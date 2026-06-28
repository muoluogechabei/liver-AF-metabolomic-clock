# ==============================================================================
# 🎯 阶段二：机器学习大逃杀 —— 构建房颤代谢风险评分 (MetRS_AF)
# ==============================================================================
library(glmnet)
library(survival)
library(dplyr)
library(doParallel)

cat("==== 1. 准备 Elastic Net 建模矩阵 (直接使用 70% 训练集) ====\n")

# 🚨 核心修正：直接使用之前物理隔离好的 merged_data_train
# 提取特征矩阵 X 和生存结局 Y (注意你的列名是 duration_updated 和 event_af_updated)
X_train_mat <- as.matrix(merged_data_train %>% select(all_of(core_metabolites)))
Y_time <- merged_data_train$duration_updated
Y_status <- merged_data_train$event_af_updated

cat(sprintf("📊 建模样本量: %d 人 | 核心池特征数: %d 个\n", nrow(merged_data_train), length(core_metabolites)))

# ==============================================================================
# 🚀 基线探路：运行全量 Cox Elastic Net (alpha = 0.5)
# ==============================================================================
cat("⏳ 正在全量训练集上运行 Cox Elastic Net，探寻最优 Lambda.1se...\n")
set.seed(2026)
# alpha = 0.5 是 Elastic Net，极其适合高共线性的脂质数据
cv_fit_af <- cv.glmnet(X_train_mat, Surv(Y_time, Y_status), 
                       family = "cox", 
                       alpha = 0.5, 
                       nfolds = 10, 
                       standardize = FALSE) # 前面已手动 scale 过了

# 提取 Lambda.1se 规则下的非零系数 (最精简、最稳健的模型)
optimal_coefs <- coef(cv_fit_af, s = "lambda.1se")
active_vars <- rownames(optimal_coefs)[which(optimal_coefs != 0)]
coef_values <- as.numeric(optimal_coefs[which(optimal_coefs != 0)])

cat("\n🎯 全量数据 Elastic Net (1se) 初步筛选完毕！入选代谢物共有:", length(active_vars), "个\n")


# ==============================================================================
# 🚀 启动 Stability Selection: Bootstrap 重抽样 (智能 1:4 采样提速)
# ==============================================================================
# 启动多核引擎 (根据你的电脑配置，可以把 3 改成你想要的核数)
n_cores <- 3
cl <- makeCluster(n_cores)
registerDoParallel(cl)
cat(sprintf("\n⚡ 已成功启动多核引擎，核心数: %d\n", n_cores))

# 准备 1:4 采样的索引
event_indices <- which(Y_status == 1)
nonevent_indices <- which(Y_status == 0)
n_events <- length(event_indices)

# ---------------------------------------------------------
# 🏁 Phase A: 前 100 次 Bootstrap 抽样
# ---------------------------------------------------------
cat("\n==== 🚀 启动 Phase A: 前 100 次 Bootstrap 抽样 ====\n")
n_boot_phase1 <- 100
results_phase1 <- list()

pb1 <- txtProgressBar(min = 0, max = n_boot_phase1, style = 3)

for (b in 1:n_boot_phase1) {
  set.seed(202600 + b) # 保持种子连续性
  
  # 1:4 智能采样逻辑 (保留所有阳性事件，随机抽取4倍的阴性对照)
  boot_idx <- c(sample(event_indices, n_events, replace = TRUE),
                sample(nonevent_indices, n_events * 4, replace = TRUE))
  
  boot_cv <- tryCatch({
    cv.glmnet(X_train_mat[boot_idx, ], Surv(Y_time[boot_idx], Y_status[boot_idx]), 
              family = "cox", alpha = 0.5, nfolds = 5, standardize = FALSE,
              parallel = TRUE, nlambda = 30)
  }, error = function(e) NULL)
  
  if (!is.null(boot_cv)) {
    boot_coefs <- coef(boot_cv, s = "lambda.1se")
    results_phase1[[b]] <- rownames(boot_coefs)[which(boot_coefs != 0)]
  }
  setTxtProgressBar(pb1, b)
}
close(pb1)

# 及时保存进度
save(results_phase1, file = "AF_Bootstrap_Checkpoint_100.RData")
cat("✅ Phase A 完成并存档。\n")

# 快速瞄一眼前 100 次的结果
phase1_summary <- as.data.frame(table(unlist(results_phase1))) %>%
  rename(Metabolite = Var1, Frequency = Freq) %>%
  mutate(Selection_Probability = Frequency / length(results_phase1) * 100) %>%
  arrange(desc(Selection_Probability))
cat("\n📊 前 100 次高频入选者预览 (Top 10):\n")
print(head(phase1_summary, 10))

# ---------------------------------------------------------
# 🏁 Phase B: 后 400 次 (101-500) 冲刺
# ---------------------------------------------------------
cat("\n==== 🚀 启动 Phase B: 后 400 次 (101-500) 冲刺 ====\n")
n_boot_phase2 <- 400
results_phase2 <- list()

pb2 <- txtProgressBar(min = 0, max = n_boot_phase2, style = 3)

for (b in 1:n_boot_phase2) {
  iter_id <- b + 100 
  set.seed(202600 + iter_id)
  
  boot_idx <- c(sample(event_indices, n_events, replace = TRUE),
                sample(nonevent_indices, n_events * 4, replace = TRUE))
  
  boot_cv <- tryCatch({
    cv.glmnet(X_train_mat[boot_idx, ], Surv(Y_time[boot_idx], Y_status[boot_idx]), 
              family = "cox", alpha = 0.5, nfolds = 5, standardize = FALSE,
              parallel = TRUE, nlambda = 30)
  }, error = function(e) NULL)
  
  if (!is.null(boot_cv)) {
    boot_coefs <- coef(boot_cv, s = "lambda.1se")
    results_phase2[[b]] <- rownames(boot_coefs)[which(boot_coefs != 0)]
  }
  setTxtProgressBar(pb2, b)
}
close(pb2)

# 释放核心
stopCluster(cl)
registerDoSEQ()

# ==============================================================================
# 🏆 最终合流：提取真金白银的核心特征
# ==============================================================================
final_boot_list <- c(results_phase1, results_phase2)
save(final_boot_list, file = "AF_Bootstrap_Final_500.RData")
cat("\n🎉 500 次抽样大功告成！全量结果已合并至 final_boot_list。\n")

# 统计最终 500 次的存活频率
final_summary <- as.data.frame(table(unlist(final_boot_list))) %>%
  rename(Metabolite = Var1, Frequency = Freq) %>%
  mutate(Selection_Probability = Frequency / length(final_boot_list)) %>%
  arrange(desc(Selection_Probability))

# 设定严苛门槛：提取出现频率 >= 80% (或 90%) 的特征
# 由于使用了 1se 和 Elastic Net，这里的保留特征会非常干净
ultra_stable_metabs <- final_summary %>% filter(Selection_Probability >= 0.80)

cat("\n======================================================\n")
cat("🏆 跨越 80% 生死线的终极王者特征共有:", nrow(ultra_stable_metabs), "个！\n")
print(ultra_stable_metabs)
cat("======================================================\n")





library(glmnet)
library(dplyr)
library(survival)

cat("\n==== 🚀 启动终极防弹版: 断点续跑 + 自动存盘 ====\n")

n_boot_phase2 <- 400
checkpoint_file <- "AF_Bootstrap_Phase2_Checkpoint.RData"

# 1. 自动检测断点：就像打游戏读档一样
if (file.exists(checkpoint_file)) {
  load(checkpoint_file)
  # 找出已经跑了多少个结果
  start_b <- length(results_phase2) + 1
  cat(sprintf("🔄 恭喜！检测到本地存档。将无缝衔接，直接从第 %d 次继续冲刺...\n", start_b))
} else {
  results_phase2 <- list()
  start_b <- 1
  cat("🆕 未检测到存档，开启全新的 400 次大逃杀...\n")
}

# 2. 开始大循环 (如果已经跑到 400，就会自动跳过)
if (start_b <= n_boot_phase2) {
  
  # 进度条起点自动适配断点
  pb2 <- txtProgressBar(min = 0, max = n_boot_phase2, initial = start_b - 1, style = 3)
  
  for (b in start_b:n_boot_phase2) {
    iter_id <- b + 100 
    set.seed(202600 + iter_id)
    
    # 1:4 抽样
    boot_idx <- c(sample(event_indices, n_events, replace = TRUE),
                  sample(nonevent_indices, n_events * 4, replace = TRUE))
    
    # 🚨 绝对保命符：关闭 parallel=TRUE，加入 maxit，保证进度条顺滑且不死锁！
    boot_cv <- tryCatch({
      cv.glmnet(X_train_mat[boot_idx, ], Surv(Y_time[boot_idx], Y_status[boot_idx]), 
                family = "cox", alpha = 0.5, nfolds = 5, standardize = FALSE,
                parallel = FALSE, nlambda = 30, maxit = 100000)
    }, error = function(e) {
      # 如果这把抽到了极端烂的数据报错了，只报个错，绝不卡死
      cat(sprintf("\n⚠️ 第 %d 次抽样遇到极端数据，已自动跳过。\n", b))
      return(NULL)
    })
    
    if (!is.null(boot_cv)) {
      boot_coefs <- coef(boot_cv, s = "lambda.1se")
      results_phase2[[b]] <- rownames(boot_coefs)[which(boot_coefs != 0)]
    } else {
      results_phase2[[b]] <- "FAILED" # 占位符，保持索引不乱
    }
    
    setTxtProgressBar(pb2, b)
    
    # 🚨 核心精髓：每跑 10 次，强制往硬盘写一次存档！天塌下来都不怕！
    if (b %% 10 == 0) {
      save(results_phase2, file = checkpoint_file)
    }
  }
  close(pb2)
  # 跑完最后保存一次
  save(results_phase2, file = checkpoint_file)
}

# ==============================================================================
# 🏆 合并与终极清算
# ==============================================================================
# 清理掉那些碰巧 FAILED 的占位符
results_phase2_clean <- results_phase2[results_phase2 != "FAILED"]
cat("\n✅ Phase B 彻底完工！有效抽样", length(results_phase2_clean), "次。\n")

# 把之前 Phase 1 (前 100次) 和现在的合在一起
load("AF_Bootstrap_Checkpoint_100.RData") 
final_boot_list <- c(results_phase1, results_phase2_clean)

# 最终保存
save(final_boot_list, file = "AF_Bootstrap_Final_500.RData")

# 计算出现频率
final_summary <- as.data.frame(table(unlist(final_boot_list))) %>%
  rename(Metabolite = Var1, Frequency = Freq) %>%
  mutate(Selection_Probability = Frequency / length(final_boot_list)) %>%
  arrange(desc(Selection_Probability))

ultra_stable_metabs <- final_summary %>% filter(Selection_Probability >= 0.80)

cat("\n======================================================\n")
cat("🏆 浴火重生的终极王者特征共有:", nrow(ultra_stable_metabs), "个！\n")
print(ultra_stable_metabs)
cat("======================================================\n")

library(glmnet)
library(dplyr)
library(survival)
library(dplyr)
library(survival)

cat("🔧 正在光速重建底层矩阵和索引...\n")
# 重建矩阵和结局
X_train_mat <- as.matrix(merged_data_train %>% select(all_of(core_metabolites)))
Y_time <- merged_data_train$duration_updated
Y_status <- merged_data_train$event_af_updated

# 重建事件索引 (这就是刚刚报错说找不到的那个家伙)
event_indices <- which(Y_status == 1)
nonevent_indices <- which(Y_status == 0)
n_events <- length(event_indices)
cat("✅ 索引重建完毕！可以开始最后冲刺了！\n")
cat("\n==== 🚀 启动终极防爆版: 断点续跑 + 暴力内存清理 ====\n")

n_boot_phase2 <- 400
checkpoint_file <- "AF_Bootstrap_Phase2_Checkpoint.RData"

# 1. 自动检测断点
if (file.exists(checkpoint_file)) {
  load(checkpoint_file)
  start_b <- length(results_phase2) + 1
  cat(sprintf("🔄 读档成功！将无缝衔接，直接从第 %d 次继续冲刺...\n", start_b))
} else {
  results_phase2 <- list()
  start_b <- 1
  cat("🆕 未检测到存档，开启全新的 400 次大逃杀...\n")
}

# 2. 开始大循环
if (start_b <= n_boot_phase2) {
  
  pb2 <- txtProgressBar(min = 0, max = n_boot_phase2, initial = start_b - 1, style = 3)
  
  for (b in start_b:n_boot_phase2) {
    iter_id <- b + 100 
    set.seed(202600 + iter_id)
    
    boot_idx <- c(sample(event_indices, n_events, replace = TRUE),
                  sample(nonevent_indices, n_events * 4, replace = TRUE))
    
    boot_cv <- tryCatch({
      cv.glmnet(X_train_mat[boot_idx, ], Surv(Y_time[boot_idx], Y_status[boot_idx]), 
                family = "cox", alpha = 0.5, nfolds = 5, standardize = FALSE,
                parallel = FALSE, nlambda = 30, maxit = 100000)
    }, error = function(e) {
      cat(sprintf("\n⚠️ 第 %d 次抽样遇到极端数据，已自动跳过。\n", b))
      return(NULL)
    })
    
    if (!is.null(boot_cv)) {
      boot_coefs <- coef(boot_cv, s = "lambda.1se")
      results_phase2[[b]] <- rownames(boot_coefs)[which(boot_coefs != 0)]
    } else {
      results_phase2[[b]] <- "FAILED" 
    }
    
    # ==========================================================
    # 🚨🚨 核心救命代码：暴力清理内存，防止撑爆电脑 🚨🚨
    # ==========================================================
    rm(boot_cv)         # 1. 彻底抹杀刚刚生成的巨大模型对象
    gc(verbose = FALSE) # 2. 拿鞭子抽 R，强制立刻回收内存给操作系统！
    # ==========================================================
    
    setTxtProgressBar(pb2, b)
    
    # 每 10 次存盘
    if (b %% 10 == 0) {
      save(results_phase2, file = checkpoint_file)
    }
  }
  close(pb2)
  save(results_phase2, file = checkpoint_file)
}

# ==============================================================================
# 🏆 合并与终极清算
# ==============================================================================
results_phase2_clean <- results_phase2[results_phase2 != "FAILED"]
cat("\n✅ Phase B 彻底完工！有效抽样", length(results_phase2_clean), "次。\n")

load("AF_Bootstrap_Checkpoint_100.RData") 
final_boot_list <- c(results_phase1, results_phase2_clean)
save(final_boot_list, file = "AF_Bootstrap_Final_500.RData")

final_summary <- as.data.frame(table(unlist(final_boot_list))) %>%
  rename(Metabolite = Var1, Frequency = Freq) %>%
  mutate(Selection_Probability = Frequency / length(final_boot_list)) %>%
  arrange(desc(Selection_Probability))

ultra_stable_metabs <- final_summary %>% filter(Selection_Probability >= 0.75)

cat("\n======================================================\n")
cat("🏆 浴火重生的终极王者特征共有:", nrow(ultra_stable_metabs), "个！\n")
print(ultra_stable_metabs)
cat("======================================================\n")
library(ggplot2)
library(dplyr)
library(stringr)
library(scales)
library(RColorBrewer) 
library(patchwork)
library(ggnewscale)

cat("⏳ 正在构建顶刊级特征权重瀑布图 (8因子带数值标出版)...\n")

# ==============================================================================
# 0. 定义全局缩写字典 (🚨 抓住了漏网之鱼)
# ==============================================================================
shorten_names <- function(x) {
  x <- gsub("\\.", " ", x) 
  # 🚨 补充的漏网之鱼：
  x <- sub("(?i).*polyunsaturated fatty acids to monounsaturated fatty acids.*", "PUFA/MUFA", x, perl=TRUE)
  
  # 原有的缩写规则
  x <- sub("(?i).*polyunsaturated fatty acids to total fatty acids.*", "PUFA% in Total FA", x, perl=TRUE)
  x <- sub("(?i).*monounsaturated fatty acids to total fatty acids.*", "MUFA% in Total FA", x, perl=TRUE)
  x <- sub("(?i).*docosahexaenoic acid to total fatty acids.*", "DHA% in Total FA", x, perl=TRUE)
  x <- gsub("(?i)Clinical LDL cholesterol", "Clinical LDL-C", x, perl=TRUE)
  x <- gsub("(?i)Glycoprotein acetyls", "GlycA", x, perl=TRUE)
  x <- gsub(" to total lipids ratio in ", "% in ", x)
  x <- gsub("Concentration of ", "", x)
  x <- gsub("Total concentration of ", "Total ", x)
  x <- gsub("Average diameter for ", "Diam: ", x)
  x <- gsub("Phospholipids", "PL", x)
  x <- gsub("Triglycerides", "TG", x)
  x <- gsub("Cholesteryl esters", "CE", x)
  x <- gsub("Free cholesterol", "FC", x)
  x <- gsub("Cholesterol", "Chol", x)
  x <- gsub("Lipoprotein", "Lipo", x)
  x <- gsub("particles", "P", x)
  x <- gsub("very large", "VL", x)
  x <- gsub("large", "L", x)
  x <- gsub("medium", "M", x)
  x <- gsub("very small", "VS", x)
  x <- gsub("small", "S", x)
  x <- gsub("Total fatty acids", "Total FA", x)
  x <- gsub("Ratio of ", "", x)
  x <- gsub("Glycolysis related metabolites", "Glycolysis", x)
  x <- gsub("Lipoprotein particle concentrations", "Lipo Concs", x)
  x <- gsub("Lipoprotein particle sizes", "Lipo Sizes", x)
  x <- gsub("Relative lipoprotein lipid concentrations", "Rel Lipo Concs", x)
  x <- gsub("Lipoprotein subclasses", "Lipo Subclasses", x)
  return(x)
}

# ==============================================================================
# 1. 构建全新数据地基 (直接使用 Cox 自由权重)
# ==============================================================================
weight_df <- data.frame(
  Metabolite = names(metRS_weights_vec),
  Cox_Beta = as.numeric(metRS_weights_vec),
  stringsAsFactors = FALSE
)

plot_weight_df <- weight_df %>%
  mutate(
    Raw_Name = gsub("\\.", " ", Metabolite),
    Display_Name = shorten_names(Raw_Name)
  ) %>%
  left_join(metab_dict %>% select(title, Group) %>% mutate(title = shorten_names(title)) %>% distinct(),
            by = c("Display_Name" = "title")) %>%
  mutate(
    Group = ifelse(is.na(Group) | Group == "", "Other Metabolites", Group),
    Group_Short = shorten_names(Group),
    Abs_Beta = abs(Cox_Beta)
  ) %>%
  arrange(Group_Short, desc(Abs_Beta)) %>%
  mutate(
    Display_Name = factor(Display_Name, levels = unique(Display_Name)),
    Group_Short = factor(Group_Short, levels = unique(Group_Short))
  )

# ==============================================================================
# 2. 同源色彩引擎 
# ==============================================================================
groups_short <- levels(plot_weight_df$Group_Short)
n_groups <- length(groups_short)
expanded_colors <- colorRampPalette(brewer.pal(8, "Set2"))(n_groups)
group_colors <- setNames(expanded_colors, groups_short)

# ==============================================================================
# 3. 双核驱动绘图：无缝色带 + 呼吸柱状图
# ==============================================================================

# 图 1：左侧色带图 (无缝衔接)
p_color_band <- ggplot(plot_weight_df, aes(x = 1, y = Display_Name)) +
  geom_tile(aes(fill = Group_Short)) +
  scale_fill_manual(values = group_colors, guide = "none") +
  facet_grid(Group_Short ~ ., scales = "free_y", space = "free_y") +
  theme_void() +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    panel.spacing = unit(0, "lines"), 
    plot.margin = margin(0, 0, 0, 0)
  )

# 图 2：右侧主图 (🚨 已添加内置具体数值标签)
p_main <- ggplot(plot_weight_df, aes(x = Cox_Beta, y = Display_Name)) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.8) +
  
  # 柱子本体
  geom_bar(aes(fill = Cox_Beta), stat = "identity", width = 0.65, color = NA) +
  
  # 🚨 新增：白色内置数值标签！
  # 动态判断正负，使得标签永远向内缩 (hjust 1.2表示向左缩，-0.2表示向右缩)
  geom_text(aes(label = sprintf("%.2f", Cox_Beta), 
                hjust = ifelse(Cox_Beta > 0, 1.2, -0.2)), 
            color = "white", size = 3.5, fontface = "bold") +
  
  scale_fill_gradient2(
    name = "Cox Coefficient (\u03b2)",
    low = "#313695", mid = "grey95", high = "#A50026", midpoint = 0,
    guide = guide_colorbar(barwidth = 1.5, barheight = 10, frame.colour = "black", ticks.colour = "black", order = 1)
  ) +
  
  ggnewscale::new_scale_fill() +
  geom_point(aes(fill = Group_Short), x = 0, shape = 22, size = 0, color = NA) +
  scale_fill_manual(
    name = "Metabolic Super-class", 
    values = group_colors,
    guide = guide_legend(override.aes = list(size = 5, color = "grey30"), order = 2)
  ) +
  
  facet_grid(Group_Short ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1)), breaks = scales::pretty_breaks(n = 5)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 12)),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x = element_line(color = "black"),
    
    axis.text.y = element_text(size = 11, color = "black", face = "bold", margin = margin(r = 5)),
    axis.title.y = element_blank(),
    
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey85", linetype = "dashed"),
    
    strip.text = element_blank(), 
    strip.background = element_blank(),
    panel.spacing = unit(0.4, "lines"), 
    
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.spacing.y = unit(0.3, "cm"), 
    
    plot.margin = margin(0, 0, 0, 0)
  )

# 终极拼合
final_plot <- p_color_band + p_main + 
  plot_layout(widths = c(0.02, 1)) + 
  plot_annotation(
    title = "Metabolic Signatures for Incident AF",
    subtitle = "Unshrunken Cox Regression coefficients (\u03b2) of the 8 core metabolites",
    theme = theme(
      plot.title = element_text(size = 18, face = "bold", hjust = 0, margin = margin(b = 6)),
      plot.subtitle = element_text(size = 12, color = "grey40", hjust = 0, margin = margin(b = 20)),
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
      plot.background = element_rect(fill = "white", color = NA)
    )
  )

# 4. 导出
ggsave("Figure_Metabolite_Weights_8Factors_With_Values.pdf", plot = final_plot, width = 8, height = 8, dpi = 400)

cat("✅ 终极完美权重图已生成：Figure_Metabolite_Weights_8Factors_With_Values.pdf\n")



library(ggplot2)
library(dplyr)
library(stringr)
library(RColorBrewer)

cat("\n==== 🎯 启动模型验证可视化双引擎 ====\n")

# ==============================================================================
# 🎨 1. 生成图 A & B：Elastic Net 系数路径与 CV 误差拼图
# ==============================================================================
cat("⏳ 正在绘制 Elastic Net 系数路径与 CV 误差图...\n")

# 提取 cv_fit_af 的底层数据
log_lambdas <- log(cv_fit_af$glmnet.fit$lambda)     
coef_matrix <- t(as.matrix(cv_fit_af$glmnet.fit$beta)) 
# 🚨 你的模型筛选用的是 1se，所以这里红线对齐 1se 最严谨
log_lambda_1se <- log(cv_fit_af$lambda.1se)         

# 设置高级配色
premium_colors <- colorRampPalette(brewer.pal(8, "Set2"))(ncol(coef_matrix))

pdf("Figure_S1_ElasticNet_Path_and_CV.pdf", width = 12, height = 6)
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2) + 0.1)

# ---------------------------------------------------------
# 图 A: 系数路径图
# ---------------------------------------------------------
matplot(log_lambdas, coef_matrix, 
        type = "l", lty = 1, lwd = 2, col = premium_colors,
        las = 1, 
        xlab = expression(bold(Log(lambda))),
        ylab = "Elastic Net Coefficients (\u03b1 = 0.5)")

mtext("A", side = 3, adj = -0.15, line = 1.5, cex = 2, font = 2)

# 在 1se 位置画红色虚线
abline(v = log_lambda_1se, col = "#E64B35FF", lty = 2, lwd = 2.5)

# 顶部添加变量个数
n_vars <- cv_fit_af$glmnet.fit$df
axis(side = 3, at = log_lambdas[seq(1, length(log_lambdas), length.out = 5)], 
     labels = n_vars[seq(1, length(log_lambdas), length.out = 5)], 
     tick = FALSE, line = -0.5)

# ---------------------------------------------------------
# 图 B: CV 误差图
# ---------------------------------------------------------
plot(cv_fit_af, 
     sign.lambda = 1, 
     las = 1, 
     lwd = 1.5, pch = 19, cex = 0.9,
     xlab = expression(bold(Log(lambda))),
     ylab = "Partial Likelihood Deviance")

mtext("B", side = 3, adj = -0.15, line = 1.5, cex = 2, font = 2)

# 强化 1se 参考线
abline(v = log_lambda_1se, col = "#E64B35FF", lty = 2, lwd = 2.5)

legend("topleft", legend = "Optimal Lambda (1se)", 
       col = "#E64B35FF", lty = 2, lwd = 2.5, bty = "n", cex = 1.1)

dev.off()
cat("✅ 图 1 完成！系数路径与 CV 误差图已保存至 Figure_S1_ElasticNet_Path_and_CV.pdf\n")


# ==============================================================================
# 🎨 2. 生成图 C：只展示 >0 的 Bootstrap 存活概率柱状图
# ==============================================================================
cat("\n⏳ 正在绘制 Bootstrap 核心变量突围柱状图...\n")

# 定义自动缩写引擎 (复用你最完美的那个)
shorten_names <- function(x) {
  x <- gsub("\\.", " ", x) 
  x <- sub("(?i).*polyunsaturated fatty acids to monounsaturated fatty acids.*", "PUFA/MUFA", x, perl=TRUE)
  x <- sub("(?i).*polyunsaturated fatty acids to total fatty acids.*", "PUFA% in Total FA", x, perl=TRUE)
  x <- sub("(?i).*monounsaturated fatty acids to total fatty acids.*", "MUFA% in Total FA", x, perl=TRUE)
  x <- sub("(?i).*docosahexaenoic acid to total fatty acids.*", "DHA% in Total FA", x, perl=TRUE)
  x <- gsub("(?i)Clinical LDL cholesterol", "Clinical LDL-C", x, perl=TRUE)
  x <- gsub("(?i)Glycoprotein acetyls", "GlycA", x, perl=TRUE)
  x <- gsub(" to total lipids ratio in ", "% in ", x)
  x <- gsub("Concentration of ", "", x)
  x <- gsub("Total concentration of ", "Total ", x)
  x <- gsub("Average diameter for ", "Diam: ", x)
  x <- gsub("Phospholipids", "PL", x)
  x <- gsub("Triglycerides", "TG", x)
  x <- gsub("Cholesteryl esters", "CE", x)
  x <- gsub("Free cholesterol", "FC", x)
  x <- gsub("Cholesterol", "Chol", x)
  x <- gsub("Lipoprotein", "Lipo", x)
  x <- gsub("particles", "P", x)
  x <- gsub("very large", "VL", x)
  x <- gsub("large", "L", x)
  x <- gsub("medium", "M", x)
  x <- gsub("very small", "VS", x)
  x <- gsub("small", "S", x)
  x <- gsub("Total fatty acids", "Total FA", x)
  x <- gsub("Ratio of ", "", x)
  x <- gsub("Glycolysis related metabolites", "Glycolysis", x)
  x <- gsub("Lipoprotein particle concentrations", "Lipo Concs", x)
  x <- gsub("Lipoprotein particle sizes", "Lipo Sizes", x)
  x <- gsub("Relative lipoprotein lipid concentrations", "Rel Lipo Concs", x)
  x <- gsub("Lipoprotein subclasses", "Lipo Subclasses", x)
  return(x)
}

# ---------------------------------------------------------
# 数据清洗与门槛划定 (只取 > 0 的，画 80% 红线)
# ---------------------------------------------------------
# 注意：你前面的 final_summary 算出来的是 0~1 的小数，所以这里 * 100 转成百分比
plot_prob_df <- final_summary %>%
  mutate(Prob_Pct = Selection_Probability * 100) %>%
  filter(Prob_Pct > 10) %>% # 🚨 核心要求：只要概率大于0的！
  mutate(
    short_name = shorten_names(as.character(Metabolite)),
    # 这里阈值设为 80% (与你上面的 ultra_stable_metabs 筛选一致)
    Status = ifelse(Prob_Pct >= 75, "Elite Retained (\u2265 75%)", "Eliminated (< 75%)") 
  ) %>%
  arrange(Prob_Pct, short_name) %>%
  mutate(short_name = factor(short_name, levels = unique(short_name)))

cat(sprintf("🔍 画图变量统计：共有 %d 个变量在 500 次重抽样中至少出现过 1 次。\n", nrow(plot_prob_df)))

# 动态字体颜色渲染
ordered_levels <- levels(plot_prob_df$short_name)
color_map <- ifelse(plot_prob_df$Status[match(ordered_levels, plot_prob_df$short_name)] == "Elite Retained (\u2265 75%)", "black", "grey60")
face_map <- ifelse(plot_prob_df$Status[match(ordered_levels, plot_prob_df$short_name)] == "Elite Retained (\u2265 75%)", "bold", "plain")

# ---------------------------------------------------------
# ggplot 绘制柱状图
# ---------------------------------------------------------
p_prob <- ggplot(plot_prob_df, aes(x = Prob_Pct, y = short_name, fill = Status)) +
  geom_col(width = 0.75, alpha = 0.9) +
  
  # 🚨 80% 的生死红线 (你可以自行改成 75)
  geom_vline(xintercept = 75, linetype = "dashed", color = "#C9184A", linewidth = 1.2) +
  
  scale_fill_manual(values = c("Elite Retained (\u2265 75%)" = "#E64B35FF", "Eliminated (< 75%)" = "#DDDDDD")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, 105), breaks = seq(0, 100, 20)) +
  
  labs(
    x = "Bootstrap Selection Probability (%)",
    y = NULL,
    title = "Stability Selection via Elastic Net",
    subtitle = "500 Bootstrap iterations (Showing only features with probability > 0)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 12, color = "grey30", margin = margin(b = 15)),
    
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    axis.text.x = element_text(size = 11, color = "black"),
    
    # 自动高亮晋级选手的文字
    axis.text.y = element_text(size = 9, color = color_map, face = face_map),
    
    legend.position = c(0.75, 0.2), # 图例放在右下角空白处
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    legend.background = element_rect(fill = alpha("white", 0.8), color = "grey80", linewidth = 0.5),
    
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7)
  )

# 根据留下来的变量数量动态调整画幅高度 (稍微紧凑一点，不至于太空)
dynamic_height <- max(5, nrow(plot_prob_df) * 0.15 + 2)

ggsave("Figure_S2_Bootstrap_Selection_Filtered.pdf", plot = p_prob, width = 9, height = dynamic_height, bg = "white")

cat(sprintf("✅ 图 2 完成！已剔除概率为 0 的变量，画布高度自适应为 %.1f，文件保存在 Figure_S2_Bootstrap_Selection_Filtered.pdf\n", dynamic_height))
library(survival)
library(dplyr)
library(broom)

cat("\n==== 🚀 最终章：构建并合并 8 因子代谢风险评分 (MetRS) ====\n")

# 1. 提取 8个王者代谢物
final_11_metabs <- ultra_stable_metabs$Metabolite

# 2. 从 NMR 大表中提取
nmr_11_data <- my_data_nmr %>%
  select(eid, all_of(final_11_metabs))

# 3. 合并到 benchmark 表
df_benchmark_scores <- df_benchmark_scores %>%
  # 如果之前跑过这一步，先去重防报错
  select(-any_of(final_11_metabs)) %>% 
  inner_join(nmr_11_data, by = "eid")

cat("✅ 8王者代谢物已成功并入临床评分主表！\n")

# ==============================================================================
# 🚨 修复核心：直接从 df_clinical_final 提取生存数据打补丁
# ==============================================================================
cat("🔧 正在打补丁：补齐生存随访时间和结局...\n")

survival_data <- df_clinical_final %>%
  select(eid = `Participant ID`, duration_updated, event_af_updated)

library(glmnet)

cat("\n==== 🚀 启动岭回归 (Ridge Cox) 获取抗共线性权重 ====\n")

# ==============================================================================
# 4. 在“训练集 (70%)”中使用标准 Cox 回归提取“自由权重”
# ==============================================================================
# 🚨 核心逻辑微调：既然特征已缩减至 0.75 频率以上，直接用 Unshrunken Cox 释放预测潜能
cat(sprintf("⏳ 正在基于 %d 名患者拟合标准 Cox 模型 (释放自由权重)...\n", nrow(df_train_final)))

# 自动构建 0.75 阈值筛选出的特征公式
f_metRS <- as.formula(paste("Surv(duration_updated, event_af_updated) ~", 
                            paste(final_11_metabs, collapse = " + ")))

# 使用普通 Cox 拟合，不再进行 L2 惩罚
fit_unshrunken <- coxph(f_metRS, data = df_train_final)

# 提取自由回归系数 (Beta)，变量名保持不变，直接衔接你后面的代码
metRS_weights_vec <- coef(fit_unshrunken)

# 打印权重表确认，你会发现系数比 Ridge 时代更具有“侵略性”
cat("\n🏆 【最终 MetRS 评分自由权重表 (Unshrunken)】 🏆\n")
print(data.frame(Metabolite = names(metRS_weights_vec), Free_Weight = metRS_weights_vec) %>% 
        arrange(desc(abs(Free_Weight))))
# ==============================================================================
# 5. 在全量数据中计算每个人的 MetRS 绝对得分
# ==============================================================================
cat("\n⏳ 正在为全队列重新计算 MetRS 得分...\n")

# 提取全队列的 11 个特征矩阵 (注意列名顺序要和权重一致)
X_all_11 <- as.matrix(df_benchmark_scores %>% select(all_of(names(metRS_weights_vec))))

# 矩阵乘法：得分 = 特征1*权重1 + 特征2*权重2 + ...
df_benchmark_scores$MetRS_Raw <- as.numeric(X_all_11 %*% metRS_weights_vec)

# 🚨 Z-score 标准化
df_benchmark_scores <- df_benchmark_scores %>%
  mutate(MetRS_Z = as.numeric(scale(MetRS_Raw)))

cat("🎉 岭回归打分完毕！现在这个 MetRS_Z 的稳健性绝对无敌了。\n")

# 6. 极速检验：重新看看分布
p_metrs_ridge <- ggplot(df_benchmark_scores, aes(x = MetRS_Z)) +
  geom_density(fill = "#1b9e77", alpha = 0.6) +
  theme_minimal() + 
  labs(title = "MetRS Distribution (Ridge-Derived Weights)",
       subtitle = "Z-score standardized, robust against collinearity",
       x = "MetRS (SD)", y = "Density")

print(p_metrs_ridge)








cat("\n==== 📦 开始进行终极数据打包与大保存 ====\n")

# ==============================================================================
# 1. 终极表合并与 7:3 物理切割
# ==============================================================================
cat("⏳ 正在将全队列划分为 Training Set (70%) 和 Testing Set (30%)...\n")

# 确保全量评分表中包含了生存随访时间 (防止之前 left_join 遗漏)
# 使用 select(-any_of(...)) 防止重复合并产生 .x / .y 后缀
df_benchmark_scores <- df_benchmark_scores %>%
  select(-any_of(c("duration_updated", "event_af_updated"))) %>% 
  left_join(survival_data, by = "eid") %>%
  # 剔除缺失随访信息的极少数样本
  drop_na(duration_updated, event_af_updated)

# 🗡️ 极其神圣的一刀：按之前锁定的 eid 切割
df_train_final <- df_benchmark_scores %>% filter(eid %in% train_eids_locked)
df_test_final  <- df_benchmark_scores %>% filter(!(eid %in% train_eids_locked))

cat("✅ 数据集物理切割完毕！\n")
cat(sprintf("   ▶ 训练集 (Training Set) 样本量: %d 人 (用于筛选和算权重)\n", nrow(df_train_final)))
cat(sprintf("   ▶ 测试集 (Testing Set)  样本量: %d 人 (极其纯洁，仅用于验证)\n", nrow(df_test_final)))

# ==============================================================================
# 2. 执行大保存 (The Grand Save)
# ==============================================================================
cat("\n⏳ 正在将所有核心资产落盘封印...\n")

# 设定终极存档文件名
final_archive_name <- "AF_MetRS_Final_Master_Data_v1.RData"

save(
  # 📊 1. 核心数据集 (重中之重)
  df_benchmark_scores,   # 全队列临床+代谢评分大表
  df_train_final,        # 70% 训练集
  df_test_final,         # 30% 测试集
  df_ready_for_scores,   # 最早清洗的纯临床数据表 (备用)
  survival_data,         # 生存结局表
  
  # 🧠 2. 模型与权重
  final_11_metabs,       # 11 个王者代谢物名单
  metRS_weights_vec,     # 岭回归平滑处理后的终极权重向量
  ultra_stable_metabs,   # 500次重抽样总结表 (可用来画特征纳入频率图)
  ridge_fit,             # 岭回归模型本体
  
  # 🔑 3. 关键索引与元数据
  train_eids_locked,     # 确保 7:3 划分绝对不可篡改的 ID 锁
  core_metabolites,      # 144 个 FIB-4 & NFS 的交集名单
  
  file = final_archive_name
)

cat("🎉 伟大的胜利！所有重要数据已成功封印至本地文件：\n")
cat(sprintf("   👉 [%s]\n", final_archive_name))
cat("======================================================\n")
cat("💡 接下来，你可以直接开一个新的空白 R 脚本，\n")
cat("   运行 load(\"", final_archive_name, "\")，即可开始画生存图和算 cNRI！\n", sep="")
cat("======================================================\n")


library(dplyr)
library(survival)
library(timeROC)
library(ggplot2)
library(tidyr)

cat("\n==== 🌟 启动终极增量价值分析 (Incremental Value Analysis) ====\n")

# ==============================================================================
# 1. 补齐 PRS 数据并清理测试集
# ==============================================================================
cat("⏳ 正在为测试集无缝拼装遗传多基因风险评分 (PRS)...\n")

# 1. 提取 PRS 数据 (保持原始列名不动)
prs_data <- df_clinical_final %>%
  select(`Participant ID`, PRS_AF = `Standard PRS for atrial fibrillation (AF)`)

# 2. 核心修正：使用 c("左表列名" = "右表列名") 的语法进行跨列名对齐
test_eval_final <- df_test_final %>%
  # 🚨 关键在这里：告诉 R 左边的 eid 对应右边的 Participant ID
  left_join(prs_data, by = c("eid" = "Participant ID")) %>%
  # 剔除缺失值
  drop_na(duration_updated, event_af_updated, PRS_AF, MetRS_Z)

cat(sprintf("✅ 物理隔离后的测试集拼装完成！\n样本量: %d 人 | 包含变量: eid, 临床评分, MetRS_Z, PRS_AF\n", nrow(test_eval_final)))
library(dplyr)
library(survival)
library(timeROC)
library(ggplot2)
library(patchwork)
library(tidyr)

cat("\n==== 🚀 拨乱反正！基于独立验证集的真正增量分析 ====\n")

# 在拼装测试集的地方，确保把 LP 字段选进去
test_eval_final <- df_test_final %>%
  select(eid, duration_updated, event_af_updated, MetRS_Z, any_of("PRS_AF"),
         CHARGE_AF_LP, ARIC_LP_centered, C2HEST_Points) %>% # 👈 确保这几行在 select 里
  left_join(prs_data, by = c("eid" = "Participant ID")) %>%
  drop_na()

train_eval_final <- df_train_final %>%
  select(eid, duration_updated, event_af_updated, MetRS_Z, any_of("PRS_AF"),
         CHARGE_AF_LP, ARIC_LP_centered, C2HEST_Points) %>% # 👈 确保这几行在 select 里
  left_join(prs_data, by = c("eid" = "Participant ID")) %>%
  drop_na()
cat(sprintf("✅ 拼装完毕！训练集: %d 人 | 测试集: %d 人\n", nrow(train_eval_final), nrow(test_eval_final)))

# ==============================================================================
# 2. 核心引擎：在 Train 拟合权重，在 Test 独立验证
# ==============================================================================
get_true_nested_roc <- function(train_data, test_data, base_var, t_eval = 5, n_boot = 500) {
  
  # A. 准备嵌套公式
  f_base <- as.formula(paste("Surv(duration_updated, event_af_updated) ~", base_var))
  f_prs  <- as.formula(paste("Surv(duration_updated, event_af_updated) ~", base_var, "+ PRS_AF"))
  f_mrs  <- as.formula(paste("Surv(duration_updated, event_af_updated) ~", base_var, "+ MetRS_Z"))
  f_full <- as.formula(paste("Surv(duration_updated, event_af_updated) ~", base_var, "+ PRS_AF + MetRS_Z"))
  
  # B. 🚨 核心纠错：在 TRAIN SET 上拟合模型，锁定最优权重
  fit_base <- coxph(f_base, data = train_data)
  fit_prs  <- coxph(f_prs,  data = train_data)
  fit_mrs  <- coxph(f_mrs,  data = train_data)
  fit_full <- coxph(f_full, data = train_data)
  
  # C. 🚨 核心纠错：在 TEST SET 上进行纯净的打分 (不再有任何作弊)
  m_base <- predict(fit_base, newdata = test_data)
  m_prs  <- predict(fit_prs,  newdata = test_data)
  m_mrs  <- predict(fit_mrs,  newdata = test_data)
  m_full <- predict(fit_full, newdata = test_data)
  
  # 计算 Test Set 真实 AUC
  roc_base <- timeROC(T=test_data$duration_updated, delta=test_data$event_af_updated, marker=m_base, cause=1, times=t_eval, iid=FALSE)
  roc_prs  <- timeROC(T=test_data$duration_updated, delta=test_data$event_af_updated, marker=m_prs,  cause=1, times=t_eval, iid=FALSE)
  roc_mrs  <- timeROC(T=test_data$duration_updated, delta=test_data$event_af_updated, marker=m_mrs,  cause=1, times=t_eval, iid=FALSE)
  roc_full <- timeROC(T=test_data$duration_updated, delta=test_data$event_af_updated, marker=m_full, cause=1, times=t_eval, iid=FALSE)
  
  # 定位时间为 t_eval (5年) 的索引
  idx <- which(roc_base$times == t_eval)
  if(length(idx)==0) idx <- length(roc_base$times)
  
  # D. 飞速 Bootstrap 抽样计算 CI 和 P 值
  # 因为直接抽预测得分，500次只需几秒钟！
  boot_res <- matrix(NA, nrow = n_boot, ncol = 4)
  colnames(boot_res) <- c("Base", "Base_PRS", "Base_MRS", "Full")
  
  cat(sprintf("   🏃 正在评估基线 [%s] (500次极速 Bootstrap)...\n", base_var))
  
  # 提前把列拿出来，极大提升 for 循环速度
  T_vec <- test_data$duration_updated
  D_vec <- test_data$event_af_updated
  
  for(b in 1:n_boot) {
    set.seed(202604 + b)
    boot_idx <- sample(1:nrow(test_data), replace = TRUE)
    
    # 极速计算 AUC (只需时间，状态，和提取好的 marker)
    b_r1 <- timeROC(T=T_vec[boot_idx], delta=D_vec[boot_idx], marker=m_base[boot_idx], cause=1, times=t_eval, iid=FALSE)$AUC[idx]
    b_r2 <- timeROC(T=T_vec[boot_idx], delta=D_vec[boot_idx], marker=m_prs[boot_idx],  cause=1, times=t_eval, iid=FALSE)$AUC[idx]
    b_r3 <- timeROC(T=T_vec[boot_idx], delta=D_vec[boot_idx], marker=m_mrs[boot_idx],  cause=1, times=t_eval, iid=FALSE)$AUC[idx]
    b_r4 <- timeROC(T=T_vec[boot_idx], delta=D_vec[boot_idx], marker=m_full[boot_idx], cause=1, times=t_eval, iid=FALSE)$AUC[idx]
    
    boot_res[b, ] <- c(b_r1, b_r2, b_r3, b_r4)
  }
  
  # E. 计算 CI 和 3个核心 P 值
  get_ci <- function(vec) {
    sprintf("(%.3f-%.3f)", quantile(vec, 0.025, na.rm=TRUE), quantile(vec, 0.975, na.rm=TRUE))
  }
  
  calc_p <- function(v_new, v_old) {
    p <- sum((v_new - v_old) <= 0, na.rm=TRUE) / n_boot
    if(p == 0) return(sprintf("P < %.3f", 1/n_boot)) else return(sprintf("P = %.3f", p))
  }
  
  p_prs_vs_base  <- calc_p(boot_res[,"Base_PRS"], boot_res[,"Base"])     # 基因 vs 临床
  p_mrs_vs_base  <- calc_p(boot_res[,"Base_MRS"], boot_res[,"Base"])     # 代谢 vs 临床
  p_full_vs_prs  <- calc_p(boot_res[,"Full"],     boot_res[,"Base_PRS"]) # 全模型 vs 基因
  
  # F. 组装作图数据 (直接把 P 值贴在字后面)
  build_df <- function(roc_obj, name, boot_vec, p_str) {
    label_str <- sprintf("%s: %.3f %s", name, roc_obj$AUC[idx], get_ci(boot_vec))
    if (p_str != "") label_str <- paste0(label_str, " | ", p_str)
    
    data.frame(
      Spec = 1 - roc_obj$FP[, idx],
      Sens = roc_obj$TP[, idx],
      Model = name,
      AUC_Label = label_str
    )
  }
  
  plot_df <- rbind(
    build_df(roc_base, "1. Clinical Alone",         boot_res[,"Base"],     ""),
    build_df(roc_prs,  "2. + Genetics (PRS)",       boot_res[,"Base_PRS"], paste0(p_prs_vs_base, " vs 1")),
    build_df(roc_mrs,  "3. + Metabolomics (MetRS)", boot_res[,"Base_MRS"], paste0(p_mrs_vs_base, " vs 1")),
    build_df(roc_full, "4. Full Model (Both)",      boot_res[,"Full"],     paste0(p_full_vs_prs, " vs 2"))
  )
  
  return(list(df = plot_df)) # 👈 注意这里不需要单独返回 p_text 了
}
library(stringr) # 确保加载了正则包

# ==============================================================================
# 3. 绘图函数 (剔除置信区间，让曲线图更清爽)
# ==============================================================================
draw_nested_roc <- function(res, title) {
  df <- res$df
  df$Model <- factor(df$Model, levels = c("4. Full Model (Both)", "3. + Metabolomics (MetRS)", 
                                          "2. + Genetics (PRS)", "1. Clinical Alone"))
  
  # 按照 Model 因子顺序提取标签
  labels_vec <- unique(df$AUC_Label[order(df$Model)])
  
  # 🚨 核心修改：使用正则表达式，精准剔除 " (0.xxx-0.xxx)" 部分，保留 AUC 和 P值
  labels_vec <- gsub(" \\([0-9.]+\\-[0-9.]+\\)", "", labels_vec)
  
  ggplot(df, aes(x = Spec, y = Sens, color = Model)) +
    geom_abline(intercept = 1, slope = 1, color = "gray80", linetype = "dashed") +
    geom_line(aes(linewidth = Model)) +
    scale_linewidth_manual(values = c("4. Full Model (Both)" = 1.2, "3. + Metabolomics (MetRS)" = 0.8, 
                                      "2. + Genetics (PRS)" = 0.8, "1. Clinical Alone" = 0.8), guide = "none") +
    scale_color_manual(values = c("4. Full Model (Both)" = "#D62728", 
                                  "3. + Metabolomics (MetRS)" = "#1F77B4", 
                                  "2. + Genetics (PRS)" = "#FF7F0E", 
                                  "1. Clinical Alone" = "#7F7F7F")) + 
    # 翻转 X 轴
    scale_x_reverse(limits = c(1, 0), breaks = seq(1, 0, -0.2), expand = c(0.01, 0.01)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
    labs(title = title, x = "1 - Specificity", y = "Sensitivity") +
    theme_bw(base_size = 14) + # 整体基础字号调大
    theme(legend.position = "none", 
          plot.title = element_text(face="bold", size=16, margin = margin(b=10)),
          panel.grid.minor = element_blank()) +
    
    # 🚨 标注文字位置对齐 (x = 0.99 保证起点靠左)
    annotate("text", x = 0.89, y = 0.32, label = labels_vec[1], color = "#D62728", fontface="bold", hjust=0, size=4.5) +
    annotate("text", x = 0.89, y = 0.24, label = labels_vec[2], color = "#1F77B4", fontface="bold", hjust=0, size=4.5) +
    annotate("text", x = 0.89, y = 0.16, label = labels_vec[3], color = "#FF7F0E", fontface="bold", hjust=0, size=4.5) +
    annotate("text", x = 0.89, y = 0.08, label = labels_vec[4], color = "#7F7F7F", fontface="bold", hjust=0, size=4.5)
}
# ==============================================================================
# 4. 自动化流水线执行
# ==============================================================================
cat("\n--- 开始计算真实增量价值 (秒级 500 次跑完) ---\n")
# ✅ 必须改为原始线性预测值 (Linear Predictor)：
res_charge <- get_true_nested_roc(train_eval_final, test_eval_final, "CHARGE_AF_LP", n_boot = 500)
res_aric   <- get_true_nested_roc(train_eval_final, test_eval_final, "ARIC_LP_centered", n_boot = 500)
res_c2hest <- get_true_nested_roc(train_eval_final, test_eval_final, "C2HEST_Points", n_boot = 500)
p1 <- draw_nested_roc(res_charge, "A. Baseline: CHARGE-AF")
p2 <- draw_nested_roc(res_aric,   "B. Baseline: ARIC")
p3 <- draw_nested_roc(res_c2hest, "C. Baseline: C2HEST")
final_fig4 <- p1 + p2 + p3 + plot_layout(ncol = 3)
ggsave("Figure_True_Nested_ROC.pdf", final_fig4, width = 18, height = 6.5, dpi = 400)

cat("\n==== 📊 正在生成附带 95% 置信区间 (误差棒+数值) 的阶梯增量柱状图 ====\n")

# 1. 提取数据的函数保持不变
extract_bar_data <- function(res_obj, baseline_name) {
  res_obj$df %>%
    select(Model, AUC_Label) %>%
    distinct() %>%
    mutate(
      Baseline = baseline_name,
      Step = Model,
      AUC = as.numeric(str_extract(AUC_Label, "(?<=: )[0-9.]+")),
      CI_Low = as.numeric(str_extract(AUC_Label, "(?<=\\()[0-9.]+(?=-)")),
      CI_High = as.numeric(str_extract(AUC_Label, "(?<=-)[0-9.]+(?=\\))"))
    )
}

# 2. 合并三大基线的数据
plot_df_bar <- bind_rows(
  extract_bar_data(res_charge, "CHARGE-AF"),
  extract_bar_data(res_aric, "ARIC"),
  extract_bar_data(res_c2hest, "C2HEST")
)

plot_df_bar$Step <- factor(plot_df_bar$Step, levels = c("1. Clinical Alone", "2. + Genetics (PRS)", "3. + Metabolomics (MetRS)", "4. Full Model (Both)"))
plot_df_bar$Baseline <- factor(plot_df_bar$Baseline, levels = c("ARIC", "C2HEST", "CHARGE-AF"))

# 3. 🚨 动态计算 Y 轴范围：把上限 + 0.05（原来是0.03），给两行文字留足天空！
ymin <- min(plot_df_bar$CI_Low, na.rm = TRUE) - 0.02
ymax <- max(plot_df_bar$CI_High, na.rm = TRUE) + 0.05

# 4. 绘制顶级期刊标准增量柱状图
p_bar_final <- ggplot(plot_df_bar, aes(x = Step, y = AUC, fill = Step)) +
  # 柱子本体
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7, alpha = 0.9, color = "black", linewidth = 0.3) +
  # 添加置信区间误差棒 (Error Bars)
  geom_errorbar(aes(ymin = CI_Low, ymax = CI_High), width = 0.15, position = position_dodge(0.9), color = "black", alpha = 0.8, linewidth = 0.6) +
  
  # 🚨 核心排版魔法：使用 \n 换行，上面是 AUC，下面是 (CI_Low-CI_High)，缩小一点字号防拥挤
  geom_text(aes(label = sprintf("%.3f\n(%.3f-%.3f)", AUC, CI_Low, CI_High), y = CI_High + 0.003), 
            position = position_dodge(0.9), vjust = 0, fontface = "bold", size = 3.2, lineheight = 0.8) +
  
  facet_wrap(~Baseline) +
  scale_fill_manual(values = c("1. Clinical Alone" = "#999999", 
                               "2. + Genetics (PRS)" = "#4DBBD5", 
                               "3. + Metabolomics (MetRS)" = "#E64B35", 
                               "4. Full Model (Both)" = "#3C5488")) +
  coord_cartesian(ylim = c(ymin, ymax)) +
  labs(
    title = "Incremental Predictive Value of MetRS Across Clinical Benchmarks",
    subtitle = "5-Year AF Prediction in Independent Test Cohort (Showing 95% CI)",
    y = "Area Under the Curve (AUC)", x = ""
  ) +
  theme_bw(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "grey20"),
    strip.text = element_text(color = "white", face = "bold", size = 13),
    axis.text.x = element_blank(), # 隐藏底部文字，靠图例
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

print(p_bar_final)
ggsave("Figure_MRS_Value_Comparison_with_CI.pdf", p_bar_final, width = 14, height = 7.5, dpi = 400) # 画幅稍微加宽加大一点

cat("✅ 柱状图生成完毕！快去查看 Figure_MRS_Value_Comparison_with_CI.pdf\n")


library(dplyr)
library(survival)
library(ggplot2)
library(patchwork)

cat("\n==== 🎯 启动终极武器：四重嵌套校准曲线 (Nested Calibration Curves) ====\n")

# ==============================================================================
# 1. 核心计算引擎：获取 4 个模型的预测概率与实际发生率
# ==============================================================================
get_4_calib_data <- function(train_data, test_data, base_var, t_eval_years = 5, n_boot = 500) {
  
  # A. 时间单位自适应
  max_t <- max(test_data$duration_updated, na.rm = TRUE)
  if (max_t > 100) {
    t_eval <- t_eval_years * 365.25
    cat(sprintf("📏 自动探测：时间单位为 [天]，校准点为 %.1f 天\n", t_eval))
  } else {
    t_eval <- t_eval_years
    cat(sprintf("📏 自动探测：时间单位为 [年]，校准点为 %d 年\n", t_eval))
  }
  
  # B. 准备公式
  f_list <- list(
    "1. Clinical Alone"         = as.formula(paste("Surv(duration_updated, event_af_updated) ~", base_var)),
    "2. + Genetics (PRS)"       = as.formula(paste("Surv(duration_updated, event_af_updated) ~", base_var, "+ PRS_AF")),
    "3. + Metabolomics (MetRS)" = as.formula(paste("Surv(duration_updated, event_af_updated) ~", base_var, "+ MetRS_Z")),
    "4. Full Model (Both)"      = as.formula(paste("Surv(duration_updated, event_af_updated) ~", base_var, "+ PRS_AF + MetRS_Z"))
  )
  
  cat(sprintf("⏳ 正在为基线 [%s] 拟合模型并计算 5 年发病概率...\n", base_var))
  
  # C. 计算单个模型校准数据的闭包函数
  calc_single_model <- function(f, model_name) {
    # 1. 在训练集拟合 Cox
    fit <- coxph(f, data = train_data)
    
    # 2. 提取基线风险 H0(t) 用于计算绝对概率
    bh <- basehaz(fit, centered = TRUE)
    bh_t <- bh[bh$time <= t_eval, ]
    h0 <- ifelse(nrow(bh_t) == 0, 0, bh_t$hazard[nrow(bh_t)])
    
    # 3. 计算测试集个体的 5 年绝对发病概率
    lp <- predict(fit, newdata = test_data, type = "lp")
    test_data$pred_prob <- 1 - exp(-h0 * exp(lp))
    
    # 4. 计算 Brier Score 和 Slope
    status_t <- ifelse(test_data$duration_updated <= t_eval & test_data$event_af_updated == 1, 1, 0)
    brier <- mean((test_data$pred_prob - status_t)^2, na.rm = TRUE)
    
    # 5. 分 10 组提取真实点
    get_obs_pred <- function(data) {
      data$group <- ntile(data$pred_prob, 10)
      data %>%
        group_by(group) %>%
        summarise(
          pred = mean(pred_prob, na.rm = TRUE),
          obs = {
            sf <- summary(survfit(Surv(duration_updated, event_af_updated) ~ 1), times = t_eval, extend = TRUE)
            if(length(sf$surv) > 0) 1 - sf$surv else NA
          },
          .groups = 'drop'
        ) %>% filter(!is.na(obs))
    }
    
    raw_pts <- get_obs_pred(test_data)
    
    # 计算 Slope
    cal_lm <- lm(qlogis(pmax(pmin(raw_pts$obs, 0.999), 0.001)) ~ qlogis(pmax(pmin(raw_pts$pred, 0.999), 0.001)))
    slope <- coef(cal_lm)[2]
    
    # 6. Bootstrap 获取阴影
    boot_list <- list()
    for(i in 1:n_boot) {
      set.seed(2026 + i)
      boot_idx <- sample(1:nrow(test_data), replace = TRUE)
      boot_pts <- get_obs_pred(test_data[boot_idx, ])
      boot_pts$boot_id <- i
      boot_list[[i]] <- boot_pts
    }
    
    df_boot <- bind_rows(boot_list)
    raw_pts$Model <- model_name
    df_boot$Model <- model_name
    
    return(list(raw = raw_pts, boot = df_boot, brier = brier, slope = slope, model = model_name))
  }
  
  # 执行四个模型
  res_all <- lapply(names(f_list), function(m_name) {
    calc_single_model(f_list[[m_name]], m_name)
  })
  
  return(res_all)
}

# ==============================================================================
# 2. 顶级期刊专属绘图函数：四色嵌套 + 动态坐标系
# ==============================================================================
draw_4_calibration <- function(res_list, title) {
  
  # 合并数据
  df_raw <- bind_rows(lapply(res_list, function(x) x$raw))
  df_boot <- bind_rows(lapply(res_list, function(x) x$boot))
  
  # 固定因子顺序，对齐 ROC 颜色隐喻
  model_levels <- c("4. Full Model (Both)", "3. + Metabolomics (MetRS)", 
                    "2. + Genetics (PRS)", "1. Clinical Alone")
  df_raw$Model <- factor(df_raw$Model, levels = model_levels)
  df_boot$Model <- factor(df_boot$Model, levels = model_levels)
  
  # 颜色映射 (与前面所有图保持绝对一致)
  color_map <- c("4. Full Model (Both)" = "#D62728", 
                 "3. + Metabolomics (MetRS)" = "#1F77B4", 
                 "2. + Genetics (PRS)" = "#FF7F0E", 
                 "1. Clinical Alone" = "#7F7F7F")
  
  # 动态确定坐标轴上限 (防止高风险点被截断)
  limit_val <- max(df_raw$pred, df_raw$obs, na.rm = TRUE) * 1.15
  limit_val <- pmin(limit_val, 0.5) # 通常房颤 5 年概率不会超过 50%，这里锁个上限
  
  # 生成指标文本 (Brier & Slope)
  anno_df <- data.frame(
    Model = sapply(res_list, function(x) x$model),
    Brier = sapply(res_list, function(x) x$brier),
    Slope = sapply(res_list, function(x) x$slope)
  )
  anno_df$Model <- factor(anno_df$Model, levels = model_levels)
  anno_df <- anno_df[order(anno_df$Model), ]
  
  # 🚨 核心修复：把 Brier Score 加回来！(保留 4 位小数，因为 Brier 通常很小)
  labels_text <- sprintf("%s | Brier: %.4f | Slope: %.2f", anno_df$Model, anno_df$Brier, anno_df$Slope)
  
  p <- ggplot() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.6) +
    stat_smooth(data = df_boot, aes(x = pred, y = obs, fill = Model, color = Model), 
                alpha = 0.1, method = "loess", linewidth = 0.1, formula = y ~ x, se = FALSE) +
    stat_smooth(data = df_raw, aes(x = pred, y = obs, color = Model), 
                method = "loess", se = FALSE, linewidth = 1.2, formula = y ~ x) +
    geom_point(data = df_raw, aes(x = pred, y = obs, color = Model), size = 2.5, alpha = 0.9) +
    coord_cartesian(xlim = c(0, limit_val), ylim = c(0, limit_val)) + 
    scale_color_manual(values = color_map) +
    scale_fill_manual(values = color_map) +
    labs(title = title, x = "Predicted 5-year AF Risk", y = "Observed 5-year AF Risk") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.minor = element_blank(), 
      legend.position = "none", 
      plot.title = element_text(face = "bold", size = 16, margin = margin(b=10))
    ) +
    
    # 🚨 精准右下角排版：展示 Brier 和 Slope，字号微调为 3.5 防止太长重叠
    annotate("text", x = limit_val * 0.98, y = limit_val * 0.28, label = labels_text[1], color = "#D62728", fontface = "bold", hjust = 1, size = 3.5) +
    annotate("text", x = limit_val * 0.98, y = limit_val * 0.20, label = labels_text[2], color = "#1F77B4", fontface = "bold", hjust = 1, size = 3.5) +
    annotate("text", x = limit_val * 0.98, y = limit_val * 0.12, label = labels_text[3], color = "#FF7F0E", fontface = "bold", hjust = 1, size = 3.5) +
    annotate("text", x = limit_val * 0.98, y = limit_val * 0.04, label = labels_text[4], color = "#7F7F7F", fontface = "bold", hjust = 1, size = 3.5)
  
  return(p)}
# ==============================================================================
# 3. 自动化流水线执行
# ==============================================================================
cat("\n--- 开始计算校准曲线 (包含 Bootstrap 阴影计算) ---\n")

# 计算数据 (n_boot = 100 足够画阴影了，速度快且平滑)
calib_charge <- get_4_calib_data(train_eval_final, test_eval_final, "CHARGE_AF_LP", n_boot = 5)
calib_aric   <- get_4_calib_data(train_eval_final, test_eval_final, "ARIC_LP_centered", n_boot = 5)
calib_c2hest <- get_4_calib_data(train_eval_final, test_eval_final, "C2HEST_Points", n_boot = 5)
# 重新画图即可，不用重新跑 Bootstrap，因为数据已经算好存在 res_list 里了
p_cal_1 <- draw_4_calibration(calib_charge, "A. Baseline: CHARGE-AF")
p_cal_2 <- draw_4_calibration(calib_aric,   "B. Baseline: ARIC")
p_cal_3 <- draw_4_calibration(calib_c2hest, "C. Baseline: C2HEST")

final_calib_fig <- p_cal_1 | p_cal_2 | p_cal_3
ggsave("Figure_Calibration_Nested_4_Models.pdf", final_calib_fig, width = 18, height = 6.5, dpi = 400)







# 安装并加载必须包
if (!requireNamespace("survIDINRI", quietly = TRUE)) install.packages("survIDINRI")
library(survIDINRI)
library(dplyr)
library(ggplot2)
library(patchwork)
library(flextable)
library(officer)

cat("\n==== 🎯 启动 NRI/IDI 终极计算 (10,000人极速测试版) ====\n")

# ==============================================================================
# 🚨 上服务器前必改区 🚨
# ==============================================================================
set.seed(2026)

# 测试模式：抽 10,000 人；正式跑：直接改成 test_sample <- test_eval_final
test_sample <- test_eval_final %>% sample_n(min(10000, n())) 

# 测试模式：5 次；正式跑：改成 500
n_boot <- 5 

cat(sprintf("✅ 当前状态：测试集抽取 %d 人，Bootstrap 设置为 %d 次\n", nrow(test_sample), n_boot))

# ==============================================================================
# 1. 自动适配时间单位
# ==============================================================================
max_t <- max(test_sample$duration_updated, na.rm = TRUE)
if (max_t > 100) {
  t_eval_actual <- 5 * 365.25
  cat(sprintf("📏 时间单位检测为 [天]，评估点自动设为: %.1f 天\n", t_eval_actual))
} else {
  t_eval_actual <- 5
  cat(sprintf("📏 时间单位检测为 [年]，评估点自动设为: %d 年\n", t_eval_actual))
}

# ==============================================================================
# 2. 自动化 survIDINRI 计算包装函数
# ==============================================================================
run_idi_nri <- function(data, cov_base, cov_new) {
  indata <- as.matrix(data[, c("duration_updated", "event_af_updated")])
  
  # 动态生成 covariate 矩阵 (去除截距项)
  f0 <- as.formula(paste("~", paste(cov_base, collapse = " + ")))
  f1 <- as.formula(paste("~", paste(c(cov_base, cov_new), collapse = " + ")))
  
  covs0 <- model.matrix(f0, data = data)[, -1, drop = FALSE]
  covs1 <- model.matrix(f1, data = data)[, -1, drop = FALSE]
  
  # 运行 IDI.INF
  res <- IDI.INF(indata, covs0, covs1, t0 = t_eval_actual, npert = n_boot)
  return(res)
}

cat("\n⏳ 正在计算 CHARGE-AF 增量...\n")
res_charge_mrs     <- run_idi_nri(test_sample, c("CHARGE_AF_LP"), "MetRS_Z")
res_charge_prs_mrs <- run_idi_nri(test_sample, c("CHARGE_AF_LP", "PRS_AF"), "MetRS_Z")

cat("⏳ 正在计算 ARIC 增量...\n")
res_aric_mrs       <- run_idi_nri(test_sample, c("ARIC_LP_centered"), "MetRS_Z")
res_aric_prs_mrs   <- run_idi_nri(test_sample, c("ARIC_LP_centered", "PRS_AF"), "MetRS_Z")

cat("⏳ 正在计算 C2HEST 增量...\n")
res_c2hest_mrs     <- run_idi_nri(test_sample, c("C2HEST_Points"), "MetRS_Z")
res_c2hest_prs_mrs <- run_idi_nri(test_sample, c("C2HEST_Points", "PRS_AF"), "MetRS_Z")


# 安装并加载包
if (!requireNamespace("nricens", quietly = TRUE)) install.packages("nricens")
library(nricens)
library(dplyr)
library(survival)
library(ggplot2)
library(scales)
library(patchwork)
library(flextable)
library(officer)

cat("\n==== 🎯 启动权威 nricens 引擎： NRI & IDI (10,000人极速测试版) ====\n")
# 1. 使用完整的十万人测试集
test_sample <- test_eval_final 

# 2. 释放顶刊级别的重抽样火力
n_iter_test <- 500

# ==============================================================================
# 2. 核心计算引擎 (整合你之前的手算 IDI 与 nricens)
# ==============================================================================
run_nricens_idi_combo <- function(data, base_vars, new_vars, comp_label) {
  cat(sprintf("⏳ 正在计算: %s ...\n", comp_label))
  
  # 1. 构建公式
  f_base <- as.formula(paste("Surv(duration_updated, event_af_updated) ~", paste(base_vars, collapse=" + ")))
  f_comb <- as.formula(paste("Surv(duration_updated, event_af_updated) ~", paste(c(base_vars, new_vars), collapse=" + ")))
  
  # 2. 拟合 Cox 模型
  fit_s <- coxph(f_base, data = data, x = TRUE)
  fit_n <- coxph(f_comb, data = data, x = TRUE)
  
  # 3. 提取生存时间与状态
  s_time <- data$duration_updated
  s_stat <- data$event_af_updated
  
  # 4. 计算 5 年绝对风险概率
  s0_5y_local <- summary(survfit(Surv(s_time, s_stat) ~ 1), times = t_eval_actual)$surv
  if(length(s0_5y_local) == 0) s0_5y_local <- min(summary(survfit(Surv(s_time, s_stat) ~ 1))$surv)
  
  ps <- as.numeric(1 - (s0_5y_local ^ exp(predict(fit_s, type="lp") - mean(predict(fit_s, type="lp"), na.rm=T))))
  pn <- as.numeric(1 - (s0_5y_local ^ exp(predict(fit_n, type="lp") - mean(predict(fit_n, type="lp"), na.rm=T))))
  
  # 5. 核心：调用 nricens 算 NRI
  set.seed(2026)
  res_nri <- tryCatch({
    nricens::nricens(p.std = ps, p.new = pn, time = as.numeric(s_time), 
                     event = as.numeric(s_stat), t0 = t_eval_actual, 
                     updown = "diff", cut = 0, niter = n_iter_test, msg = FALSE)
  }, error = function(e) { message("nricens error:", e$message); return(NULL) })
  
  n_est <- if(!is.null(res_nri)) as.numeric(res_nri$nri["NRI", 1]) else NA_real_
  n_l   <- if(!is.null(res_nri)) as.numeric(res_nri$nri["NRI", 2]) else NA_real_
  n_u   <- if(!is.null(res_nri)) as.numeric(res_nri$nri["NRI", 3]) else NA_real_
  
  n_se <- (n_u - n_l) / (2 * 1.96)
  n_p  <- if(!is.na(n_se) && n_se > 0) 2 * (1 - pnorm(abs(n_est / n_se))) else NA_real_
  
  # 6. 核心 B：终极经验法手算 IDI (精确剔除早期删失干扰，永不报错)
  diff_p <- pn - ps
  
  # 真正的 Event：在 5 年（t_eval_actual）内明确发病的人
  idx_event <- which(s_time <= t_eval_actual & s_stat == 1)
  
  # 真正的 Non-Event：随访时间明确超过 5 年，且没有发病的人
  idx_nonevent <- which(s_time > t_eval_actual)
  
  # 提取概率差值
  diff_event <- diff_p[idx_event]
  diff_nonevent <- diff_p[idx_nonevent]
  
  # 计算 IDI 及其标准误 (Pencina 公式)
  idi_est <- mean(diff_event, na.rm=TRUE) - mean(diff_nonevent, na.rm=TRUE)
  idi_se  <- sqrt(var(diff_event, na.rm=TRUE)/length(diff_event) + var(diff_nonevent, na.rm=TRUE)/length(diff_nonevent))
  
  idi_l <- idi_est - 1.96 * idi_se
  idi_u <- idi_est + 1.96 * idi_se
  idi_p <- 2 * (1 - pnorm(abs(idi_est / idi_se)))
  
  # 7. 汇总指标
  summary_stats <- data.frame(
    Comparison = comp_label,
    NRI = n_est, NRI_LCI = n_l, NRI_UCI = n_u, NRI_P = n_p,
    IDI = idi_est, IDI_LCI = idi_l, IDI_UCI = idi_u, IDI_P = idi_p,
    stringsAsFactors = FALSE
  )
  
  # 8. 散点图底表
  keep_idx <- !(s_time < t_eval_actual & s_stat == 0) # 剔除删失在评估点之前的噪音
  plot_data <- data.frame(
    Comparison = comp_label,
    Prob_Base = ps[keep_idx], 
    Prob_Comb = pn[keep_idx],
    Status = factor(ifelse(s_time[keep_idx] <= t_eval_actual & s_stat[keep_idx] == 1, "Event (AF)", "Non-Event (No AF)"),
                    levels = c("Event (AF)", "Non-Event (No AF)"))
  )
  
  return(list(summary = summary_stats, plot = plot_data))
}

# ==============================================================================
# 3. 批量执行 6 组对比
# ==============================================================================
res_1 <- run_nricens_idi_combo(test_sample, c("CHARGE_AF_LP"), "MetRS_Z", "1. CHARGE-AF vs +MetRS")
res_2 <- run_nricens_idi_combo(test_sample, c("CHARGE_AF_LP", "PRS_AF"), "MetRS_Z", "2. CHARGE-AF+PRS vs +PRS+MetRS")

res_3 <- run_nricens_idi_combo(test_sample, c("ARIC_LP_centered"), "MetRS_Z", "3. ARIC vs +MetRS")
res_4 <- run_nricens_idi_combo(test_sample, c("ARIC_LP_centered", "PRS_AF"), "MetRS_Z", "4. ARIC+PRS vs +PRS+MetRS")

res_5 <- run_nricens_idi_combo(test_sample, c("C2HEST_Points"), "MetRS_Z", "5. C2HEST vs +MetRS")
res_6 <- run_nricens_idi_combo(test_sample, c("C2HEST_Points", "PRS_AF"), "MetRS_Z", "6. C2HEST+PRS vs +PRS+MetRS")

# 合并结果
final_stats <- bind_rows(res_1$summary, res_2$summary, res_3$summary, res_4$summary, res_5$summary, res_6$summary)
final_plot  <- bind_rows(res_1$plot, res_2$plot, res_3$plot, res_4$plot, res_5$plot, res_6$plot)
# ==============================================================================
# 4. 生成顶级 Word 三线表 (精调 IDI 小数位数版)
# ==============================================================================
cat("\n📝 正在生成 Word 三线表...\n")

table_df <- final_stats %>%
  mutate(
    # cNRI 保留 3 位小数足够了
    `Continuous NRI (95% CI)` = sprintf("%.3f (%.3f-%.3f)", NRI, NRI_LCI, NRI_UCI),
    `NRI P-value` = ifelse(NRI_P < 0.001, "< 0.001", sprintf("%.3f", NRI_P)),
    
    # 🚨 核心修复：IDI 强行保留 4 位小数！
    `IDI (95% CI)` = sprintf("%.4f (%.4f-%.4f)", IDI, IDI_LCI, IDI_UCI), 
    `IDI P-value` = ifelse(IDI_P < 0.001, "< 0.001", sprintf("%.3f", IDI_P))
  ) %>%
  select(Comparison, `Continuous NRI (95% CI)`, `NRI P-value`, `IDI (95% CI)`, `IDI P-value`)

# ... 下面的 flextable 导出代码保持不变 ...

ft <- flextable(table_df) %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "all") %>%
  border_remove() %>%
  hline_top(border = fp_border(width = 2), part = "header") %>%
  hline_bottom(border = fp_border(width = 1), part = "header") %>%
  hline_bottom(border = fp_border(width = 2), part = "all") %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "all") %>%
  bold(part = "header") %>%
  autofit()

doc <- read_docx() %>%
  body_add_par("Table: Reclassification improvement (NRI & IDI) of the MetRS score.", style = "heading 1") %>%
  body_add_flextable(ft)
print(doc, target = "Table_NRI_IDI_nricens_Version.docx")

library(dplyr)
library(stringr) # 确保加载了正则包

# ==============================================================================
# 5. 绘制超清分面重分类散点图 (究极矩阵排版：一个评分一栏)
# ==============================================================================
cat("🎨 正在生成分面散点图 (二维矩阵排版)...\n")

# 1. 写一个小函数：把杂乱的 Comparison 名字拆解成“行”和“列”两个维度
process_grid_labels <- function(df) {
  df %>%
    mutate(
      # 列维度：三大临床评分 (提取纯名字)
      Baseline_Score = case_when(
        grepl("CHARGE-AF", Comparison) ~ "CHARGE-AF",
        grepl("ARIC", Comparison) ~ "ARIC",
        grepl("C2HEST", Comparison) ~ "C2HEST"
      ),
      Baseline_Score = factor(Baseline_Score, levels = c("CHARGE-AF", "ARIC", "C2HEST")),
      
      # 行维度：是否自带 PRS 基因基线
      Model_Tier = case_when(
        grepl("\\+PRS", Comparison) ~ "Base: Clinical + PRS",
        TRUE ~ "Base: Clinical Alone"
      ),
      Model_Tier = factor(Model_Tier, levels = c("Base: Clinical Alone", "Base: Clinical + PRS"))
    )
}

# 2. 对画图底表和文本标签应用双维度拆分
final_plot_grid <- process_grid_labels(final_plot) %>%
  # 🚨 核心修复1：图层顺序 (Z-order) 调整
  # 判断逻辑：是 Event 结果为 TRUE (排在数据框尾部)，不是则为 FALSE (排在头部)
  # ggplot 会最后画尾部数据，确保红点(事件)绝对浮在蓝点上方！
  arrange(Status == "Event (AF)")

ann_text_grid <- process_grid_labels(final_stats) %>%
  mutate(
    # 顺便把文字稍微加个行距，排版更好看
    stats_label = sprintf("cNRI: %.3f (P %s)\nIDI: %.3f (P %s)", 
                          NRI, ifelse(NRI_P < 0.001, "<0.001", sprintf("= %.3f", NRI_P)),
                          IDI, ifelse(IDI_P < 0.001, "<0.001", sprintf("= %.3f", IDI_P)))
  )

# 设个统一的最大坐标轴
max_p <- 0.20 

# ==============================================================================
# 3. 开始终极作图 (微调版：去顶部标签、字号翻倍、蓝点更透明)
# ==============================================================================
p_scatter_grid <- ggplot(final_plot_grid, aes(x = Prob_Base, y = Prob_Comb)) +
  # 🔵 修改点 1：蓝点透明度拉高（更透明）。将 0.3 改为了 0.15。
  # (注：如果你的“拉高”是指想让它更显眼/不那么透明，请把 0.15 改成 0.5)
  geom_point(aes(color = Status), size = 0.8, alpha = ifelse(final_plot_grid$Status == "Event (AF)", 0.8, 0.8)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.8) +
  # 保持对比配色
  scale_color_manual(values = c("Event (AF)" = "#D62728", "Non-Event (No AF)" = "#A0C4DF")) +
  
  # 🔵 修改点 2：字体放大两倍 (size 从 3.6 翻倍为 7.2)
  geom_text(data = ann_text_grid, aes(x = max_p*0.05, y = max_p*0.95, label = stats_label), 
            hjust = 0, vjust = 1, size = 7.2, fontface = "bold", color = "grey20", lineheight = 1.2) +
  
  # 十字矩阵排版！(行: 基础模型层级 ~ 列: 三大临床评分)
  facet_grid(Model_Tier ~ Baseline_Score) +
  
  scale_x_continuous(labels = scales::percent, limits = c(0, max_p)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, max_p)) +
  coord_fixed() + 
  labs(x = "Predicted Risk (Baseline Model)", y = "Predicted Risk (Baseline + MetRS)", color = "Actual 5-Year Status") +
  theme_bw(base_size = 21) +
  theme(
    legend.position = "bottom",
    
    # 🔵 修改点 3：单独去掉顶部 (x轴方向分面) 的标签框和文字
    strip.background.x = element_blank(),
    strip.text.x = element_blank(),
    
    # 保留右侧 (y轴方向分面) 标签的质感
    strip.background.y = element_rect(fill = "#f0f2f5", color = "black", linewidth = 0.8),
    strip.text.y = element_text(face = "bold", size = 18),
    
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", linewidth = 0.8)
  )

# 统一排版画幅
ggsave("Figure_Reclass_Scatter_nricens_Grid.pdf", p_scatter_grid, width = 18, height = 12, dpi = 400)

cat("✅ 散点图微调完成！顶部标签已移除，统计文字已放大，蓝点已调整透明度。\n")
cat("✅ 排版修正完成！红点已置于最上层，矩阵散点图已输出至 Figure_Reclass_Scatter_nricens_Grid.pdf\n")






library(dplyr)
library(tidyr)
library(ggplot2)

cat("\n==== 🎨 开始生成顶刊级 NRI/IDI 悬浮森林图 ====\n")

# 设定图表默认评估时间（5年）
eval_time <- 5 

# ==============================================================================
# 🔧 辅助函数：将 final_stats 大表转换为画图所需的 summary_df 格式
# ==============================================================================
prepare_forest_data <- function(stats_subset, model_names) {
  # 提取 cNRI 数据
  df_nri <- stats_subset %>%
    select(Comparison, Est = NRI, Lower = NRI_LCI, Upper = NRI_UCI) %>%
    mutate(Metric = "cNRI", Model = model_names)
  
  # 提取 IDI 数据
  df_idi <- stats_subset %>%
    select(Comparison, Est = IDI, Lower = IDI_LCI, Upper = IDI_UCI) %>%
    mutate(Metric = "IDI", Model = model_names)
  
  # 合并
  bind_rows(df_nri, df_idi) %>%
    mutate(Metric = factor(Metric, levels = c("cNRI", "IDI"))) # 锁定分面顺序
}

# ==============================================================================
# 🌲 第一张图：临床基线 vs (+ MetRS)
# ==============================================================================
cat("⏳ 正在绘制第一张图 (基线: 临床评分)...\n")

# 提取第 1, 3, 5 行（不带 PRS 的对比）
stats_base <- final_stats[c(1, 3, 5), ]
summary_df_1 <- prepare_forest_data(stats_base, c("CHARGE-AF", "ARIC", "C2HEST"))

# 1. 规范化排序与标签
summary_df_1$Model_Label <- factor(
  summary_df_1$Model, 
  levels = rev(c("CHARGE-AF", "ARIC", "C2HEST")),
  labels = rev(c("CHARGE-AF\n(+ MetRS)", "ARIC\n(+ MetRS)", "C2HEST\n(+ MetRS)"))
)

# 2. 生成超级优雅的文本标签
summary_df_1$Text_CI <- sprintf("%.4f\n(%.4f - %.4f)", summary_df_1$Est, summary_df_1$Lower, summary_df_1$Upper)

# 3. 顶刊审美绘图
p_nri_1 <- ggplot(summary_df_1, aes(x = Est, y = Model_Label, color = Metric)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#F2F7F4", alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#8C8C8C", linewidth = 1) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0, linewidth = 1.2) +
  geom_point(size = 4.5, shape = 21, fill = "white", stroke = 1.8) +
  geom_text(aes(label = Text_CI), vjust = -0.6, size = 3.5, fontface = "bold", color = "black", lineheight = 0.9) +
  facet_wrap(~Metric, scales = "free_x", strip.position = "top") +
  scale_color_manual(values = c("IDI" = "#B24745", "cNRI" = "#005083")) +
  labs(
    title = "Incremental Predictive Value of Metabolic Risk Score",
    subtitle = paste0("Reference: Clinical Baseline  |  Added: MetRS (", eval_time, "-Year AF Risk)"),
    x = "Net Estimate Value (with 95% Confidence Interval)", 
    y = "Baseline Predictive Model",
    caption = "Note: Point estimates represent the net increment in reclassification (cNRI) and discrimination (IDI)\nbrought by the 8-metabolite risk score (MetRS). Value format: Estimate (95% CI)."
  ) +
  theme_minimal(base_size = 15, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 18, margin = margin(b = 5)),
    plot.subtitle = element_text(color = "#b24745", face = "bold", size = 13, margin = margin(b = 15)),
    plot.caption = element_text(color = "grey40", size = 11, hjust = 1, face = "italic", margin = margin(t = 15)),
    axis.text.y = element_text(face = "bold", color = "black", size = 13, lineheight = 1.2),
    axis.text.x = element_text(color = "black", size = 12),
    axis.title.x = element_text(face = "bold", margin = margin(t = 12)),
    axis.title.y = element_text(face = "bold", margin = margin(r = 12)),
    strip.background = element_rect(fill = "#E8ECEF", color = NA),
    strip.text = element_text(face = "bold", size = 14, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "none"
  ) +
  scale_y_discrete(expand = expansion(mult = c(0.15, 0.25)))

ggsave("nri_idi_clinical_base.pdf", p_nri_1, width = 11, height = 6.5, device = cairo_pdf)
cat("✅ 第一张森林图已生成：nri_idi_clinical_base.pdf\n")


# ==============================================================================
# 🌲 第二张图：临床基线 + PRS vs (+ PRS + MetRS)
# ==============================================================================
cat("⏳ 正在绘制第二张图 (基线: 临床评分 + PRS)...\n")

# 提取第 2, 4, 6 行（带 PRS 的对比）
stats_prs <- final_stats[c(2, 4, 6), ]
summary_df_2 <- prepare_forest_data(stats_prs, c("CHARGE-AF", "ARIC", "C2HEST"))

# 1. 规范化排序与标签
summary_df_2$Model_Label <- factor(
  summary_df_2$Model, 
  levels = rev(c("CHARGE-AF", "ARIC", "C2HEST")),
  labels = rev(c("CHARGE-AF + PRS\n(+ MetRS)", "ARIC + PRS\n(+ MetRS)", "C2HEST + PRS\n(+ MetRS)"))
)

# 2. 文本标签
summary_df_2$Text_CI <- sprintf("%.4f\n(%.4f - %.4f)", summary_df_2$Est, summary_df_2$Lower, summary_df_2$Upper)

# 3. 绘图
p_nri_2 <- ggplot(summary_df_2, aes(x = Est, y = Model_Label, color = Metric)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = Inf, fill = "#F2F7F4", alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#8C8C8C", linewidth = 1) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0, linewidth = 1.2) +
  geom_point(size = 4.5, shape = 21, fill = "white", stroke = 1.8) +
  geom_text(aes(label = Text_CI), vjust = -0.6, size = 3.5, fontface = "bold", color = "black", lineheight = 0.9) +
  facet_wrap(~Metric, scales = "free_x", strip.position = "top") +
  scale_color_manual(values = c("IDI" = "#B24745", "cNRI" = "#005083")) +
  labs(
    title = "Incremental Predictive Value Beyond Genetics (PRS)",
    subtitle = paste0("Reference: Clinical Baseline + PRS  |  Added: MetRS (", eval_time, "-Year AF Risk)"),
    x = "Net Estimate Value (with 95% Confidence Interval)", 
    y = "Baseline Predictive Model",
    caption = "Note: Point estimates represent the net increment in reclassification (cNRI) and discrimination (IDI)\nbrought by the 8-metabolite risk score (MetRS). Value format: Estimate (95% CI)."
  ) +
  theme_minimal(base_size = 15, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 18, margin = margin(b = 5)),
    plot.subtitle = element_text(color = "#b24745", face = "bold", size = 13, margin = margin(b = 15)),
    plot.caption = element_text(color = "grey40", size = 11, hjust = 1, face = "italic", margin = margin(t = 15)),
    axis.text.y = element_text(face = "bold", color = "black", size = 13, lineheight = 1.2),
    axis.text.x = element_text(color = "black", size = 12),
    axis.title.x = element_text(face = "bold", margin = margin(t = 12)),
    axis.title.y = element_text(face = "bold", margin = margin(r = 12)),
    strip.background = element_rect(fill = "#E8ECEF", color = NA),
    strip.text = element_text(face = "bold", size = 14, color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "none"
  ) +
  scale_y_discrete(expand = expansion(mult = c(0.15, 0.25)))

ggsave("nri_idi_clinical_prs_base.pdf", p_nri_2, width = 11, height = 6.5, device = cairo_pdf)
cat("✅ 第二张森林图已生成：nri_idi_clinical_prs_base.pdf\n")
cat("🎉 两张顶级森林图全部绘制完毕，请查看文件！\n")

library(mice)
library(survival)
library(broom)
library(dplyr)
library(stringr)
library(flextable)
library(officer)
library(mice)
library(dplyr)
library(tidyr)


library(mice)
library(dplyr)
library(tidyr)
library(stringr)

imp_ready_for_cox <- mice_rediag
library(mice)
library(survival)
library(broom)
library(dplyr)
library(stringr)
library(flextable)
library(officer)
library(tidyr)

# ==============================================================================
# 0 & 1. 直接提取已经算好的连续评分 (极大简化逻辑)
# ==============================================================================
cat("⏳ 正在提取连续的 MetRS 评分并合并至插补数据...\n")

# 之前我们已经把全量队列的完美打分算好存进 df_benchmark_scores 了
# 我们直接提取它作为 MRS_Score，绝对不会报错！
score_df <- df_benchmark_scores %>% 
  select(`Participant ID` = eid, MRS_Score = MetRS_Z)

# ==============================================================================
# 2. 将分数无缝贴回全量地基 (long_dat_orig)
# ==============================================================================
long_dat_orig <- complete(imp_ready_for_cox, "long", include = TRUE)
long_dat_merged <- long_dat_orig %>% left_join(score_df, by = "Participant ID")

# 自动抓取 PRS
prs_col_name <- grep("PRS.*fibrillation", colnames(long_dat_merged), ignore.case = TRUE, value = TRUE)[1]
if(!is.na(prs_col_name)) {
  long_dat_merged$PRS_AF <- long_dat_merged[[prs_col_name]]
} else {
  fallback_prs <- grep("PRS", colnames(long_dat_merged), ignore.case = TRUE, value = TRUE)[1]
  long_dat_merged$PRS_AF <- long_dat_merged[[fallback_prs]]
}
# ==============================================================================
# 2.5 核心修改：在训练集提取阈值 (Cutoffs)，然后套用到【全量队列】
# ==============================================================================

# 1. 依然在训练集里提取绝对阈值 (严谨性拉满，防止数据泄露)
train_ids_list <- df_train_final$eid  
train_dat <- long_dat_merged %>% filter(`Participant ID` %in% train_ids_list)

prs_cutoffs <- quantile(train_dat$PRS_AF, probs = c(1/3, 2/3), na.rm = TRUE)
mrs_cutoffs <- quantile(train_dat$MRS_Score, probs = c(1/3, 2/3), na.rm = TRUE)

cat("📐 严谨标定：提取到训练集 PRS 阈值:", round(prs_cutoffs, 3), "\n")
cat("📐 严谨标定：提取到训练集 MRS 阈值:", round(mrs_cutoffs, 3), "\n")

# 2. 🚨 关键改变：不 filter 测试集了... (这句注释可以删了)
# ==============================================================================
# 3. 缩尾与联合分组 (仅在独立测试集跑交互！)
# ==============================================================================
winsor_p99 <- function(x) {
  limit <- quantile(x, 0.99, na.rm = TRUE)
  x[x > limit] <- limit
  return(x)
}

# 🚨 核心修改：通过剔除 train_ids_list，强制把全量数据变成独立的 30% 测试集
long_dat_clean <- long_dat_merged %>%
  #filter(!(`Participant ID` %in% train_ids_list)) %>%  # 👈👈👈 仅仅加了这一行！
  group_by(.imp) %>%
  mutate(
    fib4 = winsor_p99(fib4),
    # ... 后续代码完全不变
    apri = winsor_p99(apri),
    nfs = winsor_p99(nfs),
    ast_alt_ratio = winsor_p99(ast_alt_ratio),
    
    # 🚨 修复报错核心：严格使用训练集绝对阈值，先生成模型需要的数值型变量 (1, 2, 3)
    PRS_3 = case_when(
      PRS_AF <= prs_cutoffs[1] ~ 1,
      PRS_AF <= prs_cutoffs[2] ~ 2,
      TRUE                     ~ 3
    ),
    MRS_3 = case_when(
      MRS_Score <= mrs_cutoffs[1] ~ 1,
      MRS_Score <= mrs_cutoffs[2] ~ 2,
      TRUE                        ~ 3
    ),
    
    # 基于生成的数值型变量，映射对应的字符标签
    PRS_Label = case_when(
      PRS_3 == 1 ~ "Low PRS",
      PRS_3 == 2 ~ "Mid PRS",
      PRS_3 == 3 ~ "High PRS"
    ),
    MRS_Label = case_when(
      MRS_3 == 1 ~ "Low MRS",
      MRS_3 == 2 ~ "Mid MRS",
      MRS_3 == 3 ~ "High MRS"
    ),
    
    Joint_Group = if_else(is.na(PRS_Label) | is.na(MRS_Label), 
                          NA_character_, 
                          paste0(PRS_Label, " + ", MRS_Label))
  ) %>%
  ungroup()

# ==============================================================================
# 4. 提纯与重铸 mids 
# ==============================================================================
valid_ids <- long_dat_clean %>% 
  filter(.imp == 1 & !is.na(MRS_Score) & !is.na(PRS_AF)) %>% 
  pull(.id)

long_dat_balanced <- long_dat_clean %>%
  filter(.id %in% valid_ids) %>%
  arrange(.imp, .id) 

long_dat_balanced$Joint_Group <- factor(long_dat_balanced$Joint_Group)
if("Low PRS + Low MRS" %in% levels(long_dat_balanced$Joint_Group)) {
  long_dat_balanced$Joint_Group <- relevel(long_dat_balanced$Joint_Group, ref = "Low PRS + Low MRS")
}

imp_final_joint <- as.mids(long_dat_balanced)
cat("🎉 全量队列的 imp_final_joint 生成成功，MRS_3 变量已补齐！\n")
cat("🎉 全量队列的 imp_final_joint 生成成功，样本量极大，准备召唤显著的 RERI！\n")
# 保存以防万一
saveRDS(imp_final_joint, "imp_final_joint_Ridge.rds")
#imp_final_joint<-imp_final_joint_Ridge
# =======================================================
# 7. 直接运行 Cox 模型 (修复带空格列名的语法报错)
# =======================================================
cat("⏳ 正在运行联合风险的多重插补 Cox 模型...\n")

cov_m1 <- c("age_val", "sex_f", "ethnicity_f")
cov_m2 <- c(cov_m1, "bmi_val", "smoke_f", "alc_freq_f", "tdi_val")
cov_m3 <- c(cov_m2, "htn_final", "t2dm_final", "lip_f", "ckd_final", "cvd_f")

# 🚨 这里的反引号是保命符，绝对不能删！
pc_vars <- paste0("`Genetic principal components | Array ", 1:10, "`")
covs_joint <- c(cov_m3, pc_vars) 
form_joint <- as.formula(paste("Surv(duration_updated, event_af_updated) ~ Joint_Group +", paste(covs_joint, collapse = " + ")))
cat("⏳ 正在绕过环境 Bug，直接基于各插补层运行 Cox 模型...\n")

# =======================================================
# 1. 准备模型参数
# =======================================================
m_count <- imp_final_joint$m  # 获取插补次数
model_list <- list()          # 创建空列表存储模型结果

# =======================================================
# 2. 核心循环：手动提取每一层数据并运行回归
# =======================================================
for (i in 1:m_count) {
  # 🚨 关键：直接提取第 i 次插补的完整数据框
  dat_i <- complete(imp_final_joint, action = i)
  
  # 直接在这个数据框上运行 Cox
  # 这里的 data = dat_i 保证了回归器一定能找到所有变量
  model_list[[i]] <- coxph(form_joint, data = dat_i)
  
  if(i %% 2 == 0) cat(sprintf("  已完成第 %d / %d 次插补计算...\n", i, m_count))
}

# =======================================================
# 3. 汇总合并 (Pooling)
# =======================================================
# 把 model_list 喂给 pool，效果和 with(...) %>% pool() 完全一致
pooled_joint <- pool(model_list)

library(ggplot2)
library(flextable)
library(officer)
library(dplyr)
library(stringr)
# 提取回归系数

res_joint_raw <- summary(pooled_joint, conf.int = TRUE, exponentiate = TRUE) %>%
  
  filter(str_starts(term, "Joint_Group")) %>%
  
  mutate(
    
    Group_Raw = str_replace(term, "Joint_Group", ""),
    
    HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
    
    P_val = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
    
  ) %>%
  
  select(Group_Raw, HR = estimate, Lower = conf.low, Upper = conf.high, HR_CI, P_val)


long_dat_balanced<-
  # 计算各组的 No. (Total/Cases) (使用第 1 次插补数据作为代表)
  
  dat_rep <- long_dat_balanced %>% filter(.imp == 1)



counts_df <- dat_rep %>%filter(!is.na(Joint_Group)) %>%
  
  group_by(Joint_Group) %>%
  
  summarise(
    
    N_Total = n(),
    
    N_Cases = sum(event_af_updated, na.rm = TRUE)
    
  ) %>%
  
  mutate(
    
    Counts = paste0(N_Total, "/", N_Cases),
    
    Group_Raw = as.character(Joint_Group)
    
  )



# 强制创建一个 Reference 组的行 (回归 summary 中默认没有 1.0 的行)

ref_row <- data.frame(
  
  Group_Raw = "Low PRS + Low MRS",
  
  HR = 1.00, Lower = 1.00, Upper = 1.00,
  
  HR_CI = "1.00 (Reference)",
  
  P_val = "-"
  
)



# 完美合并所有信息

final_joint_table <- bind_rows(ref_row, res_joint_raw) %>%filter(!str_detect(Group_Raw, "NA")) %>% # 🚨 新增：保险起见，把任何含有 NA 的无效组踢掉
  
  left_join(counts_df, by = "Group_Raw") %>%
  
  # 将联合标签拆分成两列，方便制表
  
  mutate(
    
    PRS_Level = str_split(Group_Raw, " \\+ ", simplify = TRUE)[, 1],
    
    MRS_Level = str_split(Group_Raw, " \\+ ", simplify = TRUE)[, 2]
    
  ) %>%
  
  # 按照 PRS 和 MRS 的逻辑顺序排序
  
  arrange(factor(PRS_Level, levels = c("Low PRS", "Mid PRS", "High PRS")),
          
          factor(MRS_Level, levels = c("Low MRS", "Mid MRS", "High MRS"))) %>%
  
  select(PRS_Level, MRS_Level, Counts, HR_CI, P_val, HR, Lower, Upper)



cat("✅ 结果提取完毕！\n")


# =======================================================
# 1. 准备绘图数据 (基于你已经生成的 final_joint_table)
# =======================================================
# 确保 HR, Lower, Upper 是数值型
plot_data <- final_joint_table %>%
  mutate(
    across(c(HR, Lower, Upper), as.numeric),
    # 设置森林图颜色：根据 MRS 等级区分，保持与阶梯图颜色一致
    plot_color = case_when(
      MRS_Level == "High MRS" ~ "#E64B35",
      MRS_Level == "Mid MRS"  ~ "#F39B7F",
      TRUE                    ~ "#999999"  # Low MRS / Reference
    )
  )

# =======================================================
# 2. 生成嵌入式森林图列表 (逐行生成小图)
# =======================================================
# 定义 X 轴范围 (根据你的 HR 最大值动态调整，通常 0.5 到 10)
x_min <- 0.5
x_max <- ceiling(max(plot_data$Upper, na.rm = TRUE)) + 1

plots_list <- lapply(1:nrow(plot_data), function(i) {
  row <- plot_data[i, ]
  
  ggplot(row, aes(x = HR, y = 1)) +
    # 添加 HR=1 基准线
    geom_vline(xintercept = 1, linetype = "dotted", color = "gray50", size = 0.4) +
    # 绘制误差线
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0, 
                   color = row$plot_color, size = 0.8) +
    # 绘制 HR 点 (菱形代表联合效应)
    geom_point(size = 2.5, shape = 18, color = row$plot_color) +
    # 坐标轴对齐
    scale_x_log10(limits = c(x_min, x_max)) +
    theme_void() +
    theme(plot.margin = margin(0, 5, 0, 5))
})

# 生成 X 轴刻度尺作为脚注或底部行
axis_plot <- ggplot() +
  scale_x_log10(limits = c(x_min, x_max), breaks = c(1, 2, 4, 8)) +
  theme_minimal() +
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(),
        axis.text.y = element_blank(), 
        axis.line.x = element_line(color = "black"),
        plot.margin = margin(0, 5, 0, 5))

# =======================================================
# 3. 组装终极 Flextable (含森林图列)
# =======================================================
# 准备显示的表格列

table_to_show <- plot_data %>%
  select(PRS_Level, MRS_Level, Counts, HR_CI, P_val) %>%
  mutate(Forest = "") # 预留森林图位置

ft_forest <- flextable(table_to_show, 
                       col_keys = c("PRS_Level", "MRS_Level", "Counts", "Forest", "HR_CI", "P_val")) %>%
  # 🚨 核心步骤：将 ggplot 列表插入 Forest 列
  compose(j = "Forest", 
          value = as_paragraph(gg_chunk(value = plots_list, width = 1.5, height = 0.25))) %>%
  # 添加底部刻度轴
  add_footer_row(
    values = list(PRS_Level="", MRS_Level="", Counts="", Forest="axis", HR_CI="", P_val=""),
    colwidths = c(1, 1, 1, 1, 1, 1)
  ) %>%
  compose(i = 1, j = "Forest", 
          value = as_paragraph(gg_chunk(value = list(axis_plot), width = 1.5, height = 0.3)), 
          part = "footer") %>%
  # 格式美化
  set_header_labels(
    PRS_Level = "Genetic Risk (PRS)",
    MRS_Level = "Metabolic Risk (MRS)",
    Counts = "No. (Total/Cases)",
    Forest = "Hazard Ratio (Visual)",
    HR_CI = "HR (95% CI)",
    P_val = "P Value"
  ) %>%
  merge_v(j = "PRS_Level") %>%
  theme_booktabs() %>%
  align(align = "center", part = "all") %>%
  valign(valign = "center", j = 1:2) %>%
  autofit() %>%
  # 显著性 P 值加粗
  bold(i = ~ P_val == "<0.001" | (!is.na(as.numeric(P_val)) & as.numeric(P_val) < 0.05), j = "P_val")

# =======================================================
# 4. 导出至 Word
# =======================================================
doc_forest <- read_docx() %>%
  body_add_par("Table: Joint Effect of Genetic and Metabolic Risk with Visual Forest Plot", style = "heading 1") %>%
  body_add_par("The forest plot illustrates the hazard ratios (diamonds) and 95% confidence intervals (horizontal lines) for the risk of incident atrial fibrillation.", style = "Normal") %>%
  body_add_flextable(ft_forest)

print(doc_forest, target = "DoubleHit_ForestTable_Final.docx")

cat("✅ 带有森林图的联合风险表格已生成：DoubleHit_ForestTable_Final.docx\n")
# (运行上一条回答中的 ggplot 画图代码即可，数据输入就是这里的 res_plot)
library(survival)
library(ggplot2)
library(dplyr)

res_plot <- final_joint_table %>%
  
  rename(PRS = PRS_Level, MRS = MRS_Level) %>%
  
  mutate(
    
    PRS = factor(PRS, levels = c("Low PRS", "Mid PRS", "High PRS")),
    
    MRS = factor(MRS, levels = c("Low MRS", "Mid MRS", "High MRS"))
    
  )

# ==============================================================================
# 2. 绘制顶级期刊审美的 3x3 联合风险阶梯图
# ==============================================================================
p_joint_bar <- ggplot(res_plot, aes(x = PRS, y = HR, fill = MRS)) +
  
  # 1. 添加 HR = 1 的基准参考线 (极其专业的设计)
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.8) +
  
  # 2. 绘制主体柱子：增加极细的黑色描边，提升锐利度
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.72, 
           color = "black", linewidth = 0.3, alpha = 0.95) +
  
  # 3. 绘制误差线：线条略细，颜色加深，避免喧宾夺主
  geom_errorbar(aes(ymin = Lower, ymax = Upper), 
                position = position_dodge(0.8), width = 0.25, linewidth = 0.6, color = "grey20") +
  
  # 4. 动态文本标签：精准定位在误差线最高点 (Upper) 上方 0.3 个单位，绝不重叠
  geom_text(aes(y = Upper + 0.3, label = sprintf("%.2f", HR)), 
            position = position_dodge(0.8), size = 3.5, fontface = "bold", color = "black") +
  
  # 5. 顶级心血管大刊配色 (安全灰 -> 警示橙 -> 危险红)
  scale_fill_manual(name = "Metabolic Clock (MRS) Stratification", 
                    values = c("#E0E0E0", "#F39B7F", "#E64B35")) +
  
  # 6. Y轴动态扩展：底部紧贴X轴，顶部留出 10% 的呼吸空间给标签
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  
  # 7. 标题与标签重构
  labs(title = "Double-Hit Model: Joint Impact of Genetics and Metabolism on Incident AF",
       y = "Hazard Ratio (95% CI)", 
       x = "Genetic Predisposition (Polygenic Risk Score)") +
  
  # 8. 经典医学主题微调 (Theme Classic)
  theme_classic() +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5, margin = margin(b = 20)),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 15)),
    axis.title.y = element_text(size = 13, face = "bold", margin = margin(r = 15)),
    axis.text.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.8), # 加粗坐标轴线
    axis.ticks = element_line(color = "black", linewidth = 0.8),
    panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"), # 加入轻微的横向引导线
    legend.position = "top",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

# 保存为高清 PDF
ggsave("Figure_Joint_Risk_DoubleHit_Premium.pdf", plot = p_joint_bar, width = 9, height = 6.5)

cat("✅ 顶刊级别 3x3 联合风险阶梯图已保存：Figure_Joint_Risk_DoubleHit_Premium.pdf\n")


# ==============================================================================
# 🌟 [封神大结局 V13-Final-Adjusted] 基因-临床代谢评分 双重打击模型
# ==============================================================================
cat("⏳ 正在启动 V13 调优引擎：调整比例、坐标轴与标签位置...\n")

library(dplyr)
library(stringr)
library(ggplot2)
library(ggalluvial) 
library(patchwork)  
library(ggnewscale)

# ——————————————————————————————————————————————————————————————————————————————
# 1. 数据整理与对齐 (基于 4 倍截断逻辑)
# ——————————————————————————————————————————————————————————————————————————————
df_9_align <- final_joint_table %>%
  rename(PRS_Label = PRS_Level, MRS_Label = MRS_Level) %>%
  mutate(
    N_Total = as.numeric(str_extract(Counts, "^[0-9]+")),
    Events = as.numeric(str_extract(Counts, "[0-9]+$")),
    # 🚨 新增下面这两行：把 "<0.001" 这种字符变成数字，并计算 -log10(P)
    numeric_p = ifelse(P_val == "<0.001" | P_val == "-", 1e-4, as.numeric(P_val)),
    significance = -log10(numeric_p),
    PRS_Label = factor(PRS_Label, levels = c("Low PRS", "Mid PRS", "High PRS")),
    MRS_Label = factor(MRS_Label, levels = c("Low MRS", "Mid MRS", "High MRS")),
    
    estimate = as.numeric(HR),
    conf.low = as.numeric(Lower),
    conf.high = as.numeric(Upper)
  ) %>%
  arrange(MRS_Label, PRS_Label) %>% 
  mutate(
    weight = sqrt(N_Total), 
    ymax = cumsum(weight), 
    ymin = ymax - weight, 
    ycenter = (ymin + ymax) / 2,
    
    # 🚨 坐标轴上限改为 4
    plot_conf_high = ifelse(conf.high > 4, 4, conf.high),
    is_truncated = conf.high > 4 
  )

MAX_Y <- sum(df_9_align$weight)

# ——————————————————————————————————————————————————————————————————————————————
# 2. 🎨 色彩引擎 (保持 V13 奢华配色)
# ——————————————————————————————————————————————————————————————————————————————
node_colors <- c(
  "Low PRS" = "#7B8D9E", "Mid PRS" = "#4A5A6A", "High PRS" = "#1D2633",  
  "Low MRS" = "#7B8D9E", "Mid MRS" = "#4A5A6A", "High MRS" = "#1D2633"
)
fresh_gradient <- c("#313695", "#E0F3F8", "#FEE090", "#F48FB1", "#EC407A", "#C2185B", "#880E4F")

# ——————————————————————————————————————————————————————————————————————————————
# Panel A: 2D 桑基河流图 (占 4 成宽度)
# ——————————————————————————————————————————————————————————————————————————————
sankey_data <- df_9_align %>% select(PRS_Label, MRS_Label, freq = weight, estimate)

p_sankey <- ggplot(sankey_data, aes(y = freq, axis1 = PRS_Label, axis2 = MRS_Label)) +
  geom_alluvium(aes(fill = estimate), width = 1/6, alpha = 0.55, curve_type = "quintic", reverse = FALSE) +
  scale_fill_gradientn(colors = fresh_gradient, values = c(0, 0.1, 0.25, 0.4, 0.6, 0.8, 1), guide = "none") +
  new_scale_fill() +
  
  geom_stratum(aes(fill = after_stat(stratum)), width = 1/6, color = "grey30", linewidth = 0.2, alpha = 0.95, reverse = FALSE) +
  scale_fill_manual(values = node_colors, guide = "none") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), color = "white", fontface = "bold", size = 4.7, reverse = FALSE) + 
  
  scale_x_discrete(limits = c("Genetic Risk\n(PRS)", "Metabolic Risk\n(MRS)"), 
                   expand = expansion(mult = c(0.1, 0.02))) +
  scale_y_continuous(limits = c(0, MAX_Y), expand = c(0, 0)) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.x = element_text(face = "bold", size = 14, color = "black", vjust = 0, margin = margin(t = 15, b = 10)),
    plot.margin = margin(t = 20, r = 0, b = 10, l = 10)
  )

# ——————————————————————————————————————————————————————————————————————————————
# Panel B: 瀑布图 (占 6 成宽度 + HR 标签左挪)
# ——————————————————————————————————————————————————————————————————————————————
p_waterfall <- ggplot(df_9_align, aes(y = ycenter)) +
  # 风险色带
  geom_rect(aes(ymin = ymin, ymax = ymax, xmin = 0.1, xmax = estimate, fill = estimate), alpha = 0.6) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 1) +
  
  # 误差线 (适配截断)
  geom_linerange(aes(xmin = conf.low, xmax = plot_conf_high, color = estimate), alpha = 0.9, linewidth = 1.2) +
  
  # 截断箭头 (上限为 4)
  geom_segment(data = filter(df_9_align, is_truncated),
               aes(x = 3.6, xend = 4, y = ycenter, yend = ycenter, color = estimate),
               arrow = arrow(length = unit(0.2, "cm"), type = "closed"), linewidth = 1.2) +
  
  # 🚨 优化：白字 HR 标签通过乘以 0.85 整体向左偏移，避免挡住钻石
  geom_text(aes(x = estimate * 0.85, label = sprintf("%.2f", estimate)), 
            color = "white", fontface = "bold", size = 7) +
  
  # 钻石点 (大小映射为显著性)
  geom_point(aes(x = estimate, color = estimate, size = significance), shape = 18, alpha = 1) +
  
  scale_fill_gradientn(
    colors = fresh_gradient, 
    values = c(0, 0.1, 0.25, 0.4, 0.6, 0.8, 1), 
    name = "Adjusted HR Mapping",
    guide = guide_colorbar(order = 1, barwidth = 1.5, barheight = 8, frame.colour = "grey20", ticks.colour = "black")
  ) +
  scale_color_gradientn(colors = fresh_gradient, values = c(0, 0.1, 0.25, 0.4, 0.6, 0.8, 1), guide = "none") +
  # 🚨 图例标题改为显著性，调整 breaks 以适配 -log10(P) 常见范围 (例如 1.3 代表 P≈0.05)
  scale_size_continuous(range = c(3, 10), name = "Significance\n(-log10 P)", breaks = c(1.3, 2, 3, 4), labels = c("1.3 (P=0.05)", "2.0", "3.0", "≥4.0"), guide = guide_legend(order = 2)) +
  
  # 🚨 X轴上限设为 4
  scale_x_continuous(trans = "log10", breaks = c(0.5, 1, 2, 4)) +
  scale_y_continuous(limits = c(0, MAX_Y), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.4, 4), expand = FALSE) +
  
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
    axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.8),
    axis.ticks.x = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length.x = unit(2.5, "mm"),
    axis.text.x = element_text(face = "bold", size = 12, color = "black", margin = margin(t = 8)),
    axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 15)),
    panel.grid.major.x = element_line(color = "grey80", linetype = "dotted"),
    legend.position = "right", 
    legend.box = "vertical",
    legend.background = element_rect(fill = alpha("white", 0.9), color = "grey80", linewidth = 0.5),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 0)
  ) +
  labs(x = "Hazard Ratio (95% CI, Log Scale)")

# ——————————————————————————————————————————————————————————————————————————————
# 终极合并：黄金比例 4:6 分割
# ——————————————————————————————————————————————————————————————————————————————
final_plot <- p_sankey + p_waterfall + 
  plot_layout(widths = c(4, 6)) + # 🚨 黄金比例分割
  plot_annotation(
    title = "Double-Hit Model: Joint Impact of Genetics and Metabolic Phenotypes on AF",
    subtitle = "Multidimensional Risk Stratification encoded by HR Gradient | Y-axis scaled by Square Root of N",
    theme = theme(plot.title = element_text(face = "bold", size = 19, hjust = 0.5),
                  plot.subtitle = element_text(face = "italic", color = "#C00000", size = 14, hjust = 0.5))
  )

ggsave("Figure_Final_V13_4-6_Ratio_Adjusted.pdf", final_plot, width = 16, height = 8, device = "pdf")

cat("✅ 终极修正版完成！比例 4:6，HR上限为 4，白字标签已左移避让钻石。\n")
library(survival)
library(dplyr)
library(stringr)
library(flextable)
library(officer)

# ==============================================================================
# 核心函数：计算 MRS 与 PRS 的完整加性交互 (2x2) 与乘性交互
# ==============================================================================
run_mrs_prs_interaction_analysis <- function(imp_obj) {
  
  cat("\n====== 正在启动 MRS x PRS 全面交互分析 ======\n")
  
  # --- A. 乘性交互 P 值 (Multiplicative Interaction) ---
  cat("⏳ 正在计算乘性交互 P 值...\n")
  fit_multi <- with(imp_obj, coxph(as.formula(paste(
    "Surv(duration_updated, event_af_updated) ~ MRS_3 * PRS_3 +", 
    paste(covs_joint, collapse = " + ")
  ))))
  pooled_multi <- pool(fit_multi)
  p_multi <- summary(pooled_multi) %>% 
    filter(term == "MRS_3:PRS_3") %>% 
    pull(p.value)
  
  # --- B. 加性交互计算 (Additive Interaction) ---
  cat("⏳ 正在计算 4 种风险组合的加性交互 (RERI & AP)...\n")
  m <- imp_obj$m
  
  # 存储矩阵: [插补次数, 4种组合, 3个值(RERI_Est, RERI_Var, AP_Est)]
  res_array <- array(NA, dim = c(m, 4, 3)) 
  dimnames(res_array)[[2]] <- c("MidMRS_MidPRS", "HighMRS_MidPRS", "MidMRS_HighPRS", "HighMRS_HighPRS")
  dimnames(res_array)[[3]] <- c("RERI_Est", "RERI_Var", "AP_Est")
  
  calc_additive <- function(model_single, term_11, term_10, term_01) {
    b <- coef(model_single)
    v <- vcov(model_single)
    
    i11 <- grep(term_11, names(b), fixed = TRUE)
    i10 <- grep(term_10, names(b), fixed = TRUE)
    i01 <- grep(term_01, names(b), fixed = TRUE)
    
    if(length(i11)==0 | length(i10)==0 | length(i01)==0) return(c(NA, NA, NA))
    
    hr11 <- exp(b[i11]); hr10 <- exp(b[i10]); hr01 <- exp(b[i01])
    
    # RERI
    reri_val <- hr11 - hr10 - hr01 + 1
    # Delta Method Variance
    grad <- rep(0, length(b))
    grad[i11] <- hr11; grad[i10] <- -hr10; grad[i01] <- -hr01
    reri_var <- t(grad) %*% v %*% grad
    # AP
    ap_val <- reri_val / hr11
    
    return(c(reri_val, reri_var, ap_val))
  }
  
  for(i in 1:m) {
    curr_mod <- model_list[[i]] # 确保你的环境中已经有之前跑好的 model_list
    
    # 组合 1：Mid MRS + Mid PRS
    res_array[i, "MidMRS_MidPRS", ] <- calc_additive(curr_mod, 
                                                     "Mid PRS + Mid MRS", 
                                                     "Low PRS + Mid MRS", 
                                                     "Mid PRS + Low MRS")
    # 组合 2：High MRS + Mid PRS
    res_array[i, "HighMRS_MidPRS", ] <- calc_additive(curr_mod, 
                                                      "Mid PRS + High MRS", 
                                                      "Low PRS + High MRS", 
                                                      "Mid PRS + Low MRS")
    # 组合 3：Mid MRS + High PRS
    res_array[i, "MidMRS_HighPRS", ] <- calc_additive(curr_mod, 
                                                      "High PRS + Mid MRS", 
                                                      "Low PRS + Mid MRS", 
                                                      "High PRS + Low MRS")
    # 组合 4：High MRS + High PRS
    res_array[i, "HighMRS_HighPRS", ] <- calc_additive(curr_mod, 
                                                       "High PRS + High MRS", 
                                                       "Low PRS + High MRS", 
                                                       "High PRS + Low MRS")
  }
  
  # --- C. 池化汇总 (Rubin's Rules) ---
  pool_it <- function(ests, vars, aps) {
    Q_bar <- mean(ests, na.rm=TRUE)
    U_bar <- mean(vars, na.rm=TRUE)
    B <- var(ests, na.rm=TRUE)
    se <- sqrt(U_bar + (1 + 1/m) * B)
    
    reri_str <- sprintf("%.2f (%.2f, %.2f)", Q_bar, Q_bar - 1.96*se, Q_bar + 1.96*se)
    ap_str <- sprintf("%.2f", mean(aps, na.rm=TRUE))
    return(list(reri = reri_str, ap = ap_str))
  }
  
  p_mm <- pool_it(res_array[,"MidMRS_MidPRS",1], res_array[,"MidMRS_MidPRS",2], res_array[,"MidMRS_MidPRS",3])
  p_hm <- pool_it(res_array[,"HighMRS_MidPRS",1], res_array[,"HighMRS_MidPRS",2], res_array[,"HighMRS_MidPRS",3])
  p_mh <- pool_it(res_array[,"MidMRS_HighPRS",1], res_array[,"MidMRS_HighPRS",2], res_array[,"MidMRS_HighPRS",3])
  p_hh <- pool_it(res_array[,"HighMRS_HighPRS",1], res_array[,"HighMRS_HighPRS",2], res_array[,"HighMRS_HighPRS",3])
  
  # --- D. 组装数据框 ---
  final_df <- data.frame(
    Level = c("Moderate Risk", "High Risk"),
    Mid_RERI = c(p_mm$reri, p_hm$reri),
    Mid_AP   = c(p_mm$ap, p_hm$ap),
    High_RERI = c(p_mh$reri, p_hh$reri),
    High_AP   = c(p_mh$ap, p_hh$ap)
  )
  
  return(list(table = final_df, p_multi = p_multi))
}

# 运行分析
interaction_results <- run_mrs_prs_interaction_analysis(imp_final_joint)
cat("✅ 交互分析计算完毕！\n")


# ==============================================================================
# 创建顶刊级别 Flextable 表格 (2x2 双层表头格式)
# ==============================================================================
cat("⏳ 正在生成精美 Word 表格...\n")

ft_data <- interaction_results$table
p_val_multi <- interaction_results$p_multi

ft_inter <- flextable(ft_data) %>%
  # 1. 设置底层表头 (列名)
  set_header_labels(
    Level = "Metabolic Risk Level (MRS)",
    Mid_RERI = "RERI (95% CI)",
    Mid_AP   = "AP (95% CI)",
    High_RERI = "RERI (95% CI)",
    High_AP   = "AP (95% CI)"
  ) %>%
  
  # 2. 添加顶层表头 (Moderate PRS / High PRS)
  add_header_row(
    values = c("", "Moderate PRS", "High PRS"),
    colwidths = c(1, 2, 2) # 第一列不合并，后面每两列合并
  ) %>%
  
  # 3. 添加大标题
  add_header_lines("Table: Additive interaction between Metabolic Clock (MRS) and PRS on AF incidence") %>%
  
  # 4. 样式美化 (三线表风格)
  theme_booktabs() %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "body") %>%
  
  # 5. 调整列宽使排版更舒展
  width(j = 1, width = 1.8) %>%
  width(j = 2:5, width = 1.3) %>%
  
  # 6. 修正边框：在 Moderate/High PRS 下方加一条漂亮的横线
  hline(i = 1, j = 2:5, border = fp_border(color = "black", width = 1), part = "header") %>%
  
  # 7. 🚨 核心：在表格外部（注脚）添加乘性交互 P 值和缩写说明
  add_footer_lines(paste0("P for multiplicative interaction = ", 
                          ifelse(p_val_multi < 0.001, "<0.001", sprintf("%.3f", p_val_multi)), ".")) %>%
  add_footer_lines("RERI, Relative Excess Risk due to Interaction; AP, Attributable Proportion.") %>%
  add_footer_lines("Models adjusted for age, sex, ethnicity, lifestyle factors, comorbidities, and 10 genetic PCs.") %>%
  
  # 注脚左对齐并稍微调小字号
  align(align = "left", part = "footer") %>%
  fontsize(size = 9, part = "footer")

# 在 RStudio Viewer 中预览
print(ft_inter)

# 导出为 Word
save_as_docx(ft_inter, path = "MRS_PRS_Interaction_Table_Premium.docx")
cat("✅ 终极版交互表格已生成：MRS_PRS_Interaction_Table_Premium.docx\n")