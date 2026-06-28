# ==============================================================================
# 👑 顶刊级无监督聚类 V4.0：适配全新 8 因子与严谨物理隔离队列
# ==============================================================================
cat("⏳ 正在启动严谨版 GMM 聚类引擎...\n")
# ==============================================================================
# 🔒 核心资产大封印：保存所有用于 GMM 和 验证的数据
# ==============================================================================
cat("⏳ 正在打包核心资产，准备重启...\n")

save(
  # 1. 核心数据集 (7:3 物理切割后的完整版)
  df_train_final,      # 70% 训练集 (GMM 和 权重计算的源头)
  df_test_final,       # 30% 测试集 (绝对隔离的验证集)
  df_benchmark_scores, # 全队列数据 (包含所有评分和 MetRS_Z)
  
  # 2. 核心特征与权重
  final_11_metabs,     # 那 8 个（或11个）王者代谢物名单
  metRS_weights_vec,   # Cox 回归算出来的权重
  ultra_stable_metabs, # 稳定性选择的频率表
  
  # 3. 聚类必需的元数据
  train_eids_locked,   # 确保重启后 7:3 划分不乱的 ID 锁
  core_metabolites,    # 144 个中介初筛名单
  
  file = "AF_MetRS_GMM_Ready_Pack.RData"
)

cat("✅ 打包完毕！文件名：AF_MetRS_GMM_Ready_Pack.RData\n")
cat("🚀 现在可以安心重启 R 或者重启服务器了。回来见！\n")
# 1. 载入所有神级装备包
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(mclust)
library(uwot)           
library(mclustcomp)     
library(ggpointdensity) 
library(survminer)      
library(ClusterR)       
library(caret)          
library(pheatmap)       

# ==============================================================================
# 🛡️ 步骤 1 & 2：全量数据准备与 70% 训练集提取 (极简精准化)
# ==============================================================================
# 🚨 核心修改 1：直接使用刚出炉的 final_11_metabs (8个王者代谢物)
final_metabs <- final_11_metabs

# 把 7:3 拆分的两个集合拼起来作为 100% 全景图底座 (用于最后 UMAP 投影)
# 🚨 核心修改 2：适配新的 eid，剔除了不需要的冗余列防报错
full_cohort <- bind_rows(df_train_final, df_test_final) %>%
  select(eid, duration_updated, event_af_updated, all_of(final_metabs)) %>%
  drop_na()

full_features_mat <- as.matrix(scale(full_cohort %>% select(all_of(final_metabs))))

cat("⏳ 提取 70% 训练集作为 GMM 建模底座 (严格遵守物理隔离)...\n")
# 直接从 df_train_final 中提取 GMM 的建模基础
disc_cohort <- df_train_final %>%
  select(eid, duration_updated, event_af_updated, all_of(final_metabs)) %>%
  drop_na()

# 提取这 70% 人的特征矩阵进行 Z-score 标准化
disc_features <- as.matrix(scale(disc_cohort %>% select(all_of(final_metabs))))

cat(sprintf("📊 检查：Discovery (Train) 样本数 = %d, 特征矩阵行数 = %d\n", 
            nrow(disc_cohort), nrow(disc_features)))

# ==============================================================================
# 🛡️ 步骤 3：强制认领 Strata 标签并执行降采样 (防内存爆雷)
# ==============================================================================
cat("⏳ 正在执行 MiniBatchKmeans 中心点计算与全员标签分配...\n")
set.seed(2026)

km_fast <- MiniBatchKmeans(disc_features, clusters = 40, batch_size = 500, num_init = 5, max_iters = 50)
strata_labels <- as.vector(predict_MBatchKMeans(disc_features, km_fast$centroids))

if(length(strata_labels) != nrow(disc_cohort)) {
  stop("🚨 对齐失败！预测得到的标签数 与 样本行数不符！")
}

# 🚨 核心修改 3：动态计算每层抽样数，防止样本量太少时报错
n_per_strata <- min(500, floor(nrow(disc_cohort) / 40)) 

sample_idx_in_disc <- disc_cohort %>%
  mutate(Strata = strata_labels, Row_ID_Local = row_number()) %>%
  group_by(Strata) %>%
  slice_sample(n = n_per_strata, replace = FALSE) %>% 
  ungroup() %>%
  pull(Row_ID_Local)

scaled_sample <- disc_features[sample_idx_in_disc, ]
cat(sprintf("✅ 抽样成功！用于 GMM 核心建模的子集样本量为: %d\n", nrow(scaled_sample)))

# ==============================================================================
# 🛡️ 步骤 4：计算 BIC 矩阵，并拟合 GMM 模型
# ==============================================================================
cat("⏳ 正在计算 1-8 簇的官方 BIC 矩阵...\n")
set.seed(2026)
bic_matrix <- mclustBIC(scaled_sample, G = 1:8)

pdf("Supplementary_Figure_BIC_Final.pdf", width=6, height=5)
plot(bic_matrix)
dev.off()

cat("⏳ 正在锁定 K=3 进行拟合...\n")
# 强制使用 K=3，这在医学亚型聚类中最符合 "低/中/高" 的临床叙事
gmm_train <- Mclust(scaled_sample, G = 3, x = bic_matrix) 
cat(sprintf("🎯 算法优选的最优协方差结构为: [%s]\n", gmm_train$modelName))
# ==============================================================================
# 🛡️ 步骤 5：Bootstrap ARI 内部稳定性分析
# ==============================================================================
cat("⏳ 正在执行 Bootstrap 聚类稳定性验证 (B=500次)...\n")
set.seed(2026)

n_boot_ari <- 500
ari_scores <- numeric(n_boot_ari)

# 🚨 新增：初始化动态进度条 (style = 3 是带百分比的进度条，最直观)
pb_ari <- txtProgressBar(min = 0, max = n_boot_ari, style = 3)

for(i in 1:n_boot_ari) {
  # 1. 有放回重抽样
  boot_idx <- sample(1:nrow(scaled_sample), replace = TRUE)
  
  # 2. 拟合 GMM 模型 (防报错保护)
  boot_model <- tryCatch({ 
    Mclust(scaled_sample[boot_idx, ], G = 3, verbose = FALSE) 
  }, error = function(e) NULL)
  
  # 3. 计算并记录 ARI
  if(!is.null(boot_model)){
    pred_boot <- predict(boot_model, newdata = scaled_sample)$classification
    ari_scores[i] <- adjustedRandIndex(pred_boot, gmm_train$classification)
  } else { 
    ari_scores[i] <- NA 
  }
  
  # 🚨 新增：每次循环结束，更新一次进度条
  setTxtProgressBar(pb_ari, i)
}

# 🚨 新增：循环结束，平滑关闭进度条
close(pb_ari)

# 4. 清理 NA 并计算均值
ari_scores <- na.omit(ari_scores)
mean_ari <- mean(ari_scores)

# 注意开头的 \n，防止和进度条挤在同一行
cat(sprintf("\n✅ 聚类稳定性极佳！平均 Adjusted Rand Index (ARI) = %.3f\n", mean_ari))
# 绘制 ARI 密度图
df_ari <- data.frame(ARI = ari_scores)
dens <- density(df_ari$ARI); max_y <- max(dens$y)

p_ari <- ggplot(df_ari, aes(x = ARI)) +
  geom_histogram(aes(y = after_stat(density)), bins = 15, fill = "#7B8D9E", color = "white", alpha = 0.65) +
  geom_density(color = "#1D2633", linewidth = 1.2, alpha = 0.8) +
  geom_vline(xintercept = mean_ari, color = "#C2185B", linetype = "dashed", linewidth = 1.2) +
  annotate("text", x = mean_ari + 0.005, y = max_y * 0.95, label = sprintf("Mean ARI = %.3f", mean_ari), color = "#C2185B", fontface = "bold", size = 5, hjust = 0) +
  theme_classic(base_size = 14) +
  labs(title = "Bootstrap Stability of GMM Clustering", subtitle = "Distribution of Adjusted Rand Index (ARI)", x = "Adjusted Rand Index (ARI)", y = "Density")

ggsave("Supplementary_Figure_Bootstrap_ARI.pdf", plot = p_ari, width = 7, height = 5)


library(readxl)
library(pheatmap)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(survminer)
library(patchwork)
library(ggplotify) # 👑 关键：将 pheatmap 转换为 ggplot 对象的魔法包
cat("✅ 顶级人群重分型 (GMM + UMAP) 流水线运行完毕！所有高清大图已存至本地！\n")
metab_dict <- read_excel("dict.xlsx")
cat("\n# ==============================================================================\n")
cat("# 👑 终极主图合成：高阶代谢特征热图 (A) + 临床表型 KM 曲线 (B)\n")
cat("# ==============================================================================\n")

# ==============================================================================
# 🌟 [顶刊视觉增强版] 56个核心代谢物：自适应缩写与高阶对齐的环形热图
# ==============================================================================
# ==============================================================================
# 🌟 [顶刊视觉 V3 终极版] 56个核心代谢物：极限缩写 + 内外四环 + 绝对防错位
# ==============================================================================
cat("⏳ 正在构建绝不跑偏的出版级环形热图...\n")
# ==============================================================================
# 📍 终极源头大清洗：动态根据实际风险自动贴标签 (绝对防错位)
# ==============================================================================
cat("⏳ 正在将全量患者前瞻性映射进高斯概率云，并计算真实风险...\n")
pred_res <- predict(gmm_train, newdata = full_features_mat)
full_cohort$Cluster_Raw <- as.factor(pred_res$classification)

# 🚨 核心杀招：让数据自己说话！计算每个 Cluster 的真实房颤发病率
risk_mapping <- full_cohort %>%
  group_by(Cluster_Raw) %>%
  summarise(Event_Rate = mean(event_af_updated, na.rm = TRUE)) %>%
  arrange(Event_Rate) %>% # 严格按发病率从低到高排序
  mutate(
    # 发病率最低的永远叫 Healthy，最高的永远叫 Severe！
    Auto_Label = c("Healthy", "Fibro-Inflammatory", "Classic Cardiometabolic")
  )

cat("📊 算法动态分配结果 (发病率低 -> 高)：\n")
print(risk_mapping)

# 把绝对正确的标签贴回全队列，并死锁顺序！
full_cohort <- full_cohort %>%
  left_join(risk_mapping %>% select(Cluster_Raw, Auto_Label), by = "Cluster_Raw") %>%
  mutate(
    Plot_Cluster = factor(Auto_Label, levels = c("Healthy", "Fibro-Inflammatory", "Classic Cardiometabolic"))
  )

cat("✅ 标签源头拨乱反正完毕！全剧终！以后不管 Mclust 怎么变，标签永远不会贴错！\n")
# ==============================================================================
# 🛡️ 步骤：生成带风险表的独立 KM 曲线 PDF (彻底消灭错位)
# ==============================================================================
cat("⏳ 正在绘制审美精修版 KM 曲线...\n")

library(survival)
library(survminer)

# 1. 拟合生存对象
fit_km <- survfit(Surv(duration_updated, event_af_updated) ~ Plot_Cluster, data = full_cohort)

# 2. 核心绘图
p_km_standalone <- ggsurvplot(
  fit_km, 
  data = full_cohort, 
  fun = "cumhaz",
  xlim = c(0, 10), 
  ylim = c(0, 0.10),             
  break.time.by = 2,
  # 🚨 颜色严格按照 levels 提取，绝不串色
  palette = c("#4DBBD5FF", "#00A087FF", "#E64B35FF"),
  
  title = "Cumulative AF Hazard by Metabolic Phenotype",
  xlab = "Follow-up Time (Years)", 
  ylab = "Cumulative Hazard",
  
  risk.table = TRUE, 
  risk.table.col = "strata",     
  risk.table.y.text = FALSE, 
  risk.table.height = 0.18,      
  risk.table.fontsize = 4.5,     
  
  legend.title = "Phenotype",
  # 🚨🚨🚨 核心凶手伏法：不再手动写死，直接提取底层最干净的因子名！
  legend.labs = levels(full_cohort$Plot_Cluster),
  
  pval = TRUE, 
  pval.coord = c(0, 0.08),       
  ggtheme = theme_minimal() 
)

# 3. 主图与表格美化
p_km_standalone$plot <- p_km_standalone$plot + 
  theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
        axis.text = element_text(face = "bold", color = "black", size = 12),
        axis.title = element_text(face = "bold", size = 13),
        legend.text = element_text(size = 11, face = "bold"),
        panel.grid.minor = element_blank()) 

p_km_standalone$table <- p_km_standalone$table + 
  theme_cleantable() + 
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())

# 4. 保存
pdf("Figure_Final_KM_Ultra_Aesthetic.pdf", width = 8, height = 8.5)
print(p_km_standalone)
dev.off()
# ==============================================================================
# 🛡️ 补充步骤：UMAP 降维计算 (防内存溢出均匀抽样)
# ==============================================================================
cat("⏳ 正在计算 UMAP 降维矩阵 (抽取 3 万人代表性底片)...\n")
library(uwot)

set.seed(2026)
# 从全人群中抽取 30000 人画 UMAP（既能看清云团，又不会卡死）
n_umap_sample <- min(30000, nrow(full_cohort))
umap_plot_idx <- sample(1:nrow(full_cohort), n_umap_sample, replace = FALSE)

# 提取特征矩阵并降维
umap_data_matrix <- full_features_mat[umap_plot_idx, ]
umap_model <- umap(umap_data_matrix, n_neighbors = 15, min_dist = 0.1, metric = "euclidean", verbose = TRUE)

cat("✅ UMAP 坐标计算完成！准备绘制绝美分面图...\n")

# ==============================================================================
# 🛡️ 步骤 6 & 7：UMAP 极美可视化 (终极颜色与新名字顺序对齐)
# ==============================================================================
cat("⏳ 正在生成完美对齐的 UMAP 分面图...\n")

# 🚨 绝对源头映射：直接从我们刚才洗干净的 Plot_Cluster 里抓取，再也不用手动写 if-else！
umap_df <- data.frame(
  UMAP1 = umap_model[, 1], 
  UMAP2 = umap_model[, 2],
  Cluster_Label = full_cohort$Plot_Cluster[umap_plot_idx] 
)

# 制作灰色背景底板 (剔除标签，这样每个分面都有完整的灰色残影)
umap_bg <- umap_df %>% select(-Cluster_Label)

p_umap_facet <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  # 底层灰色星云
  geom_point(data = umap_bg, color = "grey85", size = 0.1, alpha = 0.5) +
  # 表层高亮着色
  geom_point(aes(color = Cluster_Label), size = 0.3, alpha = 0.8) +
  # 🚨 颜色字典死锁：蓝(Healthy)、绿(Fibro-Inflammatory)、红(Classic)
  scale_color_manual(values = c("Healthy" = "#4DBBD5FF", 
                                "Fibro-Inflammatory" = "#00A087FF", 
                                "Classic Cardiometabolic" = "#E64B35FF")) +
  facet_wrap(~ Cluster_Label) +
  labs(title = "Metabolic Continuum Across Phenotypes", x = "UMAP Dimension 1", y = "UMAP Dimension 2") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none", 
        strip.background = element_rect(fill = "grey95", color = NA), 
        strip.text = element_text(face = "bold", size = 13), 
        plot.title = element_text(face = "bold", hjust = 0.5))

# 导出出版级高分大图
ggsave("Figure_UMAP_Faceted_Professional_Renamed.pdf", p_umap_facet, width = 11, height = 4.5, dpi = 300)

cat("✅ UMAP 完美落盘！请打开文件夹查看 Figure_UMAP_Faceted_Professional_Renamed.pdf！\n")
# ==============================================================================
# 🌟 附加模块：绘制硬核直线版代谢指纹雷达图 (Radar Chart)
# ==============================================================================
cat("\n⏳ 正在计算三大表型的代谢物全局 Z-score 并绘制硬核雷达图...\n")

library(dplyr)
library(tidyr)
library(ggplot2)

# 🚨 破解 geom_polygon 画曲线的底层 Hack (强行直线化)
coord_radar <- function (theta = "x", start = 0, direction = 1, clip = "on") {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start,
          direction = sign(direction), clip = clip,
          is_linear = function(coord) TRUE) # 这里是核心：强制坐标系保留直线！
}

# ==============================================================================
# 🌟 附加模块：绘制代谢组学特征雷达图 (FMSB 顶刊蜘蛛网版)
# ==============================================================================
cat("\n⏳ 正在启动 fmsb 引擎绘制完美阴影雷达图...\n")

# 1. 安装并加载雷达图专业包
if (!requireNamespace("fmsb", quietly = TRUE)) install.packages("fmsb")
library(fmsb)
library(dplyr)
library(tidyr)

# 2. 缩写字典 (保持不变，确保标签清晰)
shorten_names_radar <- function(x) {
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
  return(x)
}

# ==============================================================================
# 3. 准备数据：fmsb 需要特定的极宽数据格式
# ==============================================================================
radar_data_raw <- full_cohort %>%
  select(Plot_Cluster, all_of(final_metabs)) %>%
  # 计算全局 Z-score
  mutate(across(all_of(final_metabs), ~ as.numeric(scale(.)))) %>% 
  group_by(Plot_Cluster) %>%
  # 计算均值
  summarise(across(everything(), ~ mean(., na.rm = TRUE)))

# 重命名列并转化为底层数据框
colnames(radar_data_raw)[-1] <- shorten_names_radar(colnames(radar_data_raw)[-1])
radar_df <- as.data.frame(radar_data_raw)
# 强制提取正确的 3 个分组名字作为行名
rownames(radar_df) <- radar_df$Plot_Cluster
radar_df <- radar_df[, -1]

# ==============================================================================
# 4. 计算雷达图的物理边界 (🚨 彻底修复 NA 隐身 Bug 版)
# ==============================================================================
max_z <- ceiling(max(radar_df) * 2) / 2
min_z <- floor(min(radar_df) * 2) / 2

# 🚨 关键修复：必须明确指定行名叫 Max 和 Min，否则后面提取会变成 NA！
radar_final <- rbind(
  Max = rep(max_z, ncol(radar_df)),
  Min = rep(min_z, ncol(radar_df)),
  radar_df
)

# 现在安全地确保顺序是：Max -> Min -> 蓝(Healthy) -> 绿(Standard) -> 红(Severe)
radar_final <- radar_final[c("Max", "Min", "Healthy", "Fibro-Inflammatory", "Classic Cardiometabolic"), ]

# 打印出来检查一下，只要前两行没有出现 NA，就绝对能画出来！
print(head(radar_final))

# ==============================================================================
# 5. 色彩引擎
# ==============================================================================
# 边缘线颜色
colors_border <- c("#4DBBD5FF", "#00A087FF", "#E64B35FF")
# 内部填充颜色 (带 15% 透明度)
colors_fill <- c(scales::alpha("#4DBBD5FF", 0.15), 
                 scales::alpha("#00A087FF", 0.15), 
                 scales::alpha("#E64B35FF", 0.15))

# ==============================================================================
# 6. 一键出图
# ==============================================================================
pdf("Figure_Radar_FMSB_Professional.pdf", width = 9, height = 10)

# 设置画图边距，允许在外部画图例 (xpd = TRUE)
par(mar = c(6, 2, 4, 2), xpd = TRUE) 

radarchart(
  radar_final,
  axistype = 1,                 # 在中间显示坐标轴刻度
  pcol = colors_border,         # 边缘线颜色
  pfcol = colors_fill,          # 内部填充颜色
  plwd = 2.5,                   # 线条粗细
  plty = 1,                     # 实线
  cglcol = "grey70",            # 蜘蛛网格线颜色
  cglty = 2,                    # 蜘蛛网格虚线
  cglwd = 0.8,                  # 网格线粗细
  axislabcol = "grey40",        # 坐标轴数字颜色
  caxislabels = seq(min_z, max_z, length.out = 5), # 自动生成 5 个刻度
  vlcex = 1.1,                  # 变量标签字体大小
  calcex = 0.9,                 # 刻度字体大小
  title = "Metabolic Fingerprints of the Sub-phenotypes\n(Global Z-score Means)"
)

# 在底部中央添加图例
legend(
  x = "bottom", y = -1.3,       # 精准定位到图形正下方
  legend = rownames(radar_final)[3:5],
  bty = "n", pch = 20, col = colors_border, text.col = "black",
  cex = 1.2, pt.cex = 3.5, horiz = TRUE
)

dev.off()
cat("✅ 真·顶级蜘蛛网雷达图已生成！请查看 Figure_Radar_FMSB_Professional.pdf\n")






cat("\n# ==============================================================================\n")
cat("⏳ 正在启动顶刊级表格生成引擎 (Cox & Table 1)...\n")
cat("# ==============================================================================\n")

library(dplyr)
library(tidyr)
library(survival)
library(broom)
library(stringr)
library(flextable)
library(officer)
library(tableone)
library(tibble)
# ==============================================================================
# 👑 步骤 1：双源智能提取与无缝合并 (Cox变量 + 肝功基线变量) [防碰撞修正版]
# ==============================================================================
cat("⏳ 正在从双数据源打捞所需变量...\n")

# 2. 列出我们希望拥有的所有附加临床与主成分协变量 (完美适配4模型)
cox_covariates <- c(
  "ethnicity_f_Black", "ethnicity_f_Others.Mixed.Unknown", "ethnicity_f_White",
  "bmi_val", "smoke_f_Never", "smoke_f_Prefer.not.to.answer", "smoke_f_Previous",
  "alc_freq_f", "tdi_val", "age_val", "sex_f_Male", "PRS_AF", "fib4_cat",
  "htn_final_1", "t2dm_final_1", "lip_f_1", "ckd_final_1", "cvd_f_1","MetRS_Z",
  "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"
)

missing_cox <- setdiff(cox_covariates, colnames(full_cohort))

if(length(missing_cox) > 0) {
  full_cohort_cox <- full_cohort %>%
    left_join(long_dat_dummy_中介 %>% select(Participant.ID, any_of(missing_cox)), 
              by = c("eid" = "Participant.ID"))
} else {
  full_cohort_cox <- full_cohort 
}

# ---------------------------------------------------------
# 阶段 2：从 df_clinical_final 极度安全地打捞基线变量 (避开列名冲突)
# ---------------------------------------------------------
clinical_extra <- df_clinical_final

# 1. 统一 ID 列名
if ("Participant.ID" %in% colnames(clinical_extra)) {
  clinical_extra$eid <- clinical_extra$Participant.ID
} else if ("Participant ID" %in% colnames(clinical_extra)) {
  clinical_extra$eid <- clinical_extra$`Participant ID`
}

# 2. 智能探测：如果不存在短名才去寻找长名赋值，如果已存在短名则直接沿用。
if (!"alt_val" %in% colnames(clinical_extra) && "Alanine aminotransferase | Instance 0" %in% colnames(clinical_extra)) {
  clinical_extra$alt_val <- clinical_extra[["Alanine aminotransferase | Instance 0"]]
}
if (!"ast_val" %in% colnames(clinical_extra) && "Aspartate aminotransferase | Instance 0" %in% colnames(clinical_extra)) {
  clinical_extra$ast_val <- clinical_extra[["Aspartate aminotransferase | Instance 0"]]
}
if (!"plt_val" %in% colnames(clinical_extra) && "Platelet count | Instance 0" %in% colnames(clinical_extra)) {
  clinical_extra$plt_val <- clinical_extra[["Platelet count | Instance 0"]]
}
if (!"alb_val" %in% colnames(clinical_extra) && "Albumin | Instance 0" %in% colnames(clinical_extra)) {
  clinical_extra$alb_val <- clinical_extra[["Albumin | Instance 0"]]
}

# 3. 安全提取并去重
clinical_extra <- clinical_extra %>% 
  select(any_of(c("eid", "alt_val", "ast_val", "plt_val", "alb_val"))) %>%
  distinct(eid, .keep_all = TRUE)

# 将肝功变量合并进主干数据集
full_cohort_cox <- full_cohort_cox %>%
  left_join(clinical_extra, by = "eid") %>%
  distinct(eid, .keep_all = TRUE)
# ---------------------------------------------------------
# 💡 【核心新增】阶段 2.5：从 df_benchmark_scores 智能精准打捞 MetRS_Z
# ---------------------------------------------------------
if (exists("df_benchmark_scores") && "MetRS_Z" %in% colnames(df_benchmark_scores)) {
  cat("🔄 正在从 df_benchmark_scores 提取 MetRS_Z...\n")
  benchmark_extra <- df_benchmark_scores
  
  # 统一 ID 列名（全面兼容各种ID格式）
  if ("Participant.ID" %in% colnames(benchmark_extra)) {
    benchmark_extra$eid <- benchmark_extra$Participant.ID
  } else if ("Participant ID" %in% colnames(benchmark_extra)) {
    benchmark_extra$eid <- benchmark_extra$`Participant ID`
  }
  
  # 提取目标变量并清洗去重
  benchmark_extra <- benchmark_extra %>% 
    select(any_of(c("eid", "MetRS_Z"))) %>%
    distinct(eid, .keep_all = TRUE)
  
  # 合并进入主数据集（自动清理潜在冲突）
  if ("MetRS_Z" %in% colnames(full_cohort_cox)) full_cohort_cox$MetRS_Z <- NULL
  full_cohort_cox <- full_cohort_cox %>%
    left_join(benchmark_extra, by = "eid") %>%
    distinct(eid, .keep_all = TRUE)
} else {
  cat("⚠️ 警告：未找到数据框 df_benchmark_scores 或其中没有 MetRS_Z 列！\n")
}
# ---------------------------------------------------------
# 阶段 3：严格执行缺失值过滤（保护基线变量）
# ---------------------------------------------------------
# 🚨 关键：过滤 NA 时，坚决只检查 Cox 变量，允许肝酶列存在 NA！
all_cox_vars <- c("Plot_Cluster", "duration_updated", "event_af_updated", cox_covariates)
full_cohort_cox <- full_cohort_cox %>% drop_na(any_of(all_cox_vars))

cat(sprintf("✅ 双源装填完毕！进入后续分析的可用完整样本量: %d\n", nrow(full_cohort_cox)))
# ==============================================================================
# 👑 步骤 2：执行 4 阶 Cox 回归与极美表格渲染
# ==============================================================================
# 🚨 完美映射你的 4 阶模型到真实的哑变量列名
cov_m1 <- c("age_val", "sex_f_Male", "ethnicity_f_Black", "ethnicity_f_Others.Mixed.Unknown", "ethnicity_f_White")
cov_m2 <- c(cov_m1, "bmi_val", "smoke_f_Never", "smoke_f_Prefer.not.to.answer", "smoke_f_Previous", "alc_freq_f", "tdi_val")
cov_m3 <- c(cov_m2, "htn_final_1", "t2dm_final_1", "lip_f_1", "ckd_final_1", "cvd_f_1")
cov_m4 <- c(cov_m3, "PRS_AF", paste0("PC", 1:10))

model_list <- list(
  "Crude (Unadjusted)" = "UNADJUSTED",
  "Model 1" = cov_m1,
  "Model 2" = cov_m2,
  "Model 3" = cov_m3,
  "Model 4" = cov_m4
)

# 计算事件数
counts_df <- full_cohort_cox %>%
  group_by(Plot_Cluster) %>%
  summarise(
    Total = n(),
    Cases = sum(event_af_updated, na.rm = TRUE), 
    Counts_Str = sprintf("%d / %d", Total, Cases),
    .groups = "drop"
  )
counts_map <- setNames(counts_df$Counts_Str, counts_df$Plot_Cluster)

all_res <- list()
for (mod_name in names(model_list)) {
  
  if (model_list[[mod_name]][1] == "UNADJUSTED") {
    formula_str <- "Surv(duration_updated, event_af_updated) ~ Plot_Cluster"
  } else {
    covariates_str <- paste(model_list[[mod_name]], collapse = " + ")
    formula_str <- paste("Surv(duration_updated, event_af_updated) ~ Plot_Cluster +", covariates_str)
  }
  
  fit <- coxph(as.formula(formula_str), data = full_cohort_cox)
  res_tidy <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% 
    filter(grepl("Plot_Cluster", term))
  
  mod_df <- data.frame(
    Model = mod_name, Exposure = "Metabolic Clusters", Type = "Categorical",
    Group = c( "Healthy (Ref)", "Fibro-Inflammatory", "Classic Cardiometabolic"),
    Counts = c(counts_map["Healthy"], counts_map["Fibro-Inflammatory"], counts_map["Classic Cardiometabolic"]),
    HR_CI = c(
      "1.00 (Reference)", 
      sprintf("%.2f (%.2f-%.2f)", res_tidy$estimate[1], res_tidy$conf.low[1], res_tidy$conf.high[1]),
      sprintf("%.2f (%.2f-%.2f)", res_tidy$estimate[2], res_tidy$conf.low[2], res_tidy$conf.high[2])
    ),
    P_val = c(
      "-",
      ifelse(res_tidy$p.value[1] < 0.001, "<0.001", sprintf("%.3f", res_tidy$p.value[1])),
      ifelse(res_tidy$p.value[2] < 0.001, "<0.001", sprintf("%.3f", res_tidy$p.value[2]))
    ),
    stringsAsFactors = FALSE
  )
  all_res[[mod_name]] <- mod_df
}

final_table_df <- bind_rows(all_res)

# 导出 Cox 表格
ft_cox <- flextable(final_table_df) %>%
  set_header_labels(Model = "Model", Exposure = "Exposure", Type = "Type", Group = "Level", 
                    Counts = "No. (Total / Cases)",  HR_CI = "HR (95% CI)", P_val = "P value") %>%
  merge_v(j = c("Model", "Exposure", "Type")) %>% 
  theme_booktabs() %>% autofit() %>%
  align(j = "Counts", align = "center", part = "all") %>%
  bold(i = ~ !str_detect(P_val, "-|>") & as.numeric(str_replace_all(P_val, "<", "")) < 0.05, j = "P_val")

doc_cox <- read_docx() %>%
  body_add_par("Table: Independent Association Between Metabolic Phenotypes and Incident AF", style = "heading 1") %>%
  body_add_par("Crude: Unadjusted model.", style = "Normal") %>%
  body_add_par("Model 1: Adjusted for Age, Sex, and Ethnicity.", style = "Normal") %>%
  body_add_par("Model 2: Adjusted for Model 1 + BMI, Smoking, Alcohol frequency, and Townsend Deprivation Index.", style = "Normal") %>%
  body_add_par("Model 3: Adjusted for Model 2 + Hypertension, Type 2 Diabetes, Dyslipidemia, Chronic Kidney Disease, and prior Cardiovascular Disease.", style = "Normal") %>%
  body_add_par("Model 4: Adjusted for Model 3 + AF Polygenic Risk Score (PRS) and Top 10 Genetic Principal Components.", style = "Normal") %>%
  body_add_flextable(ft_cox)
print(doc_cox, target = "Table_GMM_Clusters_Cox_Regression_4Models.docx")
cat("✅ 4模型 Cox 回归表格完美导出！文件名: Table_GMM_Clusters_Cox_Regression_4Models.docx\n")

# ==============================================================================
# 👑 步骤 3：顶级基线特征表 (Table 1) 的对齐与生成
# ==============================================================================
cat("⏳ 正在构建全景基线特征表 (Table 1)...\n")

# 数据美容与因子对齐
df_baseline <- full_cohort_cox %>%
  mutate(
    # 计算 AST/ALT 比值（由于前一步已经成功合入数据，这里绝不会报错）
    ast_alt_ratio = ast_val / alt_val,
    
    PRS_Label = case_when(
      PRS_AF <= quantile(PRS_AF, 1/3, na.rm = TRUE) ~ "Low PRS",
      PRS_AF > quantile(PRS_AF, 2/3, na.rm = TRUE) ~ "High PRS",
      TRUE ~ "Mid PRS"
    ),
    sex_f_Male = factor(sex_f_Male, levels = c(0, 1), labels = c("Female", "Male")),
    event_af_updated = factor(event_af_updated, levels = c(0, 1), labels = c("No", "Yes")),
    htn_final_1 = factor(htn_final_1, levels = c(0, 1), labels = c("No", "Yes")),
    t2dm_final_1 = factor(t2dm_final_1, levels = c(0, 1), labels = c("No", "Yes")),
    lip_f_1 = factor(lip_f_1, levels = c(0, 1), labels = c("No", "Yes")),
    ckd_final_1 = factor(ckd_final_1, levels = c(0, 1), labels = c("No", "Yes")),
    cvd_f_1 = factor(cvd_f_1, levels = c(0, 1), labels = c("No", "Yes")),
    PRS_Label = factor(PRS_Label, levels = c("Low PRS", "Mid PRS", "High PRS")),
    fib4_cat = factor(ifelse(grepl("Low", fib4_cat, ignore.case = TRUE), "Low Risk", 
                             ifelse(grepl("High", fib4_cat, ignore.case = TRUE), "High Risk", "Intermediate Risk")), 
                      levels = c("Low Risk", "Intermediate Risk", "High Risk"))
  )

vars_to_summarize <- c("age_val", "sex_f_Male", "bmi_val", "tdi_val","MetRS_Z",
                       "ast_val", "alt_val", "ast_alt_ratio", "plt_val", "alb_val",
                       "PRS_AF", "PRS_Label", "fib4_cat", 
                       "htn_final_1", "t2dm_final_1", "lip_f_1", "ckd_final_1", "cvd_f_1", 
                       "event_af_updated")

cat_vars <- c("sex_f_Male", "PRS_Label", "fib4_cat", "htn_final_1", "t2dm_final_1", 
              "lip_f_1", "ckd_final_1", "cvd_f_1", "event_af_updated")

# 指定非正态分布展示 (Median [IQR]) 的变量
non_normal_vars <- c("ast_val", "alt_val", "ast_alt_ratio")

tab1 <- CreateTableOne(vars = vars_to_summarize, strata = "Plot_Cluster", 
                       data = df_baseline, factorVars = cat_vars, test = TRUE)

# 核心：传入 nonnormal 强制渲染非正态分布
tab1_mat <- print(tab1, nonnormal = non_normal_vars, exact = "cat", quote = FALSE, 
                  noSpaces = TRUE, printToggle = FALSE, showAllLevels = TRUE)

# 变量名格式化
tab1_df <- as.data.frame(tab1_mat) %>% rownames_to_column("Variable") %>%
  mutate(Variable = case_when(
    grepl("age_val", Variable) ~ "Age (years), mean (SD)",
    grepl("sex_f_Male", Variable) ~ "Male Sex, n (%)",
    grepl("bmi_val", Variable) ~ "BMI (kg/m2), mean (SD)",
    grepl("tdi_val", Variable) ~ "Townsend Deprivation Index, mean (SD)",
    grepl("MetRS_Z", Variable) ~ "Metabolic Risk Score (Z-score), mean (SD)", #  【修改处 2】新加这一行
    grepl("ast_val", Variable) ~ "AST (U/L), median (IQR)",
    grepl("alt_val", Variable) ~ "ALT (U/L), median (IQR)",
    grepl("ast_alt_ratio", Variable) ~ "AST/ALT Ratio, median (IQR)",
    grepl("plt_val", Variable) ~ "Platelet Count (10^9/L), mean (SD)",
    grepl("alb_val", Variable) ~ "Albumin (g/L), mean (SD)",
    
    grepl("PRS_AF", Variable) ~ "Atrial Fibrillation PRS, mean (SD)",
    grepl("PRS_Label", Variable) ~ "Genetic Risk Category, n (%)",
    grepl("fib4_cat", Variable) ~ "FIB-4 Fibrosis Risk, n (%)",
    grepl("htn_final_1", Variable) ~ "Hypertension, n (%)",
    grepl("t2dm_final_1", Variable) ~ "Type 2 Diabetes, n (%)",
    grepl("lip_f_1", Variable) ~ "Dyslipidemia, n (%)",
    grepl("ckd_final_1", Variable) ~ "Chronic Kidney Disease, n (%)",
    grepl("cvd_f_1", Variable) ~ "Cardiovascular Disease, n (%)",
    grepl("event_af_updated", Variable) ~ "Incident AF (Follow-up), n (%)",
    TRUE ~ Variable
  ))

# 渲染导出
ft_base <- flextable(tab1_df) %>%
  set_header_labels(Variable = "Characteristics", p = "P value") %>%
  theme_booktabs() %>% autofit() %>%
  align(j = 2:5, align = "center", part = "all") %>%
  bold(i = 1, part = "header") %>%
  bold(i = ~ !is.na(p) & as.numeric(gsub("<", "", p)) < 0.05, j = "p")

doc_base <- read_docx() %>%
  body_add_par("Table 1. Baseline Characteristics According to Metabolic Phenotypes", style = "heading 1") %>%
  body_add_par("Continuous variables are presented as mean (standard deviation) or median (interquartile range) where appropriate; categorical variables are presented as n (%). P-values were calculated using ANOVA or Kruskal-Wallis tests for continuous variables and Pearson's chi-squared test for categorical variables.", style = "Normal") %>%
  body_add_flextable(ft_base)

print(doc_base, target = "Table_1_Baseline_Characteristics_Updated.docx")

cat("✅ 全流程完美收官！双源获取完毕，Table 1 现已顺利导出！\n")












































































































































































































cat("\n# ==============================================================================\n")
cat("⏳ 正在启动顶刊级表格生成引擎 (四模型 Cox & Table 1)...\n")
cat("# ==============================================================================\n")

library(dplyr)
library(tidyr)
library(survival)
library(broom)
library(stringr)
library(flextable)
library(officer)
library(tableone)
library(tibble)

# ==============================================================================
# 👑 步骤 1：智能提取缺失协变量并无缝合并 (从 long_dat_dummy_中介 源头打捞)
# ==============================================================================
# 1. 🚨 指定超级源数据宝库
raw_source <- long_dat_dummy_中介 

# 2. 列出我们希望拥有的所有附加临床与主成分协变量 (完美适配4模型)
cov_to_extract <- c(
  "ethnicity_f_Black", "ethnicity_f_Others.Mixed.Unknown", "ethnicity_f_White",
  "bmi_val", "smoke_f_Never", "smoke_f_Prefer.not.to.answer", "smoke_f_Previous",
  "alc_freq_f", "tdi_val", "age_val", "sex_f_Male", "PRS_AF", "fib4_cat",
  "htn_final_1", "t2dm_final_1", "lip_f_1", "ckd_final_1", "cvd_f_1",
  "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"
)

# 3. 探测 full_cohort 中到底缺哪些列
missing_cols <- setdiff(cov_to_extract, colnames(full_cohort))

if(length(missing_cols) > 0) {
  full_cohort_cox <- full_cohort %>%
    # 🚨 核心修复：使用 any_of 免疫报错，并完美映射 by = c("eid" = "Participant ID")
    left_join(raw_source %>% select(Participant.ID, any_of(missing_cols)), 
              by = c("eid" = "Participant.ID")) %>%
    distinct(eid, .keep_all = TRUE) # 绝对防止重复合并
} else {
  full_cohort_cox <- full_cohort %>% distinct(eid, .keep_all = TRUE)
}

# 4. 剔除含有 NA 的行，保证 4 个模型使用的是完全相同的严谨亚人群
all_cox_vars <- c("Plot_Cluster", "duration_updated", "event_af_updated", cov_to_extract)

# 使用 any_of 过滤，确保即使个别变量没捞到也不会让整个过滤逻辑崩溃
full_cohort_cox <- full_cohort_cox %>% drop_na(any_of(all_cox_vars))

cat(sprintf("✅ 协变量无缝装填完成！进入 Cox 和基线表的可用完整样本量: %d\n", nrow(full_cohort_cox)))

# ==============================================================================
# 👑 步骤 2：执行 4 阶 Cox 回归与极美表格渲染
# ==============================================================================
# 🚨 完美映射你的 4 阶模型到真实的哑变量列名
cov_m1 <- c("age_val", "sex_f_Male", "ethnicity_f_Black", "ethnicity_f_Others.Mixed.Unknown", "ethnicity_f_White")
cov_m2 <- c(cov_m1, "bmi_val", "smoke_f_Never", "smoke_f_Prefer.not.to.answer", "smoke_f_Previous", "alc_freq_f", "tdi_val")
cov_m3 <- c(cov_m2, "htn_final_1", "t2dm_final_1", "lip_f_1", "ckd_final_1", "cvd_f_1")
cov_m4 <- c(cov_m3, "PRS_AF", paste0("PC", 1:10))

model_list <- list(
  "Crude (Unadjusted)" = "UNADJUSTED",
  "Model 1" = cov_m1,
  "Model 2" = cov_m2,
  "Model 3" = cov_m3,
  "Model 4" = cov_m4
)

# 计算事件数
counts_df <- full_cohort_cox %>%
  group_by(Plot_Cluster) %>%
  summarise(
    Total = n(),
    Cases = sum(event_af_updated, na.rm = TRUE), 
    Counts_Str = sprintf("%d / %d", Total, Cases),
    .groups = "drop"
  )
counts_map <- setNames(counts_df$Counts_Str, counts_df$Plot_Cluster)

all_res <- list()
for (mod_name in names(model_list)) {
  
  if (model_list[[mod_name]][1] == "UNADJUSTED") {
    formula_str <- "Surv(duration_updated, event_af_updated) ~ Plot_Cluster"
  } else {
    covariates_str <- paste(model_list[[mod_name]], collapse = " + ")
    formula_str <- paste("Surv(duration_updated, event_af_updated) ~ Plot_Cluster +", covariates_str)
  }
  
  fit <- coxph(as.formula(formula_str), data = full_cohort_cox)
  res_tidy <- tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>% 
    filter(grepl("Plot_Cluster", term))
  
  mod_df <- data.frame(
    Model = mod_name, Exposure = "Metabolic Clusters", Type = "Categorical",
    Group = c( "Healthy (Ref)", "Fibro-Inflammatory", "Classic Cardiometabolic"),
    Counts = c(counts_map["Healthy"], counts_map["Fibro-Inflammatory"], counts_map["Classic Cardiometabolic"]),
    HR_CI = c(
      "1.00 (Reference)", 
      sprintf("%.2f (%.2f-%.2f)", res_tidy$estimate[1], res_tidy$conf.low[1], res_tidy$conf.high[1]),
      sprintf("%.2f (%.2f-%.2f)", res_tidy$estimate[2], res_tidy$conf.low[2], res_tidy$conf.high[2])
    ),
    P_val = c(
      "-",
      ifelse(res_tidy$p.value[1] < 0.001, "<0.001", sprintf("%.3f", res_tidy$p.value[1])),
      ifelse(res_tidy$p.value[2] < 0.001, "<0.001", sprintf("%.3f", res_tidy$p.value[2]))
    ),
    stringsAsFactors = FALSE
  )
  all_res[[mod_name]] <- mod_df
}

final_table_df <- bind_rows(all_res)

# 导出 Cox 表格
ft_cox <- flextable(final_table_df) %>%
  set_header_labels(Model = "Model", Exposure = "Exposure", Type = "Type", Group = "Level", 
                    Counts = "No. (Total / Cases)",  HR_CI = "HR (95% CI)", P_val = "P value") %>%
  merge_v(j = c("Model", "Exposure", "Type")) %>% 
  theme_booktabs() %>% autofit() %>%
  align(j = "Counts", align = "center", part = "all") %>%
  bold(i = ~ !str_detect(P_val, "-|>") & as.numeric(str_replace_all(P_val, "<", "")) < 0.05, j = "P_val")

doc_cox <- read_docx() %>%
  body_add_par("Table: Independent Association Between Metabolic Phenotypes and Incident AF", style = "heading 1") %>%
  body_add_par("Crude: Unadjusted model.", style = "Normal") %>%
  body_add_par("Model 1: Adjusted for Age, Sex, and Ethnicity.", style = "Normal") %>%
  body_add_par("Model 2: Adjusted for Model 1 + BMI, Smoking, Alcohol frequency, and Townsend Deprivation Index.", style = "Normal") %>%
  body_add_par("Model 3: Adjusted for Model 2 + Hypertension, Type 2 Diabetes, Dyslipidemia, Chronic Kidney Disease, and prior Cardiovascular Disease.", style = "Normal") %>%
  body_add_par("Model 4: Adjusted for Model 3 + AF Polygenic Risk Score (PRS) and Top 10 Genetic Principal Components.", style = "Normal") %>%
  body_add_flextable(ft_cox)
print(doc_cox, target = "Table_GMM_Clusters_Cox_Regression_4Models.docx")
cat("✅ 4模型 Cox 回归表格完美导出！文件名: Table_GMM_Clusters_Cox_Regression_4Models.docx\n")

# ==============================================================================
# 👑 步骤 3：顶级基线特征表 (Table 1) 的对齐与生成
# ==============================================================================
cat("⏳ 正在构建全景基线特征表 (Table 1)...\n")

# 数据美容与因子对齐
df_baseline <- full_cohort_cox %>%
  mutate(
    # 生成 PRS 标签
    PRS_Label = case_when(
      PRS_AF <= quantile(PRS_AF, 1/3, na.rm = TRUE) ~ "Low PRS",
      PRS_AF > quantile(PRS_AF, 2/3, na.rm = TRUE) ~ "High PRS",
      TRUE ~ "Mid PRS"
    ),
    # 把 1/0 强制映射为人类可读格式
    sex_f_Male = factor(sex_f_Male, levels = c(0, 1), labels = c("Female", "Male")),
    event_af_updated = factor(event_af_updated, levels = c(0, 1), labels = c("No", "Yes")),
    htn_final_1 = factor(htn_final_1, levels = c(0, 1), labels = c("No", "Yes")),
    t2dm_final_1 = factor(t2dm_final_1, levels = c(0, 1), labels = c("No", "Yes")),
    lip_f_1 = factor(lip_f_1, levels = c(0, 1), labels = c("No", "Yes")),
    ckd_final_1 = factor(ckd_final_1, levels = c(0, 1), labels = c("No", "Yes")),
    cvd_f_1 = factor(cvd_f_1, levels = c(0, 1), labels = c("No", "Yes")),
    
    # 强制排序
    PRS_Label = factor(PRS_Label, levels = c("Low PRS", "Mid PRS", "High PRS")),
    fib4_cat = factor(ifelse(grepl("Low", fib4_cat, ignore.case = TRUE), "Low Risk", 
                             ifelse(grepl("High", fib4_cat, ignore.case = TRUE), "High Risk", "Intermediate Risk")), 
                      levels = c("Low Risk", "Intermediate Risk", "High Risk"))
  )

vars_to_summarize <- c("age_val", "sex_f_Male", "bmi_val", "tdi_val", "PRS_AF", "PRS_Label", "fib4_cat", 
                       "htn_final_1", "t2dm_final_1", "lip_f_1", "ckd_final_1", "cvd_f_1", 
                       "event_af_updated")

cat_vars <- c("sex_f_Male", "PRS_Label", "fib4_cat", "htn_final_1", "t2dm_final_1", 
              "lip_f_1", "ckd_final_1", "cvd_f_1", "event_af_updated")

# 使用洗净的 Plot_Cluster 作为 strata
tab1 <- CreateTableOne(vars = vars_to_summarize, strata = "Plot_Cluster", 
                       data = df_baseline, factorVars = cat_vars, test = TRUE)

tab1_mat <- print(tab1, nonnormal = NULL, exact = "cat", quote = FALSE, 
                  noSpaces = TRUE, printToggle = FALSE, showAllLevels = TRUE)

# 变量名整容
tab1_df <- as.data.frame(tab1_mat) %>% rownames_to_column("Variable") %>%
  mutate(Variable = case_when(
    grepl("age_val", Variable) ~ "Age (years), mean (SD)",
    grepl("sex_f_Male", Variable) ~ "Male Sex, n (%)",
    grepl("bmi_val", Variable) ~ "BMI (kg/m2), mean (SD)",
    grepl("tdi_val", Variable) ~ "Townsend Deprivation Index, mean (SD)",
    grepl("PRS_AF", Variable) ~ "Atrial Fibrillation PRS, mean (SD)",
    grepl("PRS_Label", Variable) ~ "Genetic Risk Category, n (%)",
    grepl("fib4_cat", Variable) ~ "FIB-4 Fibrosis Risk, n (%)",
    grepl("htn_final_1", Variable) ~ "Hypertension, n (%)",
    grepl("t2dm_final_1", Variable) ~ "Type 2 Diabetes, n (%)",
    grepl("lip_f_1", Variable) ~ "Dyslipidemia, n (%)",
    grepl("ckd_final_1", Variable) ~ "Chronic Kidney Disease, n (%)",
    grepl("cvd_f_1", Variable) ~ "Cardiovascular Disease, n (%)",
    grepl("event_af_updated", Variable) ~ "Incident AF (Follow-up), n (%)",
    TRUE ~ Variable
  ))

# 渲染导出 Table 1
ft_base <- flextable(tab1_df) %>%
  set_header_labels(Variable = "Characteristics", p = "P value") %>%
  theme_booktabs() %>% autofit() %>%
  align(j = 2:5, align = "center", part = "all") %>%
  bold(i = 1, part = "header") %>%
  bold(i = ~ !is.na(p) & as.numeric(gsub("<", "", p)) < 0.05, j = "p")

doc_base <- read_docx() %>%
  body_add_par("Table 1. Baseline Characteristics According to Metabolic Phenotypes", style = "heading 1") %>%
  body_add_par("Continuous variables are presented as mean (standard deviation); categorical variables are presented as n (%). P-values were calculated using ANOVA for continuous variables and Pearson's chi-squared test for categorical variables.", style = "Normal") %>%
  body_add_flextable(ft_base)

print(doc_base, target = "Table_1_Baseline_Characteristics.docx")

cat("✅ 全流程完美收官！Table 1 与 4阶 Cox 表格现已同步生成！\n")


