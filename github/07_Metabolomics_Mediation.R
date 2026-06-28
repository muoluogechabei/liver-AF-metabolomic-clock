library(readxl)
library(dplyr)
library(tidyverse)
library(data.table)
library(lubridate)
library(stringr)
library(tidyr)
library(mice)
library(parallel)

#mice_rediag
imp_ready_for_cox<-mice_rediag
library(readxl)
library(dplyr)
library(tidyverse)
library(data.table)
library(lubridate)
library(stringr)
library(tidyr)
library(mice)
library(parallel)
library(survival)

# ==============================================================================
# 1. 读取新数据与对照表，进行 ID 转换与严格过滤
# ==============================================================================
cat("📖 1. 正在读取对照表、nmrnew 数据与字典...\n")
mapping_table <- fread("Cross_xueyu.csv")
nmr_new_raw <- fread("nmrnew.csv", na.strings = c("", "NA"))
dict <- read_excel("dict.xlsx")

cat("🔗 2. 正在通过对照表进行跨库 ID 转换...\n")
# 将 nmrnew 的 eid 与 mapping_table 的 ID_in_ckd 匹配，并将主键替换为 ID_in_participant
nmr_mapped <- nmr_new_raw %>%
  inner_join(mapping_table, by = c("eid" = "ID_in_ckd")) %>%
  mutate(eid = ID_in_participant) %>% # 重命名为 eid 确保完全兼容你后续的原始代码
  select(-ID_in_participant) 

cat("🧹 3. 正在剔除后缀为 'i1' 的列，并排除任何含有缺失值的样本...\n")
# 丢掉 i1 后缀列，随后剩下的列只要有任何 NA 就直接删掉这一行
my_data <- nmr_mapped %>%
  select(-ends_with("i1")) %>%
  drop_na()

cat("过滤后，输入管道的主数据样本量为:", nrow(my_data), "\n\n")

# ==============================================================================
# 2. 原封不动的核心数据清洗与转换代码 (完全未改动)
# ==============================================================================
# 3. 获取当前主数据的所有列名
current_colnames <- colnames(my_data)

# 4. 在字典中寻找匹配项
# match() 会返回 current_colnames 在字典 PIA_code 列中的行号，找不到的返回 NA
match_idx <- match(current_colnames, dict$id)

# 5. 安全替换列名：如果匹配到了就用 Full_name，没匹配到（是 NA）就保留原本的名字
new_colnames <- ifelse(is.na(match_idx), 
                       current_colnames, 
                       dict$title[match_idx])

clean_colnames <- make.names(new_colnames, unique = TRUE)

# 7. 将新列名赋给主数据
colnames(my_data) <- clean_colnames

# 8. 检查替换结果
head(colnames(my_data), 20)

metabolite_anno <- data.frame(
  raw_id = current_colnames,    # 关键：用最原始的列名（比如 "23400" 等数字ID）作为匹配钥匙
  clean_name = clean_colnames,  # 带着你清洗后的小数点名字
  stringsAsFactors = FALSE
) %>%
  # 使用字典里的原始 id 列（数字列）进行匹配
  inner_join(dict %>% mutate(id = as.character(id)), 
             by = c("raw_id" = "id")) %>%
  select(raw_id, clean_name, Group, Subgroup, title) %>%
  # --- 关键清洗步骤 ---
  mutate(
    # 1. 清理潜在的空白字符
    Subgroup = trimws(Subgroup),
    
    # 2. 如果 Subgroup 是 NA 或者空字符串，就用 Group 的名字填补它
    Subgroup = ifelse(is.na(Subgroup) | Subgroup == "", Group, Subgroup),
    
    # 3. 组合层级因子
    Plot_Order = paste(Group, Subgroup, sep = " - ")
  )

# 检查一下总行数，如果显示是 251（或者你之前匹配成功的数字），就说明 100% 成功了！
cat("目前包含的代谢物总数:", nrow(metabolite_anno), "\n")

# 再看一眼前 15 行，这次 VLDL 和 Cholesterol 肯定都在了
head(metabolite_anno %>% select(Group, Subgroup, Plot_Order), 15)

# 获取 251 个代谢物列名
all_metab_cols <- metabolite_anno$clean_name

# 1. 计算每个参与者（每行）在这 251 个代谢物里的 NA 数量
# (由于上方已经 drop_na()，这里的 na_count 理论上都会是 0，完美衔接旧逻辑)
my_data$na_count <- rowSums(is.na(my_data %>% select(all_of(all_metab_cols))))

# 2. 建立固定的 NMR 子队列：剔除缺失超过 10% (即 > 25 个) 的人
my_data_nmr <- my_data %>% 
  filter(na_count ==0) 

cat("过滤后，NMR 子队列最终纳入人数:", nrow(my_data_nmr), "\n")


# 4. 再次检查是否还有 NA
total_na_left <- sum(is.na(my_data_nmr %>% select(all_of(all_metab_cols))))
cat("填补后，251个代谢物中剩余 NA 数量:", total_na_left, "(必须为 0)\n")

# 5. 最后进行咱们说过的对数转换和标准化（基于这个完美的数据集）
my_data_nmr <- my_data_nmr %>%
  mutate(across(all_of(all_metab_cols), ~ {
    x_log <- log(.x + 1)
    as.numeric(scale(x_log))
  }))

cat("✅ 缺失值处理与标准化全部完成，数据准备完美！\n\n")

# ==============================================================================
# 3. 合并临床插补数据并执行 QC 核对
# ==============================================================================

# 1. 提取插补库中的第 1 个完整数据集
clinical_data_1 <- complete(imp_ready_for_cox, 1)

# 2. 将临床数据与转换好 ID 的 NMR 数据合并
merged_data <- my_data_nmr %>%
  inner_join(clinical_data_1, by = c("eid" = "Participant ID"))

cat("✅ 合并完成！最终入组并合并了临床参数的样本量为:", nrow(merged_data), "\n\n")

# 3. 运行 QC 检查，验证 ID mapping 是否精准无误
cat("🔍 开始执行 QC 校验：对比主数据 Sex 与 nmrnew 的 p31 列...\n")
if("sex_f" %in% colnames(merged_data) && "p31" %in% colnames(merged_data)) {
  # 创建列联表进行交叉比对
  qc_table <- table(Clinical_Sex = merged_data$sex_f, NMR_p31 = merged_data$p31, useNA = "ifany")
  print(qc_table)
  cat("\n💡 提示：如果矩阵的数据完全集中在对应的两个对角线位置（比如 Male对应1，Female对应0），则说明跨库嫁接的 ID 转换是完美的！\n")
} else {
  cat("⚠️ 警告：无法在合并后的数据中同时找到 'Sex' 和 'p31'。请检查它们是否在 dict 字典中被重新命名了。\n")
  cat("可能的候选列名：\n")
  print(grep("(?i)sex|p31", colnames(merged_data), value = TRUE))
}


library(dplyr)
library(survival)
imp_ready_for_cox<-mice_rediag
# 1. 提取插补库中的第 1 个完整数据集
clinical_data_1 <- complete(imp_ready_for_cox, 1)

# 2. 将临床数据与 NMR 代谢物数据合并
merged_data <- my_data_nmr %>%
  inner_join(clinical_data_1, by = c("eid" = "Participant ID"))

# ==============================================================================
# 新增：暴露指标 (肝纤维化评分) 的标准化处理
# ==============================================================================
cat("正在对 4 个肝脏暴露指标进行缩尾和标准化...\n")

# 定义缩尾函数 (1%-99%)，防止极端值拉偏标准差
winsorize_it <- function(x) {
  limits <- quantile(x, probs = c(0.01, 0.99), na.rm = TRUE)
  x[x < limits[1]] <- limits[1]
  x[x > limits[2]] <- limits[2]
  return(x)
}

# 对 4 个核心暴露指标进行处理
# 直接覆盖原变量名，这样你后面的循环代码不需要改变量名就能直接运行
merged_data <- merged_data %>%
  mutate(across(all_of(c("fib4", "nfs", "apri", "ast_alt_ratio")), 
                ~ as.numeric(scale(winsorize_it(.x)))))

cat("✅ 暴露指标标准化完成。现在所有肝脏评分单位均为 '每标准差 (SD)'。\n")
cat("合并后的总样本量:", nrow(merged_data), "\n")
cat("合并后的总列数:", ncol(merged_data), "\n")
# ==============================================================================
# 🚨 新增核心：在所有分析开始前进行 7:3 物理盲切
# ==============================================================================
set.seed(2026) 
train_size <- floor(0.7 * nrow(merged_data))
train_indices <- sample(seq_len(nrow(merged_data)), size = train_size)

# 只把这 70% 的人拿去做后续的 Cox 和中介粗筛
merged_data_train <- merged_data[train_indices, ]

# 🔑 极度关键：把这批人的 ID 锁进保险箱，保证后面的 LASSO 切分是同一拨人
train_eids_locked <- merged_data_train$eid 

cat("✅ 物理隔离完成！用于粗筛的训练集人数:", nrow(merged_data_train), "\n")
# ==============================================================================
# --- 提取所有代谢物列名 ---
# 假设 metabolite_anno 还在你的环境中
all_metab_cols <- metabolite_anno$clean_name

# --- 构造极其严谨的协变量字符串 ---
# 包含了你辛苦插补和重新诊断的所有变量
# 2. 自动抓取 10 个 PC 的列名（带有 "Genetic principal components | Array" 的列）
pc_cols <- grep("Genetic principal components \\| Array", colnames(merged_data), value = TRUE)
# 给这 10 个 PC 列名加上反引号，并用 + 号连接
pc_str <- paste0("`", pc_cols, "`", collapse = " + ")

# 3. 准备 PRS 的列名（加上反引号）
prs_str <- "`Standard PRS for atrial fibrillation (AF)`"

# 4. 拼装终极版的协变量字符串 (临床 + 遗传)
covariates_str <- paste(
  "age_val + sex_f + ethnicity_f + bmi_val + waist_val + smoke_f + alc_status_f + tdi_val + htn_final + t2dm_final + ckd_final + cvd_f + lip_f",
  prs_str,
  pc_str,
  sep = " + "
)
# --- 准备存储结果的空数据框 ---
cox_results_list <- list()

# --- 开启循环 ---
cat("开始运行 251 个代谢物的 Cox 模型...\n")
for (i in seq_along(all_metab_cols)) {
  
  metab_var <- all_metab_cols[i]
  
  # 构造动态 Cox 公式
  # Surv(随访时间, 房颤事件) ~ 代谢物 + 协变量
  formula_str <- paste0("Surv(duration_updated, event_af_updated) ~ ", metab_var, " + ", covariates_str)
  
  # 拟合 Cox 模型
  fit <- tryCatch({
    coxph(as.formula(formula_str), data = merged_data_train)
  }, error = function(e) NULL)
  
  if (!is.null(fit)) {
    sum_fit <- summary(fit)$coefficients
    conf_fit <- summary(fit)$conf.int
    
    # 提取第一行（即我们关注的那个代谢物 metab_var 的结果）
    cox_results_list[[i]] <- data.frame(
      clean_name = metab_var,
      HR = conf_fit[1, 1],
      lower_95 = conf_fit[1, 3],
      upper_95 = conf_fit[1, 4],
      P_value = sum_fit[1, 5],
      stringsAsFactors = FALSE
    )
  }
  
  if(i %% 50 == 0) cat("已完成:", i, "/ 251\n")
}

# --- 合并结果并做 FDR 校正 ---
all_cox_results <- bind_rows(cox_results_list) %>%
  mutate(P_fdr = p.adjust(P_value, method = "bonferroni")) %>%
  left_join(metabolite_anno, by = "clean_name") # 关联回 Group 和 Subgroup 分类

# 查看筛选出的显著代谢物
significant_metabs <- all_cox_results %>% filter(P_fdr < 0.05)
cat("\nFDR 校正后，与房颤显著相关的代谢物有:", nrow(significant_metabs), "个\n")

# 按照 HR 降序看看前几名是谁
head(significant_metabs %>% arrange(desc(HR)) %>% select(clean_name, HR, P_fdr, Group), 20)
#saveRDS(all_cox_results, "all_cox_results.rds")










library(dplyr)
# 1. 准备4个肝纤维化暴露指标
exposures <- c("fib4", "nfs", "apri", "ast_alt_ratio")
mediation_results_list <- list()

cat("🔥 开始计算 4个肝脏指标 -> 251个代谢物 的线性回归与粗算中介比例...\n")

row_idx <- 1
for (exp_var in exposures) {
  for (metab_var in all_metab_cols) {
    
    # 提取之前算好的 Path B (代谢物 -> 房颤) 的结果
    path_b_res <- all_cox_results %>% filter(clean_name == metab_var)
    if(nrow(path_b_res) == 0) next
    
    beta_b <- log(path_b_res$HR)
    p_b_raw <- path_b_res$P_value
    p_b_fdr <- path_b_res$P_fdr  # 🚨 关键：直接提取之前算好的 Path B FDR P值
    
    # 构造 Path A 的线性回归公式: 代谢物 ~ 肝指标 + 协变量
    # 注意：这里我们使用同一套协变量 (不加PRS和PC，因为那是心血管专用的，影响肝脏的主要是临床指标)
    formula_a <- paste0(metab_var, " ~ ", exp_var, " + age_val + sex_f + bmi_val + waist_val + smoke_f + alc_status_f + tdi_val")
    
    # 跑线性回归
    fit_a <- tryCatch(lm(as.formula(formula_a), data = merged_data_train), error=function(e) NULL)
    
    if(!is.null(fit_a)){
      sum_a <- summary(fit_a)$coefficients
      # 提取肝脏指标的系数
      beta_a <- sum_a[2, 1] 
      p_a_raw <- sum_a[2, 4]  # 🚨 提取 Path A 的原始 P 值
      
      # 粗算中介效应 (Indirect Effect = beta_a * beta_b)
      indirect_effect <- beta_a * beta_b
      
      # 记录结果 (暂不计算 Combined_P，留到循环外统一处理)
      mediation_results_list[[row_idx]] <- data.frame(
        Liver_Index = exp_var,
        clean_name = metab_var,
        Beta_A = beta_a,
        P_A_raw = p_a_raw,      # 保存原始 Path A P值
        Beta_B = beta_b,
        P_B_raw = p_b_raw,      # 保存原始 Path B P值 (备查)
        P_B_fdr = p_b_fdr,      # 保存 FDR 校正后的 Path B P值
        Indirect_Effect = indirect_effect
      )
      row_idx <- row_idx + 1
    }
  }
  cat("已完成暴露指标:", exp_var, "\n")
}

# ==============================================================================
# 🚨 核心修改：循环结束后，统一进行 Path A 的 FDR 校正并计算终极 Combined_P
# ==============================================================================
cat("🔥 正在进行 Path A 的全局 FDR 校正并计算 Combined_P...\n")

# 合并所有初步结果
all_mediation_data_raw <- bind_rows(mediation_results_list)

all_mediation_data <- all_mediation_data_raw %>%
  # 1. 对 Path A 的原始 P 值进行统一 FDR 校正
  mutate(P_A_fdr = p.adjust(P_A_raw, method = "bonferroni")) %>%
  
  # 2. 计算极度严谨的 Combined_P = P_A_fdr * P_B_fdr
  mutate(Combined_P = P_A_fdr * P_B_fdr) %>%
  
  # 3. 关联分类信息
  left_join(metabolite_anno, by = "clean_name") %>%
  
  # 4. 严厉的双向过滤：Path A 和 Path B 都要满足 FDR < 0.05
  filter(P_A_fdr < 0.05 & P_B_fdr < 0.05) %>%
  arrange(Group, Subgroup)

# 给每个代谢物一个固定的圆环X轴ID (1到N)
metab_levels <- unique(all_mediation_data$clean_name)
all_mediation_data <- all_mediation_data %>%
  mutate(x_id = as.numeric(factor(clean_name, levels = metab_levels)))

cat("✅ 双向 FDR 过滤完成！剩余具有高度显著性的组合数量为:", nrow(all_mediation_data), "\n")

library(tidyverse)
library(CMAverse)
library(patchwork)

library(fastDummies)

# ==============================================================================
# 1. 从 mice 对象提取 10 个插补集的临床数据 (长格式)
# ==============================================================================
message("正在提取 10 个插补集的临床数据...")
# mice_rediag 是你之前保存的 mids 对象
long_clinical <- complete(imp_ready_for_cox, action = "long", include = FALSE)

# ==============================================================================
# 2. 将标准化后的 NMR 代谢物数据合并进去
# ==============================================================================
# 注意：代谢物已经在 my_data_nmr 中做好了 Log 和 Z-score
# 我们只需要取出 eid 和代谢物列，合并到每一个插补集中
message("正在合并 NMR 代谢物数据...")
nmr_subset <- my_data_nmr %>% 
  select(eid, all_of(all_metab_cols))

long_dat_combined <- long_clinical %>%
  inner_join(nmr_subset, by = c("Participant ID" = "eid"))

# ==============================================================================
# 3. 暴露变量预处理 (缩尾 + 标准化)
# ==============================================================================
winsorize_it <- function(x) {
  limits <- quantile(x, probs = c(0.01, 0.99), na.rm = TRUE)
  x[x < limits[1]] <- limits[1]
  x[x > limits[2]] <- limits[2]
  return(x)
}

long_dat_processed <- long_dat_combined %>%
  mutate(
    fib4_sd = as.numeric(scale(winsorize_it(fib4))),
    nfs_sd  = as.numeric(scale(winsorize_it(nfs))),
    apri_sd = as.numeric(scale(winsorize_it(apri))),
    ratio_sd = as.numeric(scale(winsorize_it(ast_alt_ratio)))
  ) %>%
  # 统一 PRS 和 PC 的列名，防止空格报错
  rename(
    PRS_AF = `Standard PRS for atrial fibrillation (AF)`,
    PC1 = `Genetic principal components | Array 1`,
    PC2 = `Genetic principal components | Array 2`,
    PC3 = `Genetic principal components | Array 3`,
    PC4 = `Genetic principal components | Array 4`,
    PC5 = `Genetic principal components | Array 5`,
    PC6 = `Genetic principal components | Array 6`,
    PC7 = `Genetic principal components | Array 7`,
    PC8 = `Genetic principal components | Array 8`,
    PC9 = `Genetic principal components | Array 9`,
    PC10 = `Genetic principal components | Array 10`
  )

# ==============================================================================
# 4. 转换哑变量 (Dummy Variables)
# ==============================================================================
# CMAverse 处理数值型矩阵最稳定，需将因子变量转为 0/1
cat_vars <- c("sex_f", "ethnicity_f", "smoke_f", "alc_status_f", 
              "htn_final", "t2dm_final", "lip_f", "ckd_final", "cvd_f")

message("正在生成哑变量...")
long_dat_dummy <- dummy_cols(long_dat_processed, 
                             select_columns = cat_vars, 
                             remove_first_dummy = TRUE,
                             remove_selected_columns = TRUE)

# 规范化列名 (处理 dummy 产生的空格)
colnames(long_dat_dummy) <- make.names(colnames(long_dat_dummy))

# 动态获取最终协变量列表
final_dummy_cols <- colnames(long_dat_dummy)[grep(paste0("^(", paste(cat_vars, collapse="|"), ")_"), colnames(long_dat_dummy))]
covariates_final <- c("age_val", "bmi_val", "waist_val", "tdi_val", "PRS_AF", paste0("PC", 1:10), final_dummy_cols)
#saveRDS(long_dat_dummy, "long_dat_dummy_中介.rds")
#long_dat_dummy<-long_dat_dummy_中介
cat("✅ 数据预处理完成！long_dat_dummy 已就绪，样本量:", nrow(long_dat_dummy), "\n")
# ==============================================================================
# 升级版 CMAverse 函数 (自动兼容 n_imps = 1 的全景扫描 和 n_imps > 1 的池化)
# ==============================================================================
run_cmaverse_loop <- function(data, exposure, mediator, covariates, n_imps = 10) {
  res_list <- list()
  message(paste("\n⚡ 开始分析:", exposure, "->", mediator))
  
  for(i in 1:n_imps) {
    dat_i <- data %>% filter(.imp == i) %>% as.data.frame()
    analysis_cols <- c(exposure, mediator, "duration_updated", "event_af_updated", covariates)
    dat_i_clean <- dat_i %>% select(all_of(analysis_cols)) %>% drop_na()
    
    if(nrow(dat_i_clean) < 100) next
    
    valid_covs <- covariates[sapply(dat_i_clean[covariates], function(x) length(unique(x)) > 1)]
    
    tryCatch({
      m_mean <- mean(dat_i_clean[[mediator]], na.rm = TRUE)
      
      fit_cma <- cmest(data = dat_i_clean,
                       exposure = exposure, mediator = mediator,
                       outcome = "duration_updated", event = "event_af_updated",
                       basec = valid_covs, model = "rb", yreg = "coxph",
                       mreg = list("linear"), mval = list(m_mean),
                       estimation = "paramfunc", inference = "delta", EMint = FALSE)
      
      res_mat <- summary(fit_cma)$summary
      
      # 强力手动提取 Path a 和 Path b
      form_m <- as.formula(paste(mediator, "~", exposure, "+", paste(valid_covs, collapse = "+")))
      form_y <- as.formula(paste("Surv(duration_updated, event_af_updated) ~", exposure, "+", mediator, "+", paste(valid_covs, collapse = "+")))
      
      fit_m_manual <- lm(form_m, data = dat_i_clean)
      fit_y_manual <- coxph(form_y, data = dat_i_clean)
      
      m_sum <- summary(fit_m_manual)$coefficients
      y_sum <- summary(fit_y_manual)$coefficients
      
      res_list[[i]] <- data.frame(
        TE_est = res_mat["Rte", "Estimate"], TE_se = res_mat["Rte", "Std.error"],
        ACME_est = res_mat["Rtnie", "Estimate"], ACME_se = res_mat["Rtnie", "Std.error"],
        Prop_est = res_mat["pm", "Estimate"], Prop_se = res_mat["pm", "Std.error"],
        beta1_val = m_sum[exposure, "Estimate"], beta1_se = m_sum[exposure, "Std. Error"],
        beta2_val = y_sum[mediator, "coef"], beta2_se = y_sum[mediator, "se(coef)"],
        n_obs = nrow(dat_i_clean),
        n_events = sum(dat_i_clean$event_af_updated == 1)
      )
      
    }, error = function(e) message(paste("❌ 插补集", i, "报错:", e$message)))
  }
  
  if(length(res_list) == 0) return(NULL)
  
  res_df <- bind_rows(res_list) %>%
    mutate(
      log_TE = log(TE_est), log_TE_se = TE_se / TE_est,
      log_ACME = log(ACME_est), log_ACME_se = ACME_se / ACME_est,
      prop_val = Prop_est, prop_se = Prop_se
    )
  
  # =======================================================
  # 🚀 核心修复：特判 n_imps == 1 (免去方差池化，直接输出)
  # =======================================================
  if (n_imps == 1) {
    pooled_res <- res_df %>%
      summarise(
        M = 1, N_avg = n_obs[1], Events_avg = n_events[1],
        pool_log_TE = log_TE[1], se_log_TE = log_TE_se[1],
        pool_log_ACME = log_ACME[1], se_log_ACME = log_ACME_se[1],
        pool_prop = prop_val[1], se_prop = prop_se[1],
        Exposure = exposure, Mediator = mediator,
        pool_beta1 = beta1_val[1], se_beta1 = beta1_se[1],
        pool_beta2 = beta2_val[1], se_beta2 = beta2_se[1]
      )
  } else {
    # 原版的 Rubin's Rules 池化
    pooled_res <- res_df %>%
      summarise(
        M = n(), N_avg = round(mean(n_obs)), Events_avg = round(mean(n_events)),
        pool_log_TE = mean(log_TE), se_log_TE = sqrt(mean(log_TE_se^2) + (1 + 1/M) * var(log_TE)),
        pool_log_ACME = mean(log_ACME), se_log_ACME = sqrt(mean(log_ACME_se^2) + (1 + 1/M) * var(log_ACME)),
        pool_prop = mean(prop_val), se_prop = sqrt(mean(prop_se^2) + (1 + 1/M) * var(prop_val)),
        Exposure = exposure, Mediator = mediator,
        pool_beta1 = mean(beta1_val, na.rm = TRUE), se_beta1 = sqrt(mean(beta1_se^2, na.rm = TRUE) + (1 + 1/M) * var(beta1_val, na.rm = TRUE)),
        pool_beta2 = mean(beta2_val, na.rm = TRUE), se_beta2 = sqrt(mean(beta2_se^2, na.rm = TRUE) + (1 + 1/M) * var(beta2_val, na.rm = TRUE))
      )
  }
  
  # 统一计算 P 值和置信区间字符串 (单次或池化后都可以跑)
  pooled_res %>%
    mutate(
      beta1_p = 2 * (1 - pnorm(abs(pool_beta1 / se_beta1))),
      beta2_p = 2 * (1 - pnorm(abs(pool_beta2 / se_beta2))),
      
      Path_a_label = sprintf("%.3f (p=%s)", pool_beta1, ifelse(beta1_p < 0.001, "<0.001", sprintf("%.3f", beta1_p))),
      Path_b_HR = sprintf("%.3f (p=%s)", exp(pool_beta2), ifelse(beta2_p < 0.001, "<0.001", sprintf("%.3f", beta2_p))),
      
      TE_p = 2 * (1 - pnorm(abs(pool_log_TE / se_log_TE))),
      ACME_p = 2 * (1 - pnorm(abs(pool_log_ACME / se_log_ACME))),
      
      TE_HR = exp(pool_log_TE), TE_LB = exp(pool_log_TE - 1.96 * se_log_TE), TE_UB = exp(pool_log_TE + 1.96 * se_log_TE),
      ACME_HR = exp(pool_log_ACME), ACME_LB = exp(pool_log_ACME - 1.96 * se_log_ACME), ACME_UB = exp(pool_log_ACME + 1.96 * se_log_ACME),
      
      TE_CI = sprintf("%.3f (%.3f-%.3f)", TE_HR, TE_LB, TE_UB),
      ACME_CI = sprintf("%.3f (%.3f-%.3f)", ACME_HR, ACME_LB, ACME_UB),
      TE_p_val = ifelse(TE_p < 0.001, "<0.001", sprintf("%.3f", TE_p)),
      ACME_p_val = ifelse(ACME_p < 0.001, "<0.001", sprintf("%.3f", ACME_p)),
      Prop_formatted = sprintf("%.1f%% (%.1f%%, %.1f%%)", pool_prop*100, (pool_prop-1.96*se_prop)*100, (pool_prop+1.96*se_prop)*100)
    )
}
library(ggVennDiagram)
library(tidyverse)

# ==============================================================================
# 🚨 优化版策略：基于 MASLD 两大黄金指标 (FIB-4 & NFS) 的跨指标交集初筛
# ==============================================================================
cat("正在执行 FIB-4 与 NFS 的双重严苛门槛过滤...\n")

# 1. 严格过滤：锁定 FIB-4 和 NFS，应用双向显著与效应量门槛
strict_filtered_data <- all_mediation_data %>%
  # 踢掉不敏感的 APRI 和不稳定的 AST/ALT
  filter(Liver_Index %in% c("fib4", "nfs")) %>%
  
  # 门槛 1：双向 FDR < 0.05
  filter(P_A_fdr < 0.05 & P_B_fdr < 0.05) %>%
  
  # 门槛 2：效应量过滤 (Path A > 0.05 SD, Path B 风险波动 > 5%)
  filter(abs(Beta_A) > 0.05) %>%
  filter(exp(Beta_B) > 1.05 | exp(Beta_B) < 0.95)

# 2. 提取两个评分各自过线的代谢物名单
list_fib4 <- strict_filtered_data %>% filter(Liver_Index == "fib4") %>% pull(clean_name)
list_nfs  <- strict_filtered_data %>% filter(Liver_Index == "nfs") %>% pull(clean_name)

# 3. 🎨 绘制绝美的“双圈”韦恩图 (Venn Diagram)
venn_list_2 <- list(
  `FIB-4` = list_fib4,
  `NFS`   = list_nfs
)

p_venn_2 <- ggVennDiagram(venn_list_2, label_alpha = 0, set_color = "black", edge_size = 0.6) +
  scale_fill_gradient(low = "#F4FAF4", high = "#4A72B5", name = "Count") + # 优化了极简渐变色
  theme_void() +
  labs(title = "The Core Liver-Heart Metabolome",
       subtitle = "Intersection across FIB-4 and NFS\n(FDR < 0.05, |Beta_A| > 0.05, HR > 1.05 or < 0.95)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16, margin = margin(b = 10)),
        plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40", margin = margin(b = 20)),
        legend.position = "right")

# 导出这个清爽的双圈韦恩图
ggsave("Figure_Venn_Core_Metabolome_2way.pdf", plot = p_venn_2, width = 6, height = 5, dpi = 400)

# 4. 获取终极交集名单 (The Core Pool) -> 这 144 个就是你后续聚类和 LASSO 的宝库！
core_metabolites <- intersect(list_fib4, list_nfs)
cat("🔥 终极跨指标交集筛选完毕！FIB-4 与 NFS 共同的核心代谢物共有:", length(core_metabolites), "个\n")
# ==============================================================================
# 🚨 专供 CMAverse 中介画图的“代表抽样” (各大门派最强 Top 2)
# ==============================================================================
cat("正在从核心代谢池中，为因果中介分析提取各大门派的『最强代表』...\n")

top_mediators_final <- strict_filtered_data %>%
  # 1. 锁定那 144 个 Core Pool 成员
  filter(clean_name %in% core_metabolites) %>%
  
 
  
  # 4. 生成标签
  mutate(
    Exposure_SD_Name = case_when(
      Liver_Index == "fib4" ~ "fib4_sd",
      Liver_Index == "nfs"  ~ "nfs_sd"
    ),
    combo_id = paste(Exposure_SD_Name, clean_name, sep = "_")
  )

cat("✅ 代表提取完成！经过大类去重，精简为", nrow(top_mediators_final), "个最强组合。\n")

# ==============================================================================
# 2. 自动化 CMAverse 循环 (正式挂机版)
# ==============================================================================
final_cma_list <- list()
message("⚡ 正在启动 40个核心组合的 CMAverse 严谨分析 (10次插补池化)...")

for(i in 1:nrow(top_mediators_final)) {
  exp_sd <- top_mediators_final$Exposure_SD_Name[i] # 使用映射后的名字
  med    <- top_mediators_final$clean_name[i]
  
  # 这里建议直接跑 n_imps = 10，保证池化结果稳健
  res <- run_cmaverse_loop(long_dat_dummy, exp_sd, med, covariates_final, n_imps = 1)
  
  if(!is.null(res)) {
    res$Group <- top_mediators_final$Group[i]
    final_cma_list[[top_mediators_final$combo_id[i]]] <- res
  }
  cat("进度:", i, "/ 288 完成\n")
}

# ==============================================================================
# 3. 结果合并 (确保可视化标签正确)
# ==============================================================================
df_viz <- bind_rows(final_cma_list) %>%
  mutate(
    # 之前的 str_remove 改为处理映射后的名字
    Exposure_Tag = case_when(
      Exposure == "fib4_sd"  ~ "FIB-4",
      Exposure == "nfs_sd"   ~ "NFS"
      
    ),
    Plot_Label = paste0(Mediator, " (", Exposure_Tag, ")"),
    Plot_Label = reorder(Plot_Label, ACME_HR)
  )
#saveRDS(df_viz, "df_viz_4_15.rds")

library(patchwork)
library(dplyr)
library(ggnewscale) # 必须加载，用于实现双重 fill 映射
library(ggplot2)
# ==============================================================================
# 0. 定义高级配色方案
# ==============================================================================

# 2. 定义一套 Nature 风格的代谢物大类配色 (预留了足够多的颜色)
nature_pal <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", 
  "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", 
  "#7E6148", "#B09C85", "#D39200", "#00BA38"
)

# ==============================================================================
# 1. 整理绘图数据与背景数据
# ==============================================================================
# ==============================================================================
# 修正：包含 APRI 的完整绘图数据准备
# ==============================================================================
df_viz_plot <- df_viz %>%
  mutate(
    plot_p = ifelse(ACME_p == 0, 1e-300, ACME_p),
    neg_log10_p = -log10(plot_p),
    Clean_Label = Mediator,
    
    # 补全 Exposure_Tag，把 APRI 加进去
    Exposure_Tag = case_when(
      Exposure == "fib4_sd"  ~ "FIB-4",
      Exposure == "nfs_sd"   ~ "NFS",
      Exposure == "apri_sd"  ~ "APRI"       # 补上这一行
    )
  ) %>%
  # 设置分块的显示顺序
  mutate(Exposure_Tag = factor(Exposure_Tag, levels = c("FIB-4", "NFS", "APRI", "AST/ALT"))) %>%
  arrange(Exposure_Tag, ACME_HR) %>%
  mutate(Clean_Label = factor(Clean_Label, levels = unique(Clean_Label)))

# 更新背景色定义，增加 APRI 的颜色
facet_bg_colors <- c(
  "FIB-4" = "#F2F6FA", 
  "NFS"   = "#FFF7F0", 
  "APRI"  = "#F4FAF4" 
)

# 重新生成背景数据框
bg_facet <- data.frame(Exposure_Tag = factor(levels(df_viz_plot$Exposure_Tag), levels = levels(df_viz_plot$Exposure_Tag)))
# ==============================================================================
# 2. 图 A：分块气泡森林图 (带高级底色)
# ==============================================================================
# ==============================================================================
# 2. 图 A：分块气泡森林图 (找回 Y 轴行名版)
# ==============================================================================
p_forest <- ggplot(df_viz_plot, aes(y = Clean_Label)) +
  # --- 第一层：画区块底色 ---
  geom_rect(data = bg_facet, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = Exposure_Tag), inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_manual(values = facet_bg_colors, guide = "none") +
  new_scale_fill() + 
  
  # --- 第二层：画基准线与数据 ---
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ACME_LB, xmax = ACME_UB, x = ACME_HR, color = Group), height = 0.2, linewidth = 0.8) +
  geom_point(aes(x = ACME_HR, size = neg_log10_p, fill = Group), shape = 21, color = "white", stroke = 0.5) +
  
  scale_fill_manual(name = "Metabolic Super-class", values = nature_pal) +
  scale_color_manual(name = "Metabolic Super-class", values = nature_pal) +
  scale_size_continuous(name = "Significance\n(-log10 P-value)", range = c(2, 6)) +
  
  facet_grid(Exposure_Tag ~ ., scales = "free_y", space = "free_y", switch = "y") +
  
  theme_minimal() +
  labs(x = "Causal Mediation Effect (ACME HR)", y = NULL) + # 注意这里 X 轴应该是 HR
  theme(
    legend.position = "left", # 建议图例放左边或收集到右边
    axis.text.y = element_text(size = 9, face = "bold", color = "grey20"), # ✅ 核心修复：确保文字存在
    panel.grid.major.y = element_line(color = "white", linewidth = 0.6), 
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.background = element_rect(fill = "grey40", color = NA),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 11, color = "white"),
    panel.spacing = unit(0.3, "lines")
  ) +
  guides(
    color = "none", 
    # 代谢物分类图例：确保圆点大小固定且不透明
    fill = guide_legend(
      order = 1, 
      override.aes = list(size = 4, shape = 21, alpha = 1)
    ),
    # 显著性大小图例：✅ 核心修复在这里
    size = guide_legend(
      order = 2,
      override.aes = list(
        shape = 21,         # 强制使用形状 21
        fill = "grey70",    # 🚀 给图例里的泡泡一个可见的填充色
        color = "white",    # 泡泡边框颜色
        stroke = 0.5        # 边框粗细
      )
    )
  )
# ==============================================================================
# 3. 图 B：分块中介比例条形图 (带高级底色)
# ==============================================================================
p_bar <- ggplot(df_viz_plot, aes(y = Clean_Label)) +
  # --- 第一层：画区块底色 ---
  geom_rect(data = bg_facet, aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = Exposure_Tag), inherit.aes = FALSE, alpha = 0.8) +
  scale_fill_manual(values = facet_bg_colors, guide = "none") +
  new_scale_fill() +
  
  # --- 第二层：画数据 ---
  geom_col(aes(x = pool_prop * 100, fill = Group), width = 0.6, alpha = 0.9) +
  geom_text(aes(x = pool_prop * 100, label = sprintf("%.1f%%", pool_prop * 100)), hjust = -0.1, size = 3, fontface = "bold", color = "grey20") +
  
  scale_fill_manual(values = nature_pal) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) + 
  
  facet_grid(Exposure_Tag ~ ., scales = "free_y", space = "free_y") +
  
  theme_minimal() +
  labs(x = "Mediation Proportion (%)", y = NULL) +
  theme(# 🚀 重点在这里：强制不显示图例，避免 patchwork 收集到两套
    legend.position = "none",
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_line(color = "white", linewidth = 0.6),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text = element_blank(),
    panel.spacing = unit(0.3, "lines")
  )

# ==============================================================================
# 4. 完美拼图与输出
# ==============================================================================
final_figure <- p_forest + p_bar + 
  plot_layout(widths = c(1.8, 1), guides = "collect") + 
  plot_annotation(
    title = "Causal Mediation Atlas: Liver-AF Axis by Key Metabolites",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

print(final_figure)
ggsave("Figure_Mediation_Faceted_Premium.pdf", width = 14, height = 10)