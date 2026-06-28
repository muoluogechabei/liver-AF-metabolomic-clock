library(data.table)
library(dplyr)
library(mice)
# 1. 提取 data_sup 的 Participant ID 并与 data_drink 合并

df_combined <-fread("databcaa.csv", na.strings = c("", "NA"))

# 2. 输出合并后所有列的名字
names(df_combined)



library(dplyr)

library(dplyr)
library(tidyr)

library(dplyr)
library(tidyr)

df_risk <- df_combined  %>%
  # 1. 先把所有相关的饮酒摄入量列转为数值型
  # 这里建议直接列出这些 ID 变量，确保转换准确
  mutate(across(c(!!sym(red_status), !!sym(red_status2), 
                  !!sym(white_wine_status), !!sym(white_wine_status2),
                  !!sym(beer_status), !!sym(beer_status2),
                  !!sym(spirit_status), !!sym(spirit_status2),
                  !!sym(fortified_status), !!sym(fortified_status2),
                  !!sym(other_status), !!sym(other_status2)), 
                ~as.numeric(as.character(.)))) %>%
  
  # 2. 进行酒精克数计算 (Instance 0 和 Instance 2)
  mutate(
    alc_g_inst0 = 
      (replace_na(!!sym(red_status), 0) * 12) +
      (replace_na(!!sym(white_wine_status), 0) * 12) +
      (replace_na(!!sym(beer_status), 0) * 16) +
      (replace_na(!!sym(spirit_status), 0) * 8) +
      (replace_na(!!sym(fortified_status), 0) * 8) +
      (replace_na(!!sym(other_status), 0) * 8),
    
    alc_g_inst2 = 
      (replace_na(!!sym(red_status2), 0) * 12) +
      (replace_na(!!sym(white_wine_status2), 0) * 12) +
      (replace_na(!!sym(beer_status2), 0) * 16) +
      (replace_na(!!sym(spirit_status2), 0) * 8) +
      (replace_na(!!sym(fortified_status2), 0) * 8) +
      (replace_na(!!sym(other_status2), 0) * 8)
  ) %>%
  
  # 3. 判断 High Risk (结合性别)
  # 假设 Sex 列名没变，如果变了请替换为对应的 ID (如 p31)
  mutate(
    high_risk_inst0 = case_when(
      p31 %in% c("Female", 0) & alc_g_inst0 > 280 ~ 1,
      p31 %in% c("Male", 1)   & alc_g_inst0 > 400 ~ 1,
      is.na(p31) | is.na(alc_g_inst0) ~ NA_real_,
      TRUE ~ 0
    ),
    high_risk_inst2 = case_when(
      p31 %in% c("Female", 0) & alc_g_inst2 > 280 ~ 1,
      p31 %in% c("Male", 1)   & alc_g_inst2 > 400 ~ 1,
      is.na(p31) | is.na(alc_g_inst2) ~ NA_real_,
      TRUE ~ 0
    )
  )
# 1. 提取长格式数据 (包含原始数据 .imp == 0)
imp_long <- complete(mice_rediag, action = "long", include = TRUE)

# 2. 合并你刚才算的“克数法”高风险指标 (假设 df_risk 里有 Participant ID 和 high_risk_inst0)
# 如果 high_risk_inst0 已经在 imp_long 里了，可以跳过左连接
imp_long_to_filter <- imp_long %>%
  left_join(
    df_risk %>% select(eid, high_risk_inst0), 
    by = c("Participant ID" = "eid") 
  )

# 3. 执行严格剔除
imp_long_filtered <- imp_long_to_filter %>%
  filter(
    # 逻辑 A: 排除基线已有心血管疾病
    cvd_f == 0,
    
    # 逻辑 B: 排除随访前 2 年内发生房颤的（Landmark Analysis）
    !(event_af_updated == 1 & duration_updated < 2),
    
    # 逻辑 C: 排除有害饮酒者 (基于你刚才算的 280g/400g 阈值)
    high_risk_inst0 == 0, 
    
    # 逻辑 D: [新增] 排除前饮酒者 (Former Drinkers)
    # 根据 UKB 编码：0=Never, 1=Previous, 2=Current
    # 这里的 alc_status_f 只要不是 "Previous" 或 "1"
    !(as.character(alc_status_f) %in% c("Previous", "1"))
  )

# 4. 重新封装回 mids 对象 (供后续 Cox/Fine-Gray 使用)
imp_final_main <- as.mids(imp_long_filtered)

# 5. 检查剔除后的人数
cat("原始总样本量:", nrow(complete(mice_rediag, 1)), "\n")
cat("剔除后样本量:", nrow(complete(imp_final_main, 1)), "\n")

library(mice)
library(survival)
library(broom)
library(dplyr)
library(stringr)
library(flextable)
library(officer)
imp_ready_for_cox<-imp_final_main
# =======================================================
# 1. 数据预处理：99% 缩尾 (Winsorization) - 保持不变
# =======================================================
cat("正在对 FIB-4 和 APRI 进行 99% 缩尾处理...\n")

long_dat <- complete(imp_ready_for_cox, "long", include = TRUE)

winsor_p99 <- function(x) {
  limit <- quantile(x, 0.99, na.rm = TRUE)
  x[x > limit] <- limit
  return(x)
}

long_dat_clean <- long_dat %>%
  group_by(.imp) %>% 
  mutate(
    fib4 = winsor_p99(fib4),
    apri = winsor_p99(apri),
    nfs = winsor_p99(nfs),
    ast_alt_ratio = winsor_p99(ast_alt_ratio) 
  ) %>%
  ungroup()

imp_final <- as.mids(long_dat_clean)

# =======================================================
# 2. 变量与模型定义
# =======================================================
prs_var <- "`Standard PRS for atrial fibrillation (AF)`"
pc_vars <- paste0("`Genetic principal components | Array ", 1:10, "`")

cov_m1 <- c("age_val", "sex_f", "ethnicity_f")
cov_m2 <- c(cov_m1, "bmi_val", "smoke_f", "alc_freq_f", "tdi_val")
cov_m3 <- c(cov_m2, "htn_final", "t2dm_final", "lip_f", "ckd_final", "cvd_f")
cov_m4 <- c(cov_m3, prs_var, pc_vars)

model_list <- list("Model 1" = cov_m1, "Model 2" = cov_m2, "Model 3" = cov_m3, "Model 4" = cov_m4)

# =======================================================
# 3. Cox 函数 (修改：增加计数计算)
# =======================================================
run_combined_cox_sd <- function(imp_obj, exposure_name, covariates, model_label) {
  
  # --- A. 运行 Cox 模型 (保持原逻辑) ---
  fit <- with(imp_obj, {
    raw_x <- get(exposure_name)
    
    if(str_detect(exposure_name, "_cat")) {
      x <- factor(raw_x)
      lvls <- levels(x)
      # 自动寻找 Reference
      ref <- lvls[grep("Low|Normal|0", lvls, ignore.case = TRUE)[1]]
      if(!is.na(ref)) x <- relevel(x, ref = ref)
    } else {
      # 连续变量缩尾后标准化
      x <- as.numeric(scale(raw_x))
    }
    
    form <- as.formula(paste("Surv(duration_updated, event_af_updated) ~ x +", 
                             paste(covariates, collapse = " + ")))
    coxph(form)
  })
  
  pooled <- pool(fit)
  
  # --- B. 提取回归结果 ---
  res <- summary(pooled, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "x" | str_starts(term, "x")) %>% 
    mutate(
      Model = model_label,
      Exposure = toupper(str_replace(exposure_name, "_cat", "")),
      Type = if_else(str_detect(exposure_name, "_cat"), "Categorical", "Continuous"),
      # 清洗 Group 名称，移除 "x" 前缀以便后续匹配
      Group_Raw = str_replace(term, "^x", ""), 
      Group = if_else(Group_Raw == "" | is.na(Group_Raw), "Per SD increase", Group_Raw),
      HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
      P_val = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
    )
  
  # --- C. 【新增】计算 No. (Total/Cases) ---
  # 使用第1个插补数据集来计算代表性的样本量
  dat_rep <- complete(imp_obj, 1)
  
  # 确保结局变量名正确，这里假设是 event_af_updated
  # 如果你的变量名不同，请在此处修改
  event_col <- "event_af_updated" 
  
  if(str_detect(exposure_name, "_cat")) {
    # 1. 分类变量：按组计算
    counts_df <- dat_rep %>%
      rename(Exp_Var = all_of(exposure_name)) %>%
      group_by(Exp_Var) %>%
      summarise(
        N_Total = n(),
        N_Cases = sum(.data[[event_col]], na.rm = TRUE)
      ) %>%
      mutate(
        Counts = paste0(N_Total, "/", N_Cases),
        Group_Raw = as.character(Exp_Var) # 用于匹配 res 中的 Group_Raw
      )
    
    # 将计数合并回结果表
    res <- res %>%
      left_join(counts_df %>% select(Group_Raw, Counts), by = "Group_Raw")
    
  } else {
    # 2. 连续变量：计算全人群
    n_total <- nrow(dat_rep)
    n_cases <- sum(dat_rep[[event_col]], na.rm = TRUE)
    count_str <- paste0(n_total, "/", n_cases)
    
    res$Counts <- count_str
  }
  
  return(res)
}

# =======================================================
# 4. 批量运行
# =======================================================
exposures <- c("fib4", "fib4_cat", "nfs", "nfs_cat", "apri", "apri_cat", "ast_alt_ratio", "ratio_cat")
all_res <- list()

# 注意：这里需要确保 imp_final 里包含 covariates 和 outcome
for (expo in exposures) {
  for (mod_name in names(model_list)) {
    message(sprintf("Running: %s | %s", expo, mod_name))
    
    # 运行模型并获取含计数的结果
    res_df <- run_combined_cox_sd(imp_final, expo, model_list[[mod_name]], mod_name)
    
    # 简单清理一下列顺序
    all_res[[paste(expo, mod_name)]] <- res_df
  }
}

final_table_df <- bind_rows(all_res) %>%
  select(Model, Exposure, Type, Group, Counts, HR_CI, P_val)

# =======================================================
# 5. 生成表格并导出 (增加 Counts 列设置)
# =======================================================
ft <- flextable(final_table_df) %>%
  set_header_labels(
    Model = "Model", 
    Exposure = "Indicator", 
    Type = "Type", 
    Group = "Level / Change", 
    Counts = "No. (Total/Cases)",  # 新增表头
    HR_CI = "HR (95% CI)", 
    P_val = "P value"
  ) %>%
  # 合并相同单元格
  merge_v(j = c("Model", "Exposure", "Type")) %>% 
  theme_booktabs() %>%
  autofit() %>%
  # 居中对齐计数列
  align(j = "Counts", align = "center", part = "all") %>%
  # P值加粗逻辑
  bold(i = ~ !str_detect(P_val, ">") & as.numeric(str_replace_all(P_val, "<", "")) < 0.05, j = "P_val")

# 导出
doc <- read_docx() %>%
  body_add_par("Table: Liver Fibrosis Indices and Incident AF (Winsorized)", style = "heading 1") %>%
  body_add_par("Note: Continuous variables standardized after 99% winsorization.", style = "Normal") %>%
  body_add_flextable(ft)

print(doc, target = "Liver_特殊人群_Analysis_WithCounts.docx")





library(mice)
library(survival)
library(broom)
library(dplyr)
library(stringr)
library(flextable)
library(officer)
imp_ready_for_cox<-mice_sensi
# =======================================================
# 1. 数据预处理：99% 缩尾 (Winsorization) - 保持不变
# =======================================================
cat("正在对 FIB-4 和 APRI 进行 99% 缩尾处理...\n")
# 5. 检查剔除后的人数
cat("原始总样本量:", nrow(complete(mice_rediag, 1)), "\n")
cat("剔除后样本量:", nrow(complete(imp_ready_for_cox, 1)), "\n")

long_dat <- complete(imp_ready_for_cox, "long", include = TRUE)

winsor_p99 <- function(x) {
  limit <- quantile(x, 0.99, na.rm = TRUE)
  x[x > limit] <- limit
  return(x)
}

long_dat_clean <- long_dat %>%
  group_by(.imp) %>% 
  mutate(
    fib4 = winsor_p99(fib4),
    apri = winsor_p99(apri),
    nfs = winsor_p99(nfs),
    ast_alt_ratio = winsor_p99(ast_alt_ratio) 
  ) %>%
  ungroup()

imp_final <- as.mids(long_dat_clean)

# =======================================================
# 2. 变量与模型定义
# =======================================================
prs_var <- "`Standard PRS for atrial fibrillation (AF)`"
pc_vars <- paste0("`Genetic principal components | Array ", 1:10, "`")

cov_m1 <- c("age_val", "sex_f", "ethnicity_f")
cov_m2 <- c(cov_m1, "bmi_val", "smoke_f", "alc_freq_f", "tdi_val")
cov_m3 <- c(cov_m2, "htn_final", "t2dm_final", "lip_f", "ckd_final", "cvd_f")
cov_m4 <- c(cov_m3, prs_var, pc_vars)

model_list <- list("Model 1" = cov_m1, "Model 2" = cov_m2, "Model 3" = cov_m3, "Model 4" = cov_m4)

# =======================================================
# 3. Cox 函数 (修改：增加计数计算)
# =======================================================
run_combined_cox_sd <- function(imp_obj, exposure_name, covariates, model_label) {
  
  # --- A. 运行 Cox 模型 (保持原逻辑) ---
  fit <- with(imp_obj, {
    raw_x <- get(exposure_name)
    
    if(str_detect(exposure_name, "_cat")) {
      x <- factor(raw_x)
      lvls <- levels(x)
      # 自动寻找 Reference
      ref <- lvls[grep("Low|Normal|0", lvls, ignore.case = TRUE)[1]]
      if(!is.na(ref)) x <- relevel(x, ref = ref)
    } else {
      # 连续变量缩尾后标准化
      x <- as.numeric(scale(raw_x))
    }
    
    form <- as.formula(paste("Surv(duration_updated, event_af_updated) ~ x +", 
                             paste(covariates, collapse = " + ")))
    coxph(form)
  })
  
  pooled <- pool(fit)
  
  # --- B. 提取回归结果 ---
  res <- summary(pooled, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "x" | str_starts(term, "x")) %>% 
    mutate(
      Model = model_label,
      Exposure = toupper(str_replace(exposure_name, "_cat", "")),
      Type = if_else(str_detect(exposure_name, "_cat"), "Categorical", "Continuous"),
      # 清洗 Group 名称，移除 "x" 前缀以便后续匹配
      Group_Raw = str_replace(term, "^x", ""), 
      Group = if_else(Group_Raw == "" | is.na(Group_Raw), "Per SD increase", Group_Raw),
      HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
      P_val = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
    )
  
  # --- C. 【新增】计算 No. (Total/Cases) ---
  # 使用第1个插补数据集来计算代表性的样本量
  dat_rep <- complete(imp_obj, 1)
  
  # 确保结局变量名正确，这里假设是 event_af_updated
  # 如果你的变量名不同，请在此处修改
  event_col <- "event_af_updated" 
  
  if(str_detect(exposure_name, "_cat")) {
    # 1. 分类变量：按组计算
    counts_df <- dat_rep %>%
      rename(Exp_Var = all_of(exposure_name)) %>%
      group_by(Exp_Var) %>%
      summarise(
        N_Total = n(),
        N_Cases = sum(.data[[event_col]], na.rm = TRUE)
      ) %>%
      mutate(
        Counts = paste0(N_Total, "/", N_Cases),
        Group_Raw = as.character(Exp_Var) # 用于匹配 res 中的 Group_Raw
      )
    
    # 将计数合并回结果表
    res <- res %>%
      left_join(counts_df %>% select(Group_Raw, Counts), by = "Group_Raw")
    
  } else {
    # 2. 连续变量：计算全人群
    n_total <- nrow(dat_rep)
    n_cases <- sum(dat_rep[[event_col]], na.rm = TRUE)
    count_str <- paste0(n_total, "/", n_cases)
    
    res$Counts <- count_str
  }
  
  return(res)
}

# =======================================================
# 4. 批量运行
# =======================================================
exposures <- c("fib4", "fib4_cat", "nfs", "nfs_cat", "apri", "apri_cat", "ast_alt_ratio", "ratio_cat")
all_res <- list()

# 注意：这里需要确保 imp_final 里包含 covariates 和 outcome
for (expo in exposures) {
  for (mod_name in names(model_list)) {
    message(sprintf("Running: %s | %s", expo, mod_name))
    
    # 运行模型并获取含计数的结果
    res_df <- run_combined_cox_sd(imp_final, expo, model_list[[mod_name]], mod_name)
    
    # 简单清理一下列顺序
    all_res[[paste(expo, mod_name)]] <- res_df
  }
}

final_table_df <- bind_rows(all_res) %>%
  select(Model, Exposure, Type, Group, Counts, HR_CI, P_val)

# =======================================================
# 5. 生成表格并导出 (增加 Counts 列设置)
# =======================================================
ft <- flextable(final_table_df) %>%
  set_header_labels(
    Model = "Model", 
    Exposure = "Indicator", 
    Type = "Type", 
    Group = "Level / Change", 
    Counts = "No. (Total/Cases)",  # 新增表头
    HR_CI = "HR (95% CI)", 
    P_val = "P value"
  ) %>%
  # 合并相同单元格
  merge_v(j = c("Model", "Exposure", "Type")) %>% 
  theme_booktabs() %>%
  autofit() %>%
  # 居中对齐计数列
  align(j = "Counts", align = "center", part = "all") %>%
  # P值加粗逻辑
  bold(i = ~ !str_detect(P_val, ">") & as.numeric(str_replace_all(P_val, "<", "")) < 0.05, j = "P_val")

# 导出
doc <- read_docx() %>%
  body_add_par("Table: Liver Fibrosis Indices and Incident AF (Winsorized)", style = "heading 1") %>%
  body_add_par("Note: Continuous variables standardized after 99% winsorization.", style = "Normal") %>%
  body_add_flextable(ft)

print(doc, target = "Liver_血液系统肿瘤+急性肝损_Analysis_WithCounts.docx")




