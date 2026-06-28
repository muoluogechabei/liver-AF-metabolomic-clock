library(mice)
library(survival)
library(broom)
library(dplyr)
library(stringr)
library(flextable)
library(officer)

# 假设你的插补数据对象叫 imp_ready_for_cox
# imp_ready_for_cox <- mice_rediag 

# =======================================================
# 1. 数据预处理：缩尾 + 新的切点分类逻辑
# =======================================================
cat("正在进行数据预处理：99%缩尾 + 年龄分层切点 + 三分位分组...\n")

# 提取长格式数据
long_dat <- complete(imp_ready_for_cox, "long", include = TRUE)

# 定义缩尾函数
winsor_p99 <- function(x) {
  limit <- quantile(x, 0.99, na.rm = TRUE)
  x[x > limit] <- limit
  return(x)
}

long_dat_clean <- long_dat %>%
  group_by(.imp) %>% 
  mutate(
    # --- A. 连续变量缩尾 (保留操作) ---
    fib4 = winsor_p99(fib4),
    
    nfs  = winsor_p99(nfs),
    ast_alt_ratio = winsor_p99(ast_alt_ratio),
    
    # --- B. FIB-4 分类 (年龄分层) ---
    # <65岁: <1.3 / 1.3-2.67 / >=2.67
    # >=65岁: <2.0 / 2.0-2.67 / >=2.67 (你要求的修正切点)
    fib4_cat = case_when(
      # 65岁以下 (标准)
      age_val < 65 & fib4 < 1.3 ~ "Low Risk",
      age_val < 65 & fib4 >= 1.3 & fib4 < 2.67 ~ "Intermediate Risk",
      age_val < 65 & fib4 >= 2.67 ~ "High Risk",
      
      # 65岁及以上 (修正)
      age_val >= 65 & fib4 < 2.0 ~ "Low Risk",
      age_val >= 65 & fib4 >= 2.0 & fib4 < 2.67 ~ "Intermediate Risk",
      age_val >= 65 & fib4 >= 2.67 ~ "High Risk",
      
      TRUE ~ NA_character_
    ),
    
    # --- C. NFS 分类 (年龄分层) ---
    # <65岁: <-1.455 / -1.455~0.675 / >=0.675
    # >=65岁: <0.12 / 0.12~0.675 / >=0.675 (你要求的修正切点)
    nfs_cat = case_when(
      # 65岁以下 (标准)
      age_val < 65 & nfs < -1.455 ~ "Low Risk",
      age_val < 65 & nfs >= -1.455 & nfs < 0.675 ~ "Intermediate Risk",
      age_val < 65 & nfs >= 0.675 ~ "High Risk",
      
      # 65岁及以上 (修正)
      age_val >= 65 & nfs < 0.12 ~ "Low Risk",
      age_val >= 65 & nfs >= 0.12 & nfs < 0.675 ~ "Intermediate Risk",
      age_val >= 65 & nfs >= 0.675 ~ "High Risk",
      
      TRUE ~ NA_character_
    ),
    
    # --- D. APRI 分类 (两分类) ---
    # 切点: 1.0
    apri_cat = case_when(
      apri < 1.0 ~ "Low Risk",
      apri >= 1.0 ~ "High Risk",
      TRUE ~ NA_character_
    ),
    
    # --- E. AST/ALT 分类 (三分位数) ---
    # 动态计算 Low/Intermediate/High
    ast_alt_rank = ntile(ast_alt_ratio, 3), # 分为3组
    ast_alt_cat = case_when(
      ast_alt_rank == 1 ~ "Low Risk",
      ast_alt_rank == 2 ~ "Intermediate Risk",
      ast_alt_rank == 3 ~ "High Risk"
    )
  ) %>%
  ungroup() %>%
  # --- F. 设置因子水平 (确保 Low Risk 是参照组) ---
  mutate(
    fib4_cat = factor(fib4_cat, levels = c("Low Risk", "Intermediate Risk", "High Risk")),
    nfs_cat  = factor(nfs_cat,  levels = c("Low Risk", "Intermediate Risk", "High Risk")),
    apri_cat = factor(apri_cat, levels = c("Low Risk", "High Risk")), # APRI 只有两组
    ast_alt_cat = factor(ast_alt_cat, levels = c("Low Risk", "Intermediate Risk", "High Risk"))
  )

# 重新封装回 mids 对象
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
# 3. Cox 函数 (包含计数功能) - 保持逻辑不变
# =======================================================
run_combined_cox_sd <- function(imp_obj, exposure_name, covariates, model_label) {
  
  # --- A. 运行 Cox 模型 ---
  fit <- with(imp_obj, {
    # 直接获取因子变量
    x <- get(exposure_name) 
    
    # 构建公式
    form <- as.formula(paste("Surv(duration_updated, event_af_updated) ~ x +", 
                             paste(covariates, collapse = " + ")))
    coxph(form)
  })
  
  pooled <- pool(fit)
  
  # --- B. 提取回归结果 ---
  res <- summary(pooled, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(str_starts(term, "x")) %>% 
    mutate(
      Model = model_label,
      Exposure = toupper(str_replace(exposure_name, "_cat", "")),
      Type = "Categorical", # 既然只跑分类变量，这里固定为Categorical
      # 清洗 Group 名称 (移除 x 前缀)
      Group_Raw = str_replace(term, "^x", ""), 
      HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
      P_val = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
    )
  
  # --- C. 计算 No. (Total/Cases) ---
  dat_rep <- complete(imp_obj, 1)
  event_col <- "event_af_updated" 
  
  counts_df <- dat_rep %>%
    rename(Exp_Var = all_of(exposure_name)) %>%
    group_by(Exp_Var) %>%
    summarise(
      N_Total = n(),
      N_Cases = sum(.data[[event_col]], na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    mutate(
      Counts = paste0(N_Total, "/", N_Cases),
      Group_Raw = as.character(Exp_Var) 
    )
  
  # 将计数合并回结果表
  res <- res %>%
    left_join(counts_df %>% select(Group_Raw, Counts), by = "Group_Raw")
  
  return(res)
}

# =======================================================
# 4. 批量运行 (仅包含分类变量)
# =======================================================
# 修改此处：只保留 _cat 后缀的变量
exposures <- c("fib4_cat", "nfs_cat", "apri_cat", "ast_alt_cat")
all_res <- list()

for (expo in exposures) {
  for (mod_name in names(model_list)) {
    message(sprintf("Running: %s | %s", expo, mod_name))
    
    tryCatch({
      res_df <- run_combined_cox_sd(imp_final, expo, model_list[[mod_name]], mod_name)
      all_res[[paste(expo, mod_name)]] <- res_df
    }, error = function(e) {
      message(paste("Error in", expo, mod_name, ":", e$message))
    })
  }
}

final_table_df <- bind_rows(all_res) %>%
  select(Model, Exposure, Type, Group_Raw, Counts, HR_CI, P_val) %>%
  rename(Group = Group_Raw) # 重命名以便后续处理

# =======================================================
# 5. 生成表格并导出
# =======================================================
ft <- flextable(final_table_df) %>%
  set_header_labels(
    Model = "Model", 
    Exposure = "Indicator", 
    Type = "Type", 
    Group = "Risk Level", 
    Counts = "No. (Total/Cases)", 
    HR_CI = "HR (95% CI)", 
    P_val = "P value"
  ) %>%
  merge_v(j = c("Model", "Exposure", "Type")) %>% 
  theme_booktabs() %>%
  autofit() %>%
  align(j = "Counts", align = "center", part = "all") %>%
  bold(i = ~ !str_detect(P_val, ">") & as.numeric(str_replace_all(P_val, "<", "")) < 0.05, j = "P_val")

# 导出
doc <- read_docx() %>%
  body_add_par("Table: Sensitivity Analysis - Age-Adjusted Cutoffs & Modified Groupings", style = "heading 1") %>%
  body_add_par("Note: FIB-4 and NFS cutoffs adjusted for age >= 65. AST/ALT categorized by tertiles. APRI cutoff at 1.0.", style = "Normal") %>%
  body_add_flextable(ft)

print(doc, target = "Sensitivity_Analysis_NewCutoffs.docx")