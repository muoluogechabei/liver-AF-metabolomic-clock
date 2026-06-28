library(mice)
library(survival)
library(broom)
library(dplyr)
library(stringr)
library(flextable)
library(officer)
#imp_final <- imp_mri_ready_pdff
# 对每个插补数据集中的 PDFF 进行 Log + Z 转换
long_data <- complete(imp_final, "long", include = TRUE)

long_data <- long_data %>%
  group_by(.imp) %>%
  mutate(
    # 1. 加1取对数处理偏态 2. scale做Z分数转换
    pdff_z = as.numeric(scale(log(`Proton density fat fraction (PDFF) | Instance 2` + 1)))
  ) %>%
  ungroup()

imp_final <- as.mids(long_data)
#saveRDS(imp_final, "imp_final.rds")
# =======================================================
# 1. 数据准备
# =======================================================
# 假设 imp_final 是你已经包含 ct1_group 和 pdff_5_f 的插补对象
# imp_final <- imp_mri_ready_pdff 

cat("MRI 子队列分析准备就绪...\n")

# =======================================================
# 2. 变量与模型定义 
# =======================================================

# --- A. 核心变量名定义 (关键修改) ---
ct1_cont_var <- "Liver iron corrected T1 (ct1) | Instance 2"
# 【注意】这里一定要改成你新定义的 "ct1_group" (High/Int/Low)
ct1_cat_var  <- "ct1_f" 

# --- B. 协变量定义 ---
prs_var <- "`Standard PRS for atrial fibrillation (AF)`"
pc_vars <- paste0("`Genetic principal components | Array ", 1:10, "`") 
# 修改协变量定义，使其符合 EPV 原则
cov_m1 <- c("age_mri", "sex_f", "I(ethnicity_f == 'White')") # 种族在此临时二分

cov_m2 <- c(cov_m1, "tdi_val", "bmi_i2", "pdff_z")

cov_m3 <- c(cov_m2, "htn_mri_f", "t2dm_mri_f")

# 关键：手动指定 PC 1-3，不要用 contains() 全部塞进去
pc_vars_slim <- c("`Genetic principal components | Array 1`", 
                  "`Genetic principal components | Array 2`", 
                  "`Genetic principal components | Array 3`")

cov_m4 <- c(cov_m3, "`Standard PRS for atrial fibrillation (AF)`", pc_vars_slim)

# 模型列表
model_list <- list("Model 1" = cov_m1, "Model 2" = cov_m2, "Model 3" = cov_m3, "Model 4" = cov_m4)

# =======================================================
# 3. 升级版 Cox 函数 (自动计算 N 和 Cases)
# =======================================================
run_combined_cox_sd_mri <- function(imp_obj, exposure_name, covariates, model_label) {
  
  # --- 1. 提取第1个插补数据集用于计算 N 和 Cases ---
  # (多重插补中，N和结局状态通常是不变的，取第一个即可)
  df_count <- complete(imp_obj, 1)
  
  fit <- with(imp_obj, {
    raw_x <- get(exposure_name)
    
    # 判断是否为分类变量
    if(str_detect(exposure_name, "_cat|_f|group")) {
      x <- factor(raw_x)
      lvls <- levels(x)
      # 自动找 Low/Reference
      ref <- lvls[grep("Low|Normal|<|0", lvls, ignore.case = TRUE)[1]]
      if(!is.na(ref)) x <- relevel(x, ref = ref)
    } else {
      # 连续变量标准化
      x <- as.numeric(scale(raw_x))
    }
    
    form <- as.formula(paste("Surv(duration_mri, event_af_mri) ~ x +", 
                             paste(covariates, collapse = " + ")))
    coxph(form)
  })
  
  pooled <- pool(fit)
  
  # --- 2. 提取回归结果 ---
  res <- summary(pooled, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "x" | str_starts(term, "x")) %>% 
    mutate(
      Model = model_label,
      # 清理暴露名称
      Exposure = str_remove_all(exposure_name, " \\| Instance 2|_f|_cat|_group"), 
      Exposure = ifelse(str_detect(exposure_name, "ct1"), "Liver cT1", Exposure),
      Type = if_else(str_detect(exposure_name, "_cat|_f|group"), "Categorical", "Continuous"),
      # 提取组别名称
      Group = str_replace(term, "^x", ""), 
      Group = if_else(Group == "" | is.na(Group), "Per SD increase", Group),
      # 格式化 HR
      HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
      P_val = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
    )
  
  # --- 3. 计算并合并 N/Cases ---
  if(str_detect(exposure_name, "_cat|_f|group")) {
    # 分类变量：按组统计
    # 注意：这里需要确保 df_count 里的变量也是 factor 并且 levels 一致
    group_stats <- df_count %>%
      mutate(temp_group = as.character(get(exposure_name))) %>%
      group_by(temp_group) %>%
      summarise(
        N = n(),
        Cases = sum(event_af_mri == 1, na.rm = TRUE)
      ) %>%
      mutate(N_Cases = paste0(N, "/", Cases)) %>%
      select(Group = temp_group, N_Cases)
    
    # 关联回结果表
    res <- res %>% left_join(group_stats, by = "Group")
    
  } else {
    # 连续变量：统计总体
    total_n <- nrow(df_count)
    total_cases <- sum(df_count$event_af_mri == 1, na.rm = TRUE)
    res$N_Cases <- paste0(total_n, "/", total_cases)
  }
  
  return(res)
}

# =======================================================
# 4. 批量运行
# =======================================================
exposures_mri <- c(ct1_cont_var, ct1_cat_var)
all_res_mri <- list()

for (expo in exposures_mri) {
  for (mod_name in names(model_list)) {
    short_name <- ifelse(nchar(expo) > 20, substr(expo, 1, 15), expo)
    message(sprintf("Running: %s | %s", short_name, mod_name))
    
    # 运行函数
    all_res_mri[[paste(expo, mod_name)]] <- run_combined_cox_sd_mri(imp_final, expo, model_list[[mod_name]], mod_name)
  }
}

final_table_mri <- bind_rows(all_res_mri) %>%
  select(Model, Exposure, Type, Group, N_Cases, HR_CI, P_val) # 加入 N_Cases 列

# =======================================================
# 5. 生成表格 (含 N/Cases)
# =======================================================
ft_mri <- flextable(final_table_mri) %>%
  set_header_labels(
    Model="Model", Exposure="Exposure", Type="Type", 
    Group="Group", N_Cases = "No. (Total/Cases)", # 设置表头
    HR_CI="HR (95% CI)", P_val="P Value"
  ) %>%
  merge_v(j = c("Model", "Exposure", "Type")) %>% 
  theme_booktabs() %>%
  autofit() %>%
  bold(i = ~ !str_detect(P_val, ">") & as.numeric(str_replace_all(P_val, "<", "")) < 0.05, j = "P_val")

# 导出
doc_mri <- read_docx() %>%
  body_add_par("Table 2. Association of Liver cT1 with Incident AF (MRI Sub-cohort)", style = "heading 1") %>%
  body_add_par("Note: Stratification: Low (<750ms), High (>=750ms).", style = "Normal") %>%
  body_add_flextable(ft_mri)

print(doc_mri, target = "MRI_Subcohort_Results_fi.docx")

# 预览
print(final_table_mri)










# =======================================================
# 1. 数据准备
# =======================================================
# 假设 imp_final 是你已经包含 ct1_group 和 pdff_5_f 的插补对象
# imp_final <- imp_mri_ready_pdff 

cat("MRI 子队列分析准备就绪...\n")

# =======================================================
# 2. 变量与模型定义 
# =======================================================

# --- A. 核心变量名定义 (关键修改) ---
ct1_cont_var <- "Liver iron corrected T1 (ct1) | Instance 2"
# 【注意】这里一定要改成你新定义的 "ct1_group" (High/Int/Low)
ct1_cat_var  <- "ct1_f" 

# --- B. 协变量定义 ---
prs_var <- "`Standard PRS for atrial fibrillation (AF)`"
pc_vars <- paste0("`Genetic principal components | Array ", 1:10, "`") 
# 修改协变量定义，使其符合 EPV 原则
cov_m1 <- c("age_mri", "sex_f", "I(ethnicity_f == 'White')") # 种族在此临时二分
cov_m2 <- c(cov_m1, "tdi_val", "bmi_i2", "pdff_z")
cov_m3 <- c(cov_m2, "htn_mri_f", "t2dm_mri_f")

# 关键：手动指定 PC 1-3，不要用 contains() 全部塞进去
pc_vars_slim <- c("`Genetic principal components | Array 1`", 
                  "`Genetic principal components | Array 2`", 
                  "`Genetic principal components | Array 3`")

cov_m4 <- c(cov_m3, "`Standard PRS for atrial fibrillation (AF)`", pc_vars_slim)

# 🚨 修改 1：在模型列表中新增 Crude 模型
model_list <- list("Crude" = NULL, "Model 1" = cov_m1, "Model 2" = cov_m2, "Model 3" = cov_m3, "Model 4" = cov_m4)

# =======================================================
# 3. 升级版 Cox 函数 (自动计算 N 和 Cases + 强行添加 Ref)
# =======================================================
run_combined_cox_sd_mri <- function(imp_obj, exposure_name, covariates, model_label) {
  
  # --- 1. 提取第1个插补数据集用于计算 N 和 Cases ---
  df_count <- complete(imp_obj, 1)
  
  fit <- with(imp_obj, {
    raw_x <- get(exposure_name)
    
    # 判断是否为分类变量
    if(str_detect(exposure_name, "_cat|_f|group")) {
      x <- factor(raw_x)
      lvls <- levels(x)
      ref <- lvls[grep("Low|Normal|<|0", lvls, ignore.case = TRUE)[1]]
      if(!is.na(ref)) x <- relevel(x, ref = ref)
    } else {
      x <- as.numeric(scale(raw_x))
    }
    
    # 🚨 修改 2：针对 Crude 模型动态生成公式
    if(is.null(covariates) || length(covariates) == 0) {
      form <- as.formula("Surv(duration_mri, event_af_mri) ~ x")
    } else {
      form <- as.formula(paste("Surv(duration_mri, event_af_mri) ~ x +", 
                               paste(covariates, collapse = " + ")))
    }
    coxph(form)
  })
  
  pooled <- pool(fit)
  
  # 统一并清理暴露名称
  expo_clean <- str_remove_all(exposure_name, " \\| Instance 2|_f|_cat|_group") 
  expo_clean <- ifelse(str_detect(exposure_name, "ct1"), "Liver cT1", expo_clean)
  var_type <- if_else(str_detect(exposure_name, "_cat|_f|group"), "Categorical", "Continuous")
  
  # --- 2. 提取回归结果 ---
  res <- summary(pooled, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "x" | str_starts(term, "x")) %>% 
    mutate(
      Model = model_label,
      Exposure = expo_clean,
      Type = var_type,
      Group = str_replace(term, "^x", ""), 
      Group = if_else(Group == "" | is.na(Group), "Per SD increase", Group),
      HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
      P_val = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
    )
  
  # --- 3. 计算并合并 N/Cases (🚨 修改 3：强行生成参照组) ---
  if(var_type == "Categorical") {
    x_fac <- factor(df_count[[exposure_name]])
    lvls <- levels(x_fac)
    ref <- lvls[grep("Low|Normal|<|0", lvls, ignore.case = TRUE)[1]]
    if(is.na(ref)) ref <- lvls[1]
    
    group_stats <- df_count %>%
      mutate(temp_group = as.character(get(exposure_name))) %>%
      group_by(temp_group) %>%
      summarise(N = n(), Cases = sum(event_af_mri == 1, na.rm = TRUE), .groups = "drop") %>%
      mutate(N_Cases = paste0(N, "/", Cases))
    
    # 关联回结果表
    res <- res %>% left_join(group_stats %>% select(Group = temp_group, N_Cases), by = "Group")
    
    # 🌟 生成参照组 Row
    ref_count_val <- group_stats %>% filter(temp_group == ref) %>% pull(N_Cases)
    ref_row <- tibble(
      Model = model_label,
      Exposure = expo_clean,
      Type = "Categorical",
      Group = paste0(ref, " (Ref)"),
      N_Cases = ref_count_val,
      HR_CI = "1.00 (Reference)",
      P_val = "-"
    )
    
    # 合并参照组与干预组，并恢复正确的分类顺序
    res_final <- res %>%
      select(Model, Exposure, Type, Group, N_Cases, HR_CI, P_val) %>%
      bind_rows(ref_row, .) %>%
      mutate(Order_Col = match(str_remove(Group, " \\(Ref\\)"), lvls)) %>%
      arrange(Order_Col) %>%
      select(-Order_Col)
    
  } else {
    total_n <- nrow(df_count)
    total_cases <- sum(df_count$event_af_mri == 1, na.rm = TRUE)
    res_final <- res %>%
      mutate(N_Cases = paste0(total_n, "/", total_cases)) %>%
      select(Model, Exposure, Type, Group, N_Cases, HR_CI, P_val)
  }
  
  return(res_final)
}

# =======================================================
# 4. 批量运行
# =======================================================
exposures_mri <- c(ct1_cont_var, ct1_cat_var)
all_res_mri <- list()

for (expo in exposures_mri) {
  for (mod_name in names(model_list)) {
    short_name <- ifelse(nchar(expo) > 20, substr(expo, 1, 15), expo)
    message(sprintf("Running: %s | %s", short_name, mod_name))
    
    # 运行函数
    all_res_mri[[paste(expo, mod_name)]] <- run_combined_cox_sd_mri(imp_final, expo, model_list[[mod_name]], mod_name)
  }
}

final_table_mri <- bind_rows(all_res_mri) %>%
  select(Model, Exposure, Type, Group, N_Cases, HR_CI, P_val)

# =======================================================
# 5. 生成表格 (含 N/Cases)
# =======================================================
ft_mri <- flextable(final_table_mri) %>%
  set_header_labels(
    Model="Model", Exposure="Exposure", Type="Type", 
    Group="Group", N_Cases = "No. (Total/Cases)",
    HR_CI="HR (95% CI)", P_val="P Value"
  ) %>%
  merge_v(j = c("Model", "Exposure", "Type")) %>% 
  theme_booktabs() %>%
  autofit() %>%
  bold(i = ~ !str_detect(P_val, ">|-") & as.numeric(str_replace_all(P_val, "<", "")) < 0.05, j = "P_val") # 🚨 增加跳过"-"的保护

# 导出
doc_mri <- read_docx() %>%
  body_add_par("Table 2. Association of Liver cT1 with Incident AF (MRI Sub-cohort)", style = "heading 1") %>%
  body_add_par("Note: Stratification: Low (<750ms), High (>=750ms).", style = "Normal") %>%
  body_add_flextable(ft_mri)

print(doc_mri, target = "MRI_Subcohort_Results_fi.docx")

# 预览
print(final_table_mri)