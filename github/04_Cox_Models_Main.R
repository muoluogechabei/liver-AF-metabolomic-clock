library(mice)
library(survival)
library(broom)
library(dplyr)
library(stringr)
library(flextable)
library(officer)
library(tidyverse)
#imp_ready_for_cox<-mice_rediag
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
# 3. Cox 函数 (修复：增加参照组生成 + 计数计算)
# =======================================================
run_combined_cox_sd <- function(imp_obj, exposure_name, covariates, model_label) {
  
  # --- A. 运行 Cox 模型 ---
  fit <- with(imp_obj, {
    raw_x <- get(exposure_name)
    if(str_detect(exposure_name, "_cat")) {
      x <- factor(raw_x)
      lvls <- levels(x)
      ref <- lvls[grep("Low|Normal|0", lvls, ignore.case = TRUE)[1]]
      if(is.na(ref)) ref <- lvls[1]
      x <- relevel(x, ref = ref)
    } else {
      x <- as.numeric(scale(raw_x))
    }
    form <- as.formula(paste("Surv(duration_updated, event_af_updated) ~ x +", 
                             paste(covariates, collapse = " + ")))
    coxph(form)
  })
  
  pooled <- pool(fit)
  res_raw <- summary(pooled, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "x" | str_starts(term, "x"))
  
  # --- B. 提取计数并合并 ---
  dat_rep <- complete(imp_obj, 1)
  event_col <- "event_af_updated" 
  
  # 统一 Exposure 名称
  expo_base <- toupper(str_replace(exposure_name, "_cat", ""))
  if(expo_base == "RATIO") expo_base <- "AST_ALT_RATIO" 
  
  if(str_detect(exposure_name, "_cat")) {
    # 1. 分类变量
    x_fac <- factor(dat_rep[[exposure_name]])
    lvls <- levels(x_fac)
    ref <- lvls[grep("Low|Normal|0", lvls, ignore.case = TRUE)[1]]
    if(is.na(ref)) ref <- lvls[1]
    
    counts_df <- dat_rep %>%
      group_by(Exp_Var = as.character(.data[[exposure_name]])) %>%
      summarise(N_Total = n(), N_Cases = sum(.data[[event_col]], na.rm = TRUE)) %>%
      mutate(Counts = paste0(N_Total, "/", N_Cases))
    
    # 处理非参照组的回归结果
    res <- res_raw %>%
      mutate(
        Group_Raw = str_replace(term, "^x", ""),
        HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
        P_val = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
      ) %>%
      left_join(counts_df, by = c("Group_Raw" = "Exp_Var"))
    
    # 🌟 强行生成参照组 (Reference Row)
    ref_count <- counts_df %>% filter(Exp_Var == ref) %>% pull(Counts)
    ref_row <- tibble(
      Group_Raw = ref,
      HR_CI = "1.00 (Reference)",
      P_val = "-",
      Counts = ref_count
    )
    
    # 组合参照组与非参照组，并恢复原始 Level 排序
    res_final <- bind_rows(ref_row, res) %>%
      mutate(
        Model = model_label,
        Exposure = expo_base,
        Type = "Categorical",
        # 给参照组名字加上 (Ref) 标识
        Group = if_else(Group_Raw == ref, paste0(Group_Raw, " (Ref)"), Group_Raw),
        Order_Col = match(Group_Raw, lvls) # 用来排序保证 Reference 在第一行
      ) %>%
      arrange(Order_Col) %>%
      select(-Order_Col, -term, -estimate, -std.error, -statistic, -df, -p.value, -conf.low, -conf.high)
    
  } else {
    # 2. 连续变量
    n_total <- nrow(dat_rep)
    n_cases <- sum(dat_rep[[event_col]], na.rm = TRUE)
    
    res_final <- res_raw %>%
      mutate(
        Model = model_label,
        Exposure = expo_base,
        Type = "Continuous",
        Group_Raw = "Per 1 SD increase",
        Group = "Per 1 SD increase",
        HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
        P_val = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
        Counts = paste0(n_total, "/", n_cases)
      ) %>%
      select(-term, -estimate, -std.error, -statistic, -df, -p.value, -conf.low, -conf.high)
  }
  
  return(res_final)
}

# =======================================================
# 4. 批量运行与长宽表转换 (🌟 横向折叠压缩表格)
# =======================================================
exposures <- c("fib4", "fib4_cat", "nfs", "nfs_cat", "apri", "apri_cat", "ast_alt_ratio", "ratio_cat")
all_res <- list()

for (expo in exposures) {
  for (mod_name in names(model_list)) {
    message(sprintf("Running: %s | %s", expo, mod_name))
    res_df <- run_combined_cox_sd(imp_final, expo, model_list[[mod_name]], mod_name)
    all_res[[paste(expo, mod_name)]] <- res_df
  }
}

# 将长表格转换为紧凑的横向宽表
wide_table_df <- bind_rows(all_res) %>%
  mutate(
    Exposure_Clean = case_when(
      Exposure == "FIB4" ~ "FIB-4",
      Exposure == "NFS" ~ "NFS",
      Exposure == "APRI" ~ "APRI",
      Exposure == "AST_ALT_RATIO" ~ "AST/ALT Ratio",
      TRUE ~ Exposure
    ),
    # 确保排版：连续变量在上，分类变量在下
    Type = factor(Type, levels = c("Continuous", "Categorical")) 
  ) %>%
  select(Exposure_Clean, Type, Group, Counts, Model, HR_CI, P_val) %>%
  # 🌟 魔法步骤：以 Model 为列名，将 HR 和 P 值横向展开
  pivot_wider(
    names_from = Model,
    values_from = c(HR_CI, P_val),
    names_glue = "{Model}_{.value}"
  )

# =======================================================
# 5. 生成专业顶刊格式表格并导出
# =======================================================
library(flextable)
library(officer)

# 选择需要展示的列（排除 Type 列让表格更干净）
display_cols <- c("Exposure_Clean", "Group", "Counts", 
                  "Model 1_HR_CI", "Model 1_P_val", 
                  "Model 2_HR_CI", "Model 2_P_val", 
                  "Model 3_HR_CI", "Model 3_P_val", 
                  "Model 4_HR_CI", "Model 4_P_val")

ft <- flextable(wide_table_df[, display_cols]) %>%
  set_header_labels(
    Exposure_Clean = "Liver Index", 
    Group = "Level / Increment", 
    Counts = "No. (Total/Cases)", 
    `Model 1_HR_CI` = "HR (95% CI)", `Model 1_P_val` = "P value",
    `Model 2_HR_CI` = "HR (95% CI)", `Model 2_P_val` = "P value",
    `Model 3_HR_CI` = "HR (95% CI)", `Model 3_P_val` = "P value",
    `Model 4_HR_CI` = "HR (95% CI)", `Model 4_P_val` = "P value"
  ) %>%
  # 🌟 增加高阶二级表头 (跨列合并 Model 1 ~ 4)
  add_header_row(
    values = c("", "", "", "Model 1", "Model 2", "Model 3", "Model 4"), 
    colwidths = c(1, 1, 1, 2, 2, 2, 2)
  ) %>%
  merge_v(j = "Exposure_Clean") %>%          # 垂直合并相同的指标名
  valign(j = "Exposure_Clean", valign = "top") %>% 
  theme_booktabs() %>%                       # 经典三线表主题
  align(j = c("Exposure_Clean", "Group"), align = "left", part = "all") %>% # 文本左对齐
  align(j = 3:11, align = "center", part = "all") %>%                       # 数值居中
  autofit()

# 🌟 动态循环为所有的 P_val 列自动加粗显著性 P 值 (安全跳过 "-" 参照组字符)
p_cols <- grep("_P_val$", names(wide_table_df), value = TRUE)
for (p in p_cols) {
  ft <- ft %>%
    bold(i = as.formula(paste0("~ `", p, "` != '-' & as.numeric(gsub('<', '', `", p, "`)) < 0.05")), j = p)
}

# 导出文档
doc <- read_docx() %>%
  body_add_par("Table: Associations Between Liver Fibrosis Indices and Incident AF", style = "heading 1") %>%
  body_add_par("Note: Continuous variables standardized after 99% winsorization. (Ref) represents the baseline category.", style = "Normal") %>%
  body_add_flextable(ft)

print(doc, target = "Liver_AF_Analysis_CompactTable.docx")
cat("\n✅ 表格生成完毕！已横向压缩并补全参照组，保存为 Liver_AF_Analysis_CompactTable.docx\n")













library(mice)
library(survival)
library(broom)
library(dplyr)
library(stringr)
library(flextable)
library(officer)
library(tidyverse)

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
# 2. 变量与模型定义 (🚨 新增 Crude 模型)
# =======================================================
prs_var <- "`Standard PRS for atrial fibrillation (AF)`"
pc_vars <- paste0("`Genetic principal components | Array ", 1:10, "`")

cov_m1 <- c("age_val", "sex_f", "ethnicity_f")
cov_m2 <- c(cov_m1, "bmi_val", "smoke_f", "alc_freq_f", "tdi_val")
cov_m3 <- c(cov_m2, "htn_final", "t2dm_final", "lip_f", "ckd_final", "cvd_f")
cov_m4 <- c(cov_m3, prs_var, pc_vars)

# 🚨 在列表最前面加入 Crude 模型 (设为 NULL)
model_list <- list("Crude" = NULL, "Model 1" = cov_m1, "Model 2" = cov_m2, "Model 3" = cov_m3, "Model 4" = cov_m4)

# =======================================================
# 3. Cox 函数 (🚨 修复：增加无校正协变量的公式判断)
# =======================================================
run_combined_cox_sd <- function(imp_obj, exposure_name, covariates, model_label) {
  
  # --- A. 运行 Cox 模型 ---
  fit <- with(imp_obj, {
    raw_x <- get(exposure_name)
    if(str_detect(exposure_name, "_cat")) {
      x <- factor(raw_x)
      lvls <- levels(x)
      ref <- lvls[grep("Low|Normal|0", lvls, ignore.case = TRUE)[1]]
      if(is.na(ref)) ref <- lvls[1]
      x <- relevel(x, ref = ref)
    } else {
      x <- as.numeric(scale(raw_x))
    }
    
    # 🚨 针对 Crude 模型动态生成公式 (如果是 NULL 或空，则不加后面的 +)
    if(is.null(covariates) || length(covariates) == 0) {
      form <- as.formula("Surv(duration_updated, event_af_updated) ~ x")
    } else {
      form <- as.formula(paste("Surv(duration_updated, event_af_updated) ~ x +", 
                               paste(covariates, collapse = " + ")))
    }
    coxph(form)
  })
  
  pooled <- pool(fit)
  res_raw <- summary(pooled, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term == "x" | str_starts(term, "x"))
  
  # --- B. 提取计数并合并 ---
  dat_rep <- complete(imp_obj, 1)
  event_col <- "event_af_updated" 
  
  # 统一 Exposure 名称
  expo_base <- toupper(str_replace(exposure_name, "_cat", ""))
  if(expo_base == "RATIO") expo_base <- "AST_ALT_RATIO" 
  
  if(str_detect(exposure_name, "_cat")) {
    # 1. 分类变量
    x_fac <- factor(dat_rep[[exposure_name]])
    lvls <- levels(x_fac)
    ref <- lvls[grep("Low|Normal|0", lvls, ignore.case = TRUE)[1]]
    if(is.na(ref)) ref <- lvls[1]
    
    counts_df <- dat_rep %>%
      group_by(Exp_Var = as.character(.data[[exposure_name]])) %>%
      summarise(N_Total = n(), N_Cases = sum(.data[[event_col]], na.rm = TRUE)) %>%
      mutate(Counts = paste0(N_Total, "/", N_Cases))
    
    # 处理非参照组的回归结果
    res <- res_raw %>%
      mutate(
        Group_Raw = str_replace(term, "^x", ""),
        HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
        P_val = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
      ) %>%
      left_join(counts_df, by = c("Group_Raw" = "Exp_Var"))
    
    # 🌟 强行生成参照组 (Reference Row)
    ref_count <- counts_df %>% filter(Exp_Var == ref) %>% pull(Counts)
    ref_row <- tibble(
      Group_Raw = ref,
      HR_CI = "1.00 (Reference)",
      P_val = "-",
      Counts = ref_count
    )
    
    # 组合参照组与非参照组，并恢复原始 Level 排序
    res_final <- bind_rows(ref_row, res) %>%
      mutate(
        Model = model_label,
        Exposure = expo_base,
        Type = "Categorical",
        # 给参照组名字加上 (Ref) 标识
        Group = if_else(Group_Raw == ref, paste0(Group_Raw, " (Ref)"), Group_Raw),
        Order_Col = match(Group_Raw, lvls) # 用来排序保证 Reference 在第一行
      ) %>%
      arrange(Order_Col) %>%
      select(-Order_Col, -term, -estimate, -std.error, -statistic, -df, -p.value, -conf.low, -conf.high)
    
  } else {
    # 2. 连续变量
    n_total <- nrow(dat_rep)
    n_cases <- sum(dat_rep[[event_col]], na.rm = TRUE)
    
    res_final <- res_raw %>%
      mutate(
        Model = model_label,
        Exposure = expo_base,
        Type = "Continuous",
        Group_Raw = "Per 1 SD increase",
        Group = "Per 1 SD increase",
        HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
        P_val = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
        Counts = paste0(n_total, "/", n_cases)
      ) %>%
      select(-term, -estimate, -std.error, -statistic, -df, -p.value, -conf.low, -conf.high)
  }
  
  return(res_final)
}

# =======================================================
# 4. 批量运行与长宽表转换 (🌟 横向折叠压缩表格)
# =======================================================
exposures <- c("fib4", "fib4_cat", "nfs", "nfs_cat", "apri", "apri_cat", "ast_alt_ratio", "ratio_cat")
all_res <- list()

for (expo in exposures) {
  for (mod_name in names(model_list)) {
    message(sprintf("Running: %s | %s", expo, mod_name))
    res_df <- run_combined_cox_sd(imp_final, expo, model_list[[mod_name]], mod_name)
    all_res[[paste(expo, mod_name)]] <- res_df
  }
}

# 将长表格转换为紧凑的横向宽表
wide_table_df <- bind_rows(all_res) %>%
  mutate(
    Exposure_Clean = case_when(
      Exposure == "FIB4" ~ "FIB-4",
      Exposure == "NFS" ~ "NFS",
      Exposure == "APRI" ~ "APRI",
      Exposure == "AST_ALT_RATIO" ~ "AST/ALT Ratio",
      TRUE ~ Exposure
    ),
    # 确保排版：连续变量在上，分类变量在下
    Type = factor(Type, levels = c("Continuous", "Categorical")) 
  ) %>%
  select(Exposure_Clean, Type, Group, Counts, Model, HR_CI, P_val) %>%
  # 🌟 魔法步骤：以 Model 为列名，将 HR 和 P 值横向展开
  pivot_wider(
    names_from = Model,
    values_from = c(HR_CI, P_val),
    names_glue = "{Model}_{.value}"
  )

# =======================================================
# 5. 生成专业顶刊格式表格并导出 (🚨 适配 Crude 列的展示)
# =======================================================
library(flextable)
library(officer)

# 🚨 增加 Crude 列
display_cols <- c("Exposure_Clean", "Group", "Counts", 
                  "Crude_HR_CI", "Crude_P_val",
                  "Model 1_HR_CI", "Model 1_P_val", 
                  "Model 2_HR_CI", "Model 2_P_val", 
                  "Model 3_HR_CI", "Model 3_P_val", 
                  "Model 4_HR_CI", "Model 4_P_val")

ft <- flextable(wide_table_df[, display_cols]) %>%
  set_header_labels(
    Exposure_Clean = "Liver Index", 
    Group = "Level / Increment", 
    Counts = "No. (Total/Cases)", 
    `Crude_HR_CI` = "HR (95% CI)", `Crude_P_val` = "P value",
    `Model 1_HR_CI` = "HR (95% CI)", `Model 1_P_val` = "P value",
    `Model 2_HR_CI` = "HR (95% CI)", `Model 2_P_val` = "P value",
    `Model 3_HR_CI` = "HR (95% CI)", `Model 3_P_val` = "P value",
    `Model 4_HR_CI` = "HR (95% CI)", `Model 4_P_val` = "P value"
  ) %>%
  # 🌟 🚨 增加高阶二级表头 (跨列合并增加 Crude，colwidths 改为 8 项)
  add_header_row(
    values = c("", "", "", "Crude", "Model 1", "Model 2", "Model 3", "Model 4"), 
    colwidths = c(1, 1, 1, 2, 2, 2, 2, 2)
  ) %>%
  merge_v(j = "Exposure_Clean") %>%          # 垂直合并相同的指标名
  valign(j = "Exposure_Clean", valign = "top") %>% 
  theme_booktabs() %>%                       # 经典三线表主题
  align(j = c("Exposure_Clean", "Group"), align = "left", part = "all") %>% # 文本左对齐
  align(j = 3:13, align = "center", part = "all") %>%                       # 🚨 数据列对齐更新为 3:13
  autofit()

# 🌟 动态循环为所有的 P_val 列自动加粗显著性 P 值 (安全跳过 "-" 参照组字符)
p_cols <- grep("_P_val$", names(wide_table_df), value = TRUE)
for (p in p_cols) {
  ft <- ft %>%
    bold(i = as.formula(paste0("~ `", p, "` != '-' & as.numeric(gsub('<', '', `", p, "`)) < 0.05")), j = p)
}

# 导出文档
doc <- read_docx() %>%
  body_add_par("Table: Associations Between Liver Fibrosis Indices and Incident AF", style = "heading 1") %>%
  body_add_par("Note: Continuous variables standardized after 99% winsorization. (Ref) represents the baseline category.", style = "Normal") %>%
  body_add_flextable(ft)

print(doc, target = "Liver_AF_Analysis_CompactTable_WithCrude.docx")
cat("\n✅ 包含 Crude 模型的表格生成完毕！已横向压缩并补全参照组，保存为 Liver_AF_Analysis_CompactTable_WithCrude.docx\n")




# =======================================================
# 加载必备的包
# =======================================================
library(rms)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(mice)
library(survival)
library(stringr)

cat("准备绘制五拼图 RCS...\n")
imp_ready_for_cox<- mice_rediag 
# =======================================================
# 1. 准备数据 A：肝纤维化主队列数据 (Fibrosis cohort)
# =======================================================
df_rcs_fib <- complete(imp_ready_for_cox, 1)

# 清理复杂的协变量名 (避免 rms::cph 解析公式时报错)
names(df_rcs_fib)[names(df_rcs_fib) == "Standard PRS for atrial fibrillation (AF)"] <- "PRS_AF"
for(i in 1:10) {
  names(df_rcs_fib)[names(df_rcs_fib) == paste0("Genetic principal components | Array ", i)] <- paste0("PC", i)
}

# Model 4 协变量 (清理后)
cov_m1_fib <- c("age_val", "sex_f", "ethnicity_f")
cov_m2_fib <- c(cov_m1_fib, "bmi_val", "smoke_f", "alc_freq_f", "tdi_val")
cov_m3_fib <- c(cov_m2_fib, "htn_final", "t2dm_final", "lip_f", "ckd_final", "cvd_f")
cov_m4_fib <- c(cov_m3_fib, "PRS_AF", paste0("PC", 1:10))

# 数据缩尾处理
winsor_vars <- c("fib4", "apri", "ast_alt_ratio") # 99%单侧缩尾
df_rcs_fib <- df_rcs_fib %>%
  mutate(across(all_of(winsor_vars), function(x) {
    p99 <- quantile(x, 0.99, na.rm = TRUE)
    ifelse(x > p99, p99, x)
  })) %>%
  mutate(nfs = { # NFS 双侧 1%-99% 缩尾
    x <- nfs
    p01 <- quantile(x, 0.01, na.rm = TRUE)
    p99 <- quantile(x, 0.99, na.rm = TRUE)
    x[x < p01] <- p01
    x[x > p99] <- p99
    x
  })


# =======================================================
# 2. 准备数据 B：MRI 子队列数据 (MRI cohort)
# =======================================================
df_rcs_mri <- complete(imp_final, 1)

# 清理核心变量及协变量名称
df_rcs_mri <- df_rcs_mri %>%
  rename(
    ct1_val = `Liver iron corrected T1 (ct1) | Instance 2`,
    PRS_AF = `Standard PRS for atrial fibrillation (AF)`
  ) %>%
  # rms 不喜欢 I(ethnicity == 'White') 这种写法，直接转成实数分类
  mutate(white_f = as.numeric(ethnicity_f == 'White')) 

for(i in 1:3) {
  names(df_rcs_mri)[names(df_rcs_mri) == paste0("Genetic principal components | Array ", i)] <- paste0("PC", i)
}

# Model 4 协变量 (清理后)
cov_m1_mri <- c("age_mri", "sex_f", "white_f")
cov_m2_mri <- c(cov_m1_mri, "tdi_val", "bmi_i2", "pdff_z")
cov_m3_mri <- c(cov_m2_mri, "htn_mri_f", "t2dm_mri_f")
cov_m4_mri <- c(cov_m3_mri, "PRS_AF", paste0("PC", 1:3))

# cT1 数据进行 1% - 99% 双侧缩尾 (保证曲线首尾稳定)
df_rcs_mri <- df_rcs_mri %>%
  mutate(ct1_val = {
    x <- ct1_val
    p01 <- quantile(x, 0.01, na.rm = TRUE)
    p99 <- quantile(x, 0.99, na.rm = TRUE)
    x[x < p01] <- p01
    x[x > p99] <- p99
    x
  })


# =======================================================
# 3. 升级版绘制函数 (兼容多个不同数据集和时间变量)
# =======================================================
draw_my_rcs_robust <- function(data, exposure, label_name, cov_vars, time_var, event_var, line_color = "#2E508E") {
  
  # 强制赋予到全局环境，防止 rms 找不到 datadist 报错
  dd_temp <- datadist(data)
  assign("dd", dd_temp, envir = .GlobalEnv)
  options(datadist = "dd")
  
  # 构造公式 (4个节点)
  formula_str <- paste0("Surv(", time_var, ", ", event_var, ") ~ rcs(", exposure, ", 4) + ", 
                        paste(cov_vars, collapse = " + "))
  
  # 拟合模型
  fit <- cph(as.formula(formula_str), data = data, x=TRUE, y=TRUE)
  
  # 提取 P-non-linear
  av <- anova(fit)
  row_idx <- which(rownames(av) == exposure)
  p_nonlin <- av[row_idx + 1, "P"] 
  p_nonlin_text <- ifelse(p_nonlin < 0.001, "< 0.001", sprintf("%.3f", p_nonlin))
  
  # 获取绘图范围
  vals <- data[[exposure]]
  xlims <- range(vals, na.rm = TRUE)
  
  pred <- Predict(fit, name=exposure, 
                  seq(xlims[1], xlims[2], length.out = 100), 
                  fun=exp, ref.zero=TRUE)
  
  # 绘图
  p <- ggplot(pred) +
    geom_line(aes(x = .data[[exposure]], y = yhat), color = line_color, linewidth = 1.2) +
    geom_ribbon(aes(x = .data[[exposure]], ymin = lower, ymax = upper), 
                fill = line_color, alpha = 0.15) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    labs(title = label_name,
         subtitle = paste("P-non-linear:", p_nonlin_text),
         x = label_name, y = "Hazard Ratio (95% CI)") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5))
  
  return(p)
}


# =======================================================
# 4. 执行绘制 (共五张图)
# =======================================================

# A. 主队列 (FIB-4, NFS, APRI, AST/ALT) - 使用 cov_m4_fib
p1 <- draw_my_rcs_robust(df_rcs_fib, "fib4", "FIB-4 Score", cov_m4_fib, "duration_updated", "event_af_updated")
p2 <- draw_my_rcs_robust(df_rcs_fib, "nfs", "NFS Score", cov_m4_fib, "duration_updated", "event_af_updated")
p3 <- draw_my_rcs_robust(df_rcs_fib, "apri", "APRI Score", cov_m4_fib, "duration_updated", "event_af_updated")
p4 <- draw_my_rcs_robust(df_rcs_fib, "ast_alt_ratio", "AST/ALT Ratio", cov_m4_fib, "duration_updated", "event_af_updated")

# B. MRI队列 (Liver cT1) - 使用 cov_m4_mri，稍微换个颜色(如暗红)以作区分(如果不需要区分，直接删掉line_color参数)
p5 <- draw_my_rcs_robust(df_rcs_mri, "ct1_val", "Liver cT1 (ms)", cov_m4_mri, "duration_mri", "event_af_mri", line_color = "#992224")

# =======================================================
# 5. 合并五拼图与导出 (A4 横向尺寸) - 已修复 Layout Bug
# =======================================================

# 巧妙的自定义网格设计 (6列):
# 11 22 33  (前三张图各占2列，排满)
# # 44 55 # (# 代表官方的留白区域，4和5各占2列，刚好居中对齐)
custom_layout <- "
  112233
  #4455#
"

final_fig <- p1 + p2 + p3 + p4 + p5 + 
  plot_layout(design = custom_layout) +
  plot_annotation(tag_levels = 'A', 
                  title = "Association of Liver Fibrosis Scores and Liver cT1 with AF Risk",
                  subtitle = "Restricted Cubic Splines (Fully Adjusted Model 4, Values Winsorized)") &
  theme(plot.tag = element_text(face = 'bold', size = 16))

# 在控制台预览
print(final_fig)

# 导出为横向 A4 尺寸 (宽 11.69 英寸, 高 8.27 英寸)
ggsave("RCS_Combined_Model4_5Panels_A4.pdf", final_fig, 
       width = 11.69, height = 8.27, dpi = 300)

cat("导出成功！请查看当前目录下的 RCS_Combined_Model4_5Panels_A4.pdf \n")