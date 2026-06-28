library(mice)
library(survival)
library(broom)
library(dplyr)
library(stringr)
library(flextable)
library(officer)
#imp_ready_for_cox<-mice_rediag
# =======================================================
# 1. 数据预处理：99% 缩尾 (Winsorization) - 保持不变
# =======================================================
cat("正在对 FIB-4 和 APRI 进行 99% 缩尾处理...\n")
prs_var <- "`Standard PRS for atrial fibrillation (AF)`"
pc_vars <- paste0("`Genetic principal components | Array ", 1:10, "`")

cov_m1 <- c("age_val", "sex_f", "ethnicity_f")
cov_m2 <- c(cov_m1, "bmi_val", "smoke_f", "alc_freq_f", "tdi_val")
cov_m3 <- c(cov_m2, "htn_final", "t2dm_final", "lip_f", "ckd_final", "cvd_f")
cov_m4 <- c(cov_m3, prs_var, pc_vars)


winsor_p99 <- function(x) {
  limit <- quantile(x, 0.99, na.rm = TRUE)
  x[x > limit] <- limit
  return(x)
}

library(mice)
library(survival)
library(broom)
library(dplyr)
library(stringr)
library(flextable)
library(officer)


# 假设 admin_end_date 已在环境里，或从数据中推断
# 定义竞争风险状态列
# =======================================================
# 1. 预处理：从原始表提取死因等必需变量并补全
# =======================================================
# =======================================================
# 1. 预处理：定义竞争风险状态 (修正版)
# =======================================================
extra_vars <- df_clinical_final %>%
  select(`Participant ID`, date_death, follow_up_end_updated) %>%
  distinct()

df_long_imp <- complete(imp_ready_for_cox, action = "long")

df_long_ready <- df_long_imp %>%
  left_join(extra_vars, by = "Participant ID") %>%
  # 关键：剔除随访时间异常的数据，防止 finegray 报错 y(0,-1)
  filter(!is.na(duration_updated) & duration_updated > 0) %>% 
  group_by(.imp) %>%
  mutate(across(c(fib4, apri, nfs, ast_alt_ratio), winsor_p99)) %>%
  ungroup() %>%
  mutate(
    # 第一步：定义数值
    status_num = case_when(
      event_af_updated == 1 ~ 1,
      event_af_updated == 0 & !is.na(date_death) & date_death <= follow_up_end_updated ~ 2,
      TRUE ~ 0
    ),
    # 第二步：转为因子（必须在 case_when 外面）
    status_cr_raw = factor(status_num, levels = c(0, 1, 2))
  )
# =======================================================
# 2. Fine-Gray 批量运行函数
# =======================================================
run_fine_gray_mice <- function(imp_long_df, exposure_name, covariates, model_label) {
  
  # 存储每个插补集的模型结果
  model_fits <- list()
  m <- max(imp_long_df$.imp)
  
  # 循环每个插补数据集 (Fine-Gray 必须对每个数据集单独变换)
  for(i in 1:m) {
    # 提取第 i 个数据集
    temp_dat <- imp_long_df %>% filter(.imp == i)
    
    # 准备暴露变量处理
    if(str_detect(exposure_name, "_cat")) {
      temp_dat$x_var <- factor(temp_dat[[exposure_name]])
      lvls <- levels(temp_dat$x_var)
      ref <- lvls[grep("Low|Normal|0", lvls, ignore.case = TRUE)[1]]
      if(!is.na(ref)) temp_dat$x_var <- relevel(temp_dat$x_var, ref = ref)
    } else {
      temp_dat$x_var <- as.numeric(scale(temp_dat[[exposure_name]]))
    }
    
    # --- 修正后的 Fine-Gray 变换 ---
    # 确保 status_cr_raw 是 factor
    # etype = 1 表示我们要关注的状态是 level '1' (房颤)
    fg_dat <- finegray(Surv(duration_updated, status_cr_raw) ~ ., 
                       data = temp_dat, 
                       etype = 1)
    
    # 此时 fg_dat 会包含新变量：fgstart, fgstop, fgstatus, fgwt
    form_str <- paste("Surv(fgstart, fgstop, fgstatus) ~ x_var +", 
                      paste(covariates, collapse = " + "))
    
    model_fits[[i]] <- coxph(as.formula(form_str), data = fg_dat, weights = fgwt)
    
    rm(temp_dat, fg_dat); if (i %% 2 == 0) gc()
    if (i %% 2 == 0) gc() # 每跑两个数据集进行一次深度垃圾回收
  }
  
  # 使用 mice::pool 汇总 m 个模型
  pooled <- pool(as.mira(model_fits))
  
  # --- 提取结果 ---
  res <- summary(pooled, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(str_starts(term, "x_var")) %>%
    mutate(
      Model = model_label,
      Exposure = toupper(str_replace(exposure_name, "_cat", "")),
      Group = if_else(term == "x_var", "Per SD increase", str_replace(term, "x_var", "")),
      HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
      P_val = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
    )
  
  # --- 计算计数 (基于第一个数据集) ---
  dat_rep <- imp_long_df %>% filter(.imp == 1)
  if(str_detect(exposure_name, "_cat")) {
    counts_df <- dat_rep %>%
      group_by(Group_Raw = as.character(.data[[exposure_name]])) %>%
      summarise(Counts = paste0(n(), "/", sum(status_cr_raw == 1)))
    res <- res %>% left_join(counts_df, by = c("Group" = "Group_Raw"))
  } else {
    res$Counts <- paste0(nrow(dat_rep), "/", sum(dat_rep$status_cr_raw == 1))
  }
  
  return(res)
}

# =======================================================
# 3. 执行分析
# =======================================================
# 模型定义保持不变
model_list <- list("Model 1" = cov_m1, "Model 2" = cov_m2, "Model 3" = cov_m3, "Model 4" = cov_m4)
exposures <- c("fib4", "fib4_cat", "nfs", "nfs_cat", "apri", "apri_cat", "ast_alt_ratio", "ratio_cat")

all_fg_results <- list()

for (expo in exposures) {
  for (mod_name in names(model_list)) {
    message(sprintf("Running Fine-Gray: %s | %s", expo, mod_name))
    all_fg_results[[paste(expo, mod_name)]] <- run_fine_gray_mice(df_long_ready, expo, model_list[[mod_name]], mod_name)
  }
}

# =======================================================
# 4. 结果汇总与列筛选 (SELECT)
# =======================================================

# 将列表合并为数据框
final_fg_table_df <- bind_rows(all_fg_results) %>%
  # 使用 select 挑选分析需要的核心列并重命名部分列
  select(
    Model, 
    Indicator = Exposure, 
    Group, 
    `No. (Total/Cases)` = Counts, 
    `sHR (95% CI)` = HR_CI,   # 注意：Fine-Gray 输出的是 subdistribution HR
    P_val
  )

# =======================================================
# 5. 生成最终表格 (Flextable)
# =======================================================

# 创建美化表格
ft_fg <- flextable(final_fg_table_df) %>%
  theme_booktabs() %>%
  # 合并同类项，让表格更整洁
  merge_v(j = c("Model", "Indicator")) %>% 
  # 设置表头
  set_header_labels(
    Model = "Model",
    Indicator = "Liver Fibrosis Index",
    Group = "Level / Change",
    `No. (Total/Cases)` = "No. (Total/Cases)",
    `sHR (95% CI)` = "sHR (95% CI)",
    P_val = "P value"
  ) %>%
  # 格式美化
  autofit() %>%
  align(j = c("No. (Total/Cases)", "sHR (95% CI)", "P_val"), align = "center", part = "all") %>%
  # P值加粗逻辑 (显著性 < 0.05)
  bold(i = ~ !str_detect(P_val, ">") & as.numeric(str_replace_all(P_val, "<", "")) < 0.05, j = "P_val") %>%
  # 给表格加个注脚解释 sHR
  add_footer_lines("sHR: Subdistribution Hazard Ratio calculated by Fine-Gray model (competing risk: Death).")

# =======================================================
# 6. 导出到 Word
# =======================================================

doc_fg <- read_docx() %>%
  body_add_par("Table: Fine-Gray Subdistribution Hazard Models for Incident AF", style = "heading 1") %>%
  body_add_par("This table accounts for the competing risk of death before AF diagnosis.", style = "Normal") %>%
  body_add_flextable(ft_fg)

print(doc_fg, target = "Liver_AF_FineGray_Results.docx")

cat("🎉 Fine-Gray 分析结果已成功导出至 Liver_AF_FineGray_Results.docx\n")
