library(mice)
library(survival)
library(broom)
library(dplyr)
library(stringr)
library(flextable)
library(officer)
library(data.table)
#imp_ready_for_cox<-mice_rediag

# 1. 提取 data_sup 的 Participant ID 并与 data_drink 合并


df_combined_1 <- fread("databcaa.csv", na.strings = c("", "NA"))

# 2. 输出合并后所有列的名字
names(df_combined_1)
library(dplyr)

# 假设 data 表格中的性别列名为 "Sex" 或 "sex_f"
# 1. 从 data 中提取 ID 和性别，然后合并到之前的 df_combined 中


library(dplyr)

library(dplyr)
library(tidyr)

df_risk <- df_combined_1  %>%
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

# --- 修正 NA 并准备最终 mids 对象 ---
long_sub_fixed <- complete(imp_ready_for_cox, "long", include = TRUE) %>%
  left_join(df_final %>% select(`Participant ID`, final_masld, metald_status), by = "Participant ID") %>%
  # 处理连接后的 NA 值
  mutate(
    final_masld = replace_na(final_masld, "No"),
    metald_status = replace_na( metald_status, "No")
  ) %>%
  group_by(.imp) %>%
  mutate(
    # 1. 人口学
    sub_age = factor(if_else(age_val < 65, "<65", ">=65"), levels = c("<65", ">=65")),
    sub_sex = factor(sex_f), # 假设原始 level 是 Female, Male
    
    # 2. 代谢与生活方式
    sub_bmi = factor(if_else(bmi_val < 40, "<40", ">=40"), levels = c("<40", ">=40")),
    sub_dm  = factor(t2dm_final, labels = c("No DM", "DM")),
    sub_alc_daily = factor(if_else(as.character(alc_freq_f) == "Daily or almost daily", "Daily", "Non-daily"), 
                           levels = c("Non-daily", "Daily")),
    
    # 3. 临床病因
    sub_masld = factor(if_else(final_masld == "Yes", "MASLD", "Non-MASLD"), 
                       levels = c("Non-MASLD", "MASLD")),
    
    sub_metald = factor(case_when(
      metald_status == "Yes"  ~ "Met-ALD",
      TRUE ~ "Non-MetALD"
    ), levels = c("Non-MetALD", "Met-ALD"))
  ) %>%
  ungroup()

imp_sub <- as.mids(long_sub_fixed)
#saveRDS(imp_sub, "imp_sub.rds")
# =======================================================
# 亚组分布检验 (基于第 1 个插补数据集)
# =======================================================

# 1. 提取代表性数据集
dat_check <- complete(imp_sub, 1)

# 2. 定义所有亚组变量
subgroup_vars <- c("sub_age", "sub_sex", "sub_bmi", "sub_dm", "sub_alc_daily", "sub_masld", "sub_metald")

# 3. 循环计算每个亚组的 Total 和 Cases
subgroup_dist <- lapply(subgroup_vars, function(var_name) {
  dat_check %>%
    group_by(Group = var_name, Level = get(var_name)) %>%
    summarise(
      Total = n(),
      Cases = sum(event_af_updated, na.rm = TRUE),
      Incidence_Rate = sprintf("%.2f%%", (Cases / Total) * 100)
    ) %>%
    ungroup()
}) %>% bind_rows()

# 4. 美化展示
cat("--- 各亚组样本量及房颤事件分布 (m=1) ---\n")
print(subgroup_dist, n = 40)

# 5. 预警提示：检查是否有极小样本量组
small_groups <- subgroup_dist %>% filter(Cases < 20)
if(nrow(small_groups) > 0) {
  warning("注意：以下亚组的事件数(Cases)过少，可能会导致回归模型不收敛或 CI 极宽：\n")
  print(small_groups)
}






library(mice)
library(survival)
library(broom)
library(dplyr)
library(stringr)
library(flextable)
library(officer)

# =======================================================
# 1. 准备工作：定义变量与协变量
# =======================================================

# 确保使用你已经生成的包含亚组信息的 mids 对象
# imp_ready <- imp_sub

# 定义 Model 4 协变量 (包含 PRS 和 PCs)
cov_m4_base <- c("age_val", "sex_f", "ethnicity_f",
                 "bmi_val", "smoke_f", "alc_freq_f", "tdi_val",
                 "htn_final", "t2dm_final", "lip_f", "ckd_final", "cvd_f",
                 "`Standard PRS for atrial fibrillation (AF)`",
                 paste0("`Genetic principal components | Array ", 1:10, "`"))

# 定义要分析的 4 个分类指标
cat_exposures <- c("fib4_cat", "nfs_cat", "apri_cat", "ratio_cat")

# 定义亚组变量
subgroups <- c("sub_age", "sub_sex", "sub_bmi", "sub_dm", "sub_alc_daily", "sub_masld", "sub_metald")

# 辅助函数：根据亚组变量自动剔除冲突协变量
get_adj_covs <- function(sub_var, all_covs) {
  exclude <- c()
  if (sub_var == "sub_age") exclude <- c("age_val")
  if (sub_var == "sub_sex") exclude <- c("sex_f")
  if (sub_var == "sub_bmi") exclude <- c("bmi_val")
  if (sub_var == "sub_dm") exclude <- c("t2dm_final")
  if (sub_var == "sub_alc_daily") exclude <- c("alc_freq_f")
  # Met-ALD 和 MASLD 通常不需要剔除基础协变量，除非共线性极其严重
  setdiff(all_covs, exclude)
}

# =======================================================
# 替换代码点 1: 修正后的分析函数
# =======================================================
analyze_subgroup_interaction <- function(imp_data, exposure, sub_var, covariates) {
  
  # --- Step 1: 交互作用 P 值计算 ---
  form_int_str <- paste("Surv(duration_updated, event_af_updated) ~", 
                        exposure, "*", sub_var, "+", 
                        paste(covariates, collapse = "+"))
  
  # 强制硬编码执行，避开 mice 环境变量坑
  cmd_int <- sprintf("with(imp_data, coxph(as.formula('%s')))", form_int_str)
  fit_int <- eval(parse(text = cmd_int))
  
  p_interaction <- summary(pool(fit_int)) %>%
    filter(str_detect(term, ":")) %>%
    slice(n()) %>%
    pull(p.value)
  if(length(p_interaction) == 0) p_interaction <- NA
  
  # --- Step 2: 分层计算 HR ---
  levels_sub <- levels(as.factor(imp_data$data[[sub_var]]))
  res_list <- list()
  
  form_sub_str <- paste("Surv(duration_updated, event_af_updated) ~", 
                        exposure, "+", paste(covariates, collapse = "+"))
  
  for(lvl in levels_sub) {
    # 核心修复：用 sprintf 强行把子组名称写成死代码传入 subset
    cmd_sub <- sprintf("with(imp_data, coxph(as.formula('%s'), subset = %s == '%s'))", 
                       form_sub_str, sub_var, lvl)
    fit_sub <- eval(parse(text = cmd_sub))
    
    summ <- summary(pool(fit_sub), conf.int = TRUE, exponentiate = TRUE)
    
    target_row <- summ %>%
      filter(str_detect(term, exposure)) %>%
      filter(str_detect(term, "High|Risk|2|3")) %>%
      slice(n())
    
    # 获取真实 HR
    if(nrow(target_row) == 0) {
      est <- NA; low <- NA; high <- NA; pval <- NA
    } else {
      est <- target_row$estimate; low <- target_row$conf.low; high <- target_row$conf.high; pval <- target_row$p.value
    }
    
    # 计算计数 (使用第一个插补集)
    dat_rep <- complete(imp_data, 1)
    n_tot <- sum(dat_rep[[sub_var]] == lvl, na.rm = TRUE)
    n_evt <- sum(dat_rep[[sub_var]] == lvl & dat_rep$event_af_updated == 1, na.rm = TRUE)
    
    res_list[[lvl]] <- data.frame(
      Subgroup_Var = sub_var, Level = lvl,
      Counts = paste0(n_tot, "/", n_evt),
      HR = est, LCI = low, UCI = high, P_val = pval, P_int = NA
    )
  }
  
  res_df <- bind_rows(res_list)
  res_df$P_int[1] <- p_interaction
  return(res_df)
}

# =======================================================
# 替换代码点 2: 修正后的第 3 部分批量执行循环
# =======================================================
final_tables <- list()

for (expo in cat_exposures) {
  message(paste("Analyzing Exposure:", expo))
  sub_res_list <- list()
  
  # --- 3.1 全人群 (All Patients) ---
  # 同样采用强行写入公式法，确保模型抓取真实的 expo 列
  form_all_str <- paste("Surv(duration_updated, event_af_updated) ~", expo, "+", paste(cov_m4_base, collapse = "+"))
  cmd_all <- sprintf("with(imp_sub, coxph(as.formula('%s')))", form_all_str)
  fit_all <- eval(parse(text = cmd_all))
  
  summ_all <- summary(pool(fit_all), conf.int = TRUE, exponentiate = TRUE)
  
  target_all <- summ_all %>% 
    filter(str_detect(term, expo) & str_detect(term, "High|Risk|2|3")) %>% 
    slice(n())
  
  dat_rep <- complete(imp_sub, 1)
  all_row <- data.frame(
    Subgroup_Var = "Overall", Level = "All Participants",
    Counts = paste0(nrow(dat_rep), "/", sum(dat_rep$event_af_updated, na.rm=TRUE)),
    HR = target_all$estimate, LCI = target_all$conf.low, UCI = target_all$conf.high,
    P_val = target_all$p.value, P_int = NA
  )
  sub_res_list[["All"]] <- all_row
  
  # --- 3.2 循环亚组 ---
  for (sub in subgroups) {
    adj_covs <- get_adj_covs(sub, cov_m4_base)
    sub_res <- analyze_subgroup_interaction(imp_sub, expo, sub, adj_covs)
    sub_res_list[[sub]] <- sub_res
  }
  
  # --- 3.3 合并并进行 FDR 校正 ---
  table_df <- bind_rows(sub_res_list) %>%
    mutate(Indicator = toupper(str_replace(expo, "_cat", "")))
  
  # 提取有效的 P_int 并做 BH(FDR) 校正
  valid_idx <- !is.na(table_df$P_int)
  table_df$FDR_P_int <- NA
  table_df$FDR_P_int[valid_idx] <- p.adjust(table_df$P_int[valid_idx], method = "BH")
  
  # 格式化输出
  table_df <- table_df %>%
    mutate(
      HR_CI = sprintf("%.2f (%.2f-%.2f)", HR, LCI, UCI),
      P_val_fmt = ifelse(P_val < 0.001, "<0.001", sprintf("%.3f", P_val)),
      P_int_fmt = ifelse(is.na(P_int), "", ifelse(P_int < 0.001, "<0.001", sprintf("%.3f", P_int))),
      FDR_P_int_fmt = ifelse(is.na(FDR_P_int), "", ifelse(FDR_P_int < 0.001, "<0.001", sprintf("%.3f", FDR_P_int)))
    ) %>%
    select(Indicator, Subgroup_Var, Level, Counts, HR_CI, P_val_fmt, P_int_fmt, FDR_P_int_fmt)
  
  final_tables[[expo]] <- table_df
}

#saveRDS(final_tables, "final_tables")
# =======================================================
# 4. 导出漂亮的 Word 表格 (带森林图列)
# =======================================================

library(ggplot2)
library(flextable)
library(dplyr)
library(stringr)

make_forest_flextable_with_gg <- function(df, title_str) {
  
  # --- 1. 数据预处理 (增强正则兼容性) ---
  df_plot <- df %>%
    mutate(
      HR_val = as.numeric(str_extract(HR_CI, "^[0-9.]+")),
      # 兼容不同类型的横杠 (短横线, 中横线)
      LCI_val = as.numeric(str_extract(HR_CI, "(?<=\\()[0-9.]+")),
      UCI_val = as.numeric(str_extract(HR_CI, "(?<=[\\-–])[0-9.]+(?=\\))")),
      
      Subgroup_Display = case_when(
        Subgroup_Var == "Overall" ~ "All Participants",
        Subgroup_Var == "sub_age" ~ "Age",
        Subgroup_Var == "sub_sex" ~ "Sex",
        Subgroup_Var == "sub_bmi" ~ "BMI Status",
        Subgroup_Var == "sub_dm"  ~ "Diabetes",
        Subgroup_Var == "sub_alc_daily" ~ "Alcohol Frequency",
        Subgroup_Var == "sub_masld" ~ "MASLD Status",
        Subgroup_Var == "sub_metald" ~ "Met-ALD Subgroup",
        TRUE ~ Subgroup_Var
      )
    ) %>%
    # 如果提取失败 (比如 Reference 组)，填充 1 避免 ggplot 崩溃
    mutate(across(c(HR_val, LCI_val, UCI_val), ~ifelse(is.na(.), 1, .)))
  
  # --- 2. 为每一行生成 ggplot 对象 ---
  x_limits <- c(0.4, 6.0) 
  
  # 关键修复 1: 使用 seq_along 确保纯净索引，并显式返回对象
  plots_list <- lapply(seq_len(nrow(df_plot)), function(i) {
    row_data <- df_plot[i, ]
    p_color <- ifelse(row_data$Subgroup_Var == "Overall", "black", "#2E5A88")
    
    p <- ggplot(row_data, aes(x = HR_val, y = 1)) +
      geom_vline(xintercept = 1, linetype = "dotted", color = "gray60", size = 0.4) +
      geom_errorbarh(aes(xmin = LCI_val, xmax = UCI_val), height = 0, color = p_color, size = 0.7) +
      geom_point(size = 2, shape = 18, color = p_color) +
      scale_x_log10(limits = x_limits, expand = c(0, 0)) +
      theme_void()
    return(p)
  })
  
  # 关键修复 2: 彻底清除列表的名称属性，这是导致 "list matrix to function" 报错的主因
  plots_list <- unname(plots_list)
  
  # --- 3. 生成底部坐标轴 ---
  axis_plot <- ggplot() +
    scale_x_log10(limits = x_limits, expand = c(0, 0), 
                  breaks = c(0.5, 1, 2, 4, 6), labels = c("0.5", "1", "2", "4", "6")) +
    theme_minimal() +
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(),
          axis.text.y = element_blank(), 
          axis.ticks.length.x = unit(.1, "cm"),
          axis.line.x = element_line(color = "black", size = 0.5),
          axis.text.x = element_text(size = 8))
  
  # --- 4. 构建 Flextable ---
  # 增加 FDR_P_int_fmt 列
  ft <- flextable(df_plot, col_keys = c("Subgroup_Display", "Level", "Counts", "Forest", "HR_CI", "P_val_fmt", "P_int_fmt", "FDR_P_int_fmt")) %>%
    flextable::compose(j = "Forest", value = as_paragraph(
      gg_chunk(value = plots_list, width = 1.3, height = 0.25)
    )) %>%
    add_footer_row(
      # 这里加了一个空的 FDR_P_int_fmt=""
      values = list(Subgroup_Display="", Level="", Counts="", Forest="axis", HR_CI="", P_val_fmt="", P_int_fmt="", FDR_P_int_fmt=""),
      colwidths = rep(1, 8) # 列数变成 8
    ) %>%
    flextable::compose(i = 1, j = "Forest", value = as_paragraph(
      gg_chunk(value = list(axis_plot), width = 1.3, height = 0.3)
    ), part = "footer") %>%
    set_header_labels(
      Subgroup_Display = "Subgroup", Level = "Level",
      Counts = "No. (Total/Cases)", Forest = "Forest Plot",
      HR_CI = "HR (95% CI)", P_val_fmt = "P Value", 
      P_int_fmt = "P for Interaction", FDR_P_int_fmt = "FDR-adjusted P" # 新增表头
    ) %>%
    merge_v(j = "Subgroup_Display") %>%
    theme_vanilla() %>%
    # 把 FDR 列加进居中对齐里
    align(j = c("Counts", "Forest", "HR_CI", "P_val_fmt", "P_int_fmt", "FDR_P_int_fmt"), align = "center", part = "all") %>%
    bold(i = ~ str_detect(P_val_fmt, "<0.05") | (as.numeric(P_val_fmt) < 0.05), j = "P_val_fmt") %>%
    # 也可以让 FDR P < 0.05 的加粗，引起注意
    bold(i = ~ str_detect(FDR_P_int_fmt, "<0.05") | (as.numeric(FDR_P_int_fmt) < 0.05), j = "FDR_P_int_fmt") %>%
    add_header_lines(values = title_str) %>%
    fontsize(size = 9, part = "all") %>%
    width(j = "Forest", width = 1.5) %>%
    autofit()
  
  return(ft)
}
# 导出到 Word
doc <- read_docx()

for (expo in cat_exposures) {
  df_res <- final_tables[[expo]]
  title <- paste0("Table: Subgroup Analysis of ", toupper(str_replace(expo, "_cat", "")), 
                  " (High Risk) on Incident AF (Model 4)")
  
  # 调用带森林图的新函数
  ft <- make_forest_flextable_with_gg(df_res, title)
  
  doc <- doc %>% 
    body_add_flextable(ft) %>%
    body_add_break() 
}

print(doc, target = "Subgroup_Analysis_Forest_Plots_Final.docx")