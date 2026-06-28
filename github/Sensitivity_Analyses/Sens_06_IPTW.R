library(mice)
library(WeightIt)
library(survey)
library(broom)
library(cobalt)
library(dplyr)
library(purrr)
library(stringr)
library(flextable)
library(officer)
#imp_ready_for_cox<-mice_rediag
run_full_iptw_pipeline_v4 <- function(mids_obj, exposure_name) {
  cat("\n>>> [IPTW Analysis] 正在处理:", exposure_name, "<<<\n")
  
  dat_long <- complete(mids_obj, action = "long")
  dat_long$treatment <- as.factor(dat_long[[exposure_name]])
  
  diag_list <- list()
  fit_list <- list()
  imp_groups <- dat_long %>% group_by(.imp) %>% group_split()
  
  for (i in seq_along(imp_groups)) {
    d <- imp_groups[[i]]
    
    formula_str <- paste0(
      "treatment ~ poly(age_val, 2) + sex_f + ethnicity_f + tdi_val + ",
      "poly(bmi_val, 2) + smoke_f + alc_freq_f + htn_final + t2dm_final + ",
      "ckd_final + cvd_f + lip_f + age_val:bmi_val + age_val:sex_f"
    )
    
    w <- weightit(as.formula(formula_str), data = d, method = "ps", estimand = "ATO")
    
    # --- 修正后的变量提取逻辑 ---
    # 1. 提取 ESS
    ess_summary <- summary(w)$effective.sample.size
    adj_ess <- ess_summary[nrow(ess_summary), ] # 提取 Adjusted 行
    
    # 2. 映射 APRI 等二分类指标的名字
    if (length(levels(d$treatment)) == 2) {
      names(adj_ess) <- levels(d$treatment)
    }
    
    
    
    # --- 抓取 Max ASMD ---
    b <- bal.tab(w, stats = "m", binary = "std", abs = TRUE)
    max_asmd <- if (!is.null(b[["Balance.Across.Pairs"]])) {
      max(b$Balance.Across.Pairs$Max.Diff.Adj, na.rm = TRUE)
    } else {
      max(b$Balance$Diff.Adj, na.rm = TRUE)
    }
    
    # 3. 统计原始计数 (先统计，再 mutate)
    counts_raw <- d %>%
      group_by(treatment) %>%
      summarise(N_raw = n(), Events_raw = sum(event_af_updated == 1, na.rm = TRUE), .groups = "drop")
    
    # 4. 合并诊断信息 (确保这里使用的是 adj_ess)
    diag_list[[i]] <- counts_raw %>%
      mutate(
        .imp = d$.imp[1], 
        Max_ASMD = max_asmd,
        # 修正：将 ess_vals 改为 adj_ess
        ESS = as.numeric(adj_ess[as.character(treatment)]) 
      )
    # --- 计算加权 Cox ---
    d$weights <- w$weights
    # 针对极端权重进行截断 (解决 NFS 警告)
    q <- quantile(d$weights, c(0.01, 0.99))
    d$weights <- pmin(pmax(d$weights, q[1]), q[2])
    
    des <- svydesign(ids = ~1, weights = ~weights, data = d)
    pc_vars <- paste0("`Genetic principal components | Array ", 1:10, "`", collapse = " + ")
    cox_formula <- paste0("Surv(duration_updated, event_af_updated) ~ treatment + ",
                          "`Standard PRS for atrial fibrillation (AF)` + ", pc_vars)
    
    fit_list[[i]] <- svycoxph(as.formula(cox_formula), design = des)
  }
  
  # 池化与合并
  res_model <- summary(pool(fit_list), conf.int = TRUE, exponentiate = TRUE)
  res_diag <- bind_rows(diag_list) %>%
    group_by(treatment) %>%
    summarise(Mean_N = mean(N_raw), Mean_E = mean(Events_raw), 
              Mean_ESS = mean(ESS), Mean_ASMD = mean(Max_ASMD), .groups = "drop")
  
  return(list(model = res_model, diag = res_diag))
}
# =======================================================
# 2. 增强版清理函数：手动添加 Reference 行
# =======================================================
clean_iptw_res_v4 <- function(res_obj, label) {
  diag_df <- res_obj$diag
  model_df <- res_obj$model
  
  # 1. 确定 Reference Level (通常是第一行)
  ref_level <- as.character(diag_df$treatment[1])
  
  # 2. 清理模型输出
  model_clean <- model_df %>%
    filter(str_detect(term, "treatment")) %>%
    mutate(
      Level = str_replace(term, "treatment", ""),
      HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
      P = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
    )
  
  # 3. 拼装最终行
  diag_df %>%
    mutate(Level = as.character(treatment)) %>%
    left_join(model_clean, by = "Level") %>%
    mutate(
      Exposure = label,
      `No. (Events/Total)` = sprintf("%.0f/%.0f", Mean_E, Mean_N),
      `ESS (Weighted)` = sprintf("%.0f", Mean_ESS),
      `Max ASMD` = sprintf("%.3f", Mean_ASMD),
      `HR (95% CI)` = ifelse(Level == ref_level, "1.00 (Reference)", HR_CI),
      `P value` = ifelse(Level == ref_level, "-", P)
    ) %>%
    select(Exposure, Level, `No. (Events/Total)`, `ESS (Weighted)`, `Max ASMD`, `HR (95% CI)`, `P value`)
}

# 运行分析
out_fib4 <- run_full_iptw_pipeline_v4(imp_ready_for_cox, "fib4_cat")
out_nfs  <- run_full_iptw_pipeline_v4(imp_ready_for_cox, "nfs_cat")
out_apri <- run_full_iptw_pipeline_v4(imp_ready_for_cox, "apri_cat")
out_ratio <- run_full_iptw_pipeline_v4(imp_ready_for_cox, "ratio_cat")
# 合并结果表格
final_iptw_table <- bind_rows(
  clean_iptw_res_v4(out_fib4, "FIB-4"),
  clean_iptw_res_v4(out_nfs, "NFS"),
  clean_iptw_res_v4(out_apri, "APRI"),
  clean_iptw_res_v4(out_ratio, "AST_ALT")
)

# Flextable 导出
ft <- flextable(final_iptw_table) %>%
  merge_v(j = "Exposure") %>%
  theme_booktabs() %>%
  autofit() %>%
  align(align = "center", part = "all") %>%
  bold(part = "header") %>%
  add_footer_lines("ATO: Overlap weights; ESS: Effective Sample Size; ASMD: Absolute Standardized Mean Difference.")

# 保存到 Word
doc <- read_docx() %>%
  body_add_par("Table: Sensitivity Analysis using IPTW with ATO Weights", style = "heading 1") %>%
  body_add_flextable(ft)
print(doc, target = "Final_Diagnostic_IPTW.docx")




library(cobalt)
library(ggplot2)

# 1. 定义需要绘图的指标及其展示标签
exposure_map <- list(
  "fib4_cat"  = "FIB-4",
  "nfs_cat"   = "NFS",
  "apri_cat"  = "APRI",
  "ratio_cat" = "AST/ALT Ratio"
)

# 2. 遍历生成图片
for (expo_name in names(exposure_map)) {
  
  cat("\n>>> 正在生成 Love Plot:", exposure_map[[expo_name]], "<<<\n")
  
  # 提取第一个插补集进行诊断 (IPTW 诊断惯例)
  d_plot <- complete(imp_ready_for_cox, 1)
  d_plot$treatment <- as.factor(d_plot[[expo_name]])
  
  # 构造与 v4 pipeline 完全一致的公式 (包含多项式和交互项)
  formula_plot <- as.formula(paste0(
    "treatment ~ poly(age_val, 2) + sex_f + ethnicity_f + tdi_val + ",
    "poly(bmi_val, 2) + smoke_f + alc_freq_f + htn_final + t2dm_final + ",
    "ckd_final + cvd_f + lip_f + age_val:bmi_val + age_val:sex_f"
  ))
  
  # 计算 ATO 权重 (与分析端逻辑对齐)
  w_out <- weightit(formula_plot, 
                    data = d_plot, 
                    method = "ps", 
                    estimand = "ATO")
  
  # 准备绘图对象
  # 使用 which.treat = .all 确保多分类也能计算
  b_tab <- bal.tab(w_out, un = TRUE, stats = c("m"), abs = TRUE)
  
  # 3. 绘制 Love Plot
  p <- love.plot(
    b_tab,
    threshold = 0.1,        # 统计学公认平衡线
    abs = TRUE,             # 显示绝对差值
    agg.fun = "max",        # 多分类下展示最大偏差 (最严谨)
    var.order = "unadjusted", # 按原始偏差排序，方便观察改善程度
    line = TRUE,            # 连线方便观察
    stars = "std",          # 标注标准化差值
    colors = c("#E41A1C", "#377EB8"), # 红色未调整，蓝色 ATO
    shapes = c(21, 19),
    title = paste("Covariate Balance (ATO):", exposure_map[[expo_name]]),
    sample.names = c("Unadjusted", "Balanced (ATO)")
  ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  # 4. 保存文件
  file_name <- paste0("LovePlot_", expo_name, ".png")
  ggsave(filename = file_name, plot = p, width = 8, height = 10, dpi = 300)
  
  cat("已保存:", file_name, "\n")
}
library(WeightIt)
library(survey)
library(cobalt)
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
run_mri_iptw_pipeline <- function(mids_obj, exposure_name) {
  cat("\n>>> [MRI IPTW Analysis] 正在处理:", exposure_name, "<<<\n")
  
  dat_long <- complete(mids_obj, action = "long")
  dat_long$treatment <- as.factor(dat_long[[exposure_name]])
  
  diag_list <- list()
  fit_list <- list()
  imp_groups <- dat_long %>% group_by(.imp) %>% group_split()
  
  for (i in seq_along(imp_groups)) {
    d <- imp_groups[[i]]
    
    ps_formula <- as.formula(paste0(
      "treatment ~ poly(age_mri, 2) + sex_f + ethnicity_f + tdi_val + ",
      "poly(bmi_i2, 2) + pdff_z + htn_mri_f + t2dm_mri_f + ",
      "`Standard PRS for atrial fibrillation (AF)`"
    ))
    
    # 核心：计算 ATO 权重
    w <- weightit(ps_formula, data = d, method = "ps", estimand = "ATO")
    
    # --- 完全对齐主队列 APRI 的 ESS 提取逻辑 ---
    ess_summary <- summary(w)$effective.sample.size
    adj_ess <- ess_summary[nrow(ess_summary), ] # 提取 Adjusted 行
    
    # 关键点：对于二分类指标，WeightIt 默认列名是 Control/Treated
    # 我们强制将其映射为数据真实的 Levels (Low/High)
    if (length(levels(d$treatment)) == 2) {
      names(adj_ess) <- levels(d$treatment)
    }
    
    # --- 抓取 Max ASMD ---
    b <- bal.tab(w, stats = "m", binary = "std", abs = TRUE)
    max_asmd <- max(b$Balance$Diff.Adj, na.rm = TRUE)
    
    # 统计原始计数
    counts_raw <- d %>%
      group_by(treatment) %>%
      summarise(N_raw = n(), 
                Events_raw = sum(event_af_mri == 1, na.rm = TRUE), 
                .groups = "drop")
    
    # 合并诊断信息
    diag_list[[i]] <- counts_raw %>%
      mutate(
        .imp = d$.imp[1], 
        Max_ASMD = max_asmd,
        # 此时 adj_ess 的 Names 已经是 Low/High，匹配成功
        ESS = as.numeric(adj_ess[as.character(treatment)]) 
      )
    
    # --- 计算加权 Cox ---
    d$weights <- w$weights
    q <- quantile(d$weights, c(0.005, 0.995))
    d$weights <- pmin(pmax(d$weights, q[1]), q[2])
    
    des <- svydesign(ids = ~1, weights = ~weights, data = d)
    pc_vars <- paste0("`Genetic principal components | Array ", 1:3, "`", collapse = " + ")
    cox_formula <- as.formula(paste0("Surv(duration_mri, event_af_mri) ~ treatment + ", pc_vars))
    
    fit_list[[i]] <- svycoxph(cox_formula, design = des)
  }
  
  # 池化与合并
  res_model <- summary(pool(fit_list), conf.int = TRUE, exponentiate = TRUE)
  res_diag <- bind_rows(diag_list) %>%
    group_by(treatment) %>%
    summarise(Mean_N = mean(N_raw), Mean_E = mean(Events_raw), 
              Mean_ESS = mean(ESS), Mean_ASMD = mean(Max_ASMD), .groups = "drop")
  
  return(list(model = res_model, diag = res_diag))
}
# 运行子队列分析 (cT1 分类变量)
out_ct1_iptw <- run_mri_iptw_pipeline(imp_final, "ct1_f")

clean_iptw_res_v4 <- function(res_obj, label) {
  diag_df <- res_obj$diag
  model_df <- res_obj$model
  
  # 1. 确定 Reference Level (通常是第一行)
  ref_level <- as.character(diag_df$treatment[1])
  
  # 2. 清理模型输出
  model_clean <- model_df %>%
    filter(str_detect(term, "treatment")) %>%
    mutate(
      Level = str_replace(term, "treatment", ""),
      HR_CI = sprintf("%.2f (%.2f-%.2f)", estimate, conf.low, conf.high),
      P = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
    )
  
  # 3. 拼装最终行
  diag_df %>%
    mutate(Level = as.character(treatment)) %>%
    left_join(model_clean, by = "Level") %>%
    mutate(
      Exposure = label,
      `No. (Events/Total)` = sprintf("%.0f/%.0f", Mean_E, Mean_N),
      `ESS (Weighted)` = sprintf("%.0f", Mean_ESS),
      `Max ASMD` = sprintf("%.3f", Mean_ASMD),
      `HR (95% CI)` = ifelse(Level == ref_level, "1.00 (Reference)", HR_CI),
      `P value` = ifelse(Level == ref_level, "-", P)
    ) %>%
    select(Exposure, Level, `No. (Events/Total)`, `ESS (Weighted)`, `Max ASMD`, `HR (95% CI)`, `P value`)
}

# 使用你主队列的一致清理逻辑
final_mri_iptw_table <- clean_iptw_res_v4(out_ct1_iptw, "Liver cT1 (Categorical)")

# Flextable 导出
ft_iptw_mri <- flextable(final_mri_iptw_table) %>%
  theme_booktabs() %>% autofit() %>%
  align(align = "center", part = "all") %>%
  add_footer_lines("Weighted using ATO (Overlap Weights). Covariates adjusted: age, sex, ethnicity, TDI, BMI, PDFF, HTN, T2DM, PRS, and PC1-3.")

doc_iptw <- read_docx() %>%
  body_add_par("Table: Sensitivity Analysis of Liver cT1 using IPTW (MRI Sub-cohort)", style = "heading 1") %>%
  body_add_flextable(ft_iptw_mri)
print(doc_iptw, target = "MRI_Subcohort_IPTW_Results.docx")
# 提取一个插补集进行 Love Plot 演示
d_plot <- complete(imp_final, 1)
d_plot$treatment <- as.factor(d_plot$ct1_f)

# 重复计算一次权重用于绘图
w_plot <- weightit(
  treatment ~ poly(age_mri, 2) + sex_f + ethnicity_f + tdi_val + 
    poly(bmi_i2, 2) + pdff_z + htn_mri_f + t2dm_mri_f + 
    `Standard PRS for atrial fibrillation (AF)`,
  data = d_plot, method = "ps", estimand = "ATO"
)

# --- 1. 定义变量名映射表 (让图表更学术) ---
var_map <- data.frame(
  old = c("poly(age_mri, 2)1", "poly(age_mri, 2)2", "sex_f_Male", 
          "ethnicity_f_White", "tdi_val", "poly(bmi_i2, 2)1", 
          "poly(bmi_i2, 2)2", "pdff_z", "htn_mri_f_1", "t2dm_mri_f_1",
          "Standard PRS for atrial fibrillation (AF)"),
  new = c("Age (Linear)", "Age (Quadratic)", "Sex (Male)", 
          "Ethnicity (White)", "TDI", "BMI (Linear)", 
          "BMI (Quadratic)", "Liver Fat (PDFF, log-Z)", "Hypertension", "T2DM",
          "Genetic Risk (PRS)")
)

# --- 2. 重新绘图 ---
p_love <- love.plot(
  w_plot, 
  threshold = 0.1, 
  abs = TRUE,
  # 核心修正：强制所有变量（包含二分类）使用标准化差值
  binary = "std", 
  # 使用 stars 标注标准化变量（消除警告的关键）
  stars = "std",
  var.order = "unadjusted", 
  line = TRUE,
  # 替换变量名
  var.names = var_map,
  colors = c("#E41A1C", "#377EB8"), 
  shapes = c(21, 19),
  title = "Covariate Balance: Liver cT1 High vs Low",
  sample.names = c("Original", "Weighted (ATO)")
) + 
  theme_bw() + 
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title.x = element_text(face = "bold")
  ) +
  xlab("Standardized Mean Differences") # 明确标注 X 轴

ggsave("LovePlot_MRI_cT1_v2.pdf", plot = p_love, width = 8, height = 9, dpi = 300)