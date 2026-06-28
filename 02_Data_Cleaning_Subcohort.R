library(tidyverse)
library(data.table)
library(mice)
library(survival)
# --- A. ICD 10/9 正则模式 ---
# 注意：CVD 比较复杂，包含冠心病、心衰、卒中、外周血管、瓣膜病
# 1. 房颤 (AF/AFL) - 用于排除基线和定义结局
pat_af_icd <- "^I48" # ICD-10: I48.x; ICD-9: 4273
pat_htn_icd    <- "^I10|^I11|^I12|^I13|^I15|^40[1-5]"
pat_lipids_icd <- "^E78|^272"
pat_ckd_icd    <- "^N18|^585"
pat_cvd_icd    <- paste0(
  "^I2[0-5]|^I50|",              # 冠心病 & 心衰 (ICD-10)
  "^I0[5-8]|^I3[4-9]|",          # 瓣膜病
  "^I7[0-3]|",                   # 外周血管病
  "^I6[0-4]|",                   # 卒中
  "^41[0-4]|^428|",              # 冠心病 & 心衰 (ICD-9)
  "^39[4-6]|^424|^44[0-3]|",     # 瓣膜病 & 外周
  "^43[01346]"                   # 卒中
)
pat_t1dm_icd   <- "^E10|^250\\.?[13579]" 
pat_t2dm_icd   <- "^E11|^250\\.?[02468]"
pat_gdm_icd    <- "^O24[489]|^6488"

# --- B. 自述疾病编码 (Field 20002 Codes) ---
# 建议直接使用数值编码，比文本匹配更准。
# 1065=HTN, 1473=High chol, 1220/1223=Diabetes, 1075=Heart attack, 1081=Stroke, 1079=Angina
codes_sr_htn  <- c("1065", "hypertension", "essential hypertension")
codes_sr_lip  <- c("1473", "high cholesterol")
codes_sr_ckd  <- c("1192", "1193", "1194", "renal failure", "kidney failure")
codes_sr_cvd  <- c("1075", "1074", "1076", "1077", "1078", "1079", "1081", "1082", "1086", "1471", "1496", "1583", # 包含心梗, 心绞痛, 心衰, 卒中
                   "heart attack", "myocardial infarction", "angina", "stroke", "heart failure") 
codes_sr_t2dm <- c("1220", "1223", "diabetes", "type 2 diabetes")
codes_sr_t1dm <- c("1222", "type 1 diabetes")

# --- C. 药物正则 (Instance 0) ---
pat_htn_meds <- "atenolol|bisoprolol|metoprolol|propranolol|amlodipine|nifedipine|felodipine|ramipril|lisinopril|enalapril|losartan|candesartan|valsartan|bendroflumethiazide|furosemide|spironolactone"
pat_lipid_meds <- "atorvastatin|simvastatin|rosuvastatin|pravastatin|fluvastatin|beclofibrate|fenofibrate|ezetimibe"
pat_t2dm_meds <- "metformin|gliclazide|glimepiride|glipizide|tolbutamide|pioglitazone|rosiglitazone|sitagliptin|vildagliptin|saxagliptin|linagliptin|dapagliflozin|empagliflozin|canagliflozin|acarbose|repaglinide|nateglinide|exenatide|liraglutide|semaglutide"




# ==============================================================================
# 第一部分：MRI 子队列构建与基线平移 (Time-shift)
# ==============================================================================

# 1. 提取有 MRI 数据的人群，并确定每个人的 MRI 检查日期
df_mri_raw <- df_clinical_final %>% 
  filter(!is.na(`Liver iron corrected T1 (ct1) | Instance 2`)) %>%
  mutate(date_mri = as.Date(`Date of attending assessment centre | Instance 2`))

# 2. 核心：利用 ICD 全记录表，抓取 MRI 之前的共病 (解决 ICD10 日期在 I2 之前的问题)
# 假设 df_all_icd_long 是包含所有 ICD 记录（不限日期）的长表
mri_date_lookup <- df_mri_raw %>% select(`Participant ID`, date_mri)

long_records_mri <- df_all_icd_long %>% 
  inner_join(mri_date_lookup, by = "Participant ID") %>%
  filter(event_date <= date_mri) # 关键：截止到 MRI 当天

# 提取各类疾病在 MRI 之前的 ICD 确诊 ID 列表
ids_htn_mri_icd  <- long_records_mri %>% filter(str_detect(code, pat_htn_icd)) %>% pull(`Participant ID`) %>% unique()
ids_lip_mri_icd  <- long_records_mri %>% filter(str_detect(code, pat_lipids_icd)) %>% pull(`Participant ID`) %>% unique()
ids_t2dm_mri_icd <- long_records_mri %>% filter(str_detect(code, pat_t2dm_icd)) %>% pull(`Participant ID`) %>% unique()
ids_ckd_mri_icd  <- long_records_mri %>% filter(str_detect(code, pat_ckd_icd)) %>% pull(`Participant ID`) %>% unique()
ids_cvd_mri_icd  <- long_records_mri %>% filter(str_detect(code, pat_cvd_icd)) %>% pull(`Participant ID`) %>% unique()
ids_af_mri_icd   <- long_records_mri %>% filter(str_detect(code, pat_af_icd)) %>% pull(`Participant ID`) %>% unique()

# ==============================================================================
# 第二部分：全面定义 Instance 2 状态的共病 (含 Self-report & CKD/高脂)
# ==============================================================================

df_mri_vars_full <- df_mri_raw %>%
  mutate(
    # --- 1. 房颤排除逻辑 (Prevalent AF at MRI) ---
    sr_af_inst2 = if_any(contains("Non-cancer illness code, self-reported | Instance 2"), 
                         ~as.character(.) %in% c("atrial fibrillation", "atrial flutter")),
    is_prevalent_af_mri = if_else((`Participant ID` %in% ids_af_mri_icd) | sr_af_inst2, 1, 0, missing = 0),
    
    # --- 2. 整合 Instance 2 的药物与自述 (Self-report Instance 2) ---
    # 药物
    med_htn_mri  = if_any(contains("Treatment/medication code | Instance 2"), ~str_detect(tolower(as.character(.)), pat_htn_meds) %in% TRUE),
    med_lip_mri  = if_any(contains("Treatment/medication code | Instance 2"), ~str_detect(tolower(as.character(.)), pat_lipid_meds) %in% TRUE),
    med_t2dm_mri = if_any(contains("Treatment/medication code | Instance 2"), ~str_detect(tolower(as.character(.)), pat_t2dm_meds) %in% TRUE),
    # 自述
    sr_htn_mri  = if_any(contains("Non-cancer illness code, self-reported | Instance 2"), ~as.character(.) %in% codes_sr_htn),
    sr_lip_mri  = if_any(contains("Non-cancer illness code, self-reported | Instance 2"), ~as.character(.) %in% codes_sr_lip),
    sr_t2dm_mri = if_any(contains("Non-cancer illness code, self-reported | Instance 2"), ~as.character(.) %in% codes_sr_t2dm),
    sr_ckd_mri  = if_any(contains("Non-cancer illness code, self-reported | Instance 2"), ~as.character(.) %in% codes_sr_ckd),
    sr_cvd_mri  = if_any(contains("Non-cancer illness code, self-reported | Instance 2"), ~as.character(.) %in% codes_sr_cvd),
    
    # --- 3. 终极共病定义 (Logic Assembly) ---
    # 高血压 (ICD OR 自述 OR 药物 OR MRI当天实测高)
    sbp_mri_val = rowMeans(select(., matches("Systolic blood pressure.*Instance 2")), na.rm = TRUE),
    dbp_mri_val = rowMeans(select(., matches("Diastolic blood pressure.*Instance 2")), na.rm = TRUE),
    htn_mri_f = if_else((`Participant ID` %in% ids_htn_mri_icd) | sr_htn_mri | med_htn_mri | (sbp_mri_val >= 140 | dbp_mri_val >= 90), 1, 0, missing = 0),
    
    # 高脂血症 (ICD OR 自述 OR 药物 OR 继承基线异常)
    lip_mri_f = if_else((`Participant ID` %in% ids_lip_mri_icd) | sr_lip_mri | med_lip_mri | baseline_lipids == 1, 1, 0, missing = 0),
    
    # 2型糖尿病 (ICD OR 自述 OR 药物 OR 继承基线确诊)
    t2dm_mri_f = if_else((`Participant ID` %in% ids_t2dm_mri_icd) | sr_t2dm_mri | med_t2dm_mri | baseline_t2dm == 1, 1, 0, missing = 0),
    
    # 慢性肾病 CKD (ICD OR 自述 OR 继承基线eGFR<60)
    ckd_mri_f = if_else((`Participant ID` %in% ids_ckd_mri_icd) | sr_ckd_mri | baseline_ckd == 1, 1, 0, missing = 0),
    
    # 心血管疾病 CVD (ICD OR 自述 OR 继承基线)
    cvd_mri_f = if_else((`Participant ID` %in% ids_cvd_mri_icd) | sr_cvd_mri | baseline_cvd == 1, 1, 0, missing = 0)
  ) %>%
  # 排除 MRI 之前的房颤流行病例
  filter(is_prevalent_af_mri == 0) %>%
  # 更新随访终点
  mutate(
    duration_mri = as.numeric(follow_up_end_updated - date_mri) / 365.25,
    event_af_mri = if_else(event_af_updated == 1 & af_date > date_mri, 1, 0, missing = 0)
  ) %>%
  filter(duration_mri > 0)
#saveRDS(df_mri_vars_full, "df_mri_vars_full.rds")

covariates_to_check_i2 <- c(
  "Townsend deprivation index at recruitment",
  "Smoking status | Instance 2",
  "Systolic blood pressure, automated reading | Instance 2 | Array 0",
  "Systolic blood pressure, automated reading | Instance 2 | Array 1",
  "Diastolic blood pressure, automated reading | Instance 2 | Array 0",
  "Diastolic blood pressure, automated reading | Instance 2 | Array 1",
  "Age at recruitment",
  "Alcohol drinker status | Instance 2",
  "Alcohol intake frequency. | Instance 2",
  "Sex",
  "Waist circumference | Instance 2",
  "Proton density fat fraction (PDFF) | Instance 2"
)

# 2. 检查缺失情况
missing_summary <- df_mri_vars_full %>%
  select(all_of(covariates_to_check_i2)) %>%
  summarise(across(everything(), list(
    na_count = ~sum(is.na(.)),
    na_percent = ~round(sum(is.na(.)) / n() * 100, 2)
  ))) %>%
  pivot_longer(everything(), names_to = c("variable", "stat"), names_sep = "_(?!.*_)") %>% # 使用正则匹配最后一个下划线
  pivot_wider(names_from = stat, values_from = value)

print(missing_summary)
# 1. 定义 cT1 分组
df_mri_vars_full <- df_mri_vars_full %>%
  mutate(ct1_group = case_when(`Liver iron corrected T1 (ct1) | Instance 2` >= 750 ~ "High Risk", 
                               `Liver iron corrected T1 (ct1) | Instance 2` < 750 ~ "Low Risk"
                               ),
                               ct1_group = factor(ct1_group, levels = c("Low Risk", "High Risk")))


# 2. 计算各组事件数和总随访时间
ir_table <- df_mri_vars_full %>%
  group_by(ct1_group) %>%
  summarise(
    n = n(),
    events = sum(event_af_mri),
    total_person_years = sum(duration_mri),
    .groups = 'drop'
  ) %>%
  mutate(
    ir_1000py = (events / total_person_years) * 1000,
    lower_ci = (qchisq(0.05/2, 2*events)/2) / total_person_years * 1000, # 泊松分布置信区间
    upper_ci = (qchisq(1-0.05/2, 2*events+2)/2) / total_person_years * 1000
  )

print(ir_table)
#saveRDS(df_mri_vars_full, "df_mri_vars_full.rds")
# ==============================================================================
# 第三部分：多重插补 (MICE) 准备与执行
# ==============================================================================

# 选择插补变量 (Instance 2 测值为优先)
df_mi_mri <- df_mri_vars_full %>%
  mutate(ct1_f=factor(ct1_group),
    sex_f = factor(Sex),
    smoke_i2_f = factor(`Smoking status | Instance 2`),
    alc_stati2_f = factor(`Alcohol drinker status | Instance 2`),
    alc_freqi2_f = factor(`Alcohol intake frequency. | Instance 2`),
    ethnicity_f = case_when(
      str_detect(as.character(`Ethnic background | Instance 0`), "White|British|Irish") ~ "White",
      str_detect(as.character(`Ethnic background | Instance 0`), "Asian|Indian|Pakistani|Chinese|Bangladesh") ~ "Asian",
      str_detect(as.character(`Ethnic background | Instance 0`), "Black|African|Caribbean") ~ "Black",
      TRUE ~ "Others/Mixed/Unknown"
    ) %>% factor(), # 使用之前洗好的种族
    # 连续变量
    bmi_i2 = `Body mass index (BMI) | Instance 2`,
    waist_i2 = `Waist circumference | Instance 2`,
    age_mri = `Age when attended assessment centre | Instance 2`,
    tdi_val = `Townsend deprivation index at recruitment`
  ) %>%
  select(
    `Participant ID`, event_af_mri, duration_mri,ct1_f,
    `Liver iron corrected T1 (ct1) | Instance 2`, # 暴露
    age_mri, sex_f,tdi_val,ethnicity_f, bmi_i2, waist_i2, sbp_mri_val, dbp_mri_val,
    smoke_i2_f, alc_stati2_f,alc_freqi2_f,
    htn_mri_f, t2dm_mri_f, lip_mri_f, ckd_mri_f, cvd_mri_f,
    `Mean estimate of area of pericardial fat | Instance 2`,
    `LA maximum volume | Instance 2`,
    `Standard PRS for atrial fibrillation (AF)`,
    contains("Genetic principal components | Array"),
    `Proton density fat fraction (PDFF) | Instance 2`
    
  )

# MICE 设置
init_mri <- mice(df_mi_mri, maxit = 0)
meth_mri <- init_mri$method
pred <- init_mri$predictorMatrix
# 不插补 ID、结局、核心暴露和锁定变量
meth_mri[c("Participant ID","event_af_mri", "duration_mri", "Liver iron corrected T1 (ct1) | Instance 2", "Standard PRS for atrial fibrillation (AF)")] <- ""
pred_optimized_i2 <- quickpred(df_mi_mri, 
                            mincor = 0.1, 
                            exclude = "Participant ID", 
                            include = c("event_af_mri", "duration_mri", "Liver iron corrected T1 (ct1) | Instance 2",
                                        "age_mri", "sex_f", "sbp_mri_val","dbp_mri_val","ct1_f",
                                        "ethnicity_f", "bmi_i2","waist_i2", "alc_stati2_f","alc_freqi2_f", 
                                        "htn_mri_f", "t2dm_mri_f", "cvd_mri_f", "ckd_mri_f", "lip_mri_f","tdi_val",
                                        "Genetic principal components | Array 1", 
                                        "Genetic principal components | Array 2", 
                                        "Genetic principal components | Array 3",
                                        "Proton density fat fraction (PDFF) | Instance 2","Standard PRS for atrial fibrillation (AF)"))
# 执行插补
imp_mri_final<- mice(df_mi_mri, 
                     m = 10, 
                     maxit = 10, 
                     method = meth_mri, 
                     predictorMatrix = pred_optimized_i2, 
                     seed = 123,
                     printFlag = TRUE)

# ==============================================================================
# 第四部分：生成分析数据集 (mids 对象)
# ==============================================================================
# 1. 从插补对象中提取长表
imp_mri_long_new <- complete(imp_mri_final, action = "long", include = TRUE)

# 2. 在长表中划分 PDFF 分类 + 【新增】根据插补血压修正高血压定义
imp_mri_long_new <- imp_mri_long_new %>%
  mutate(
    # --- A. PDFF 5% 分类 ---
    pdff_5_f = if_else(`Proton density fat fraction (PDFF) | Instance 2` >= 5, ">=5%", "<5%"),
    pdff_5_f = factor(pdff_5_f, levels = c("<5%", ">=5%")), # 以 <5% 为参照组
    
    # --- B. 【新增】根据插补后的血压值修正高血压诊断 ---
    # 逻辑：原有诊断(1) OR 插补后收缩压>=140 OR 插补后舒张压>=90
    # 注意：先转为 character 或 numeric 比较，避免 factor 等级混乱
    htn_mri_f = if_else(
      as.character(htn_mri_f) == "1" | sbp_mri_val >= 140 | dbp_mri_val >= 90, 
      1, 
      0
    ),
    htn_mri_f = factor(htn_mri_f) # 重新封装为因子
  )

# 3. 将更新后的长表重新包装回 mids 对象
imp_mri_ready_pdff <- as.mids(imp_mri_long_new)
#saveRDS(imp_mri_ready_pdff, "imp_mri_ready_pdff.rds")













# ==============================================================================
# MRI 子队列数据质量验证 (QC Report)
# ==============================================================================

# 1. 样本量与房颤结局验证
qc_outcome <- df_mri_vars_full %>%
  summarise(
    Total_N = n(),
    Excl_Prevalent_AF = sum(is_prevalent_af_mri),
    Incident_AF_Events = sum(event_af_mri == 1),
    Median_FollowUp = median(duration_mri, na.rm = TRUE),
    Max_FollowUp = max(duration_mri, na.rm = TRUE),
    Min_FollowUp = min(duration_mri, na.rm = TRUE)
  )

# 2. 临床共病流行率验证 (Instance 2 截面)
qc_comorb <- df_mri_vars_full %>%
  summarise(
    HTN_Prev = mean(htn_mri_f == 1) * 100,
    T2DM_Prev = mean(t2dm_mri_f == 1) * 100,
    Lipid_Prev = mean(lip_mri_f == 1) * 100,
    CKD_Prev = mean(ckd_mri_f == 1) * 100,
    CVD_Prev = mean(cvd_mri_f == 1) * 100
  )

# 3. 时间逻辑自洽性检查 (核心：确保结局发生在 MRI 之后)
time_check <- df_mri_vars_full %>%
  filter(event_af_mri == 1) %>%
  mutate(time_diff = as.numeric(af_date - date_mri)) %>%
  summarise(
    Negative_Time_Count = sum(time_diff < 0), # 理论上应为 0
    Zero_Time_Count = sum(time_diff == 0)      # 理论上应为 0
  )

print("--- 1. 结局与随访统计 ---")
print(qc_outcome)
print("--- 2. MRI基线共病患病率 (%) ---")
print(round(qc_comorb, 2))
print("--- 3. 时间逻辑错误检查 ---")
print(time_check)










# ==============================================================================
# 质量检查报告
# ==============================================================================
cat("\n=== MRI 子队列构建报告 ===\n")
cat("1. 样本总量:", nrow(df_mri_vars_full), "\n")
cat("2. 发生房颤事件数:", sum(df_mri_vars_full$event_af_mri), "\n")
cat("3. 平均随访时间:", round(mean(df_mri_vars_full$duration_mri), 2), "年\n")
cat("4. cT1 缺失数 (应为0):", sum(is.na(df_mri_vars_full$ct1_val)), "\n")
cat("5. 关键中介变量缺失检查:\n")
print(colSums(is.na(df_mri_vars_full %>% select("Proton density fat fraction (PDFF) | Instance 2"))))
library(dplyr)
library(tidyr)
library(stringr)
library(officer)     
library(tableone)    
library(flextable)   
library(magrittr)

df_combined_1 <- fread("databcaa.csv", na.strings = c("", "NA"))

# 使用你定义的变量名列表进行选择
drink_cols_raw <- df_combined_1 %>%
  select(
    `eid`,
    # 动态引用你定义的所有酒精摄入量相关列名
    !!sym(red_status), !!sym(red_status2),
    !!sym(white_wine_status), !!sym(white_wine_status2),
    !!sym(beer_status), !!sym(beer_status2),
    !!sym(spirit_status), !!sym(spirit_status2),
    !!sym(fortified_status), !!sym(fortified_status2),
    !!sym(other_status), !!sym(other_status2),
    # 别忘了频率列，因为后面 Met-ALD 判断要用到
    !!sym(a_frequency), !!sym(a_frequency2)
  )
code_viral_icd <- "B15|B16|B17|B18|B19|070"
code_viral_sr <- "infectious/viral hepatitis|hepatitis a|hepatitis b|hepatitis c|hepatitis d|hepatitis e"
# 从原始 df 中提取 Instance 2 的评估日期作为 MRI 日期
mri_date_lookup <- df_mri_vars_full %>%
  select(`Participant ID`, 
         date_mri = `Date of attending assessment centre | Instance 2`) %>%
  mutate(date_mri = as.Date(date_mri)) %>%
  filter(!is.na(date_mri))

# 检查一下日期是否成功转换
summary(mri_date_lookup$date_mri)
# 这里的 ids_viral_mri_icd 逻辑就会引用上面定义的 date_mri
ids_viral_mri_icd <- df_all_icd_long %>%
  # 1. 关联 MRI 日期
  inner_join(mri_date_lookup, by = "Participant ID") %>%
  # 2. 匹配肝炎代码 (ICD10: B15-B19, 070)
  filter(str_detect(code, code_viral_icd)) %>%
  # 3. 关键逻辑：诊断日期必须在 MRI 日期之前或当天
  filter(as.Date(event_date) <= date_mri) %>% 
  pull(`Participant ID`) %>% 
  unique()
# ==============================================================================
# 1. 数据准备与变量清洗 (基于 MRI 子队列 df_mri_vars_full)
# ==============================================================================

# 假设 df_mri_vars_full 是上一轮代码生成的经过清洗的 MRI 数据集
# 必须包含：Participant ID, ct1_group, pdff_z, 饮酒量(Inst2), 以及各共病 _mri_f

df_table_prep <- df_mri_vars_full %>%
  left_join(drink_cols_raw, by = c("Participant ID" = "eid")) %>%
  mutate(
    # --- A. 时间相关变量平移 ---
    # 计算 MRI 时的年龄 (如果数据中没有直接提供 age_inst2)
    age_at_mri = `Age at recruitment` + as.numeric(as.Date(`Date of attending assessment centre | Instance 2`) - as.Date(`Date of attending assessment centre | Instance 0`))/365.25,
    age_at_mri = as.numeric(age_at_mri),
    # --- A. 提取 MRI 时的自我报告 (Instance 2) ---
    # 检查 Instance 2 的所有 self-reported 字段是否包含肝炎
    sr_viral_mri = if_any(contains("Non-cancer illness code, self-reported | Instance 2"), 
                          ~str_detect(tolower(as.character(.)), code_viral_sr)),
    
    # --- B. 综合判定：历史 ICD (ids_viral_mri_icd) OR 当下自我报告 ---
    viral_hep_status = if_else(
      (`Participant ID` %in% ids_viral_mri_icd) | replace_na(sr_viral_mri, FALSE), 
      1, 0
    ),
    # --- 2. 重新定义 MASLD ---
    pdff_val = `Proton density fat fraction (PDFF) | Instance 2`,
    # --- E. 补充缺失变量 (Smoking & PRS) ---
    # 吸烟状态 (通常使用 Instance 0 作为基线特征，或同步到 Instance 2)
    smoking_status_i2 = factor(`Smoking status | Instance 2`), 
    
    # PRS Tertile (同步主队列计算逻辑)
    prs_val = as.numeric(`Standard PRS for atrial fibrillation (AF)`),
    prs_tertile = ntile(prs_val, 3),
    prs_group_mri = factor(prs_tertile, levels = c(1, 2, 3), 
                           labels = c("Low PRS", "Intermediate PRS", "High PRS")),
    
    # 核心逻辑：PDFF >= 5% 且 排除病毒性肝炎
    masld_mri_status = if_else(pdff_val >= 5 & viral_hep_status == 0, 1, 0, missing = 0),
    
    # --- C. Met-ALD (基于 MRI PDFF + Instance 2 酒精) ---
    # 重新计算 Instance 2 的酒精克数 (为了保险起见，再次确认)
    
    # 1. 批量将酒精摄入列转化为数值，处理潜在的字符型问题
    across(c(!!sym(red_status2), !!sym(white_wine_status2), 
             !!sym(beer_status2), !!sym(spirit_status2), 
             !!sym(fortified_status2), !!sym(other_status2)), 
           ~as.numeric(as.character(.))),
    
    # 2. 现在计算就不会报错了，直接使用 replace_na(., 0)
    alc_g_inst2 = 
      (replace_na(!!sym(red_status2), 0) * 12) +
      (replace_na(!!sym(white_wine_status2), 0) * 12) +
      (replace_na(!!sym(beer_status2), 0) * 16) +
      (replace_na(!!sym(spirit_status2), 0) * 8) +
      (replace_na(!!sym(fortified_status2), 0) * 8) +
      (replace_na(!!sym(other_status2), 0) * 8),
    
    # ... 后续逻辑保持不变 ...
    
    
    metald_mri_status = case_when(
      masld_mri_status == 1 & (
        (Sex %in% c("Female", "0") & alc_g_inst2 >= 140 & alc_g_inst2 <= 350) |
          (Sex %in% c("Male", "1")   & alc_g_inst2 >= 210 & alc_g_inst2 <= 420) |
          (str_detect(as.character(`Alcohol intake frequency. | Instance 2`), "Daily"))
      ) ~ 1,
      TRUE ~ 0
    ),
    
    # --- D. 饮酒频率 (Instance 2) ---
    alc_raw_str_i2 = as.character(`Alcohol intake frequency. | Instance 2`),
    alcohol_freq_i2 = factor(case_when(
      alc_raw_str_i2 == "Daily or almost daily" ~ "Daily",
      alc_raw_str_i2 == "Never" ~ "Never",
      alc_raw_str_i2 %in% c("Three or four times a week", 
                            "Once or twice a week", 
                            "One to three times a month", 
                            "Special occasions only") ~ "Others",
      TRUE ~ as.character(NA)
    ), levels = c("Daily", "Others", "Never")),
    
    # --- E. 种族 (沿用基线) ---
    ethnicity_simplified = factor(case_when(
      grepl("White|British|Irish", as.character(`Ethnic background | Instance 0`)) ~ "White",
      grepl("Asian|Indian|Pakistani|Chinese", as.character(`Ethnic background | Instance 0`)) ~ "Asian",
      grepl("Black|African|Caribbean", as.character(`Ethnic background | Instance 0`)) ~ "Black",
      TRUE ~ "Other/Mixed"
    ), levels = c("White", "Asian", "Black", "Other/Mixed")),
    
    # --- F. 临床指标因子化 (0/1 -> No/Yes) ---
    # 注意：这里使用的是上一轮代码生成的 _mri_f 后缀变量
    across(any_of(c("htn_mri_f", "lip_mri_f", "t2dm_mri_f", 
                    "ckd_mri_f", "cvd_mri_f", "masld_mri_status", "metald_mri_status","viral_hep_status")), 
           ~ factor(., levels = c(0, 1), labels = c("No", "Yes"))),
    # 确保分层变量存在且是因子
    ct1_group = case_when(
      `Liver iron corrected T1 (ct1) | Instance 2` >= 750 ~ "High Risk",
      `Liver iron corrected T1 (ct1) | Instance 2` < 750 ~ "Low Risk"
    ),
    ct1_group = factor(ct1_group, levels = c("Low Risk", "High Risk"))
  )

# ==============================================================================
# 2. 映射表 (Rename Map) - 针对 MRI 变量
# ==============================================================================

rename_map_mri <- c(
  "Liver cT1 Group"              = "ct1_group",
  "Age at MRI scan, years"       = "age_at_mri",
  "Sex"                          = "Sex",
  "Race/Ethnicity"               = "ethnicity_simplified",
  "Body mass index (MRI), kg/m2" = "Body mass index (BMI) | Instance 2",
  "Townsend Deprivation Index"   = "Townsend deprivation index at recruitment",
  "Smoking status"               = "smoking_status_i2",
  
  "Diabetes mellitus"            = "t2dm_mri_f",
  "Hypertension"                 = "htn_mri_f",
  "Hyperlipidemia"               = "lip_mri_f",
  "Chronic kidney disease"       = "ckd_mri_f",
  "CVD History (excl. AF)"       = "cvd_mri_f",
  
  "NAFLD/MASLD"                  = "masld_mri_status",
  "Met-ALD"                      = "metald_mri_status",
  "Viral hepatitis"              = "viral_hep_status",
  
  "Liver PDFF, %"                = "pdff_val",
  "Alcohol Consumption (Inst. 2)" = "alcohol_freq_i2",
  "PRS Tertile Group"            = "prs_group_mri",
  
  "Incident AF (Post-MRI)"       = "event_af_mri",
  "Follow-up duration, years"       = "duration_mri"
)

# 准备用于 TableOne 的最终数据框
df_table_mri_final <- df_table_prep %>%
  # 将 event_af_mri 转为因子用于计数
  mutate(event_af_mri = factor(event_af_mri, levels=c(0,1), labels=c("No", "Yes"))) %>%
  rename(any_of(rename_map_mri)) %>%
  select(any_of(names(rename_map_mri)))

# ==============================================================================
# 3. 配置变量类型
# ==============================================================================

# 所有列名
all_vars_mri <- names(rename_map_mri)
# 排除分组变量本身
all_vars_mri <- setdiff(all_vars_mri, "Liver cT1 Group")

cat_vars_mri <- c(
  "Sex", "Race/Ethnicity", "Smoking status",
  "Diabetes mellitus", "Hypertension", "Hyperlipidemia", 
  "Chronic kidney disease", "CVD History (excl. AF)",
  "Incident AF (Post-MRI)", "NAFLD/MASLD", "Met-ALD", 
  "Viral hepatitis", "Alcohol Consumption (Inst. 2)", "PRS Tertile Group"
)

non_normal_vars_mri <- c(
  "Age at MRI scan, years", "Body mass index (MRI), kg/m2", 
  "Townsend Deprivation Index", "Follow-up duration, years",
  "Liver PDFF, %"
)

# ==============================================================================
# 4. 生成 TableOne 并导出
# ==============================================================================

message("正在生成 MRI 子队列基线表...")

# A. 创建对象
tab_mri <- CreateTableOne(
  vars = all_vars_mri,
  strata = "Liver cT1 Group",  # 按 cT1 分组
  data = df_table_mri_final,
  factorVars = cat_vars_mri,
  addOverall = TRUE
)

# B. 打印矩阵
tab_mat_mri <- print(
  tab_mri,
  nonnormal = non_normal_vars_mri,
  showAllLevels = TRUE,
  noSpaces = TRUE,
  printToggle = FALSE
)

# C. 格式美化
df_temp_mri <- as.data.frame(tab_mat_mri)
new_rownames <- rownames(tab_mat_mri)

# 优化行名显示
for (i in seq_along(new_rownames)) {
  if (grepl(" = ", new_rownames[i])) {
    level_name <- str_split(new_rownames[i], " = ", simplify = TRUE)[2]
    new_rownames[i] <- paste0("    ", str_to_sentence(level_name))
  }
}
df_temp_mri <- cbind(Characteristic = new_rownames, df_temp_mri)

# D. Flextable 渲染
doc_mri_base <- read_docx()

ft_mri <- flextable(df_temp_mri) %>%
  theme_booktabs() %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 9, part = "all") %>%
  # 表头加粗
  bold(part = "header") %>%
  # 变量名加粗，层级名变斜体
  bold(i = ~ !grepl("^    ", Characteristic), j = 1) %>% 
  italic(i = ~ grepl("^    ", Characteristic), j = 1) %>% 
  # 内容居中
  align(j = 2:ncol(df_temp_mri), align = "center", part = "all") %>%
  autofit() %>%
  set_caption(caption = "Table 1. Baseline Characteristics of the MRI Sub-cohort Stratified by Liver cT1")

# 写入 Word
doc_mri_base <- doc_mri_base %>% 
  body_add_flextable(value = ft_mri)

print(doc_mri_base, target = "Table1_MRI_Subcohort_Final.docx")

message("MRI 基线表已导出: Table1_MRI_Subcohort_Final.docx")