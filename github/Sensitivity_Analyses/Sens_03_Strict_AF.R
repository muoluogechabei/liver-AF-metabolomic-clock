#（一）洗数据
library(tidyverse)
library(data.table)
library(lubridate)
df <- fread("data.csv", na.strings = c("", "NA"))
df_raw <- df
df_codes_long <- df_raw %>%
  select(`Participant ID`, diag_icd10 = `Diagnoses - ICD10`) %>% # 假设 p41270 是你那一长串 ICD
  # 按 "|" 拆分成多行
  separate_rows(diag_icd10, sep = "\\|") %>%
  # 给拆出来的每一行打上序号 (a0, a1, a2...)，以便和日期列对齐
  group_by(`Participant ID`) %>%
  mutate(array_idx = paste0("a", row_number() - 1)) %>% 
  ungroup()

df_dates_long <- df_raw %>%
  select(`Participant ID`, contains("Date of first in-patient diagnosis - ICD10 | Array ")) %>%
  pivot_longer(
    cols = -`Participant ID`, 
    names_to = "array_idx", 
    values_to = "event_date"
  ) %>%
  # 修正：从列名中提取数字，并加上 "a" 前缀，确保和 codes 表完全一致
  mutate(
    array_idx = paste0("a", str_extract(array_idx, "\\d+$")),
    event_date = as.Date(event_date)
  ) %>%
  filter(!is.na(event_date))

# --- 第三步：合并两者，并筛选房颤 (I48) ---
# 合并编码和日期
df_combined <- df_codes_long %>%
  left_join(df_dates_long, by = c("Participant ID", "array_idx"))

# 提取纯房颤病人（排除心房扑动 I48.3, I48.4）
df_af_final <- df_combined %>%
  # 匹配 I48 但排除 I483 和 I484
  # 正则解释：^I48 后面不跟着 3 或 4
  filter(str_detect(diag_icd10, "I48") & !str_detect(diag_icd10, "I483|I484")) %>%
  group_by(`Participant ID`) %>%
  summarise(
    af_date = min(event_date, na.rm = TRUE),
    af_status = 1
  ) %>%
  ungroup()

# --- 处理 ICD9 编码 ---
df_codes9_long <- df_raw %>%
  # 假设 ICD9 的列名类似于 "Diagnoses - ICD9"
  select(`Participant ID`, diag_icd9 = `Diagnoses - ICD9`) %>% 
  separate_rows(diag_icd9, sep = "\\|") %>%
  group_by(`Participant ID`) %>%
  mutate(array_idx = paste0("a", row_number() - 1)) %>% 
  ungroup()

# --- 处理 ICD9 日期 ---
df_dates9_long <- df_raw %>%
  # 假设日期列名包含 "Date of first in-patient diagnosis - ICD9"
  select(`Participant ID`, contains("Date of first in-patient diagnosis - ICD9 | Array ")) %>%
  pivot_longer(
    cols = -`Participant ID`, 
    names_to = "array_idx", 
    values_to = "event_date"
  ) %>%
  mutate(
    array_idx = paste0("a", str_extract(array_idx, "\\d+$")),
    event_date = as.Date(event_date)
  ) %>%
  filter(!is.na(event_date))

# --- 合并并筛选 ICD9 纯房颤 (42731) ---
df_af_icd9 <- df_codes9_long %>%
  left_join(df_dates9_long, by = c("Participant ID", "array_idx")) %>%
  # 精确匹配 42731 (AF)，排除 42732 (AFL)
  # 使用 ^42731 确保是以该代码开头
  filter(str_detect(diag_icd9, "^42731")) %>% 
  group_by(`Participant ID`) %>%
  summarise(
    af_date_icd9 = min(event_date, na.rm = TRUE)
  ) %>%
  ungroup()
# --- 处理自述房颤 (Field 20002) ---
# 房颤在 20002 字段中的 Coding ID 是 1471

df_self_long <- df_raw %>%
  # 1. 选中所有自述疾病列 (Instance 0)
  # 使用 contains 或 matches 匹配 "Non-cancer illness code, self-reported | Instance 0"
  select(`Participant ID`, contains("Non-cancer illness code, self-reported | Instance 0")) %>%
  
  # 2. 宽表转长表，方便搜索
  pivot_longer(
    cols = -`Participant ID`, 
    names_to = "array_index", 
    values_to = "illness_code"
  ) %>%
  
  # 3. 筛选出代码为 1471，1483 的记录
  filter(illness_code %in% c("atrial fibrillation")) %>%
  
  # 4. 只要出现过一次，就标记该患者基线已有房颤
  distinct(`Participant ID`) %>%
  mutate(self_reported_af = 1)

# --- 终极大合并 ---
df_af_all_sources <- df_raw %>%
  # 1. 提取基线日期
  select(`Participant ID`, date_baseline = `Date of attending assessment centre | Instance 0`) %>% 
  mutate(date_baseline = as.Date(date_baseline)) %>%
  
  # 2. 关联三路数据
  left_join(df_af_final, by = "Participant ID") %>%  # ICD10 (提供 af_date)
  left_join(df_af_icd9, by = "Participant ID") %>%   # ICD9 (提供 af_date_icd9)
  left_join(df_self_long, by = "Participant ID") %>% # 自述 (提供 self_reported_af)
  
  # 3. 计算最早日期和最终状态
  rowwise() %>%
  mutate(
    # 计算 ICD 记录的最早日期（忽略 NA）
    # pmin 在处理两列对比时比 min 更简洁，且 handle NA 更好
    af_earliest_date = pmin(af_date, af_date_icd9, na.rm = TRUE),
    
    # 房颤状态：只要有一路信号就算（自述、ICD10 或 ICD9）
    af_status = if_else(!is.na(af_earliest_date) | (!is.na(self_reported_af) && self_reported_af == 1), 1, 0)
  ) %>%
  ungroup()

# 1. 将三路合并的房颤信息连回 df_raw
df_cleaned <- df_raw %>%
  # 关联 ICD10/9 的最早日期
  left_join(df_af_final %>% select(`Participant ID`, af_date), by = "Participant ID") %>%
  left_join(df_af_icd9 %>% select(`Participant ID`, af_date_icd9), by = "Participant ID") %>%
  # 关联自述标志
  left_join(df_self_long, by = "Participant ID") %>%
  
  # 2. 在主表里进行统一判定
  mutate(
    # 统一转换基线日期格式
    baseline_dt = as.Date(`Date of attending assessment centre | Instance 0`),
    
    # 取 ICD 系统记录的最早房颤日期
    earliest_icd_dt = pmin(af_date, af_date_icd9, na.rm = TRUE),
    
    # 【定义基线房颤 (Prevalent AF)】
    # 只要符合以下任一条件，即为 1：
    # (a) 自述过房颤 (self_reported_af == 1)
    # (b) ICD 记录日期在基线当天或之前
    is_prevalent_af = if_else(
      (!is.na(self_reported_af) & self_reported_af == 1) | 
        (!is.na(earliest_icd_dt) & earliest_icd_dt <= baseline_dt),
      1, 0, missing = 0
    )
  ) %>%
  
  # 3. 【执行排除】
  # 只留下基线时“心脏健康”的人
  filter(is_prevalent_af == 0)

# 4. 检查排除后的样本量
print(paste("排除基线房颤后，剩余样本量为:", nrow(df_cleaned)))

library(tidyverse)
library(mice)

# 1. 定义核心变量名（基于你提供的列名）
core_vars <- c(
  "Alanine aminotransferase | Instance 0", 
  "Aspartate aminotransferase | Instance 0", 
  "Platelet count | Instance 0", 
  "Albumin | Instance 0",
  "Body mass index (BMI) | Instance 0",
  "Age at recruitment",
  "Standard PRS for atrial fibrillation (AF)"
)

# 10个遗传主成分的列名
pc_vars <- paste0("Genetic principal components | Array ", 1:10)

# 2. 执行排除逻辑
df_subs <- df_cleaned %>%
  # 排除生化、BMI、年龄、糖尿病标志缺失的人
  filter(across(all_of(core_vars), ~ !is.na(.))) %>%
  # 排除 PRS 和 10 个 PC 缺失的人
  filter(across(all_of(pc_vars), ~ !is.na(.)))

print(paste("排除核心指标、PRS及PC缺失后，剩余样本量：", nrow(df_subs)))


# --- 提取糖尿病相关 ID ---

# 1. 先建立一个 ID 和 基线日期 的对照表
baseline_lookup <- df_subs %>%
  select(`Participant ID`, baseline_dt)

# 2. 提取 ICD10 糖尿病 ID (确保在基线前)
dm_icd10_ids <- df_combined %>%
  # 关联基线日期
  left_join(baseline_lookup, by = "Participant ID") %>%
  # 筛选：代码匹配 E10-E14 且 日期在基线前
  filter(str_detect(diag_icd10, "^E1[0-4]"), event_date <= baseline_dt) %>%
  pull(`Participant ID`) %>% 
  unique()

# --- 修正后的 ICD9 糖尿病提取 ---
dm_icd9_ids <- df_codes9_long %>%
  # 只要匹配到 250 开头的 ICD9 编码，就认为是基线糖尿病
  # 因为 ICD9 是英国早期的住院记录系统，出现在系统里的时间必然早于 UKB 基线 (2006-2010)
  filter(str_detect(diag_icd9, "^250")) %>%
  pull(`Participant ID`) %>% 
  unique()

# 3. 自述 (20002 编码 1220, 1221, 1222, 1223)
dm_self_ids <- df_raw %>%
  select(`Participant ID`, contains("Non-cancer illness code, self-reported | Instance 0")) %>%
  pivot_longer(-`Participant ID`, values_to = "code") %>%
  filter(code %in% c("diabetes", "gestational diabetes", "type 1 diabetes", "type 2 diabetes")) %>%
  pull(`Participant ID`) %>% unique()

# --- 在 df_subs 中标记整合后的糖尿病 ---
df_ready <- df_subs %>%
  mutate(
    is_dm_clinical = if_else(`Participant ID` %in% c(dm_icd10_ids, dm_icd9_ids, dm_self_ids), 1, 0),
    is_dm_biochem = if_else(
      (`Glycated haemoglobin (HbA1c) | Instance 0` >= 48) | 
        (`Glucose | Instance 0` >= 11.1), 1, 0, missing = 0
    ),
    # 最终糖尿病状态 (临床诊断 OR 生化异常 OR 医生口头诊断)
    final_diabetes = if_else(
      is_dm_clinical == 1 | is_dm_biochem == 1 | `Diabetes diagnosed by doctor | Instance 0` == 1, 1, 0, missing = 0
    )
  )
print(paste("整合判定出的糖尿病总人数:", sum(df_ready$final_diabetes)))

# 1. 设定行政截断日期
admin_end_date <- as.Date("2022-05-31")

# 2. 准备死亡日期
death_dates <- df_raw %>%
  select(`Participant ID`, date_death = `Date of death | Instance 0`) %>%
  mutate(date_death = as.Date(date_death))

# 3. 计算随访指标
df_survival <- df_ready %>%
  left_join(death_dates, by = "Participant ID") %>%
  # 这里的 af_date 是你之前洗出的 ICD10/9 综合最早日期
  mutate(
    # 随访终点：发病、死亡、行政截断三者取最早
    follow_up_end = pmin(af_date, date_death, admin_end_date, na.rm = TRUE),
    
    # 结局状态：只有在行政截断前发病才算 1
    event_af = if_else(!is.na(af_date) & af_date <= pmin(date_death, admin_end_date, na.rm = TRUE), 1, 0),
    
    # 计算时间 (年)
    duration = as.numeric(follow_up_end - baseline_dt) / 365.25
  ) 
# --- 计算核心肝纤维化评分 ---
df_with_scores <- df_survival %>%
  mutate(
    # 提取并重命名基础指标，方便计算
    age_val = `Age at recruitment`,
    alt_val = `Alanine aminotransferase | Instance 0`,
    ast_val = `Aspartate aminotransferase | Instance 0`,
    plt_val = `Platelet count | Instance 0`,
    alb_val = `Albumin | Instance 0`/10,
    bmi_val = `Body mass index (BMI) | Instance 0`,
    
    # 1. FIB-4 Index
    # 公式: (Age * AST) / (PLT * sqrt(ALT))
    fib4 = (age_val * ast_val) / (plt_val * sqrt(alt_val)),
    
    # 2. APRI (AST to Platelet Ratio Index)
    # 公式: ((AST / AST_ULN) / PLT) * 100  (通常 AST 上限取 40 U/L)
    apri = ((ast_val / 40) / plt_val) * 100,
    
    # 3. NFS (NAFLD Fibrosis Score)
    # 公式: -1.675 + 0.037*age + 0.094*BMI + 1.13*diabetes(yes=1) + 
    #       0.99*(AST/ALT) - 0.013*PLT - 0.66*albumin
    nfs = -1.675 + (0.037 * age_val) + (0.094 * bmi_val) + 
      (1.13 * final_diabetes) + (0.99 * (ast_val / alt_val)) - 
      (0.013 * plt_val) - (0.66 * alb_val)
  )

# 检查评分分布，看看有没有异常大值
summary(df_with_scores %>% select(fib4, apri, nfs))
# 检查是否有评分计算失败 (NA)
colSums(is.na(df_with_scores %>% select(fib4, apri, nfs)))

# --- 整合计算与临床分组 ---
df_final_clean <- df_with_scores %>%
  mutate(
    # 1. 计算 AST/ALT 比值 (De Ritis Ratio)
    ast_alt_ratio = ast_val / alt_val,
    
    # 2. FIB-4 分组: <1.3 (低), 1.3-2.67 (中), >=2.67 (高)
    fib4_cat = case_when(
      fib4 < 1.3 ~ "Low Risk",
      fib4 >= 1.3 & fib4 < 2.67 ~ "Intermediate Risk",
      fib4 >= 2.67 ~ "High Risk"
    ),
    fib4_cat = factor(fib4_cat, levels = c("Low Risk", "Intermediate Risk", "High Risk")),
    
    # 3. NFS 分组: < -1.455 (低), -1.455 to 0.675 (中), >= 0.675 (高)
    nfs_cat = case_when(
      nfs < -1.455 ~ "Low Risk",
      nfs >= -1.455 & nfs < 0.675 ~ "Intermediate Risk",
      nfs >= 0.675 ~ "High Risk"
    ),
    nfs_cat = factor(nfs_cat, levels = c("Low Risk", "Intermediate Risk", "High Risk")),
    
    # 4. APRI 分组: <1 (低), >=1 (高)
    apri_cat = case_when(
      apri < 0.5 ~ "Low Risk",
      apri >= 0.5 ~ "High Risk"
    ),
    apri_cat = factor(apri_cat, levels = c("Low Risk", "High Risk")),
    
    # 5. AST/ALT 比值分组: <1 (低), 1-2 (中), >2 (高)
    ratio_cat = case_when(
      ast_alt_ratio < 1 ~ "Low (<1)",
      ast_alt_ratio >= 1 & ast_alt_ratio <= 2 ~ "Intermediate (1-2)",
      ast_alt_ratio > 2 ~ "High (>2)"
    ),
    ratio_cat = factor(ratio_cat, levels = c("Low (<1)", "Intermediate (1-2)", "High (>2)"))
  )


# 1. 确保 data_sip 已经加载
df_sup <- fread("data_sup.csv") 

# 2. 执行合并
# 我们使用 left_join，以你洗好的 40万+ 样本为准
df_final_clean_v2 <- df_final_clean %>%
  left_join(df_sup, by = "Participant ID")

# 3. 核心质量检查 (Audit)
if(nrow(df_final_clean_v2) == nrow(df_final_clean)) {
  message("✅ 样本量对齐：合并成功，样本总数依然是 ", nrow(df_final_clean_v2))
} else {
  warning("⚠️ 样本量异常：合并后行数变了，可能存在重复 ID！")
}

# 4. 检查新并入的字段（以 191 为例）
# 看看失访日期列（Field 191）并进来后的非空记录有多少
new_cols <- setdiff(names(df_sup), names(df_final_clean))
message("新并入的变量包括: ", paste(new_cols, collapse = ", "))

# 5. 更新你的“分析终极版”文件
#saveRDS(df_final_clean_v2, "df_final_clean_v2.rds")





# 如果你的列名不同，请修改下面的 matches 部分
df_final_clean_v3 <- df_final_clean_v2 %>%
  mutate(
    # 1. 转换为日期格式
    date_lost = as.Date(`Date lost to follow-up`),
    
    # 2. 更新随访终点：发病、死亡、失访、行政截断 四者取最早
    # 注意：na.rm = TRUE 必须加，否则只要有一个 NA 结果就是 NA
    follow_up_end_updated = pmin(af_date, date_death, date_lost, admin_end_date, na.rm = TRUE),
    
    # 3. 更新结局状态：只有在“更新后的终点”确实是由于房颤引起时，才算 1
    
    # 只有在死亡、失访或截断之前（或当天）发生的房颤，才计入结局
    event_af_updated = if_else(!is.na(af_date) & af_date <= follow_up_end_updated, 1, 0),
    # 4. 重新计算随访年限
    duration_updated = as.numeric(follow_up_end_updated - baseline_dt) / 365.25
  )
library(dplyr)
library(stringr)
library(tidyr)
library(data.table)

# =======================================================
# 1. 准备工作：建立 ID 和 基线日期的对照
# =======================================================
# 提取基线日期对照表 (如果 df_final_clean_v3 里有 baseline_dt)
baseline_lookup <- df_final_clean_v3 %>%
  select(`Participant ID`, baseline_dt) %>%
  distinct()

# =======================================================
# 2. 定义 ICD 和 自述的匹配模式 (Regex & Codes)
# =======================================================

# --- A. ICD 10/9 正则模式 ---
# 注意：CVD 比较复杂，包含冠心病、心衰、卒中、外周血管、瓣膜病
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


# =======================================================
# 3. 核心步骤：整合 ICD10 & ICD9 并提取基线前确诊 ID
# =======================================================

# 3.1 统一 ICD10 和 ICD9 为一个长表格
df_all_icd_long <- bind_rows(
  # 处理 ICD10
  df_combined %>% 
    select(`Participant ID`, array_idx, event_date, code = diag_icd10) %>%
    mutate(icd_version = 10),
  # 处理 ICD9 (关联日期)
  df_codes9_long %>%
    left_join(df_dates9_long, by = c("Participant ID", "array_idx")) %>%
    select(`Participant ID`, array_idx, event_date, code = diag_icd9) %>%
    mutate(icd_version = 9)
)
# 3.2 关联基线日期并过滤
long_records_clean <- df_all_icd_long %>%
  left_join(baseline_lookup, by = "Participant ID") %>%
  # [关键]：过滤掉缺失日期和基线后发生的记录
  filter(!is.na(event_date), event_date <= baseline_dt) %>%
  filter(!is.na(code))

# 3.3 分别提取各疾病的 ID 列表
# 现在直接在统一的 code 列中匹配
ids_htn_icd_base <- long_records_clean %>% filter(str_detect(code, pat_htn_icd)) %>% pull(`Participant ID`) %>% unique()
ids_lip_icd_base <- long_records_clean %>% filter(str_detect(code, pat_lipids_icd)) %>% pull(`Participant ID`) %>% unique()
ids_ckd_icd_base <- long_records_clean %>% filter(str_detect(code, pat_ckd_icd)) %>% pull(`Participant ID`) %>% unique()
ids_cvd_icd_base <- long_records_clean %>% filter(str_detect(code, pat_cvd_icd)) %>% pull(`Participant ID`) %>% unique()

# 糖尿病特异性 ID (基线前)
ids_t1dm_icd_base <- long_records_clean %>% filter(str_detect(code, pat_t1dm_icd)) %>% pull(`Participant ID`) %>% unique()
ids_t2dm_icd_base <- long_records_clean %>% filter(str_detect(code, pat_t2dm_icd)) %>% pull(`Participant ID`) %>% unique()
ids_gdm_icd_base  <- long_records_clean %>% filter(str_detect(code, pat_gdm_icd)) %>% pull(`Participant ID`) %>% unique()

# =======================================================
# 4. 主表整合 (Instance 0 数据 + 上面提取的 ID)
# =======================================================

df_clinical_final <- df_final_clean_v3 %>%
  mutate(
    # --- A. 药物判定 (Instance 0) ---
    # 只要在基线药单里出现，就是基线用药
    med_htn  = if_any(contains("Treatment/medication code | Instance 0"), ~str_detect(tolower(as.character(.)), pat_htn_meds) %in% TRUE),
    med_lip  = if_any(contains("Treatment/medication code | Instance 0"), ~str_detect(tolower(as.character(.)), pat_lipid_meds) %in% TRUE),
    med_t2dm = if_any(contains("Treatment/medication code | Instance 0"), ~str_detect(tolower(as.character(.)), pat_t2dm_meds) %in% TRUE),
    
    # --- B. 自述判定 (Instance 0) ---
    # 只要在基线问卷里说了，就是基线患病
    sr_htn  = if_any(contains("Non-cancer illness code, self-reported | Instance 0"), ~as.character(.) %in% codes_sr_htn),
    sr_lip  = if_any(contains("Non-cancer illness code, self-reported | Instance 0"), ~as.character(.) %in% codes_sr_lip),
    sr_ckd  = if_any(contains("Non-cancer illness code, self-reported | Instance 0"), ~as.character(.) %in% codes_sr_ckd),
    sr_cvd  = if_any(contains("Non-cancer illness code, self-reported | Instance 0"), ~as.character(.) %in% codes_sr_cvd),
    sr_t2dm = if_any(contains("Non-cancer illness code, self-reported | Instance 0"), ~as.character(.) %in% codes_sr_t2dm),
    sr_t1dm = if_any(contains("Non-cancer illness code, self-reported | Instance 0"), ~as.character(.) %in% codes_sr_t1dm),
    
    # --- C. 生化/体测指标 (Instance 0) ---
    sbp_mean = rowMeans(select(., contains("Systolic blood pressure, automated reading | Instance 0")), na.rm = TRUE),
    dbp_mean = rowMeans(select(., contains("Diastolic blood pressure, automated reading | Instance 0")), na.rm = TRUE),
    
    # eGFR 计算 (Instance 0 肌酐)
    scr_mgdl = `Creatinine | Instance 0` / 88.4,
    kappa = if_else(Sex == "Female", 0.7, 0.9),
    alpha = if_else(Sex == "Female", -0.241, -0.302),
    egfr = 142 * pmin(scr_mgdl/kappa, 1)^alpha * pmax(scr_mgdl/kappa, 1)^-1.200 * 0.9938^`Age at recruitment` * (if_else(Sex == "Female", 1.012, 1.0)),
    
    # ----------------------------------------------------------------
    # --- D. 最终基线判定 (Logic Assembly) ---
    # ----------------------------------------------------------------
    
    # 1. 高血压: ID在基线前ICD名单 OR 自述 OR 药物 OR 测值高
    baseline_htn = if_else(
      (`Participant ID` %in% ids_htn_icd_base) | 
        sr_htn | 
        med_htn | 
        (sbp_mean >= 140 | dbp_mean >= 90), 1, 0, missing = 0),
    
    # 2. 血脂异常: ID在基线前ICD名单 OR 自述 OR 药物 OR 测值异常
    baseline_lipids = if_else(
      (`Participant ID` %in% ids_lip_icd_base) | 
        sr_lip | 
        med_lip |
        (Sex == "Female" & `HDL cholesterol | Instance 0` < 1.29) |
        (Sex == "Male" & `HDL cholesterol | Instance 0` < 1.03) |
        (`Triglycerides | Instance 0` >= 1.7), 1, 0, missing = 0),
    
    # 3. CKD: ID在基线前ICD名单 OR 自述 OR eGFR<60
    baseline_ckd = if_else(
      (`Participant ID` %in% ids_ckd_icd_base) | 
        sr_ckd | 
        egfr < 60, 1, 0, missing = 0),
    
    # 4. CVD: ID在基线前ICD名单 OR 自述
    baseline_cvd = if_else(
      (`Participant ID` %in% ids_cvd_icd_base) | 
        sr_cvd, 1, 0, missing = 0)
  ) %>%
  
  # --- E. 糖尿病复杂判定 (T1DM 排除法) ---
  mutate(
    # 判断是否为1型 (自述1型 OR ICD有1型记录)
    is_t1dm_evidence = (`Participant ID` %in% ids_t1dm_icd_base) | sr_t1dm | (`Participant ID` %in% ids_gdm_icd_base),
    
    # 判断是否有2型证据 (ICD 2型 OR 药物 OR 自述2型 OR 生化指标)
    has_t2dm_evidence = (
      (`Participant ID` %in% ids_t2dm_icd_base) | 
        med_t2dm | 
        sr_t2dm | 
        (`Glycated haemoglobin (HbA1c) | Instance 0` >= 48) | 
        (`Glucose | Instance 0` >= 11.1)
    ),
    
    # 最终判定: 有T2DM证据 且 没有T1DM证据
    baseline_t2dm = if_else(has_t2dm_evidence & !is_t1dm_evidence, 1, 0, missing = 0)
  )

# =======================================================
# 5. 质量检查 (QC)
# =======================================================
cat("--- 基线患病率检查 (应符合一般人群常识) ---\n")
print(prop.table(table(df_clinical_final$baseline_htn)) * 100)
print(prop.table(table(df_clinical_final$baseline_t2dm)) * 100)
print(prop.table(table(df_clinical_final$baseline_ckd)) * 100)
print(prop.table(table(df_clinical_final$baseline_cvd)) * 100)
#saveRDS(df_clinical_final, "df_clinical_final.rds")


# 1. 定义你数据集中确实存在的协变量名
covariates_to_check <- c(
  "Townsend deprivation index at recruitment",
  "Smoking status | Instance 0",
  "Systolic blood pressure, automated reading | Instance 0 | Array 0",
  "Systolic blood pressure, automated reading | Instance 0 | Array 1",
  "Diastolic blood pressure, automated reading | Instance 0 | Array 0",
  "Diastolic blood pressure, automated reading | Instance 0 | Array 1",
  "HDL cholesterol | Instance 0",
  "Triglycerides | Instance 0",
  "Creatinine | Instance 0",
  "Age at recruitment",
  "Alcohol drinker status | Instance 0",
  "Alcohol intake frequency. | Instance 0",
  "Sex",
  "Waist circumference | Instance 0",
  "Gamma glutamyltransferase | Instance 0",
  "Glycated haemoglobin (HbA1c) | Instance 0",
  "Platelet count | Instance 0",
  "Albumin | Instance 0",
  "Glucose | Instance 0"
)

# 2. 检查缺失情况
missing_summary <- df_clinical_final %>%
  select(all_of(covariates_to_check)) %>%
  summarise(across(everything(), list(
    na_count = ~sum(is.na(.)),
    na_percent = ~round(sum(is.na(.)) / n() * 100, 2)
  ))) %>%
  pivot_longer(everything(), names_to = c("variable", "stat"), names_sep = "_(?!.*_)") %>% # 使用正则匹配最后一个下划线
  pivot_wider(names_from = stat, values_from = value)

print(missing_summary)

# =======================================================
# 第二部分：准备 MICE 数据集 (Data Prep)
# =======================================================
library(dplyr)
library(stringr)
library(tidyr)
library(mice)
library(parallel)

# =======================================================
# 1. 数据准备：加入种族、腰围、饮酒状态
# =======================================================

df_mi_strict <- df_clinical_final %>%
  mutate(
    # --- A. 显式转换为 Factor ---
    sex_f      = factor(Sex),
    smoke_f    = factor(`Smoking status | Instance 0`), 
    alc_status_f = factor(`Alcohol drinker status | Instance 0`),
    alc_freq_f = factor(`Alcohol intake frequency. | Instance 0`), 
    ethnicity_f = case_when(
      str_detect(as.character(`Ethnic background | Instance 0`), "White|British|Irish") ~ "White",
      str_detect(as.character(`Ethnic background | Instance 0`), "Asian|Indian|Pakistani|Chinese|Bangladesh") ~ "Asian",
      str_detect(as.character(`Ethnic background | Instance 0`), "Black|African|Caribbean") ~ "Black",
      TRUE ~ "Others/Mixed/Unknown"
    ) %>% factor(),
    
    # 共病因子 (Hard Facts)
    htn_f      = factor(baseline_htn),
    t2dm_f     = factor(baseline_t2dm),
    cvd_f      = factor(baseline_cvd),
    ckd_f      = factor(baseline_ckd),
    lip_f      = factor(baseline_lipids),
    
    # --- B. 连续变量预处理 ---
    bmi_val   = `Body mass index (BMI) | Instance 0`,
    sbp_avg = rowMeans(across(c(`Systolic blood pressure, automated reading | Instance 0 | Array 0`, 
                                `Systolic blood pressure, automated reading | Instance 0 | Array 1`)), na.rm = TRUE),
    dbp_avg = rowMeans(across(c(`Diastolic blood pressure, automated reading | Instance 0 | Array 0`, 
                                `Diastolic blood pressure, automated reading | Instance 0 | Array 1`)), na.rm = TRUE),
    # 将 NaN 转回 NA（rowMeans 在全缺失时可能返回 NaN）
    sbp_val = ifelse(is.nan(sbp_avg), NA, sbp_avg),
    dbp_val = ifelse(is.nan(dbp_avg), NA, dbp_avg),
    tdi_val = `Townsend deprivation index at recruitment`,
    # [新加入] 腰围
    waist_val = `Waist circumference | Instance 0`
  ) %>%
  select(
    `Participant ID`, 
    event_af_updated, duration_updated,
    
    # 核心暴露
    fib4, nfs, apri, ast_alt_ratio,
    fib4_cat, nfs_cat, apri_cat, ratio_cat,
    # 预测因子/待插补变量 (加入新变量)
    age_val, sex_f, ethnicity_f, # +种族
    bmi_val, waist_val,          # +腰围
    sbp_val, dbp_val, 
    smoke_f, alc_freq_f, alc_status_f, # +饮酒状态
    `Alanine aminotransferase | Instance 0`,
    `Aspartate aminotransferase | Instance 0`,
    `Platelet count | Instance 0`,
    `Albumin | Instance 0`,
    tdi_val,
    `HDL cholesterol | Instance 0`,
    `Triglycerides | Instance 0`,
    `Creatinine | Instance 0`,
    `Glycated haemoglobin (HbA1c) | Instance 0`,
    `Glucose | Instance 0`,
    `Gamma glutamyltransferase | Instance 0`,
    # 共病因子
    htn_f, t2dm_f, cvd_f, ckd_f, lip_f,
    
    # 锁定变量
    `Standard PRS for atrial fibrillation (AF)`,
    contains("Genetic principal components | Array")
  )


# =======================================================
# 2. MICE 配置 (Method & Predictor Matrix)
# =======================================================

# 初始试运行
init <- mice(df_mi_strict, maxit = 0) 
meth <- init$method
pred <- init$predictorMatrix

# --- Method 设置 ---
# 锁定不需要插补的列
meth["Standard PRS for atrial fibrillation (AF)"] <- ""
meth["event_af_updated"] <- ""
meth["duration_updated"] <- ""
meth["Participant ID"] <- ""

# --- Predictor Matrix 设置 ---
# A. 排除 ID
pred[, "Participant ID"] <- 0

# B. 优化预测关系 (Quickpred)
# [关键] 把新加入的 ethnicity_f, waist_val, alc_status_f 加入 include
pred_optimized <- quickpred(df_mi_strict, 
                            mincor = 0.1, 
                            exclude = "Participant ID", 
                            include = c("event_af_updated", "duration_updated", 
                                        "age_val", "sex_f", "sbp_val","dbp_val",
                                        "ethnicity_f", "bmi_val","waist_val", "alc_status_f","alc_freq_f", 
                                        "Alanine aminotransferase | Instance 0",
                                        "Aspartate aminotransferase | Instance 0",
                                        "Platelet count | Instance 0",
                                        "Glucose | Instance 0",
                                        "Glycated haemoglobin (HbA1c) | Instance 0",
                                        "Albumin | Instance 0",# <--- 新变量作为强预测因子
                                        "htn_f", "t2dm_f", "cvd_f", "ckd_f", "lip_f","tdi_val",
                                        "Genetic principal components | Array 1", 
                                        "Genetic principal components | Array 2", 
                                        "Genetic principal components | Array 3"))

# C. 锁定 Target (Hard Facts 不被插补)
vars_comorb <- c("htn_f", "t2dm_f", "cvd_f", "ckd_f", "lip_f")


vars_no_na <- c(vars_comorb, "Standard PRS for atrial fibrillation (AF)")
meth[vars_no_na] <- "" 

vars_cats <- c("fib4_cat", "nfs_cat", "apri_cat", "ratio_cat")
meth[vars_cats] <- ""

# 注意：ethnicity_f 和 alc_status_f 如果有缺失值，它们会被插补。
# 如果你认为它们也是"Hard Facts" (不应该变)，也可以加入上面的 vars_comorb 列表锁定。
# 但通常 UKB 中这两个可能有少量缺失，建议允许插补。

# =======================================================
# 3. 执行插补 (Final Run)
# =======================================================
# 变量变多了，建议跑的时候去喝杯咖啡
imp_final <- mice(df_mi_strict, 
                  m = 10, 
                  maxit = 10, 
                  method = meth, 
                  predictorMatrix = pred_optimized, 
                  seed = 123,
                  printFlag = TRUE)

# =======================================================
# 4. 保存与检查
# =======================================================
#saveRDS(imp_final, "mice_imputation_result_full_v2.rds")








# 检查
df_check <- complete(imp_final, 1)
cat("检查关键变量是否补全:\n")
colSums(is.na(df_check[, c("waist_val", "ethnicity_f", "alc_status_f", "bmi_val")]))

#saveRDS(df_clinical_final, "df_clinical_final_pureAF.rds")
#plot(mice_imputation_result_full_v2, c("waist_val", "sbp_val", "dbp_val", "tdi_val","`Glycated haemoglobin (HbA1c) | Instance 0`","`Glucose | Instance 0`","`Gamma glutamyltransferase | Instance 0`"))
# 观察不同插补集（.imp）中 smoke_f 的分布
# 如果某一个插补集的分布显著异于其他集，说明该变量插补不稳定



cat("新发房颤平均随访年限:", round(mean(df_clinical_final$duration_updated, na.rm=TRUE), 2), "年\n")
cat("新发房颤中位随访年限:", round(median(df_clinical_final$duration_updated, na.rm=TRUE), 2), "年\n")
cat("随访四分位区间 (IQR):", 
    round(quantile(df_clinical_final$duration_updated, 0.25, na.rm=TRUE), 2), "-", 
    round(quantile(df_clinical_final$duration_updated, 0.75, na.rm=TRUE), 2), "年\n")
cat("新发房颤事件数:", sum(df_clinical_final$event_af_updated, na.rm=TRUE), "例\n")



# 1. 重新提取长格式，但这次包含原始数据 (include = TRUE)
# 这一步很快，不涉及重新插补
imp_long_with_orig <- complete(imp_final, action = "long", include = TRUE)

# 2. 把你刚才算的疾病逻辑应用到这个带原始数据的大表上
# 注意：mutate 之后，.imp = 0 (原始) 和 .imp = 1:10 (插补) 都会被计算
# 2. 重新计算诊断逻辑
imp_long_diagnosed <- imp_long_with_orig %>%
  mutate(
    # --- [关键补全]：先计算基于插补值的 eGFR ---
    scr_mgdl = `Creatinine | Instance 0` / 88.4,
    kappa = if_else(sex_f == "Female", 0.7, 0.9),
    alpha = if_else(sex_f == "Female", -0.241, -0.302),
    # CKD-EPI 2021 公式
    egfr = 142 * pmin(scr_mgdl/kappa, 1)^alpha * pmax(scr_mgdl/kappa, 1)^-1.200 * 0.9938^age_val * (if_else(sex_f == "Female", 1.012, 1.0)),
    
    # --- 现在的诊断逻辑就完美了 ---
    
    # 高血压判定 (安全对比)
    htn_final = factor(if_else(as.character(htn_f) %in% "1" | sbp_val >= 140 | dbp_val >= 90, 1, 0)),
    
    # T2DM 判定 (加入随机血糖)
    t2dm_final = factor(if_else(as.character(t2dm_f) %in% "1" | 
                                  `Glycated haemoglobin (HbA1c) | Instance 0` >= 48 | 
                                  `Glucose | Instance 0` >= 11.1, 1, 0)),
    
    # CKD 判定
    ckd_final = factor(if_else(as.character(ckd_f) %in% "1" | egfr < 60, 1, 0))
  )

# 3. 这次再封装，就不会报错了喵！
imp_ready_for_cox <- as.mids(imp_long_diagnosed)
saveRDS(imp_ready_for_cox, "mice_rediag_pureAF.rds")






# =======================================================
# 验证：插补前后（.imp=0 vs .imp>0）病例数对比报告
# =======================================================

library(tidyr)

# 1. 汇总统计
diagnosis_verification <- imp_long_diagnosed %>%
  group_by(.imp) %>%
  summarise(
    # 样本总数
    Total_N = n(),
    # 高血压人数及比例
    HTN_Count = sum(htn_final == "1", na.rm = TRUE),
    HTN_NA = sum(is.na(htn_final)),
    # 糖尿病人数及比例
    T2DM_Count = sum(t2dm_final == "1", na.rm = TRUE),
    T2DM_NA = sum(is.na(t2dm_final)),
    # CKD人数及比例
    CKD_Count = sum(ckd_final == "1", na.rm = TRUE),
    CKD_NA = sum(is.na(ckd_final))
  ) %>%
  mutate(
    # 标记哪些是插补后的
    Type = if_else(.imp == 0, "Original (Raw)", "Imputed (m)")
  )

# 2. 打印对比表格
cat("--- 诊断结果对比表 --- \n")
print(as.data.frame(diagnosis_verification))

# 3. 计算“增量”：看看插补帮你多抓住了多少病人
cat("\n--- 插补带来的‘增量’分析 (以第1个插补集为例) --- \n")
diff_check <- diagnosis_verification %>%
  filter(.imp %in% c(0, 1)) %>%
  summarise(
    Additional_HTN = HTN_Count[2] - HTN_Count[1],
    Additional_T2DM = T2DM_Count[2] - T2DM_Count[1],
    Additional_CKD = CKD_Count[2] - CKD_Count[1]
  )
print(diff_check)
#writeLines(names(df_clinical_final), "df_clinical_final_columns.txt")
library(tidyverse)
library(mice)
imp_ready_for_cox<-mice_rediag_pureAF
# =======================================================
# 1. 定义常量与编码
# =======================================================
follow_up_end_dt <- as.Date("2022-10-31") # 随访截止日期
pat_heme_cancer  <- "^C8[1-9]|^C9[0-6]|^20[0-8]" # 血液肿瘤
# [关键] 纯房颤正则：匹配 I480-I482, I489 (支持带点或不带点) 以及 ICD9 427.31
pat_pure_af      <- "I48[0129]|I48\\.[0129]|42731|427\\.31" 

# =======================================================
# 2. 准备基础数据 (ID名单与日期)
# =======================================================

# 2.1 提取“纯房颤”最早发病日期
# 注意：这里从原始 ICD 长表中提取，确保来源最纯净
df_pure_af_dates <- df_all_icd_long %>%
  # 清洗 code 中的点，统一格式
  mutate(code_clean = str_replace_all(as.character(code), "\\.", "")) %>%
  filter(str_detect(code_clean, "I48[0129]|42731")) %>%
  group_by(`Participant ID`) %>%
  summarise(pure_af_date = min(event_date, na.rm = TRUE)) %>%
  ungroup() %>%
  # 防止 min 返回 Inf
  mutate(pure_af_date = if_else(is.infinite(as.numeric(pure_af_date)), as.Date(NA), as.Date(pure_af_date)))

# 2.2 准备日期补丁 (防止 imp_long 缺列导致计算报错)
# 从你清洗好的 df_clinical_final 或 df_raw 中提取基线和死亡日期
date_patch <- df_clinical_final_pureAF %>% 
  select(`Participant ID`, baseline_dt, date_death) %>% 
  distinct() %>%
  mutate(
    baseline_dt = as.Date(baseline_dt),
    date_death = as.Date(date_death)
  )


# =======================================================
# 4. 处理插补数据：剔除 + 覆盖变量
# =======================================================

# 4.1 提取长格式并补齐日期
imp_long <- complete(imp_ready_for_cox, action = "long", include = TRUE) %>%
  # 先移除可能存在的残缺日期列，防止 join 冲突
  select(-any_of(c("baseline_dt", "date_death"))) %>%
  # 重新补入完整的日期基准
  left_join(date_patch, by = "Participant ID")

# 4.2 识别急性肝损 (基于 .imp==0 的原始值)
ids_acute_liver <- imp_long %>%
  filter(.imp == 0) %>%
  filter(`Alanine aminotransferase | Instance 0` > 400 | 
           `Aspartate aminotransferase | Instance 0` > 400) %>%
  pull(`Participant ID`) %>% unique()

# 汇总剔除名单
ids_to_exclude <- unique(c(ids_heme_final, ids_acute_liver))

# 4.3 [核心步骤] 过滤并重写结局变量
imp_long_final <- imp_long %>%
  
  
  # B. 关联纯房颤日期
  left_join(df_pure_af_dates, by = "Participant ID") %>%
  
  mutate(
    # --- C. 覆盖原有的 event_af_updated ---
    # 定义：有纯房颤日期 且 日期晚于基线 且 日期早于截止
    event_af_updated = if_else(
      !is.na(pure_af_date) & pure_af_date > baseline_dt & pure_af_date <= follow_up_end_dt, 
      1, 0, missing = 0
    ),
    
    # --- D. 覆盖原有的 duration_updated ---
    # 逻辑：发病算到发病日，死掉算到死亡日，其余算到截止日
    duration_updated = case_when(
      event_af_updated == 1 ~ as.numeric(pure_af_date - baseline_dt) / 365.25,
      !is.na(date_death) & date_death < follow_up_end_dt ~ as.numeric(date_death - baseline_dt) / 365.25,
      TRUE ~ as.numeric(follow_up_end_dt - baseline_dt) / 365.25
    )
  ) %>%
  
  # E. 数据清洗：移除负随访时间或极短时间 (反向因果：排除2年内发病)
  filter(duration_updated > 0) 

# =======================================================
# 5. 验证与保存
# =======================================================

# 5.1 检查事件数 (确保不是 0)
n_events <- sum(imp_long_final$event_af_updated[imp_long_final$.imp == 1], na.rm = TRUE)
n_sample <- nrow(imp_long_final %>% filter(.imp == 1))

cat("--- 最终数据集构建完成 ---\n")
cat("分析总样本量 (单套):", n_sample, "\n")
cat("新定义(纯房颤)事件数:", n_events, "\n")


# 5.2 转换为 mids 对象
imp_expo_clean <- as.mids(imp_long_final)

# 5.3 保存 (这个文件现在可以直接进 Cox 模型，变量名没变)
saveRDS(imp_expo_clean, "mice_sensi_pureAF.rds")


library(mice)
library(survival)
library(broom)
library(dplyr)
library(stringr)
library(flextable)
library(officer)
imp_ready_for_cox<-mice_sensi_pureAF
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

print(doc, target = "Liver_纯房颤_Analysis_WithCounts.docx")




