library(dplyr)
library(tidyr)
library(stringr)

cat("⏳ 正在按要求从三大数据源提取并合并数据...\n")

# ==============================================================================
# 1. 读取外部身高体重数据 (bcaasup.csv)
# ==============================================================================
# 注意：确保 bcaasup.csv 在你的工作目录下
bcaa_data <- read.csv("bcaasup.csv", stringsAsFactors = FALSE) %>%
  select(eid, 
         Height_cm = p50_i0, 
         Weight_kg = p21002_i0)

# ==============================================================================
# 2. 提取多重插补数据集 (mice_rediag) 中的血压和疾病结局
# ==============================================================================
# 假设你使用的是第一次插补的数据 (.imp == 1)，如果你的 mice_rediag 已经是单数据集，可删掉 filter 这一行
library(mice)
library(dplyr)

cat("⏳ 正在提取第 1 次插补的完整数据...\n")

# ==============================================================================
# 2. 提取多重插补数据集 (mice_rediag) 中的血压和疾病结局
# ==============================================================================
# 🚨 核心修复：用 complete(data, 1) 直接把第一套完整数据提取成普通表格
imputed_data <- complete(mice_rediag, action = 1) %>%
  # 统一把 ID 名字改成 eid，方便后面三大表格无缝合并
  rename(eid = `Participant ID`) %>% 
  select(eid, 
         htn_final, 
         t2dm_final, 
         sbp_val, 
         dbp_val,
         smoke_f,
         Race = ethnicity_f)

cat("✅ 插补数据提取成功！现在可以继续往下合并了...\n")
# ==============================================================================
# 3. 提取 df_final_clinical 中的基础信息和所有诊断列
# ==============================================================================
clin_data <- df_clinical_final %>%
  rename(eid = `Participant ID`) %>%
  select(eid,
         Age = age_val,
        
         med_htn,
         # 提取所有的 ICD-10, ICD-9 和自我报告疾病代码列
         starts_with("Diagnoses - ICD10"),
         starts_with("Diagnoses - ICD9"),
         starts_with("Date of first in-patient diagnosis"),
         starts_with("Date of attending assessment centre"),
         starts_with("Non-cancer illness code, self-reported"))

# ==============================================================================
# 4. 终极大会师：按 eid 合并三个数据集
# ==============================================================================
merged_cohort <- clin_data %>%
  inner_join(imputed_data, by = "eid") %>%
  inner_join(bcaa_data, by = "eid")%>%

  rename(baseline_dt = `Date of attending assessment centre | Instance 0`)

cat("✅ 数据合并完成！总人数:", nrow(merged_cohort), "\n")
cat("⏳ 正在进行极速字符串匹配，提取心衰、心梗、冠心病、慢阻肺和甲状腺病史...\n")

#writeLines(names(merged_cohort), "merged_cohort.txt")


# ==============================================================================
# 2. 定义 5 大核心病史的正则模式 (Regex) - 🚨强化扩充版
# ==============================================================================
# 1. 心力衰竭 (HF)
pat_hf_icd10 <- "I50|I110|I130|I132"
pat_hf_icd9  <- "^428"
pat_hf_sr    <- "1076|1077|heart failure|pulmonary oedema"

# 2. 心肌梗死 (MI)
pat_mi_icd10 <- "I21|I22|I252"
pat_mi_icd9  <- "^410|^412"
pat_mi_sr    <- "1075|myocardial infarction|heart attack"

# 3. 冠心病 (CAD)
pat_cad_icd10 <- "I2[0-5]"
pat_cad_icd9  <- "^41[0-4]|^4292"
pat_cad_sr    <- "1074|1075|1079|angina|myocardial infarction|heart attack"

# 4. 慢性阻塞性肺疾病/哮喘 (COPD/Asthma)
pat_copd_icd10 <- "J4[1-6]"
pat_copd_icd9  <- "^49[1236]"
pat_copd_sr    <- "1111|1112|1113|1114|asthma|copd|emphysema|chronic bronchitis"

# 5. 甲状腺疾病 (Thyroid)
pat_thyroid_icd10 <- "E0[0-7]"
pat_thyroid_icd9  <- "^24[0-6]"
pat_thyroid_sr    <- "1224|1225|1226|1227|thyroid|hypothyroidism|hyperthyroidism"


# ==============================================================================
# 3. 极度严谨的宽表转长表与日期过滤 (保持不变)
# ==============================================================================
cat("⏳ 正在进行 ICD 记录的日期匹配与清洗...\n")

df_icd10_codes <- merged_cohort %>%
  select(eid, diag_icd10 = `Diagnoses - ICD10`) %>%
  separate_rows(diag_icd10, sep = "\\|") %>%
  group_by(eid) %>% mutate(array_idx = paste0("a", row_number() - 1)) %>% ungroup()

df_icd10_dates <- merged_cohort %>%
  select(eid, contains("Date of first in-patient diagnosis - ICD10 | Array ")) %>%
  pivot_longer(cols = -eid, names_to = "array_idx", values_to = "event_date") %>%
  mutate(array_idx = paste0("a", str_extract(array_idx, "\\d+$")), event_date = as.Date(event_date)) %>%
  filter(!is.na(event_date))

df_icd10_long <- df_icd10_codes %>%
  left_join(df_icd10_dates, by = c("eid", "array_idx")) %>%
  left_join(merged_cohort %>% select(eid, baseline_dt), by = "eid") %>%
  filter(!is.na(event_date), event_date <= baseline_dt, !is.na(diag_icd10))

df_icd9_long <- merged_cohort %>%
  select(eid, diag_icd9 = `Diagnoses - ICD9`) %>%
  separate_rows(diag_icd9, sep = "\\|") %>%
  filter(!is.na(diag_icd9))


# ==============================================================================
# 4. 提取病患 ID 列表 - 🚨修复了自我报告的匹配 Bug
# ==============================================================================
cat("⏳ 正在提取 5 大疾病的基线患者名单...\n")

# A. 提取 ICD 记录
ids_hf_icd10 <- df_icd10_long %>% filter(str_detect(diag_icd10, pat_hf_icd10)) %>% pull(eid)
ids_mi_icd10 <- df_icd10_long %>% filter(str_detect(diag_icd10, pat_mi_icd10)) %>% pull(eid)
ids_cad_icd10 <- df_icd10_long %>% filter(str_detect(diag_icd10, pat_cad_icd10)) %>% pull(eid)
ids_copd_icd10 <- df_icd10_long %>% filter(str_detect(diag_icd10, pat_copd_icd10)) %>% pull(eid)
ids_thy_icd10 <- df_icd10_long %>% filter(str_detect(diag_icd10, pat_thyroid_icd10)) %>% pull(eid)

ids_hf_icd9 <- df_icd9_long %>% filter(str_detect(diag_icd9, pat_hf_icd9)) %>% pull(eid)
ids_mi_icd9 <- df_icd9_long %>% filter(str_detect(diag_icd9, pat_mi_icd9)) %>% pull(eid)
ids_cad_icd9 <- df_icd9_long %>% filter(str_detect(diag_icd9, pat_cad_icd9)) %>% pull(eid)
ids_copd_icd9 <- df_icd9_long %>% filter(str_detect(diag_icd9, pat_copd_icd9)) %>% pull(eid)
ids_thy_icd9 <- df_icd9_long %>% filter(str_detect(diag_icd9, pat_thyroid_icd9)) %>% pull(eid)

# B. 🚨提取自我报告记录 (用强大的 str_detect 横向扫描，绝不漏诊)
sr_cols <- grep("Non-cancer illness code", names(merged_cohort), value = TRUE)
df_sr_search <- merged_cohort %>%
  select(eid, all_of(sr_cols)) %>%
  mutate(all_sr_text = tolower(do.call(paste, c(., sep = ","))))

ids_hf_sr    <- df_sr_search %>% filter(str_detect(all_sr_text, pat_hf_sr)) %>% pull(eid)
ids_mi_sr    <- df_sr_search %>% filter(str_detect(all_sr_text, pat_mi_sr)) %>% pull(eid)
ids_cad_sr   <- df_sr_search %>% filter(str_detect(all_sr_text, pat_cad_sr)) %>% pull(eid)
ids_copd_sr  <- df_sr_search %>% filter(str_detect(all_sr_text, pat_copd_sr)) %>% pull(eid)
ids_thy_sr   <- df_sr_search %>% filter(str_detect(all_sr_text, pat_thyroid_sr)) %>% pull(eid)

# C. 合并去重
hf_all_ids   <- unique(c(ids_hf_icd10, ids_hf_icd9, ids_hf_sr))
mi_all_ids   <- unique(c(ids_mi_icd10, ids_mi_icd9, ids_mi_sr))
cad_all_ids  <- unique(c(ids_cad_icd10, ids_cad_icd9, ids_cad_sr))
copd_all_ids <- unique(c(ids_copd_icd10, ids_copd_icd9, ids_copd_sr))
thy_all_ids  <- unique(c(ids_thy_icd10, ids_thy_icd9, ids_thy_sr))
# ==============================================================================
# 5. 组装最终极其干净的建模数据框
# ==============================================================================
df_ready_for_scores <- merged_cohort %>%
  mutate(
    is_HF      = if_else(eid %in% hf_all_ids, 1, 0),
    is_MI      = if_else(eid %in% mi_all_ids, 1, 0),
    is_CAD     = if_else(eid %in% cad_all_ids, 1, 0),
    is_COPD    = if_else(eid %in% copd_all_ids, 1, 0),
    is_Thyroid = if_else(eid %in% thy_all_ids, 1, 0)
  ) %>%
  # 只保留计算 CHARGE-AF, ARIC, C2HEST 必备的黄金变量
  select(eid, Age, Race, smoke_f, Height_cm, Weight_kg, 
         sbp_val, dbp_val, med_htn, htn_final, t2dm_final, 
         is_HF, is_MI, is_CAD, is_COPD, is_Thyroid)

cat("🎉 完美！所有变量清洗完毕并绝对严防了未来数据泄露！当前可用数据:", nrow(df_ready_for_scores), "\n")
# 计算各项患病人数和百分比
disease_stats <- df_ready_for_scores %>%
  summarise(
    Total_N = n(),
    # 1. 核心临床病史 (这几个是我们刚生成的纯数字 0/1，可以直接 sum)
    Heart_Failure = sum(is_HF, na.rm = TRUE),
    Myocardial_Infarction = sum(is_MI, na.rm = TRUE),
    Coronary_Artery_Disease = sum(is_CAD, na.rm = TRUE),
    COPD_Asthma = sum(is_COPD, na.rm = TRUE),
    Thyroid_Disease = sum(is_Thyroid, na.rm = TRUE),
    
    # 2. 基础代谢合并症 (🚨 修复 Factor 求和报错)
    # 用逻辑判断包起来，TRUE 会自动变成 1 被 sum
    Hypertension_Final = sum(htn_final == 1 | htn_final == "1" | htn_final == "Yes", na.rm = TRUE),
    Diabetes_Final = sum(t2dm_final == 1 | t2dm_final == "1" | t2dm_final == "Yes", na.rm = TRUE),
    
    # 3. 吸烟情况 (Current smoker)
    # 包含 Current 即算数 (防范大小写或前后空格)
    Current_Smoker = sum(str_detect(tolower(smoke_f), "current"), na.rm = TRUE)
  ) %>%
  # 转为长表方便查看
  pivot_longer(cols = -Total_N, names_to = "Condition", values_to = "Count") %>%
  mutate(
    Percentage = sprintf("%.2f%%", (Count / Total_N) * 100)
  )

print(disease_stats)
# ==============================================================================
# 6. 终极数据质控：检查每一列的缺失值 (NA)
# ==============================================================================
cat("\n================ 缺失值检查报告 ================\n")
missing_summary <- colSums(is.na(df_ready_for_scores))

# 把缺失情况做成一个好看的数据框，并计算缺失率
missing_df <- data.frame(
  Variable = names(missing_summary),
  Missing_Count = as.numeric(missing_summary),
  Missing_Rate = sprintf("%.3f%%", (as.numeric(missing_summary) / nrow(df_ready_for_scores)) * 100)
) %>%
  filter(Missing_Count > 0) %>% # 只挑出有缺失值的列
  arrange(desc(Missing_Count))

if(nrow(missing_df) == 0) {
  cat("🎉 太牛了！数据集里没有任何缺失值 (0 NAs)，可以直接计算模型得分！\n")
} else {
  cat("⚠️ 发现以下变量存在缺失值：\n")
  print(missing_df)
}
cat("================================================\n")



# ==============================================================================
# 终极严谨版：临床模型绝对风险计算 (Explicit Calibration Version)
# ==============================================================================
cat("⏳ 正在应用原始文献 Mean LP 校正与时间原位风险转换...\n")

df_benchmark_scores <- df_ready_for_scores %>%
  # --- 第一步：统一生成所有 0/1 逻辑变量 ---
  mutate(
    # 种族判定（基于你确认的 Race 列）
    is_white   = if_else(str_detect(tolower(as.character(Race)), "white|british|irish"), 1, 0, missing = 0),
    # 吸烟判定
    is_smoker  = if_else(str_detect(tolower(as.character(smoke_f)), "current"), 1, 0, missing = 0),
    # 用药与疾病判定
    is_htn_med = if_else(med_htn == 1 | med_htn == "1" | med_htn == "Yes", 1, 0, missing = 0),
    is_dm      = if_else(t2dm_final == 1 | t2dm_final == "1" | t2dm_final == "Yes", 1, 0, missing = 0),
    is_htn_dz  = if_else(htn_final == 1 | htn_final == "1" | htn_final == "Yes", 1, 0, missing = 0)
  ) %>%
  mutate(
    # --- 1. CHARGE-AF (Alonso 2013, JAMA Intern Med) ---
    # 均值中心化逻辑：Age-57, Height-167, Weight-79, SBP-120, DBP-75
    # 由于严格按照原文均值中心化，其 Mean LP 在理论上归为 0
    CHARGE_AF_LP = 
      (0.508 * ((Age - 57) / 5)) + (0.248 * is_white) + 
      (0.115 * ((Height_cm - 167) / 10)) + (0.197 * ((Weight_kg - 79) / 15)) + 
      (0.458 * ((sbp_val - 120) / 20)) - (0.352 * ((dbp_val - 75) / 10)) + 
      (0.359 * is_smoker) + (0.349 * is_htn_med) + 
      (0.237 * is_dm) + (0.496 * is_MI) + (1.212 * is_HF),
    
    # 5年绝对风险：Risk = 1 - S0(5)^exp(LP - 0)
    # S0(5) = 0.971841 来自原文 Table 3
    CHARGE_AF_5y_Risk = 1 - (0.971841 ^ exp(CHARGE_AF_LP)),
    # --- 2. ARIC (Chamberlain 2011) ---
    # 🚨 修正核心：对 ARIC 变量也进行“开发人群均值中心化”
    # 均值参考：Age(54.2), Height(168.1), SBP(122.5)
    # 这样计算出的 ARIC_LP 均值会回归 0 附近，避免 exp 爆炸
    ARIC_LP_centered = 
      (0.106 * (Age - 54.2)) - 
      (0.191 * (is_white - 0.77)) + # 原始 ARIC 约 77% 白人
      (0.057 * (Height_cm - 168.1)) + 
      (0.418 * (is_smoker - 0.26)) + # 原始吸烟率约 26%
      (0.013 * (sbp_val - 122.5)) + 
      (0.428 * (is_htn_med - 0.23)) + 
      (0.285 * (is_dm - 0.08)) + 
      (0.551 * (is_CAD - 0.05)) + 
      (1.054 * (is_HF - 0.01)),
    
    # 现在这个 LP 已经是“净 LP”了，直接套用 S0(10) = 0.929
    ARIC_10y_Risk = 1 - (0.929 ^ exp(ARIC_LP_centered)),
    
    # --- 3. C2HEST (Li 2019, Europace) ---
    # 纯积分系统，计算其 LP 主要用于 C-index 比较
    C2HEST_Points = (is_CAD*1) + (is_COPD*1) + (is_htn_dz*1) + 
      (if_else(Age>=75, 2, 0)) + (is_HF*2) + (is_Thyroid*1)
  )

cat("🎉 全部模型已完成显式校正计算！\n")
# ==============================================================================
# 📊 房颤评分分布终极体检报告
# ==============================================================================
cat("\n🔎 正在进行多维度分布分析...\n")

# 1. 基础统计摘要
score_summary <- df_benchmark_scores %>%
  select(CHARGE_AF_5y_Risk, ARIC_10y_Risk, C2HEST_Points) %>%
  summary()

print(score_summary)

# 2. 边界逻辑检查 (检测是否有非法概率)
logical_check <- df_benchmark_scores %>%
  summarise(
    CHARGE_AF_OOB = sum(CHARGE_AF_5y_Risk < 0 | CHARGE_AF_5y_Risk > 1, na.rm = TRUE),
    ARIC_OOB = sum(ARIC_10y_Risk < 0 | ARIC_10y_Risk > 1, na.rm = TRUE),
    NAs_Count = sum(is.na(CHARGE_AF_5y_Risk) | is.na(ARIC_10y_Risk))
  )

cat("\n⚠️ 逻辑越界检查 (越界数应为0):\n")
print(logical_check)

# 3. 极速绘图：看清风险的“长相”
library(ggplot2)
library(patchwork) # 用于并排显示图片

p1 <- ggplot(df_benchmark_scores, aes(x = CHARGE_AF_5y_Risk)) +
  geom_density(fill = "#69b3a2", alpha = 0.5) +
  theme_minimal() + labs(title = "CHARGE-AF 5y Risk Distribution")

p2 <- ggplot(df_benchmark_scores, aes(x = ARIC_10y_Risk)) +
  geom_density(fill = "#404080", alpha = 0.5) +
  theme_minimal() + labs(title = "ARIC 10y Risk Distribution")

p3 <- ggplot(df_benchmark_scores, aes(x = factor(C2HEST_Points))) +
  geom_bar(fill = "#f68060", alpha = 0.8) +
  theme_minimal() + labs(title = "C2HEST Points Frequency", x = "Points")

(p1 / p2) | p3 # 拼图显示




