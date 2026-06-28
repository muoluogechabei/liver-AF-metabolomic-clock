


# 1. 载入所有神级装备包
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(mclust)
library(uwot)           
library(mclustcomp)     
library(ggpointdensity) 
library(survminer)      
library(ClusterR)       
library(caret)          
library(pheatmap)       

# ==============================================================================
# 🛡️ 步骤 1 & 2：全量数据准备与 70% 训练集提取 (极简精准化)
# ==============================================================================
# 🚨 核心修改 1：直接使用刚出炉的 final_11_metabs (8个王者代谢物)
final_metabs <- final_11_metabs

# 把 7:3 拆分的两个集合拼起来作为 100% 全景图底座 (用于最后 UMAP 投影)
# 🚨 核心修改 2：适配新的 eid，剔除了不需要的冗余列防报错
full_cohort <- bind_rows(df_train_final, df_test_final) %>%
  select(eid, duration_updated, event_af_updated, all_of(final_metabs)) %>%
  drop_na()


# 推荐加载 data.table 和 dplyr，处理 UKB 大数据的黄金组合
library(data.table)
library(dplyr)

# 1. 定义你要从大表里抓取的所有原始列名（加上了身高和体重）
target_cols <- c(
  "eid",
  "p50_i2", "p21002_i2",    # Instance 2 的身高和体重
  "p22420_i2", "p22420_i3", # LV ejection fraction
  "p24105_i2", "p24105_i3", # LV myocardial mass
  "p24110_i2", "p24110_i3", # LA maximum volume
  "p24181_i2", "p24181_i3", # LV longitudinal strain global
  "p24157_i2", "p24157_i3", # LV circumferential strain global
  "p24106_i2", "p24106_i3", # RV end diastolic volume
  "p24113_i2", "p24113_i3" , # LA ejection fraction
  "p21003_i2"
)

# 2. 高效读取（只读需要的列，极大节省内存）
df_raw <- fread("databcaa.csv", select = target_cols)

# 3. 临床变量重命名 (The Magic Step)
# 把反人类的编号变成极具临床意义的简写，方便后续写公式
df_cmr_base <- df_raw %>%
  rename(
    Height_i2 = p50_i2,
    Weight_i2 = p21002_i2,
    
    LVEF_i2 = p22420_i2,         LVEF_i3 = p22420_i3,
    LVMass_i2 = p24105_i2,       LVMass_i3 = p24105_i3,
    LA_MaxVol_i2 = p24110_i2,    LA_MaxVol_i3 = p24110_i3,
    LV_GLS_i2 = p24181_i2,       LV_GLS_i3 = p24181_i3,
    LV_GCS_i2 = p24157_i2,       LV_GCS_i3 = p24157_i3,
    RVEDV_i2 = p24106_i2,        RVEDV_i3 = p24106_i3,
    LAAEF_i2 = p24113_i2,        LAAEF_i3 = p24113_i3
  ) %>%
  # 4. 顺手把 BSA 和临床核心指数直接算出来！
  mutate(
    # 使用 Mosteller 公式计算体表面积 (m2)，身高单位需为cm，体重为kg
    BSA_i2 = sqrt((Height_i2 * Weight_i2) / 3600),
    
    # 计算心内科最看重的两大指数：LAVI 和 LVMI
    LAVI_i2 = LA_MaxVol_i2 / BSA_i2,
    LVMI_i2 = LVMass_i2 / BSA_i2
  )

# 5. 查看完美成型的底座数据
print(dim(df_cmr_base))
glimpse(df_cmr_base)
library(data.table)
library(dplyr)

# =========================================================================
# 步骤 1：读取师兄的数据，精准匹配带点号的列名
# =========================================================================
# 注意：UKB 导出时，'|' 往往变成了 '...'，空格变成了 '.'
shixiong_cols <- c(
  "Participant ID", 
  "Body surface area | Instance 2", 
  "Date of attending assessment centre | Instance 2",
  "LV end diastolic volume | Instance 2", 
  "LV end systolic volume | Instance 2",  
  "LV stroke volume | Instance 2",
  "RA maximum volume | Instance 2"
)

# 读取并立即重命名为干净的变量名
df_shixiong <- fread("cardiac mri pertinent index.csv", select = shixiong_cols) %>%
  rename(
    eid = `Participant ID`, # 将师兄的 ID 对齐到你的 eid
    BSA_official_i2 =`Body surface area | Instance 2`,
    Date_MRI_i2 = `Date of attending assessment centre | Instance 2`,
    LVEDV_i2 =`LV end diastolic volume | Instance 2`,
    LVESV_i2 = `LV end systolic volume | Instance 2`,
    LVSV_i2 = `LV stroke volume | Instance 2`,
    RAMV_i2 = `RA maximum volume | Instance 2`
  )

# =========================================================================
# 步骤 2：合并并计算 Indexed 临床指标
# =========================================================================
df_mechanism_structural <- df_cmr_base %>%
  left_join(df_shixiong, by = "eid") %>%
  mutate(
    # 1. 处理日期格式
    Date_MRI_i2 = as.Date(Date_MRI_i2),
    
    # 2. 计算心内科核心指数 (使用师兄提供的官方 BSA)
    # LA_MaxVol_i2 是你自己从 databcaa 提出来的 p24110_i2
    # LVMass_i2 是你提出来的 p24105_i2
    LAVI_i2 = LA_MaxVol_i2 / BSA_official_i2,     # 左房容积指数
    LVMI_i2 = LVMass_i2 / BSA_official_i2,
    RAVI_i2 = RAMV_i2 / BSA_official_i2,# 左室质量指数
    LVEDVI_i2 = LVEDV_i2 / BSA_official_i2     # 左室舒张末期容积指数
  )

# =========================================================================
# 步骤 3：验证合并是否成功
# =========================================================================
# 看看非空值 (NA) 的比例，如果数量对得上说明 join 成功了
summary(df_mechanism_structural$LAVI_i2)


library(data.table)
library(dplyr)

# =========================================================================
# 步骤 1：读取 xdt.csv 并处理列名
# =========================================================================
# 假设 xdt.csv 里的 ID 也是 "Participant ID"，我们先读取并改名
# 如果它里面本来就叫 "eid"，请把 rename 里的 `Participant ID` 改成 eid
df_xdt <- fread("ECG.csv") 

# 自动兼容 ID 命名（如果是 Participant ID 就改名，如果已经是 eid 就保持）
if ("Participant ID" %in% colnames(df_xdt)) {
  df_xdt <- df_xdt %>% rename(eid = `Participant ID`)
}

# =========================================================================
# 步骤 2：全量大合并 (Final Join)
# =========================================================================
# 将之前算好的结构数据与 xdt 中的电生理等数据合并
df_final_all <- df_mechanism_structural %>%
  left_join(df_xdt, by = "eid")

# =========================================================================
# 步骤 3：输出所有列名到 TXT 文件，方便“查漏补缺”
# =========================================================================
# 获取所有列名
all_colnames <- colnames(df_final_all)

# 写入 txt 文件
writeLines(all_colnames, "all_column_names_check.txt")

# 在控制台打印摘要
cat(sprintf("✅ 大合拢完成！最终数据集共包含 %d 行, %d 个变量。\n", 
            nrow(df_final_all), ncol(df_final_all)))
cat("📊 所有列名已导出至：all_column_names_check.txt，请打开该文件核对。\n")

# =========================================================================
# 步骤 4：核心变量存活检查 (快速预览)
# =========================================================================
# 你可以运行下面这一行，看看最核心的几块拼图在不在
core_check <- c("eid", "LAVI_i2", "LVMI_i2", "BSA_official_i2", "Date_MRI_i2")
print(df_final_all %>% select(any_of(core_check)) %>% head())
saveRDS(df_final_all,"ECG.rds")

library(data.table)
fwrite(df_final_all, "ECG.csv")









library(dplyr)
library(stringr)
library(tidyr)
library(data.table)

cat("\n=======================================================\n")
cat("🫀 极速版：构建以 Instance 2 (MRI) 为锚点的纯净人群与用药矩阵\n")
cat("=======================================================\n")

# 假设 df_clinical_final 已经读入环境中
# 先将 ID 标准化为 eid，并提取核磁日期
df_target_base <- df_clinical_final %>%
  rename(
    eid = `Participant ID`,
    Date_MRI_i2 = `Date of attending assessment centre | Instance 2`
  ) %>%
  # 必须要有核磁日期的人才有资格进入我们的分析
  filter(!is.na(Date_MRI_i2)) %>%
  mutate(Date_MRI_i2 = as.Date(Date_MRI_i2))

# =========================================================================
# 💊 步骤 1：构建用药矩阵 (扫描 Instance 0 和 Instance 2 的所有用药列)
# =========================================================================
cat("正在扫描 Instance 0 和 Instance 2 的服药记录...\n")

# 定义正则表达式
regex_htn_meds <- "amlodipine|ramipril|bisoprolol|1140866148|1140860806|1140883446|1140860696|1140888746|1140860808|1140867878"
regex_lip_meds <- "statin|atorvastatin|simvastatin|1140888594|1140861958|1141146234|1140861970|1141192410|1140881748"
regex_dm_meds  <- "metformin|gliclazide|insulin|1140884600|1140883066|1140884646|1141152732|1140888648|1140874686"

# 将所有用药列合并成一个长字符串，方便正则匹配 (极其提速的写法)
df_meds <- df_target_base %>%
  select(eid, starts_with("Treatment/medication code | Instance 0")) %>%
  # 将一个人的所有用药代码拼接在一起，忽略 NA
  unite("all_meds", -eid, sep = "|", na.rm = TRUE) %>%
  mutate(
    meds_htn_i0 = ifelse(str_detect(all_meds, regex_htn_meds), 1, 0),
    meds_lip_i0 = ifelse(str_detect(all_meds, regex_lip_meds), 1, 0),
    meds_dm_i0  = ifelse(str_detect(all_meds, regex_dm_meds), 1, 0)
  ) %>%
  select(eid, meds_htn_i0, meds_lip_i0, meds_dm_i0)

# =========================================================================
# 🛡️ 步骤 2：致命疾病排阴 (剔除 Instance 2 之前的房颤与心衰/冠心病)
# =========================================================================
cat("正在根据核磁日期 (Instance 2) 执行严格排阴逻辑...\n")

df_exclusion <- df_target_base %>%
  select(eid, Date_MRI_i2, is_prevalent_af, af_date, baseline_cvd, earliest_icd_dt) %>%
  mutate(
    af_date = as.Date(af_date),
    earliest_icd_dt = as.Date(earliest_icd_dt),
    
    # 1. 判定在做核磁前是否已有房颤 (AF)
    # 如果 baseline 时已有 (is_prevalent_af == 1) 或 af_date 早于核磁当天，则标记为 TRUE
    has_af_before_mri = case_when(
      is_prevalent_af == 1 ~ TRUE,
      !is.na(af_date) & af_date <= Date_MRI_i2 ~ TRUE,
      TRUE ~ FALSE
    ),
    
    # 2. 判定在做核磁前是否已有严重心血管并发症 (CVD/HF)
    # 利用你现成的 baseline_cvd，以及如果 Instance 0~2 之间新发 CVD 也算
    has_cvd_before_mri = case_when(
      baseline_cvd == 1 ~ TRUE,
      # 假设你原始的 earliest_icd_dt 捕捉了 CVD 首次确诊时间
      !is.na(earliest_icd_dt) & earliest_icd_dt <= Date_MRI_i2 ~ TRUE, 
      TRUE ~ FALSE
    )
  )

# =========================================================================
# 🎯 步骤 3：输出最终纯净的人群 ID 与用药矩阵
# =========================================================================
df_clean_cohort_i2 <- df_target_base %>%
  select(eid, Date_MRI_i2) %>%
  left_join(df_meds, by = "eid") %>%
  left_join(df_exclusion %>% select(eid, has_af_before_mri, has_cvd_before_mri), by = "eid") %>%
  
  # 执行核心剔除：只要核磁当天干干净净的人！
  filter(has_af_before_mri == FALSE & has_cvd_before_mri == FALSE) %>%
  
  # 丢掉排阴标记列，只保留我们需要的：ID, 核磁日期，用药矩阵
  select(eid, Date_MRI_i2, meds_htn_i0, meds_lip_i0, meds_dm_i0)

cat("\n✅ 清洗完成！\n")
cat(sprintf("经过严格的 Instance 2 日期排阴，共锁定 %d 名纯净受试者。\n", nrow(df_clean_cohort_i2)))
cat("生成的 df_clean_cohort_i2 包含列名：\n")
print(colnames(df_clean_cohort_i2))

# 下一步：拿这个 df_clean_cohort_i2 的 eid 去你的多重插补数据集里提取连续变量（血压等）
library(dplyr)
library(mice)

cat("\n=======================================================\n")
cat("🎯 终极拼图：利用纯净 ID 从插补数据集中捞取协变量底座\n")
cat("=======================================================\n")

# 假设你之前跑出来的纯净人群表叫做 df_clean_cohort_i2，里面有 eid 和 meds_xxx_i2
# 假设你已经 readRDS 载入了 mice_rediag (Instance 0 插补) 和 imp_ready_for_cox (里面有 i2 年龄)

# =========================================================================
# 步骤 1：从 mice_rediag 中提取 Instance 0 协变量 (取第 1 次插补结果)
# =========================================================================
# 提取插补集 1 (如果想更严谨，可以用 pool 逻辑，但在这种单次机制验证中，用 dataset 1 足矣)
df_cov_i0 <- complete(mice_rediag, 1) %>%
  # 将你原本用来跑 mice 的 ID 名字对齐为 eid
  rename(eid = `Participant ID`) %>%
  select(
    eid,
    sex_f,                      # 性别 (恒定)
    bmi_val,                    # Instance 0 BMI (用于不带 BSA 的模型)
    sbp_val,                    # Instance 0 血压 (插补后的干净血压)
    smoke_f,                    # 抽烟
    alc_status_f                # 饮酒状态
    # 注意：你说降压药问题你处理好了，所以我这里没捞 htn_final，
    # 假设你会在回归代码里直接用你之前生成的 meds_htn_i0 (如果不在里面，请在这里加上 htn_final)
  )

# =========================================================================
# 步骤 2：从另一个插补集 imp_ready_for_cox 中提取 Instance 2 的年龄
# =========================================================================
# 提取插补集 1
df_cov_i2_age <- df_raw %>%
  # 这里你明确说过 ID 叫 eid，年龄叫 age
  select(
    eid,
    age_i2 = p21003_i2  # 重命名为 age_i2 避免与基线年龄混淆
  )

# =========================================================================
# 步骤 3：三表合一 (The Final Merge)
# =========================================================================
# 以你筛选出的纯净 ID (df_clean_cohort_i2) 为基准，向左合并所有协变量
df_mechanism_analysis <- df_clean_cohort_i2 %>%
  # 1. 拼入 I0 底座
  left_join(df_cov_i0, by = "eid") %>%
  # 2. 拼入 I2 年龄
  left_join(df_cov_i2_age, by = "eid") 

# =========================================================================
# 步骤 4：检查终极分析矩阵
# =========================================================================
cat("\n✅ 终极分析矩阵组装完毕！\n")
cat(sprintf("总样本量：%d 人\n", nrow(df_mechanism_analysis)))
cat("包含的核心协变量：\n")
print(colnames(df_mechanism_analysis %>% select(age_i2, sex_f, bmi_val, sbp_val, meds_htn_i0)))

# 下一步：全部 scale() 并启动批量回归！

library(dplyr)
# =========================================================================
# 🚀 步骤 5：填充代谢物“子弹” (对接筛选后的最新名单)
# =========================================================================
# 1. 从 full_cohort（之前算好 MetRS 的表）中提取最新的王者代谢物
# 🚨 注意：请检查 full_cohort 里的 ID 列名，如果是 eid 就用 eid，如果是 Participant.ID 就改一下
df_metabs_to_join <- full_cohort %>%
  select(eid, all_of(final_metabs)) # 使用刚才定义的最新 final_metabs

# 2. 与分析底座合并
df_mechanism_step5 <- df_mechanism_analysis %>%
  left_join(df_metabs_to_join, by = "eid") %>%
  drop_na(all_of(final_metabs))

cat(sprintf("✅ 代谢物装填完毕！当前分析矩阵包含 %d 人, %d 个代谢物。\n", 
            nrow(df_mechanism_step5), length(final_metabs)))
cat(sprintf("当前剩余样本量：%d 人\n", nrow(df_mechanism_step5)))
cat(sprintf("已装入代谢物数量：%d 个\n", length(final_metabs)))

# 4. 预览一眼（看看前 5 个代谢物和核心协变量）
glimpse(df_mechanism_step5 %>% select(eid, age_i2, sbp_val, any_of(final_metabs[1:5])))

library(dplyr)
library(data.table)

cat("\n=======================================================\n")
cat("🏁 终极合龙：合并机制分析底座与心脏指标\n")
cat("=======================================================\n")

# 1. 执行最终合并
# 以 df_mechanism_step5 为主（因为它已经选好了纯净 ID、协变量和代谢物）
df_mechanism_final <- df_mechanism_step5 %>%
  left_join(df_final_all, by = "eid")

# 2. 导出所有列名到 txt 文件，方便你最后一次“总点名”
final_columns <- colnames(df_mechanism_final)
writeLines(final_columns, "final_mechanism_database_columns.txt")

# 3. 统计一下最终有多少人、多少个指标
n_obs <- nrow(df_mechanism_final)
n_vars <- ncol(df_mechanism_final)

cat(sprintf("✅ 合并成功！\n"))
cat(sprintf("📊 最终分析样本量：%d 人\n", n_obs))
cat(sprintf("📊 最终总变量数：%d 个\n", n_vars))
cat("📝 完整清单已保存至：final_mechanism_database_columns.txt\n")

# 4. 预览核心模块是否都在
cat("\n--- 核心模块检查 ---\n")
check_list <- list(
  Covariates = c("age_i2", "sex_f", "sbp_val", "meds_htn_i0"),
  Metabolites = head(final_metabs, 3), # 选前 3 个代谢物
  Cardiac_MRI = c("LAVI_i2", "LVMI_i2", "LVEF_i2","RAVI_i2"),
  Cardiac_ECG = c("PQ_i2", "QTc_i2", "HR_i2")
)
print(check_list)


library(dplyr)
library(tidyr)
library(broom)

cat("\n=======================================================\n")
cat("🔥 启动机制验证引擎：全自动分轨批量线性回归\n")
cat("=======================================================\n")

# =========================================================================
# 1. 清洗电生理列名 (防止带有空格和特殊符号的列名导致 lm() 报错)
# =========================================================================
df_regression <- df_mechanism_final %>%
  rename(
    HR_i2  = `Ventricular rate | Instance 2`,
    PR_i2  = `PQ interval | Instance 2`,
    QRS_i2 = `QRS duration | Instance 2`,
    QTc_i2 = `QTC interval | Instance 2`,
    pd_i2 = `P duration | Instance 2`
  )

# 修改为：
# 🚨 核心修改：直接对接 500 次 Bootstrap 选出的王者名单
final_metabs <- ultra_stable_metabs$Metabolite 
metabolites <- final_metabs

cat(sprintf("🚀 机制分析已同步：正在对 %d 个高频稳定代谢物进行回归...\n", length(metabolites)))
# B. 靶点分组 (根据你的策略，严格分流)
# 组别 1：已用 BSA 索引，坚决不加 BMI
targets_bsa <- c("LAVI_i2.x", "LVMI_i2.x","RAVI_i2.x", "LVEDVI_i2.x")

# 组别 2：未用 BSA 索引，必须加 BMI
targets_non_bsa <- c("LVEF_i2.x", "LAAEF_i2.x", "LV_GLS_i2.x", "HR_i2", "PR_i2", "QRS_i2", "QTc_i2","pd_i2")

# =========================================================================
# 3. 执行批量回归循环
# =========================================================================
results_list <- list()

for (metab in metabolites) {
  # 遍历所有靶点
  for (target in c(targets_bsa, targets_non_bsa)) {
    
    # 根据靶点类型，动态生成公式
    if (target %in% targets_bsa) {
      # 不加 BMI
      formula_str <- paste0("scale(", target, ") ~ scale(", metab, ") + age_i2 + sex_f + scale(sbp_val) + meds_htn_i0")
    } else {
      # 加 BMI
      formula_str <- paste0("scale(", target, ") ~ scale(", metab, ") + age_i2 + sex_f + scale(sbp_val) + meds_htn_i0 + scale(bmi_val)")
    }
    
    # 拟合模型
    fit <- lm(as.formula(formula_str), data = df_regression)
    
    # 提取代谢物的统计结果并加上【实际样本量】
    res <- tidy(fit) %>%
      filter(term == paste0("scale(", metab, ")")) %>%
      mutate(
        Metabolite = metab,
        Target = target,
        Model_Type = ifelse(target %in% targets_bsa, "No_BMI", "With_BMI"),
        Sample_Size = nobs(fit)  # 🌟 新增：提取该模型实际纳入的无缺失样本量！
      ) %>%
      select(Metabolite, Target, Model_Type, Sample_Size, estimate, std.error, p.value)
    
    results_list[[paste(metab, target, sep = "_")]] <- res
  }
}

# 合并所有结果
df_results <- bind_rows(results_list)

# =========================================================================
# 4. 后处理：方向翻转与 FDR 校正 (画图前必须做的预处理)
# =========================================================================
df_results_final <- df_results %>%
  mutate(
    # 🚨 极其关键的一步：
    # LVEF, LAAEF (射血分数) 越低越差；GLS (应变) 通常绝对值越低越差 (取决于UKB的正负号约定，通常需要翻转)。
    # 这里我们将 LVEF 和 LAAEF 的 Beta 乘以 -1。
    Final_Beta = ifelse(Target %in% c("LVEF_i2.x", "LAAEF_i2.x", "LVEF_i2", "LAAEF_i2"), -estimate, estimate),
    # 计算全局 FDR (Benjamini-Hochberg)
    FDR = p.adjust(p.value, method = "fdr"),
    
    # 按照显著性打星号 (画图用)
    Significance = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE ~ ""
    )
  ) %>%
  arrange(Target, FDR)

# =========================================================================
# 5. 查看战斗成果
# =========================================================================
cat("\n✅ 批量回归跑完了！\n")
cat(sprintf("共完成 %d 次线性回归模型拟合。\n", nrow(df_results_final)))
cat("\n前 10 个最显著的关联：\n")
print(head(df_results_final %>% select(Metabolite, Target, Final_Beta, FDR, Significance), 10))

library(dplyr)

cat("\n=======================================================\n")
cat("📊 各心脏靶点实际纳入回归的有效样本量审查\n")
cat("=======================================================\n")

df_sample_check <- df_results_final %>%
  group_by(Target, Model_Type) %>%
  summarise(
    # 因为我们之前已经把代谢物 (X) 强制 drop_na 补齐了，
    # 所以同一个 Target 对不同代谢物跑回归时，样本量应该是一模一样的
    Sample_Size = unique(Sample_Size), 
    .groups = "drop"
  ) %>%
  arrange(desc(Sample_Size))

print(df_sample_check)
library(readxl)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(grid)
library(RColorBrewer)
metab_dict<-read_excel("dict.xlsx")
# ——————————————————————————————————————————————
# 1. 保留缩写字典
# ——————————————————————————————————————————————
shorten_names <- function(x) {
  x <- gsub("\\.", " ", x) 
  # 🚨 补充的漏网之鱼：
  x <- sub("(?i).*polyunsaturated fatty acids to monounsaturated fatty acids.*", "PUFA/MUFA", x, perl=TRUE)
  
  x <- sub("(?i).*polyunsaturated fatty acids to total fatty acids.*", "PUFA% in Total FA", x, perl=TRUE)
  x <- sub("(?i).*monounsaturated fatty acids to total fatty acids.*", "MUFA% in Total FA", x, perl=TRUE)
  x <- sub("(?i).*docosahexaenoic acid to total fatty acids.*", "DHA% in Total FA", x, perl=TRUE)
  x <- gsub("(?i)Clinical LDL cholesterol", "Clinical LDL-C", x, perl=TRUE)
  x <- gsub("(?i)Glycoprotein acetyls", "GlycA", x, perl=TRUE)
  x <- gsub(" to total lipids ratio in ", "% in ", x)
  x <- gsub("Concentration of ", "", x)
  x <- gsub("Total concentration of ", "Total ", x)
  x <- gsub("Average diameter for ", "Diam: ", x)
  x <- gsub("Phospholipids", "PL", x)
  x <- gsub("Triglycerides", "TG", x)
  x <- gsub("Cholesteryl esters", "CE", x)
  x <- gsub("Free cholesterol", "FC", x)
  x <- gsub("Cholesterol", "Chol", x)
  x <- gsub("Lipoprotein", "Lipo", x)
  x <- gsub("particles", "P", x)
  x <- gsub("very large", "VL", x)
  x <- gsub("large", "L", x)
  x <- gsub("medium", "M", x)
  x <- gsub("very small", "VS", x)
  x <- gsub("small", "S", x)
  x <- gsub("Total fatty acids", "Total FA", x)
  x <- gsub("Ratio of ", "", x)
  x <- gsub("Glycolysis related metabolites", "Glycolysis", x)
  x <- gsub("Lipoprotein particle concentrations", "Lipo Concs", x)
  x <- gsub("Lipoprotein particle sizes", "Lipo Sizes", x)
  x <- gsub("Relative lipoprotein lipid concentrations", "Rel Lipo Concs", x)
  x <- gsub("Lipoprotein subclasses", "Lipo Subclasses", x)
  return(x)
}

# ——————————————————————————————————————————————
# 1. 字典映射与靶点重命名 (修正匹配逻辑版)
# ——————————————————————————————————————————————
df_plot <- df_results_final %>%
  mutate(
    # 🚨 关键：先将回归结果里的点号还原成空格，才能匹配字典
    Match_Name = gsub("\\.", " ", Metabolite)
  ) %>%
  # 关联字典：用还原后的名字去匹配 title
  left_join(metab_dict %>% select(title, Group), by = c("Match_Name" = "title")) %>%
  mutate(
    Group = ifelse(is.na(Group) | Group == "", "Other Metabolites", Group),
    # 统一应用缩写魔法
    Metabolite_Short = shorten_names(Match_Name),
    Group_Short = shorten_names(Group)
  )

# 1.2 心脏靶点重命名 (保持原样)
target_info <- data.frame(
  Target = c("LAVI_i2.x", "LAAEF_i2.x", "RAVI_i2.x", 
             "LVMI_i2.x", "LVEDVI_i2.x", "LVEF_i2.x", "LV_GLS_i2.x", 
             "HR_i2", "pd_i2", "PR_i2", "QRS_i2", "QTc_i2"),
  Pub_Name = c("LA Volume Index", "LA Ejection Fraction (Inv)", "RA Volume Index",
               "LV Mass Index", "LV EDV Index", "LV Ejection Fraction (Inv)", "LV Global Long. Strain",
               "Heart Rate", "P-wave Duration", "PR Interval", "QRS Duration", "QTc Interval"),
  Category = factor(c(rep("Atrial Remodeling", 3), 
                      rep("Ventricular Remodeling", 4), 
                      rep("Electrophysiology", 5)), 
                    levels = c("Atrial Remodeling", "Ventricular Remodeling", "Electrophysiology"))
)
df_plot <- df_plot %>% left_join(target_info, by = "Target")

# ——————————————————————————————————————————————
# 2. 构造 ComplexHeatmap 需要的三大宽矩阵 (Wide Matrices)
# ——————————————————————————————————————————————
# 我们需要 3 个矩阵，维度一模一样 (Rows = Targets, Cols = Metabolites)

# 确保行列顺序的排列逻辑
row_order_names <- target_info$Pub_Name
col_order_df <- df_plot %>% select(Metabolite_Short, Group) %>% distinct() %>% arrange(Group, Metabolite_Short)
col_order_names <- col_order_df$Metabolite_Short

# A. Beta 矩阵 (决定颜色)
mat_beta <- df_plot %>% select(Pub_Name, Metabolite_Short, Final_Beta) %>%
  pivot_wider(names_from = Metabolite_Short, values_from = Final_Beta) %>%
  tibble::column_to_rownames("Pub_Name")
mat_beta <- as.matrix(mat_beta[row_order_names, col_order_names])

# B. FDR 矩阵 (决定大小)
mat_fdr <- df_plot %>% select(Pub_Name, Metabolite_Short, FDR) %>%
  pivot_wider(names_from = Metabolite_Short, values_from = FDR) %>%
  tibble::column_to_rownames("Pub_Name")
mat_fdr <- as.matrix(mat_fdr[row_order_names, col_order_names])

# C. 星号矩阵 (文字标注)
mat_sig <- df_plot %>% select(Pub_Name, Metabolite_Short, Significance) %>%
  pivot_wider(names_from = Metabolite_Short, values_from = Significance) %>%
  tibble::column_to_rownames("Pub_Name")
mat_sig <- as.matrix(mat_sig[row_order_names, col_order_names])

# ——————————————————————————————————————————————
# 1. 颜色与基础图例准备
# ——————————————————————————————————————————————
target_cat_colors = c(
  "Atrial Remodeling" = "#E64B35FF",      
  "Ventricular Remodeling" = "#8431E0",   
  "Electrophysiology" = "#3C5488FF"       
)

max_b <- quantile(abs(mat_beta), 0.98, na.rm = TRUE)
col_fun_heat <- colorRamp2(c(-max_b, 0, max_b), c("navy", "white", "firebrick3"))

lgd_sig = Legend(
  labels = c("p = 0.05", "p = 1e-3", "p <= 1e-5"), 
  type = "points", pch = 16, 
  legend_gp = gpar(col = "grey60"),
  size = unit(c(3, 5.5, 8), "mm"), 
  title = "Significance (FDR)",
  title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_gap = unit(2.5, "mm") 
)

lgd_cardiac = Legend(
  labels = names(target_cat_colors),
  title = "Cardiac Target Class",
  legend_gp = gpar(fill = target_cat_colors),
  title_gp = gpar(fontsize = 10, fontface = "bold")
)
# ——————————————————————————————————————————————
# 1.5 同步权重图的颜色引擎
# ——————————————————————————————————————————————
# 1. 重新锁定横轴顺序：先按大类，再按代谢物缩写
col_order_df <- df_plot %>% 
  select(Metabolite_Short, Group_Short) %>% 
  distinct() %>% 
  arrange(Group_Short, Metabolite_Short)

col_order_names <- col_order_df$Metabolite_Short
groups_short_vec <- col_order_df$Group_Short 

# 2. 颜色定义：必须在这里手动生成和权重图一致的命名向量
unique_groups <- unique(groups_short_vec)
n_groups <- length(unique_groups)
# 🚨 确保颜色生成逻辑与权重图 100% 相同
group_colors <- setNames(colorRampPalette(brewer.pal(8, "Set2"))(n_groups), unique_groups)

# 打印检查一下，现在应该不止两类了
print(group_colors)

# 找到 2. 极致纯净的横纵色带部分，修改底部色带：
bottom_ann_v20 = HeatmapAnnotation(
  Category = groups_short_vec, # 🚨 使用缩写后的组名
  col = list(Category = group_colors), # 🚨 使用同步后的颜色向量
  show_annotation_name = FALSE,
  simple_anno_size = unit(4, "mm"),
  border = FALSE
)
# 重新制作三大矩阵 (这部分只需运行即可)
row_order_names <- target_info$Pub_Name
mat_beta <- df_plot %>% select(Pub_Name, Metabolite_Short, Final_Beta) %>%
  pivot_wider(names_from = Metabolite_Short, values_from = Final_Beta) %>%
  tibble::column_to_rownames("Pub_Name")
mat_beta <- as.matrix(mat_beta[row_order_names, col_order_names])

mat_fdr <- df_plot %>% select(Pub_Name, Metabolite_Short, FDR) %>%
  pivot_wider(names_from = Metabolite_Short, values_from = FDR) %>%
  tibble::column_to_rownames("Pub_Name")
mat_fdr <- as.matrix(mat_fdr[row_order_names, col_order_names])

# ——————————————————————————————————————————————
# 3. 最终组装 Heatmap (已修正 column_names_gp 冲突)
# ——————————————————————————————————————————————
ht_v20 = Heatmap(
  mat_beta, 
  name = "Std. Beta",
  rect_gp = gpar(type = "none"), 
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    fdr_val = as.numeric(mat_fdr[i, j])
    beta_val = as.numeric(mat_beta[i, j])
    if (!is.na(fdr_val) && fdr_val < 0.05) {
      w_pt = convertWidth(width, "pt", valueOnly = TRUE)
      h_pt = convertHeight(height, "pt", valueOnly = TRUE)
      max_safe_r = min(w_pt, h_pt) * 0.45 
      log_p = pmin(pmax(-log10(fdr_val), 1.3), 5)
      size_factor = 0.35 + (log_p - 1.3) / (5 - 1.3) * 0.65
      final_r = max_safe_r * size_factor
      grid.circle(x, y, r = unit(final_r, "pt"), 
                  gp = gpar(fill = col_fun_heat(beta_val), col = "white", lwd = 0.4))
    }
  },
  
  # 底部色带使用同步后的 Group_Short
  bottom_annotation = HeatmapAnnotation(
    Category = groups_short_vec,
    col = list(Category = group_colors),
    show_annotation_name = FALSE,
    simple_anno_size = unit(4, "mm"),
    border = FALSE
  ),
  left_annotation = left_ann_v20,
  
  # 🌟 分割逻辑同步
  column_split = factor(groups_short_vec, levels = unique(groups_short_vec)),
  column_gap = unit(0, "mm"),
  row_split = target_info$Category,
  row_gap = unit(0, "mm"), 
  
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, row_names_side = "left", 
  row_names_gp = gpar(fontsize = 11, fontface = "bold", col = target_cat_colors[target_info$Category]),
  
  show_column_names = TRUE, column_names_side = "bottom", column_names_rot = 90,
  # 🚨 使用同步后的 group_colors，颜色将与权重图完美匹配
  column_names_gp = gpar(fontsize = 10, fontface = "bold", 
                         col = group_colors[groups_short_vec]), 
  
  col = col_fun_heat
)


# 渲染保存
pdf("Mechanism_Bubble_Heatmap_Synced.pdf", width = 8, height = 7)
draw(ht_v20, annotation_legend_list = list(lgd_sig, lgd_cardiac), merge_legend = TRUE)
dev.off()

cat("✅ 报错已修复！同步了权重图大类颜色的气泡图已生成：Mechanism_Bubble_Heatmap_Synced.pdf\n")
library(dplyr)
library(readr)

# 专为大模型生成“摘要版”文字结果
df_for_llm <- df_results_final %>%
  # 1. 只看显著的，没必要让大模型读垃圾数据
  filter(FDR < 0.05) %>%
  # 2. 清理一下靶点名字（把恶心的 .x 去掉，方便大模型阅读）
  mutate(Target = gsub("\\.x$", "", Target)) %>%
  # 3. 按心脏靶点分组，并且把 Beta 绝对值最大的排在前面
  arrange(Target, desc(abs(Final_Beta))) %>%
  # 4. 只保留最核心的列
  select(Target, Metabolite, Final_Beta, FDR, Significance)

# 导出为 CSV 文件
write_csv(df_for_llm, "Results_for_LLM_Clean.csv")

cat("✅ 已成功导出大模型专用版结果！请用记事本或 Excel 打开并复制给它。\n")
