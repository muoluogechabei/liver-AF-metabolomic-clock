# AI-Driven Metabolomic Clock for Atrial Fibrillation

This repository contains the customized R scripts used for the statistical analyses and machine learning modelling in the manuscript: **"An AI-Driven Metabolomic Clock Reveals a Dual-Track Metabolic Axis Connecting Hepatic Fibro-Inflammation and Atrial Fibrillation"** (*Nature Communications*, 2026).

## Environment and Dependencies
All analyses were performed in **R version 4.5.1**. 
Key packages used in this pipeline include:
* **Survival Analysis:** `survival`, `cmprsk` (Fine-Gray models)
* **Machine Learning & Mediation:** `glmnet` (Elastic Net), `CMAverse` (Causal Mediation)
* **Clustering & Dimensionality Reduction:** `mclust` (Gaussian Mixture Models), `uwot` (UMAP)
* **Prediction & Calibration:** `riskRegression`
* **Data Imputation:** `mice`

## Repository Structure

The scripts are numbered sequentially to reflect the analytical workflow:

* `01_Data_Cleaning_Main.R`: Data extraction and preprocessing for the primary UK Biobank cohort.
* `02_Data_Cleaning_Subcohort.R`: Preprocessing for the liver-imaging and NMR metabolomics subcohorts.
* `03_Clinical_Scores_Calculation.R`: Computation of established clinical models (CHARGE-AF, ARIC, C2HEST).
* `04_Cox_Models_Main.R`: Multivariable Cox proportional hazards modelling for biochemical liver fibrosis indices.
* `05_Cox_Models_Subcohort.R`: Cox modelling for liver cT1 adjusted for PDFF.
* `06_Subgroup_Analysis.R`: Stratified analyses across prespecified baseline subgroups.
* `07_Metabolomics_Mediation.R`: Metabolome-wide screening and formal causal mediation analysis for 249 NMR features.
* `08_ML_Feature_Selection_and_MRS.R`: Elastic Net feature prioritization (500 bootstraps), Metabolomic Risk Score (MRS) derivation, and evaluation of incremental predictive value (AUC, cNRI, IDI).
* `09_GMM_Clustering_and_Survival.R`: Unsupervised phenotype clustering, internal stability validation (ARI), and UMAP visualization.
* `10_Target_Organ_Validation.R`: Phenomapping of core metabolites to CMR and ECG metrics.
* `/Sensitivity_Analyses/`: Folder containing six distinct scripts for robustness checks (e.g., competing risk models, IPTW, alternative thresholds).

## Data Availability
Individual-level data are available through the UK Biobank (Application No. 366548). To protect participant privacy, raw data are not hosted in this repository. Access requires standard UK Biobank approval (http://www.ukbiobank.ac.uk/).
