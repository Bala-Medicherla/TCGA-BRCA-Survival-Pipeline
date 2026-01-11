message("Running TCGA-BRCA survival pipeline...")

source("scripts/01_download_tcga.R")
source("scripts/02_clean_clinical.R")
source("scripts/03_survival_analysis.R")
source("scripts/04_model_validation.R")

message("Pipeline completed. Check the results/ folder.")