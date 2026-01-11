# ------------------------------------------------------------
# Script: 01_download_tcga.R
# Purpose:
#   Download publicly available TCGA-BRCA clinical data
#   and store it locally for downstream survival analysis.
#
# Notes:
#   - Uses TCGAbiolinks to access de-identified, open-access
#     clinical data from the NCI Genomic Data Commons (GDC).
#   - No protected or controlled-access data are required.
#   - This script is intended to be run once per environment
#     unless the raw data need to be refreshed.
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(readr)
})

# Create required directories if they do not already exist
dir.create("data_processed", showWarnings = FALSE, recursive = TRUE)
dir.create("results/logs", showWarnings = FALSE, recursive = TRUE)

# Set up a log file to capture messages from this script
log_file <- file("results/logs/download_log.txt", open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

cat("Starting download of TCGA-BRCA clinical data...\n")
cat("Source: NCI Genomic Data Commons (public clinical dataset)\n\n")

# Query and retrieve clinical-level data for the TCGA breast cancer cohort
clinical_raw <- GDCquery_clinic(
  project = "TCGA-BRCA",
  type = "clinical"
)

# Basic sanity checks and logging
cat("Download completed successfully.\n")
cat("Number of patient records:", nrow(clinical_raw), "\n")
cat("Number of clinical variables:", ncol(clinical_raw), "\n\n")

# Save raw clinical data in RDS format for reproducibility
saveRDS(
  clinical_raw,
  file = "data_processed/clinical_raw.rds"
)

cat("Raw clinical data saved to data_processed/clinical_raw.rds\n")
cat("This dataset will be used as input for endpoint derivation\n")
cat("and downstream survival analyses.\n")

# Close log sinks
sink(type = "output")
sink(type = "message")
close(log_file)