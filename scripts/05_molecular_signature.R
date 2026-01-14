# ------------------------------------------------------------
# Script: 05_molecular_signature.R
# Purpose:
#   Optional molecular extension for TCGA-BRCA:
#     1) Download gene expression (STAR - Counts)
#     2) Compute a simple proliferation signature
#     3) Merge with clinical survival data
#     4) Fit a Cox model (age + stage + proliferation score)
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tibble)
  library(survival)
})

dir.create("data_processed", showWarnings = FALSE, recursive = TRUE)
results_tables_dir <- file.path("results", "tables")
legacy_tables_dir <- file.path("results", "Tables")
results_logs_dir <- file.path("results", "logs")

if (dir.exists(legacy_tables_dir) && !dir.exists(results_tables_dir)) {
  results_tables_dir <- legacy_tables_dir
}

dir.create("results/Tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/logs", showWarnings = FALSE, recursive = TRUE)

log_file <- file("results/logs/molecular_log.txt", open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")
on.exit({
  sink(type = "output")
  sink(type = "message")
  close(log_file)
}, add = TRUE)

cat("Starting optional molecular extension (TCGA-BRCA gene expression)...\n")
cat("This step may take a while depending on download speed.\n\n")

# ------------------------------------------------------------
# 1) Load clinical survival dataset
# ------------------------------------------------------------
clinical <- readRDS("data_processed/clinical_surv.rds")
cat("Loaded clinical survival dataset.\n")
cat("Clinical dataset rows:", nrow(clinical), "\n\n")

pick_first_existing <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(nm)
  return(NA_character_)
}

patient_id_col <- pick_first_existing(
  clinical,
  c("bcr_patient_barcode", "submitter_id", "case_id", "patient_id")
)

if (is.na(patient_id_col)) {
  stop("Patient ID column not found in clinical dataset.")
}

clinical <- clinical %>%
  mutate(patient_id = str_sub(as.character(.data[[patient_id_col]]), 1, 12))

# ------------------------------------------------------------
# 2) Download gene expression (STAR - Counts)
# ------------------------------------------------------------
query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)


manifest <- tryCatch({
  getResults(query)
}, error = function(e) {
  query$results[[1]]
})

if (nrow(manifest) == 0) {
  cat("No gene expression files returned for the query. Skipping molecular extension.\n")
  cat("Molecular extension completed with no data.\n")
} else {
  download_attempt <- tryCatch({
    GDCdownload(query)
    TRUE
  }, error = function(e) {
    cat("GDCdownload failed. Skipping molecular extension.\n")
    cat("Error:", conditionMessage(e), "\n")
    FALSE
  })
}

# ------------------------------------------------------------
# 3) Build a simple proliferation signature (Memory-Efficient)
# ------------------------------------------------------------
prolif_genes <- c("MKI67", "CCNB1", "CDK1", "TOP2A", "BUB1")

cat("Retrieving file manifest from query...\n")
# Get the manifest which links file names to sample barcodes
manifest <- tryCatch({
  getResults(query)
}, error = function(e) {
  # Fallback for older TCGAbiolinks versions
  query$results[[1]]
})

# Filter for the files we just downloaded
# Using the project directory structure created by GDCdownload
source_dir <- file.path("GDCdata", "TCGA-BRCA", "Transcriptome_Profiling", "Gene_Expression_Quantification")

# We need to map file names to patient IDs. 
# The manifest usually has 'file_name' and 'cases' (which is the sample barcode)
files_to_process <- manifest %>%
  select(file_name, cases) %>%
  distinct()

cat("Found", nrow(files_to_process), "files to process.\n")

# Function to read a single file and extract only the proliferation genes
read_prolif_genes <- function(f_name, case_id, required_genes) {
  # Build full path (GDCdownload organizes by UUID folders usually, but let's try recursive search or direct path if possible)
  # Actually, GDCdownload creates: GDCdata/Project/Category/Type/UUID/filename
  # But we can verify with list.files
  
  full_path <- list.files(
    path = "GDCdata", 
    pattern = f_name, 
    recursive = TRUE, 
    full.names = TRUE
  )
  
  if (length(full_path) == 0) {
    warning(paste("File not found:", f_name))
    return(NULL)
  }
  
  # Read only necessary columns
  # STAR-Counts file usually has: gene_id, gene_name, gene_type, unstranded, ..
  # We skip lines (often 2-6 header lines in some formats, but STAR-Counts TSV usually has a header row or comments)
  # using read_tsv with comment argument is safest.
  
  d <- suppressMessages(read_tsv(
    full_path[1], 
    col_types = cols_only(
      gene_name = col_character(), 
      unstranded = col_double()
    ),
    comment = "#",
    show_col_types = FALSE
  ))
  
  # Filter for our genes
  d_sub <- d %>%
    filter(gene_name %in% required_genes) %>%
    mutate(sample_barcode = case_id)
  
  return(d_sub)
}

cat("Reading files iteratively (low memory mode)...\n")

# Use lapply to read files and bind them (much lighter than loading a giant matrix)
# Requires 'purrr' or 'dplyr' bind_rows
expr_list <- vector("list", nrow(files_to_process))

# Progress bar could be nice, but simple print is safer for logs
for (i in seq_len(nrow(files_to_process))) {
  if (i %% 100 == 0) cat("Processed", i, "/", nrow(files_to_process), "samples...\n")
  
  expr_list[[i]] <- read_prolif_genes(
    files_to_process$file_name[i], 
    files_to_process$cases[i], 
    prolif_genes
  )
}

expr_df <- bind_rows(expr_list)
cat("Finished reading expression data.\n")

if (nrow(expr_df) == 0) {
  stop("No proliferation gene data extracted. Check gene names or file formats.")
}

# Calculate proliferation score per sample -> then average per patient
# Score = mean of Z-scores of the 5 genes
# First, pivot to wide format to calculate Z-scores across samples
expr_wide <- expr_df %>%
  select(sample_barcode, gene_name, unstranded) %>%
  tidyr::pivot_wider(names_from = gene_name, values_from = unstranded) %>%
  column_to_rownames("sample_barcode")

# Calculate Z-scores
# scale() operates on columns (genes)
expr_z <- scale(expr_wide)

# Calculate mean Z-score (Proliferation Score) for each sample
prolif_score_vec <- rowMeans(expr_z, na.rm = TRUE)

expr_scores <- tibble(
  sample_barcode = names(prolif_score_vec),
  patient_id = str_sub(sample_barcode, 1, 12),
  proliferation_score = as.numeric(prolif_score_vec)
) %>%
  group_by(patient_id) %>%
  summarize(proliferation_score = mean(proliferation_score, na.rm = TRUE), .groups = "drop")

# ------------------------------------------------------------
# 4) Merge with clinical data and fit Cox model
# ------------------------------------------------------------
merged <- clinical %>%
  inner_join(expr_scores, by = "patient_id") %>%
  filter(!is.na(stage_group), !is.na(age_years), !is.na(proliferation_score))

cat("Merged clinical + expression data.\n")
cat("Merged dataset rows:", nrow(merged), "\n\n")

if (nrow(merged) == 0 || sum(merged$os_event, na.rm = TRUE) == 0) {
  stop("Insufficient data/events after merging molecular scores.")
}


cat("Stage distribution:\n")
print(table(merged$stage_group, merged$os_event, dnn = c("Stage", "Event")))

# Collapse stages into binary groups to ensure model stability (fix convergence warning)
merged <- merged %>%
  mutate(stage_binary = case_when(
    stage_group %in% c("I", "II") ~ "Early",
    stage_group %in% c("III", "IV") ~ "Late",
    TRUE ~ NA_character_
  ))

cat("Binary stage distribution:\n")
print(table(merged$stage_binary, merged$os_event, dnn = c("Stage_Binary", "Event")))

cox_fit <- coxph(
  Surv(os_time_days, os_event) ~ age_years + stage_binary + proliferation_score,
  data = merged
)

cox_sum <- summary(cox_fit)
cox_table <- as.data.frame(cox_sum$coefficients) %>%
  rownames_to_column("term")

write_csv(cox_table, "results/Tables/cox_molecular_summary.csv")
cat("Saved molecular Cox summary to results/Tables/cox_molecular_summary.csv\n")
cat("Tip: compare this model to the clinical-only Cox model for context.\n")

cat("Molecular extension completed successfully.\n")