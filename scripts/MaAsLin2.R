# ============================================================
# Title:   MaAsLin2 Microbiome Association Model
# Author:  Lou Langhammer
# Purpose: Test associations between gut microbial taxa and
#          MRI score, adjusting for age, disease duration,
#          weight, and hypertension. Uses MaAsLin2 defaults
#          for normalisation and transformation with BH
#          multiple-testing correction and a 33% minimum
#          prevalence filter to reduce low-prevalence noise.
# ============================================================


# -- 1. Package management -----------------------------------

library(Maaslin2)

# -- 2. User settings ----------------------------------------

# Project working directory. Adjust before running.
workdir <- "/Users/.../study_1/08-MaasLin2"

# Input files
infile_data     <- "Maaslin_data_ASV_input.txt"       # feature table (samples x taxa or taxa x samples)
infile_metadata <- "study_1_metadata_4_including_sample_removed.txt"  # sample metadata

# Metadata column names. Must match the metadata file header exactly.
outcome_var  <- "MRI_Score"    # continuous outcome
cov_age      <- "Age"          # continuous covariate
cov_dd       <- "DD"           # disease duration (continuous)
cov_weight   <- "Weight"       # body weight (continuous)
cov_htn      <- "Hypertension" # binary covariate; must be coded 0/1

# MaAsLin2 model settings.
# Normalisation and transformation are left at MaAsLin2 defaults (TSS + LOG).
# standardize = TRUE z-scores continuous metadata columns before fitting,
# which puts regression coefficients on a comparable scale across covariates.
correction    <- "BH"    # multiple-testing correction method
max_sig       <- 0.05    # q-value threshold for significance
min_prev      <- 0.33    # minimum sample prevalence to retain a feature (33%)

# Output directory. Name encodes the key model choices for traceability.
outdir <- "251108_ASV_MRI_Score_genus_cov_Age_DD_Weight_Hypertension_MinPrev33pct"


# -- 3. Derived settings -------------------------------------

setwd(workdir)

# Collect all model variable names for downstream validation steps.
model_vars <- c(outcome_var, cov_age, cov_dd, cov_weight, cov_htn)


# -- 4. I/O helper -------------------------------------------

#' Read a tab-separated file, stopping with a clear message if it is missing.
read_table_tsv <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  read.table(
    path, header = TRUE, sep = "\t",
    row.names = 1, check.names = FALSE,
    stringsAsFactors = FALSE
  )
}


# -- 5. Data loading -----------------------------------------

message("Loading feature table: ",     infile_data)
message("Loading sample metadata: ",   infile_metadata)

df_data <- read_table_tsv(infile_data)
df_meta <- read_table_tsv(infile_metadata)


# -- 6. Sample ID cleaning -----------------------------------

# Strip invisible whitespace from row and column names;
# mismatches here are a common source of hard-to-debug errors.
rownames(df_meta) <- trimws(rownames(df_meta))
rownames(df_data) <- trimws(rownames(df_data))
colnames(df_data) <- trimws(colnames(df_data))


# -- 7. Feature-table orientation ----------------------------

# MaAsLin2 expects samples in rows and features in columns.
# If sample IDs are currently in the columns, transpose the table.
if (all(colnames(df_data) %in% rownames(df_meta))) {
  message("Samples detected in columns of feature table -> transposing.")
  df_data <- t(df_data)
}


# -- 8. Metadata column validation ---------------------------

missing_cols <- setdiff(model_vars, colnames(df_meta))
if (length(missing_cols) > 0) {
  stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}


# -- 9. Coerce model variables to numeric --------------------

#' Convert a metadata column to numeric, stopping if non-numeric
#' values are found (excluding NA and empty strings).
coerce_numeric <- function(x, col_name) {
  if (is.factor(x))    x <- as.character(x)
  if (is.character(x)) x <- trimws(x)
  xn  <- suppressWarnings(as.numeric(x))
  bad <- is.na(xn) & !is.na(x) & x != ""
  if (any(bad)) stop("Non-numeric values in column: ", col_name)
  xn
}

df_meta[[outcome_var]] <- coerce_numeric(df_meta[[outcome_var]], outcome_var)
df_meta[[cov_age]]     <- coerce_numeric(df_meta[[cov_age]],     cov_age)
df_meta[[cov_dd]]      <- coerce_numeric(df_meta[[cov_dd]],      cov_dd)
df_meta[[cov_weight]]  <- coerce_numeric(df_meta[[cov_weight]],  cov_weight)
df_meta[[cov_htn]]     <- coerce_numeric(df_meta[[cov_htn]],     cov_htn)

# Confirm hypertension is binary; any value outside {0, 1} is an error.
htn_levels <- sort(unique(df_meta[[cov_htn]]))
if (!all(htn_levels %in% c(0, 1))) {
  stop(cov_htn, " must be coded 0/1 only. Found: ", paste(htn_levels, collapse = ", "))
}


# -- 10. Drop samples with missing model variables -----------

complete_idx  <- complete.cases(df_meta[, model_vars])
df_meta       <- df_meta[complete_idx, , drop = FALSE]
message(sprintf("Samples retained after removing incomplete cases: %d", nrow(df_meta)))


# -- 11. Match samples between metadata and feature table ----

common_ids <- intersect(rownames(df_meta), rownames(df_data))

if (length(common_ids) < 5) {
  stop(
    "Too few overlapping sample IDs (", length(common_ids), "). ",
    "Check that row names in the metadata and feature table use the same format.\n",
    "  Example metadata IDs: ", paste(head(rownames(df_meta)), collapse = ", "), "\n",
    "  Example feature IDs:  ", paste(head(rownames(df_data)), collapse = ", ")
  )
}

df_meta <- df_meta[common_ids, , drop = FALSE]
df_data <- df_data[common_ids, , drop = FALSE]
message("Matched samples: ", length(common_ids))


# -- 12. Feature-table sanity checks -------------------------

feat_range <- range(as.matrix(df_data), na.rm = TRUE)
message(sprintf("Feature table range: [%.4f, %.4f]", feat_range[1], feat_range[2]))

# Negative values are incompatible with MaAsLin2's default TSS + LOG pipeline.
# If the input is CLR-transformed, set normalization = "NONE", transform = "NONE".
if (feat_range[1] < 0) {
  stop(
    "Feature table contains negative values, which are incompatible with the ",
    "default TSS + LOG normalisation. If the table is already CLR-transformed, ",
    "rerun with normalization = 'NONE' and transform = 'NONE'."
  )
}


# -- 13. Run MaAsLin2 ----------------------------------------

message("Running MaAsLin2 ...")

fit <- Maaslin2(
  input_data     = df_data,
  input_metadata = df_meta,
  output         = outdir,
  fixed_effects  = model_vars,      # outcome + all covariates enter as fixed effects
  correction     = correction,      # BH false discovery rate
  max_significance = max_sig,       # q-value cutoff for reported associations
  min_prevalence = min_prev,        # drop features present in < 33% of samples
  standardize    = TRUE             # z-score continuous predictors before fitting
)

message("Done. Results written to: ", normalizePath(outdir))
