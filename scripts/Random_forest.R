# ============================================================
# Title:   Microbiome + Clinical Random Forest Pipeline 
# Author:  Lou Langhammer
# Purpose: Predict a continuous clinical outcome (MRI_Score)
#          from CLR-transformed microbial taxa + clinical
#          covariates. Uses repeated train/test splits with
#          Boruta feature selection strictly inside each fold
#          to avoid data leakage. Includes a permutation test
#          comparing full vs. covariate-only (reduced) models.
# ============================================================


# -- 1. Package management -----------------------------------

library(caret)
library(Boruta)
library(randomForest)
library(dplyr)
library(tidyr)

# -- 2. User settings ----------------------------------------

infile       <- "data/RF_input_totalscore_CLR.txt"  # path to input table
outdir       <- "results"                            # output directory (created if absent)
outcome_var  <- "MRI_Score"                          # column name of the response variable

# Clinical covariate column names (exist in infile)
cont_covars  <- c("Age", "weight", "DD")   # continuous covariates
bin_covars   <- c("hypertension")          # binary covariates (0/1) -> coerced to factor
fac_covars   <- c()                        # factor/categorical covariates (coerced to factor)

# Regex to identify taxon columns in the input table.
# NOTE: this is GTDB-style; change for SILVA / NCBI taxonomy.
taxa_pattern <- "^d__(Bacteria|Archaea)"

# Split & resampling settings
train_prop   <- 0.60   # fraction of samples used for training in each split
n_splits     <- 5      # number of train/test splits  (publication: >= 50)
n_boruta     <- 20     # Boruta runs per split        (publication: >= 50)
n_iter       <- 20     # RF iterations per split      (publication: >= 50)
base_seed    <- 42     # base seed; split s uses base_seed + s

# Boruta settings
boruta_min_freq <- 0.25   # minimum WITHIN-SPLIT selection frequency to retain a feature
boruta_maxRuns  <- 300    # maximum Boruta iterations (higher = more stable, slower)

# Random Forest settings
subsample_prop <- 0.80   # fraction of training data used in each RF iteration
rf_ntree       <- 500    # number of trees per forest
cv_folds       <- 5      # k for cross-validation within caret

# Permutation / sign-flip test settings
do_scramble  <- TRUE     # also run a scrambled-label negative-control RF
perm_B       <- 10000    # number of permutations for the sign-flip test
perm_stat    <- "median" # summary statistic for permutation test: "median" or "mean"

# Filtering: if TRUE, drops any taxon with non-finite values;
#            if FALSE, replaces non-finite values with 0 before modelling.
strict_taxa_finite <- FALSE


# -- 3. Derived output paths ---------------------------------

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

out_metrics  <- file.path(outdir, "rf_validation_metrics.tsv")
out_boruta   <- file.path(outdir, "boruta_selection_frequencies.tsv")
out_imp      <- file.path(outdir, "rf_feature_importance.tsv")
out_deltas   <- file.path(outdir, "full_vs_reduced_per_split_deltas.tsv")
out_test     <- file.path(outdir, "full_vs_reduced_signflip_test.tsv")
out_scramble <- file.path(outdir, "rf_scrambled_negative_control.tsv")
out_session  <- file.path(outdir, "session_info.txt")


# -- 4. I/O and metric helpers -------------------------------

#' Read a tab-separated file, stopping with a clear message if it is missing.
read_table_tsv <- function(path) {
  if (!file.exists(path)) stop("Input file not found: ", path)
  read.table(
    path, header = TRUE, sep = "\t",
    stringsAsFactors = FALSE, check.names = FALSE,
    comment.char = "", quote = ""
  )
}

#' Write a data frame as a tab-separated file (no row names, no quoting).
write_tsv_local <- function(x, path) {
  write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE)
}

#' If the first column looks like a unique sample-ID column, promote it to
#' rownames and drop it from the data frame; otherwise return df unchanged.
use_first_col_as_ids <- function(df) {
  first <- df[[1]]
  is_id_col <- !is.numeric(first) &&
    length(unique(first)) == nrow(df) &&
    !anyDuplicated(first)

  if (is_id_col) {
    rownames(df) <- trimws(as.character(first))
    df[[1]] <- NULL
  }
  df
}

#' Standard coefficient of determination: 1 - SS_res / SS_tot.
#' Unlike cor(pred, obs)^2, this penalises systematic bias / scale errors
#' and can be negative if the model performs worse than the mean.
r2_standard <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs)
  if (sum(ok) < 2) return(NA_real_)
  ss_res <- sum((obs[ok] - pred[ok])^2)
  ss_tot <- sum((obs[ok] - mean(obs[ok]))^2)
  if (ss_tot == 0) return(NA_real_)
  1 - ss_res / ss_tot
}

#' Squared Pearson correlation (matches the original script's R^2).
r2_cor <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs)
  if (sum(ok) < 2) return(NA_real_)
  suppressWarnings(cor(pred[ok], obs[ok]))^2
}


# -- 5. Column resolution helper -----------------------------

#' Find a column in `df` by name, with case-insensitive fallback.
#' Stops with an informative message if the column is absent or ambiguous.
resolve_col <- function(df, target) {
  if (target %in% names(df)) return(target)

  hits <- names(df)[tolower(names(df)) == tolower(target)]

  if (length(hits) == 1) {
    message("Note: column '", target, "' matched case-insensitively as '", hits, "'.")
    return(hits)
  }
  if (length(hits) > 1) {
    stop(
      "Ambiguous column name for '", target, "'. ",
      "Candidates: ", paste(hits, collapse = ", ")
    )
  }

  stop("Column not found: '", target, "'. Check spelling and the input file headers.")
}


# -- 6. Data loading & validation ----------------------------

message("Loading data from: ", infile)
raw <- read_table_tsv(infile)
raw <- use_first_col_as_ids(raw)

# Resolve and validate all user-specified columns
outcome_col  <- resolve_col(raw, outcome_var)
cont_cols    <- sapply(cont_covars, resolve_col, df = raw)
bin_cols     <- if (length(bin_covars) > 0) sapply(bin_covars, resolve_col, df = raw) else character(0)
fac_cols     <- if (length(fac_covars) > 0) sapply(fac_covars, resolve_col, df = raw) else character(0)

# Coerce binary and factor covariates to factor so RF / Boruta treat them
# as categorical rather than continuous.
for (col in c(bin_cols, fac_cols)) {
  raw[[col]] <- as.factor(raw[[col]])
}

# Identify taxon columns
taxa_cols <- grep(taxa_pattern, names(raw), value = TRUE)
if (length(taxa_cols) == 0) {
  stop("No taxon columns matched the pattern: ", taxa_pattern)
}
message(sprintf("Found %d taxon columns matching pattern '%s'.", length(taxa_cols), taxa_pattern))

# Handle non-finite taxon values
n_nonfinite <- sum(!is.finite(as.matrix(raw[, taxa_cols])))
if (n_nonfinite > 0) {
  if (strict_taxa_finite) {
    bad_taxa <- taxa_cols[apply(raw[, taxa_cols], 2, function(x) any(!is.finite(x)))]
    message(sprintf(
      "strict_taxa_finite = TRUE: dropping %d taxa with non-finite values.", length(bad_taxa)
    ))
    taxa_cols <- setdiff(taxa_cols, bad_taxa)
  } else {
    message(sprintf(
      "strict_taxa_finite = FALSE: replacing %d non-finite taxon values with 0.", n_nonfinite
    ))
    raw[, taxa_cols][!is.finite(as.matrix(raw[, taxa_cols]))] <- 0
  }
}

all_feature_cols <- c(taxa_cols, cont_cols, bin_cols, fac_cols)
model_data       <- raw[, c(outcome_col, all_feature_cols), drop = FALSE]
model_data       <- model_data[complete.cases(model_data), ]
message(sprintf("Samples after removing rows with missing values: %d", nrow(model_data)))


# -- 7. Main loop: repeated train/test splits ----------------
# Boruta feature selection runs entirely within the training fold to prevent
# information from the test set influencing which features are chosen.
# PATCH: the threshold uses a split-local Boruta counter; a separate global
# counter accumulates across splits purely for the final report.

message(sprintf("\nStarting %d train/test splits ...", n_splits))

metrics_list <- vector("list", n_splits)
boruta_freq_global <- setNames(numeric(length(all_feature_cols)), all_feature_cols)
imp_list     <- vector("list", n_splits)
delta_list   <- vector("list", n_splits)
scramble_list <- vector("list", n_splits)

covar_cols <- c(cont_cols, bin_cols, fac_cols)

for (s in seq_len(n_splits)) {

  message(sprintf("  Split %d / %d", s, n_splits))
  set.seed(base_seed + s)

  # 7a. Partition -----------------------------------------------------------
  train_idx <- createDataPartition(model_data[[outcome_col]],
                                    p = train_prop, list = FALSE)
  train_df  <- model_data[ train_idx, ]
  test_df   <- model_data[-train_idx, ]

  # 7b. Boruta on training data only ----------------------------------------
  # PATCH: use a split-local counter for thresholding.
  boruta_freq_local <- setNames(numeric(length(all_feature_cols)), all_feature_cols)

  for (b in seq_len(n_boruta)) {
    sub_idx <- sample(nrow(train_df), size = floor(subsample_prop * nrow(train_df)))
    sub_df  <- train_df[sub_idx, ]

    bor <- tryCatch(
      Boruta(
        x       = sub_df[, all_feature_cols],
        y       = sub_df[[outcome_col]],
        maxRuns = boruta_maxRuns,
        doTrace = 0
      ),
      error = function(e) {
        message("    Boruta run ", b, " failed: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(bor)) {
      confirmed <- names(bor$finalDecision[bor$finalDecision == "Confirmed"])
      boruta_freq_local[confirmed]  <- boruta_freq_local[confirmed]  + 1
      boruta_freq_global[confirmed] <- boruta_freq_global[confirmed] + 1
    }
  }

  # Retain only features selected in >= boruta_min_freq of THIS split's runs
  local_rate        <- boruta_freq_local / n_boruta
  selected_features <- names(local_rate[local_rate >= boruta_min_freq])

  if (length(selected_features) == 0) {
    message("    No features met the Boruta frequency threshold; skipping split.")
    next
  }

  # 7c. Paired RF iterations: full and reduced share subsample --------------
  # PATCH: same sub_idx for full and reduced => paired delta-R^2.
  split_r2_full     <- rep(NA_real_, n_iter)   # squared correlation (legacy)
  split_r2_full_std <- rep(NA_real_, n_iter)   # standard R^2
  split_r2_red      <- rep(NA_real_, n_iter)
  split_r2_red_std  <- rep(NA_real_, n_iter)
  split_r2_scram    <- rep(NA_real_, n_iter)
  split_r2_scram_std <- rep(NA_real_, n_iter)
  split_imp         <- vector("list", n_iter)

  ctrl <- trainControl(method = "cv", number = cv_folds)

  full_features <- union(selected_features, covar_cols)

  for (i in seq_len(n_iter)) {
    sub_idx <- sample(nrow(train_df), size = floor(subsample_prop * nrow(train_df)))
    sub_df  <- train_df[sub_idx, ]

    # --- Full model: selected taxa + clinical covariates ---
    fit_full <- tryCatch(
      train(
        x         = sub_df[, full_features, drop = FALSE],
        y         = sub_df[[outcome_col]],
        method    = "rf",
        ntree     = rf_ntree,
        trControl = ctrl
      ),
      error = function(e) NULL
    )

    if (!is.null(fit_full)) {
      preds_full            <- predict(fit_full, newdata = test_df[, full_features, drop = FALSE])
      split_r2_full[i]      <- r2_cor(preds_full,      test_df[[outcome_col]])
      split_r2_full_std[i]  <- r2_standard(preds_full, test_df[[outcome_col]])
      split_imp[[i]]        <- varImp(fit_full)$importance
    }

    # --- Reduced model: clinical covariates only, SAME sub_idx ---
    if (length(covar_cols) > 0) {
      fit_red <- tryCatch(
        train(
          x         = sub_df[, covar_cols, drop = FALSE],
          y         = sub_df[[outcome_col]],
          method    = "rf",
          ntree     = rf_ntree,
          trControl = ctrl
        ),
        error = function(e) NULL
      )

      if (!is.null(fit_red)) {
        preds_red             <- predict(fit_red, newdata = test_df[, covar_cols, drop = FALSE])
        split_r2_red[i]       <- r2_cor(preds_red,      test_df[[outcome_col]])
        split_r2_red_std[i]   <- r2_standard(preds_red, test_df[[outcome_col]])
      }
    }

    # --- Scrambled-label negative control, SAME sub_idx ---
    if (do_scramble) {
      sub_df_scram <- sub_df
      sub_df_scram[[outcome_col]] <- sample(sub_df[[outcome_col]])

      fit_scram <- tryCatch(
        train(
          x         = sub_df_scram[, full_features, drop = FALSE],
          y         = sub_df_scram[[outcome_col]],
          method    = "rf",
          ntree     = rf_ntree,
          trControl = ctrl
        ),
        error = function(e) NULL
      )

      if (!is.null(fit_scram)) {
        preds_scram            <- predict(fit_scram, newdata = test_df[, full_features, drop = FALSE])
        split_r2_scram[i]      <- r2_cor(preds_scram,      test_df[[outcome_col]])
        split_r2_scram_std[i]  <- r2_standard(preds_scram, test_df[[outcome_col]])
      }
    }
  }

  # 7d. Per-split summaries -------------------------------------------------
  metrics_list[[s]] <- data.frame(
    split                 = s,
    n_selected            = length(selected_features),
    mean_R2cor_full       = mean(split_r2_full,      na.rm = TRUE),
    median_R2cor_full     = median(split_r2_full,    na.rm = TRUE),
    mean_R2std_full       = mean(split_r2_full_std,  na.rm = TRUE),
    median_R2std_full     = median(split_r2_full_std, na.rm = TRUE),
    mean_R2cor_reduced    = mean(split_r2_red,       na.rm = TRUE),
    median_R2cor_reduced  = median(split_r2_red,     na.rm = TRUE),
    mean_R2std_reduced    = mean(split_r2_red_std,   na.rm = TRUE),
    median_R2std_reduced  = median(split_r2_red_std, na.rm = TRUE),
    mean_R2std_scrambled  = if (do_scramble) mean(split_r2_scram_std,   na.rm = TRUE) else NA_real_,
    median_R2std_scrambled = if (do_scramble) median(split_r2_scram_std, na.rm = TRUE) else NA_real_
  )

  # Paired deltas (standard R^2)
  delta_list[[s]] <- data.frame(
    split         = s,
    iter          = seq_len(n_iter),
    R2_full_std   = split_r2_full_std,
    R2_red_std    = split_r2_red_std,
    delta_R2_std  = split_r2_full_std - split_r2_red_std
  )

  if (do_scramble) {
    scramble_list[[s]] <- data.frame(
      split          = s,
      iter           = seq_len(n_iter),
      R2_scram_cor   = split_r2_scram,
      R2_scram_std   = split_r2_scram_std
    )
  }

  # Pool feature importances across iterations
  if (length(split_imp) > 0) {
    imp_combined <- do.call(rbind, lapply(seq_along(split_imp), function(i) {
      if (is.null(split_imp[[i]])) return(NULL)
      df <- as.data.frame(split_imp[[i]])
      df$feature <- rownames(df)
      df$split   <- s
      df$iter    <- i
      df
    }))
    imp_list[[s]] <- imp_combined
  }
}


# -- 8. Write outputs ---------------------------------------

message("\nWriting outputs ...")

# Validation metrics
metrics_df <- do.call(rbind, Filter(Negate(is.null), metrics_list))
write_tsv_local(metrics_df, out_metrics)

# Boruta selection frequencies (global, normalised to [0, 1])
boruta_df <- data.frame(
  feature   = names(boruta_freq_global),
  frequency = boruta_freq_global / (n_splits * n_boruta)
)
boruta_df <- boruta_df[order(-boruta_df$frequency), ]
write_tsv_local(boruta_df, out_boruta)

# Feature importances
imp_df <- do.call(rbind, Filter(Negate(is.null), imp_list))
if (!is.null(imp_df)) write_tsv_local(imp_df, out_imp)

# Per-split, per-iter paired deltas
delta_df <- do.call(rbind, Filter(Negate(is.null), delta_list))
write_tsv_local(delta_df, out_deltas)

# Scrambled negative control
if (do_scramble) {
  scramble_df <- do.call(rbind, Filter(Negate(is.null), scramble_list))
  if (!is.null(scramble_df) && nrow(scramble_df) > 0) {
    write_tsv_local(scramble_df, out_scramble)
  }
}


# -- 9. Sign-flip permutation test on paired deltas ----------
# Uses standard R^2 because that's the metric whose differences are
# directly interpretable as variance-explained gains.

if (nrow(delta_df) > 0) {
  delta_vec <- delta_df$delta_R2_std
  delta_vec <- delta_vec[is.finite(delta_vec)]

  message(sprintf(
    "Running sign-flip permutation test (B = %d, stat = %s, n_pairs = %d) ...",
    perm_B, perm_stat, length(delta_vec)
  ))

  stat_fn  <- if (perm_stat == "median") median else mean
  obs_stat <- stat_fn(delta_vec)

  perm_stats <- replicate(perm_B, {
    flipped <- delta_vec * sample(c(-1, 1), length(delta_vec), replace = TRUE)
    stat_fn(flipped)
  })

  # Two-sided p-value with +1 numerator/denominator for a valid Monte-Carlo bound
  p_val <- (sum(abs(perm_stats) >= abs(obs_stat)) + 1) / (perm_B + 1)

  test_df <- data.frame(
    statistic      = perm_stat,
    metric         = "standard_R2",
    observed       = obs_stat,
    p_value        = p_val,
    n_pairs        = length(delta_vec),
    n_permutations = perm_B
  )
  write_tsv_local(test_df, out_test)

  message(sprintf(
    "  Observed delta-%s standard R^2: %.4f  |  p-value: %.4f",
    perm_stat, obs_stat, p_val
  ))
}


# -- 10. Session info for reproducibility --------------------

message("\nSession information (for reproducibility):")
print(sessionInfo())

# Also write sessionInfo to a file in outdir
con <- file(out_session, open = "wt")
writeLines(capture.output(print(sessionInfo())), con)
close(con)

message("\nPipeline complete. Results written to: ", normalizePath(outdir))
