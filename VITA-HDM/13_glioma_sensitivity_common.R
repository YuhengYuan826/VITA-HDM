# ==============================================================================
# 13_glioma_sensitivity_common.R
# Common utilities for target-anchored glioma VITA-HDM sensitivity analyses.
# Aligned with the repository case-study workflow.
# Main input:
#   results/case_study/discovery/
#     glioma_final_selection_objects.rds
# Key design:
#   - Target-anchored reporting rule:
#       source evidence annotates transferability but does not veto target-only
#       stable signals.
#   - Single final source:
#       TCGA_LGG
#   - Sensitivity analyses track:
#       target-anchored reported set
#       transfer-supported set
#       source-supported subset
#       high-stability core
#       target-only stable genes
#       target-only source-discordant genes
#       target-enhanced transfer-retained genes
# Outputs:
#   results/case_study/sensitivity/
# No sink() logging is used.
# ==============================================================================
if (!exists("PROJECT_ROOT", inherits = FALSE) ||
    !exists("PATHS", inherits = FALSE)) {
  stop(
    "13_glioma_sensitivity_common.R must be loaded by a numbered sensitivity script.",
    call. = FALSE
  )
}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(survival)
  library(ggplot2)
})
source(file.path(PROJECT_ROOT, "08_glioma_vita_functions.R"))
# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
PROJECT_DIR <- PROJECT_ROOT
RESULTS_ROOT <- PATHS$case_results
DATA_DIR <- PATHS$case_data
DISC_DIR <- PATHS$discovery
CV_DIR <- PATHS$fixed_signature_cv
EFF_DIR <- PATHS$conditional_effects
FINAL_SELECTION_RDS <- file.path(
  DISC_DIR,
  "glioma_final_selection_objects.rds"
)
if (!file.exists(FINAL_SELECTION_RDS)) {
  stop(
    "Cannot find final selection object: ", FINAL_SELECTION_RDS,
    "\nRun 09_run_glioma_discovery.R before sensitivity analyses."
  )
}
SENS_DIR <- ensure_directory(PATHS$sensitivity)
SENS_FIG_DIR <- ensure_directory(PATHS$sensitivity_figures)
SENS_TAB_DIR <- ensure_directory(PATHS$sensitivity_tables)
SENS_PRIVATE_DIR <- ensure_directory(PATHS$sensitivity_private)
SENS_FINAL_DIR <- ensure_directory(PATHS$sensitivity_final)
# ------------------------------------------------------------------------------
# Load final object and parameters
# ------------------------------------------------------------------------------
obj <- readRDS(FINAL_SELECTION_RDS)
if (is.null(obj$parameters)) {
  stop("Final selection object is missing obj$parameters.")
}
params <- obj$parameters
get_param <- function(name, default) {
  if (!is.null(params[[name]])) params[[name]] else default
}
TARGET_NAME <- get_param("TARGET_NAME", "CGGA_mRNA_array_301")
SOURCE_NAME <- get_param("SOURCE_NAME", "TCGA_LGG")
COVARIATES <- get_param("COVARIATES", c("Age", "Gender"))
BINARY_COVARIATES <- get_param("BINARY_COVARIATES", c("Gender"))
N_GENES_KEEP <- get_param("N_GENES_KEEP", length(obj$keep_genes))
TARGET_BFDR <- get_param("TARGET_BFDR", 0.20)
BASE_ETA <- get_param("BASE_ETA", 0.10)
SOURCE_ETA <- get_param("SOURCE_ETA", 0.20)
BASE_TAU <- get_param("BASE_TAU", 0.50)
FINAL_STR <- get_param("FINAL_STR", 1.50)
TRANSFER_BETA_STRENGTH <- get_param("TRANSFER_BETA_STRENGTH", 0.25)
WEIGHT_MODE <- get_param("WEIGHT_MODE", "single_source_lpml_reported")
LPML_TEMPERATURE <- get_param("LPML_TEMPERATURE", 5)
MAX_ITER <- get_param("MAX_ITER", 100)
BOOT_MAIN <- get_param("BOOT_B", 200)
BOOT_SEED <- get_param("BOOT_SEED", 4100)
EXTENDED_STABILITY_CUTOFF <- get_param("EXTENDED_STABILITY_CUTOFF", 0.50)
CORE_STABILITY_CUTOFF <- get_param("CORE_STABILITY_CUTOFF", 0.70)
SOURCE_SUPPORT_CUTOFF <- get_param("SOURCE_SUPPORT_CUTOFF", 0.50)
TARGET_ENHANCED_STABILITY_CUTOFF <- get_param("TARGET_ENHANCED_STABILITY_CUTOFF", 0.50)
TARGET_ONLY_SOURCE_DISCORDANT_CUTOFF <- get_param(
  "TARGET_ONLY_SOURCE_DISCORDANT_CUTOFF",
  TARGET_ENHANCED_STABILITY_CUTOFF
)
if (is.null(obj$target)) stop("Final object is missing obj$target.")
if (is.null(obj$source)) stop("Final object is missing obj$source.")
if (is.null(obj$keep_genes) || length(obj$keep_genes) == 0) {
  stop("Final object is missing non-empty obj$keep_genes.")
}
target <- obj$target
source_cohort <- obj$source
keep_genes <- unique(as.character(obj$keep_genes))
# ------------------------------------------------------------------------------
# Basic helpers
# ------------------------------------------------------------------------------
strip_M <- function(x) {
  sub("^M_", "", as.character(x))
}
add_M <- function(x) {
  x <- as.character(x)
  ifelse(grepl("^M_", x), x, paste0("M_", x))
}
bt <- function(x) {
  paste0("`", gsub("`", "", x), "`")
}
safe_mean_local <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}
safe_sd_local <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  sd(x)
}
safe_quantile_local <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, prob, na.rm = TRUE, names = FALSE))
}
safe_ratio_local <- function(num, den) {
  if (!is.finite(num) || !is.finite(den) || abs(den) < 1e-8) return(NA_real_)
  num / den
}
jaccard_index <- function(a, b) {
  a <- unique(as.character(a))
  b <- unique(as.character(b))
  u <- union(a, b)
  if (length(u) == 0) return(NA_real_)
  length(intersect(a, b)) / length(u)
}
write_csv_out <- function(x, filename, subdir = SENS_TAB_DIR) {
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(subdir, filename)
  write.csv(x, out, row.names = FALSE, na = "")
  message(">>> Saved: ", out)
  invisible(out)
}
save_plot_out <- function(plot_obj,
                          filename_base,
                          width = 7,
                          height = 5,
                          dpi = 300,
                          subdir = SENS_FIG_DIR) {
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  png_path <- file.path(subdir, paste0(filename_base, ".png"))
  pdf_path <- file.path(subdir, paste0(filename_base, ".pdf"))
  ggplot2::ggsave(
    png_path,
    plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  ggplot2::ggsave(
    pdf_path,
    plot_obj,
    width = width,
    height = height,
    bg = "white"
  )
  message(">>> Saved: ", png_path)
  message(">>> Saved: ", pdf_path)
  invisible(c(png = png_path, pdf = pdf_path))
}
# ------------------------------------------------------------------------------
# Cohort consistency checks
# ------------------------------------------------------------------------------
check_cohort_orientation <- function(cohort, cohort_name = "cohort") {
  if (is.null(cohort$clinical) || is.null(cohort$expr)) {
    stop(cohort_name, " must contain $clinical and $expr.")
  }
  if (!("Sample_ID" %in% colnames(cohort$clinical))) {
    stop(cohort_name, "$clinical is missing Sample_ID.")
  }
  sample_ids <- as.character(cohort$clinical$Sample_ID)
  expr_cols <- colnames(cohort$expr)
  if (is.null(expr_cols)) {
    stop(cohort_name, "$expr must have sample IDs as column names.")
  }
  missing_samples <- setdiff(sample_ids, expr_cols)
  if (length(missing_samples) > 0) {
    stop(
      cohort_name,
      " has clinical samples missing from expression columns. First missing: ",
      paste(head(missing_samples, 10), collapse = ", ")
    )
  }
  invisible(TRUE)
}
check_cohort_orientation(target, "target")
check_cohort_orientation(source_cohort, "source")
missing_keep_target <- setdiff(keep_genes, rownames(target$expr))
if (length(missing_keep_target) > 0) {
  stop(
    "Some keep_genes are missing from target expression rows. First missing: ",
    paste(head(missing_keep_target, 10), collapse = ", ")
  )
}
missing_keep_source <- setdiff(keep_genes, rownames(source_cohort$expr))
if (length(missing_keep_source) > 0) {
  stop(
    "Some keep_genes are missing from source expression rows. First missing: ",
    paste(head(missing_keep_source, 10), collapse = ", ")
  )
}
# ------------------------------------------------------------------------------
# Read final discovery tables
# ------------------------------------------------------------------------------
read_gene_table <- function(path, required = FALSE, label = "gene table") {
  if (!file.exists(path)) {
    if (required) stop("Cannot find ", label, ": ", path)
    return(data.frame())
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  if (nrow(df) == 0) return(df)
  if (!("Gene" %in% colnames(df)) && "Raw_Gene" %in% colnames(df)) {
    df$Gene <- add_M(df$Raw_Gene)
  }
  if (!("Raw_Gene" %in% colnames(df)) && "Gene" %in% colnames(df)) {
    df$Raw_Gene <- strip_M(df$Gene)
  }
  df
}
selection_report_file <- file.path(DISC_DIR, "glioma_bootstrap_selection_report.csv")
target_anchored_file <- if (file.exists(file.path(DISC_DIR, "glioma_final_selected_genes_target_anchored.csv"))) {
  file.path(DISC_DIR, "glioma_final_selected_genes_target_anchored.csv")
} else {
  file.path(DISC_DIR, "glioma_final_selected_genes_manuscript.csv")
}
transfer_supported_file <- file.path(DISC_DIR, "glioma_final_selected_genes_transfer_supported.csv")
source_supported_file <- file.path(DISC_DIR, "glioma_final_selected_genes_source_supported.csv")
core_file <- file.path(DISC_DIR, "glioma_final_selected_genes_core.csv")
targetonly_file <- file.path(DISC_DIR, "glioma_targetonly_stable_genes.csv")
target_discordant_file <- file.path(DISC_DIR, "glioma_targetonly_source_discordant_genes.csv")
selection_report <- read_gene_table(selection_report_file, required = TRUE, label = "selection report")
selected_target_anchored <- read_gene_table(
  target_anchored_file,
  required = TRUE,
  label = "target-anchored reported set"
)
selected_transfer_supported <- read_gene_table(transfer_supported_file)
selected_source_supported <- read_gene_table(source_supported_file)
selected_core <- read_gene_table(core_file)
selected_targetonly <- read_gene_table(targetonly_file)
selected_target_discordant <- read_gene_table(target_discordant_file)
get_gene_vec <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0 || !("Gene" %in% colnames(df))) {
    return(character(0))
  }
  unique(as.character(df$Gene[df$Gene != ""]))
}
available_gene_names <- paste0("M_", keep_genes)
target_anchored_signature <- intersect(add_M(get_gene_vec(selected_target_anchored)), available_gene_names)
transfer_supported_signature <- intersect(add_M(get_gene_vec(selected_transfer_supported)), available_gene_names)
source_supported_signature <- intersect(add_M(get_gene_vec(selected_source_supported)), available_gene_names)
core_signature <- intersect(add_M(get_gene_vec(selected_core)), available_gene_names)
target_only_signature <- intersect(add_M(get_gene_vec(selected_targetonly)), available_gene_names)
target_discordant_signature <- intersect(add_M(get_gene_vec(selected_target_discordant)), available_gene_names)
target_enhanced_signature <- if ("Target_Enhanced" %in% colnames(selection_report)) {
  selection_report$Gene[selection_report$Target_Enhanced %in% TRUE]
} else {
  selection_report$Gene[grepl("Target-enhanced", selection_report$Category)]
}
target_enhanced_signature <- intersect(add_M(target_enhanced_signature), available_gene_names)
if (length(target_anchored_signature) == 0) {
  stop("Target-anchored reported set is empty or has no overlap with keep_genes.")
}
if (length(transfer_supported_signature) == 0) {
  warning("Transfer-supported signature is empty. Some sensitivity summaries will show NA.")
}
message(">>> Loaded final target-anchored glioma object.")
message(">>> Target cohort: ", TARGET_NAME)
message(">>> Source cohort: ", SOURCE_NAME)
message(">>> Candidate genes: ", length(keep_genes))
message(">>> Target-anchored reported set: ", paste(target_anchored_signature, collapse = ", "))
message(">>> Transfer-supported set: ",
        ifelse(length(transfer_supported_signature) == 0, "none", paste(transfer_supported_signature, collapse = ", ")))
message(">>> Source-supported subset: ",
        ifelse(length(source_supported_signature) == 0, "none", paste(source_supported_signature, collapse = ", ")))
message(">>> Core subset: ",
        ifelse(length(core_signature) == 0, "none", paste(core_signature, collapse = ", ")))
message(">>> Target-only stable genes: ",
        ifelse(length(target_only_signature) == 0, "none", paste(target_only_signature, collapse = ", ")))
message(">>> Target-only source-discordant genes: ",
        ifelse(length(target_discordant_signature) == 0, "none", paste(target_discordant_signature, collapse = ", ")))
# ------------------------------------------------------------------------------
# Baseline source PIP and target stability vectors
# ------------------------------------------------------------------------------
get_named_vector_from_report <- function(df, value_col, default = NA_real_) {
  out <- rep(default, length(available_gene_names))
  names(out) <- available_gene_names
  if (!("Gene" %in% colnames(df)) || !(value_col %in% colnames(df))) {
    return(out)
  }
  common <- intersect(available_gene_names, as.character(df$Gene))
  out[common] <- as.numeric(df[[value_col]][match(common, df$Gene)])
  out
}
baseline_target_stability <- get_named_vector_from_report(
  selection_report,
  "TargetOnly_Stability",
  default = 0
)
baseline_vita_stability <- get_named_vector_from_report(
  selection_report,
  "VITA_Stability",
  default = 0
)
baseline_source_pip <- get_named_vector_from_report(
  selection_report,
  "Source_PIP_Weighted",
  default = NA_real_
)
baseline_source_beta <- get_named_vector_from_report(
  selection_report,
  "Source_Beta_Weighted",
  default = NA_real_
)
baseline_target_stability[!is.finite(baseline_target_stability)] <- 0
baseline_vita_stability[!is.finite(baseline_vita_stability)] <- 0
# ------------------------------------------------------------------------------
# Target-anchored rule
# ------------------------------------------------------------------------------
apply_target_anchored_rule <- function(gene_names,
                                       vita_stability,
                                       target_stability = baseline_target_stability,
                                       source_pip = baseline_source_pip,
                                       source_support_cutoff = SOURCE_SUPPORT_CUTOFF,
                                       extended_cutoff = EXTENDED_STABILITY_CUTOFF,
                                       core_cutoff = CORE_STABILITY_CUTOFF,
                                       targetonly_cutoff = TARGET_ONLY_SOURCE_DISCORDANT_CUTOFF) {
  gene_names <- unique(as.character(gene_names))
  vita_stability <- vita_stability[gene_names]
  target_stability <- target_stability[gene_names]
  source_pip <- source_pip[gene_names]
  vita_stability[!is.finite(vita_stability)] <- 0
  target_stability[!is.finite(target_stability)] <- 0
  source_supported <- is.finite(source_pip) & source_pip >= source_support_cutoff
  vita_extended <- vita_stability >= extended_cutoff
  vita_core <- vita_stability >= core_cutoff
  targetonly_stable <- target_stability >= targetonly_cutoff
  source_supported_stable <- vita_extended & source_supported
  target_enhanced <- targetonly_stable & vita_extended & !source_supported
  target_discordant <- targetonly_stable & !source_supported & !vita_extended
  transfer_supported <- source_supported_stable | target_enhanced
  target_anchored <- transfer_supported | target_discordant
  category <- rep("Other", length(gene_names))
  category[source_supported & !vita_extended] <- "Filtered source-supported signal"
  category[target_discordant] <- "Target-only source-discordant"
  category[target_enhanced] <- "Target-enhanced transfer-retained"
  category[source_supported_stable] <- "Source-supported transfer-stabilized"
  category[vita_core & source_supported] <- "High-stability transfer core"
  data.frame(
    Gene = gene_names,
    Raw_Gene = strip_M(gene_names),
    VITA_Stability = as.numeric(vita_stability),
    TargetOnly_Stability = as.numeric(target_stability),
    Source_PIP_Weighted = as.numeric(source_pip),
    Source_Supported = source_supported,
    VITA_Extended = vita_extended,
    VITA_Core = vita_core,
    TargetOnly_Stable = targetonly_stable,
    SourceSupported_Stable = source_supported_stable,
    Target_Enhanced = target_enhanced,
    TargetOnly_SourceDiscordant = target_discordant,
    TransferSupported_Signature = transfer_supported,
    TargetAnchored_Reported = target_anchored,
    Category = category,
    stringsAsFactors = FALSE
  )
}
summarise_classification_recovery <- function(class_df,
                                              label = "analysis",
                                              cindex = NA_real_,
                                              cindex_sd = NA_real_) {
  target_anchored_selected <- class_df$Gene[class_df$TargetAnchored_Reported]
  transfer_selected <- class_df$Gene[class_df$TransferSupported_Signature]
  source_supported_selected <- class_df$Gene[class_df$SourceSupported_Stable]
  core_selected <- class_df$Gene[class_df$VITA_Core & class_df$Source_Supported]
  targetonly_selected <- class_df$Gene[class_df$TargetOnly_Stable]
  discordant_selected <- class_df$Gene[class_df$TargetOnly_SourceDiscordant]
  data.frame(
    Analysis_Label = label,
    N_TargetAnchored_Selected = length(target_anchored_selected),
    N_TransferSupported_Selected = length(transfer_selected),
    N_SourceSupported_Selected = length(source_supported_selected),
    N_Core_Selected = length(core_selected),
    N_TargetOnlyStable_Selected = length(targetonly_selected),
    N_TargetDiscordant_Selected = length(discordant_selected),
    TargetAnchored_Recall = ifelse(
      length(target_anchored_signature) > 0,
      length(intersect(target_anchored_signature, target_anchored_selected)) / length(target_anchored_signature),
      NA_real_
    ),
    TransferSupported_Recall = ifelse(
      length(transfer_supported_signature) > 0,
      length(intersect(transfer_supported_signature, transfer_selected)) / length(transfer_supported_signature),
      NA_real_
    ),
    SourceSupported_Recall = ifelse(
      length(source_supported_signature) > 0,
      length(intersect(source_supported_signature, source_supported_selected)) / length(source_supported_signature),
      NA_real_
    ),
    Core_Recall = ifelse(
      length(core_signature) > 0,
      length(intersect(core_signature, core_selected)) / length(core_signature),
      NA_real_
    ),
    TargetOnlyStable_Recall = ifelse(
      length(target_only_signature) > 0,
      length(intersect(target_only_signature, targetonly_selected)) / length(target_only_signature),
      NA_real_
    ),
    TargetDiscordant_Retention = ifelse(
      length(target_discordant_signature) > 0,
      length(intersect(target_discordant_signature, target_anchored_selected)) / length(target_discordant_signature),
      NA_real_
    ),
    Jaccard_vs_TargetAnchored = jaccard_index(target_anchored_signature, target_anchored_selected),
    Jaccard_vs_TransferSupported = jaccard_index(transfer_supported_signature, transfer_selected),
    SourceSupported_Fraction_TargetAnchored = ifelse(
      length(target_anchored_selected) > 0,
      mean(class_df$Source_Supported[match(target_anchored_selected, class_df$Gene)], na.rm = TRUE),
      NA_real_
    ),
    Mean_C_Index = cindex,
    SD_C_Index = cindex_sd,
    TargetAnchored_Genes = ifelse(
      length(target_anchored_selected) > 0,
      paste(target_anchored_selected, collapse = " | "),
      "None"
    ),
    TransferSupported_Genes = ifelse(
      length(transfer_selected) > 0,
      paste(transfer_selected, collapse = " | "),
      "None"
    ),
    TargetDiscordant_Genes = ifelse(
      length(discordant_selected) > 0,
      paste(discordant_selected, collapse = " | "),
      "None"
    ),
    stringsAsFactors = FALSE
  )
}
# ------------------------------------------------------------------------------
# Cohort manipulation helpers
# ------------------------------------------------------------------------------
subset_cohort_by_rows <- function(cohort, idx, suffix = "subset") {
  clin <- cohort$clinical[idx, , drop = FALSE]
  sample_ids <- as.character(clin$Sample_ID)
  expr <- cohort$expr[, sample_ids, drop = FALSE]
  list(
    clinical = clin,
    expr = expr,
    cohort = paste0(cohort$cohort, "_", suffix)
  )
}
subsample_source <- function(source_local,
                             dose,
                             seed = 1,
                             min_n = 20) {
  if (!is.finite(dose) || dose <= 0 || dose > 1) {
    stop("dose must be in (0, 1].")
  }
  set.seed(seed)
  n <- nrow(source_local$clinical)
  n_sample <- max(min_n, floor(n * dose))
  n_sample <- min(n_sample, n)
  idx <- sample(seq_len(n), size = n_sample, replace = FALSE)
  subset_cohort_by_rows(
    source_local,
    idx = idx,
    suffix = paste0("dose_", round(dose * 100))
  )
}
add_noise_genes_to_cohort <- function(cohort,
                                      fake_names,
                                      seed = 1) {
  set.seed(seed)
  expr <- as.matrix(cohort$expr)
  sample_ids <- colnames(expr)
  fake_mat <- matrix(
    rnorm(length(fake_names) * length(sample_ids)),
    nrow = length(fake_names),
    ncol = length(sample_ids)
  )
  rownames(fake_mat) <- fake_names
  colnames(fake_mat) <- sample_ids
  cohort$expr <- rbind(expr, fake_mat)
  cohort
}
# ------------------------------------------------------------------------------
# Model refit helpers: single-source final workflow
# ------------------------------------------------------------------------------
fit_source_for_current <- function(source_local,
                                   keep_genes_local,
                                   source_eta = SOURCE_ETA,
                                   tau_sq = BASE_TAU,
                                   max_iter = MAX_ITER) {
  source_design_current <- prepare_source_design(
    source_cohort = source_local,
    keep_genes = keep_genes_local,
    covariates = COVARIATES,
    binary_covariates = BINARY_COVARIATES
  )
  source_fit_current <- fit_source_model(
    source_design = source_design_current,
    source_name = SOURCE_NAME,
    source_eta0 = source_eta,
    tau_sq = tau_sq,
    max_iter = max_iter,
    verbose = FALSE
  )
  list(
    source_design = source_design_current,
    source_fit = source_fit_current
  )
}
fit_vita_full_current <- function(target_local,
                                  keep_genes_local,
                                  source_design_current,
                                  source_fit_current,
                                  eta0 = BASE_ETA,
                                  source_eta0 = SOURCE_ETA,
                                  strength = FINAL_STR,
                                  tau_sq = BASE_TAU,
                                  transfer_beta_strength = TRANSFER_BETA_STRENGTH,
                                  weight_mode = WEIGHT_MODE,
                                  lpml_temperature = LPML_TEMPERATURE,
                                  max_iter = MAX_ITER) {
  target_design <- prepare_design(
    cohort = target_local,
    keep_genes = keep_genes_local,
    covariates = COVARIATES,
    binary_covariates = BINARY_COVARIATES,
    fit_scalers = TRUE
  )
  fit <- fit_target_with_transfer(
    target_design = target_design,
    source_design = source_design_current,
    source_name = SOURCE_NAME,
    source_fit = source_fit_current,
    transfer = TRUE,
    eta0 = eta0,
    source_eta0 = source_eta0,
    strength = strength,
    tau_sq = tau_sq,
    transfer_beta_strength = transfer_beta_strength,
    weight_mode = weight_mode,
    lpml_temperature = lpml_temperature,
    max_iter = max_iter,
    verbose = FALSE
  )
  
  list(
    target_design = target_design,
    fit = fit
  )
}
run_vita_bootstrap_selection <- function(target_local,
                                         keep_genes_local,
                                         source_design_current,
                                         source_fit_current,
                                         boot_b = 100,
                                         target_bfdr = TARGET_BFDR,
                                         eta0 = BASE_ETA,
                                         source_eta0 = SOURCE_ETA,
                                         strength = FINAL_STR,
                                         tau_sq = BASE_TAU,
                                         transfer_beta_strength = TRANSFER_BETA_STRENGTH,
                                         weight_mode = WEIGHT_MODE,
                                         lpml_temperature = LPML_TEMPERATURE,
                                         seed_base = 7000,
                                         store_pip = TRUE,
                                         max_iter = MAX_ITER,
                                         show_progress = TRUE) {
  n_target <- nrow(target_local$clinical)
  full_idx <- seq_len(n_target)
  reference_genes <- paste0("M_", keep_genes_local)
  n_genes <- length(reference_genes)
  boot_mat <- matrix(0, nrow = boot_b, ncol = n_genes)
  colnames(boot_mat) <- reference_genes
  pip_mat <- NULL
  if (store_pip) {
    pip_mat <- matrix(NA_real_, nrow = boot_b, ncol = n_genes)
    colnames(pip_mat) <- reference_genes
  }
  diag_rows <- list()
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = boot_b, style = 3)
  }
  for (b in seq_len(boot_b)) {
    set.seed(seed_base + b)
    boot_idx <- sample(full_idx, replace = TRUE)
    target_b <- make_bootstrap_cohort(
      cohort = target_local,
      boot_idx = boot_idx,
      prefix = paste0("sens_boot_", b)
    )
    d_tab <- table(target_b$clinical$IDH)
    d_has_both <- all(c("0", "1") %in% names(d_tab))
    fit_success <- FALSE
    selected <- character(0)
    err_msg <- ""
    if (d_has_both) {
      target_design_b <- try(
        prepare_design(
          cohort = target_b,
          keep_genes = keep_genes_local,
          covariates = COVARIATES,
          binary_covariates = BINARY_COVARIATES,
          fit_scalers = TRUE
        ),
        silent = TRUE
      )
      if (!inherits(target_design_b, "try-error")) {
        fit_b <- try(
          fit_target_with_transfer(
            target_design = target_design_b,
            source_design = source_design_current,
            source_name = SOURCE_NAME,
            source_fit = source_fit_current,
            transfer = TRUE,
            eta0 = eta0,
            source_eta0 = source_eta0,
            strength = strength,
            tau_sq = tau_sq,
            transfer_beta_strength = transfer_beta_strength,
            weight_mode = weight_mode,
            lpml_temperature = lpml_temperature,
            max_iter = max_iter,
            verbose = FALSE
          ),
          silent = TRUE
        )
        fit_success <- !inherits(fit_b, "try-error")
        if (fit_success) {
          sel_b <- get_selected_genes(
            fit_b$pip_gene,
            target_bfdr = target_bfdr
          )
          selected <- intersect(sel_b$selected, reference_genes)
          if (length(selected) > 0) {
            boot_mat[b, selected] <- 1
          }
          if (store_pip) {
            common <- intersect(reference_genes, names(fit_b$pip_gene))
            pip_mat[b, common] <- fit_b$pip_gene[common]
          }
        } else {
          err_msg <- as.character(attr(fit_b, "condition")$message)
        }
      } else {
        err_msg <- as.character(attr(target_design_b, "condition")$message)
      }
    } else {
      err_msg <- "Bootstrap sample lacks both IDH classes."
    }
    diag_rows[[b]] <- data.frame(
      Bootstrap_ID = b,
      N_Unique_Target = length(unique(boot_idx)),
      D_Has_Both_Levels = d_has_both,
      Fit_Success = fit_success,
      N_Selected = length(selected),
      Error_Message = err_msg,
      stringsAsFactors = FALSE
    )
    if (show_progress) {
      setTxtProgressBar(pb, b)
    }
  }
  if (show_progress) {
    close(pb)
  }
  stability <- colMeans(boot_mat, na.rm = TRUE)
  names(stability) <- reference_genes
  list(
    boot_mat = boot_mat,
    pip_mat = pip_mat,
    stability = stability,
    diagnostics = bind_rows(diag_rows)
  )
}
source_pip_from_fit <- function(fit_obj, gene_names) {
  out <- rep(NA_real_, length(gene_names))
  names(out) <- gene_names
  pip_bar <- NULL
  if (!is.null(fit_obj$priors$source_pip_bar)) {
    pip_bar <- fit_obj$priors$source_pip_bar
  }
  if (is.null(pip_bar)) return(out)
  common <- intersect(gene_names, names(pip_bar))
  out[common] <- pip_bar[common]
  out
}
# ------------------------------------------------------------------------------
# Prediction helpers
# ------------------------------------------------------------------------------
cindex_survival_score <- function(time, status, survival_score) {
  df <- data.frame(
    Time = as.numeric(time),
    Status = as.numeric(status),
    Score = as.numeric(survival_score)
  )
  df <- df[
    is.finite(df$Time) &
      df$Time > 0 &
      df$Status %in% c(0, 1) &
      is.finite(df$Score),
    ,
    drop = FALSE
  ]
  if (nrow(df) < 10) return(NA_real_)
  if (length(unique(df$Status)) < 2) return(NA_real_)
  if (length(unique(df$Score)) < 2) return(NA_real_)
  out <- try(
    survival::concordance(
      survival::Surv(Time, Status) ~ Score,
      data = df
    )$concordance,
    silent = TRUE
  )
  if (inherits(out, "try-error")) return(NA_real_)
  as.numeric(out)
}
fit_cox_cindex_safe <- function(train_df, test_df, vars) {
  vars <- intersect(vars, colnames(train_df))
  vars <- vars[vars %in% colnames(test_df)]
  if (length(vars) == 0) return(NA_real_)
  rhs <- paste(bt(vars), collapse = " + ")
  form <- as.formula(paste("Surv(Time, Status) ~", rhs))
  fit <- try(
    survival::coxph(
      form,
      data = train_df,
      singular.ok = TRUE,
      control = survival::coxph.control(iter.max = 50)
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) return(NA_real_)
  lp <- try(
    as.numeric(predict(fit, newdata = test_df, type = "lp")),
    silent = TRUE
  )
  if (inherits(lp, "try-error") || any(!is.finite(lp))) return(NA_real_)
  cindex_survival_score(test_df$Time, test_df$Status, -lp)
}
evaluate_multiple_fixed_signatures_cox <- function(target_cohort,
                                                   signature_list,
                                                   covariates = COVARIATES,
                                                   binary_covariates = BINARY_COVARIATES,
                                                   repeats = 10,
                                                   folds = 5,
                                                   seed = 2030) {
  signature_list <- signature_list[lengths(signature_list) > 0]
  all_genes <- unique(strip_M(unlist(signature_list, use.names = FALSE)))
  all_genes <- intersect(all_genes, rownames(target_cohort$expr))
  clin_vars <- c("D", covariates)
  rows <- list()
  row_id <- 1
  set.seed(seed)
  for (r in seq_len(repeats)) {
    fold_list <- make_stratified_folds(
      target_cohort$clinical$OS_status,
      k = folds,
      seed = seed + r
    )
    for (k in seq_len(folds)) {
      test_idx <- fold_list[[k]]
      train_idx <- setdiff(seq_len(nrow(target_cohort$clinical)), test_idx)
      dat <- make_cox_fold_data_target(
        target_cohort = target_cohort,
        idx_train = train_idx,
        idx_test = test_idx,
        raw_genes = all_genes,
        covariates = covariates,
        binary_covariates = binary_covariates
      )
      train_df <- dat$train
      test_df <- dat$test
      c_clin <- fit_cox_cindex_safe(
        train_df = train_df,
        test_df = test_df,
        vars = clin_vars
      )
      rows[[row_id]] <- data.frame(
        Repeat = r,
        Fold = k,
        Model = "Clinical-Only",
        N_Genes = 0,
        Training_Domain = "Target",
        C_Index = c_clin,
        Delta_vs_Clinical = 0,
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1
      for (nm in names(signature_list)) {
        raw_genes <- strip_M(signature_list[[nm]])
        raw_genes <- intersect(raw_genes, rownames(target_cohort$expr))
        vars <- c(clin_vars, paste0("M_", raw_genes))
        c_sig <- fit_cox_cindex_safe(
          train_df = train_df,
          test_df = test_df,
          vars = vars
        )
        rows[[row_id]] <- data.frame(
          Repeat = r,
          Fold = k,
          Model = nm,
          N_Genes = length(raw_genes),
          Training_Domain = "Target",
          C_Index = c_sig,
          Delta_vs_Clinical = c_sig - c_clin,
          stringsAsFactors = FALSE
        )
        row_id <- row_id + 1
      }
    }
  }
  repeated_cv <- bind_rows(rows)
  summary <- repeated_cv %>%
    group_by(Model, N_Genes, Training_Domain) %>%
    summarise(
      RepeatedCV_N = sum(is.finite(C_Index)),
      RepeatedCV_Mean_C_Index = mean(C_Index, na.rm = TRUE),
      RepeatedCV_SD_C_Index = sd(C_Index, na.rm = TRUE),
      Mean_Delta_vs_Clinical = mean(Delta_vs_Clinical, na.rm = TRUE),
      Positive_Delta_N = sum(Delta_vs_Clinical > 0, na.rm = TRUE),
      Positive_Delta_Rate = mean(Delta_vs_Clinical > 0, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(RepeatedCV_Mean_C_Index))
  list(
    repeated_cv = repeated_cv,
    summary = summary
  )
}
extract_complete_beta <- function(fit_obj) {
  beta <- NULL
  if (!is.null(fit_obj$beta_all)) {
    beta <- fit_obj$beta_all
  } else if (!is.null(fit_obj$fit) && !is.null(fit_obj$fit$beta_hat)) {
    beta <- fit_obj$fit$beta_hat
  } else if (!is.null(fit_obj$fit_target) &&
             !is.null(fit_obj$fit_target$beta_hat)) {
    beta <- fit_obj$fit_target$beta_hat
  } else if (!is.null(fit_obj$beta_hat)) {
    beta <- fit_obj$beta_hat
  }
  if (is.null(beta)) return(NULL)
  beta <- as.numeric(beta)
  if (any(!is.finite(beta))) return(NULL)
  beta
}
predict_aft_cindex_safe <- function(fit_obj, test_design) {
  if (exists("predict_cindex_aft", mode = "function")) {
    out <- try(
      predict_cindex_aft(fit_obj, test_design),
      silent = TRUE
    )
    if (!inherits(out, "try-error") && is.finite(out)) {
      return(as.numeric(out))
    }
  }
  if (is.null(test_design$X)) return(NA_real_)
  x_test <- try(as.matrix(test_design$X), silent = TRUE)
  if (inherits(x_test, "try-error")) return(NA_real_)
  suppressWarnings(storage.mode(x_test) <- "double")
  x_test[!is.finite(x_test)] <- 0
  beta <- extract_complete_beta(fit_obj)
  if (is.null(beta) || ncol(x_test) != length(beta)) return(NA_real_)
  lp <- try(as.numeric(x_test %*% beta), silent = TRUE)
  if (inherits(lp, "try-error") || any(!is.finite(lp))) return(NA_real_)
  aft_cindex(test_design$Time, test_design$Status, lp)
}
evaluate_vita_aft_cv_current <- function(target_local,
                                         keep_genes_local,
                                         source_design_current,
                                         source_fit_current,
                                         n_folds = 5,
                                         seed = 8000,
                                         eta0 = BASE_ETA,
                                         source_eta0 = SOURCE_ETA,
                                         strength = FINAL_STR,
                                         tau_sq = BASE_TAU,
                                         transfer_beta_strength = TRANSFER_BETA_STRENGTH,
                                         weight_mode = WEIGHT_MODE,
                                         lpml_temperature = LPML_TEMPERATURE,
                                         max_iter = MAX_ITER) {
  folds <- make_stratified_folds(
    target_local$clinical$OS_status,
    k = n_folds,
    seed = seed
  )
  c_vita <- rep(NA_real_, n_folds)
  keep_genes_local <- intersect(keep_genes_local, rownames(target_local$expr))
  for (k in seq_len(n_folds)) {
    test_idx <- folds[[k]]
    train_idx <- setdiff(seq_len(nrow(target_local$clinical)), test_idx)
    dat <- try(
      prepare_target_train_test(
        target_cohort = target_local,
        train_idx = train_idx,
        test_idx = test_idx,
        keep_genes = keep_genes_local,
        covariates = COVARIATES,
        binary_covariates = BINARY_COVARIATES
      ),
      silent = TRUE
    )
    if (inherits(dat, "try-error")) next
    fit_k <- try(
      fit_target_with_transfer(
        target_design = dat$train,
        source_design = source_design_current,
        source_name = SOURCE_NAME,
        source_fit = source_fit_current,
        transfer = TRUE,
        eta0 = eta0,
        source_eta0 = source_eta0,
        strength = strength,
        tau_sq = tau_sq,
        transfer_beta_strength = transfer_beta_strength,
        weight_mode = weight_mode,
        lpml_temperature = lpml_temperature,
        max_iter = max_iter,
        verbose = FALSE
      ),
      silent = TRUE
    )
    if (inherits(fit_k, "try-error")) next
    c_vita[k] <- predict_aft_cindex_safe(
      fit_obj = fit_k,
      test_design = dat$test
    )
  }
  c_vita
}
build_signature_analysis_data <- function(target_cohort,
                                          signature_genes,
                                          covariates = COVARIATES,
                                          binary_covariates = BINARY_COVARIATES) {
  raw_genes <- strip_M(signature_genes)
  raw_genes <- intersect(raw_genes, rownames(target_cohort$expr))
  if (length(raw_genes) == 0) {
    stop("No signature genes are available in target expression matrix.")
  }
  full_idx <- seq_len(nrow(target_cohort$clinical))
  dat <- make_cox_fold_data_target(
    target_cohort = target_cohort,
    idx_train = full_idx,
    idx_test = full_idx,
    raw_genes = raw_genes,
    covariates = covariates,
    binary_covariates = binary_covariates
  )
  list(
    data = dat$train,
    clin_vars = c("D", covariates),
    signature_genes = paste0("M_", raw_genes)
  )
}
# ------------------------------------------------------------------------------
# Source-target feature diagnostics
# ------------------------------------------------------------------------------
compute_source_target_pca <- function(raw_genes = keep_genes,
                                      scale_mode = c(
                                        "pooled_zscore",
                                        "within_cohort_zscore"
                                      )) {
  scale_mode <- match.arg(scale_mode)
  raw_genes <- intersect(raw_genes, rownames(target$expr))
  raw_genes <- intersect(raw_genes, rownames(source_cohort$expr))
  if (length(raw_genes) < 2) {
    stop("At least two common genes are required for PCA diagnostics.")
  }
  target_expr <- t(target$expr[raw_genes, target$clinical$Sample_ID, drop = FALSE])
  source_expr <- t(source_cohort$expr[raw_genes, source_cohort$clinical$Sample_ID, drop = FALSE])
  target_expr <- as.matrix(target_expr)
  source_expr <- as.matrix(source_expr)
  storage.mode(target_expr) <- "double"
  storage.mode(source_expr) <- "double"
  dataset_label <- c(
    rep(paste0("Target (", TARGET_NAME, ")"), nrow(target_expr)),
    rep(paste0("Source (", SOURCE_NAME, ")"), nrow(source_expr))
  )
  if (scale_mode == "pooled_zscore") {
    # Raw feature-space diagnostic:
    # target and source are pooled first, then gene-wise z-scored together.
    # This intentionally displays platform-level separation if it exists.
    combined_X <- rbind(target_expr, source_expr)
    combined_sc <- safe_scale_train(combined_X)
    X_for_pca <- combined_sc$X
    scale_label <- "Raw pooled z-score"
  } else {
    # Within-cohort standardized diagnostic:
    # target and source are z-scored separately before PCA. This removes
    # first-order platform mean/scale differences and better reflects the
    # domain-specific standardization used before target-stage modeling.
    target_sc <- safe_scale_train(target_expr)$X
    source_sc <- safe_scale_train(source_expr)$X
    X_for_pca <- rbind(target_sc, source_sc)
    X_for_pca[!is.finite(X_for_pca)] <- 0
    scale_label <- "Within-cohort z-score"
  }
  pca_res <- prcomp(X_for_pca, center = FALSE, scale. = FALSE)
  var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 2)
  data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    Dataset = dataset_label,
    PC1_Var = var_exp[1],
    PC2_Var = var_exp[2],
    N_Genes = length(raw_genes),
    Scaling_Mode = scale_label,
    stringsAsFactors = FALSE
  )
}
# ------------------------------------------------------------------------------
# Baseline classification
# ------------------------------------------------------------------------------
baseline_classification <- apply_target_anchored_rule(
  gene_names = available_gene_names,
  vita_stability = baseline_vita_stability,
  target_stability = baseline_target_stability,
  source_pip = baseline_source_pip,
  source_support_cutoff = SOURCE_SUPPORT_CUTOFF,
  extended_cutoff = EXTENDED_STABILITY_CUTOFF,
  core_cutoff = CORE_STABILITY_CUTOFF,
  targetonly_cutoff = TARGET_ONLY_SOURCE_DISCORDANT_CUTOFF
)
write_csv_out(
  baseline_classification,
  "Baseline_TargetAnchored_Classification.csv"
)
baseline_recovery <- summarise_classification_recovery(
  baseline_classification,
  label = "Baseline final discovery"
)
write_csv_out(
  baseline_recovery,
  "Baseline_TargetAnchored_Recovery.csv"
)