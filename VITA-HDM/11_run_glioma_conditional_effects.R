# ==============================================================================
# 11_run_glioma_conditional_effects.R
# Run after file 09: Rscript 11_run_glioma_conditional_effects.R
# Conditional bootstrap mediation-effect quantification for the glioma
# VITA-HDM example.
# Analysis specification:
#   - Reads the target-anchored reported set from
#     09_run_glioma_discovery.R.
#   - Computes conditional bootstrap effects for:
#       1. Target-anchored reported set
#       2. Transfer-supported subset
#       3. Source-supported subset
#       4. Target-enhanced subset
#       5. Target-only stable signature
#       6. Target-only source-discordant subset
#   - Gene-level effects are reported for all genes in the union of these sets.
# Important:
#   This is a conditional bootstrap for fixed selected mediator sets.
#   It is NOT a full post-selection bootstrap of the entire discovery pipeline.
# Outputs:
#   results/case_study/conditional_effects/
# No sink() logging is used.
# ==============================================================================
script_args <- commandArgs(trailingOnly = FALSE)
script_file_arg <- grep("^--file=", script_args, value = TRUE)
source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
SCRIPT_DIR <- if (length(script_file_arg) > 0) {
  dirname(normalizePath(
    sub("^--file=", "", script_file_arg[[1]]),
    winslash = "/",
    mustWork = FALSE
  ))
} else if (!is.null(source_file)) {
  dirname(normalizePath(
    source_file,
    winslash = "/",
    mustWork = FALSE
  ))
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}
source(file.path(SCRIPT_DIR, "00_project_setup.R"))
PROJECT_DIR <- find_project_root(SCRIPT_DIR)
PATHS <- project_paths(PROJECT_DIR)
suppressPackageStartupMessages({
  library(dplyr)
})
source(file.path(PROJECT_DIR, "08_glioma_vita_functions.R"))
# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
RESULTS_ROOT <- PATHS$case_results
SELECTION_RDS <- file.path(
  PATHS$discovery, "glioma_final_selection_objects.rds"
)
TARGET_ANCHORED_CSV <- file.path(
  PATHS$discovery, "glioma_final_selected_genes_target_anchored.csv"
)
MANUSCRIPT_CSV <- file.path(
  PATHS$discovery, "glioma_final_selected_genes_manuscript.csv"
)
TRANSFER_SUPPORTED_CSV <- file.path(
  PATHS$discovery, "glioma_final_selected_genes_transfer_supported.csv"
)
TARGETONLY_CSV <- file.path(
  PATHS$discovery, "glioma_targetonly_stable_genes.csv"
)
TARGET_DISCORDANT_CSV <- file.path(
  PATHS$discovery, "glioma_targetonly_source_discordant_genes.csv"
)
SELECTION_REPORT_CSV <- file.path(
  PATHS$discovery, "glioma_bootstrap_selection_report.csv"
)
OUT_DIR <- ensure_directory(PATHS$conditional_effects)
if (!file.exists(SELECTION_RDS)) {
  stop(
    "Cannot find selection object. Run 09_run_glioma_discovery.R first: ",
    SELECTION_RDS
  )
}
if (!file.exists(TARGET_ANCHORED_CSV) && !file.exists(MANUSCRIPT_CSV)) {
  stop(
    "Cannot find target-anchored or manuscript signature CSV. ",
    "Run 09_run_glioma_discovery.R first."
  )
}
# ------------------------------------------------------------------------------
# Settings
# ------------------------------------------------------------------------------
BOOT_EFFECT_B <- 200
BOOT_SEED <- 5000
# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------
read_selected_gene_df <- function(path, fallback_category = "Selected") {
  if (!file.exists(path)) {
    return(data.frame())
  }
  df <- read.csv(path, stringsAsFactors = FALSE)
  if (nrow(df) == 0) {
    return(data.frame())
  }
  if (!("Gene" %in% colnames(df))) {
    stop("Selection file must contain Gene column: ", path)
  }
  if (!("Raw_Gene" %in% colnames(df))) {
    df$Raw_Gene <- sub("^M_", "", df$Gene)
  }
  if (!("Category" %in% colnames(df))) {
    df$Category <- fallback_category
  }
  df
}
gene_vec_from_df <- function(df) {
  if (is.null(df) || nrow(df) == 0 || !("Gene" %in% colnames(df))) {
    return(character(0))
  }
  unique(as.character(df$Gene[df$Gene != ""]))
}
safe_gene_intersect <- function(x, available) {
  x <- unique(as.character(x))
  x <- x[x != ""]
  intersect(x, available)
}
# ------------------------------------------------------------------------------
# Load objects and selected sets
# ------------------------------------------------------------------------------
obj <- readRDS(SELECTION_RDS)
target <- obj$target
keep_genes <- obj$keep_genes
source_design <- obj$source_design
source_fit <- obj$source_fit
params <- obj$parameters
COVARIATES <- params$COVARIATES
BINARY_COVARIATES <- params$BINARY_COVARIATES
target_anchored_path <- if (file.exists(TARGET_ANCHORED_CSV)) {
  TARGET_ANCHORED_CSV
} else {
  MANUSCRIPT_CSV
}
target_anchored_df <- read_selected_gene_df(
  target_anchored_path,
  fallback_category = "Target-anchored reported"
)
transfer_supported_df <- read_selected_gene_df(
  TRANSFER_SUPPORTED_CSV,
  fallback_category = "Transfer-supported"
)
targetonly_df <- read_selected_gene_df(
  TARGETONLY_CSV,
  fallback_category = "Target-only stable"
)
target_discordant_df <- read_selected_gene_df(
  TARGET_DISCORDANT_CSV,
  fallback_category = "Target-only source-discordant"
)
selection_report <- if (file.exists(SELECTION_REPORT_CSV)) {
  read.csv(SELECTION_REPORT_CSV, stringsAsFactors = FALSE)
} else {
  data.frame()
}
target_anchored_genes <- gene_vec_from_df(target_anchored_df)
transfer_supported_genes <- gene_vec_from_df(transfer_supported_df)
targetonly_genes <- gene_vec_from_df(targetonly_df)
target_discordant_genes <- gene_vec_from_df(target_discordant_df)
if (length(target_anchored_genes) == 0) {
  stop("No target-anchored selected genes found.")
}
selected_union <- unique(c(
  target_anchored_genes,
  transfer_supported_genes,
  targetonly_genes,
  target_discordant_genes
))
cat("\n============================================================\n")
cat("Conditional mediation-effect bootstrap: target-anchored version\n")
cat("============================================================\n")
cat("Target-anchored reported genes:",
    paste(sub("^M_", "", target_anchored_genes), collapse = ", "), "\n")
cat("Transfer-supported genes:",
    ifelse(length(transfer_supported_genes) == 0, "none",
           paste(sub("^M_", "", transfer_supported_genes), collapse = ", ")), "\n")
cat("Target-only stable genes:",
    ifelse(length(targetonly_genes) == 0, "none",
           paste(sub("^M_", "", targetonly_genes), collapse = ", ")), "\n")
cat("Target-only source-discordant genes:",
    ifelse(length(target_discordant_genes) == 0, "none",
           paste(sub("^M_", "", target_discordant_genes), collapse = ", ")), "\n")
# ------------------------------------------------------------------------------
# Full-data fit and point estimates
# ------------------------------------------------------------------------------
target_full <- prepare_design(
  cohort = target,
  keep_genes = keep_genes,
  covariates = COVARIATES,
  binary_covariates = BINARY_COVARIATES,
  fit_scalers = TRUE
)
fit_full <- fit_target_with_transfer(
  target_design = target_full,
  source_design = source_design,
  source_name = params$SOURCE_NAME,
  source_fit = source_fit,
  transfer = TRUE,
  eta0 = params$BASE_ETA,
  source_eta0 = params$SOURCE_ETA,
  strength = params$FINAL_STR,
  tau_sq = params$BASE_TAU,
  transfer_beta_strength = params$TRANSFER_BETA_STRENGTH,
  weight_mode = params$WEIGHT_MODE,
  lpml_temperature = params$LPML_TEMPERATURE,
  max_iter = params$MAX_ITER,
  verbose = FALSE
)
selected_union <- safe_gene_intersect(
  selected_union,
  names(fit_full$beta_gene)
)
if (length(selected_union) == 0) {
  stop("Selected genes are not present in final beta vector.")
}
alpha_full <- compute_alpha_paths_safe(
  target_design = target_full,
  gene_names = selected_union,
  covariates = COVARIATES
)
beta_full <- fit_full$beta_gene[selected_union]
pip_full <- fit_full$pip_gene[selected_union]
ide_full <- alpha_full[selected_union] * beta_full[selected_union]
direct_effect <- as.numeric(fit_full$fixed_beta["D"])
# ------------------------------------------------------------------------------
# Conditional bootstrap
# ------------------------------------------------------------------------------
n_target <- nrow(target$clinical)
full_idx <- seq_len(n_target)
boot_ide <- matrix(
  NA_real_,
  nrow = BOOT_EFFECT_B,
  ncol = length(selected_union)
)
colnames(boot_ide) <- selected_union
boot_direct <- rep(NA_real_, BOOT_EFFECT_B)
boot_diag_rows <- list()
cat("\nRunning conditional bootstrap: B =", BOOT_EFFECT_B, "\n")
pb <- txtProgressBar(min = 0, max = BOOT_EFFECT_B, style = 3)
for (b in seq_len(BOOT_EFFECT_B)) {
  set.seed(BOOT_SEED + b)
  boot_idx <- sample(full_idx, replace = TRUE)
  target_b <- make_bootstrap_cohort(
    target,
    boot_idx,
    prefix = paste0("effect_boot", b)
  )
  d_tab <- table(target_b$clinical$IDH)
  d_has_both <- all(c("0", "1") %in% names(d_tab))
  fit_success <- FALSE
  n_estimated <- 0L
  err_msg <- ""
  if (d_has_both) {
    target_b_design <- try(
      prepare_design(
        cohort = target_b,
        keep_genes = keep_genes,
        covariates = COVARIATES,
        binary_covariates = BINARY_COVARIATES,
        fit_scalers = TRUE
      ),
      silent = TRUE
    )
    if (!inherits(target_b_design, "try-error")) {
      fit_b <- try(
        fit_target_with_transfer(
          target_design = target_b_design,
          source_design = source_design,
          source_name = params$SOURCE_NAME,
          source_fit = source_fit,
          transfer = TRUE,
          eta0 = params$BASE_ETA,
          source_eta0 = params$SOURCE_ETA,
          strength = params$FINAL_STR,
          tau_sq = params$BASE_TAU,
          transfer_beta_strength = params$TRANSFER_BETA_STRENGTH,
          weight_mode = params$WEIGHT_MODE,
          lpml_temperature = params$LPML_TEMPERATURE,
          max_iter = params$MAX_ITER,
          verbose = FALSE
        ),
        silent = TRUE
      )
      if (!inherits(fit_b, "try-error")) {
        fit_success <- TRUE
        common_selected <- intersect(selected_union, names(fit_b$beta_gene))
        if (length(common_selected) > 0) {
          alpha_b <- compute_alpha_paths_safe(
            target_design = target_b_design,
            gene_names = common_selected,
            covariates = COVARIATES
          )
          beta_b <- fit_b$beta_gene[common_selected]
          boot_ide[b, common_selected] <-
            alpha_b[common_selected] * beta_b[common_selected]
          n_estimated <- sum(is.finite(boot_ide[b, common_selected]))
        }
        boot_direct[b] <- as.numeric(fit_b$fixed_beta["D"])
      } else {
        err_msg <- as.character(attr(fit_b, "condition")$message)
      }
    } else {
      err_msg <- as.character(attr(target_b_design, "condition")$message)
    }
  } else {
    err_msg <- "Bootstrap sample lacks both IDH classes."
  }
  boot_diag_rows[[b]] <- data.frame(
    Bootstrap = b,
    Has_Both_IDH_Classes = d_has_both,
    Fit_Success = fit_success,
    N_Selected_Genes_Estimated = n_estimated,
    Error_Message = err_msg,
    stringsAsFactors = FALSE
  )
  setTxtProgressBar(pb, b)
}
close(pb)
boot_diag <- bind_rows(boot_diag_rows)
# ------------------------------------------------------------------------------
# Define effect sets
# ------------------------------------------------------------------------------
source_supported_genes <- character(0)
target_enhanced_genes <- character(0)
if (nrow(selection_report) > 0) {
  if ("SourceSupported_Stable" %in% colnames(selection_report)) {
    source_supported_genes <- selection_report$Gene[
      selection_report$Gene %in% selected_union &
        selection_report$SourceSupported_Stable %in% TRUE
    ]
  } else if ("Category" %in% colnames(selection_report)) {
    source_supported_genes <- selection_report$Gene[
      selection_report$Gene %in% selected_union &
        grepl("Source-supported|High-stability", selection_report$Category)
    ]
  }
  if ("Target_Enhanced" %in% colnames(selection_report)) {
    target_enhanced_genes <- selection_report$Gene[
      selection_report$Gene %in% selected_union &
        selection_report$Target_Enhanced %in% TRUE
    ]
  } else if ("Category" %in% colnames(selection_report)) {
    target_enhanced_genes <- selection_report$Gene[
      selection_report$Gene %in% selected_union &
        grepl("Target-enhanced", selection_report$Category)
    ]
  }
}
if (length(transfer_supported_genes) == 0) {
  transfer_supported_genes <- unique(c(source_supported_genes, target_enhanced_genes))
}
target_anchored_genes <- safe_gene_intersect(target_anchored_genes, selected_union)
transfer_supported_genes <- safe_gene_intersect(transfer_supported_genes, selected_union)
source_supported_genes <- safe_gene_intersect(source_supported_genes, selected_union)
target_enhanced_genes <- safe_gene_intersect(target_enhanced_genes, selected_union)
targetonly_genes <- safe_gene_intersect(targetonly_genes, selected_union)
target_discordant_genes <- safe_gene_intersect(target_discordant_genes, selected_union)
effect_sets <- list(
  TargetAnchored_Reported_Set = target_anchored_genes,
  TransferSupported_Set = transfer_supported_genes,
  Source_Supported_Subset = source_supported_genes,
  Target_Enhanced_Subset = target_enhanced_genes,
  TargetOnly_Stable_Signature = targetonly_genes,
  TargetOnly_SourceDiscordant_Subset = target_discordant_genes
)
effect_sets <- effect_sets[lengths(effect_sets) > 0]
# ------------------------------------------------------------------------------
# Summaries
# ------------------------------------------------------------------------------
summarise_overall <- function(group_name, gene_set) {
  gene_set <- safe_gene_intersect(gene_set, selected_union)
  
  if (length(gene_set) == 0) {
    return(NULL)
  }
  total_ide_point <- sum_finite_or_na(ide_full[gene_set])
  total_effect_point <- direct_effect + total_ide_point
  mp_point <- safe_ratio(total_ide_point, total_effect_point)
  boot_total_ide <- apply(
    boot_ide[, gene_set, drop = FALSE],
    1,
    sum_finite_or_na
  )
  boot_total_effect <- boot_direct + boot_total_ide
  boot_mp <- safe_ratio(boot_total_ide, boot_total_effect)
  data.frame(
    Effect_Set = group_name,
    N_Genes = length(gene_set),
    Parameter = c(
      "Direct_Effect",
      "Total_Indirect_Effect",
      "Total_Effect",
      "Mediation_Proportion"
    ),
    Estimate = c(
      direct_effect,
      total_ide_point,
      total_effect_point,
      mp_point
    ),
    Bootstrap_Mean = c(
      safe_mean(boot_direct),
      safe_mean(boot_total_ide),
      safe_mean(boot_total_effect),
      safe_mean(boot_mp)
    ),
    Bootstrap_SD = c(
      safe_sd(boot_direct),
      safe_sd(boot_total_ide),
      safe_sd(boot_total_effect),
      safe_sd(boot_mp)
    ),
    CI_95_Low = c(
      safe_quantile(boot_direct, 0.025),
      safe_quantile(boot_total_ide, 0.025),
      safe_quantile(boot_total_effect, 0.025),
      safe_quantile(boot_mp, 0.025)
    ),
    CI_95_Upp = c(
      safe_quantile(boot_direct, 0.975),
      safe_quantile(boot_total_ide, 0.975),
      safe_quantile(boot_total_effect, 0.975),
      safe_quantile(boot_mp, 0.975)
    ),
    N_Boot_Used = c(
      sum(is.finite(boot_direct)),
      sum(is.finite(boot_total_ide)),
      sum(is.finite(boot_total_effect)),
      sum(is.finite(boot_mp))
    ),
    stringsAsFactors = FALSE
  )
}
overall_effects <- bind_rows(
  lapply(
    names(effect_sets),
    function(nm) summarise_overall(nm, effect_sets[[nm]])
  )
)
overall_effects$CI_Crosses_Zero <- with(
  overall_effects,
  CI_95_Low <= 0 & CI_95_Upp >= 0
)
# ------------------------------------------------------------------------------
# Gene-level effects
# ------------------------------------------------------------------------------
get_gene_category <- function(g) {
  if (nrow(selection_report) > 0 &&
      "Gene" %in% colnames(selection_report) &&
      "Category" %in% colnames(selection_report)) {
    idx <- match(g, selection_report$Gene)
    if (!is.na(idx)) return(as.character(selection_report$Category[idx]))
  }
  if (g %in% target_discordant_genes) return("Target-only source-discordant")
  if (g %in% target_enhanced_genes) return("Target-enhanced transfer-retained")
  if (g %in% source_supported_genes) return("Source-supported transfer-stabilized")
  "Selected"
}
get_membership <- function(g, gene_set) {
  g %in% gene_set
}
gene_effect_rows <- lapply(selected_union, function(g) {
  vals <- boot_ide[, g]
  vals_finite <- vals[is.finite(vals)]
  data.frame(
    Gene = g,
    Raw_Gene = sub("^M_", "", g),
    Category = get_gene_category(g),
    In_TargetAnchored_Reported_Set = get_membership(g, target_anchored_genes),
    In_TransferSupported_Set = get_membership(g, transfer_supported_genes),
    In_Source_Supported_Subset = get_membership(g, source_supported_genes),
    In_Target_Enhanced_Subset = get_membership(g, target_enhanced_genes),
    In_TargetOnly_Stable_Signature = get_membership(g, targetonly_genes),
    In_TargetOnly_SourceDiscordant_Subset = get_membership(g, target_discordant_genes),
    Alpha_D_to_M = as.numeric(alpha_full[g]),
    Beta_M_to_Y = as.numeric(beta_full[g]),
    PIP_Full = as.numeric(pip_full[g]),
    IDE_Point = as.numeric(ide_full[g]),
    IDE_Boot_Mean = safe_mean(vals),
    IDE_Boot_SD = safe_sd(vals),
    IDE_95_Low = safe_quantile(vals, 0.025),
    IDE_95_Upp = safe_quantile(vals, 0.975),
    IDE_Crosses_Zero =
      safe_quantile(vals, 0.025) <= 0 &
      safe_quantile(vals, 0.975) >= 0,
    IDE_Sign_Stability = ifelse(
      length(vals_finite) == 0,
      NA_real_,
      max(mean(vals_finite > 0), mean(vals_finite < 0))
    ),
    N_Boot_Used = length(vals_finite),
    stringsAsFactors = FALSE
  )
})
gene_effects <- bind_rows(gene_effect_rows) %>%
  arrange(
    desc(In_TargetAnchored_Reported_Set),
    desc(In_TransferSupported_Set),
    desc(IDE_Sign_Stability),
    Raw_Gene
  )
# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------
run_meta <- data.frame(
  Item = c(
    "BOOT_EFFECT_B",
    "BOOT_SEED",
    "N_TargetAnchored_Reported_Set",
    "N_TransferSupported_Set",
    "N_Source_Supported_Subset",
    "N_Target_Enhanced_Subset",
    "N_TargetOnly_Stable_Signature",
    "N_TargetOnly_SourceDiscordant_Subset",
    "Bootstrap_Fit_Success",
    "Bootstrap_Has_Both_IDH_Classes",
    "Interpretation"
  ),
  Value = c(
    BOOT_EFFECT_B,
    BOOT_SEED,
    length(target_anchored_genes),
    length(transfer_supported_genes),
    length(source_supported_genes),
    length(target_enhanced_genes),
    length(targetonly_genes),
    length(target_discordant_genes),
    sum(boot_diag$Fit_Success),
    sum(boot_diag$Has_Both_IDH_Classes),
    "Conditional bootstrap for fixed target-anchored and transfer-supported mediator sets; not full post-selection bootstrap"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  boot_diag,
  file.path(OUT_DIR, "glioma_conditional_effect_bootstrap_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  run_meta,
  file.path(OUT_DIR, "glioma_conditional_effect_run_metadata.csv"),
  row.names = FALSE
)
write.csv(
  gene_effects,
  file.path(OUT_DIR, "glioma_gene_mediation_effects_conditional_bootstrap.csv"),
  row.names = FALSE
)
write.csv(
  overall_effects,
  file.path(OUT_DIR, "glioma_overall_effects_conditional_bootstrap.csv"),
  row.names = FALSE
)
cat("\nConditional effects complete.\n")
cat("Overall effects:\n")
print(overall_effects)
cat("\nGene-level effects:\n")
print(
  gene_effects[, c(
    "Raw_Gene",
    "Category",
    "IDE_Point",
    "IDE_95_Low",
    "IDE_95_Upp",
    "IDE_Crosses_Zero",
    "IDE_Sign_Stability",
    "In_TargetOnly_SourceDiscordant_Subset"
  )]
)
cat("\nOutput directory:\n", OUT_DIR, "\n", sep = "")