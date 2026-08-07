# ==============================================================================
# 09_run_glioma_discovery.R
# Run after file 07: Rscript 09_run_glioma_discovery.R
# Glioma VITA-HDM discovery analysis.
# Analysis specification:
#   - Uses the two-cohort case-study design:
#       Target: CGGA_mRNA_array_301
#       Source: TCGA_LGG
#       Exposure: IDH mutation status
#       Outcome: overall survival
#       Covariates: Age + Gender
#   - Uses the fixed audit-selected VITA-HDM hyperparameters.
#   - Changes the final reporting rule:
#       Source evidence is used for annotation and prioritization, not as a veto.
#   - Reports:
#       1. VITA-HDM transfer-supported set
#       2. Target-only source-discordant genes
#       3. Target-anchored final reported set
# Important:
#   - This is a discovery analysis.
#   - Bootstrap stability is target-resampling based.
#   - The final reported set is target-anchored:
#       source-supported stable genes
#       union target-enhanced transfer-retained genes
#       union target-only source-discordant genes.
# Outputs:
#   results/case_study/discovery/
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
ALIGNED_RDS <- file.path(
  PATHS$case_data, "glioma_example_aligned_object.rds"
)
OUT_DIR <- ensure_directory(PATHS$discovery)
if (!file.exists(ALIGNED_RDS)) {
  stop(
    "Cannot find aligned RDS. Run 07_prepare_glioma_data.R first: ",
    ALIGNED_RDS
  )
}
# ------------------------------------------------------------------------------
# Fixed main-analysis settings
# ------------------------------------------------------------------------------
TARGET_NAME <- "CGGA_mRNA_array_301"
SOURCE_NAME <- "TCGA_LGG"
COVARIATES <- c("Age", "Gender")
BINARY_COVARIATES <- c("Gender")
N_GENES_KEEP <- 3000
TARGET_BFDR <- 0.20
BASE_ETA <- 0.10
SOURCE_ETA <- 0.20
BASE_TAU <- 0.50
FINAL_STR <- 1.50
TRANSFER_BETA_STRENGTH <- 0.25
WEIGHT_MODE <- "single_source_lpml_reported"
LPML_TEMPERATURE <- 5
MAX_ITER <- 100
BOOT_B <- 200
BOOT_SEED <- 4100
EXTENDED_STABILITY_CUTOFF <- 0.50
CORE_STABILITY_CUTOFF <- 0.70
SOURCE_SUPPORT_CUTOFF <- 0.50
TARGET_ENHANCED_STABILITY_CUTOFF <- 0.50
# For target-anchored reporting:
# target-only stable genes are not discarded simply because source support is weak.
TARGET_ONLY_SOURCE_DISCORDANT_CUTOFF <- TARGET_ENHANCED_STABILITY_CUTOFF
# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
aligned_obj <- readRDS(ALIGNED_RDS)
target <- get_complete_cohort(
  aligned_obj,
  "target",
  covariates = COVARIATES
)
source_cohort <- get_complete_cohort(
  aligned_obj,
  "source",
  covariates = COVARIATES
)
cat("\n============================================================\n")
cat("Main glioma discovery analysis: target-anchored reporting\n")
cat("============================================================\n")
cat("Target:", target$cohort, "N =", nrow(target$clinical),
    "Events =", sum(target$clinical$OS_status == 1), "\n")
cat("Source:", source_cohort$cohort, "N =", nrow(source_cohort$clinical),
    "Events =", sum(source_cohort$clinical$OS_status == 1), "\n")
cat("Gene universe:", aligned_obj$meta$gene_universe_definition, "\n")
keep_genes <- select_top_variable_genes(
  target_cohort = target,
  n_genes_keep = N_GENES_KEEP,
  gene_universe = aligned_obj$common_genes
)
cat("Candidate genes retained by target variance:", length(keep_genes), "\n")
# ------------------------------------------------------------------------------
# Full-data target/source designs and source fit
# ------------------------------------------------------------------------------
target_full <- prepare_design(
  cohort = target,
  keep_genes = keep_genes,
  covariates = COVARIATES,
  binary_covariates = BINARY_COVARIATES,
  fit_scalers = TRUE
)
source_design <- prepare_source_design(
  source_cohort = source_cohort,
  keep_genes = keep_genes,
  covariates = COVARIATES,
  binary_covariates = BINARY_COVARIATES
)
source_fit <- fit_source_model(
  source_design = source_design,
  source_name = SOURCE_NAME,
  source_eta0 = SOURCE_ETA,
  tau_sq = BASE_TAU,
  max_iter = MAX_ITER,
  verbose = FALSE
)
# ------------------------------------------------------------------------------
# Full-data target fits: target-only and VITA-HDM
# ------------------------------------------------------------------------------
fit_target_full <- fit_target_with_transfer(
  target_design = target_full,
  transfer = FALSE,
  eta0 = BASE_ETA,
  tau_sq = BASE_TAU,
  max_iter = MAX_ITER,
  verbose = FALSE
)
fit_vita_full <- fit_target_with_transfer(
  target_design = target_full,
  source_design = source_design,
  source_name = SOURCE_NAME,
  source_fit = source_fit,
  transfer = TRUE,
  eta0 = BASE_ETA,
  source_eta0 = SOURCE_ETA,
  strength = FINAL_STR,
  tau_sq = BASE_TAU,
  transfer_beta_strength = TRANSFER_BETA_STRENGTH,
  weight_mode = WEIGHT_MODE,
  lpml_temperature = LPML_TEMPERATURE,
  max_iter = MAX_ITER,
  verbose = FALSE
)
sel_source_full <- get_selected_genes(
  source_fit$pip_gene,
  target_bfdr = TARGET_BFDR
)
sel_target_full <- get_selected_genes(
  fit_target_full$pip_gene,
  target_bfdr = TARGET_BFDR
)
sel_vita_full <- get_selected_genes(
  fit_vita_full$pip_gene,
  target_bfdr = TARGET_BFDR
)
cat("\nFull-data BFDR-selected genes:\n")
cat("  Source-only:", length(sel_source_full$selected), "\n")
cat("  Target-only:", length(sel_target_full$selected), "\n")
cat("  VITA-HDM:", length(sel_vita_full$selected), "\n")
# ------------------------------------------------------------------------------
# Bootstrap stability selection
# ------------------------------------------------------------------------------
all_gene_names <- target_full$GeneNames
boot_select_vita <- matrix(FALSE, nrow = BOOT_B, ncol = length(all_gene_names))
boot_select_target <- matrix(FALSE, nrow = BOOT_B, ncol = length(all_gene_names))
colnames(boot_select_vita) <- all_gene_names
colnames(boot_select_target) <- all_gene_names
boot_diag_rows <- list()
boot_weight_rows <- list()
cat("\nRunning bootstrap stability selection: B =", BOOT_B, "\n")
pb <- txtProgressBar(min = 0, max = BOOT_B, style = 3)
for (b in seq_len(BOOT_B)) {
  set.seed(BOOT_SEED + b)
  idx <- sample(seq_len(nrow(target$clinical)), replace = TRUE)
  target_b <- make_bootstrap_cohort(
    target,
    idx,
    prefix = paste0("boot", b)
  )
  d_tab <- table(target_b$clinical$IDH)
  d_has_both <- all(c("0", "1") %in% names(d_tab))
  fit_target_success <- FALSE
  fit_vita_success <- FALSE
  n_target_sel <- NA_integer_
  n_vita_sel <- NA_integer_
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
      fit_t <- try(
        fit_target_with_transfer(
          target_design = target_b_design,
          transfer = FALSE,
          eta0 = BASE_ETA,
          tau_sq = BASE_TAU,
          max_iter = MAX_ITER,
          verbose = FALSE
        ),
        silent = TRUE
      )
      if (!inherits(fit_t, "try-error")) {
        fit_target_success <- TRUE
        sel_t <- get_selected_genes(
          fit_t$pip_gene,
          target_bfdr = TARGET_BFDR
        )
        if (length(sel_t$selected) > 0) {
          boot_select_target[b, sel_t$selected] <- TRUE
        }
        n_target_sel <- length(sel_t$selected)
      } else {
        err_msg <- paste(
          err_msg,
          "target-fit failed:",
          as.character(attr(fit_t, "condition")$message)
        )
      }
      fit_v <- try(
        fit_target_with_transfer(
          target_design = target_b_design,
          source_design = source_design,
          source_name = SOURCE_NAME,
          source_fit = source_fit,
          transfer = TRUE,
          eta0 = BASE_ETA,
          source_eta0 = SOURCE_ETA,
          strength = FINAL_STR,
          tau_sq = BASE_TAU,
          transfer_beta_strength = TRANSFER_BETA_STRENGTH,
          weight_mode = WEIGHT_MODE,
          lpml_temperature = LPML_TEMPERATURE,
          max_iter = MAX_ITER,
          verbose = FALSE
        ),
        silent = TRUE
      )
      if (!inherits(fit_v, "try-error")) {
        fit_vita_success <- TRUE
        sel_v <- get_selected_genes(
          fit_v$pip_gene,
          target_bfdr = TARGET_BFDR
        )
        if (length(sel_v$selected) > 0) {
          boot_select_vita[b, sel_v$selected] <- TRUE
        }
        n_vita_sel <- length(sel_v$selected)
        if (!is.null(fit_v$priors$source_info) &&
            nrow(fit_v$priors$source_info) > 0) {
          boot_weight_rows[[length(boot_weight_rows) + 1]] <- data.frame(
            Bootstrap = b,
            fit_v$priors$source_info,
            stringsAsFactors = FALSE
          )
        }
      } else {
        err_msg <- paste(
          err_msg,
          "vita-fit failed:",
          as.character(attr(fit_v, "condition")$message)
        )
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
    Target_Fit_Success = fit_target_success,
    VITA_Fit_Success = fit_vita_success,
    Target_N_Selected = n_target_sel,
    VITA_N_Selected = n_vita_sel,
    Error_Message = err_msg,
    stringsAsFactors = FALSE
  )
  setTxtProgressBar(pb, b)
}
close(pb)
boot_diag <- bind_rows(boot_diag_rows)
boot_weights <- if (length(boot_weight_rows) > 0) {
  bind_rows(boot_weight_rows)
} else {
  data.frame()
}
valid_vita <- boot_diag$VITA_Fit_Success
valid_target <- boot_diag$Target_Fit_Success
if (sum(valid_vita) == 0) stop("All VITA bootstrap fits failed.")
if (sum(valid_target) == 0) warning("All target-only bootstrap fits failed.")
vita_stability <- colMeans(
  boot_select_vita[valid_vita, , drop = FALSE],
  na.rm = TRUE
)
target_stability <- colMeans(
  boot_select_target[valid_target, , drop = FALSE],
  na.rm = TRUE
)
# ------------------------------------------------------------------------------
# Selection report and target-anchored reporting rule
# ------------------------------------------------------------------------------
source_pip_bar <- fit_vita_full$priors$source_pip_bar[all_gene_names]
source_beta_bar <- fit_vita_full$priors$source_beta_bar[all_gene_names]
source_pip_bar[!is.finite(source_pip_bar)] <- NA_real_
source_beta_bar[!is.finite(source_beta_bar)] <- NA_real_
selection_report <- data.frame(
  Gene = all_gene_names,
  Raw_Gene = sub("^M_", "", all_gene_names),
  Final_Target_PIP_VITA = as.numeric(fit_vita_full$pip_gene[all_gene_names]),
  Final_Target_Beta_VITA = as.numeric(fit_vita_full$beta_gene[all_gene_names]),
  Final_Target_PIP_TargetOnly = as.numeric(fit_target_full$pip_gene[all_gene_names]),
  Final_Target_Beta_TargetOnly = as.numeric(fit_target_full$beta_gene[all_gene_names]),
  Source_PIP_Weighted = as.numeric(source_pip_bar),
  Source_Beta_Weighted = as.numeric(source_beta_bar),
  VITA_Stability = as.numeric(vita_stability[all_gene_names]),
  TargetOnly_Stability = as.numeric(target_stability[all_gene_names]),
  FullData_Selected_Source = all_gene_names %in% sel_source_full$selected,
  FullData_Selected_TargetOnly = all_gene_names %in% sel_target_full$selected,
  FullData_Selected_VITA = all_gene_names %in% sel_vita_full$selected,
  stringsAsFactors = FALSE
)
selection_report$Source_Supported <- with(
  selection_report,
  is.finite(Source_PIP_Weighted) & Source_PIP_Weighted >= SOURCE_SUPPORT_CUTOFF
)
selection_report$VITA_Extended <- with(
  selection_report,
  is.finite(VITA_Stability) & VITA_Stability >= EXTENDED_STABILITY_CUTOFF
)
selection_report$VITA_Core <- with(
  selection_report,
  is.finite(VITA_Stability) & VITA_Stability >= CORE_STABILITY_CUTOFF
)
selection_report$TargetOnly_Stable <- with(
  selection_report,
  is.finite(TargetOnly_Stability) &
    TargetOnly_Stability >= TARGET_ONLY_SOURCE_DISCORDANT_CUTOFF
)
selection_report$SourceSupported_Stable <- with(
  selection_report,
  VITA_Extended & Source_Supported
)
selection_report$Target_Enhanced <- with(
  selection_report,
  TargetOnly_Stable & VITA_Extended & !Source_Supported
)
selection_report$TargetOnly_SourceDiscordant <- with(
  selection_report,
  TargetOnly_Stable & !Source_Supported & !VITA_Extended
)
# Transfer-supported set keeps the previous VITA-HDM emphasis.
selection_report$TransferSupported_Signature <- with(
  selection_report,
  SourceSupported_Stable | Target_Enhanced
)
# Target-anchored reported set prevents source-domain evidence from acting as a veto.
selection_report$TargetAnchored_Reported <- with(
  selection_report,
  TransferSupported_Signature | TargetOnly_SourceDiscordant
)
# Backward-compatible column name.
# In this workflow, Manuscript_Signature means the target-anchored reported set.
selection_report$Manuscript_Signature <- selection_report$TargetAnchored_Reported
selection_report$Category <- "Other"
selection_report$Category[
  selection_report$Source_Supported & !selection_report$VITA_Extended
] <- "Filtered source-supported signal"
selection_report$Category[
  selection_report$TargetOnly_SourceDiscordant
] <- "Target-only source-discordant"
selection_report$Category[
  selection_report$Target_Enhanced
] <- "Target-enhanced transfer-retained"
selection_report$Category[
  selection_report$SourceSupported_Stable
] <- "Source-supported transfer-stabilized"
selection_report$Category[
  selection_report$VITA_Core & selection_report$Source_Supported
] <- "High-stability transfer core"
selection_report <- selection_report %>%
  arrange(
    desc(TargetAnchored_Reported),
    desc(TransferSupported_Signature),
    desc(VITA_Core),
    desc(VITA_Stability),
    desc(TargetOnly_Stability),
    desc(Source_PIP_Weighted),
    Raw_Gene
  )
selected_all_extended <- selection_report %>%
  filter(VITA_Extended)
selected_source_supported <- selection_report %>%
  filter(SourceSupported_Stable)
selected_core <- selection_report %>%
  filter(VITA_Core)
selected_transfer_supported <- selection_report %>%
  filter(TransferSupported_Signature)
selected_target_discordant <- selection_report %>%
  filter(TargetOnly_SourceDiscordant)
selected_target_anchored <- selection_report %>%
  filter(TargetAnchored_Reported)
targetonly_stable <- selection_report %>%
  filter(TargetOnly_Stable)
target_anchored_source_frac <- if (nrow(selected_target_anchored) > 0) {
  mean(selected_target_anchored$Source_Supported, na.rm = TRUE)
} else {
  NA_real_
}
transfer_supported_source_frac <- if (nrow(selected_transfer_supported) > 0) {
  mean(selected_transfer_supported$Source_Supported, na.rm = TRUE)
} else {
  NA_real_
}
# ------------------------------------------------------------------------------
# Summary and rule metadata
# ------------------------------------------------------------------------------
selection_summary <- data.frame(
  Parameter = c(
    "Target_N",
    "Target_Events",
    "Source_N",
    "Source_Events",
    "Candidate_Genes_Kept",
    "Source_FullData_BFDR_Selected",
    "Target_FullData_BFDR_Selected",
    "VITA_FullData_BFDR_Selected",
    "TargetOnly_Stable_Selected",
    "VITA_Extended_Stable_Selected",
    "VITA_SourceSupported_Stable_Selected",
    "VITA_Core_Stable_Selected",
    "TransferSupported_Signature_Selected",
    "TargetOnly_SourceDiscordant_Selected",
    "TargetAnchored_Reported_Selected",
    "TransferSupported_SourceSupported_Fraction",
    "TargetAnchored_SourceSupported_Fraction",
    "Bootstrap_B",
    "Bootstrap_Target_Success",
    "Bootstrap_VITA_Success"
  ),
  Value = c(
    nrow(target$clinical),
    sum(target$clinical$OS_status == 1),
    nrow(source_cohort$clinical),
    sum(source_cohort$clinical$OS_status == 1),
    length(keep_genes),
    length(sel_source_full$selected),
    length(sel_target_full$selected),
    length(sel_vita_full$selected),
    nrow(targetonly_stable),
    nrow(selected_all_extended),
    nrow(selected_source_supported),
    nrow(selected_core),
    nrow(selected_transfer_supported),
    nrow(selected_target_discordant),
    nrow(selected_target_anchored),
    transfer_supported_source_frac,
    target_anchored_source_frac,
    BOOT_B,
    sum(valid_target),
    sum(valid_vita)
  ),
  stringsAsFactors = FALSE
)
manuscript_rule <- data.frame(
  Rule_Item = c(
    "Gene universe",
    "Primary source",
    "N_GENES_KEEP",
    "TARGET_BFDR",
    "Extended stability cutoff",
    "Core stability cutoff",
    "Source support cutoff",
    "Target-only stability cutoff",
    "Source-supported stable definition",
    "Target-enhanced transfer-retained definition",
    "Target-only source-discordant definition",
    "Transfer-supported signature rule",
    "Target-anchored reported-set rule",
    "Interpretation"
  ),
  Value = c(
    "intersection(CGGA_mRNA_array_301, TCGA_LGG)",
    SOURCE_NAME,
    as.character(N_GENES_KEEP),
    as.character(TARGET_BFDR),
    as.character(EXTENDED_STABILITY_CUTOFF),
    as.character(CORE_STABILITY_CUTOFF),
    paste0("weighted source PIP >= ", SOURCE_SUPPORT_CUTOFF),
    paste0("target-only bootstrap stability >= ", TARGET_ONLY_SOURCE_DISCORDANT_CUTOFF),
    paste0("VITA stability >= ", EXTENDED_STABILITY_CUTOFF,
           " and source-supported"),
    paste0("target-only stability >= ", TARGET_ONLY_SOURCE_DISCORDANT_CUTOFF,
           " and VITA stability >= ", EXTENDED_STABILITY_CUTOFF,
           " and not source-supported"),
    paste0("target-only stability >= ", TARGET_ONLY_SOURCE_DISCORDANT_CUTOFF,
           " and not source-supported and VITA stability < ", EXTENDED_STABILITY_CUTOFF),
    "source-supported stable genes union target-enhanced transfer-retained genes",
    "transfer-supported signature union target-only source-discordant genes",
    "target-domain stability is retained and annotated; source evidence prioritizes transferability but is not used as a veto"
  ),
  stringsAsFactors = FALSE
)
source_weights_full <- fit_vita_full$priors$source_info
# ------------------------------------------------------------------------------
# Save outputs
# ------------------------------------------------------------------------------
write.csv(
  selection_report,
  file.path(OUT_DIR, "glioma_bootstrap_selection_report.csv"),
  row.names = FALSE
)
write.csv(
  selected_all_extended,
  file.path(OUT_DIR, "glioma_final_selected_genes_all_extended.csv"),
  row.names = FALSE
)
write.csv(
  selected_source_supported,
  file.path(OUT_DIR, "glioma_final_selected_genes_source_supported.csv"),
  row.names = FALSE
)
write.csv(
  selected_core,
  file.path(OUT_DIR, "glioma_final_selected_genes_core.csv"),
  row.names = FALSE
)
write.csv(
  selected_transfer_supported,
  file.path(OUT_DIR, "glioma_final_selected_genes_transfer_supported.csv"),
  row.names = FALSE
)
write.csv(
  selected_target_discordant,
  file.path(OUT_DIR, "glioma_targetonly_source_discordant_genes.csv"),
  row.names = FALSE
)
write.csv(
  selected_target_anchored,
  file.path(OUT_DIR, "glioma_final_selected_genes_target_anchored.csv"),
  row.names = FALSE
)
# Backward-compatible manuscript filename.
# This alias contains the target-anchored reported set.
write.csv(
  selected_target_anchored,
  file.path(OUT_DIR, "glioma_final_selected_genes_manuscript.csv"),
  row.names = FALSE
)
write.csv(
  targetonly_stable,
  file.path(OUT_DIR, "glioma_targetonly_stable_genes.csv"),
  row.names = FALSE
)
write.csv(
  selection_summary,
  file.path(OUT_DIR, "glioma_selection_summary.csv"),
  row.names = FALSE
)
write.csv(
  manuscript_rule,
  file.path(OUT_DIR, "glioma_manuscript_signature_rule.csv"),
  row.names = FALSE
)
write.csv(
  source_weights_full,
  file.path(OUT_DIR, "glioma_source_weights_full.csv"),
  row.names = FALSE
)
write.csv(
  boot_diag,
  file.path(OUT_DIR, "glioma_bootstrap_diagnostics.csv"),
  row.names = FALSE
)
write.csv(
  boot_weights,
  file.path(OUT_DIR, "glioma_bootstrap_source_weights.csv"),
  row.names = FALSE
)
selection_objects <- list(
  target = target,
  source = source_cohort,
  keep_genes = keep_genes,
  target_full = target_full,
  source_design = source_design,
  source_fit = source_fit,
  fit_target_full = fit_target_full,
  fit_vita_full = fit_vita_full,
  selection_report = selection_report,
  selected_manuscript = selected_target_anchored,
  selected_target_anchored = selected_target_anchored,
  selected_transfer_supported = selected_transfer_supported,
  selected_all_extended = selected_all_extended,
  selected_source_supported = selected_source_supported,
  selected_core = selected_core,
  selected_target_discordant = selected_target_discordant,
  targetonly_stable = targetonly_stable,
  boot_diag = boot_diag,
  boot_weights = boot_weights,
  parameters = list(
    TARGET_NAME = TARGET_NAME,
    SOURCE_NAME = SOURCE_NAME,
    COVARIATES = COVARIATES,
    BINARY_COVARIATES = BINARY_COVARIATES,
    N_GENES_KEEP = N_GENES_KEEP,
    TARGET_BFDR = TARGET_BFDR,
    BASE_ETA = BASE_ETA,
    SOURCE_ETA = SOURCE_ETA,
    BASE_TAU = BASE_TAU,
    FINAL_STR = FINAL_STR,
    TRANSFER_BETA_STRENGTH = TRANSFER_BETA_STRENGTH,
    WEIGHT_MODE = WEIGHT_MODE,
    LPML_TEMPERATURE = LPML_TEMPERATURE,
    MAX_ITER = MAX_ITER,
    BOOT_B = BOOT_B,
    BOOT_SEED = BOOT_SEED,
    EXTENDED_STABILITY_CUTOFF = EXTENDED_STABILITY_CUTOFF,
    CORE_STABILITY_CUTOFF = CORE_STABILITY_CUTOFF,
    SOURCE_SUPPORT_CUTOFF = SOURCE_SUPPORT_CUTOFF,
    TARGET_ENHANCED_STABILITY_CUTOFF = TARGET_ENHANCED_STABILITY_CUTOFF,
    TARGET_ONLY_SOURCE_DISCORDANT_CUTOFF = TARGET_ONLY_SOURCE_DISCORDANT_CUTOFF,
    FINAL_REPORTING_RULE = "target-anchored; source evidence annotates transferability but does not veto target-only stable signals"
  )
)
saveRDS(
  selection_objects,
  file.path(OUT_DIR, "glioma_final_selection_objects.rds")
)
cat("\nMain discovery complete.\n")
cat("Target-anchored final reported genes:\n")
if (nrow(selected_target_anchored) == 0) {
  cat("  None selected under the target-anchored rule.\n")
} else {
  print(
    selected_target_anchored[, c(
      "Raw_Gene",
      "Category",
      "VITA_Stability",
      "TargetOnly_Stability",
      "Source_PIP_Weighted",
      "Final_Target_PIP_VITA",
      "TransferSupported_Signature",
      "TargetOnly_SourceDiscordant"
    )]
  )
}
cat("\nTransfer-supported subset:\n")
if (nrow(selected_transfer_supported) == 0) {
  cat("  None.\n")
} else {
  print(
    selected_transfer_supported[, c(
      "Raw_Gene",
      "Category",
      "VITA_Stability",
      "TargetOnly_Stability",
      "Source_PIP_Weighted"
    )]
  )
}
cat("\nTarget-only source-discordant genes:\n")
if (nrow(selected_target_discordant) == 0) {
  cat("  None.\n")
} else {
  print(
    selected_target_discordant[, c(
      "Raw_Gene",
      "Category",
      "VITA_Stability",
      "TargetOnly_Stability",
      "Source_PIP_Weighted"
    )]
  )
}
cat("\nOutput directory:\n", OUT_DIR, "\n", sep = "")