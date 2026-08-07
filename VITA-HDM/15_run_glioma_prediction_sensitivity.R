# ==============================================================================
# 15_run_glioma_prediction_sensitivity.R
# Run after file 09: Rscript 15_run_glioma_prediction_sensitivity.R
# Target-anchored glioma diagnostics:
# SA4: Fixed-signature repeated-CV prediction diagnostics
# SA5: Cox-Snell residual diagnostics for log-normal AFT assumption
# SA6: Source-target transferability diagnostics
# Uses 13_glioma_sensitivity_common.R. No sink() logging is used.
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
PROJECT_ROOT <- find_project_root(SCRIPT_DIR)
PATHS <- project_paths(PROJECT_ROOT)
source(file.path(PROJECT_ROOT, "13_glioma_sensitivity_common.R"))
N_REPEATS <- 10
N_FOLDS <- 5
# ==============================================================================
# SA4: Fixed-signature repeated-CV prediction diagnostics
# ==============================================================================
message("\n================================================================")
message("SA4: Fixed-signature repeated-CV prediction diagnostics")
message("================================================================")
signature_sets <- list(
  VITA_HDM_TargetAnchored_Reported_Set = target_anchored_signature
)
if (length(transfer_supported_signature) > 0) {
  signature_sets$VITA_HDM_TransferSupported_Set <- transfer_supported_signature
}
if (length(source_supported_signature) > 0) {
  signature_sets$VITA_HDM_SourceSupported_Subset <- source_supported_signature
}
if (length(core_signature) > 0) {
  signature_sets$VITA_HDM_Core <- core_signature
}
if (length(target_only_signature) > 0) {
  signature_sets$TargetOnly_Stable <- target_only_signature
}
if (length(target_discordant_signature) > 0) {
  signature_sets$TargetOnly_SourceDiscordant <- target_discordant_signature
}
signature_sets <- signature_sets[lengths(signature_sets) > 0]
if (length(signature_sets) == 0) {
  stop("No non-empty signature sets are available for SA4.")
}
pred <- evaluate_multiple_fixed_signatures_cox(
  target_cohort = target,
  signature_list = signature_sets,
  covariates = COVARIATES,
  binary_covariates = BINARY_COVARIATES,
  repeats = N_REPEATS,
  folds = N_FOLDS,
  seed = 2030
)
prediction_all <- pred$repeated_cv
prediction_summary <- pred$summary
if (is.null(prediction_all) || nrow(prediction_all) == 0) {
  stop("SA4 prediction_all is empty. Check evaluate_multiple_fixed_signatures_cox().")
}
if (is.null(prediction_summary) || nrow(prediction_summary) == 0) {
  stop("SA4 prediction_summary is empty. Check evaluate_multiple_fixed_signatures_cox().")
}
required_pred_cols <- c("Repeat", "Fold", "Model", "C_Index")
missing_pred_cols <- setdiff(required_pred_cols, colnames(prediction_all))
if (length(missing_pred_cols) > 0) {
  stop(
    "prediction_all is missing required columns: ",
    paste(missing_pred_cols, collapse = ", ")
  )
}
wide_cv <- prediction_all %>%
  dplyr::select(Repeat, Fold, Model, C_Index) %>%
  tidyr::pivot_wider(
    names_from = Model,
    values_from = C_Index
  )
make_contrast <- function(df, model_a, model_b, label) {
  if (!(model_a %in% colnames(df)) || !(model_b %in% colnames(df))) {
    return(NULL)
  }
  diff_vec <- df[[model_a]] - df[[model_b]]
  diff_vec <- diff_vec[is.finite(diff_vec)]
  if (length(diff_vec) == 0) {
    return(NULL)
  }
  data.frame(
    Contrast = label,
    Model_A = model_a,
    Model_B = model_b,
    N = length(diff_vec),
    Mean_Difference = mean(diff_vec),
    SD_Difference = ifelse(length(diff_vec) > 1, stats::sd(diff_vec), NA_real_),
    Positive_N = sum(diff_vec > 0),
    Positive_Rate = mean(diff_vec > 0),
    stringsAsFactors = FALSE
  )
}
contrast_list <- list(
  anchored_vs_transfer = make_contrast(
    wide_cv,
    "VITA_HDM_TargetAnchored_Reported_Set",
    "VITA_HDM_TransferSupported_Set",
    "Target-anchored reported set minus transfer-supported set"
  ),
  anchored_vs_targetonly = make_contrast(
    wide_cv,
    "VITA_HDM_TargetAnchored_Reported_Set",
    "TargetOnly_Stable",
    "Target-anchored reported set minus target-only stable signature"
  ),
  anchored_vs_clinical = make_contrast(
    wide_cv,
    "VITA_HDM_TargetAnchored_Reported_Set",
    "Clinical-Only",
    "Target-anchored reported set minus clinical-only"
  ),
  transfer_vs_clinical = make_contrast(
    wide_cv,
    "VITA_HDM_TransferSupported_Set",
    "Clinical-Only",
    "Transfer-supported set minus clinical-only"
  ),
  discordant_vs_clinical = make_contrast(
    wide_cv,
    "TargetOnly_SourceDiscordant",
    "Clinical-Only",
    "Target-only source-discordant subset minus clinical-only"
  )
)
contrast_list <- contrast_list[!vapply(contrast_list, is.null, logical(1))]
prediction_contrasts <- if (length(contrast_list) > 0) {
  dplyr::bind_rows(contrast_list)
} else {
  data.frame(
    Contrast = character(0),
    Model_A = character(0),
    Model_B = character(0),
    N = integer(0),
    Mean_Difference = numeric(0),
    SD_Difference = numeric(0),
    Positive_N = integer(0),
    Positive_Rate = numeric(0),
    stringsAsFactors = FALSE
  )
}
write_csv_out(prediction_all, "SA4_FixedSignature_RepeatedCV_All.csv")
write_csv_out(prediction_summary, "SA4_FixedSignature_Prediction_Summary.csv")
write_csv_out(prediction_contrasts, "SA4_FixedSignature_Pairwise_Contrasts.csv")
print(prediction_summary)
if (nrow(prediction_contrasts) > 0) {
  print(prediction_contrasts)
} else {
  message(">>> SA4 pairwise contrast table is empty because one or more compared models were unavailable.")
}
plot_sa4_df <- prediction_summary %>%
  dplyr::filter(Model != "Clinical-Only") %>%
  dplyr::mutate(
    Model = factor(Model, levels = Model[order(RepeatedCV_Mean_C_Index)])
  )
if (nrow(plot_sa4_df) > 0) {
  p_sa4 <- ggplot2::ggplot(
    plot_sa4_df,
    ggplot2::aes(x = Model, y = RepeatedCV_Mean_C_Index)
  ) +
    ggplot2::geom_col(width = 0.65, alpha = 0.85) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = RepeatedCV_Mean_C_Index - RepeatedCV_SD_C_Index,
        ymax = RepeatedCV_Mean_C_Index + RepeatedCV_SD_C_Index
      ),
      width = 0.2
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Fixed-signature prediction",
      x = NULL,
      y = "Repeated-CV C-index"
    ) +
    ggplot2::theme_classic(base_size = 12)
  save_plot_out(
    p_sa4,
    "SA4_FixedSignature_Prediction",
    width = 8,
    height = 5
  )
}
# ==============================================================================
# SA5: Cox-Snell residual diagnostics
# ==============================================================================
message("\n================================================================")
message("SA5: Cox-Snell residual diagnostics")
message("================================================================")
run_coxsnell_for_signature <- function(signature_name, signature_genes) {
  signature_genes <- intersect(add_M(signature_genes), paste0("M_", keep_genes))
  if (length(signature_genes) == 0) {
    warning("Skipping ", signature_name, ": empty signature.")
    return(NULL)
  }
  dat_obj <- try(
    build_signature_analysis_data(
      target_cohort = target,
      signature_genes = signature_genes,
      covariates = COVARIATES,
      binary_covariates = BINARY_COVARIATES
    ),
    silent = TRUE
  )
  if (inherits(dat_obj, "try-error")) {
    warning("Skipping ", signature_name, ": failed to build analysis data.")
    return(NULL)
  }
  analysis_data <- dat_obj$data
  clin_vars <- dat_obj$clin_vars
  sig_vars <- dat_obj$signature_genes
  all_vars <- c(clin_vars, sig_vars)
  all_vars <- intersect(all_vars, colnames(analysis_data))
  if (length(all_vars) == 0) {
    warning("Skipping ", signature_name, ": no usable variables.")
    return(NULL)
  }
  form_aft <- stats::as.formula(
    paste(
      "Surv(Time, Status) ~",
      paste(bt(all_vars), collapse = " + ")
    )
  )
  fit_aft <- try(
    survival::survreg(
      form_aft,
      data = analysis_data,
      dist = "lognormal"
    ),
    silent = TRUE
  )
  if (inherits(fit_aft, "try-error")) {
    warning("Skipping ", signature_name, ": survreg failed.")
    return(NULL)
  }
  mu_hat <- try(
    stats::predict(fit_aft, type = "lp"),
    silent = TRUE
  )
  if (inherits(mu_hat, "try-error") || any(!is.finite(mu_hat))) {
    warning("Skipping ", signature_name, ": invalid AFT linear predictor.")
    return(NULL)
  }
  sigma_hat <- fit_aft$scale
  if (!is.finite(sigma_hat) || sigma_hat <= 0) {
    warning("Skipping ", signature_name, ": invalid AFT scale.")
    return(NULL)
  }
  surv_prob <- 1 - stats::pnorm((log(analysis_data$Time) - mu_hat) / sigma_hat)
  surv_prob <- pmax(surv_prob, 1e-12)
  cs_residuals <- -log(surv_prob)
  fit_cs <- try(
    survival::survfit(
      survival::Surv(cs_residuals, analysis_data$Status) ~ 1
    ),
    silent = TRUE
  )
  if (inherits(fit_cs, "try-error")) {
    warning("Skipping ", signature_name, ": Cox-Snell survfit failed.")
    return(NULL)
  }
  cs_plot_data <- data.frame(
    Signature_Set = signature_name,
    Residual_Time = fit_cs$time,
    Cumulative_Hazard = -log(fit_cs$surv),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(
      is.finite(Residual_Time),
      is.finite(Cumulative_Hazard)
    )
  if (nrow(cs_plot_data) == 0) {
    warning("Skipping ", signature_name, ": empty Cox-Snell plot data.")
    return(NULL)
  }
  cs_summary <- data.frame(
    Signature_Set = signature_name,
    N_Genes = length(sig_vars),
    N_Points = nrow(cs_plot_data),
    RMSE_from_UnitSlope = sqrt(
      mean(
        (cs_plot_data$Residual_Time - cs_plot_data$Cumulative_Hazard)^2,
        na.rm = TRUE
      )
    ),
    MaxAbs_Deviation = max(
      abs(cs_plot_data$Residual_Time - cs_plot_data$Cumulative_Hazard),
      na.rm = TRUE
    ),
    Tail_Residual_Max = max(cs_plot_data$Residual_Time, na.rm = TRUE),
    Tail_CumHazard_Max = max(cs_plot_data$Cumulative_Hazard, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  list(
    data = cs_plot_data,
    summary = cs_summary
  )
}
cs_sets <- list(
  TargetAnchored_Reported_Set = target_anchored_signature
)
if (length(transfer_supported_signature) > 0) {
  cs_sets$TransferSupported_Set <- transfer_supported_signature
}
if (length(target_only_signature) > 0) {
  cs_sets$TargetOnly_Stable <- target_only_signature
}
cs_sets <- cs_sets[lengths(cs_sets) > 0]
cs_results <- lapply(
  names(cs_sets),
  function(nm) run_coxsnell_for_signature(nm, cs_sets[[nm]])
)
cs_results <- cs_results[!vapply(cs_results, is.null, logical(1))]
if (length(cs_results) == 0) {
  warning("No Cox-Snell diagnostics could be computed.")
  cs_plot_all <- data.frame()
  cs_summary_all <- data.frame()
} else {
  cs_plot_all <- dplyr::bind_rows(lapply(cs_results, `[[`, "data"))
  cs_summary_all <- dplyr::bind_rows(lapply(cs_results, `[[`, "summary"))
}
# Publication labels for Cox-Snell facets.  The raw Signature_Set column is
# preserved in the output CSV for reproducibility.
cs_signature_label_map <- c(
  "TargetAnchored_Reported_Set" = "Target-anchored reported set",
  "TargetOnly_Stable" = "Target-only stable signature",
  "TransferSupported_Set" = "Transfer-supported set"
)
if (nrow(cs_plot_all) > 0 && "Signature_Set" %in% colnames(cs_plot_all)) {
  cs_plot_all <- cs_plot_all %>%
    dplyr::mutate(
      Signature_Set_Label = dplyr::recode(
        Signature_Set,
        !!!cs_signature_label_map,
        .default = Signature_Set
      )
    )
}
if (nrow(cs_summary_all) > 0 && "Signature_Set" %in% colnames(cs_summary_all)) {
  cs_summary_all <- cs_summary_all %>%
    dplyr::mutate(
      Signature_Set_Label = dplyr::recode(
        Signature_Set,
        !!!cs_signature_label_map,
        .default = Signature_Set
      )
    )
}
write_csv_out(cs_plot_all, "SA5_CoxSnell_Residual_Data.csv")
write_csv_out(cs_summary_all, "SA5_CoxSnell_Residual_Summary.csv")
if (nrow(cs_plot_all) > 0) {
  max_axis <- max(
    c(cs_plot_all$Residual_Time, cs_plot_all$Cumulative_Hazard),
    na.rm = TRUE
  )
  p_cs <- ggplot2::ggplot(
    cs_plot_all,
    ggplot2::aes(x = Residual_Time, y = Cumulative_Hazard)
  ) +
    ggplot2::geom_step(linewidth = 1.0) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      linewidth = 0.8
    ) +
    ggplot2::facet_wrap(~ Signature_Set_Label) +
    ggplot2::labs(
      title = "Cox-Snell residual diagnostics",
      x = "Cox-Snell residuals",
      y = "Estimated cumulative hazard"
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::coord_fixed(
      ratio = 1,
      xlim = c(0, max_axis),
      ylim = c(0, max_axis)
    )
  save_plot_out(
    p_cs,
    "SA5_CoxSnell_Residual_Diagnostic",
    width = 9,
    height = 5
  )
}
print(cs_summary_all)
# ============================================================================== 
# SA6: Source-target transferability diagnostics
# ============================================================================== 
message("\n================================================================")
message("SA6: Source-target transferability diagnostics")
message("================================================================")
# ------------------------------------------------------------------------------
# SA6A: Source LPML / source weight
# ------------------------------------------------------------------------------
source_info <- if (!is.null(obj$fit_vita_full$priors$source_info)) {
  obj$fit_vita_full$priors$source_info
} else if (!is.null(obj$source_weights_full)) {
  obj$source_weights_full
} else {
  data.frame()
}
write_csv_out(source_info, "SA6_SourceWeights_FullData.csv")
# ------------------------------------------------------------------------------
# SA6B: PCA source-target feature distribution
# ------------------------------------------------------------------------------
# Two PCA diagnostics are reported:
#   1. Raw pooled feature-space PCA: evaluates cross-platform separation in the
#      observed candidate-gene expression space.
#   2. Within-cohort z-score PCA: target and source are standardized separately
#      gene-by-gene before pooling. This removes first-order cohort/platform mean
#      and scale differences and better reflects the fact that VITA-HDM transfers
#      source information through model-derived priors rather than direct pooling
#      of raw expression values.
compute_source_target_pca_local <- function(raw_genes,
                                            scale_mode = c("raw_pooled_zscore", "within_cohort_zscore"),
                                            max_genes = NULL) {
  scale_mode <- match.arg(scale_mode)
  raw_genes <- unique(as.character(raw_genes))
  raw_genes <- intersect(raw_genes, rownames(target$expr))
  raw_genes <- intersect(raw_genes, rownames(source_cohort$expr))
  if (!is.null(max_genes) && length(raw_genes) > max_genes) {
    raw_genes <- raw_genes[seq_len(max_genes)]
  }
  if (length(raw_genes) < 2) {
    stop("PCA requires at least two shared genes.")
  }
  target_x <- t(target$expr[raw_genes, target$clinical$Sample_ID, drop = FALSE])
  source_x <- t(source_cohort$expr[raw_genes, source_cohort$clinical$Sample_ID, drop = FALSE])
  colnames(target_x) <- raw_genes
  colnames(source_x) <- raw_genes
  target_x <- as.matrix(target_x)
  source_x <- as.matrix(source_x)
  mode(target_x) <- "numeric"
  mode(source_x) <- "numeric"
  target_x[!is.finite(target_x)] <- NA_real_
  source_x[!is.finite(source_x)] <- NA_real_
  if (scale_mode == "within_cohort_zscore") {
    target_sc <- safe_scale_train(target_x)$X
    source_sc <- safe_scale_train(source_x)$X
    x_all <- rbind(target_sc, source_sc)
    pca_label <- "Within-cohort z-score PCA"
  } else {
    x_all <- rbind(target_x, source_x)
    x_all <- safe_scale_train(x_all)$X
    pca_label <- "Raw pooled feature-space PCA"
  }
  x_all[!is.finite(x_all)] <- 0
  dataset <- c(
    rep(TARGET_NAME, nrow(target_x)),
    rep(SOURCE_NAME, nrow(source_x))
  )
  pca <- stats::prcomp(x_all, center = FALSE, scale. = FALSE)
  var_exp <- 100 * (pca$sdev^2) / sum(pca$sdev^2)
  data.frame(
    Sample_ID = c(rownames(target_x), rownames(source_x)),
    Dataset = dataset,
    PCA_Mode = pca_label,
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    PC1_Var = round(var_exp[1], 2),
    PC2_Var = round(var_exp[2], 2),
    N_Genes = length(raw_genes),
    stringsAsFactors = FALSE
  )
}
plot_pca_local <- function(pca_df, title, subtitle, filename_base) {
  if (is.null(pca_df) || nrow(pca_df) == 0) return(invisible(NULL))
  var_exp <- c(
    unique(pca_df$PC1_Var)[1],
    unique(pca_df$PC2_Var)[1]
  )
  p_pca <- ggplot2::ggplot(
    pca_df,
    ggplot2::aes(x = PC1, y = PC2, color = Dataset)
  ) +
    ggplot2::geom_point(alpha = 0.45, size = 1.4) +
    ggplot2::stat_ellipse(
      level = 0.95,
      linetype = "dashed",
      linewidth = 0.8,
      na.rm = TRUE
    ) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = paste0("PC1 (", var_exp[1], "% variance)"),
      y = paste0("PC2 (", var_exp[2], "% variance)")
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_blank()
    )
  save_plot_out(
    p_pca,
    filename_base,
    width = 7,
    height = 5
  )
}
pca_raw <- try(
  compute_source_target_pca_local(
    raw_genes = keep_genes,
    scale_mode = "raw_pooled_zscore"
  ),
  silent = TRUE
)
if (inherits(pca_raw, "try-error")) {
  warning("SA6 raw feature-space PCA failed: ", as.character(attr(pca_raw, "condition")$message))
  pca_raw <- data.frame()
}
pca_within <- try(
  compute_source_target_pca_local(
    raw_genes = keep_genes,
    scale_mode = "within_cohort_zscore"
  ),
  silent = TRUE
)
if (inherits(pca_within, "try-error")) {
  warning("SA6 within-cohort z-score PCA failed: ", as.character(attr(pca_within, "condition")$message))
  pca_within <- data.frame()
}
# Participant-level PCA coordinates retain Sample_ID and are intentionally
# isolated from publication tables. The figures and aggregate diagnostics remain
# available under the standard sensitivity output directories.
write_csv_out(
  pca_raw,
  "SA6_PCA_SourceTarget_Data.csv",
  subdir = SENS_PRIVATE_DIR
)
write_csv_out(
  pca_raw,
  "SA6_PCA_SourceTarget_RawPooled_Data.csv",
  subdir = SENS_PRIVATE_DIR
)
write_csv_out(
  pca_within,
  "SA6_PCA_SourceTarget_WithinCohortZ_Data.csv",
  subdir = SENS_PRIVATE_DIR
)
if (nrow(pca_raw) > 0) {
  plot_pca_local(
    pca_raw,
    title = "Raw source-target feature-space separation",
    subtitle = paste0("Raw pooled candidate-gene space: ", length(keep_genes), " genes"),
    filename_base = "SA6_PCA_SourceTarget_RawFeatureSpace"
  )
  # Backward-compatible output name; content is the raw feature-space diagnostic.
  plot_pca_local(
    pca_raw,
    title = "Raw source-target feature-space separation",
    subtitle = paste0("Raw pooled candidate-gene space: ", length(keep_genes), " genes"),
    filename_base = "SA6_PCA_SourceTarget_Overlap"
  )
}
if (nrow(pca_within) > 0) {
  plot_pca_local(
    pca_within,
    title = "Within-cohort standardized source-target PCA",
    subtitle = "Each cohort was gene-wise z-scored before pooled PCA",
    filename_base = "SA6_PCA_SourceTarget_WithinCohortZ"
  )
}
# ------------------------------------------------------------------------------
# SA6C: Source PIP vs target-domain stability
# ------------------------------------------------------------------------------
sa6_stab <- baseline_classification %>%
  dplyr::mutate(
    In_TargetAnchored_Reference = Gene %in% target_anchored_signature,
    In_TransferSupported_Reference = Gene %in% transfer_supported_signature,
    In_SourceSupported_Reference = Gene %in% source_supported_signature,
    In_Core_Reference = Gene %in% core_signature,
    In_TargetOnly_Reference = Gene %in% target_only_signature,
    In_TargetDiscordant_Reference = Gene %in% target_discordant_signature,
    In_TargetEnhanced_Reference = Gene %in% target_enhanced_signature,
    SourceTarget_Quadrant = dplyr::case_when(
      Source_Supported & TargetOnly_Stable ~ "High source PIP + high target stability",
      Source_Supported & !TargetOnly_Stable ~ "High source PIP + low target stability",
      !Source_Supported & TargetOnly_Stable ~ "Low source PIP + high target stability",
      TRUE ~ "Low source PIP + low target stability"
    ),
    Plot_Class = dplyr::case_when(
      In_Core_Reference ~ "High-stability transfer core",
      SourceSupported_Stable ~ "Source-supported transfer-stabilized",
      Target_Enhanced ~ "Target-enhanced transfer-retained",
      TargetOnly_SourceDiscordant ~ "Target-only source-discordant",
      Source_Supported ~ "Source-supported background",
      TRUE ~ "Background"
    )
  )
safe_cor_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(
    stats::cor(x[ok], y[ok], method = "spearman")
  )
}
sa6_corr <- data.frame(
  Metric = c(
    "Spearman_SourcePIP_vs_VITAStability",
    "Spearman_SourcePIP_vs_TargetOnlyStability",
    "Spearman_VITAStability_vs_TargetOnlyStability"
  ),
  Estimate = c(
    safe_cor_spearman(sa6_stab$Source_PIP_Weighted, sa6_stab$VITA_Stability),
    safe_cor_spearman(sa6_stab$Source_PIP_Weighted, sa6_stab$TargetOnly_Stability),
    safe_cor_spearman(sa6_stab$VITA_Stability, sa6_stab$TargetOnly_Stability)
  ),
  stringsAsFactors = FALSE
)
quadrant_levels <- c(
  "High source PIP + high target stability",
  "High source PIP + low target stability",
  "Low source PIP + high target stability",
  "Low source PIP + low target stability"
)
sa6_quadrant <- sa6_stab %>%
  dplyr::mutate(
    SourceTarget_Quadrant = factor(SourceTarget_Quadrant, levels = quadrant_levels)
  ) %>%
  dplyr::group_by(SourceTarget_Quadrant) %>%
  dplyr::summarise(
    N_Genes = dplyr::n(),
    N_TargetAnchored = sum(TargetAnchored_Reported, na.rm = TRUE),
    N_TransferSupported = sum(TransferSupported_Signature, na.rm = TRUE),
    N_TargetDiscordant = sum(TargetOnly_SourceDiscordant, na.rm = TRUE),
    Genes_TargetAnchored = paste(Gene[TargetAnchored_Reported], collapse = " | "),
    .groups = "drop"
  ) %>%
  tidyr::complete(
    SourceTarget_Quadrant = factor(quadrant_levels, levels = quadrant_levels),
    fill = list(
      N_Genes = 0L,
      N_TargetAnchored = 0L,
      N_TransferSupported = 0L,
      N_TargetDiscordant = 0L,
      Genes_TargetAnchored = ""
    )
  ) %>%
  dplyr::arrange(SourceTarget_Quadrant)
write_csv_out(sa6_stab, "SA6_SourcePIP_TargetStability_Data.csv")
write_csv_out(sa6_corr, "SA6_SourcePIP_TargetStability_Correlation.csv")
write_csv_out(sa6_quadrant, "SA6_SourceTarget_Quadrant_Summary.csv")
plot_class_short <- c(
  "High-stability transfer core" = "Core",
  "Source-supported transfer-stabilized" = "Source-supported",
  "Target-enhanced transfer-retained" = "Target-enhanced",
  "Target-only source-discordant" = "Source-discordant",
  "Source-supported background" = "Source-supported background",
  "Background" = "Background"
)
shape_short <- c(
  "Core" = 16,
  "Source-supported" = 17,
  "Target-enhanced" = 15,
  "Source-discordant" = 3,
  "Source-supported background" = 16,
  "Background" = 16
)
plot_df <- sa6_stab %>%
  dplyr::mutate(
    TargetOnly_Stability = suppressWarnings(as.numeric(TargetOnly_Stability)),
    Source_PIP_Weighted = suppressWarnings(as.numeric(Source_PIP_Weighted)),
    Plot_Class_Short = dplyr::recode(
      Plot_Class,
      !!!plot_class_short,
      .default = Plot_Class
    )
  ) %>%
  dplyr::filter(
    is.finite(TargetOnly_Stability),
    is.finite(Source_PIP_Weighted)
  )
background_df <- plot_df %>% dplyr::filter(!TargetAnchored_Reported)
label_df <- plot_df %>% dplyr::filter(TargetAnchored_Reported)
p_pip <- ggplot2::ggplot() +
  ggplot2::geom_hline(
    yintercept = SOURCE_SUPPORT_CUTOFF,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  ggplot2::geom_vline(
    xintercept = TARGET_ONLY_SOURCE_DISCORDANT_CUTOFF,
    linetype = "dashed",
    linewidth = 0.8
  ) +
  ggplot2::geom_point(
    data = background_df,
    ggplot2::aes(
      x = TargetOnly_Stability,
      y = Source_PIP_Weighted
    ),
    alpha = 0.18,
    size = 1.4
  ) +
  ggplot2::geom_point(
    data = label_df,
    ggplot2::aes(
      x = TargetOnly_Stability,
      y = Source_PIP_Weighted,
      shape = Plot_Class_Short
    ),
    alpha = 0.95,
    size = 2.8
  ) +
  ggplot2::scale_shape_manual(
    values = shape_short,
    breaks = c("Core", "Source-supported", "Target-enhanced", "Source-discordant")
  ) +
  ggplot2::guides(
    shape = ggplot2::guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(size = 3.2)
    )
  ) +
  ggplot2::labs(
    title = "Source support vs target-domain stability",
    subtitle = "Background genes are faded; labels show the target-anchored reported set",
    x = "Target-only bootstrap stability",
    y = "Weighted source PIP",
    shape = NULL
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    legend.box.margin = ggplot2::margin(t = 4, b = 4),
    legend.key.width = grid::unit(1.1, "cm")
  )
if (nrow(label_df) > 0 && requireNamespace("ggrepel", quietly = TRUE)) {
  p_pip <- p_pip +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(
        x = TargetOnly_Stability,
        y = Source_PIP_Weighted,
        label = Raw_Gene
      ),
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE
    )
} else if (nrow(label_df) > 0) {
  p_pip <- p_pip +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(
        x = TargetOnly_Stability,
        y = Source_PIP_Weighted,
        label = Raw_Gene
      ),
      vjust = -0.8,
      size = 3,
      show.legend = FALSE
    )
}
save_plot_out(
  p_pip,
  "SA6_SourcePIP_vs_TargetStability",
  width = 8.8,
  height = 6.4
)
print(sa6_corr)
print(sa6_quadrant)
message("\nPrediction and transferability diagnostics complete.")