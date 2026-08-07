# ==============================================================================
# 14_run_glioma_selection_sensitivity.R
# Run after files 09-11: Rscript 14_run_glioma_selection_sensitivity.R
# Target-anchored glioma sensitivity analyses:
# SA1: Noise knock-in
# SA2: Source dose-response
# SA3: BFDR x stability threshold grid
# SA7: Hyperparameter robustness grid
# SA8: Target-anchoring rule audit
# Uses 13_glioma_sensitivity_common.R. Failed AFT-CV folds are recorded as NA.
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
PROJECT_ROOT <- find_project_root(SCRIPT_DIR)
PATHS <- project_paths(PROJECT_ROOT)
source(file.path(PROJECT_ROOT, "13_glioma_sensitivity_common.R"))
# ------------------------------------------------------------------------------
# Runtime settings
# ------------------------------------------------------------------------------
SENS_BOOT_B <- min(BOOT_MAIN, 200)
# SA1
N_FAKE_MAIN <- 50
RUN_SA1_FAKE_COUNT_GRID <- FALSE
SA1_FAKE_COUNT_GRID <- c(25, 50, 100)
SA1_GRID_BOOT_B <- min(80, SENS_BOOT_B)
# SA2
DOSE_REPS <- 20
DOSE_LEVELS <- c(0.30, 0.60)
DOSE_BOOT_B <- min(100, SENS_BOOT_B)
DOSE_CV_FOLDS <- 5
# SA3
GRID_BFDR <- c(0.10, 0.20, 0.30)
GRID_STABILITY <- c(0.40, 0.50, 0.60)
SA3_BOOT_B <- SENS_BOOT_B
# SA7
RUN_SA7_HYPERPARAM_GRID <- TRUE
HYPER_BOOT_B <- 60
BASE_ETA_GRID <- sort(unique(c(BASE_ETA, 0.10, 0.15)))
FINAL_STR_GRID <- sort(unique(c(FINAL_STR - 0.50, FINAL_STR, FINAL_STR + 0.50)))
FINAL_STR_GRID <- FINAL_STR_GRID[is.finite(FINAL_STR_GRID) & FINAL_STR_GRID > 0]
BASE_TAU_GRID <- sort(unique(c(BASE_TAU, 0.25, 0.50, 1.00)))
TRANSFER_BETA_STRENGTH_GRID <- sort(unique(c(TRANSFER_BETA_STRENGTH, 0.00, 0.25, 0.50)))
# ==============================================================================
# SA1: Noise knock-in
# ==============================================================================
message("\n================================================================")
message("SA1: Noise knock-in sensitivity")
message("================================================================")
run_noise_knockin_once <- function(n_fake,
                                   boot_b,
                                   fake_seed = 9999,
                                   boot_seed_base = 4025,
                                   label = "main") {
  fake_raw_names <- paste0("Fake_Noise_", seq_len(n_fake))
  fake_gene_names <- paste0("M_", fake_raw_names)
  target_noise <- add_noise_genes_to_cohort(
    cohort = target,
    fake_names = fake_raw_names,
    seed = fake_seed + n_fake
  )
  source_noise <- add_noise_genes_to_cohort(
    cohort = source_cohort,
    fake_names = fake_raw_names,
    seed = fake_seed + n_fake + 100
  )
  keep_genes_noise <- c(keep_genes, fake_raw_names)
  gene_names_noise <- paste0("M_", keep_genes_noise)
  message(sprintf(
    ">>> SA1 %s | injected %d synthetic noise features | bootstrap B=%d",
    label,
    n_fake,
    boot_b
  ))
  src_noise <- fit_source_for_current(
    source_local = source_noise,
    keep_genes_local = keep_genes_noise,
    source_eta = SOURCE_ETA,
    tau_sq = BASE_TAU,
    max_iter = MAX_ITER
  )
  noise_full <- fit_vita_full_current(
    target_local = target_noise,
    keep_genes_local = keep_genes_noise,
    source_design_current = src_noise$source_design,
    source_fit_current = src_noise$source_fit,
    eta0 = BASE_ETA,
    source_eta0 = SOURCE_ETA,
    strength = FINAL_STR,
    tau_sq = BASE_TAU,
    transfer_beta_strength = TRANSFER_BETA_STRENGTH,
    weight_mode = WEIGHT_MODE,
    lpml_temperature = LPML_TEMPERATURE,
    max_iter = MAX_ITER
  )
  noise_boot <- run_vita_bootstrap_selection(
    target_local = target_noise,
    keep_genes_local = keep_genes_noise,
    source_design_current = src_noise$source_design,
    source_fit_current = src_noise$source_fit,
    boot_b = boot_b,
    target_bfdr = TARGET_BFDR,
    eta0 = BASE_ETA,
    source_eta0 = SOURCE_ETA,
    strength = FINAL_STR,
    tau_sq = BASE_TAU,
    transfer_beta_strength = TRANSFER_BETA_STRENGTH,
    weight_mode = WEIGHT_MODE,
    lpml_temperature = LPML_TEMPERATURE,
    seed_base = boot_seed_base + n_fake * 10,
    store_pip = TRUE,
    max_iter = MAX_ITER,
    show_progress = TRUE
  )
  target_stability_noise <- rep(0, length(gene_names_noise))
  names(target_stability_noise) <- gene_names_noise
  common_base <- intersect(names(baseline_target_stability), gene_names_noise)
  target_stability_noise[common_base] <- baseline_target_stability[common_base]
  source_pip_noise <- source_pip_from_fit(
    fit_obj = noise_full$fit,
    gene_names = gene_names_noise
  )
  class_noise <- apply_target_anchored_rule(
    gene_names = gene_names_noise,
    vita_stability = noise_boot$stability,
    target_stability = target_stability_noise,
    source_pip = source_pip_noise,
    source_support_cutoff = SOURCE_SUPPORT_CUTOFF,
    extended_cutoff = EXTENDED_STABILITY_CUTOFF,
    core_cutoff = CORE_STABILITY_CUTOFF,
    targetonly_cutoff = TARGET_ONLY_SOURCE_DISCORDANT_CUTOFF
  )
  recovery <- summarise_classification_recovery(
    class_df = class_noise,
    label = label
  )
  selected_target_anchored <- class_noise$Gene[class_noise$TargetAnchored_Reported]
  selected_fake <- intersect(selected_target_anchored, fake_gene_names)
  recovery$Injected_Noise_Features <- n_fake
  recovery$Bootstrap_B <- boot_b
  recovery$Noise_Selected <- length(selected_fake)
  recovery$Noise_Selection_Rate <- length(selected_fake) / n_fake
  recovery$Background_Selected_Not_In_TargetAnchored <- length(
    setdiff(selected_target_anchored, c(target_anchored_signature, fake_gene_names))
  )
  gene_report <- class_noise %>%
    mutate(
      Noise_Test_Label = label,
      Injected_Noise_Features = n_fake,
      Is_Fake = Gene %in% fake_gene_names,
      In_TargetAnchored_Reference = Gene %in% target_anchored_signature,
      In_TransferSupported_Reference = Gene %in% transfer_supported_signature,
      In_SourceSupported_Reference = Gene %in% source_supported_signature,
      In_Core_Reference = Gene %in% core_signature,
      In_TargetOnly_Reference = Gene %in% target_only_signature,
      In_TargetDiscordant_Reference = Gene %in% target_discordant_signature,
      Selected_TargetAnchored = TargetAnchored_Reported,
      Selected_TransferSupported = TransferSupported_Signature,
      Status = case_when(
        Is_Fake & TargetAnchored_Reported ~ "Synthetic noise selected",
        Is_Fake & !TargetAnchored_Reported ~ "Synthetic noise rejected",
        In_Core_Reference & TargetAnchored_Reported ~ "Core retained",
        In_TransferSupported_Reference & TargetAnchored_Reported ~ "Transfer-supported retained",
        In_TargetDiscordant_Reference & TargetAnchored_Reported ~ "Target-discordant retained",
        In_TargetAnchored_Reference & TargetAnchored_Reported ~ "Target-anchored retained",
        In_TargetAnchored_Reference & !TargetAnchored_Reported ~ "Target-anchored not retained",
        TargetAnchored_Reported ~ "Additional background selected",
        TRUE ~ "Background"
      )
    ) %>%
    arrange(
      desc(Selected_TargetAnchored),
      desc(Is_Fake),
      desc(VITA_Stability),
      Raw_Gene
    )
  list(
    summary = recovery,
    gene_report = gene_report,
    diagnostics = noise_boot$diagnostics
  )
}
sa1_main <- run_noise_knockin_once(
  n_fake = N_FAKE_MAIN,
  boot_b = SENS_BOOT_B,
  fake_seed = 9999,
  boot_seed_base = 4025,
  label = paste0("Main_", N_FAKE_MAIN, "_pseudo_genes")
)
sa1_summary_long <- data.frame(
  Metric = names(sa1_main$summary),
  Value = as.character(unlist(sa1_main$summary, use.names = FALSE)),
  stringsAsFactors = FALSE
)
write_csv_out(sa1_summary_long, "SA1_NoiseKnockin_Summary.csv")
write_csv_out(sa1_main$summary, "SA1_NoiseKnockin_MainSummary.csv")
write_csv_out(sa1_main$gene_report, "SA1_NoiseKnockin_GeneLevel.csv")
write_csv_out(sa1_main$diagnostics, "SA1_NoiseKnockin_Bootstrap_Diagnostics.csv")
if (isTRUE(RUN_SA1_FAKE_COUNT_GRID)) {
  grid_summaries <- list()
  grid_gene_reports <- list()
  for (nf in SA1_FAKE_COUNT_GRID) {
    if (nf == N_FAKE_MAIN) {
      grid_summaries[[as.character(nf)]] <- sa1_main$summary
      grid_gene_reports[[as.character(nf)]] <- sa1_main$gene_report
    } else {
      tmp <- run_noise_knockin_once(
        n_fake = nf,
        boot_b = SA1_GRID_BOOT_B,
        fake_seed = 9999,
        boot_seed_base = 4025,
        label = paste0("Grid_", nf, "_pseudo_genes")
      )
      grid_summaries[[as.character(nf)]] <- tmp$summary
      grid_gene_reports[[as.character(nf)]] <- tmp$gene_report
    }
  }
  write_csv_out(
    bind_rows(grid_summaries),
    "SA1_NoiseKnockin_FakeCountGrid_Summary.csv"
  )
  write_csv_out(
    bind_rows(grid_gene_reports),
    "SA1_NoiseKnockin_GeneLevel_AllLevels.csv"
  )
}
print(sa1_main$summary)
# ==============================================================================
# SA2: Source dose-response
# ==============================================================================
message("\n================================================================")
message("SA2: Source dose-response sensitivity")
message("================================================================")
dose_rows <- list()
dose_id <- 1
# ------------------------------------------------------------------------------
# Full-source reference.
# The 100% source condition is the main full-source analysis benchmark rather
# than a repeated downsampling analysis. 30% and 60% are repeatedly subsampled
# because their selected source samples vary by random seed. The 100% row uses
# the final full TCGA_LGG source object from the discovery analysis, and it
# also computes a
# full-source AFT-CV C-index so that the prediction panel has a valid 100% point.
# ------------------------------------------------------------------------------
message(">>> Computing 100% full-source reference for SA2")
source_design_ref <- if (!is.null(obj$source_design)) {
  obj$source_design
} else {
  prepare_source_design(
    source_cohort = source_cohort,
    keep_genes = keep_genes,
    covariates = COVARIATES,
    binary_covariates = BINARY_COVARIATES
  )
}
source_fit_ref <- if (!is.null(obj$source_fit)) {
  obj$source_fit
} else {
  fit_source_model(
    source_design = source_design_ref,
    source_name = SOURCE_NAME,
    source_eta0 = SOURCE_ETA,
    tau_sq = BASE_TAU,
    max_iter = MAX_ITER,
    verbose = FALSE
  )
}
c_vita_ref <- evaluate_vita_aft_cv_current(
  target_local = target,
  keep_genes_local = keep_genes,
  source_design_current = source_design_ref,
  source_fit_current = source_fit_ref,
  n_folds = DOSE_CV_FOLDS,
  seed = 9100,
  eta0 = BASE_ETA,
  source_eta0 = SOURCE_ETA,
  strength = FINAL_STR,
  tau_sq = BASE_TAU,
  transfer_beta_strength = TRANSFER_BETA_STRENGTH,
  weight_mode = WEIGHT_MODE,
  lpml_temperature = LPML_TEMPERATURE,
  max_iter = MAX_ITER
)
baseline_ref <- summarise_classification_recovery(
  class_df = baseline_classification,
  label = "100% full-source reference",
  cindex = mean(c_vita_ref, na.rm = TRUE),
  cindex_sd = safe_sd_local(c_vita_ref)
)
baseline_ref$Dose <- "100%"
baseline_ref$Dose_Proportion <- 1.00
baseline_ref$Repetition <- 0L
baseline_ref$N_Source <- nrow(source_cohort$clinical)
baseline_ref$Analysis_Type <- "Full-source reference"
baseline_ref$SA2_Fit_Success <- all(is.finite(c_vita_ref))
baseline_ref$SA2_Error_Message <- ifelse(
  baseline_ref$SA2_Fit_Success,
  "",
  "One or more full-source AFT-CV folds returned NA."
)
dose_rows[[dose_id]] <- baseline_ref
dose_id <- dose_id + 1
for (dose in DOSE_LEVELS) {
  message(sprintf(
    ">>> Reduced source dose %.0f%% | repetitions=%d | bootstrap B=%d",
    dose * 100,
    DOSE_REPS,
    DOSE_BOOT_B
  ))
  for (r in seq_len(DOSE_REPS)) {
    seed_now <- 5200 + round(dose * 100) * 10 + r
    source_sub <- subsample_source(
      source_local = source_cohort,
      dose = dose,
      seed = seed_now,
      min_n = 20
    )
    n_source_sub <- nrow(source_sub$clinical)
    message(sprintf("  -> Dose %.0f%% repetition %d/%d", dose * 100, r, DOSE_REPS))
    rec_sub <- try({
      src_sub <- fit_source_for_current(
        source_local = source_sub,
        keep_genes_local = keep_genes,
        source_eta = SOURCE_ETA,
        tau_sq = BASE_TAU,
        max_iter = MAX_ITER
      )
      full_sub <- fit_vita_full_current(
        target_local = target,
        keep_genes_local = keep_genes,
        source_design_current = src_sub$source_design,
        source_fit_current = src_sub$source_fit,
        eta0 = BASE_ETA,
        source_eta0 = SOURCE_ETA,
        strength = FINAL_STR,
        tau_sq = BASE_TAU,
        transfer_beta_strength = TRANSFER_BETA_STRENGTH,
        weight_mode = WEIGHT_MODE,
        lpml_temperature = LPML_TEMPERATURE,
        max_iter = MAX_ITER
      )
      dose_boot <- run_vita_bootstrap_selection(
        target_local = target,
        keep_genes_local = keep_genes,
        source_design_current = src_sub$source_design,
        source_fit_current = src_sub$source_fit,
        boot_b = DOSE_BOOT_B,
        target_bfdr = TARGET_BFDR,
        eta0 = BASE_ETA,
        source_eta0 = SOURCE_ETA,
        strength = FINAL_STR,
        tau_sq = BASE_TAU,
        transfer_beta_strength = TRANSFER_BETA_STRENGTH,
        weight_mode = WEIGHT_MODE,
        lpml_temperature = LPML_TEMPERATURE,
        seed_base = 6200 + round(dose * 100) * 10 + r * 1000,
        store_pip = FALSE,
        max_iter = MAX_ITER,
        show_progress = TRUE
      )
      source_pip_sub <- source_pip_from_fit(
        fit_obj = full_sub$fit,
        gene_names = available_gene_names
      )
      class_sub <- apply_target_anchored_rule(
        gene_names = available_gene_names,
        vita_stability = dose_boot$stability,
        target_stability = baseline_target_stability,
        source_pip = source_pip_sub,
        source_support_cutoff = SOURCE_SUPPORT_CUTOFF,
        extended_cutoff = EXTENDED_STABILITY_CUTOFF,
        core_cutoff = CORE_STABILITY_CUTOFF,
        targetonly_cutoff = TARGET_ONLY_SOURCE_DISCORDANT_CUTOFF
      )
      c_vita <- evaluate_vita_aft_cv_current(
        target_local = target,
        keep_genes_local = keep_genes,
        source_design_current = src_sub$source_design,
        source_fit_current = src_sub$source_fit,
        n_folds = DOSE_CV_FOLDS,
        seed = 8100 + r + round(dose * 100) * 10,
        eta0 = BASE_ETA,
        source_eta0 = SOURCE_ETA,
        strength = FINAL_STR,
        tau_sq = BASE_TAU,
        transfer_beta_strength = TRANSFER_BETA_STRENGTH,
        weight_mode = WEIGHT_MODE,
        lpml_temperature = LPML_TEMPERATURE,
        max_iter = MAX_ITER
      )
      rec <- summarise_classification_recovery(
        class_df = class_sub,
        label = paste0("Dose_", round(dose * 100), "_rep_", r),
        cindex = mean(c_vita, na.rm = TRUE),
        cindex_sd = safe_sd_local(c_vita)
      )
      rec$Dose <- sprintf("%d%%", round(dose * 100))
      rec$Dose_Proportion <- dose
      rec$Repetition <- r
      rec$N_Source <- n_source_sub
      rec$Analysis_Type <- "Reduced-source repeated subsampling"
      rec$SA2_Fit_Success <- TRUE
      rec$SA2_Error_Message <- ""
      rec
    }, silent = TRUE)
    if (inherits(rec_sub, "try-error")) {
      err_msg <- as.character(attr(rec_sub, "condition")$message)
      rec_sub <- data.frame(
        Analysis_Label = paste0("Dose_", round(dose * 100), "_rep_", r),
        N_TargetAnchored_Selected = NA_integer_,
        N_TransferSupported_Selected = NA_integer_,
        N_SourceSupported_Selected = NA_integer_,
        N_Core_Selected = NA_integer_,
        N_TargetOnlyStable_Selected = NA_integer_,
        N_TargetDiscordant_Selected = NA_integer_,
        TargetAnchored_Recall = NA_real_,
        TransferSupported_Recall = NA_real_,
        SourceSupported_Recall = NA_real_,
        Core_Recall = NA_real_,
        TargetOnlyStable_Recall = NA_real_,
        TargetDiscordant_Retention = NA_real_,
        Jaccard_vs_TargetAnchored = NA_real_,
        Jaccard_vs_TransferSupported = NA_real_,
        SourceSupported_Fraction_TargetAnchored = NA_real_,
        Mean_C_Index = NA_real_,
        SD_C_Index = NA_real_,
        TargetAnchored_Genes = NA_character_,
        TransferSupported_Genes = NA_character_,
        TargetDiscordant_Genes = NA_character_,
        Dose = sprintf("%d%%", round(dose * 100)),
        Dose_Proportion = dose,
        Repetition = r,
        N_Source = n_source_sub,
        Analysis_Type = "Reduced-source repeated subsampling",
        SA2_Fit_Success = FALSE,
        SA2_Error_Message = err_msg,
        stringsAsFactors = FALSE
      )
    }
    dose_rows[[dose_id]] <- rec_sub
    dose_id <- dose_id + 1
  }
}
dose_detailed <- bind_rows(dose_rows)
dose_summary_reduced <- dose_detailed %>%
  filter(Analysis_Type == "Reduced-source repeated subsampling") %>%
  group_by(Dose, Dose_Proportion, Analysis_Type) %>%
  summarise(
    MC_Runs = n(),
    Successful_Runs = sum(SA2_Fit_Success %in% TRUE, na.rm = TRUE),
    Avg_N_Source = mean(N_Source, na.rm = TRUE),
    Avg_C_Index = mean(Mean_C_Index, na.rm = TRUE),
    SD_C_Index = safe_sd_local(Mean_C_Index),
    Avg_TargetAnchored_Recall = mean(TargetAnchored_Recall, na.rm = TRUE),
    SD_TargetAnchored_Recall = safe_sd_local(TargetAnchored_Recall),
    Avg_TransferSupported_Recall = mean(TransferSupported_Recall, na.rm = TRUE),
    SD_TransferSupported_Recall = safe_sd_local(TransferSupported_Recall),
    Avg_SourceSupported_Recall = mean(SourceSupported_Recall, na.rm = TRUE),
    SD_SourceSupported_Recall = safe_sd_local(SourceSupported_Recall),
    Avg_Core_Recall = mean(Core_Recall, na.rm = TRUE),
    SD_Core_Recall = safe_sd_local(Core_Recall),
    Avg_TargetOnlyStable_Recall = mean(TargetOnlyStable_Recall, na.rm = TRUE),
    SD_TargetOnlyStable_Recall = safe_sd_local(TargetOnlyStable_Recall),
    Avg_TargetDiscordant_Retention = mean(TargetDiscordant_Retention, na.rm = TRUE),
    SD_TargetDiscordant_Retention = safe_sd_local(TargetDiscordant_Retention),
    Avg_TargetAnchored_Size = mean(N_TargetAnchored_Selected, na.rm = TRUE),
    Avg_TransferSupported_Size = mean(N_TransferSupported_Selected, na.rm = TRUE),
    Avg_SourceSupported_Fraction_TargetAnchored =
      mean(SourceSupported_Fraction_TargetAnchored, na.rm = TRUE),
    .groups = "drop"
  )
dose_summary_reference <- dose_detailed %>%
  filter(Analysis_Type == "Full-source reference") %>%
  transmute(
    Dose,
    Dose_Proportion,
    Analysis_Type,
    MC_Runs = 1L,
    Successful_Runs = as.integer(SA2_Fit_Success %in% TRUE),
    Avg_N_Source = as.numeric(N_Source),
    Avg_C_Index = Mean_C_Index,
    SD_C_Index = NA_real_,
    Avg_TargetAnchored_Recall = TargetAnchored_Recall,
    SD_TargetAnchored_Recall = NA_real_,
    Avg_TransferSupported_Recall = TransferSupported_Recall,
    SD_TransferSupported_Recall = NA_real_,
    Avg_SourceSupported_Recall = SourceSupported_Recall,
    SD_SourceSupported_Recall = NA_real_,
    Avg_Core_Recall = Core_Recall,
    SD_Core_Recall = NA_real_,
    Avg_TargetOnlyStable_Recall = TargetOnlyStable_Recall,
    SD_TargetOnlyStable_Recall = NA_real_,
    Avg_TargetDiscordant_Retention = TargetDiscordant_Retention,
    SD_TargetDiscordant_Retention = NA_real_,
    Avg_TargetAnchored_Size = as.numeric(N_TargetAnchored_Selected),
    Avg_TransferSupported_Size = as.numeric(N_TransferSupported_Selected),
    Avg_SourceSupported_Fraction_TargetAnchored = SourceSupported_Fraction_TargetAnchored
  )
dose_summary <- bind_rows(
  dose_summary_reduced,
  dose_summary_reference
) %>%
  arrange(Dose_Proportion)
write_csv_out(dose_detailed, "SA2_SourceDose_Detailed.csv")
write_csv_out(dose_summary, "SA2_SourceDose_Summary.csv")
print(dose_summary)
# ==============================================================================
# SA3: BFDR x stability threshold grid
# ==============================================================================
message("\n================================================================")
message("SA3: BFDR x stability threshold grid")
message("================================================================")
src_main <- fit_source_for_current(
  source_local = source_cohort,
  keep_genes_local = keep_genes,
  source_eta = SOURCE_ETA,
  tau_sq = BASE_TAU,
  max_iter = MAX_ITER
)
full_main_refit <- fit_vita_full_current(
  target_local = target,
  keep_genes_local = keep_genes,
  source_design_current = src_main$source_design,
  source_fit_current = src_main$source_fit,
  eta0 = BASE_ETA,
  source_eta0 = SOURCE_ETA,
  strength = FINAL_STR,
  tau_sq = BASE_TAU,
  transfer_beta_strength = TRANSFER_BETA_STRENGTH,
  weight_mode = WEIGHT_MODE,
  lpml_temperature = LPML_TEMPERATURE,
  max_iter = MAX_ITER
)
sa3_source_pip <- source_pip_from_fit(
  fit_obj = full_main_refit$fit,
  gene_names = available_gene_names
)
sa3_boot <- run_vita_bootstrap_selection(
  target_local = target,
  keep_genes_local = keep_genes,
  source_design_current = src_main$source_design,
  source_fit_current = src_main$source_fit,
  boot_b = SA3_BOOT_B,
  target_bfdr = TARGET_BFDR,
  eta0 = BASE_ETA,
  source_eta0 = SOURCE_ETA,
  strength = FINAL_STR,
  tau_sq = BASE_TAU,
  transfer_beta_strength = TRANSFER_BETA_STRENGTH,
  weight_mode = WEIGHT_MODE,
  lpml_temperature = LPML_TEMPERATURE,
  seed_base = 9300,
  store_pip = TRUE,
  max_iter = MAX_ITER,
  show_progress = TRUE
)
if (is.null(sa3_boot$pip_mat)) {
  stop("SA3 requires pip_mat, but pip_mat is NULL.")
}
select_genes_from_pip_bfdr <- function(pip_vec, bfdr) {
  pip_vec <- pip_vec[is.finite(pip_vec)]
  if (length(pip_vec) == 0) return(character(0))
  sorted <- sort(pip_vec, decreasing = TRUE)
  fdr_curve <- cumsum(1 - sorted) / seq_along(sorted)
  valid <- which(fdr_curve <= bfdr)
  if (length(valid) == 0) return(character(0))
  names(sorted)[seq_len(max(valid))]
}
sa3_rows <- list()
sa3_gene_rows <- list()
row_id <- 1
gene_id <- 1
for (bfdr in GRID_BFDR) {
  sel_mat <- matrix(
    0,
    nrow = nrow(sa3_boot$pip_mat),
    ncol = ncol(sa3_boot$pip_mat)
  )
  colnames(sel_mat) <- colnames(sa3_boot$pip_mat)
  for (b in seq_len(nrow(sa3_boot$pip_mat))) {
    pip_b <- sa3_boot$pip_mat[b, ]
    names(pip_b) <- colnames(sa3_boot$pip_mat)
    sel_b <- select_genes_from_pip_bfdr(pip_b, bfdr)
    sel_b <- intersect(sel_b, colnames(sel_mat))
    if (length(sel_b) > 0) {
      sel_mat[b, sel_b] <- 1
    }
  }
  stability_bfdr <- colMeans(sel_mat, na.rm = TRUE)
  names(stability_bfdr) <- colnames(sel_mat)
  for (stab in GRID_STABILITY) {
    class_grid <- apply_target_anchored_rule(
      gene_names = available_gene_names,
      vita_stability = stability_bfdr,
      target_stability = baseline_target_stability,
      source_pip = sa3_source_pip,
      source_support_cutoff = SOURCE_SUPPORT_CUTOFF,
      extended_cutoff = stab,
      core_cutoff = CORE_STABILITY_CUTOFF,
      targetonly_cutoff = stab
    )
    rec_grid <- summarise_classification_recovery(
      class_df = class_grid,
      label = paste0("BFDR_", bfdr, "_STAB_", stab)
    )
    rec_grid$TARGET_BFDR <- bfdr
    rec_grid$STABILITY_CUTOFF <- stab
    rec_grid$Bootstrap_B <- SA3_BOOT_B
    sa3_rows[[row_id]] <- rec_grid
    row_id <- row_id + 1
    gene_table <- class_grid %>%
      mutate(
        TARGET_BFDR = bfdr,
        STABILITY_CUTOFF = stab,
        Bootstrap_B = SA3_BOOT_B,
        In_TargetAnchored_Reference = Gene %in% target_anchored_signature,
        In_TransferSupported_Reference = Gene %in% transfer_supported_signature,
        In_TargetDiscordant_Reference = Gene %in% target_discordant_signature
      )
    sa3_gene_rows[[gene_id]] <- gene_table
    gene_id <- gene_id + 1
  }
}
sa3_summary <- bind_rows(sa3_rows) %>%
  arrange(TARGET_BFDR, STABILITY_CUTOFF)
sa3_genelevel <- bind_rows(sa3_gene_rows)
write_csv_out(sa3_summary, "SA3_ThresholdGrid_Summary.csv")
write_csv_out(sa3_genelevel, "SA3_ThresholdGrid_GeneLevel.csv")
write_csv_out(sa3_boot$diagnostics, "SA3_ThresholdGrid_Bootstrap_Diagnostics.csv")
print(sa3_summary)
# ==============================================================================
# SA7: Hyperparameter robustness grid
# ==============================================================================
message("\n================================================================")
message("SA7: Hyperparameter robustness grid")
message("================================================================")
if (isTRUE(RUN_SA7_HYPERPARAM_GRID)) {
  hyper_grid <- expand.grid(
    BASE_ETA = BASE_ETA_GRID,
    FINAL_STR = FINAL_STR_GRID,
    BASE_TAU = BASE_TAU_GRID,
    TRANSFER_BETA_STRENGTH = TRANSFER_BETA_STRENGTH_GRID,
    stringsAsFactors = FALSE
  ) %>%
    distinct()
  hyper_rows <- list()
  hyper_gene_rows <- list()
  row_id <- 1
  gene_id <- 1
  for (ii in seq_len(nrow(hyper_grid))) {
    hp <- hyper_grid[ii, ]
    message(sprintf(
      ">>> SA7 cell %d/%d | BASE_ETA=%.3f | FINAL_STR=%.3f | TAU=%.3f | BETA_STR=%.3f",
      ii,
      nrow(hyper_grid),
      hp$BASE_ETA,
      hp$FINAL_STR,
      hp$BASE_TAU,
      hp$TRANSFER_BETA_STRENGTH
    ))
    tmp <- try({
      src_hp <- fit_source_for_current(
        source_local = source_cohort,
        keep_genes_local = keep_genes,
        source_eta = SOURCE_ETA,
        tau_sq = hp$BASE_TAU,
        max_iter = MAX_ITER
      )
      full_hp <- fit_vita_full_current(
        target_local = target,
        keep_genes_local = keep_genes,
        source_design_current = src_hp$source_design,
        source_fit_current = src_hp$source_fit,
        eta0 = hp$BASE_ETA,
        source_eta0 = SOURCE_ETA,
        strength = hp$FINAL_STR,
        tau_sq = hp$BASE_TAU,
        transfer_beta_strength = hp$TRANSFER_BETA_STRENGTH,
        weight_mode = WEIGHT_MODE,
        lpml_temperature = LPML_TEMPERATURE,
        max_iter = MAX_ITER
      )
      boot_hp <- run_vita_bootstrap_selection(
        target_local = target,
        keep_genes_local = keep_genes,
        source_design_current = src_hp$source_design,
        source_fit_current = src_hp$source_fit,
        boot_b = HYPER_BOOT_B,
        target_bfdr = TARGET_BFDR,
        eta0 = hp$BASE_ETA,
        source_eta0 = SOURCE_ETA,
        strength = hp$FINAL_STR,
        tau_sq = hp$BASE_TAU,
        transfer_beta_strength = hp$TRANSFER_BETA_STRENGTH,
        weight_mode = WEIGHT_MODE,
        lpml_temperature = LPML_TEMPERATURE,
        seed_base = 11000 + ii * 100,
        store_pip = FALSE,
        max_iter = MAX_ITER,
        show_progress = TRUE
      )
      source_pip_hp <- source_pip_from_fit(
        fit_obj = full_hp$fit,
        gene_names = available_gene_names
      )
      class_hp <- apply_target_anchored_rule(
        gene_names = available_gene_names,
        vita_stability = boot_hp$stability,
        target_stability = baseline_target_stability,
        source_pip = source_pip_hp,
        source_support_cutoff = SOURCE_SUPPORT_CUTOFF,
        extended_cutoff = EXTENDED_STABILITY_CUTOFF,
        core_cutoff = CORE_STABILITY_CUTOFF,
        targetonly_cutoff = TARGET_ONLY_SOURCE_DISCORDANT_CUTOFF
      )
      rec_hp <- summarise_classification_recovery(
        class_df = class_hp,
        label = paste0("HP_Cell_", ii)
      )
      rec_hp <- cbind(hp, rec_hp)
      rec_hp$Fit_Success <- TRUE
      rec_hp$Error_Message <- ""
      gene_hp <- class_hp %>%
        mutate(
          Hyper_Cell = ii,
          BASE_ETA = hp$BASE_ETA,
          FINAL_STR = hp$FINAL_STR,
          BASE_TAU = hp$BASE_TAU,
          TRANSFER_BETA_STRENGTH = hp$TRANSFER_BETA_STRENGTH,
          In_TargetAnchored_Reference = Gene %in% target_anchored_signature,
          In_TransferSupported_Reference = Gene %in% transfer_supported_signature,
          In_TargetDiscordant_Reference = Gene %in% target_discordant_signature
        )
      list(summary = rec_hp, gene = gene_hp)
    }, silent = TRUE)
    if (inherits(tmp, "try-error")) {
      err_msg <- as.character(attr(tmp, "condition")$message)
      rec_hp <- cbind(
        hp,
        data.frame(
          Analysis_Label = paste0("HP_Cell_", ii),
          N_TargetAnchored_Selected = NA_integer_,
          N_TransferSupported_Selected = NA_integer_,
          N_SourceSupported_Selected = NA_integer_,
          N_Core_Selected = NA_integer_,
          N_TargetOnlyStable_Selected = NA_integer_,
          N_TargetDiscordant_Selected = NA_integer_,
          TargetAnchored_Recall = NA_real_,
          TransferSupported_Recall = NA_real_,
          SourceSupported_Recall = NA_real_,
          Core_Recall = NA_real_,
          TargetOnlyStable_Recall = NA_real_,
          TargetDiscordant_Retention = NA_real_,
          Jaccard_vs_TargetAnchored = NA_real_,
          Jaccard_vs_TransferSupported = NA_real_,
          SourceSupported_Fraction_TargetAnchored = NA_real_,
          Mean_C_Index = NA_real_,
          SD_C_Index = NA_real_,
          TargetAnchored_Genes = NA_character_,
          TransferSupported_Genes = NA_character_,
          TargetDiscordant_Genes = NA_character_,
          Fit_Success = FALSE,
          Error_Message = err_msg,
          stringsAsFactors = FALSE
        )
      )
      hyper_rows[[row_id]] <- rec_hp
      row_id <- row_id + 1
    } else {
      hyper_rows[[row_id]] <- tmp$summary
      row_id <- row_id + 1
      hyper_gene_rows[[gene_id]] <- tmp$gene
      gene_id <- gene_id + 1
    }
  }
  sa7_summary <- bind_rows(hyper_rows)
  sa7_genelevel <- if (length(hyper_gene_rows) > 0) {
    bind_rows(hyper_gene_rows)
  } else {
    data.frame()
  }
  write_csv_out(sa7_summary, "SA7_HyperparameterGrid_Summary.csv")
  write_csv_out(sa7_genelevel, "SA7_HyperparameterGrid_GeneLevel.csv")
  print(sa7_summary)
}
# ==============================================================================
# SA8: Target-anchoring rule audit
# ==============================================================================
message("\n================================================================")
message("SA8: Target-anchoring rule audit")
message("================================================================")
sa8_gene_membership <- baseline_classification %>%
  mutate(
    In_TargetAnchored_Reference = Gene %in% target_anchored_signature,
    In_TransferSupported_Reference = Gene %in% transfer_supported_signature,
    In_SourceSupported_Reference = Gene %in% source_supported_signature,
    In_Core_Reference = Gene %in% core_signature,
    In_TargetOnly_Reference = Gene %in% target_only_signature,
    In_TargetDiscordant_Reference = Gene %in% target_discordant_signature,
    In_TargetEnhanced_Reference = Gene %in% target_enhanced_signature
  ) %>%
  filter(
    In_TargetAnchored_Reference |
      In_TransferSupported_Reference |
      In_SourceSupported_Reference |
      In_Core_Reference |
      In_TargetOnly_Reference |
      In_TargetDiscordant_Reference |
      In_TargetEnhanced_Reference
  ) %>%
  arrange(
    desc(In_TargetAnchored_Reference),
    desc(In_TransferSupported_Reference),
    desc(In_TargetOnly_Reference),
    desc(VITA_Stability),
    Raw_Gene
  )
sa8_rule_summary <- data.frame(
  Rule_Set = c(
    "Transfer-supported set",
    "Target-only stable set",
    "Target-only source-discordant subset",
    "Target-anchored reported set"
  ),
  N_Genes = c(
    length(transfer_supported_signature),
    length(target_only_signature),
    length(target_discordant_signature),
    length(target_anchored_signature)
  ),
  Genes = c(
    ifelse(
      length(transfer_supported_signature) > 0,
      paste(transfer_supported_signature, collapse = " | "),
      "None"
    ),
    ifelse(
      length(target_only_signature) > 0,
      paste(target_only_signature, collapse = " | "),
      "None"
    ),
    ifelse(
      length(target_discordant_signature) > 0,
      paste(target_discordant_signature, collapse = " | "),
      "None"
    ),
    ifelse(
      length(target_anchored_signature) > 0,
      paste(target_anchored_signature, collapse = " | "),
      "None"
    )
  ),
  Interpretation = c(
    "Genes retained by transfer-supported or target-enhanced transfer-retained evidence.",
    "Genes reproducibly selected by target-only analysis.",
    "Target-only stable genes with weak source support and insufficient VITA stability.",
    "Final reported set under target-anchored rule; source evidence annotates transferability but does not veto target-only stability."
  ),
  stringsAsFactors = FALSE
)
cv_summary_file <- file.path(CV_DIR, "glioma_fixed_signature_cindex_summary.csv")
effect_summary_file <- file.path(EFF_DIR, "glioma_overall_effects_conditional_bootstrap.csv")
sa8_cv <- if (file.exists(cv_summary_file)) {
  read.csv(cv_summary_file, stringsAsFactors = FALSE)
} else {
  data.frame()
}
sa8_effects <- if (file.exists(effect_summary_file)) {
  read.csv(effect_summary_file, stringsAsFactors = FALSE)
} else {
  data.frame()
}
write_csv_out(sa8_gene_membership, "SA8_TargetAnchoring_GeneMembership.csv")
write_csv_out(sa8_rule_summary, "SA8_TargetAnchoring_RuleSummary.csv")
write_csv_out(sa8_cv, "SA8_TargetAnchoring_CIndexComparison.csv")
write_csv_out(sa8_effects, "SA8_TargetAnchoring_EffectComparison.csv")
print(sa8_rule_summary)
message("\nSelection robustness sensitivity analyses complete.")