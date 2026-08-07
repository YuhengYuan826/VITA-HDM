# ==============================================================================
# 05_run_simulation_study.R
# Run: Rscript 05_run_simulation_study.R
# Run cross-platform VITA-HDM simulation study
# Analysis design:
#   1. Homogeneous:
#      Source and target share the same underlying mediation architecture,
#      but all source domains have platform-level mediator measurement distortion.
#
#   2. Heterogeneous:
#      Source domains have platform-level measurement distortion plus
#      structural perturbations: effect-size shifts, covariate-distribution shift,
#      and source-specific spurious mediators.
# Optional:
#   Structure-only transfer sensitivity is run on a reduced difficult grid.
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
SIMULATION_OUTPUT_DIR <- ensure_directory(PATHS$simulation_results)
suppressPackageStartupMessages({
  library(doParallel)
  library(foreach)
  library(dplyr)
})
SIMULATION_UTILS_FILE <- file.path(PROJECT_ROOT, "01_simulation_utils.R")
SIMULATION_DATA_FILE <- file.path(PROJECT_ROOT, "02_simulation_data.R")
SIMULATION_CAVI_FILE <- file.path(PROJECT_ROOT, "03_simulation_aft_cavi.R")
SIMULATION_MODEL_FILE <- file.path(PROJECT_ROOT, "04_simulation_vita_hdm.R")
source(SIMULATION_UTILS_FILE)
source(SIMULATION_DATA_FILE)
source(SIMULATION_CAVI_FILE)
source(SIMULATION_MODEL_FILE)
n_simu <- 100
n_cores <- 20
n_target <- 100
n_source <- 500
p_x <- 5
p_m_list <- c(1000, 2000)
rho_list <- c(0, 0.3, 0.6)
cens_rate_list <- c(0.2, 0.4)
h_hetero <- 5
# Main simulation design.
design_list <- c("Homogeneous", "Heterogeneous")
# Main simulation does not include structure-only rows.
# Structure-only is handled separately on a reduced difficult grid.
run_structure_only_in_main <- FALSE
run_structure_only_reduced_sensitivity <- TRUE
# Moderate platform-shift setting.
# This is intentionally stronger than the previous mild test,
# but not extreme enough to destroy source informativeness.
platform_params <- list(
  platform_frac = 0.35,
  platform_mean_sd = 0.50,
  platform_scale_sd = 0.30,
  platform_noise_sd = 0.20
)
run_one_setting <- function(design = c("Homogeneous", "Heterogeneous"),
                            p_m,
                            rho,
                            cens_rate,
                            n_simu_local = n_simu,
                            n_target_local = n_target,
                            n_source_local = n_source,
                            p_x_local = p_x,
                            h_hetero_local = h_hetero,
                            n_cores_local = n_cores,
                            include_structure_only = FALSE,
                            output_prefix = "Sim_Result") {
  design <- match.arg(design)
  # Platform shift is now a basic property of all source-domain transfer settings.
  # Structural heterogeneity is additionally imposed only in the Heterogeneous design.
  structural_heterogeneous <- design == "Heterogeneous"
  source_platform_shift <- TRUE
  cl <- parallel::makeCluster(n_cores_local)
  registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  res <- foreach(
    s = seq_len(n_simu_local),
    .combine = rbind,
    .packages = c("MASS", "survival"),
    .export = c(
      "n_target_local", "n_source_local", "p_x_local", "h_hetero_local",
      "generate_truth_params", "make_source_params", "generate_aft_dataset",
      "apply_platform_shift",
      "VITA_HDM_sim", "make_ar1_cov", "calibrate_censoring_rate",
      "prepare_simulation_data", "synthesize_priors_sim",
      "fit_aft_cavi", "get_bfdr_threshold", "aft_cindex",
      "safe_scale_train", "safe_scale_apply", "validate_sim_design",
      "platform_params",
      "SIMULATION_UTILS_FILE", "SIMULATION_DATA_FILE",
      "SIMULATION_CAVI_FILE", "SIMULATION_MODEL_FILE"
    )
  ) %dopar% {
    source(SIMULATION_UTILS_FILE)
    source(SIMULATION_DATA_FILE)
    source(SIMULATION_CAVI_FILE)
    source(SIMULATION_MODEL_FILE)
    params <- generate_truth_params(
      p_x = p_x_local,
      p_m = p_m,
      seed = 10000 + s
    )
    target_train <- generate_aft_dataset(
      n = n_target_local,
      params = params,
      rho = rho,
      censoring_rate = cens_rate,
      source_shift = FALSE,
      platform_shift = FALSE,
      seed = 20000 + s
    )
    target_test <- generate_aft_dataset(
      n = n_target_local,
      params = params,
      rho = rho,
      censoring_rate = cens_rate,
      source_shift = FALSE,
      platform_shift = FALSE,
      seed = 30000 + s
    )
    source_pool <- list()
    for (k in 1:3) {
      src_params <- make_source_params(
        params = params,
        h = ifelse(structural_heterogeneous, h_hetero_local, 0),
        heterogeneous = structural_heterogeneous,
        seed = 40000 + s * 10 + k
      )
      source_pool[[k]] <- generate_aft_dataset(
        n = n_source_local,
        params = src_params,
        rho = rho,
        censoring_rate = cens_rate,
        source_shift = structural_heterogeneous,
        platform_shift = source_platform_shift,
        platform_frac = platform_params$platform_frac,
        platform_mean_sd = platform_params$platform_mean_sd,
        platform_scale_sd = platform_params$platform_scale_sd,
        platform_noise_sd = platform_params$platform_noise_sd,
        seed = 50000 + s * 10 + k
      )$data
    }
    true_vars <- paste0("M", params$active_idx)
    get_alpha_estimates <- function(fit_obj) {
      X_train <- fit_obj$prep$train$X
      med_cols <- fit_obj$prep$mediator_cols
      D_vec <- X_train[, "D"]
      x_cov_cols <- grep("^X", colnames(X_train), value = TRUE)
      X_cov <- X_train[, x_cov_cols, drop = FALSE]
      M_mat <- X_train[, med_cols, drop = FALSE]
      if (ncol(X_cov) > 0) {
        resid_D <- lm.fit(cbind(1, X_cov), D_vec)$residuals
        resid_M <- lm.fit(cbind(1, X_cov), M_mat)$residuals
      } else {
        resid_D <- D_vec
        resid_M <- M_mat
      }
      denom <- sum(resid_D^2)
      if (is.na(denom) || denom < 1e-8) {
        alpha_hat <- rep(0, length(med_cols))
      } else {
        alpha_hat <- as.vector(t(resid_D) %*% resid_M / denom)
      }
      names(alpha_hat) <- med_cols
      alpha_hat
    }
    get_metrics <- function(fit_obj, method, n_src) {
      selected <- fit_obj$selected_mediators
      tp <- length(intersect(selected, true_vars))
      fp <- length(setdiff(selected, true_vars))
      recall <- tp / length(true_vars)
      fdr <- ifelse(tp + fp > 0, fp / (tp + fp), 0)
      precision <- 1 - fdr
      f1 <- ifelse(
        precision + recall > 0,
        2 * precision * recall / (precision + recall),
        0
      )
      scale_sd <- fit_obj$prep$scale_sd
      # ------------------------------------------------------------------------
      # Simulation convention:
      #   Intercept is explicit and unscaled.
      #   D, X, and M are standardized using target-training scaling.
      #   Y is not standardized.
      # Therefore:
      #   beta_D_scaled = beta1_raw * sd(D)
      #   beta_M_scaled = beta2_raw * sd(M)
      #   alpha_scaled  = alpha_raw * sd(D) / sd(M)
      #   IE_scaled_j   = alpha_raw * beta2_raw * sd(D)
      # ------------------------------------------------------------------------
      sd_D <- scale_sd["D"]
      m_names <- paste0("M", seq_along(params$beta2))
      sd_M <- scale_sd[m_names]
      true_beta_scaled <- params$beta2 * sd_M
      names(true_beta_scaled) <- m_names
      true_alpha_scaled <- params$alpha1 * sd_D / sd_M
      names(true_alpha_scaled) <- m_names
      true_de <- params$beta1 * sd_D
      true_ie <- sum(true_alpha_scaled * true_beta_scaled, na.rm = TRUE)
      true_te <- true_de + true_ie
      true_mp <- ifelse(abs(true_te) > 1e-8, true_ie / true_te, 0)
      inf <- fit_obj$inference
      beta_hat <- inf$Effect
      names(beta_hat) <- inf$Variable
      de_est <- beta_hat["D"]
      if (is.na(de_est)) de_est <- 0
      alpha_hat <- get_alpha_estimates(fit_obj)
      if (length(selected) > 0) {
        selected <- intersect(selected, names(alpha_hat))
        beta_selected <- beta_hat[selected]
        alpha_selected <- alpha_hat[selected]
        ie_est <- sum(alpha_selected * beta_selected, na.rm = TRUE)
      } else {
        ie_est <- 0
      }
      te_est <- de_est + ie_est
      mp_est <- ifelse(abs(te_est) > 1e-8, ie_est / te_est, 0)
      gene_inf <- inf[grepl("^M", inf$Variable), ]
      true_gene_beta <- true_beta_scaled[gene_inf$Variable]
      covered <- true_gene_beta >= gene_inf$Mean - 1.96 * gene_inf$SD &
        true_gene_beta <= gene_inf$Mean + 1.96 * gene_inf$SD
      active_mask <- true_gene_beta != 0
      null_mask <- true_gene_beta == 0
      coverage_active <- mean(covered[active_mask], na.rm = TRUE)
      coverage_null <- mean(covered[null_mask], na.rm = TRUE)
      data.frame(
        Design = design,
        p_m = p_m,
        rho = rho,
        Cens_Rate = cens_rate,
        Actual_Censor_Train = target_train$actual_censoring,
        Method = method,
        N_Source = n_src,
        Sim_ID = s,
        Recall = recall,
        FDR = fdr,
        F1 = f1,
        Cindex = fit_obj$test_cindex,
        DE_Est = as.numeric(de_est),
        IE_Est = as.numeric(ie_est),
        MP_Est = as.numeric(mp_est),
        True_DE = as.numeric(true_de),
        True_IE = as.numeric(true_ie),
        True_MP = as.numeric(true_mp),
        DE_Error = as.numeric(de_est - true_de),
        IE_Error = as.numeric(ie_est - true_ie),
        MP_Error = as.numeric(mp_est - true_mp),
        Coverage_Active = coverage_active,
        Coverage_Null = coverage_null
      )
    }
    rows <- list()
    row_id <- 1
    fit_hdm <- VITA_HDM_sim(
      target_train = target_train$data,
      source_list = NULL,
      target_test = target_test$data,
      transfer = FALSE,
      verbose = FALSE
    )
    rows[[row_id]] <- get_metrics(
      fit_obj = fit_hdm,
      method = "HDM",
      n_src = 0
    )
    row_id <- row_id + 1
    for (k in 1:3) {
      fit_vita <- VITA_HDM_sim(
        target_train = target_train$data,
        source_list = source_pool[seq_len(k)],
        target_test = target_test$data,
        transfer = TRUE,
        transfer_beta_strength = 1.0,
        verbose = FALSE
      )
      rows[[row_id]] <- get_metrics(
        fit_obj = fit_vita,
        method = "VITA-HDM",
        n_src = k
      )
      row_id <- row_id + 1
      if (include_structure_only) {
        fit_vita_struct <- VITA_HDM_sim(
          target_train = target_train$data,
          source_list = source_pool[seq_len(k)],
          target_test = target_test$data,
          transfer = TRUE,
          transfer_beta_strength = 0.0,
          verbose = FALSE
        )
        rows[[row_id]] <- get_metrics(
          fit_obj = fit_vita_struct,
          method = "VITA-HDM-StructOnly",
          n_src = k
        )
        row_id <- row_id + 1
      }
    }
    do.call(rbind, rows)
  }
  out_name <- sprintf(
    "%s_%s_p%d_rho%s_cens%s.csv",
    output_prefix,
    design,
    p_m,
    rho,
    cens_rate
  )
  out_path <- file.path(SIMULATION_OUTPUT_DIR, out_name)
  write.csv(res, out_path, row.names = FALSE)
  message("Saved: ", out_path)
  res
}
# ==============================================================================
# Main simulation
# ==============================================================================
all_results <- list()
for (design in design_list) {
  for (p_m in p_m_list) {
    for (rho in rho_list) {
      for (cens in cens_rate_list) {
        message(sprintf(
          "Running MAIN %s | p=%d | rho=%s | censor=%s",
          design,
          p_m,
          rho,
          cens
        ))
        all_results[[length(all_results) + 1]] <- run_one_setting(
          design = design,
          p_m = p_m,
          rho = rho,
          cens_rate = cens,
          include_structure_only = run_structure_only_in_main,
          output_prefix = "Sim_Result"
        )
      }
    }
  }
}
final_raw <- dplyr::bind_rows(all_results)
write.csv(
  final_raw,
  file.path(SIMULATION_OUTPUT_DIR, "Simulation_Raw_All.csv"),
  row.names = FALSE
)
message("========== Main simulation settings completed ==========")
message("Saved: Simulation_Raw_All.csv")
# ==============================================================================
# Optional reduced structure-only sensitivity
# This is intentionally limited to the difficult setting to reduce runtime.
# It should be summarized briefly, not plotted as a main figure.
# ==============================================================================
if (run_structure_only_reduced_sensitivity) {
  sens_results <- list()
  sens_design_list <- c("Homogeneous", "Heterogeneous")
  sens_p_m_list <- c(1000, 2000)
  sens_rho_list <- c(0.6)
  sens_cens_rate_list <- c(0.4)
  for (design in sens_design_list) {
    for (p_m in sens_p_m_list) {
      for (rho in sens_rho_list) {
        for (cens in sens_cens_rate_list) {
          message(sprintf(
            "Running STRUCTURE-ONLY SENSITIVITY %s | p=%d | rho=%s | censor=%s",
            design,
            p_m,
            rho,
            cens
          ))
          sens_results[[length(sens_results) + 1]] <- run_one_setting(
            design = design,
            p_m = p_m,
            rho = rho,
            cens_rate = cens,
            include_structure_only = TRUE,
            output_prefix = "Sensitivity_StructOnly"
          )
        }
      }
    }
  }
  sens_raw <- dplyr::bind_rows(sens_results)
  write.csv(
    sens_raw,
    file.path(SIMULATION_OUTPUT_DIR, "StructureOnly_Sensitivity_Raw.csv"),
    row.names = FALSE
  )
  message("========== Structure-only sensitivity completed ==========")
  message("Saved: StructureOnly_Sensitivity_Raw.csv")
}