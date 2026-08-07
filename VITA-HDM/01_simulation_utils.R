# ==============================================================================
# 01_simulation_utils.R
# Utility functions for VITA-HDM simulation study
# ==============================================================================
suppressPackageStartupMessages({
  library(MASS)
  library(survival)
})
make_ar1_cov <- function(p, rho) {
  toeplitz(rho^(0:(p - 1)))
}
safe_scale_train <- function(X_train) {
  X_train <- as.matrix(X_train)
  center <- colMeans(X_train, na.rm = TRUE)
  scalev <- apply(X_train, 2, sd, na.rm = TRUE)
  scalev[is.na(scalev) | scalev < 1e-8] <- 1
  X_scaled <- sweep(sweep(X_train, 2, center, "-"), 2, scalev, "/")
  X_scaled[is.na(X_scaled)] <- 0
  list(X = X_scaled, center = center, scale = scalev)
}
safe_scale_apply <- function(X, center, scalev) {
  X <- as.matrix(X)
  X_scaled <- sweep(sweep(X, 2, center, "-"), 2, scalev, "/")
  X_scaled[is.na(X_scaled)] <- 0
  X_scaled
}
# ------------------------------------------------------------------------------
# Deterministic censoring calibration
# Fixes one set of uniforms and solves a deterministic censoring objective.
# ------------------------------------------------------------------------------
calibrate_censoring_rate <- function(T_true, target_rate) {
  if (target_rate <= 0) return(rep(Inf, length(T_true)))
  if (target_rate >= 1) stop("target_rate must be < 1.")
  U <- runif(length(T_true))
  f <- function(rate) {
    C_time <- -log(U) / rate
    mean(T_true > C_time) - target_rate
  }
  lower <- 1e-8
  upper <- 1
  while (f(upper) < 0 && upper < 1e8) {
    upper <- upper * 2
  }
  if (upper >= 1e8 && f(upper) < 0) {
    warning("Could not fully bracket censoring-rate root; using upper bound.")
    rate_star <- upper
  } else {
    rate_star <- uniroot(f, lower = lower, upper = upper)$root
  }
  -log(U) / rate_star
}
get_bfdr_threshold <- function(pip_vec, target_bfdr = 0.20) {
  pip_vec <- pip_vec[is.finite(pip_vec)]
  if (length(pip_vec) == 0) return(Inf)
  sorted_pip <- sort(pip_vec, decreasing = TRUE)
  fdr_curve <- cumsum(1 - sorted_pip) / seq_along(sorted_pip)
  valid_idx <- which(fdr_curve <= target_bfdr)
  if (length(valid_idx) == 0) return(Inf)
  sorted_pip[max(valid_idx)]
}
# ------------------------------------------------------------------------------
# C-index helper
# Higher score = longer predicted survival.
# ------------------------------------------------------------------------------
aft_cindex <- function(time, status, lp_survival) {
  survival::concordance(
    survival::Surv(time, status) ~ as.numeric(lp_survival)
  )$concordance
}
# ------------------------------------------------------------------------------
# Sanity checks for simulation design matrices
# ------------------------------------------------------------------------------
validate_sim_design <- function(X, n_fixed, require_intercept = TRUE) {
  X <- as.matrix(X)
  if (require_intercept) {
    if (!("Intercept" %in% colnames(X))) {
      stop("Design matrix does not contain an Intercept column.")
    }
    if (!all(abs(X[, "Intercept"] - 1) < 1e-12)) {
      stop("Intercept column is not all ones.")
    }
  }
  if (!("D" %in% colnames(X))) {
    stop("Design matrix does not contain exposure column D.")
  }
  if (n_fixed < 2) {
    stop("n_fixed should include at least Intercept and D.")
  }
  if (n_fixed > ncol(X)) {
    stop("n_fixed exceeds number of columns in X.")
  }
  if (any(!is.finite(X))) {
    stop("Design matrix contains non-finite values.")
  }
  invisible(TRUE)
}