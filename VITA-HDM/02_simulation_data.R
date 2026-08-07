# ==============================================================================
# 02_simulation_data.R
# Data-generation functions for VITA-HDM simulation study
# ==============================================================================
generate_truth_params <- function(p_x = 5, p_m = 1000, n_active = 4, seed = 1) {
  set.seed(seed)
  alpha1 <- rep(0, p_m)
  beta2  <- rep(0, p_m)
  alpha1[1:n_active] <- runif(n_active, 0.4, 0.6)
  beta2[1:n_active]  <- runif(n_active, 0.4, 0.6)
  list(
    gamma = rep(0.3, p_x),
    alpha1 = alpha1,
    alpha2 = matrix(0.1, nrow = p_x, ncol = p_m),
    beta0 = 0.1,
    beta1 = 0.5,
    beta2 = beta2,
    beta4 = rep(0.1, p_x),
    active_idx = 1:n_active
  )
}
make_source_params <- function(params,
                               h = 0,
                               magnitude_shift = 0.15,
                               heterogeneous = FALSE,
                               seed = 1) {
  out <- params
  if (!heterogeneous) return(out)
  set.seed(seed)
  nz_a <- which(out$alpha1 != 0)
  nz_b <- which(out$beta2 != 0)
  out$alpha1[nz_a] <- out$alpha1[nz_a] +
    runif(length(nz_a), -magnitude_shift, magnitude_shift)
  out$beta2[nz_b] <- out$beta2[nz_b] +
    runif(length(nz_b), -magnitude_shift, magnitude_shift)
  if (h > 0) {
    zeros <- which(params$alpha1 == 0 & params$beta2 == 0)
    trap <- sample(zeros, h)
    out$alpha1[trap] <- runif(h, 0.4, 0.6)
    out$beta2[trap]  <- runif(h, 0.4, 0.6)
  }
  out
}
# ------------------------------------------------------------------------------
# Moderate platform-shift perturbation
# This simulates source-domain transcriptomic measurement distortion:
#   1. gene-wise mean shift
#   2. gene-wise scale distortion
#   3. additional measurement noise
# Important:
#   Survival is generated from the latent biological mediator matrix M_bio.
#   The fitted model observes M_obs, which may be platform-distorted.
# ------------------------------------------------------------------------------
apply_platform_shift <- function(M,
                                 frac = 0.35,
                                 mean_shift_sd = 0.50,
                                 scale_log_sd = 0.30,
                                 noise_sd = 0.20,
                                 seed = 1) {
  set.seed(seed)
  M <- as.matrix(M)
  n <- nrow(M)
  p <- ncol(M)
  if (p == 0) return(M)
  affected <- sample(seq_len(p), size = max(1, floor(frac * p)))
  gene_shift <- rnorm(length(affected), mean = 0, sd = mean_shift_sd)
  gene_scale <- exp(rnorm(length(affected), mean = 0, sd = scale_log_sd))
  M[, affected] <- sweep(
    M[, affected, drop = FALSE],
    2,
    gene_shift,
    "+"
  )
  M[, affected] <- sweep(
    M[, affected, drop = FALSE],
    2,
    gene_scale,
    "*"
  )
  if (noise_sd > 0) {
    M <- M + matrix(
      rnorm(n * p, mean = 0, sd = noise_sd),
      nrow = n,
      ncol = p
    )
  }
  M
}
generate_aft_dataset <- function(n,
                                 params,
                                 rho = 0,
                                 censoring_rate = 0.2,
                                 source_shift = FALSE,
                                 platform_shift = FALSE,
                                 platform_frac = 0.35,
                                 platform_mean_sd = 0.50,
                                 platform_scale_sd = 0.30,
                                 platform_noise_sd = 0.20,
                                 seed = 1) {
  set.seed(seed)
  p_x <- length(params$gamma)
  p_m <- length(params$alpha1)
  Sigma_X <- make_ar1_cov(p_x, 0.5)
  if (source_shift) {
    A <- matrix(
      sample(c(0, 0.2), p_x * p_x, replace = TRUE, prob = c(0.8, 0.2)),
      p_x,
      p_x
    )
    Sigma_X <- Sigma_X + t(A) %*% A
  }
  X <- MASS::mvrnorm(n, rep(0, p_x), Sigma_X)
  D <- as.numeric(X %*% params$gamma + rnorm(n))
  Sigma_M <- make_ar1_cov(p_m, rho)
  E_M <- MASS::mvrnorm(n, rep(0, p_m), Sigma_M)
  M_bio <- matrix(0, n, p_m)
  for (j in seq_len(p_m)) {
    M_bio[, j] <- params$alpha1[j] * D +
      X %*% params$alpha2[, j] +
      E_M[, j]
  }
  Z_true <- as.numeric(
    params$beta0 +
      params$beta1 * D +
      M_bio %*% params$beta2 +
      X %*% params$beta4 +
      rnorm(n)
  )
  T_true <- exp(Z_true)
  C_time <- calibrate_censoring_rate(T_true, censoring_rate)
  Time <- pmin(T_true, C_time)
  Status <- as.numeric(T_true <= C_time)
  M_obs <- M_bio
  if (platform_shift) {
    M_obs <- apply_platform_shift(
      M_obs,
      frac = platform_frac,
      mean_shift_sd = platform_mean_sd,
      scale_log_sd = platform_scale_sd,
      noise_sd = platform_noise_sd,
      seed = seed + 7000
    )
  }
  df <- data.frame(
    Time = Time,
    Status = Status,
    Y = log(Time),
    D = D,
    M_obs,
    X
  )
  colnames(df) <- c(
    "Time",
    "Status",
    "Y",
    "D",
    paste0("M", seq_len(p_m)),
    paste0("X", seq_len(p_x))
  )
  list(
    data = df,
    params = params,
    actual_censoring = mean(Status == 0)
  )
}