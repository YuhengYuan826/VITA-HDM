# ==============================================================================
# 03_simulation_aft_cavi.R
# Coordinate Ascent Variational Inference engine for log-normal AFT model
# ==============================================================================
fit_aft_cavi <- function(X,
                         Y,
                         Status,
                         n_fixed,
                         priors = list(Lambda = NULL, mu_prior = NULL, tau_sq = 1),
                         max_iter = 200,
                         tol = 1e-5,
                         verbose = FALSE) {
  X <- as.matrix(X)
  Y <- as.numeric(Y)
  Status <- as.numeric(Status)
  X[is.na(X)] <- 0
  Y[is.na(Y)] <- 0
  n <- nrow(X)
  p <- ncol(X)
  if (length(Y) != n) stop("Length of Y does not match nrow(X).")
  if (length(Status) != n) stop("Length of Status does not match nrow(X).")
  if (n_fixed < 0 || n_fixed > p) stop("Invalid n_fixed.")
  Lambda <- priors$Lambda
  if (is.null(Lambda)) {
    Lambda <- rep(log(0.05 / 0.95), p)
  }
  if (length(Lambda) != p) stop("Length of Lambda must equal ncol(X).")
  mu_0 <- priors$mu_prior
  if (is.null(mu_0)) {
    mu_0 <- rep(0, p)
  }
  if (length(mu_0) != p) stop("Length of mu_prior must equal ncol(X).")
  tau_sq <- ifelse(is.null(priors$tau_sq), 1, priors$tau_sq)
  if (!is.finite(tau_sq) || tau_sq <= 0) stop("tau_sq must be positive.")
  alpha <- rep(0.5, p)
  mu <- rep(0, p)
  s_sq <- rep(1, p)
  sigma_sq <- max(var(Y), 1e-4)
  if (n_fixed > 0) {
    alpha[seq_len(n_fixed)] <- 1
    init_fit <- try(
      lm.fit(X[, seq_len(n_fixed), drop = FALSE], Y),
      silent = TRUE
    )
    if (!inherits(init_fit, "try-error")) {
      coef_init <- init_fit$coefficients
      coef_init[is.na(coef_init)] <- 0
      mu[seq_len(n_fixed)] <- coef_init
    }
  }
  Z <- Y
  censored <- which(Status == 0)
  if (length(censored) > 0) {
    Z[censored] <- Y[censored] + sqrt(sigma_sq)
  }
  beta_old <- alpha * mu
  objective_history <- numeric(max_iter)
  diff_val <- Inf
  for (iter in seq_len(max_iter)) {
    beta_curr <- alpha * mu
    lp <- as.numeric(X %*% beta_curr)
    sigma <- sqrt(max(sigma_sq, 1e-8))
    if (length(censored) > 0) {
      a <- (Y[censored] - lp[censored]) / sigma
      log_mills <- dnorm(a, log = TRUE) -
        pnorm(a, lower.tail = FALSE, log.p = TRUE)
      mills <- exp(log_mills)
      bad <- !is.finite(mills)
      if (any(bad)) {
        aa <- abs(a[bad])
        mills[bad] <- aa + 1 / pmax(aa, 1e-8)
      }
      Z[censored] <- lp[censored] + sigma * mills
    }
    residual <- as.numeric(Z - lp)
    for (j in seq_len(p)) {
      xj <- X[, j]
      beta_old_j <- alpha[j] * mu[j]
      r_j <- residual + xj * beta_old_j
      x2 <- sum(xj^2)
      xr <- sum(xj * r_j)
      if (j <= n_fixed) {
        s_sq[j] <- sigma_sq / (x2 + 1e-8)
        mu[j] <- xr / (x2 + 1e-8)
        alpha[j] <- 1
      } else {
        s_sq[j] <- 1 / (x2 / sigma_sq + 1 / tau_sq)
        mu[j] <- s_sq[j] * (xr / sigma_sq + mu_0[j] / tau_sq)
        log_odds <- Lambda[j] +
          0.5 * log(s_sq[j] / tau_sq) +
          (mu[j]^2) / (2 * s_sq[j]) -
          (mu_0[j]^2) / (2 * tau_sq)
        if (!is.finite(log_odds)) log_odds <- -30
        alpha[j] <- plogis(pmax(pmin(log_odds, 30), -30))
      }
      beta_new_j <- alpha[j] * mu[j]
      residual <- r_j - xj * beta_new_j
    }
    beta_mean <- alpha * mu
    beta_var <- alpha * (s_sq + mu^2) - beta_mean^2
    beta_var[beta_var < 0] <- 0
    sse <- sum(residual^2)
    variance_term <- sum(colSums(X^2) * beta_var)
    sigma_sq <- max((sse + variance_term) / n, 1e-6)
    loglik_part <- -0.5 * n * log(sigma_sq) -
      0.5 * (sse + variance_term) / sigma_sq
    pi_prior <- plogis(Lambda)
    prior_part <- 0
    entropy_part <- 0
    if (n_fixed < p) {
      idx <- (n_fixed + 1):p
      prior_part <- sum(
        alpha[idx] *
          (
            log(pi_prior[idx] + 1e-12) -
              0.5 * log(2 * pi * tau_sq) -
              (s_sq[idx] + (mu[idx] - mu_0[idx])^2) / (2 * tau_sq)
          ) +
          (1 - alpha[idx]) * log(1 - pi_prior[idx] + 1e-12)
      )
      entropy_part <- sum(
        -alpha[idx] * log(alpha[idx] + 1e-12) -
          (1 - alpha[idx]) * log(1 - alpha[idx] + 1e-12) +
          alpha[idx] * 0.5 * log(2 * pi * exp(1) * s_sq[idx])
      )
    }
    objective_history[iter] <- loglik_part + prior_part + entropy_part
    diff_val <- max(abs(beta_mean - beta_old), na.rm = TRUE)
    if (diff_val < tol) break
    beta_old <- beta_mean
  }
  if (verbose) {
    message(sprintf("CAVI finished at iter %d; diff = %.6g", iter, diff_val))
  }
  list(
    beta_hat = alpha * mu,
    PIP = alpha,
    mu = mu,
    s_sq = s_sq,
    sigma_sq = sigma_sq,
    objective_history = objective_history[seq_len(iter)]
  )
}