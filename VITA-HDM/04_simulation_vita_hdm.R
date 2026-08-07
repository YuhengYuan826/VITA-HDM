# ==============================================================================
# 04_simulation_vita_hdm.R
# Simulation wrapper for VITA-HDM with explicit Intercept
# ==============================================================================
prepare_simulation_data <- function(target_train,
                                    source_list = NULL,
                                    target_test = NULL,
                                    center_y = FALSE,
                                    run_sanity_checks = TRUE) {
  raw_design_cols <- c(
    "D",
    grep("^X", names(target_train), value = TRUE),
    grep("^M", names(target_train), value = TRUE)
  )
  m_cols <- grep("^M", raw_design_cols, value = TRUE)
  x_cols <- grep("^X", raw_design_cols, value = TRUE)
  # Current simulation convention:
  #   Intercept is explicit and unscaled.
  #   D, X, and M are standardized using target-training scaling.
  scale_obj <- safe_scale_train(as.matrix(target_train[, raw_design_cols]))
  scaled_train <- scale_obj$X
  colnames(scaled_train) <- raw_design_cols
  intercept_train <- matrix(1, nrow = nrow(target_train), ncol = 1)
  colnames(intercept_train) <- "Intercept"
  train_X <- cbind(intercept_train, scaled_train)
  y_mean <- if (center_y) mean(target_train$Y, na.rm = TRUE) else 0
  train_Y <- target_train$Y - y_mean
  n_fixed <- 1 + 1 + length(x_cols)  # Intercept + D + X covariates
  out <- list(
    train = list(
      X = train_X,
      Y = train_Y,
      Status = target_train$Status,
      Time = target_train$Time
    ),
    n_fixed = n_fixed,
    design_cols = colnames(train_X),
    raw_scaled_cols = raw_design_cols,
    mediator_cols = m_cols,
    covariate_cols = c("Intercept", "D", x_cols),
    scale_center = scale_obj$center,
    scale_sd = scale_obj$scale,
    Y_mean = y_mean
  )
  if (run_sanity_checks) {
    validate_sim_design(out$train$X, out$n_fixed)
  }
  if (!is.null(target_test)) {
    scaled_test <- safe_scale_apply(
      as.matrix(target_test[, raw_design_cols]),
      scale_obj$center,
      scale_obj$scale
    )
    colnames(scaled_test) <- raw_design_cols
    intercept_test <- matrix(1, nrow = nrow(target_test), ncol = 1)
    colnames(intercept_test) <- "Intercept"
    test_X <- cbind(intercept_test, scaled_test)
    out$test <- list(
      X = test_X,
      Y = target_test$Y - y_mean,
      Status = target_test$Status,
      Time = target_test$Time
    )
    if (run_sanity_checks) {
      validate_sim_design(out$test$X, out$n_fixed)
    }
  }
  if (!is.null(source_list)) {
    out$sources <- lapply(source_list, function(src) {
      scaled_src <- safe_scale_apply(
        as.matrix(src[, raw_design_cols]),
        scale_obj$center,
        scale_obj$scale
      )
      colnames(scaled_src) <- raw_design_cols
      intercept_src <- matrix(1, nrow = nrow(src), ncol = 1)
      colnames(intercept_src) <- "Intercept"
      src_X <- cbind(intercept_src, scaled_src)
      src_out <- list(
        X = src_X,
        Y = src$Y - y_mean,
        Status = src$Status,
        Time = src$Time
      )
      if (run_sanity_checks) {
        validate_sim_design(src_out$X, out$n_fixed)
      }
      src_out
    })
  } else {
    out$sources <- list()
  }
  out
}
synthesize_priors_sim <- function(prep,
                                  eta0 = log(0.05 / 0.95),
                                  eta1 = 1,
                                  clamp_pip = 0.90,
                                  tau_sq = 1,
                                  transfer_beta_strength = 1.0) {
  K <- length(prep$sources)
  p <- ncol(prep$train$X)
  if (K == 0) {
    return(list(
      Lambda = c(rep(10, prep$n_fixed), rep(eta0, p - prep$n_fixed)),
      mu_prior = rep(0, p),
      tau_sq = tau_sq,
      weights = numeric(0),
      source_lpml = numeric(0)
    ))
  }
  fits <- vector("list", K)
  lpml <- numeric(K)
  for (k in seq_len(K)) {
    src <- prep$sources[[k]]
    pri <- list(
      Lambda = c(rep(10, prep$n_fixed), rep(log(0.05 / 0.95), p - prep$n_fixed)),
      mu_prior = rep(0, p),
      tau_sq = tau_sq
    )
    fit_k <- fit_aft_cavi(
      X = src$X,
      Y = src$Y,
      Status = src$Status,
      n_fixed = prep$n_fixed,
      priors = pri,
      max_iter = 200,
      verbose = FALSE
    )
    fits[[k]] <- fit_k
    pred <- as.numeric(prep$train$X %*% fit_k$beta_hat)
    sig <- sqrt(fit_k$sigma_sq)
    z <- (prep$train$Y - pred) / sig
    log_pdf <- dnorm(z, log = TRUE) - log(sig)
    log_surv <- pnorm(z, lower.tail = FALSE, log.p = TRUE)
    lpml[k] <- sum(
      prep$train$Status * log_pdf +
        (1 - prep$train$Status) * log_surv,
      na.rm = TRUE
    )
  }
  w <- exp(lpml - max(lpml))
  w <- w / sum(w)
  logit_min <- qlogis(1 - clamp_pip)
  logit_max <- qlogis(clamp_pip)
  pip_logit_bar <- rep(0, p)
  beta_bar <- rep(0, p)
  for (k in seq_len(K)) {
    pip_k <- pmin(pmax(fits[[k]]$PIP, 1e-6), 1 - 1e-6)
    logit_k <- pmin(pmax(qlogis(pip_k), logit_min), logit_max)
    pip_logit_bar <- pip_logit_bar + w[k] * logit_k
    beta_bar <- beta_bar + w[k] * fits[[k]]$beta_hat
  }
  Lambda <- eta0 + eta1 * pip_logit_bar
  Lambda[seq_len(prep$n_fixed)] <- 10
  mu_prior <- transfer_beta_strength * beta_bar
  mu_prior[seq_len(prep$n_fixed)] <- 0
  list(
    Lambda = Lambda,
    mu_prior = mu_prior,
    tau_sq = tau_sq,
    weights = w,
    source_lpml = lpml
  )
}
VITA_HDM_sim <- function(target_train,
                         source_list = NULL,
                         target_test = NULL,
                         transfer = TRUE,
                         fdr_target = 0.20,
                         eta0 = log(0.05 / 0.95),
                         eta1 = 1,
                         clamp_pip = 0.90,
                         tau_sq = 1,
                         transfer_beta_strength = 1.0,
                         center_y = FALSE,
                         run_sanity_checks = TRUE,
                         verbose = FALSE) {
  prep <- prepare_simulation_data(
    target_train = target_train,
    source_list = source_list,
    target_test = target_test,
    center_y = center_y,
    run_sanity_checks = run_sanity_checks
  )
  p <- ncol(prep$train$X)
  if (transfer && length(prep$sources) > 0) {
    priors <- synthesize_priors_sim(
      prep = prep,
      eta0 = eta0,
      eta1 = eta1,
      clamp_pip = clamp_pip,
      tau_sq = tau_sq,
      transfer_beta_strength = transfer_beta_strength
    )
  } else {
    priors <- list(
      Lambda = c(rep(10, prep$n_fixed), rep(eta0, p - prep$n_fixed)),
      mu_prior = rep(0, p),
      tau_sq = tau_sq,
      weights = numeric(0),
      source_lpml = numeric(0)
    )
  }
  fit <- fit_aft_cavi(
    X = prep$train$X,
    Y = prep$train$Y,
    Status = prep$train$Status,
    n_fixed = prep$n_fixed,
    priors = priors,
    max_iter = 200,
    verbose = verbose
  )
  vars <- colnames(prep$train$X)
  inf <- data.frame(
    Variable = vars,
    Effect = fit$beta_hat,
    PIP = fit$PIP,
    Mean = fit$mu,
    SD = sqrt(fit$s_sq),
    stringsAsFactors = FALSE
  )
  med_tab <- inf[inf$Variable %in% prep$mediator_cols, ]
  th <- get_bfdr_threshold(med_tab$PIP, fdr_target)
  selected <- med_tab$Variable[med_tab$PIP >= th]
  cindex <- NA_real_
  if (!is.null(prep$test)) {
    lp_test <- as.numeric(prep$test$X %*% fit$beta_hat)
    cindex <- aft_cindex(prep$test$Time, prep$test$Status, lp_test)
  }
  list(
    fit = fit,
    inference = inf,
    selected_mediators = selected,
    threshold = th,
    prep = prep,
    priors = priors,
    test_cindex = cindex
  )
}