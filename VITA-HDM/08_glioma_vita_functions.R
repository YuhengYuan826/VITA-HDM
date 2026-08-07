# ==============================================================================
# 08_glioma_vita_functions.R
# Core utility and model functions for the final cleaned glioma VITA-HDM example.
# This script defines functions only and does not execute the analysis.
# No sink() logging is used by default, so errors remain visible in the R console.
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(survival)
})
# ==============================================================================
# Basic utilities
# ==============================================================================
as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
safe_plogis <- function(x) plogis(pmax(pmin(x, 30), -30))
safe_logit <- function(p) {
  p <- pmin(pmax(p, 1e-8), 1 - 1e-8)
  qlogis(p)
}
safe_scale_train <- function(X) {
  X <- as.matrix(X)
  center <- colMeans(X, na.rm = TRUE)
  scalev <- apply(X, 2, sd, na.rm = TRUE)
  scalev[is.na(scalev) | scalev < 1e-8] <- 1
  Xs <- sweep(sweep(X, 2, center, "-"), 2, scalev, "/")
  Xs[is.na(Xs)] <- 0
  list(X = Xs, center = center, scale = scalev)
}
safe_scale_apply <- function(X, center, scalev) {
  X <- as.matrix(X)
  Xs <- sweep(sweep(X, 2, center, "-"), 2, scalev, "/")
  Xs[is.na(Xs)] <- 0
  Xs
}
aft_cindex <- function(time, status, lp_survival) {
  survival::concordance(
    survival::Surv(time, status) ~ as.numeric(lp_survival)
  )$concordance
}
make_stratified_folds <- function(status, k = 5, seed = 1) {
  set.seed(seed)
  status <- as.numeric(status)
  idx_event <- sample(which(status == 1))
  idx_cens  <- sample(which(status == 0))
  folds <- vector("list", k)
  for (i in seq_len(k)) folds[[i]] <- integer(0)
  for (i in seq_along(idx_event)) {
    folds[[((i - 1) %% k) + 1]] <- c(folds[[((i - 1) %% k) + 1]], idx_event[i])
  }
  for (i in seq_along(idx_cens)) {
    folds[[((i - 1) %% k) + 1]] <- c(folds[[((i - 1) %% k) + 1]], idx_cens[i])
  }
  lapply(folds, sort)
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
get_selected_genes <- function(pip_gene, target_bfdr = 0.20) {
  pip_gene <- pip_gene[is.finite(pip_gene)]
  if (length(pip_gene) == 0) {
    return(list(selected = character(0), threshold = Inf, bfdr = NA_real_))
  }
  th <- get_bfdr_threshold(pip_gene, target_bfdr)
  selected <- names(pip_gene)[pip_gene >= th]
  if (length(selected) == 0 || !is.finite(th)) {
    return(list(selected = character(0), threshold = th, bfdr = NA_real_))
  }
  bfdr <- mean(1 - pip_gene[selected], na.rm = TRUE)
  list(selected = selected, threshold = th, bfdr = bfdr)
}
safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  mean(x)
}
safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) return(NA_real_)
  sd(x)
}
safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(quantile(x, probs = prob, na.rm = TRUE, names = FALSE))
}
safe_ratio <- function(num, den) {
  out <- rep(NA_real_, length(num))
  ok <- is.finite(num) & is.finite(den) & abs(den) > 1e-8
  out[ok] <- num[ok] / den[ok]
  out
}
sum_finite_or_na <- function(x) {
  if (all(!is.finite(x))) return(NA_real_)
  sum(x, na.rm = TRUE)
}
# ==============================================================================
# Design construction
# ==============================================================================
get_complete_cohort <- function(aligned_obj,
                                cohort_name = c("target", "source"),
                                covariates = c("Age", "Gender")) {
  cohort_name <- match.arg(cohort_name)
  obj <- aligned_obj[[cohort_name]]
  if (is.null(obj)) stop("Cannot find cohort in aligned object: ", cohort_name)
  clin <- obj$clinical
  expr <- obj$expr
  clin$Sample_ID <- as.character(clin$Sample_ID)
  colnames(expr) <- as.character(colnames(expr))
  common_samples <- intersect(clin$Sample_ID, colnames(expr))
  clin <- clin[match(common_samples, clin$Sample_ID), , drop = FALSE]
  expr <- expr[, common_samples, drop = FALSE]
  req <- c("OS_time", "OS_status", "IDH", covariates)
  missing_req <- setdiff(req, colnames(clin))
  if (length(missing_req) > 0) {
    stop(cohort_name, " is missing required clinical columns: ",
         paste(missing_req, collapse = ", "))
  }
  keep <- !is.na(clin$OS_time) & clin$OS_time > 0 &
    !is.na(clin$OS_status) & clin$OS_status %in% c(0, 1) &
    !is.na(clin$IDH) & clin$IDH %in% c(0, 1)
  for (cc in covariates) keep <- keep & !is.na(clin[[cc]])
  clin <- clin[keep, , drop = FALSE]
  expr <- expr[, clin$Sample_ID, drop = FALSE]
  list(clinical = clin, expr = expr, cohort = obj$cohort)
}
select_top_variable_genes <- function(target_cohort,
                                      n_genes_keep = 3000,
                                      gene_universe = NULL) {
  expr <- target_cohort$expr
  if (!is.null(gene_universe)) {
    gene_universe <- intersect(gene_universe, rownames(expr))
    expr <- expr[gene_universe, , drop = FALSE]
  }
  vars <- apply(expr, 1, var, na.rm = TRUE)
  vars <- vars[is.finite(vars)]
  if (length(vars) == 0) stop("No finite gene variances in target cohort.")
  names(sort(vars, decreasing = TRUE))[seq_len(min(n_genes_keep, length(vars)))]
}
validate_design <- function(design) {
  X <- as.matrix(design$X)
  if (!("Intercept" %in% colnames(X))) stop("Missing Intercept.")
  if (!all(abs(X[, "Intercept"] - 1) < 1e-12)) stop("Intercept is not all ones.")
  if (!("D" %in% colnames(X))) stop("Missing D exposure column.")
  dvals <- unique(as.numeric(X[, "D"]))
  if (!all(dvals %in% c(0, 1))) stop("D is not binary 0/1.")
  if (any(!is.finite(X))) stop("X contains non-finite values.")
  if (any(!is.finite(design$Y))) stop("Y contains non-finite values.")
  if (length(design$Y) != nrow(X)) stop("Y length mismatch.")
  if (length(design$Status) != nrow(X)) stop("Status length mismatch.")
  invisible(TRUE)
}
prepare_design <- function(cohort,
                           keep_genes,
                           covariates = c("Age", "Gender"),
                           binary_covariates = c("Gender"),
                           cov_scaler = NULL,
                           gene_scaler = NULL,
                           y_mean = NULL,
                           fit_scalers = TRUE) {
  clin <- cohort$clinical
  expr <- cohort$expr
  keep_genes <- intersect(keep_genes, rownames(expr))
  if (length(keep_genes) == 0) stop("No keep_genes available in expression matrix.")
  clin <- clin[match(colnames(expr), clin$Sample_ID), , drop = FALSE]
  if (!all(clin$Sample_ID == colnames(expr))) stop("Clinical/expression order mismatch.")
  n <- nrow(clin)
  intercept <- matrix(1, nrow = n, ncol = 1)
  colnames(intercept) <- "Intercept"
  D <- matrix(as.numeric(clin$IDH), ncol = 1)
  colnames(D) <- "D"
  cov_parts <- list()
  cov_scaler_out <- cov_scaler
  if (is.null(cov_scaler_out)) cov_scaler_out <- list()
  for (cc in covariates) {
    v <- as.numeric(clin[[cc]])
    if (cc %in% binary_covariates) {
      mat <- matrix(v, ncol = 1)
      colnames(mat) <- cc
      cov_parts[[cc]] <- mat
    } else {
      raw <- matrix(v, ncol = 1)
      colnames(raw) <- cc
      if (fit_scalers || is.null(cov_scaler) || is.null(cov_scaler[[cc]])) {
        sc <- safe_scale_train(raw)
        cov_scaler_out[[cc]] <- list(center = sc$center, scale = sc$scale)
        mat <- sc$X
      } else {
        mat <- safe_scale_apply(raw, cov_scaler[[cc]]$center, cov_scaler[[cc]]$scale)
      }
      colnames(mat) <- cc
      cov_parts[[cc]] <- mat
    }
  }
  cov_mat <- if (length(cov_parts) > 0) do.call(cbind, cov_parts) else NULL
  gene_raw <- t(expr[keep_genes, , drop = FALSE])
  colnames(gene_raw) <- paste0("M_", keep_genes)
  if (fit_scalers || is.null(gene_scaler)) {
    gene_sc <- safe_scale_train(gene_raw)
    gene_scaler_out <- list(center = gene_sc$center, scale = gene_sc$scale)
    gene_mat <- gene_sc$X
  } else {
    gene_mat <- safe_scale_apply(gene_raw, gene_scaler$center, gene_scaler$scale)
    gene_scaler_out <- gene_scaler
  }
  colnames(gene_mat) <- paste0("M_", keep_genes)
  X <- cbind(intercept, D, cov_mat, gene_mat)
  if (is.null(y_mean)) y_mean <- mean(log(clin$OS_time), na.rm = TRUE)
  Y <- log(clin$OS_time) - y_mean
  n_fixed <- 2 + length(covariates)
  out <- list(
    X = X,
    Y = Y,
    Status = as.numeric(clin$OS_status),
    Time = as.numeric(clin$OS_time),
    clinical = clin,
    GeneNames = colnames(gene_mat),
    RawGeneNames = keep_genes,
    CovNames = colnames(X)[seq_len(n_fixed)],
    n_fixed = n_fixed,
    y_mean = y_mean,
    cov_scaler = cov_scaler_out,
    gene_scaler = gene_scaler_out
  )
  validate_design(out)
  out
}
prepare_target_train_test <- function(target_cohort,
                                      train_idx,
                                      test_idx,
                                      keep_genes,
                                      covariates = c("Age", "Gender"),
                                      binary_covariates = c("Gender")) {
  train_cohort <- list(
    clinical = target_cohort$clinical[train_idx, , drop = FALSE],
    expr = target_cohort$expr[, target_cohort$clinical$Sample_ID[train_idx], drop = FALSE],
    cohort = paste0(target_cohort$cohort, "_train")
  )
  test_cohort <- list(
    clinical = target_cohort$clinical[test_idx, , drop = FALSE],
    expr = target_cohort$expr[, target_cohort$clinical$Sample_ID[test_idx], drop = FALSE],
    cohort = paste0(target_cohort$cohort, "_test")
  )
  train_design <- prepare_design(
    cohort = train_cohort,
    keep_genes = keep_genes,
    covariates = covariates,
    binary_covariates = binary_covariates,
    fit_scalers = TRUE
  )
  test_design <- prepare_design(
    cohort = test_cohort,
    keep_genes = keep_genes,
    covariates = covariates,
    binary_covariates = binary_covariates,
    cov_scaler = train_design$cov_scaler,
    gene_scaler = train_design$gene_scaler,
    y_mean = train_design$y_mean,
    fit_scalers = FALSE
  )
  list(train = train_design, test = test_design)
}
prepare_source_design <- function(source_cohort,
                                  keep_genes,
                                  covariates = c("Age", "Gender"),
                                  binary_covariates = c("Gender")) {
  prepare_design(
    cohort = source_cohort,
    keep_genes = keep_genes,
    covariates = covariates,
    binary_covariates = binary_covariates,
    fit_scalers = TRUE
  )
}
make_bootstrap_cohort <- function(cohort, boot_idx, prefix = "boot") {
  clin <- cohort$clinical[boot_idx, , drop = FALSE]
  old_ids <- as.character(clin$Sample_ID)
  new_ids <- paste0(prefix, "_", seq_along(old_ids), "_", old_ids)
  clin$Sample_ID <- new_ids
  clin$Patient_ID <- new_ids
  expr <- cohort$expr[, old_ids, drop = FALSE]
  colnames(expr) <- new_ids
  list(
    clinical = clin,
    expr = expr,
    cohort = paste0(cohort$cohort, "_", prefix)
  )
}
# ==============================================================================
# CAVI engine for log-normal AFT model
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
  if (is.null(Lambda)) Lambda <- rep(log(0.05 / 0.95), p)
  if (length(Lambda) != p) stop("Length of Lambda must equal ncol(X).")
  mu_0 <- priors$mu_prior
  if (is.null(mu_0)) mu_0 <- rep(0, p)
  if (length(mu_0) != p) stop("Length of mu_prior must equal ncol(X).")
  tau_sq <- ifelse(is.null(priors$tau_sq), 1, priors$tau_sq)
  if (!is.finite(tau_sq) || tau_sq <= 0) stop("tau_sq must be positive.")
  alpha <- rep(0.5, p)
  mu <- rep(0, p)
  s_sq <- rep(1, p)
  sigma_sq <- max(var(Y), 1e-4)
  if (n_fixed > 0) {
    alpha[seq_len(n_fixed)] <- 1
    init_fit <- try(lm.fit(X[, seq_len(n_fixed), drop = FALSE], Y), silent = TRUE)
    if (!inherits(init_fit, "try-error")) {
      coef_init <- init_fit$coefficients
      coef_init[is.na(coef_init)] <- 0
      mu[seq_len(n_fixed)] <- coef_init
    }
  }
  Z <- Y
  censored <- which(Status == 0)
  if (length(censored) > 0) Z[censored] <- Y[censored] + sqrt(sigma_sq)
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
        alpha[j] <- safe_plogis(log_odds)
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
    objective_history = objective_history[seq_len(iter)],
    iter = iter,
    diff = diff_val
  )
}
# ==============================================================================
# Prior synthesis and model wrappers
# ==============================================================================
make_base_priors <- function(p,
                             n_fixed,
                             eta0 = 0.10,
                             tau_sq = 0.50) {
  if (!is.finite(eta0) || eta0 <= 0 || eta0 >= 1) stop("eta0 must be between 0 and 1.")
  if (!is.finite(tau_sq) || tau_sq <= 0) stop("tau_sq must be positive.")
  if (n_fixed < 0 || n_fixed > p) stop("Invalid n_fixed.")
  list(
    Lambda = c(rep(10, n_fixed), rep(log(eta0 / (1 - eta0)), p - n_fixed)),
    mu_prior = rep(0, p),
    tau_sq = tau_sq,
    source_info = data.frame(),
    source_pip_bar = numeric(0),
    source_logit_bar = numeric(0),
    source_beta_bar = numeric(0)
  )
}
extract_gene_vectors <- function(fit, design) {
  vars <- colnames(design$X)
  n_fixed <- design$n_fixed
  gene_vars <- vars[(n_fixed + 1):length(vars)]
  
  beta_gene <- fit$beta_hat[match(gene_vars, vars)]
  pip_gene <- fit$PIP[match(gene_vars, vars)]
  names(beta_gene) <- gene_vars
  names(pip_gene) <- gene_vars
  
  list(beta_gene = beta_gene, pip_gene = pip_gene)
}
fit_source_model <- function(source_design,
                             source_name = "source",
                             source_eta0 = 0.20,
                             tau_sq = 0.50,
                             max_iter = 100,
                             verbose = FALSE) {
  p <- ncol(source_design$X)
  pri <- make_base_priors(
    p = p,
    n_fixed = source_design$n_fixed,
    eta0 = source_eta0,
    tau_sq = tau_sq
  )
  fit <- fit_aft_cavi(
    X = source_design$X,
    Y = source_design$Y,
    Status = source_design$Status,
    n_fixed = source_design$n_fixed,
    priors = pri,
    max_iter = max_iter,
    verbose = verbose
  )
  gv <- extract_gene_vectors(fit, source_design)
  list(
    source_name = source_name,
    fit = fit,
    beta_gene = gv$beta_gene,
    pip_gene = gv$pip_gene,
    n = nrow(source_design$X),
    events = sum(source_design$Status == 1, na.rm = TRUE)
  )
}
compute_aft_lpml <- function(design, beta_hat, sigma_sq) {
  pred <- as.numeric(design$X %*% beta_hat)
  sig <- sqrt(max(sigma_sq, 1e-8))
  z <- (design$Y - pred) / sig
  log_pdf <- dnorm(z, log = TRUE) - log(sig)
  log_surv <- pnorm(z, lower.tail = FALSE, log.p = TRUE)
  sum(
    design$Status * log_pdf +
      (1 - design$Status) * log_surv,
    na.rm = TRUE
  )
}
synthesize_single_source_priors <- function(target_design,
                                            source_design,
                                            source_fit,
                                            eta0 = 0.10,
                                            source_eta0 = 0.20,
                                            strength = 1.50,
                                            tau_sq = 0.50,
                                            transfer_beta_strength = 0.25,
                                            clamp_pip = 0.90,
                                            weight_mode = "mean_lpml_temperature",
                                            lpml_temperature = 5) {
  p <- ncol(target_design$X)
  n_fixed <- target_design$n_fixed
  base <- make_base_priors(
    p = p,
    n_fixed = n_fixed,
    eta0 = eta0,
    tau_sq = tau_sq
  )
  all_vars <- colnames(target_design$X)
  gene_vars <- all_vars[(n_fixed + 1):length(all_vars)]
  # Single-source transfer: LPML is still reported, but softmax weight degenerates to 1.
  source_lpml <- compute_aft_lpml(
    design = target_design,
    beta_hat = source_fit$fit$beta_hat,
    sigma_sq = source_fit$fit$sigma_sq
  )
  source_pip <- source_fit$pip_gene[gene_vars]
  source_beta <- source_fit$beta_gene[gene_vars]
  source_pip[!is.finite(source_pip)] <- 0
  source_beta[!is.finite(source_beta)] <- 0
  logit_min <- qlogis(1 - clamp_pip)
  logit_max <- qlogis(clamp_pip)
  source_logit <- pmin(pmax(safe_logit(source_pip), logit_min), logit_max)
  gene_idx <- match(gene_vars, all_vars)
  base$Lambda[gene_idx] <- log(eta0 / (1 - eta0)) + strength * source_logit
  base$mu_prior[gene_idx] <- transfer_beta_strength * source_beta
  base$source_info <- data.frame(
    Source = source_fit$source_name,
    LPML = source_lpml,
    Weight = 1,
    Weight_Mode = weight_mode,
    LPML_Temperature = lpml_temperature,
    stringsAsFactors = FALSE
  )
  base$source_pip_bar <- source_pip
  base$source_logit_bar <- source_logit
  base$source_beta_bar <- source_beta
  names(base$source_pip_bar) <- gene_vars
  names(base$source_logit_bar) <- gene_vars
  names(base$source_beta_bar) <- gene_vars
  base
}
fit_target_with_transfer <- function(target_design,
                                     source_design = NULL,
                                     source_name = "source",
                                     source_fit = NULL,
                                     transfer = TRUE,
                                     eta0 = 0.10,
                                     source_eta0 = 0.20,
                                     strength = 1.50,
                                     tau_sq = 0.50,
                                     transfer_beta_strength = 0.25,
                                     weight_mode = "mean_lpml_temperature",
                                     lpml_temperature = 5,
                                     max_iter = 100,
                                     verbose = FALSE) {
  p <- ncol(target_design$X)
  if (transfer) {
    if (is.null(source_design) || is.null(source_fit)) {
      stop("source_design and source_fit must be provided when transfer=TRUE.")
    }
    pri <- synthesize_single_source_priors(
      target_design = target_design,
      source_design = source_design,
      source_fit = source_fit,
      eta0 = eta0,
      source_eta0 = source_eta0,
      strength = strength,
      tau_sq = tau_sq,
      transfer_beta_strength = transfer_beta_strength,
      weight_mode = weight_mode,
      lpml_temperature = lpml_temperature
    )
  } else {
    pri <- make_base_priors(
      p = p,
      n_fixed = target_design$n_fixed,
      eta0 = eta0,
      tau_sq = tau_sq
    )
  }
  fit <- fit_aft_cavi(
    X = target_design$X,
    Y = target_design$Y,
    Status = target_design$Status,
    n_fixed = target_design$n_fixed,
    priors = pri,
    max_iter = max_iter,
    verbose = verbose
  )
  gv <- extract_gene_vectors(fit, target_design)
  fixed_vars <- colnames(target_design$X)[seq_len(target_design$n_fixed)]
  fixed_beta <- fit$beta_hat[seq_len(target_design$n_fixed)]
  names(fixed_beta) <- fixed_vars
  list(
    fit = fit,
    priors = pri,
    beta_gene = gv$beta_gene,
    pip_gene = gv$pip_gene,
    fixed_beta = fixed_beta,
    beta_all = fit$beta_hat,
    target_design = target_design
  )
}
predict_cindex_aft <- function(fit_obj, test_design) {
  lp <- as.numeric(test_design$X %*% fit_obj$fit$beta_hat)
  aft_cindex(test_design$Time, test_design$Status, lp)
}
# ==============================================================================
# Cox helper functions for fixed-signature predictive assessment
# ==============================================================================
fit_cox_cindex <- function(train_df, test_df, vars) {
  vars <- intersect(vars, colnames(train_df))
  if (length(vars) == 0) return(NA_real_)
  rhs <- paste(sprintf("`%s`", vars), collapse = " + ")
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
  lp <- try(as.numeric(predict(fit, newdata = test_df, type = "lp")), silent = TRUE)
  if (inherits(lp, "try-error") || any(!is.finite(lp))) return(NA_real_)
  # Cox lp is higher hazard, so use -lp as a survival score.
  aft_cindex(test_df$Time, test_df$Status, -lp)
}
make_cox_fold_data_target <- function(target_cohort,
                                      idx_train,
                                      idx_test,
                                      raw_genes,
                                      covariates = c("Age", "Gender"),
                                      binary_covariates = c("Gender")) {
  clin <- target_cohort$clinical
  raw_genes <- intersect(raw_genes, rownames(target_cohort$expr))
  clin_train <- clin[idx_train, , drop = FALSE]
  clin_test <- clin[idx_test, , drop = FALSE]
  train_df <- data.frame(
    Time = clin_train$OS_time,
    Status = clin_train$OS_status,
    D = clin_train$IDH,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  test_df <- data.frame(
    Time = clin_test$OS_time,
    Status = clin_test$OS_status,
    D = clin_test$IDH,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (cc in covariates) {
    if (cc %in% binary_covariates) {
      train_df[[cc]] <- as.numeric(clin_train[[cc]])
      test_df[[cc]] <- as.numeric(clin_test[[cc]])
    } else {
      sc <- safe_scale_train(matrix(clin_train[[cc]], ncol = 1))
      train_df[[cc]] <- as.numeric(sc$X[, 1])
      test_df[[cc]] <- as.numeric(
        safe_scale_apply(matrix(clin_test[[cc]], ncol = 1), sc$center, sc$scale)[, 1]
      )
    }
  }
  if (length(raw_genes) > 0) {
    expr_train_raw <- t(target_cohort$expr[raw_genes, clin_train$Sample_ID, drop = FALSE])
    expr_test_raw  <- t(target_cohort$expr[raw_genes, clin_test$Sample_ID, drop = FALSE])
    colnames(expr_train_raw) <- paste0("M_", raw_genes)
    colnames(expr_test_raw) <- paste0("M_", raw_genes)
    gene_sc <- safe_scale_train(expr_train_raw)
    expr_train <- gene_sc$X
    expr_test <- safe_scale_apply(expr_test_raw, gene_sc$center, gene_sc$scale)
    colnames(expr_test) <- colnames(expr_train)
    train_df <- cbind(train_df, expr_train)
    test_df <- cbind(test_df, expr_test)
  }
  list(train = train_df, test = test_df)
}
make_cox_source_train_target_test <- function(source_cohort,
                                              target_test_clin,
                                              target_expr,
                                              raw_genes,
                                              covariates = c("Age", "Gender"),
                                              binary_covariates = c("Gender")) {
  raw_genes <- intersect(raw_genes, intersect(rownames(source_cohort$expr), rownames(target_expr)))
  source_clin <- source_cohort$clinical
  train_df <- data.frame(
    Time = source_clin$OS_time,
    Status = source_clin$OS_status,
    D = source_clin$IDH,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  test_df <- data.frame(
    Time = target_test_clin$OS_time,
    Status = target_test_clin$OS_status,
    D = target_test_clin$IDH,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (cc in covariates) {
    if (cc %in% binary_covariates) {
      train_df[[cc]] <- as.numeric(source_clin[[cc]])
      test_df[[cc]] <- as.numeric(target_test_clin[[cc]])
    } else {
      sc <- safe_scale_train(matrix(source_clin[[cc]], ncol = 1))
      train_df[[cc]] <- as.numeric(sc$X[, 1])
      test_df[[cc]] <- as.numeric(
        safe_scale_apply(matrix(target_test_clin[[cc]], ncol = 1), sc$center, sc$scale)[, 1]
      )
    }
  }
  if (length(raw_genes) > 0) {
    expr_source_raw <- t(source_cohort$expr[raw_genes, source_clin$Sample_ID, drop = FALSE])
    expr_target_raw <- t(target_expr[raw_genes, target_test_clin$Sample_ID, drop = FALSE])
    colnames(expr_source_raw) <- paste0("M_", raw_genes)
    colnames(expr_target_raw) <- paste0("M_", raw_genes)
    gene_sc <- safe_scale_train(expr_source_raw)
    expr_source <- gene_sc$X
    expr_target <- safe_scale_apply(expr_target_raw, gene_sc$center, gene_sc$scale)
    colnames(expr_target) <- colnames(expr_source)
    train_df <- cbind(train_df, expr_source)
    test_df <- cbind(test_df, expr_target)
  }
  list(train = train_df, test = test_df)
}
# ==============================================================================
# Mediation-effect helpers
# ==============================================================================
compute_alpha_paths_safe <- function(target_design,
                                     gene_names,
                                     covariates = c("Age", "Gender")) {
  X <- target_design$X
  gene_names <- intersect(gene_names, colnames(X))
  alpha <- rep(NA_real_, length(gene_names))
  names(alpha) <- gene_names
  rhs_vars <- c("D", covariates)
  rhs_vars <- intersect(rhs_vars, colnames(X))
  if (!("D" %in% rhs_vars)) return(alpha)
  rhs <- paste(sprintf("`%s`", rhs_vars), collapse = " + ")
  for (g in gene_names) {
    dat <- as.data.frame(X[, rhs_vars, drop = FALSE], check.names = FALSE)
    dat$M <- X[, g]
    fit <- try(lm(as.formula(paste("M ~", rhs)), data = dat), silent = TRUE)
    if (!inherits(fit, "try-error")) {
      co <- coef(fit)
      alpha[g] <- ifelse("D" %in% names(co), as.numeric(co["D"]), NA_real_)
    }
  }
  alpha
}