################################################################################
# Tables_generating_code.R
#
# Replication code for Tables 3, 4, and 5 of:
# "Is There an AI Bubble? Robust Date-Stamping for Periods of Exuberance"
#
# Table 3:
#   Monte Carlo identification-rate comparison between:
#     - SV-ADF procedure: origination cutoff log(ns)/10,
#                         collapse cutoff log(ns)/2
#     - PWY benchmark:    origination/collapse cutoff log(log(ns))/100
#
# Table 4:
#   Identification-rate comparison across volatility designs:
#     - Homoskedastic innovations
#     - Persistent stochastic volatility
#     - GARCH volatility
#
# Table 5:
#   Bias and MSE comparison for estimated origination and collapse dates.
#
# Notes:
# - Package installation is intentionally omitted; see README.md.
# - This script uses only base R.
# - The Monte Carlo sample size R = 1000 is used to match the paper.
# - Set SAVE_TABLES <- TRUE to write CSV outputs.
################################################################################

# ------------------------------------------------------------------------------
# Global settings
# ------------------------------------------------------------------------------

SAVE_TABLES <- TRUE

TABLE_DIR <- "tables"

if (!dir.exists(TABLE_DIR)) {
  dir.create(TABLE_DIR, recursive = TRUE)
}

set.seed(123)

# ------------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------------

TABLE3_CSV <- file.path(TABLE_DIR, "table3_identification_rates_svadf_pwy.csv")
TABLE4_CSV <- file.path(TABLE_DIR, "table4_identification_rates_volatility_conditions.csv")
TABLE5_CSV <- file.path(TABLE_DIR, "table5_bias_mse_svadf_pwy.csv")

################################################################################
# Common helper functions
################################################################################

safe_mean <- function(x) {
  out <- mean(x, na.rm = TRUE)
  if (is.nan(out)) NA_real_ else out
}

round_numeric_columns <- function(df, digits = 4) {
  out <- df
  numeric_cols <- sapply(out, is.numeric)
  out[numeric_cols] <- lapply(out[numeric_cols], round, digits = digits)
  out
}

# ------------------------------------------------------------------------------
# Data-generating process
# ------------------------------------------------------------------------------

# Generate X_0, ..., X_n with one bubble episode.
#
# Model:
#   X_t = a_t X_{t-1} + u_t,
#
# where:
#   a_t = 1 before origination,
#   a_t = delta_n during the bubble,
#   a_t = 1 after collapse,
#   delta_n = 1 + c / n^alpha.
#
# Volatility choices:
#   homosk: u_t = eps_t
#   sv:     log(sigma_t^2) = phi_n log(sigma_{t-1}^2) + eta_t
#   garch:  sigma_t^2 = omega + alpha_g u_{t-1}^2 + beta_g sigma_{t-1}^2
generate_bubble_series <- function(
    n,
    r_e,
    r_f,
    c_par,
    alpha,
    vol_model = c("sv", "homosk", "garch"),
    d_par = 0.1,
    eta_par = 0.1,
    omega_g = 0.01,
    alpha_g = 0.05,
    beta_g = 0.94,
    x0 = 5,
    sigma2_0 = 1
) {
  vol_model <- match.arg(vol_model)
  
  tau_e <- floor(n * r_e)
  tau_f <- floor(n * r_f)
  
  stopifnot(tau_e >= 1, tau_f > tau_e, tau_f <= n)
  
  delta_n <- 1 + c_par / (n^alpha)
  
  X <- numeric(n + 1)
  X[1] <- x0
  
  log_sigma2_prev <- log(sigma2_0)
  sigma2_prev <- sigma2_0
  u_prev <- 0
  
  if (vol_model == "sv") {
    phi_n <- 1 - d_par / log(log(n))
  }
  
  for (t in seq_len(n)) {
    
    if (vol_model == "homosk") {
      sigma_t <- 1
      
    } else if (vol_model == "sv") {
      eta_t <- rnorm(1, mean = 0, sd = eta_par)
      log_sigma2_t <- phi_n * log_sigma2_prev + eta_t
      sigma_t <- sqrt(exp(log_sigma2_t))
      
    } else if (vol_model == "garch") {
      sigma2_t <- omega_g + alpha_g * u_prev^2 + beta_g * sigma2_prev
      sigma_t <- sqrt(sigma2_t)
    }
    
    u_t <- sigma_t * rnorm(1)
    
    a_t <- if (t < tau_e) {
      1
    } else if (t < tau_f) {
      delta_n
    } else {
      1
    }
    
    X[t + 1] <- a_t * X[t] + u_t
    
    if (vol_model == "sv") {
      log_sigma2_prev <- log_sigma2_t
    } else if (vol_model == "garch") {
      sigma2_prev <- sigma2_t
      u_prev <- u_t
    }
  }
  
  X
}

# Modified post-collapse DGP used for PWY collapse calculations.
#
# Up to tau_f, the path follows the same bubble dynamics. After collapse, the
# post-collapse level is reset to X_{tau_e} plus a positive random jump:
#
#   X_{tau_f+1} = X_{tau_e} + U,
#   U ~ Uniform(sqrt(n), 10 sqrt(n)).
#
# The process then continues as a unit root with the same volatility recursion.
generate_pwy_modified_collapse_series <- function(
    n,
    r_e,
    r_f,
    c_par,
    alpha,
    vol_model = c("sv", "homosk", "garch"),
    d_par = 0.1,
    eta_par = 0.1,
    omega_g = 0.01,
    alpha_g = 0.05,
    beta_g = 0.94,
    x0 = 5,
    sigma2_0 = 1
) {
  vol_model <- match.arg(vol_model)
  
  tau_e <- floor(n * r_e)
  tau_f <- floor(n * r_f)
  
  stopifnot(tau_e >= 1, tau_f > tau_e, tau_f < n)
  
  delta_n <- 1 + c_par / (n^alpha)
  
  X <- numeric(n + 1)
  X[1] <- x0
  
  log_sigma2_prev <- log(sigma2_0)
  sigma2_prev <- sigma2_0
  u_prev <- 0
  
  if (vol_model == "sv") {
    phi_n <- 1 - d_par / log(log(n))
  }
  
  # Simulate through the explosive episode.
  for (t in seq_len(tau_f)) {
    
    if (vol_model == "homosk") {
      sigma_t <- 1
      
    } else if (vol_model == "sv") {
      eta_t <- rnorm(1, mean = 0, sd = eta_par)
      log_sigma2_t <- phi_n * log_sigma2_prev + eta_t
      sigma_t <- sqrt(exp(log_sigma2_t))
      
    } else if (vol_model == "garch") {
      sigma2_t <- omega_g + alpha_g * u_prev^2 + beta_g * sigma2_prev
      sigma_t <- sqrt(sigma2_t)
    }
    
    u_t <- sigma_t * rnorm(1)
    
    a_t <- if (t < tau_e) 1 else delta_n
    X[t + 1] <- a_t * X[t] + u_t
    
    if (vol_model == "sv") {
      log_sigma2_prev <- log_sigma2_t
    } else if (vol_model == "garch") {
      sigma2_prev <- sigma2_t
      u_prev <- u_t
    }
  }
  
  # In R indexing, X[tau_f + 2] corresponds to X_{tau_f+1}.
  jump_u <- runif(1, min = sqrt(n), max = 10 * sqrt(n))
  X[tau_f + 2] <- X[tau_e + 1] + jump_u
  
  # Continue as a unit root after the modified collapse.
  if (tau_f + 2 <= n) {
    for (t in (tau_f + 2):n) {
      
      if (vol_model == "homosk") {
        sigma_t <- 1
        
      } else if (vol_model == "sv") {
        eta_t <- rnorm(1, mean = 0, sd = eta_par)
        log_sigma2_t <- phi_n * log_sigma2_prev + eta_t
        sigma_t <- sqrt(exp(log_sigma2_t))
        
      } else if (vol_model == "garch") {
        sigma2_t <- omega_g + alpha_g * u_prev^2 + beta_g * sigma2_prev
        sigma_t <- sqrt(sigma2_t)
      }
      
      u_t <- sigma_t * rnorm(1)
      X[t + 1] <- X[t] + u_t
      
      if (vol_model == "sv") {
        log_sigma2_prev <- log_sigma2_t
      } else if (vol_model == "garch") {
        sigma2_prev <- sigma2_t
        u_prev <- u_t
      }
    }
  }
  
  X
}

# ------------------------------------------------------------------------------
# Recursive ADF statistic and date-stamping rules
# ------------------------------------------------------------------------------

# Compute recursive coefficient-based DF statistic:
#
#   DF_delta(s) = tau * (delta_hat_tau - 1),
#   tau = floor(ns).
adf_for_sgrid <- function(
    X_full,
    n,
    s_seq = seq(0.10, 1.00, by = 0.001)
) {
  nX <- length(X_full) - 1
  stopifnot(nX >= max(floor(n * s_seq)))
  
  csX <- cumsum(X_full)
  csX2 <- cumsum(X_full^2)
  xy <- X_full[2:(n + 1)] * X_full[1:n]
  csXY <- cumsum(xy)
  
  tau_seq <- pmax(2, floor(n * s_seq))
  DF_delta <- rep(NA_real_, length(tau_seq))
  
  for (k in seq_along(tau_seq)) {
    tau <- tau_seq[k]
    
    # Regression: X_t on a constant and X_{t-1}, t = 1, ..., tau.
    S_reg <- csX[tau]
    S_y <- csX[tau + 1] - X_full[1]
    S_reg2 <- csX2[tau]
    S_y2 <- csX2[tau + 1] - X_full[1]^2
    S_cross <- csXY[tau]
    
    Xbar <- S_y / tau
    SSR_t <- S_reg2 - 2 * Xbar * S_reg + tau * Xbar^2
    
    if (!is.finite(SSR_t) || SSR_t <= 0) {
      next
    }
    
    Sxtyt <- S_cross - Xbar * (S_reg + S_y) + tau * Xbar^2
    SSY_t <- S_y2 - 2 * Xbar * S_y + tau * Xbar^2
    
    delta_hat <- Sxtyt / SSR_t
    sigma2_hat <- (SSY_t - 2 * delta_hat * Sxtyt + delta_hat^2 * SSR_t) / tau
    
    if (!is.finite(sigma2_hat) || sigma2_hat <= 0) {
      next
    }
    
    DF_delta[k] <- tau * (delta_hat - 1)
  }
  
  data.frame(
    s = tau_seq / n,
    DF_delta = DF_delta
  )
}

# SV-ADF critical values.
cv_svadf_orig <- function(n, s) log(n * s) / 10
cv_svadf_coll <- function(n, s) log(n * s) / 2

# PWY-style benchmark critical values.
cv_pwy_orig <- function(n, s) log(log(n * s)) / 100
cv_pwy_coll <- function(n, s) log(log(n * s)) / 100

# Generic date-stamping rule.
#
# Origination:
#   first s such that DF_delta(s) > origination cutoff.
#
# Collapse:
#   first s after rhat_e + min_gap such that DF_delta(s) < collapse cutoff.
estimate_breaks_generic <- function(
    X_full,
    n,
    s_seq = seq(0.10, 1.00, by = 0.001),
    min_gap = 0.05,
    cv_orig_fun,
    cv_coll_fun
) {
  adf_df <- adf_for_sgrid(X_full, n, s_seq)
  
  cv_orig <- cv_orig_fun(n, adf_df$s)
  cv_coll <- cv_coll_fun(n, adf_df$s)
  
  idx_e <- which(adf_df$DF_delta > cv_orig)[1]
  rhat_e <- if (!is.na(idx_e)) adf_df$s[idx_e] else NA_real_
  
  idx_f <- NA_integer_
  rhat_f <- NA_real_
  
  if (!is.na(rhat_e)) {
    idx_start_f <- which(adf_df$s >= rhat_e + min_gap)[1]
    
    if (!is.na(idx_start_f)) {
      idx_f_rel <- which(
        adf_df$DF_delta[idx_start_f:nrow(adf_df)] <
          cv_coll[idx_start_f:nrow(adf_df)]
      )[1]
      
      if (!is.na(idx_f_rel)) {
        idx_f <- idx_start_f + idx_f_rel - 1
        rhat_f <- adf_df$s[idx_f]
      }
    }
  }
  
  c(rhat_e = rhat_e, rhat_f = rhat_f)
}

# Estimate PWY collapse on the modified collapse path, while keeping PWY
# origination from the original path.
estimate_pwy_with_modified_collapse <- function(
    X_original,
    X_modified,
    n,
    s_seq = seq(0.10, 1.00, by = 0.001),
    min_gap = 0.05
) {
  pwy_original <- estimate_breaks_generic(
    X_full = X_original,
    n = n,
    s_seq = s_seq,
    min_gap = min_gap,
    cv_orig_fun = cv_pwy_orig,
    cv_coll_fun = cv_pwy_coll
  )
  
  rhat_e <- pwy_original["rhat_e"]
  rhat_f <- NA_real_
  
  if (!is.na(rhat_e)) {
    adf_df_modified <- adf_for_sgrid(X_modified, n, s_seq)
    cv_coll <- cv_pwy_coll(n, adf_df_modified$s)
    
    idx_start_f <- which(adf_df_modified$s >= rhat_e + min_gap)[1]
    
    if (!is.na(idx_start_f)) {
      idx_f_rel <- which(
        adf_df_modified$DF_delta[idx_start_f:nrow(adf_df_modified)] <
          cv_coll[idx_start_f:nrow(adf_df_modified)]
      )[1]
      
      if (!is.na(idx_f_rel)) {
        idx_f <- idx_start_f + idx_f_rel - 1
        rhat_f <- adf_df_modified$s[idx_f]
      }
    }
  }
  
  c(rhat_e = rhat_e, rhat_f = rhat_f)
}

################################################################################
# Monte Carlo engine
################################################################################

run_identification_design <- function(
    n,
    r_e,
    r_f,
    c_par,
    alpha,
    R,
    tol,
    min_gap,
    vol_model = "sv",
    d_par = 0.1,
    eta_par = 0.1,
    omega_g = 0.01,
    alpha_g = 0.05,
    beta_g = 0.94,
    s_seq = seq(0.10, 1.00, by = 0.001),
    use_modified_pwy_collapse = TRUE
) {
  est_svadf <- matrix(NA_real_, nrow = R, ncol = 2)
  est_pwy <- matrix(NA_real_, nrow = R, ncol = 2)
  
  colnames(est_svadf) <- c("rhat_e", "rhat_f")
  colnames(est_pwy) <- c("rhat_e", "rhat_f")
  
  for (rep in seq_len(R)) {
    
    X_original <- generate_bubble_series(
      n = n,
      r_e = r_e,
      r_f = r_f,
      c_par = c_par,
      alpha = alpha,
      vol_model = vol_model,
      d_par = d_par,
      eta_par = eta_par,
      omega_g = omega_g,
      alpha_g = alpha_g,
      beta_g = beta_g,
      x0 = 5
    )
    
    est_svadf[rep, ] <- estimate_breaks_generic(
      X_full = X_original,
      n = n,
      s_seq = s_seq,
      min_gap = min_gap,
      cv_orig_fun = cv_svadf_orig,
      cv_coll_fun = cv_svadf_coll
    )
    
    if (isTRUE(use_modified_pwy_collapse)) {
      X_modified <- generate_pwy_modified_collapse_series(
        n = n,
        r_e = r_e,
        r_f = r_f,
        c_par = c_par,
        alpha = alpha,
        vol_model = vol_model,
        d_par = d_par,
        eta_par = eta_par,
        omega_g = omega_g,
        alpha_g = alpha_g,
        beta_g = beta_g,
        x0 = 5
      )
      
      est_pwy[rep, ] <- estimate_pwy_with_modified_collapse(
        X_original = X_original,
        X_modified = X_modified,
        n = n,
        s_seq = s_seq,
        min_gap = min_gap
      )
      
    } else {
      est_pwy[rep, ] <- estimate_breaks_generic(
        X_full = X_original,
        n = n,
        s_seq = s_seq,
        min_gap = min_gap,
        cv_orig_fun = cv_pwy_orig,
        cv_coll_fun = cv_pwy_coll
      )
    }
  }
  
  orig_id_svadf <- !is.na(est_svadf[, "rhat_e"]) &
    abs(est_svadf[, "rhat_e"] - r_e) < tol
  
  coll_id_svadf <- !is.na(est_svadf[, "rhat_f"]) &
    abs(est_svadf[, "rhat_f"] - r_f) < tol
  
  orig_id_pwy <- !is.na(est_pwy[, "rhat_e"]) &
    abs(est_pwy[, "rhat_e"] - r_e) < tol
  
  coll_id_pwy <- !is.na(est_pwy[, "rhat_f"]) &
    abs(est_pwy[, "rhat_f"] - r_f) < tol
  
  data.frame(
    n = n,
    r_e = r_e,
    r_f = r_f,
    c = c_par,
    alpha = alpha,
    R = R,
    tol = tol,
    vol_model = vol_model,
    d = d_par,
    eta = eta_par,
    omega_g = omega_g,
    alpha_g = alpha_g,
    beta_g = beta_g,
    origination_id_rate_svadf = mean(orig_id_svadf),
    collapse_id_rate_svadf = mean(coll_id_svadf),
    origination_id_rate_pwy = mean(orig_id_pwy),
    collapse_id_rate_pwy = mean(coll_id_pwy),
    stringsAsFactors = FALSE
  )
}

run_bias_mse_design <- function(
    n,
    r_e,
    r_f,
    c_par,
    alpha,
    R,
    min_gap,
    vol_model = "sv",
    d_par = 1,
    eta_par = 0.1,
    omega_g = 0.01,
    alpha_g = 0.05,
    beta_g = 0.94,
    s_seq = seq(0.10, 1.00, by = 0.001),
    use_modified_pwy_collapse = TRUE
) {
  est_svadf <- matrix(NA_real_, nrow = R, ncol = 2)
  est_pwy <- matrix(NA_real_, nrow = R, ncol = 2)
  
  colnames(est_svadf) <- c("rhat_e", "rhat_f")
  colnames(est_pwy) <- c("rhat_e", "rhat_f")
  
  for (rep in seq_len(R)) {
    
    X_original <- generate_bubble_series(
      n = n,
      r_e = r_e,
      r_f = r_f,
      c_par = c_par,
      alpha = alpha,
      vol_model = vol_model,
      d_par = d_par,
      eta_par = eta_par,
      omega_g = omega_g,
      alpha_g = alpha_g,
      beta_g = beta_g,
      x0 = 5
    )
    
    est_svadf[rep, ] <- estimate_breaks_generic(
      X_full = X_original,
      n = n,
      s_seq = s_seq,
      min_gap = min_gap,
      cv_orig_fun = cv_svadf_orig,
      cv_coll_fun = cv_svadf_coll
    )
    
    if (isTRUE(use_modified_pwy_collapse)) {
      X_modified <- generate_pwy_modified_collapse_series(
        n = n,
        r_e = r_e,
        r_f = r_f,
        c_par = c_par,
        alpha = alpha,
        vol_model = vol_model,
        d_par = d_par,
        eta_par = eta_par,
        omega_g = omega_g,
        alpha_g = alpha_g,
        beta_g = beta_g,
        x0 = 5
      )
      
      est_pwy[rep, ] <- estimate_pwy_with_modified_collapse(
        X_original = X_original,
        X_modified = X_modified,
        n = n,
        s_seq = s_seq,
        min_gap = min_gap
      )
      
    } else {
      est_pwy[rep, ] <- estimate_breaks_generic(
        X_full = X_original,
        n = n,
        s_seq = s_seq,
        min_gap = min_gap,
        cv_orig_fun = cv_pwy_orig,
        cv_coll_fun = cv_pwy_coll
      )
    }
  }
  
  err_svadf_e <- est_svadf[, "rhat_e"] - r_e
  err_svadf_f <- est_svadf[, "rhat_f"] - r_f
  
  err_pwy_e <- est_pwy[, "rhat_e"] - r_e
  err_pwy_f <- est_pwy[, "rhat_f"] - r_f
  
  data.frame(
    n = n,
    r_e = r_e,
    r_f = r_f,
    c = c_par,
    alpha = alpha,
    R = R,
    estimated_r_e_svadf = safe_mean(est_svadf[, "rhat_e"]),
    estimated_r_f_svadf = safe_mean(est_svadf[, "rhat_f"]),
    bias_rhat_e_svadf = safe_mean(err_svadf_e),
    mse_rhat_e_svadf = safe_mean(err_svadf_e^2),
    bias_rhat_f_svadf = safe_mean(err_svadf_f),
    mse_rhat_f_svadf = safe_mean(err_svadf_f^2),
    estimated_r_e_pwy = safe_mean(est_pwy[, "rhat_e"]),
    estimated_r_f_pwy = safe_mean(est_pwy[, "rhat_f"]),
    bias_rhat_e_pwy = safe_mean(err_pwy_e),
    mse_rhat_e_pwy = safe_mean(err_pwy_e^2),
    bias_rhat_f_pwy = safe_mean(err_pwy_f),
    mse_rhat_f_pwy = safe_mean(err_pwy_f^2),
    stringsAsFactors = FALSE
  )
}

################################################################################
# Table 3
# Identification-rate comparison: SV-ADF versus PWY
################################################################################

make_table3 <- function() {
  
  n <- 1000
  c_par <- 1
  alpha <- 0.3
  d_par <- 0.1
  eta_par <- 0.1
  R <- 1000
  tol <- 0.1
  min_gap <- 0.1
  
  design_grid <- expand.grid(
    r_e = c(0.20, 0.30, 0.40, 0.50),
    r_f = c(0.50, 0.65, 0.75),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  design_grid <- subset(design_grid, r_f > r_e + min_gap)
  
  results_list <- vector("list", nrow(design_grid))
  
  for (i in seq_len(nrow(design_grid))) {
    message(sprintf("Table 3: design %d of %d", i, nrow(design_grid)))
    
    results_list[[i]] <- run_identification_design(
      n = n,
      r_e = design_grid$r_e[i],
      r_f = design_grid$r_f[i],
      c_par = c_par,
      alpha = alpha,
      R = R,
      tol = tol,
      min_gap = min_gap,
      vol_model = "sv",
      d_par = d_par,
      eta_par = eta_par,
      use_modified_pwy_collapse = TRUE
    )
  }
  
  results <- do.call(rbind, results_list)
  results <- results[order(results$r_e, results$r_f), ]
  rownames(results) <- NULL
  
  results
}

################################################################################
# Table 4
# Identification rates under alternative volatility designs
################################################################################

make_table4 <- function() {
  
  n <- 1000
  alpha <- 0.3
  R <- 1000
  tol <- 0.1
  min_gap <- 0.1
  
  c_grid <- c(0.3, 0.5, 1.0)
  
  design_grid <- expand.grid(
    r_e = c(0.20, 0.30, 0.40),
    r_f = c(0.55, 0.65, 0.75),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  design_grid <- subset(design_grid, r_f > r_e + min_gap)
  
  vol_grid <- data.frame(
    scenario_id = 1:6,
    vol_model = c("homosk", "sv", "sv", "sv", "garch", "garch"),
    vol_label = c(
      "Homoskedastic",
      "SV (d=0.1, eta=0.1)",
      "SV (d=0.1, eta=0.5)",
      "SV (d=0.1, eta=1.0)",
      "GARCH (omega=0.01, alpha_g=0.05, beta_g=0.94)",
      "GARCH (omega=0.01, alpha_g=0.10, beta_g=0.89)"
    ),
    d = c(0.0, 0.1, 0.1, 0.1, NA, NA),
    eta = c(0.0, 0.1, 0.5, 1.0, NA, NA),
    omega_g = c(NA, NA, NA, NA, 0.01, 0.01),
    alpha_g = c(NA, NA, NA, NA, 0.05, 0.10),
    beta_g = c(NA, NA, NA, NA, 0.94, 0.89),
    stringsAsFactors = FALSE
  )
  
  results_list <- list()
  counter <- 1
  
  for (j in seq_len(nrow(vol_grid))) {
    for (cc in c_grid) {
      for (i in seq_len(nrow(design_grid))) {
        
        message(sprintf(
          "Table 4: volatility scenario %d of %d, c = %.1f, design %d of %d",
          j, nrow(vol_grid), cc, i, nrow(design_grid)
        ))
        
        results_list[[counter]] <- run_identification_design(
          n = n,
          r_e = design_grid$r_e[i],
          r_f = design_grid$r_f[i],
          c_par = cc,
          alpha = alpha,
          R = R,
          tol = tol,
          min_gap = min_gap,
          vol_model = vol_grid$vol_model[j],
          d_par = vol_grid$d[j],
          eta_par = vol_grid$eta[j],
          omega_g = vol_grid$omega_g[j],
          alpha_g = vol_grid$alpha_g[j],
          beta_g = vol_grid$beta_g[j],
          use_modified_pwy_collapse = TRUE
        )
        
        results_list[[counter]]$scenario_id <- vol_grid$scenario_id[j]
        results_list[[counter]]$vol_label <- vol_grid$vol_label[j]
        
        counter <- counter + 1
      }
    }
  }
  
  results <- do.call(rbind, results_list)
  
  results <- results[
    order(
      results$scenario_id,
      results$c,
      results$r_e,
      results$r_f
    ),
  ]
  
  rownames(results) <- NULL
  
  results
}

################################################################################
# Table 5
# Bias and MSE comparison
################################################################################

make_table5 <- function() {
  
  n <- 1000
  c_par <- 0.5
  d_par <- 1
  eta_par <- 0.1
  R <- 1000
  min_gap <- 0.05
  
  design_grid <- expand.grid(
    r_e = c(0.20, 0.30, 0.40, 0.50),
    r_f = c(0.50, 0.65, 0.75, 0.85),
    alpha = c(0.30, 0.50, 0.70, 1.00),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  
  design_grid <- subset(design_grid, r_f > r_e + min_gap)
  
  results_list <- vector("list", nrow(design_grid))
  
  for (i in seq_len(nrow(design_grid))) {
    message(sprintf("Table 5: design %d of %d", i, nrow(design_grid)))
    
    results_list[[i]] <- run_bias_mse_design(
      n = n,
      r_e = design_grid$r_e[i],
      r_f = design_grid$r_f[i],
      c_par = c_par,
      alpha = design_grid$alpha[i],
      R = R,
      min_gap = min_gap,
      vol_model = "sv",
      d_par = d_par,
      eta_par = eta_par,
      use_modified_pwy_collapse = TRUE
    )
  }
  
  results <- do.call(rbind, results_list)
  
  results <- results[
    order(
      results$alpha,
      results$r_e,
      results$r_f
    ),
  ]
  
  rownames(results) <- NULL
  
  results
}

################################################################################
# Run all tables
################################################################################

set.seed(123)
table3_results <- make_table3()
table3_print <- round_numeric_columns(table3_results, digits = 4)

print(table3_print)

if (isTRUE(SAVE_TABLES)) {
  write.csv(table3_results, TABLE3_CSV, row.names = FALSE)
}

set.seed(123)
table4_results <- make_table4()
table4_print <- round_numeric_columns(table4_results, digits = 4)

print(table4_print)

if (isTRUE(SAVE_TABLES)) {
  write.csv(table4_results, TABLE4_CSV, row.names = FALSE)
}

set.seed(123)
table5_results <- make_table5()
table5_print <- round_numeric_columns(table5_results, digits = 4)

print(table5_print)

if (isTRUE(SAVE_TABLES)) {
  write.csv(table5_results, TABLE5_CSV, row.names = FALSE)
}

################################################################################
# End of script
################################################################################
