#Simulation Power comparison with PWY (updated DGP for PWY)

n       <- 1000
c_par   <- 1
alpha   <- 0.3
d_par   <- 0.1
eta_par <- 0.1
R       <- 1000
tol     <- 0.1
min_gap <- 0.1

design_grid <- expand.grid(
  r_e = c(0.20, 0.30, 0.40, 0.50),
  r_f = c(0.50, 0.65, 0.75),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

design_grid <- subset(design_grid, r_f > r_e + min_gap)

# Original DGP
generate_series <- function(n, r_e, r_f, c, alpha, d, eta,
                            x0 = 5, sigma2_0 = 1) {
  tau_e <- floor(n * r_e)
  tau_f <- floor(n * r_f)
  stopifnot(tau_e >= 1, tau_f > tau_e, tau_f <= n)
  
  delta_n <- 1 + c / (n^alpha)
  phi_n   <- 1 - d / log(log(n))
  
  X <- numeric(n + 1)
  X[1] <- x0
  
  logsig2_prev <- log(sigma2_0)
  
  for (t in 1:n) {
    eta_t     <- rnorm(1, 0, eta)
    logsig2_t <- phi_n * logsig2_prev + eta_t
    sigma_t   <- sqrt(exp(logsig2_t))
    u_t       <- sigma_t * rnorm(1)
    
    a_t <- if (t < tau_e) {
      1
    } else if (t < tau_f) {
      delta_n
    } else {
      1
    }
    
    X[t + 1] <- a_t * X[t] + u_t
    logsig2_prev <- logsig2_t
  }
  
  X
}

# Modified DGP used only for Algorithm 2 collapse:
# up to tau_f same as original,
# then X_{tau_f+1} = X_{tau_e+1} + U, U ~ Unif(sqrt(n), 10 sqrt(n)),
# then continue as unit root with same SV volatility recursion.
generate_series_alg2_collapse <- function(n, r_e, r_f, c, alpha, d, eta,
                                          x0 = 5, sigma2_0 = 1) {
  tau_e <- floor(n * r_e)
  tau_f <- floor(n * r_f)
  stopifnot(tau_e >= 1, tau_f > tau_e, tau_f < n)
  
  delta_n <- 1 + c / (n^alpha)
  phi_n   <- 1 - d / log(log(n))
  
  X <- numeric(n + 1)
  X[1] <- x0
  
  logsig2_prev <- log(sigma2_0)
  
  # Simulate up to X_{tau_f}
  for (t in 1:tau_f) {
    eta_t     <- rnorm(1, 0, eta)
    logsig2_t <- phi_n * logsig2_prev + eta_t
    sigma_t   <- sqrt(exp(logsig2_t))
    u_t       <- sigma_t * rnorm(1)
    
    a_t <- if (t < tau_e) 1 else delta_n
    X[t + 1] <- a_t * X[t] + u_t
    logsig2_prev <- logsig2_t
  }
  
  # Set X_{tau_f + 1} = X_{tau_e} + U
  jump_u <- runif(1, min = sqrt(n), max = 10 * sqrt(n))
  X[tau_f + 2] <- X[tau_e + 1] + jump_u
  
  # Continue as unit root from time tau_f + 2 onward
  if (tau_f + 2 <= n) {
    for (t in (tau_f + 2):n) {
      eta_t     <- rnorm(1, 0, eta)
      logsig2_t <- phi_n * logsig2_prev + eta_t
      sigma_t   <- sqrt(exp(logsig2_t))
      u_t       <- sigma_t * rnorm(1)
      
      X[t + 1] <- X[t] + u_t
      logsig2_prev <- logsig2_t
    }
  }
  
  X
}

adf_for_sgrid <- function(X_full, n, s_seq) {
  nX <- length(X_full) - 1
  stopifnot(nX >= max(floor(n * s_seq)))
  
  csX  <- cumsum(X_full)
  csX2 <- cumsum(X_full^2)
  xy   <- X_full[2:(n + 1)] * X_full[1:n]
  csXY <- cumsum(xy)
  
  tau_seq  <- pmax(2, floor(n * s_seq))
  DF_delta <- rep(NA_real_, length(tau_seq))
  
  for (k in seq_along(tau_seq)) {
    tau <- tau_seq[k]
    
    S_reg   <- csX[tau]
    S_y     <- csX[tau + 1] - X_full[1]
    S_reg2  <- csX2[tau]
    S_y2    <- csX2[tau + 1] - X_full[1]^2
    S_cross <- csXY[tau]
    
    Xbar  <- S_y / tau
    SSR_t <- S_reg2 - 2 * Xbar * S_reg + tau * Xbar^2
    if (!is.finite(SSR_t) || SSR_t <= 0) next
    
    Sxtyt <- S_cross - Xbar * (S_reg + S_y) + tau * Xbar^2
    SSY_t <- S_y2    - 2 * Xbar * S_y       + tau * Xbar^2
    
    delta_hat  <- Sxtyt / SSR_t
    sigma2_hat <- (SSY_t - 2 * delta_hat * Sxtyt + (delta_hat^2) * SSR_t) / tau
    if (!is.finite(sigma2_hat) || sigma2_hat <= 0) next
    
    DF_delta[k] <- tau * (delta_hat - 1)
  }
  
  data.frame(s = tau_seq / n, DF_delta = DF_delta)
}

estimate_breaks_generic <- function(X_full, n,
                                    s_seq = seq(0.10, 1.00, by = 0.001),
                                    min_gap = 0.05,
                                    cv_orig_fun,
                                    cv_coll_fun) {
  adf_df <- adf_for_sgrid(X_full, n, s_seq)
  
  cv_orig <- cv_orig_fun(n, adf_df$s)
  cv_coll <- cv_coll_fun(n, adf_df$s)
  
  idx_e  <- which(adf_df$DF_delta > cv_orig)[1]
  rhat_e <- if (!is.na(idx_e)) adf_df$s[idx_e] else NA_real_
  
  idx_f  <- NA_integer_
  rhat_f <- NA_real_
  
  if (!is.na(rhat_e)) {
    idx_start_f <- which(adf_df$s >= (rhat_e + min_gap))[1]
    
    if (!is.na(idx_start_f)) {
      idx_f_rel <- which(adf_df$DF_delta[idx_start_f:nrow(adf_df)] <
                           cv_coll[idx_start_f:nrow(adf_df)])[1]
      
      if (!is.na(idx_f_rel)) {
        idx_f  <- idx_start_f + idx_f_rel - 1
        rhat_f <- adf_df$s[idx_f]
      }
    }
  }
  
  c(rhat_e = rhat_e, rhat_f = rhat_f)
}

cv1_orig_fun <- function(n, s) log(n * s) / 10
cv1_coll_fun <- function(n, s) log(n * s) / 2

cv2_orig_fun <- function(n, s) log(log(n * s)) / 100
cv2_coll_fun <- function(n, s) log(log(n * s)) / 100

run_design_power <- function(n, r_e, r_f, c_par, alpha, d_par, eta_par, R,
                             tol = 0.05,
                             s_seq = seq(0.10, 1.00, by = 0.001),
                             min_gap = 0.05) {
  est_alg1 <- matrix(NA_real_, nrow = R, ncol = 2)
  est_alg2 <- matrix(NA_real_, nrow = R, ncol = 2)
  colnames(est_alg1) <- c("rhat_e", "rhat_f")
  colnames(est_alg2) <- c("rhat_e", "rhat_f")
  
  for (rep in 1:R) {
    # Original path: used for Algorithm 1 and Algorithm 2 origination
    X_full_orig <- generate_series(
      n = n, r_e = r_e, r_f = r_f,
      c = c_par, alpha = alpha, d = d_par, eta = eta_par,
      x0 = 5
    )
    
    est_alg1[rep, ] <- estimate_breaks_generic(
      X_full      = X_full_orig,
      n           = n,
      s_seq       = s_seq,
      min_gap     = min_gap,
      cv_orig_fun = cv1_orig_fun,
      cv_coll_fun = cv1_coll_fun
    )
    
    # For Algorithm 2:
    # keep origination estimate based on the original path
    est_alg2_orig <- estimate_breaks_generic(
      X_full      = X_full_orig,
      n           = n,
      s_seq       = s_seq,
      min_gap     = min_gap,
      cv_orig_fun = cv2_orig_fun,
      cv_coll_fun = cv2_coll_fun
    )
    
    # Modified path: used only for Algorithm 2 collapse
    X_full_mod <- generate_series_alg2_collapse(
      n = n, r_e = r_e, r_f = r_f,
      c = c_par, alpha = alpha, d = d_par, eta = eta_par,
      x0 = 5
    )
    
    adf_df_mod <- adf_for_sgrid(X_full_mod, n, s_seq)
    cv_coll2   <- cv2_coll_fun(n, adf_df_mod$s)
    
    rhat_e_alg2 <- est_alg2_orig["rhat_e"]
    rhat_f_alg2 <- NA_real_
    
    if (!is.na(rhat_e_alg2)) {
      idx_start_f <- which(adf_df_mod$s >= (rhat_e_alg2 + min_gap))[1]
      if (!is.na(idx_start_f)) {
        idx_f_rel <- which(adf_df_mod$DF_delta[idx_start_f:nrow(adf_df_mod)] <
                             cv_coll2[idx_start_f:nrow(adf_df_mod)])[1]
        if (!is.na(idx_f_rel)) {
          idx_f_alg2 <- idx_start_f + idx_f_rel - 1
          rhat_f_alg2 <- adf_df_mod$s[idx_f_alg2]
        }
      }
    }
    
    est_alg2[rep, "rhat_e"] <- rhat_e_alg2
    est_alg2[rep, "rhat_f"] <- rhat_f_alg2
  }
  
  orig_id_alg1 <- !is.na(est_alg1[, "rhat_e"]) & (abs(est_alg1[, "rhat_e"] - r_e) < tol)
  coll_id_alg1 <- !is.na(est_alg1[, "rhat_f"]) & (abs(est_alg1[, "rhat_f"] - r_f) < tol)
  
  orig_id_alg2 <- !is.na(est_alg2[, "rhat_e"]) & (abs(est_alg2[, "rhat_e"] - r_e) < tol)
  coll_id_alg2 <- !is.na(est_alg2[, "rhat_f"]) & (abs(est_alg2[, "rhat_f"] - r_f) < tol)
  
  data.frame(
    n   = n,
    r_e = r_e,
    r_f = r_f,
    alpha = alpha,
    R   = R,
    tol = tol,
    origination_id_rate_alg1 = mean(orig_id_alg1),
    collapse_id_rate_alg1    = mean(coll_id_alg1),
    origination_id_rate_alg2 = mean(orig_id_alg2),
    collapse_id_rate_alg2    = mean(coll_id_alg2),
    stringsAsFactors = FALSE
  )
}

set.seed(123)

results_list <- vector("list", nrow(design_grid))

for (i in seq_len(nrow(design_grid))) {
  results_list[[i]] <- run_design_power(
    n       = n,
    r_e     = design_grid$r_e[i],
    r_f     = design_grid$r_f[i],
    c_par   = c_par,
    alpha   = alpha,
    d_par   = d_par,
    eta_par = eta_par,
    R       = R,
    tol     = tol,
    min_gap = min_gap
  )
}

results_table <- do.call(rbind, results_list)
results_table <- results_table[order(results_table$r_e, results_table$r_f), ]
row.names(results_table) <- NULL

results_table_print <- results_table
num_cols <- sapply(results_table_print, is.numeric)
results_table_print[num_cols] <- lapply(
  results_table_print[num_cols],
  function(x) round(x, 4)
)

print(results_table_print)
###################################################################################

# Different volatility models considered


#Comparison with other volatility processes (GARCH)
# ------------------------------------------------------------
# Power table under different volatility conditions
#
# Algorithm 1:
#   original DGP only
#
# Algorithm 2:
#   origination on original DGP
#   collapse on modified post-tau_f jump DGP
#
# Reports the same four columns:
#   origination_id_rate_alg1
#   collapse_id_rate_alg1
#   origination_id_rate_alg2
#   collapse_id_rate_alg2
# ------------------------------------------------------------

# -----------------------------
# Fixed parameters
# -----------------------------
n       <- 1000
alpha   <- 0.3
R       <- 1000
tol     <- 0.1
min_gap <- 0.1

c_grid <- c(0.3, 0.5, 1.0)

design_grid <- expand.grid(
  r_e = c(0.20, 0.30, 0.40),
  r_f = c(0.55, 0.65, 0.75),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
design_grid <- subset(design_grid, r_f > r_e + min_gap)

# -----------------------------
# Volatility scenarios
# -----------------------------
include_garch <- TRUE   # set TRUE if you also want GARCH rows

vol_grid <- data.frame(
  scenario_id = 1:4,
  vol_model   = c("homosk", "sv", "sv", "sv"),
  vol_label   = c("Homoskedastic (d=0, eta=0)",
                  "SV (d=0.1, eta=0.1)",
                  "SV (d=0.1, eta=0.5)",
                  "SV (d=0.1, eta=1.0)"),
  d           = c(0.0, 0.1, 0.1, 0.1),
  eta         = c(0.0, 0.1, 0.5, 1.0),
  omega_g     = NA_real_,
  alpha_g     = NA_real_,
  beta_g      = NA_real_,
  stringsAsFactors = FALSE
)

if (include_garch) {
  garch_grid <- data.frame(
    scenario_id = 5:6,
    vol_model   = c("garch", "garch"),
    vol_label   = c("GARCH (omega=0.01, a=0.05, b=0.94)",
                    "GARCH (omega=0.01, a=0.10, b=0.89)"),
    d           = NA_real_,
    eta         = NA_real_,
    omega_g     = c(0.01, 0.01),
    alpha_g     = c(0.05, 0.10),
    beta_g      = c(0.94, 0.89),
    stringsAsFactors = FALSE
  )
  vol_grid <- rbind(vol_grid, garch_grid)
}

# -----------------------------
# Original DGP
# -----------------------------
generate_series <- function(n, r_e, r_f, c, alpha,
                            vol_model = c("homosk", "sv", "garch"),
                            d = 0.1, eta = 0.1,
                            omega_g = 0.01, alpha_g = 0.05, beta_g = 0.94,
                            x0 = 5, sigma2_0 = 1) {
  vol_model <- match.arg(vol_model)
  
  tau_e <- floor(n * r_e)
  tau_f <- floor(n * r_f)
  stopifnot(tau_e >= 1, tau_f > tau_e, tau_f <= n)
  
  delta_n <- 1 + c / (n^alpha)
  
  X <- numeric(n + 1)
  X[1] <- x0
  
  logsig2_prev <- log(sigma2_0)
  sigma2_prev  <- sigma2_0
  u_prev       <- 0
  
  if (vol_model == "sv") {
    phi_n <- 1 - d / log(log(n))
  }
  
  for (t in 1:n) {
    if (vol_model == "homosk") {
      sigma_t <- 1
      
    } else if (vol_model == "sv") {
      eta_t     <- rnorm(1, 0, eta)
      logsig2_t <- phi_n * logsig2_prev + eta_t
      sigma_t   <- sqrt(exp(logsig2_t))
      
    } else if (vol_model == "garch") {
      sigma2_t <- omega_g + alpha_g * (u_prev^2) + beta_g * sigma2_prev
      sigma_t  <- sqrt(sigma2_t)
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
      logsig2_prev <- logsig2_t
    } else if (vol_model == "garch") {
      sigma2_prev <- sigma2_t
      u_prev      <- u_t
    }
  }
  
  X
}

# -----------------------------
# Modified DGP used only for Algorithm 2 collapse
# Up to tau_f: same as the user's updated code
# Then X_{tau_f+1} = X_{tau_e} + U, U ~ Unif(sqrt(n), 10 sqrt(n))
# Then continue as unit root with the same volatility recursion
# -----------------------------
generate_series_alg2_collapse <- function(n, r_e, r_f, c, alpha,
                                          vol_model = c("homosk", "sv", "garch"),
                                          d = 0.1, eta = 0.1,
                                          omega_g = 0.01, alpha_g = 0.05, beta_g = 0.94,
                                          x0 = 5, sigma2_0 = 1) {
  vol_model <- match.arg(vol_model)
  
  tau_e <- floor(n * r_e)
  tau_f <- floor(n * r_f)
  stopifnot(tau_e >= 1, tau_f > tau_e, tau_f < n)
  
  delta_n <- 1 + c / (n^alpha)
  
  X <- numeric(n + 1)
  X[1] <- x0
  
  logsig2_prev <- log(sigma2_0)
  sigma2_prev  <- sigma2_0
  u_prev       <- 0
  
  if (vol_model == "sv") {
    phi_n <- 1 - d / log(log(n))
  }
  
  # Simulate through time tau_f using the updated pre-jump rule
  for (t in 1:tau_f) {
    if (vol_model == "homosk") {
      sigma_t <- 1
      
    } else if (vol_model == "sv") {
      eta_t     <- rnorm(1, 0, eta)
      logsig2_t <- phi_n * logsig2_prev + eta_t
      sigma_t   <- sqrt(exp(logsig2_t))
      
    } else if (vol_model == "garch") {
      sigma2_t <- omega_g + alpha_g * (u_prev^2) + beta_g * sigma2_prev
      sigma_t  <- sqrt(sigma2_t)
    }
    
    u_t <- sigma_t * rnorm(1)
    
    a_t <- if (t < tau_e) 1 else delta_n
    X[t + 1] <- a_t * X[t] + u_t
    
    if (vol_model == "sv") {
      logsig2_prev <- logsig2_t
    } else if (vol_model == "garch") {
      sigma2_prev <- sigma2_t
      u_prev      <- u_t
    }
  }
  
  # Set X_{tau_f + 1} = X_{tau_e} + U
  jump_u <- runif(1, min = sqrt(n), max = 10 * sqrt(n))
  X[tau_f + 2] <- X[tau_e + 1] + jump_u
  
  # Continue as unit root from time tau_f + 2 onward
  if (tau_f + 2 <= n) {
    for (t in (tau_f + 2):n) {
      if (vol_model == "homosk") {
        sigma_t <- 1
        
      } else if (vol_model == "sv") {
        eta_t     <- rnorm(1, 0, eta)
        logsig2_t <- phi_n * logsig2_prev + eta_t
        sigma_t   <- sqrt(exp(logsig2_t))
        
      } else if (vol_model == "garch") {
        sigma2_t <- omega_g + alpha_g * (u_prev^2) + beta_g * sigma2_prev
        sigma_t  <- sqrt(sigma2_t)
      }
      
      u_t <- sigma_t * rnorm(1)
      X[t + 1] <- X[t] + u_t
      
      if (vol_model == "sv") {
        logsig2_prev <- logsig2_t
      } else if (vol_model == "garch") {
        sigma2_prev <- sigma2_t
        u_prev      <- u_t
      }
    }
  }
  
  X
}

# -----------------------------
# Recursive ADF path
# -----------------------------
adf_for_sgrid <- function(X_full, n, s_seq) {
  nX <- length(X_full) - 1
  stopifnot(nX >= max(floor(n * s_seq)))
  
  csX  <- cumsum(X_full)
  csX2 <- cumsum(X_full^2)
  xy   <- X_full[2:(n + 1)] * X_full[1:n]
  csXY <- cumsum(xy)
  
  tau_seq  <- pmax(2, floor(n * s_seq))
  DF_delta <- rep(NA_real_, length(tau_seq))
  
  for (k in seq_along(tau_seq)) {
    tau <- tau_seq[k]
    
    S_reg   <- csX[tau]
    S_y     <- csX[tau + 1] - X_full[1]
    S_reg2  <- csX2[tau]
    S_y2    <- csX2[tau + 1] - X_full[1]^2
    S_cross <- csXY[tau]
    
    Xbar  <- S_y / tau
    SSR_t <- S_reg2 - 2 * Xbar * S_reg + tau * Xbar^2
    if (!is.finite(SSR_t) || SSR_t <= 0) next
    
    Sxtyt <- S_cross - Xbar * (S_reg + S_y) + tau * Xbar^2
    SSY_t <- S_y2    - 2 * Xbar * S_y       + tau * Xbar^2
    
    delta_hat  <- Sxtyt / SSR_t
    sigma2_hat <- (SSY_t - 2 * delta_hat * Sxtyt + (delta_hat^2) * SSR_t) / tau
    if (!is.finite(sigma2_hat) || sigma2_hat <= 0) next
    
    DF_delta[k] <- tau * (delta_hat - 1)
  }
  
  data.frame(s = tau_seq / n, DF_delta = DF_delta)
}

# -----------------------------
# Generic break estimator
# -----------------------------
estimate_breaks_generic <- function(X_full, n,
                                    s_seq = seq(0.10, 1.00, by = 0.001),
                                    min_gap = 0.05,
                                    cv_orig_fun,
                                    cv_coll_fun) {
  adf_df <- adf_for_sgrid(X_full, n, s_seq)
  
  cv_orig <- cv_orig_fun(n, adf_df$s)
  cv_coll <- cv_coll_fun(n, adf_df$s)
  
  idx_e  <- which(adf_df$DF_delta > cv_orig)[1]
  rhat_e <- if (!is.na(idx_e)) adf_df$s[idx_e] else NA_real_
  
  idx_f  <- NA_integer_
  rhat_f <- NA_real_
  
  if (!is.na(rhat_e)) {
    idx_start_f <- which(adf_df$s >= (rhat_e + min_gap))[1]
    
    if (!is.na(idx_start_f)) {
      idx_f_rel <- which(adf_df$DF_delta[idx_start_f:nrow(adf_df)] <
                           cv_coll[idx_start_f:nrow(adf_df)])[1]
      
      if (!is.na(idx_f_rel)) {
        idx_f  <- idx_start_f + idx_f_rel - 1
        rhat_f <- adf_df$s[idx_f]
      }
    }
  }
  
  c(rhat_e = rhat_e, rhat_f = rhat_f)
}

# -----------------------------
# Cutoffs
# -----------------------------
cv1_orig_fun <- function(n, s) log(n * s) / 10
cv1_coll_fun <- function(n, s) log(n * s) / 2

cv2_orig_fun <- function(n, s) log(log(n * s)) / 100
cv2_coll_fun <- function(n, s) log(log(n * s)) / 100

# -----------------------------
# One design x one volatility scenario x one c
# -----------------------------
run_design_power <- function(n, r_e, r_f, c_par, alpha, R,
                             vol_model, d, eta, omega_g, alpha_g, beta_g,
                             tol = 0.05,
                             s_seq = seq(0.10, 1.00, by = 0.001),
                             min_gap = 0.05) {
  est_alg1 <- matrix(NA_real_, nrow = R, ncol = 2)
  est_alg2 <- matrix(NA_real_, nrow = R, ncol = 2)
  colnames(est_alg1) <- c("rhat_e", "rhat_f")
  colnames(est_alg2) <- c("rhat_e", "rhat_f")
  
  for (rep in 1:R) {
    # Original path: used for Algorithm 1 and Algorithm 2 origination
    X_full_orig <- generate_series(
      n         = n,
      r_e       = r_e,
      r_f       = r_f,
      c         = c_par,
      alpha     = alpha,
      vol_model = vol_model,
      d         = d,
      eta       = eta,
      omega_g   = omega_g,
      alpha_g   = alpha_g,
      beta_g    = beta_g,
      x0        = 5
    )
    
    est_alg1[rep, ] <- estimate_breaks_generic(
      X_full      = X_full_orig,
      n           = n,
      s_seq       = s_seq,
      min_gap     = min_gap,
      cv_orig_fun = cv1_orig_fun,
      cv_coll_fun = cv1_coll_fun
    )
    
    # Algorithm 2 origination on original path
    est_alg2_orig <- estimate_breaks_generic(
      X_full      = X_full_orig,
      n           = n,
      s_seq       = s_seq,
      min_gap     = min_gap,
      cv_orig_fun = cv2_orig_fun,
      cv_coll_fun = cv2_coll_fun
    )
    
    # Algorithm 2 collapse on modified post-tau_f jump path
    X_full_mod <- generate_series_alg2_collapse(
      n         = n,
      r_e       = r_e,
      r_f       = r_f,
      c         = c_par,
      alpha     = alpha,
      vol_model = vol_model,
      d         = d,
      eta       = eta,
      omega_g   = omega_g,
      alpha_g   = alpha_g,
      beta_g    = beta_g,
      x0        = 5
    )
    
    adf_df_mod <- adf_for_sgrid(X_full_mod, n, s_seq)
    cv_coll2   <- cv2_coll_fun(n, adf_df_mod$s)
    
    rhat_e_alg2 <- est_alg2_orig["rhat_e"]
    rhat_f_alg2 <- NA_real_
    
    if (!is.na(rhat_e_alg2)) {
      idx_start_f <- which(adf_df_mod$s >= (rhat_e_alg2 + min_gap))[1]
      
      if (!is.na(idx_start_f)) {
        idx_f_rel <- which(adf_df_mod$DF_delta[idx_start_f:nrow(adf_df_mod)] <
                             cv_coll2[idx_start_f:nrow(adf_df_mod)])[1]
        
        if (!is.na(idx_f_rel)) {
          idx_f_alg2 <- idx_start_f + idx_f_rel - 1
          rhat_f_alg2 <- adf_df_mod$s[idx_f_alg2]
        }
      }
    }
    
    est_alg2[rep, "rhat_e"] <- rhat_e_alg2
    est_alg2[rep, "rhat_f"] <- rhat_f_alg2
  }
  
  orig_id_alg1 <- !is.na(est_alg1[, "rhat_e"]) & (abs(est_alg1[, "rhat_e"] - r_e) < tol)
  coll_id_alg1 <- !is.na(est_alg1[, "rhat_f"]) & (abs(est_alg1[, "rhat_f"] - r_f) < tol)
  
  orig_id_alg2 <- !is.na(est_alg2[, "rhat_e"]) & (abs(est_alg2[, "rhat_e"] - r_e) < tol)
  coll_id_alg2 <- !is.na(est_alg2[, "rhat_f"]) & (abs(est_alg2[, "rhat_f"] - r_f) < tol)
  
  data.frame(
    vol_model = vol_model,
    d         = d,
    eta       = eta,
    omega_g   = omega_g,
    alpha_g   = alpha_g,
    beta_g    = beta_g,
    c         = c_par,
    n         = n,
    r_e       = r_e,
    r_f       = r_f,
    alpha     = alpha,
    R         = R,
    tol       = tol,
    origination_id_rate_alg1 = mean(orig_id_alg1),
    collapse_id_rate_alg1    = mean(coll_id_alg1),
    origination_id_rate_alg2 = mean(orig_id_alg2),
    collapse_id_rate_alg2    = mean(coll_id_alg2),
    stringsAsFactors = FALSE
  )
}

# -----------------------------
# Run all scenarios
# -----------------------------
set.seed(123)

results_list <- list()
counter <- 1

for (j in seq_len(nrow(vol_grid))) {
  for (cc in c_grid) {
    for (i in seq_len(nrow(design_grid))) {
      results_list[[counter]] <- run_design_power(
        n         = n,
        r_e       = design_grid$r_e[i],
        r_f       = design_grid$r_f[i],
        c_par     = cc,
        alpha     = alpha,
        R         = R,
        vol_model = vol_grid$vol_model[j],
        d         = vol_grid$d[j],
        eta       = vol_grid$eta[j],
        omega_g   = vol_grid$omega_g[j],
        alpha_g   = vol_grid$alpha_g[j],
        beta_g    = vol_grid$beta_g[j],
        tol       = tol,
        min_gap   = min_gap
      )
      counter <- counter + 1
    }
  }
}

results_table <- do.call(rbind, results_list)

# Add readable volatility label
results_table$scenario_id <- NA_integer_
for (j in seq_len(nrow(vol_grid))) {
  pick <- results_table$vol_model == vol_grid$vol_model[j] &
    ((is.na(results_table$d) & is.na(vol_grid$d[j])) | (!is.na(results_table$d) & !is.na(vol_grid$d[j]) & results_table$d == vol_grid$d[j])) &
    ((is.na(results_table$eta) & is.na(vol_grid$eta[j])) | (!is.na(results_table$eta) & !is.na(vol_grid$eta[j]) & results_table$eta == vol_grid$eta[j])) &
    ((is.na(results_table$omega_g) & is.na(vol_grid$omega_g[j])) | (!is.na(results_table$omega_g) & !is.na(vol_grid$omega_g[j]) & results_table$omega_g == vol_grid$omega_g[j])) &
    ((is.na(results_table$alpha_g) & is.na(vol_grid$alpha_g[j])) | (!is.na(results_table$alpha_g) & !is.na(vol_grid$alpha_g[j]) & results_table$alpha_g == vol_grid$alpha_g[j])) &
    ((is.na(results_table$beta_g) & is.na(vol_grid$beta_g[j])) | (!is.na(results_table$beta_g) & !is.na(vol_grid$beta_g[j]) & results_table$beta_g == vol_grid$beta_g[j]))
  
  results_table$scenario_id[pick] <- vol_grid$scenario_id[j]
}

results_table <- merge(
  results_table,
  vol_grid[, c("scenario_id", "vol_label")],
  by = "scenario_id",
  all.x = TRUE,
  sort = FALSE
)

results_table <- results_table[
  order(results_table$scenario_id,
        results_table$c,
        results_table$r_e,
        results_table$r_f),
]

row.names(results_table) <- NULL

# Display table
results_table_print <- results_table[, c(
  "vol_label", "c", "n", "r_e", "r_f", "alpha", "R", "tol",
  "origination_id_rate_alg1", "collapse_id_rate_alg1",
  "origination_id_rate_alg2", "collapse_id_rate_alg2"
)]

num_cols <- sapply(results_table_print, is.numeric)
results_table_print[num_cols] <- lapply(
  results_table_print[num_cols],
  function(x) round(x, 4)
)

print(results_table_print)

# Optional: save full table
write.csv(
  results_table,
  file = "power_identification_rates_volatility_conditions_updated_alg2_dgp.csv",
  row.names = FALSE
)


###PWY approach comparison (BIAS MSE COMPARISON (WITHOUT NEW DGP))
# ------------------------------------------------------------
# Monte Carlo code:
# For each (r_e, r_f, alpha) design, report
#   - average estimated r_e and r_f
#   - bias and MSE of rhat_e and rhat_f
# under two algorithms:
#
# Algorithm 1:
#   origination cutoff = log(ns)/10
#   collapse cutoff    = log(ns)/2
#
# Algorithm 2:
#   origination cutoff = log(log(ns))/100
#   collapse cutoff    = log(log(ns))/100
#
# In both algorithms, rhat_f is searched only after rhat_e + 0.05.
# No plots; output is a table only.
# ------------------------------------------------------------

# -----------------------------
# Fixed parameters
# -----------------------------
n       <- 1000
c_par   <- 0.5
d_par   <- 1
eta_par <- 0.1
R       <- 1000   # number of Monte Carlo replications

# Design grid: edit as needed
design_grid <- expand.grid(
  r_e   = c(0.20, 0.30, 0.40, 0.50),
  r_f   = c(0.50, 0.65, 0.75, 0.85),
  alpha = c(0.30, 0.50, 0.70, 1.00),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# Keep only admissible designs with r_f > r_e + 0.05
design_grid <- subset(design_grid, r_f > r_e + 0.05)

# -----------------------------
# Data generating process
# -----------------------------
generate_series <- function(n, r_e, r_f, c, alpha, d, eta,
                            x0 = 5, sigma2_0 = 1) {
  tau_e <- floor(n * r_e)
  tau_f <- floor(n * r_f)
  stopifnot(tau_e >= 1, tau_f > tau_e, tau_f <= n)
  
  delta_n <- 1 + c / (n^alpha)
  phi_n   <- 1 - d / log(log(n))
  
  X <- numeric(n + 1)
  X[1] <- x0
  
  logsig2_prev <- log(sigma2_0)
  
  for (t in 1:n) {
    eta_t     <- rnorm(1, 0, eta)
    logsig2_t <- phi_n * logsig2_prev + eta_t
    sigma_t   <- sqrt(exp(logsig2_t))
    u_t       <- sigma_t * rnorm(1)
    
    a_t <- if (t < tau_e) {
      1
    } else if (t < tau_f) {
      delta_n
    } else {
      1
    }
    
    X[t + 1] <- a_t * X[t] + u_t
    logsig2_prev <- logsig2_t
  }
  
  X
}

# -----------------------------
# Recursive ADF path
# -----------------------------
adf_for_sgrid <- function(X_full, n, s_seq) {
  nX <- length(X_full) - 1
  stopifnot(nX >= max(floor(n * s_seq)))
  
  csX  <- cumsum(X_full)
  csX2 <- cumsum(X_full^2)
  xy   <- X_full[2:(n + 1)] * X_full[1:n]
  csXY <- cumsum(xy)
  
  tau_seq  <- pmax(2, floor(n * s_seq))
  DF_delta <- rep(NA_real_, length(tau_seq))
  
  for (k in seq_along(tau_seq)) {
    tau <- tau_seq[k]
    
    # Regression: X_t on (1, X_{t-1}), t = 1,...,tau
    S_reg   <- csX[tau]
    S_y     <- csX[tau + 1] - X_full[1]
    S_reg2  <- csX2[tau]
    S_y2    <- csX2[tau + 1] - X_full[1]^2
    S_cross <- csXY[tau]
    
    Xbar  <- S_y / tau
    SSR_t <- S_reg2 - 2 * Xbar * S_reg + tau * Xbar^2
    if (!is.finite(SSR_t) || SSR_t <= 0) next
    
    Sxtyt <- S_cross - Xbar * (S_reg + S_y) + tau * Xbar^2
    SSY_t <- S_y2    - 2 * Xbar * S_y       + tau * Xbar^2
    
    delta_hat  <- Sxtyt / SSR_t
    sigma2_hat <- (SSY_t - 2 * delta_hat * Sxtyt + (delta_hat^2) * SSR_t) / tau
    if (!is.finite(sigma2_hat) || sigma2_hat <= 0) next
    
    DF_delta[k] <- tau * (delta_hat - 1)
  }
  
  data.frame(s = tau_seq / n, DF_delta = DF_delta)
}

# -----------------------------
# Generic estimator
#   rhat_e = first s with DF_delta > cv_orig(s)
#   rhat_f = first s after rhat_e + min_gap with DF_delta < cv_coll(s)
# -----------------------------
estimate_breaks_generic <- function(X_full, n,
                                    s_seq = seq(0.10, 1.00, by = 0.001),
                                    min_gap = 0.05,
                                    cv_orig_fun,
                                    cv_coll_fun) {
  adf_df <- adf_for_sgrid(X_full, n, s_seq)
  
  cv_orig <- cv_orig_fun(n, adf_df$s)
  cv_coll <- cv_coll_fun(n, adf_df$s)
  
  idx_e  <- which(adf_df$DF_delta > cv_orig)[1]
  rhat_e <- if (!is.na(idx_e)) adf_df$s[idx_e] else NA_real_
  
  idx_f  <- NA_integer_
  rhat_f <- NA_real_
  
  if (!is.na(rhat_e)) {
    idx_start_f <- which(adf_df$s >= (rhat_e + min_gap))[1]
    if (!is.na(idx_start_f)) {
      idx_f_rel <- which(adf_df$DF_delta[idx_start_f:nrow(adf_df)] <
                           cv_coll[idx_start_f:nrow(adf_df)])[1]
      if (!is.na(idx_f_rel)) {
        idx_f  <- idx_start_f + idx_f_rel - 1
        rhat_f <- adf_df$s[idx_f]
      }
    }
  }
  
  c(rhat_e = rhat_e, rhat_f = rhat_f)
}

# -----------------------------
# Cutoff functions for the two algorithms
# -----------------------------
# Algorithm 1
cv1_orig_fun <- function(n, s) log(n * s) / 10
cv1_coll_fun <- function(n, s) log(n * s) / 2

# Algorithm 2
cv2_orig_fun <- function(n, s) log(log(n * s)) / 100
cv2_coll_fun <- function(n, s) log(log(n * s)) / 100

# -----------------------------
# One Monte Carlo design
# For each algorithm, report:
#   - average estimated r_e and r_f over non-NA estimates
#   - bias and MSE over non-NA estimates
# -----------------------------
run_design <- function(n, r_e, r_f, c_par, alpha, d_par, eta_par, R,
                       s_seq = seq(0.10, 1.00, by = 0.001),
                       min_gap = 0.05) {
  est_alg1 <- matrix(NA_real_, nrow = R, ncol = 2)
  est_alg2 <- matrix(NA_real_, nrow = R, ncol = 2)
  colnames(est_alg1) <- c("rhat_e", "rhat_f")
  colnames(est_alg2) <- c("rhat_e", "rhat_f")
  
  for (rep in 1:R) {
    X_full <- generate_series(
      n = n, r_e = r_e, r_f = r_f,
      c = c_par, alpha = alpha, d = d_par, eta = eta_par,
      x0 = 5
    )
    
    est_alg1[rep, ] <- estimate_breaks_generic(
      X_full      = X_full,
      n           = n,
      s_seq       = s_seq,
      min_gap     = min_gap,
      cv_orig_fun = cv1_orig_fun,
      cv_coll_fun = cv1_coll_fun
    )
    
    est_alg2[rep, ] <- estimate_breaks_generic(
      X_full      = X_full,
      n           = n,
      s_seq       = s_seq,
      min_gap     = min_gap,
      cv_orig_fun = cv2_orig_fun,
      cv_coll_fun = cv2_coll_fun
    )
  }
  
  # Errors: Algorithm 1
  err1_e <- est_alg1[, "rhat_e"] - r_e
  err1_f <- est_alg1[, "rhat_f"] - r_f
  
  # Errors: Algorithm 2
  err2_e <- est_alg2[, "rhat_e"] - r_e
  err2_f <- est_alg2[, "rhat_f"] - r_f
  
  out <- data.frame(
    n = n,
    r_e = r_e,
    r_f = r_f,
    alpha = alpha,
    R = R,
    
    # Algorithm 1: log(ns)/10 and log(ns)/2
    estimated_r_e_alg1 = mean(est_alg1[, "rhat_e"], na.rm = TRUE),
    estimated_r_f_alg1 = mean(est_alg1[, "rhat_f"], na.rm = TRUE),
    bias_rhat_e_alg1   = mean(err1_e, na.rm = TRUE),
    mse_rhat_e_alg1    = mean(err1_e^2, na.rm = TRUE),
    bias_rhat_f_alg1   = mean(err1_f, na.rm = TRUE),
    mse_rhat_f_alg1    = mean(err1_f^2, na.rm = TRUE),
    
    # Algorithm 2: loglog(ns)/100 for both
    estimated_r_e_alg2 = mean(est_alg2[, "rhat_e"], na.rm = TRUE),
    estimated_r_f_alg2 = mean(est_alg2[, "rhat_f"], na.rm = TRUE),
    bias_rhat_e_alg2   = mean(err2_e, na.rm = TRUE),
    mse_rhat_e_alg2    = mean(err2_e^2, na.rm = TRUE),
    bias_rhat_f_alg2   = mean(err2_f, na.rm = TRUE),
    mse_rhat_f_alg2    = mean(err2_f^2, na.rm = TRUE),
    
    stringsAsFactors = FALSE
  )
  
  # Replace NaN by NA when no estimate is ever found
  out[] <- lapply(out, function(x) {
    if (is.numeric(x)) x[is.nan(x)] <- NA_real_
    x
  })
  
  out
}

# -----------------------------
# Run all designs
# -----------------------------
set.seed(123)

results_list <- vector("list", nrow(design_grid))

for (i in seq_len(nrow(design_grid))) {
  results_list[[i]] <- run_design(
    n       = n,
    r_e     = design_grid$r_e[i],
    r_f     = design_grid$r_f[i],
    c_par   = c_par,
    alpha   = design_grid$alpha[i],
    d_par   = d_par,
    eta_par = eta_par,
    R       = R
  )
}

results_table <- do.call(rbind, results_list)

# Sort nicely
results_table <- results_table[order(results_table$alpha,
                                     results_table$r_e,
                                     results_table$r_f), ]

row.names(results_table) <- NULL

# Round for presentation
results_table_print <- results_table
num_cols <- sapply(results_table_print, is.numeric)
results_table_print[num_cols] <- lapply(
  results_table_print[num_cols],
  function(x) round(x, 4)
)

# Print table
print(results_table_print)




# -----------------------------
n       <- 1000
c_par   <- 0.5
d_par   <- 1
eta_par <- 0.1
R       <- 1000   # number of Monte Carlo replications

# Design grid: edit as needed
design_grid <- expand.grid(
  r_e   = c(0.20, 0.30, 0.40, 0.50),
  r_f   = c(0.50, 0.65, 0.75, 0.85),
  alpha = c(0.30, 0.50, 0.70, 1.00),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# Keep only admissible designs with r_f > r_e + 0.05
design_grid <- subset(design_grid, r_f > r_e + 0.05)

# -----------------------------
# Original data generating process
# -----------------------------
generate_series <- function(n, r_e, r_f, c, alpha, d, eta,
                            x0 = 5, sigma2_0 = 1) {
  tau_e <- floor(n * r_e)
  tau_f <- floor(n * r_f)
  stopifnot(tau_e >= 1, tau_f > tau_e, tau_f <= n)
  
  delta_n <- 1 + c / (n^alpha)
  phi_n   <- 1 - d / log(log(n))
  
  X <- numeric(n + 1)
  X[1] <- x0
  
  logsig2_prev <- log(sigma2_0)
  
  for (t in 1:n) {
    eta_t     <- rnorm(1, 0, eta)
    logsig2_t <- phi_n * logsig2_prev + eta_t
    sigma_t   <- sqrt(exp(logsig2_t))
    u_t       <- sigma_t * rnorm(1)
    
    a_t <- if (t < tau_e) {
      1
    } else if (t < tau_f) {
      delta_n
    } else {
      1
    }
    
    X[t + 1] <- a_t * X[t] + u_t
    logsig2_prev <- logsig2_t
  }
  
  X
}

# -----------------------------
# Modified DGP used only for Algorithm 2 collapse
# Up to tau_f same as original,
# then X_{tau_f+1} = X_{tau_e+1} + U, U ~ Unif(sqrt(n), 10*sqrt(n)),
# then continue as unit root with same SV volatility recursion.
# -----------------------------
generate_series_alg2_collapse <- function(n, r_e, r_f, c, alpha, d, eta,
                                          x0 = 5, sigma2_0 = 1) {
  tau_e <- floor(n * r_e)
  tau_f <- floor(n * r_f)
  stopifnot(tau_e >= 1, tau_f > tau_e, tau_f < n)
  
  delta_n <- 1 + c / (n^alpha)
  phi_n   <- 1 - d / log(log(n))
  
  X <- numeric(n + 1)
  X[1] <- x0
  
  logsig2_prev <- log(sigma2_0)
  
  # Simulate through t = tau_f
  for (t in 1:tau_f) {
    eta_t     <- rnorm(1, 0, eta)
    logsig2_t <- phi_n * logsig2_prev + eta_t
    sigma_t   <- sqrt(exp(logsig2_t))
    u_t       <- sigma_t * rnorm(1)
    
    a_t <- if (t < tau_e) 1 else delta_n
    X[t + 1] <- a_t * X[t] + u_t
    logsig2_prev <- logsig2_t
  }
  
  # Set X_{tau_f + 1} = X_{tau_e + 1} + U
  jump_u <- runif(1, min = sqrt(n), max = 10 * sqrt(n))
  X[tau_f + 2] <- X[tau_e + 1] + jump_u
  
  # Continue as unit root from time tau_f + 2 onward
  if (tau_f + 2 <= n) {
    for (t in (tau_f + 2):n) {
      eta_t     <- rnorm(1, 0, eta)
      logsig2_t <- phi_n * logsig2_prev + eta_t
      sigma_t   <- sqrt(exp(logsig2_t))
      u_t       <- sigma_t * rnorm(1)
      
      X[t + 1] <- X[t] + u_t
      logsig2_prev <- logsig2_t
    }
  }
  
  X
}

# -----------------------------
# Recursive ADF path
# -----------------------------
adf_for_sgrid <- function(X_full, n, s_seq) {
  nX <- length(X_full) - 1
  stopifnot(nX >= max(floor(n * s_seq)))
  
  csX  <- cumsum(X_full)
  csX2 <- cumsum(X_full^2)
  xy   <- X_full[2:(n + 1)] * X_full[1:n]
  csXY <- cumsum(xy)
  
  tau_seq  <- pmax(2, floor(n * s_seq))
  DF_delta <- rep(NA_real_, length(tau_seq))
  
  for (k in seq_along(tau_seq)) {
    tau <- tau_seq[k]
    
    # Regression: X_t on (1, X_{t-1}), t = 1,...,tau
    S_reg   <- csX[tau]
    S_y     <- csX[tau + 1] - X_full[1]
    S_reg2  <- csX2[tau]
    S_y2    <- csX2[tau + 1] - X_full[1]^2
    S_cross <- csXY[tau]
    
    Xbar  <- S_y / tau
    SSR_t <- S_reg2 - 2 * Xbar * S_reg + tau * Xbar^2
    if (!is.finite(SSR_t) || SSR_t <= 0) next
    
    Sxtyt <- S_cross - Xbar * (S_reg + S_y) + tau * Xbar^2
    SSY_t <- S_y2    - 2 * Xbar * S_y       + tau * Xbar^2
    
    delta_hat  <- Sxtyt / SSR_t
    sigma2_hat <- (SSY_t - 2 * delta_hat * Sxtyt + (delta_hat^2) * SSR_t) / tau
    if (!is.finite(sigma2_hat) || sigma2_hat <= 0) next
    
    DF_delta[k] <- tau * (delta_hat - 1)
  }
  
  data.frame(s = tau_seq / n, DF_delta = DF_delta)
}

# -----------------------------
# Generic estimator
#   rhat_e = first s with DF_delta > cv_orig(s)
#   rhat_f = first s after rhat_e + min_gap with DF_delta < cv_coll(s)
# -----------------------------
estimate_breaks_generic <- function(X_full, n,
                                    s_seq = seq(0.10, 1.00, by = 0.001),
                                    min_gap = 0.05,
                                    cv_orig_fun,
                                    cv_coll_fun) {
  adf_df <- adf_for_sgrid(X_full, n, s_seq)
  
  cv_orig <- cv_orig_fun(n, adf_df$s)
  cv_coll <- cv_coll_fun(n, adf_df$s)
  
  idx_e  <- which(adf_df$DF_delta > cv_orig)[1]
  rhat_e <- if (!is.na(idx_e)) adf_df$s[idx_e] else NA_real_
  
  idx_f  <- NA_integer_
  rhat_f <- NA_real_
  
  if (!is.na(rhat_e)) {
    idx_start_f <- which(adf_df$s >= (rhat_e + min_gap))[1]
    if (!is.na(idx_start_f)) {
      idx_f_rel <- which(adf_df$DF_delta[idx_start_f:nrow(adf_df)] <
                           cv_coll[idx_start_f:nrow(adf_df)])[1]
      if (!is.na(idx_f_rel)) {
        idx_f  <- idx_start_f + idx_f_rel - 1
        rhat_f <- adf_df$s[idx_f]
      }
    }
  }
  
  c(rhat_e = rhat_e, rhat_f = rhat_f)
}

# -----------------------------
# Cutoff functions for the two algorithms
# -----------------------------
# Algorithm 1
cv1_orig_fun <- function(n, s) log(n * s) / 10
cv1_coll_fun <- function(n, s) log(n * s) / 2

# Algorithm 2
cv2_orig_fun <- function(n, s) log(log(n * s)) / 100
cv2_coll_fun <- function(n, s) log(log(n * s)) / 100

# -----------------------------
# One Monte Carlo design
# For Algorithm 1: everything from original DGP
# For Algorithm 2:
#   - rhat_e comes from the original DGP
#   - rhat_f comes from the modified collapse DGP,
#     starting search after the original-path rhat_e + min_gap
# Hence only estimated_r_f_alg2, bias_rhat_f_alg2, mse_rhat_f_alg2 change.
# -----------------------------
run_design <- function(n, r_e, r_f, c_par, alpha, d_par, eta_par, R,
                       s_seq = seq(0.10, 1.00, by = 0.001),
                       min_gap = 0.05) {
  est_alg1 <- matrix(NA_real_, nrow = R, ncol = 2)
  est_alg2 <- matrix(NA_real_, nrow = R, ncol = 2)
  colnames(est_alg1) <- c("rhat_e", "rhat_f")
  colnames(est_alg2) <- c("rhat_e", "rhat_f")
  
  for (rep in 1:R) {
    # Original path: used for Algorithm 1 and Algorithm 2 origination
    X_full_orig <- generate_series(
      n = n, r_e = r_e, r_f = r_f,
      c = c_par, alpha = alpha, d = d_par, eta = eta_par,
      x0 = 5
    )
    
    # Algorithm 1 unchanged
    est_alg1[rep, ] <- estimate_breaks_generic(
      X_full      = X_full_orig,
      n           = n,
      s_seq       = s_seq,
      min_gap     = min_gap,
      cv_orig_fun = cv1_orig_fun,
      cv_coll_fun = cv1_coll_fun
    )
    
    # Algorithm 2 origination stays on the original path
    est_alg2_orig <- estimate_breaks_generic(
      X_full      = X_full_orig,
      n           = n,
      s_seq       = s_seq,
      min_gap     = min_gap,
      cv_orig_fun = cv2_orig_fun,
      cv_coll_fun = cv2_coll_fun
    )
    
    # Modified path: used only for Algorithm 2 collapse
    X_full_mod <- generate_series_alg2_collapse(
      n = n, r_e = r_e, r_f = r_f,
      c = c_par, alpha = alpha, d = d_par, eta = eta_par,
      x0 = 5
    )
    
    adf_df_mod <- adf_for_sgrid(X_full_mod, n, s_seq)
    cv_coll2   <- cv2_coll_fun(n, adf_df_mod$s)
    
    rhat_e_alg2 <- est_alg2_orig["rhat_e"]
    rhat_f_alg2 <- NA_real_
    
    if (!is.na(rhat_e_alg2)) {
      idx_start_f <- which(adf_df_mod$s >= (rhat_e_alg2 + min_gap))[1]
      
      if (!is.na(idx_start_f)) {
        idx_f_rel <- which(adf_df_mod$DF_delta[idx_start_f:nrow(adf_df_mod)] <
                             cv_coll2[idx_start_f:nrow(adf_df_mod)])[1]
        
        if (!is.na(idx_f_rel)) {
          idx_f_alg2 <- idx_start_f + idx_f_rel - 1
          rhat_f_alg2 <- adf_df_mod$s[idx_f_alg2]
        }
      }
    }
    
    est_alg2[rep, "rhat_e"] <- rhat_e_alg2
    est_alg2[rep, "rhat_f"] <- rhat_f_alg2
  }
  
  # Errors: Algorithm 1
  err1_e <- est_alg1[, "rhat_e"] - r_e
  err1_f <- est_alg1[, "rhat_f"] - r_f
  
  # Errors: Algorithm 2
  # rhat_e uses original DGP, rhat_f uses modified collapse DGP
  err2_e <- est_alg2[, "rhat_e"] - r_e
  err2_f <- est_alg2[, "rhat_f"] - r_f
  
  out <- data.frame(
    n = n,
    r_e = r_e,
    r_f = r_f,
    alpha = alpha,
    R = R,
    
    # Algorithm 1: unchanged
    estimated_r_e_alg1 = mean(est_alg1[, "rhat_e"], na.rm = TRUE),
    estimated_r_f_alg1 = mean(est_alg1[, "rhat_f"], na.rm = TRUE),
    bias_rhat_e_alg1   = mean(err1_e, na.rm = TRUE),
    mse_rhat_e_alg1    = mean(err1_e^2, na.rm = TRUE),
    bias_rhat_f_alg1   = mean(err1_f, na.rm = TRUE),
    mse_rhat_f_alg1    = mean(err1_f^2, na.rm = TRUE),
    
    # Algorithm 2:
    # r_e components unchanged; r_f components updated
    estimated_r_e_alg2 = mean(est_alg2[, "rhat_e"], na.rm = TRUE),
    estimated_r_f_alg2 = mean(est_alg2[, "rhat_f"], na.rm = TRUE),
    bias_rhat_e_alg2   = mean(err2_e, na.rm = TRUE),
    mse_rhat_e_alg2    = mean(err2_e^2, na.rm = TRUE),
    bias_rhat_f_alg2   = mean(err2_f, na.rm = TRUE),
    mse_rhat_f_alg2    = mean(err2_f^2, na.rm = TRUE),
    
    stringsAsFactors = FALSE
  )
  
  # Replace NaN by NA when no estimate is ever found
  out[] <- lapply(out, function(x) {
    if (is.numeric(x)) x[is.nan(x)] <- NA_real_
    x
  })
  
  out
}

# -----------------------------
# Run all designs
# -----------------------------
set.seed(123)

results_list <- vector("list", nrow(design_grid))

for (i in seq_len(nrow(design_grid))) {
  results_list[[i]] <- run_design(
    n       = n,
    r_e     = design_grid$r_e[i],
    r_f     = design_grid$r_f[i],
    c_par   = c_par,
    alpha   = design_grid$alpha[i],
    d_par   = d_par,
    eta_par = eta_par,
    R       = R
  )
}

results_table <- do.call(rbind, results_list)

# Sort nicely
results_table <- results_table[order(results_table$alpha,
                                     results_table$r_e,
                                     results_table$r_f), ]

row.names(results_table) <- NULL

# Round for presentation
results_table_print <- results_table
num_cols <- sapply(results_table_print, is.numeric)
results_table_print[num_cols] <- lapply(
  results_table_print[num_cols],
  function(x) round(x, 4)
)

# Print table
print(results_table_print)

# Optional: save to CSV
write.csv(results_table_print,
          file = "bias_mse_estimated_re_rf_two_algorithms_DGP.csv",
          row.names = FALSE)