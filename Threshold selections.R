################################################################################
# Threshold selections.R
#
# Replication code for Figure 12 of the main paper and Figure A.3 / Table A.1
# of the Online Appendix for:
#
# "Is There an AI Bubble? Robust Date-Stamping for Periods of Exuberance"
#
# Figure 12:
#   One simulated recursive SV-ADF path with:
#     - true bubble origination fraction r_e,
#     - true bubble collapse fraction r_f,
#     - estimated origination fraction \hat r_e,
#     - estimated collapse fraction \hat r_f.
#
# Online Appendix Figure A.3:
#   Simulated coefficient-based critical values under H0 and H1.
#
# Online Appendix Table A.1:
#   Simulated coefficient-based critical values for n = 500, 550, ..., 1000.
#
# Notes:
# - Package installation is intentionally omitted; see README.md.
# - This script uses only base R.
# - Set SAVE_FIGURES <- TRUE to save PDF/PNG figures.
# - Set RECOMPUTE_CRITICAL_VALUES <- TRUE to rerun the Monte Carlo simulation
#   for Table A.1 and Figure A.3. Otherwise, the paper values are used directly.
################################################################################

# ------------------------------------------------------------------------------
# Global settings
# ------------------------------------------------------------------------------

SAVE_FIGURES <- TRUE
SAVE_TABLES <- TRUE
SAVE_PNG_PREVIEW <- TRUE

RECOMPUTE_CRITICAL_VALUES <- FALSE

OUTPUT_DIR <- "figures"
TABLE_DIR <- "tables"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

if (!dir.exists(TABLE_DIR)) {
  dir.create(TABLE_DIR, recursive = TRUE)
}

FIGURE12_PDF <- file.path(OUTPUT_DIR, "sv_adf_plot_simulated_one.pdf")
FIGURE12_PNG <- file.path(OUTPUT_DIR, "sv_adf_plot_simulated_one.png")

APPENDIX_FIGURE_A3_PDF <- file.path(OUTPUT_DIR, "critical_values_plot.pdf")
APPENDIX_FIGURE_A3_PNG <- file.path(OUTPUT_DIR, "critical_values_plot.png")

APPENDIX_TABLE_A1_CSV <- file.path(TABLE_DIR, "tableA1_critical_values.csv")

################################################################################
# Helper functions
################################################################################

# Axis-label helper.
# Large values are displayed in scientific notation, e.g. 10^5 or 2 * 10^5.
sci_tick_labels <- function(vals) {
  out <- character(length(vals))
  
  for (i in seq_along(vals)) {
    v <- vals[i]
    
    if (!is.finite(v)) {
      out[i] <- "''"
      next
    }
    
    if (abs(v) < .Machine$double.eps) {
      out[i] <- "0"
      next
    }
    
    abs_v <- abs(v)
    exponent <- floor(log10(abs_v))
    
    if (exponent >= 4 || exponent <= -3) {
      multiplier <- round(v / (10^exponent), 2)
      
      if (abs(multiplier - round(multiplier)) < 1e-8) {
        multiplier <- round(multiplier)
      }
      
      if (abs(multiplier - 1) < 1e-8) {
        out[i] <- paste0("10^", exponent)
      } else if (abs(multiplier + 1) < 1e-8) {
        out[i] <- paste0("-10^", exponent)
      } else {
        out[i] <- paste0(multiplier, "%*%10^", exponent)
      }
    } else {
      out[i] <- format(round(v, 2), trim = TRUE, scientific = FALSE)
    }
  }
  
  parse(text = out)
}

# Generate X_0, ..., X_n under stochastic volatility with one bubble episode.
#
# Model:
#   X_t = a_t X_{t-1} + u_t,
#   u_t = sigma_t epsilon_t,
#   log(sigma_t^2) = phi_n log(sigma_{t-1}^2) + eta_t,
#
# where:
#   a_t = 1 before the bubble,
#   a_t = delta_n during the bubble,
#   a_t = 1 after collapse.
generate_series_bubble_sv <- function(
    n,
    r_e,
    r_f,
    c,
    alpha,
    d,
    eta,
    x0 = 5,
    sigma2_0 = 1
) {
  tau_e <- floor(n * r_e)
  tau_f <- floor(n * r_f)
  
  stopifnot(tau_e >= 1, tau_f > tau_e, tau_f <= n)
  
  delta_n <- 1 + c / (n^alpha)
  phi_n <- 1 - d / log(log(n))
  
  X <- numeric(n + 1)
  X[1] <- x0
  
  log_sigma2_prev <- log(sigma2_0)
  
  for (t in seq_len(n)) {
    eta_t <- rnorm(1, mean = 0, sd = eta)
    log_sigma2_t <- phi_n * log_sigma2_prev + eta_t
    sigma_t <- sqrt(exp(log_sigma2_t))
    u_t <- sigma_t * rnorm(1)
    
    a_t <- if (t < tau_e) {
      1
    } else if (t < tau_f) {
      delta_n
    } else {
      1
    }
    
    X[t + 1] <- a_t * X[t] + u_t
    log_sigma2_prev <- log_sigma2_t
  }
  
  X
}

# Generate X_0, ..., X_n under the unit-root null with stochastic volatility.
#
# Model:
#   X_t = X_{t-1} + u_t.
generate_series_h0_sv <- function(
    n,
    d,
    eta,
    x0 = 100,
    sigma2_0 = 1
) {
  phi_n <- 1 - d / log(log(n))
  
  X <- numeric(n + 1)
  X[1] <- x0
  
  log_sigma2_prev <- log(sigma2_0)
  
  for (t in seq_len(n)) {
    eta_t <- rnorm(1, mean = 0, sd = eta)
    log_sigma2_t <- phi_n * log_sigma2_prev + eta_t
    sigma_t <- sqrt(exp(log_sigma2_t))
    u_t <- sigma_t * rnorm(1)
    
    X[t + 1] <- X[t] + u_t
    log_sigma2_prev <- log_sigma2_t
  }
  
  X
}

# Compute recursive coefficient-based and t-type ADF statistics over an s-grid.
#
# The coefficient-based statistic is:
#   DF_delta(s) = tau * (delta_hat_tau - 1),
# where tau = floor(ns).
adf_for_sgrid <- function(X_full, n, s_seq) {
  max_tau <- max(floor(n * s_seq))
  nX <- length(X_full) - 1
  
  stopifnot(nX >= max_tau)
  
  csX <- cumsum(X_full)
  csX2 <- cumsum(X_full^2)
  xy <- X_full[2:(n + 1)] * X_full[1:n]
  csXY <- cumsum(xy)
  
  tau_seq <- pmax(2, floor(n * s_seq))
  DF_delta <- rep(NA_real_, length(tau_seq))
  DF_t <- rep(NA_real_, length(tau_seq))
  
  for (k in seq_along(tau_seq)) {
    tau <- tau_seq[k]
    
    # Regression: X_t on a constant and X_{t-1}, for t = 1, ..., tau.
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
    DF_t[k] <- (delta_hat - 1) * sqrt(SSR_t) / sqrt(sigma2_hat)
  }
  
  data.frame(
    s = tau_seq / n,
    DF_delta = DF_delta,
    DF_t = DF_t
  )
}

# Full-sample coefficient-based DF statistic without intercept:
#   DF_delta = n * (delta_hat - 1).
df_delta_fullsample <- function(X_full) {
  n <- length(X_full) - 1
  
  Xlag <- X_full[1:n]
  Y <- X_full[2:(n + 1)]
  
  denom <- sum(Xlag^2)
  
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }
  
  delta_hat <- sum(Y * Xlag) / denom
  
  n * (delta_hat - 1)
}

# Coefficient-based DF statistic at a fixed recursive endpoint tau.
adf_delta_at_tau <- function(X_full, tau) {
  Xlag <- X_full[1:tau]
  Y <- X_full[2:(tau + 1)]
  
  Ybar <- mean(Y)
  SSR <- sum((Xlag - Ybar)^2)
  
  if (!is.finite(SSR) || SSR <= 0) {
    return(NA_real_)
  }
  
  Sxy <- sum((Xlag - Ybar) * (Y - Ybar))
  delta_hat <- Sxy / SSR
  
  tau * (delta_hat - 1)
}

################################################################################
# Figure 12
# One simulated recursive SV-ADF path
################################################################################

make_simulated_svadf_path_plot <- function(
    n = 1000,
    r_e = 0.3,
    r_f = 0.6,
    c_par = 0.5,
    alpha = 0.5,
    d_par = 1,
    eta_par = 0.1,
    seed = 123,
    save_figure = SAVE_FIGURES
) {
  set.seed(seed)
  
  X_full <- generate_series_bubble_sv(
    n = n,
    r_e = r_e,
    r_f = r_f,
    c = c_par,
    alpha = alpha,
    d = d_par,
    eta = eta_par,
    x0 = 10,
    sigma2_0 = 1
  )
  
  s_seq <- seq(0.1, 1.00, by = 0.001)
  adf_df <- adf_for_sgrid(X_full, n, s_seq)
  
  # Origination and collapse threshold curves.
  cv_orig <- log(n * adf_df$s) / 10
  cv_coll <- log(n * adf_df$s) / 2
  
  # Estimated origination:
  # first s where DF_delta exceeds log(ns)/10.
  idx_e <- which(adf_df$DF_delta > cv_orig)[1]
  rhat_e <- if (!is.na(idx_e)) adf_df$s[idx_e] else NA_real_
  
  # Estimated collapse:
  # first s after rhat_e + 0.05 where DF_delta falls below log(ns)/2.
  idx_f <- NA_integer_
  rhat_f <- NA_real_
  
  if (!is.na(rhat_e)) {
    idx_start_f <- which(adf_df$s >= rhat_e + 0.05)[1]
    
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
  
  cat(sprintf("Estimated r_e (rhat_e) = %.3f\n", rhat_e))
  cat(sprintf("Estimated r_f (rhat_f) = %.3f\n", rhat_f))
  
  # Prepare |X_t| for the right y-axis overlay.
  t_idx <- 0:n
  s_x <- t_idx / n
  X_vals <- abs(X_full)
  
  keep <- s_x >= min(adf_df$s) & s_x <= max(adf_df$s)
  s_x <- s_x[keep]
  X_vals <- X_vals[keep]
  
  x_ylim <- range(X_vals, finite = TRUE)
  
  left_ticks <- pretty(range(c(adf_df$DF_delta, cv_orig, cv_coll), finite = TRUE))
  right_ticks <- pretty(x_ylim)
  
  draw_plot <- function() {
    par(mfrow = c(1, 1), mar = c(5.5, 5, 1.5, 5) + 0.1)
    
    plot(
      adf_df$s,
      adf_df$DF_delta,
      type = "l",
      col = "steelblue",
      lwd = 2,
      xlab = "Proportion of samples for recursive SV-ADF test",
      ylab = "SV-ADF Statistic",
      main = "",
      xaxt = "s",
      yaxt = "n"
    )
    
    axis(1)
    axis(2, at = left_ticks, labels = sci_tick_labels(left_ticks), las = 1)
    
    lines(adf_df$s, cv_orig, col = "orange", lwd = 2, lty = 1)
    lines(adf_df$s, cv_coll, col = "purple", lwd = 2, lty = 1)
    
    abline(v = r_e, col = "red", lty = 1, lwd = 2)
    abline(v = r_f, col = "brown", lty = 1, lwd = 2)
    
    if (!is.na(rhat_e)) {
      abline(v = rhat_e, col = "blue", lty = 1, lwd = 2)
    }
    
    if (!is.na(rhat_f)) {
      abline(v = rhat_f, col = "darkgreen", lty = 1, lwd = 2)
    }
    
    par(new = TRUE)
    
    plot(
      s_x,
      X_vals,
      type = "l",
      col = rgb(0.3, 0.3, 0.3, 0.7),
      lwd = 1.5,
      axes = FALSE,
      xlab = "",
      ylab = "",
      xlim = range(adf_df$s),
      ylim = x_ylim
    )
    
    axis(4, at = right_ticks, labels = sci_tick_labels(right_ticks), las = 1)
    mtext("Price process Time Series", side = 4, line = 3)
    
    legend(
      "topleft",
      legend = c(
        expression("True " * r[e]),
        expression("True " * r[f]),
        expression(hat(r)[e]),
        expression(hat(r)[f])
      ),
      col = c("red", "brown", "blue", "darkgreen"),
      lty = c(1, 1, 1, 1),
      lwd = c(2, 2, 2, 2),
      bty = "n"
    )
  }
  
  draw_plot()
  
  if (isTRUE(save_figure)) {
    grDevices::pdf(FIGURE12_PDF, width = 10, height = 6)
    draw_plot()
    grDevices::dev.off()
    
    if (isTRUE(SAVE_PNG_PREVIEW)) {
      grDevices::png(FIGURE12_PNG, width = 10, height = 6, units = "in", res = 300)
      draw_plot()
      grDevices::dev.off()
    }
  }
  
  invisible(list(
    X_full = X_full,
    adf_df = adf_df,
    rhat_e = rhat_e,
    rhat_f = rhat_f
  ))
}

################################################################################
# Online Appendix Table A.1 and Figure A.3
# Critical values under H0 and H1
################################################################################

# Simulate the 90% upper critical value under H0 for one sample size n.
critical_value_h0_upper <- function(
    n,
    B,
    d,
    eta,
    x0 = 1,
    sigma2_0 = 1
) {
  vals <- numeric(B)
  
  for (b in seq_len(B)) {
    X_full <- generate_series_h0_sv(
      n = n,
      d = d,
      eta = eta,
      x0 = x0,
      sigma2_0 = sigma2_0
    )
    
    vals[b] <- df_delta_fullsample(X_full)
  }
  
  stats::quantile(
    vals,
    probs = 0.90,
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )
}

# Simulate the 10% lower critical value under H1 for one sample size n.
critical_value_h1_lower <- function(
    n,
    B,
    r_e,
    r_f,
    c,
    alpha,
    d,
    eta,
    x0 = 1,
    sigma2_0 = 1
) {
  tau_after_rf <- floor(n * r_f)
  
  if (tau_after_rf > n) {
    stop("tau_after_rf exceeds n")
  }
  
  vals <- numeric(B)
  
  for (b in seq_len(B)) {
    X_full <- generate_series_bubble_sv(
      n = n,
      r_e = r_e,
      r_f = r_f,
      c = c,
      alpha = alpha,
      d = d,
      eta = eta,
      x0 = x0,
      sigma2_0 = sigma2_0
    )
    
    vals[b] <- adf_delta_at_tau(X_full, tau = tau_after_rf)
  }
  
  stats::quantile(
    vals,
    probs = 0.10,
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )
}

# Recompute critical values by Monte Carlo.
# This can be time-consuming, so the script also allows use of the paper values.
simulate_critical_value_table <- function(
    B = 1000,
    n_grid = seq(500, 1000, by = 50),
    d_par = 0.01,
    eta_par = 0.5,
    r_e = 0.3,
    r_f = 0.4,
    c_par = 1,
    alpha = 0.5,
    seed = 123
) {
  set.seed(seed)
  
  crit_90 <- numeric(length(n_grid))
  crit_10 <- numeric(length(n_grid))
  
  for (i in seq_along(n_grid)) {
    n <- n_grid[i]
    
    crit_90[i] <- critical_value_h0_upper(
      n = n,
      B = B,
      d = d_par,
      eta = eta_par,
      x0 = 1,
      sigma2_0 = 1
    )
    
    crit_10[i] <- critical_value_h1_lower(
      n = n,
      B = B,
      r_e = r_e,
      r_f = r_f,
      c = c_par,
      alpha = alpha,
      d = d_par,
      eta = eta_par,
      x0 = 1,
      sigma2_0 = 1
    )
    
    cat(sprintf(
      "Done n = %d, cv_90_H0 = %.4f, cv_10_H1 = %.4f\n",
      n,
      crit_90[i],
      crit_10[i]
    ))
  }
  
  data.frame(
    n = n_grid,
    cv_90_upper_H0 = crit_90,
    cv_10_lower_H1 = crit_10
  )
}

# Paper values used in Online Appendix Table A.1.
get_paper_critical_value_table <- function() {
  data.frame(
    n = c(500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000),
    cv_90_upper_H0 = c(
      0.6463, 0.5563, 0.6304, 0.7157, 0.7097, 0.7037,
      0.7823, 0.7888, 0.8174, 0.8374, 0.8403
    ),
    cv_10_lower_H1 = c(
      -0.2566, -0.9916, 1.2746, 1.9418, 1.0183, 2.6310,
      3.2163, 3.8318, 4.6070, 3.2394, 4.1479
    )
  )
}

make_critical_value_plot <- function(
    cv_table,
    save_figure = SAVE_FIGURES
) {
  y_limits <- range(
    c(cv_table$cv_90_upper_H0, cv_table$cv_10_lower_H1),
    finite = TRUE
  )
  
  draw_plot <- function() {
    plot(
      cv_table$n,
      cv_table$cv_90_upper_H0,
      type = "b",
      pch = 19,
      lwd = 2,
      col = "firebrick",
      ylim = y_limits,
      xlab = "Sample Size",
      ylab = "Critical value",
      main = expression(paste("Coefficient-based critical values"))
    )
    
    lines(
      cv_table$n,
      cv_table$cv_10_lower_H1,
      type = "b",
      pch = 17,
      lwd = 2,
      col = "steelblue"
    )
    
    legend(
      "topleft",
      legend = c(
        expression(cv[0.10,H[0]]^delta),
        expression(cv[0.10,H[1]]^delta)
      ),
      col = c("firebrick", "steelblue"),
      lty = 1,
      pch = c(19, 17),
      lwd = 2,
      bty = "n"
    )
  }
  
  draw_plot()
  
  if (isTRUE(save_figure)) {
    grDevices::pdf(APPENDIX_FIGURE_A3_PDF, width = 8, height = 5)
    draw_plot()
    grDevices::dev.off()
    
    if (isTRUE(SAVE_PNG_PREVIEW)) {
      grDevices::png(APPENDIX_FIGURE_A3_PNG, width = 8, height = 5, units = "in", res = 300)
      draw_plot()
      grDevices::dev.off()
    }
  }
  
  invisible(cv_table)
}

################################################################################
# Run all outputs
################################################################################

# Figure 12 of the main paper.
figure12_output <- make_simulated_svadf_path_plot(
  n = 1000,
  r_e = 0.3,
  r_f = 0.6,
  c_par = 0.5,
  alpha = 0.5,
  d_par = 1,
  eta_par = 0.1,
  seed = 123,
  save_figure = SAVE_FIGURES
)

# Table A.1 and Figure A.3 of the Online Appendix.
critical_value_table <- if (isTRUE(RECOMPUTE_CRITICAL_VALUES)) {
  simulate_critical_value_table(
    B = 1000,
    n_grid = seq(500, 1000, by = 50),
    d_par = 0.01,
    eta_par = 0.5,
    r_e = 0.3,
    r_f = 0.4,
    c_par = 1,
    alpha = 0.5,
    seed = 123
  )
} else {
  get_paper_critical_value_table()
}

print(critical_value_table)

if (isTRUE(SAVE_TABLES)) {
  write.csv(
    critical_value_table,
    file = APPENDIX_TABLE_A1_CSV,
    row.names = FALSE
  )
}

make_critical_value_plot(
  cv_table = critical_value_table,
  save_figure = SAVE_FIGURES
)

################################################################################
# End of script
################################################################################
