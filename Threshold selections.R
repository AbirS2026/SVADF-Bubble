##Simulations
# One path shown as an example 
n       <- 1000
r_e     <- 0.3
r_f     <- 0.6
c_par   <- 0.5
alpha   <- 0.5
d_par   <- 1
eta_par <- 0.1

# Generate X_0,...,X_n (length n+1) with bubble origination at tau_e
# and collapse at tau_f:
#   t < tau_e            : a_t = 1
#   tau_e <= t < tau_f   : a_t = delta_n
#   t >= tau_f           : a_t = 1
generate_series <- function(n, r_e, r_f, c, alpha, d, eta,
                            x0 = 5, sigma2_0 = 1) {
  tau_e   <- floor(n * r_e)
  tau_f   <- floor(n * r_f)
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

# Compute ADF paths for a given s-grid (tau = floor(n s))
adf_for_sgrid <- function(X_full, n, s_seq) {
  nX <- length(X_full) - 1
  stopifnot(nX >= max(floor(n * s_seq)))
  
  csX  <- cumsum(X_full)
  csX2 <- cumsum(X_full^2)
  xy   <- X_full[2:(n + 1)] * X_full[1:n]
  csXY <- cumsum(xy)
  
  tau_seq  <- pmax(2, floor(n * s_seq))
  DF_delta <- DF_t <- rep(NA_real_, length(tau_seq))
  
  for (k in seq_along(tau_seq)) {
    tau <- tau_seq[k]
    
    # Regression: X_t on (1, X_{t-1}), t = 1,...,tau
    S_reg   <- csX[tau]                       # sum X_0,...,X_{tau-1}
    S_y     <- csX[tau + 1] - X_full[1]       # sum X_1,...,X_tau
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
    DF_t[k]     <- (delta_hat - 1) * sqrt(SSR_t) / sqrt(sigma2_hat)
  }
  
  data.frame(s = tau_seq / n, DF_delta = DF_delta, DF_t = DF_t)
}

# Axis label helper: show large values as 10^5, 2 x 10^5, etc.
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
    
    av <- abs(v)
    e  <- floor(log10(av))
    
    if (e >= 4 || e <= -3) {
      m <- round(v / (10^e), 2)
      
      if (abs(m - round(m)) < 1e-8) {
        m <- round(m)
      }
      
      if (abs(m - 1) < 1e-8) {
        out[i] <- paste0("10^", e)
      } else if (abs(m + 1) < 1e-8) {
        out[i] <- paste0("-10^", e)
      } else {
        out[i] <- paste0(m, "%*%10^", e)
      }
    } else {
      out[i] <- format(round(v, 2), trim = TRUE, scientific = FALSE)
    }
  }
  
  parse(text = out)
}

# -----------------------------
# Run one replicate
# -----------------------------
set.seed(123)

X_full <- generate_series(n, r_e, r_f, c_par, alpha, d_par, eta_par, x0 = 10)
s_seq  <- seq(0.1, 1.00, by = 0.001)
adf_df <- adf_for_sgrid(X_full, n, s_seq)

# Critical value curves
cv_orig <- log(n * adf_df$s) / 10
cv_coll <- log(n * adf_df$s) / 2

# Estimated origination time:
# first s where DF_delta exceeds log(ns)/10
idx_e  <- which(adf_df$DF_delta > cv_orig)[1]
rhat_e <- if (!is.na(idx_e)) adf_df$s[idx_e] else NA_real_

# Estimated collapse time:
# first s after rhat_e + 0.05 where DF_delta falls below log(ns)/2
idx_f  <- NA_integer_
rhat_f <- NA_real_

if (!is.na(rhat_e)) {
  idx_start_f <- which(adf_df$s >= (rhat_e + 0.05))[1]
  if (!is.na(idx_start_f)) {
    idx_f_rel <- which(adf_df$DF_delta[idx_start_f:nrow(adf_df)] <
                         cv_coll[idx_start_f:nrow(adf_df)])[1]
    if (!is.na(idx_f_rel)) {
      idx_f  <- idx_start_f + idx_f_rel - 1
      rhat_f <- adf_df$s[idx_f]
    }
  }
}

cat(sprintf("estimated r_e (rhat_e) = %.3f\n", rhat_e))
cat(sprintf("estimated r_f (rhat_f) = %.3f\n", rhat_f))

# Prepare |X_t| series for right y-axis overlay
t_idx  <- 0:n
s_x    <- t_idx / n
X_vals <- abs(X_full)

keep   <- s_x >= min(adf_df$s) & s_x <= max(adf_df$s)
s_x    <- s_x[keep]
X_vals <- X_vals[keep]
x_ylim <- range(X_vals, finite = TRUE)

# Custom axis ticks/labels
left_ticks  <- pretty(range(c(adf_df$DF_delta, cv_orig, cv_coll), finite = TRUE))
right_ticks <- pretty(x_ylim)


pdf("sv_adf_plot_simulated_one.pdf", width = 10, height = 6)
par(mfrow = c(1, 1), mar = c(5.5, 5, 1.5, 5) + 0.1)

plot(adf_df$s, adf_df$DF_delta, type = "l", col = "steelblue", lwd = 2,
     xlab = "Proportion of samples for recursive SV-ADF test",
     ylab = "SV-ADF Statistic",
     main = "",
     xaxt = "s", yaxt = "n")

axis(1)
axis(2, at = left_ticks, labels = sci_tick_labels(left_ticks), las = 1)

lines(adf_df$s, cv_orig, col = "orange", lwd = 2, lty = 1)
lines(adf_df$s, cv_coll, col = "purple", lwd = 2, lty = 1)

abline(v = r_e, col = "red", lty = 1, lwd = 2)
abline(v = r_f, col = "brown", lty = 1, lwd = 2)
if (!is.na(rhat_e)) abline(v = rhat_e, col = "blue", lty = 1, lwd = 2)
if (!is.na(rhat_f)) abline(v = rhat_f, col = "darkgreen", lty = 1, lwd = 2)

par(new = TRUE)
plot(s_x, X_vals, type = "l", col = rgb(0.3, 0.3, 0.3, 0.7), lwd = 1.5,
     axes = FALSE, xlab = "", ylab = "",
     xlim = range(adf_df$s), ylim = x_ylim)

axis(4, at = right_ticks, labels = sci_tick_labels(right_ticks), las = 1)
mtext("Price process Time Series", side = 4, line = 3)

legend("topleft",
       legend = c(expression("True " * r[e]),
                  expression("True " * r[f]),
                  expression(hat(r)[e]),
                  expression(hat(r)[f])),
       col = c("red", "brown", "blue", "darkgreen"),
       lty = c(1, 1, 1, 1),
       lwd = c(2, 2, 2, 2),
       bty = "n")
dev.off()
##################################################################

# Threshold selections citical values


# ============================================================
# 90% upper critical value under H0: delta = 1
# 10% lower critical value under H1: bubble model
# coefficient-based DF statistic
# ============================================================

B       <- 1000
n_grid  <- seq(500, 1000, by = 50)

# Common volatility parameters
d_par   <- 0.01
eta_par <- 0.5

# H1 bubble parameters
r_e     <- 0.3
r_f     <- 0.4
c_par   <- 1
alpha   <- 0.5

# ------------------------------------------------------------
# Generate X_0,...,X_n under H0: delta = 1 for all t
# ------------------------------------------------------------
generate_series_h0 <- function(n, d, eta, x0 = 100, sigma2_0 = 1) {
  phi_n <- 1 - d / log(log(n))
  
  X <- numeric(n + 1)
  X[1] <- x0
  
  logsig2_prev <- log(sigma2_0)
  
  for (t in 1:n) {
    eta_t     <- rnorm(1, 0, eta)
    logsig2_t <- phi_n * logsig2_prev + eta_t
    sigma_t   <- sqrt(exp(logsig2_t))
    u_t       <- sigma_t * rnorm(1)
    
    X[t + 1] <- X[t] + u_t
    logsig2_prev <- logsig2_t
  }
  
  X
}

# ------------------------------------------------------------
# Generate X_0,...,X_n with bubble window [tau_e, tau_f]
# ------------------------------------------------------------
generate_series <- function(n, r_e, r_f, c, alpha, d, eta,
                            x0 = 100, sigma2_0 = 1) {
  tau_e   <- floor(n * r_e)
  tau_f   <- floor(n * r_f)
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
    } else if (t <= tau_f) {
      delta_n
    } else {
      1
    }
    
    X[t + 1] <- a_t * X[t] + u_t
    logsig2_prev <- logsig2_t
  }
  
  X
}

# ------------------------------------------------------------
# Compute coefficient-based DF statistic on full sample
# Regression: X_t on X_{t-1}, t = 1,...,n
# DF_delta = n * (delta_hat - 1)
# ------------------------------------------------------------
df_delta_fullsample <- function(X_full) {
  n <- length(X_full) - 1
  
  Xlag <- X_full[1:n]
  Y    <- X_full[2:(n + 1)]
  
  denom <- sum(Xlag^2)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)
  
  delta_hat <- sum(Y * Xlag) / denom
  
  n * (delta_hat - 1)
}

# ------------------------------------------------------------
# Compute coefficient-based DF statistic at fixed tau
# Regression: X_t on constant + X_{t-1}, t = 1,...,tau
# DF_delta = tau * (delta_hat - 1)
# ------------------------------------------------------------
adf_delta_at_tau <- function(X_full, tau) {
  Xlag <- X_full[1:tau]
  Y    <- X_full[2:(tau + 1)]
  
  Ybar <- mean(Y)
  SSR  <- sum((Xlag - Ybar)^2)
  if (!is.finite(SSR) || SSR <= 0) return(NA_real_)
  
  Sxy <- sum((Xlag - Ybar) * (Y - Ybar))
  delta_hat <- Sxy / SSR
  
  tau * (delta_hat - 1)
}

# ------------------------------------------------------------
# For one n:
# simulate B paths under H0
# compute full-sample DF_delta
# return 90% upper critical value
# ------------------------------------------------------------
critical_value_h0_upper <- function(n, B, d, eta, x0 = 1, sigma2_0 = 1) {
  vals <- numeric(B)
  
  for (b in seq_len(B)) {
    X_full <- generate_series_h0(
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

# ------------------------------------------------------------
# For one n:
# simulate B paths under H1
# evaluate DF_delta at tau = floor(n*r_f)
# return 10% lower critical value
# ------------------------------------------------------------
critical_value_h1_lower <- function(n, B, r_e, r_f, c, alpha, d, eta,
                                    x0 = 1, sigma2_0 = 1) {
  tau_after_rf <- floor(n * r_f)
  if (tau_after_rf > n) stop("tau_after_rf exceeds n")
  
  vals <- numeric(B)
  
  for (b in seq_len(B)) {
    X_full <- generate_series(
      n = n, r_e = r_e, r_f = r_f,
      c = c, alpha = alpha,
      d = d, eta = eta,
      x0 = x0, sigma2_0 = sigma2_0
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

# ------------------------------------------------------------
# Run over n-grid
# ------------------------------------------------------------
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
    "Done n = %d, cv_90 = %.4f, cv_10 = %.4f\n",
    n, crit_90[i], crit_10[i]
  ))
}

# ------------------------------------------------------------
# Table
# ------------------------------------------------------------
cv_table <- data.frame(
  n = n_grid,
  cv_90_upper_H0 = crit_90,
  cv_10_lower_H1 = crit_10
)

print(cv_table)

n <- c(500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000)

cv_h0 <- c(0.7463, 0.8563, 0.8304, 0.0597, 0.8997, 0.7037, 0.9823, 0.7888, 0.9474, 0.8374, 0.7803)
cv_h1 <- c(-0.2566, -0.9916, 1.2746, 1.9418, 1.0183, 2.6310, 3.2163, 3.8318, 4.6070, 3.2394, 4.1479)

ylims <- range(c(cv_h0, cv_h1))

plot(
  n, cv_h0,
  type = "b", pch = 19, lwd = 2,
  col = "firebrick",
  ylim = ylims,
  xlab = "n",
  ylab = "Critical value",
  main = expression(paste("Coefficient-based critical values"))
)

lines(
  n, cv_h1,
  type = "b", pch = 17, lwd = 2,
  col = "steelblue"
)

legend(
  "topleft",
  legend = c(expression(cv[0.10,H[0]]^delta),
             expression(cv[0.10,H[1]]^delta)),
  col = c("firebrick", "steelblue"),
  lty = 1,
  pch = c(19, 17),
  lwd = 2,
  bty = "n"
)

pdf("critical_values_plot.pdf", width = 8, height = 5)

plot(
  n, cv_h0,
  type = "b", pch = 19, lwd = 2,
  col = "firebrick",
  ylim = range(c(cv_h0, cv_h1)),
  xlab = "Sample Size",
  ylab = "Critical value",
  main = expression(paste("Coefficient-based critical values"))
)

lines(
  n, cv_h1,
  type = "b", pch = 17, lwd = 2,
  col = "steelblue"
)



dev.off()
