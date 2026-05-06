################################################################################
# Appendix Crypto.R
#
# Replication code for Online Appendix Figures A.1 and A.2 of:
# "Is There an AI Bubble? Robust Date-Stamping for Periods of Exuberance"
#
# Online Appendix Figure A.1:
#   Simulation evidence on the magnitude of X_{tau_f} - X_{tau_e}, comparing
#   homoskedastic and stochastic-volatility designs.
#
# Online Appendix Figure A.2:
#   Additional cryptocurrency bubble-detection results for:
#     - Cardano
#     - Solana
#     - Binance Coin
#     - Dogecoin
#
# Data source for crypto prices:
# Yahoo Finance daily adjusted closing prices.
#
# Notes:
# - Package installation is intentionally omitted; see README.md.
# - Set SAVE_FIGURES <- TRUE to save PDF and PNG versions.
# - PDF files are suitable for the paper / Online Appendix.
# - PNG files are included for GitHub preview.
################################################################################

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------

library(quantmod)
library(zoo)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)
library(scales)

# ------------------------------------------------------------------------------
# Global settings
# ------------------------------------------------------------------------------

SAVE_FIGURES <- TRUE
SAVE_PNG_PREVIEW <- TRUE

OUTPUT_DIR <- "figures"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

FIGURE_A1_PDF <- file.path(OUTPUT_DIR, "log_gap_const_vs_sv.pdf")
FIGURE_A1_PNG <- file.path(OUTPUT_DIR, "log_gap_const_vs_sv.png")

FIGURE_A2_PDF <- file.path(OUTPUT_DIR, "crypto_four_panel_bubbles.pdf")
FIGURE_A2_PNG <- file.path(OUTPUT_DIR, "crypto_four_panel_bubbles.png")

################################################################################
# Online Appendix Figure A.1
# Simulation evidence for X_{tau_f} - X_{tau_e}
################################################################################

# Simulate the gap X_{tau_f} - X_{tau_e} under homoskedastic innovations.
gap_X_tf_te_const <- function(
    n,
    X0,
    r_e,
    r_f,
    c_par,
    alpha,
    sigma
) {
  tau_e <- floor(n * r_e)
  tau_f <- floor(n * r_f)
  delta_n <- 1 + c_par / (n^alpha)
  
  X <- numeric(n + 1)
  X[1] <- X0
  
  for (t in seq_len(n)) {
    u_t <- sigma * rnorm(1)
    
    a_t <- if (t < tau_e) {
      1
    } else if (t <= tau_f) {
      delta_n
    } else {
      1
    }
    
    X[t + 1] <- a_t * X[t] + u_t
  }
  
  X[tau_f + 1] - X[tau_e + 1]
}

# Simulate the gap X_{tau_f} - X_{tau_e} under persistent stochastic volatility.
gap_X_tf_te_sv <- function(
    n,
    X0,
    r_e,
    r_f,
    c_par,
    alpha,
    d,
    eta,
    sigma2_0 = 1
) {
  tau_e <- floor(n * r_e)
  tau_f <- floor(n * r_f)
  delta_n <- 1 + c_par / (n^alpha)
  phi_n <- 1 - d / log(log(n))
  
  X <- numeric(n + 1)
  X[1] <- X0
  
  log_sigma2_prev <- log(sigma2_0)
  
  for (t in seq_len(n)) {
    eta_t <- rnorm(1, mean = 0, sd = eta)
    log_sigma2_t <- phi_n * log_sigma2_prev + eta_t
    sigma_t <- sqrt(exp(log_sigma2_t))
    u_t <- sigma_t * rnorm(1)
    
    a_t <- if (t < tau_e) {
      1
    } else if (t <= tau_f) {
      delta_n
    } else {
      1
    }
    
    X[t + 1] <- a_t * X[t] + u_t
    log_sigma2_prev <- log_sigma2_t
  }
  
  X[tau_f + 1] - X[tau_e + 1]
}

make_log_gap_comparison_plot <- function(
    seed = 123,
    n_grid = seq(50, 500, by = 50),
    B = 1000,
    X0 = 1,
    r_e = 0.4,
    r_f = 0.6,
    c_par = 1,
    alpha = 0.5,
    sigma = 1,
    d_par = 1,
    eta_par = 0.5,
    sigma2_0 = 1,
    save_figure = SAVE_FIGURES
) {
  set.seed(seed)
  
  # Monte Carlo mean absolute gap under homoskedasticity.
  gap_mean_abs_const <- sapply(n_grid, function(n) {
    mean(abs(replicate(
      B,
      gap_X_tf_te_const(
        n = n,
        X0 = X0,
        r_e = r_e,
        r_f = r_f,
        c_par = c_par,
        alpha = alpha,
        sigma = sigma
      )
    )))
  })
  
  # Monte Carlo mean absolute gap under stochastic volatility.
  gap_mean_abs_sv <- sapply(n_grid, function(n) {
    mean(abs(replicate(
      B,
      gap_X_tf_te_sv(
        n = n,
        X0 = X0,
        r_e = r_e,
        r_f = r_f,
        c_par = c_par,
        alpha = alpha,
        d = d_par,
        eta = eta_par,
        sigma2_0 = sigma2_0
      )
    )))
  })
  
  log_gap_mean_abs_const <- log(pmax(gap_mean_abs_const, 1e-8))
  log_gap_mean_abs_sv <- log(pmax(gap_mean_abs_sv, 1e-8))
  
  plot_data <- list(
    n_grid = n_grid,
    log_gap_mean_abs_const = log_gap_mean_abs_const,
    log_gap_mean_abs_sv = log_gap_mean_abs_sv
  )
  
  y_limits <- range(
    c(log_gap_mean_abs_const, log_gap_mean_abs_sv),
    na.rm = TRUE
  )
  
  draw_plot <- function() {
    old_par <- par(mfrow = c(1, 2), mar = c(4, 6, 3, 1))
    on.exit(par(old_par), add = TRUE)
    
    plot(
      n_grid,
      log_gap_mean_abs_const,
      type = "b",
      pch = 19,
      lwd = 2,
      col = "firebrick",
      ylim = y_limits,
      xlab = "Sample size",
      ylab = expression(log(E * abs(X[tau[f]] - X[tau[e]]))),
      main = "Homoskedasticity"
    )
    
    plot(
      n_grid,
      log_gap_mean_abs_sv,
      type = "b",
      pch = 19,
      lwd = 2,
      col = "steelblue",
      ylim = y_limits,
      xlab = "Sample size",
      ylab = expression(log(E * abs(X[tau[f]] - X[tau[e]]))),
      main = "Stochastic volatility"
    )
  }
  
  draw_plot()
  
  if (isTRUE(save_figure)) {
    grDevices::pdf(FIGURE_A1_PDF, width = 12, height = 5)
    draw_plot()
    grDevices::dev.off()
    
    if (isTRUE(SAVE_PNG_PREVIEW)) {
      grDevices::png(FIGURE_A1_PNG, width = 12, height = 5, units = "in", res = 300)
      draw_plot()
      grDevices::dev.off()
    }
  }
  
  invisible(plot_data)
}

################################################################################
# Online Appendix Figure A.2
# Additional cryptocurrency SV-ADF date-stamping results
################################################################################

# Convert an event date into a month label.
# If the event occurs after the middle of the month, the label is rounded forward.
collapse_month_label <- function(date_value) {
  day_value <- as.integer(format(date_value, "%d"))
  
  if (day_value < 15) {
    format(date_value, "%b-%Y")
  } else {
    format(seq(date_value, by = "month", length.out = 2)[2], "%b-%Y")
  }
}

# Fetch adjusted closing prices from Yahoo Finance.
fetch_adjusted_prices <- function(ticker, start_date, end_date) {
  xts_obj <- quantmod::getSymbols(
    ticker,
    src = "yahoo",
    from = start_date,
    to = end_date,
    auto.assign = FALSE
  )
  
  price_xts <- quantmod::Ad(xts_obj)
  price_xts <- na.omit(price_xts)
  
  if (NROW(price_xts) < 20) {
    stop(paste("Not enough observations for", ticker))
  }
  
  price_xts
}

# Compute the recursive coefficient-based ADF statistic:
#   DF_delta(r) = tau * (delta_hat_tau - 1),
# where tau = floor(nr).
compute_recursive_adf <- function(price_vector, date_vector, min_fraction = 0.10) {
  n <- length(price_vector) - 1
  
  s_grid <- seq(min_fraction, 1, by = 1 / n)
  tau_grid <- floor(n * s_grid)
  
  grid_df <- unique(data.frame(
    s = s_grid,
    tau = tau_grid
  ))
  
  grid_df <- grid_df[grid_df$tau >= 2, ]
  
  if (nrow(grid_df) == 0) {
    stop("No valid recursive windows.")
  }
  
  delta_hat_vec <- sapply(grid_df$tau, function(tau) {
    y_tau <- price_vector[2:(tau + 1)]
    xlag_tau <- price_vector[1:tau]
    
    y_tilde_tau <- y_tau - mean(y_tau)
    xlag_tilde_tau <- xlag_tau - mean(xlag_tau)
    
    denom <- sum(xlag_tilde_tau^2)
    
    if (denom == 0) {
      return(NA_real_)
    }
    
    sum(xlag_tilde_tau * y_tilde_tau) / denom
  })
  
  grid_df$delta_hat <- delta_hat_vec
  grid_df$ADF_delta <- grid_df$tau * (grid_df$delta_hat - 1)
  grid_df$Date <- as.Date(date_vector[grid_df$tau + 1])
  
  ns_vec <- n * grid_df$s
  
  grid_df$boundary_orig <- log(ns_vec) / 10
  grid_df$boundary_screen <- log(ns_vec)
  grid_df$boundary_collapse <- log(ns_vec) / 2
  
  keep <- is.finite(grid_df$ADF_delta) &
    is.finite(grid_df$boundary_orig) &
    is.finite(grid_df$boundary_screen) &
    is.finite(grid_df$boundary_collapse)
  
  grid_df <- grid_df[keep, ]
  
  if (nrow(grid_df) == 0) {
    stop("All recursive ADF statistics are missing.")
  }
  
  grid_df
}

# Detect first bubble origination.
# Origination occurs when the statistic exceeds log(ns)/10 for M consecutive
# recursive windows.
detect_origination <- function(grid_df, M) {
  above_orig <- grid_df$ADF_delta > grid_df$boundary_orig
  N <- nrow(grid_df)
  
  if (N < M) {
    return(NA_integer_)
  }
  
  for (j in seq_len(N - M + 1)) {
    if (all(above_orig[j:(j + M - 1)])) {
      return(j)
    }
  }
  
  NA_integer_
}

# Detect post-bridge collapse.
# Candidate start j must satisfy ADF_delta[j] < log(ns), and then the statistic
# must remain below log(ns)/2 for R consecutive recursive windows.
find_post_bridge_collapse <- function(grid_df, idx_start, R) {
  N <- nrow(grid_df)
  
  if (is.na(idx_start) || idx_start > N) {
    return(list(
      idx = NA_integer_,
      date = as.Date(NA),
      label = NA_character_
    ))
  }
  
  last_possible_start <- N - R + 1
  
  if (last_possible_start < idx_start) {
    return(list(
      idx = NA_integer_,
      date = as.Date(NA),
      label = NA_character_
    ))
  }
  
  for (j in idx_start:last_possible_start) {
    if (grid_df$ADF_delta[j] < grid_df$boundary_screen[j]) {
      run_idx <- j:(j + R - 1)
      
      if (all(grid_df$ADF_delta[run_idx] < grid_df$boundary_collapse[run_idx])) {
        collapse_date <- grid_df$Date[j]
        
        return(list(
          idx = j,
          date = collapse_date,
          label = collapse_month_label(collapse_date)
        ))
      }
    }
  }
  
  list(
    idx = NA_integer_,
    date = as.Date(NA),
    label = NA_character_
  )
}

# Display-only adjustment for short-lived pre-origination spikes.
#
# Important:
# This does NOT affect bubble detection or date-stamping. It only rescales
# transient pre-origination spikes in the plotted statistic to make the final
# figure visually readable.
make_display_adf_series <- function(grid_df, orig_idx, compression_factor = 6) {
  grid_df$ADF_delta_plot <- grid_df$ADF_delta
  
  if (is.na(orig_idx) || orig_idx <= 1) {
    return(grid_df)
  }
  
  pre_idx <- seq_len(orig_idx - 1)
  above_orig_pre <- grid_df$ADF_delta[pre_idx] > grid_df$boundary_orig[pre_idx]
  above_noise_pre <- grid_df$ADF_delta[pre_idx] > grid_df$boundary_orig[pre_idx]
  
  if (!any(above_orig_pre)) {
    return(grid_df)
  }
  
  runs <- rle(above_orig_pre)
  run_end <- cumsum(runs$lengths)
  run_start <- run_end - runs$lengths + 1
  true_runs <- which(runs$values)
  
  compress_idx <- integer(0)
  
  for (rr in true_runs) {
    start_i <- run_start[rr]
    end_i <- run_end[rr]
    
    expanded_start <- start_i
    while (expanded_start > 1 && above_noise_pre[expanded_start - 1]) {
      expanded_start <- expanded_start - 1
    }
    
    expanded_end <- end_i
    while (expanded_end < length(pre_idx) && above_noise_pre[expanded_end + 1]) {
      expanded_end <- expanded_end + 1
    }
    
    compress_idx <- union(compress_idx, pre_idx[expanded_start:expanded_end])
  }
  
  if (length(compress_idx) > 0) {
    grid_df$ADF_delta_plot[compress_idx] <-
      grid_df$ADF_delta_plot[compress_idx] / compression_factor
  }
  
  grid_df
}

# Detect all bubble episodes after the first one.
# This is used only to print an interpretable console summary.
detect_bubble_episodes <- function(grid_df, M, R, bridge_days) {
  above <- grid_df$ADF_delta > grid_df$boundary_orig
  N <- nrow(grid_df)
  
  episodes <- list()
  i <- 1
  episode_id <- 1
  
  while (i <= N) {
    orig_idx <- NA_integer_
    
    if (N >= M && i <= N - M + 1) {
      for (j in i:(N - M + 1)) {
        if (all(above[j:(j + M - 1)])) {
          orig_idx <- j
          break
        }
      }
    }
    
    if (is.na(orig_idx)) {
      break
    }
    
    orig_date <- grid_df$Date[orig_idx]
    bridge_end_date <- orig_date + bridge_days
    idx_bridge_end <- which(grid_df$Date >= bridge_end_date)[1]
    
    if (is.na(idx_bridge_end)) {
      idx_bridge_end <- N
    }
    
    collapse_idx <- NA_integer_
    collapse_date <- as.Date(NA)
    collapse_label <- NA_character_
    collapse_regime <- NA_character_
    
    # Collapse within bridge: ADF < log(ns)/10.
    j_start <- orig_idx + 1
    j_end <- idx_bridge_end - 1
    
    if (j_start <= j_end) {
      for (jj in j_start:j_end) {
        if (grid_df$ADF_delta[jj] < grid_df$boundary_orig[jj]) {
          collapse_idx <- jj
          collapse_date <- grid_df$Date[jj]
          collapse_label <- collapse_month_label(collapse_date)
          collapse_regime <- "bridge"
          break
        }
      }
    }
    
    # Collapse post-bridge: log(ns) screen and R-run below log(ns)/2.
    if (is.na(collapse_idx) && idx_bridge_end <= N) {
      tmp <- find_post_bridge_collapse(grid_df, idx_bridge_end, R)
      
      if (!is.na(tmp$idx)) {
        collapse_idx <- tmp$idx
        collapse_date <- tmp$date
        collapse_label <- tmp$label
        collapse_regime <- "post_bridge"
      }
    }
    
    episodes[[episode_id]] <- data.frame(
      episode = episode_id,
      orig_idx = orig_idx,
      orig_date = grid_df$Date[orig_idx],
      orig_label = format(grid_df$Date[orig_idx], "%b-%Y"),
      orig_s = grid_df$s[orig_idx],
      orig_tau = grid_df$tau[orig_idx],
      orig_adf = grid_df$ADF_delta[orig_idx],
      orig_boundary = grid_df$boundary_orig[orig_idx],
      collapse_idx = ifelse(is.na(collapse_idx), NA, collapse_idx),
      collapse_date = if (is.na(collapse_idx)) as.character(NA) else as.character(collapse_date),
      collapse_label = if (is.na(collapse_idx)) NA_character_ else collapse_label,
      collapse_s = ifelse(is.na(collapse_idx), NA, grid_df$s[collapse_idx]),
      collapse_tau = ifelse(is.na(collapse_idx), NA, grid_df$tau[collapse_idx]),
      collapse_adf = ifelse(is.na(collapse_idx), NA, grid_df$ADF_delta[collapse_idx]),
      collapse_boundary = ifelse(
        is.na(collapse_idx),
        NA,
        if (identical(collapse_regime, "bridge")) {
          grid_df$boundary_orig[collapse_idx]
        } else {
          grid_df$boundary_collapse[collapse_idx]
        }
      ),
      collapse_regime = ifelse(is.na(collapse_idx), NA, collapse_regime),
      stringsAsFactors = FALSE
    )
    
    if (is.na(collapse_idx)) {
      break
    }
    
    i <- collapse_idx + 1
    episode_id <- episode_id + 1
  }
  
  if (length(episodes) == 0) {
    return(NULL)
  }
  
  do.call(rbind, episodes)
}

# Print a compact summary of detected episodes.
print_episode_summary <- function(company_name, episodes_df, M, R, bridge_days) {
  if (is.null(episodes_df) || nrow(episodes_df) == 0) {
    cat("No bubble episodes found for", company_name, "\n\n")
    return(invisible(NULL))
  }
  
  cat(
    "Detected bubble episodes for", company_name,
    "with M =", M,
    ", R =", R,
    "and bridge_days =", bridge_days, "\n\n"
  )
  
  for (k in seq_len(nrow(episodes_df))) {
    cat("Episode", episodes_df$episode[k], "\n")
    cat("  Origination Month :", episodes_df$orig_label[k], "\n")
    cat("  Origination date  :", as.character(episodes_df$orig_date[k]), "\n")
    
    if (!is.na(episodes_df$collapse_idx[k])) {
      cat("  Collapse regime   :", episodes_df$collapse_regime[k], "\n")
      cat("  Collapse Month    :", episodes_df$collapse_label[k], "\n")
      cat("  Collapse date     :", as.character(episodes_df$collapse_date[k]), "\n")
    } else {
      cat("  Collapse          : not found\n")
    }
    
    cat("\n")
  }
  
  invisible(NULL)
}

# Main plotting function for one cryptocurrency.
plot_bubble_episodes_horizontal <- function(
    ticker,
    company_name,
    start_date,
    end_date,
    L = 20,
    M = 5,
    R = 30,
    bridge_days = 60
) {
  # L is retained for compatibility with earlier scripts.
  # The current implementation uses the coefficient-based recursive statistic.
  
  price_xts <- fetch_adjusted_prices(ticker, start_date, end_date)
  
  X <- as.numeric(price_xts)
  dates <- zoo::index(price_xts)
  
  grid_df <- compute_recursive_adf(X, dates, min_fraction = 0.10)
  
  # --------------------------------------------------------------------------
  # Date-stamping step
  # --------------------------------------------------------------------------
  
  orig_idx <- detect_origination(grid_df, M = M)
  
  collapse_within_bridge <- FALSE
  collapse_bridge_date <- as.Date(NA)
  collapse_bridge_label <- NA_character_
  collapse_idx <- NA_integer_
  collapse_date <- as.Date(NA)
  collapse_label <- NA_character_
  collapse_regime <- NA_character_
  
  if (!is.na(orig_idx)) {
    bridge_start_date <- grid_df$Date[orig_idx]
    bridge_end_date <- bridge_start_date + bridge_days
    
    idx_bridge_end <- which(grid_df$Date >= bridge_end_date)[1]
    
    if (is.na(idx_bridge_end)) {
      idx_bridge_end <- nrow(grid_df)
    }
    
    tau_plot <- 1:grid_df$tau[orig_idx]
    orig_curve_df <- data.frame(
      Date = as.Date(dates[tau_plot + 1]),
      y = log(tau_plot) / 10
    )
    
    collapse_curve_start_date <- bridge_start_date + bridge_days + R
    idx_collapse_curve_start <- which(grid_df$Date >= collapse_curve_start_date)[1]
    
    collapse_curve_df <- if (!is.na(idx_collapse_curve_start) &&
                             idx_collapse_curve_start <= nrow(grid_df)) {
      data.frame(
        Date = grid_df$Date[idx_collapse_curve_start:nrow(grid_df)],
        y = grid_df$boundary_collapse[idx_collapse_curve_start:nrow(grid_df)]
      )
    } else {
      data.frame(Date = as.Date(character()), y = numeric())
    }
    
    # Collapse within bridge: first crossing below log(ns)/10.
    j_start <- orig_idx + 1
    j_end <- idx_bridge_end - 1
    
    if (j_start <= j_end) {
      for (jj in j_start:j_end) {
        if (grid_df$ADF_delta[jj] < grid_df$boundary_orig[jj]) {
          collapse_within_bridge <- TRUE
          collapse_idx <- jj
          collapse_date <- grid_df$Date[jj]
          collapse_bridge_date <- collapse_date
          collapse_bridge_label <- collapse_month_label(collapse_date)
          collapse_label <- collapse_bridge_label
          collapse_regime <- "bridge"
          break
        }
      }
    }
    
    # Collapse post-bridge: log(ns) screen and R-run below log(ns)/2.
    if (is.na(collapse_idx) && idx_bridge_end <= nrow(grid_df)) {
      post_bridge_res <- find_post_bridge_collapse(grid_df, idx_bridge_end, R)
      
      if (!is.na(post_bridge_res$idx)) {
        collapse_idx <- post_bridge_res$idx
        collapse_date <- post_bridge_res$date
        collapse_label <- post_bridge_res$label
        collapse_regime <- "post_bridge"
      }
    }
  } else {
    bridge_start_date <- as.Date(NA)
    bridge_end_date <- as.Date(NA)
    
    tau_plot <- 1:max(grid_df$tau)
    orig_curve_df <- data.frame(
      Date = as.Date(dates[tau_plot + 1]),
      y = log(tau_plot) / 10
    )
    
    collapse_curve_df <- data.frame(Date = as.Date(character()), y = numeric())
  }
  
  # Display-only transformation.
  grid_df <- make_display_adf_series(
    grid_df = grid_df,
    orig_idx = orig_idx,
    compression_factor = 6
  )
  
  orig_curve_df$y_plot <- orig_curve_df$y
  collapse_curve_df$y_plot <- collapse_curve_df$y
  
  price_df <- data.frame(
    Date = as.Date(dates),
    Price = X,
    stringsAsFactors = FALSE
  )
  
  # Colors match the main paper / appendix convention.
  adf_col <- "#1f77b4"
  boundary_col <- "#ff7f0e"
  boundary2_col <- "#8c564b"
  price_col <- "#d62728"
  orig_col <- "#2ca02c"
  collapse_col <- "#9467bd"
  
  left_min_raw <- min(c(
    grid_df$ADF_delta_plot,
    orig_curve_df$y_plot,
    if (nrow(collapse_curve_df) > 0) collapse_curve_df$y_plot else Inf
  ), na.rm = TRUE)
  
  left_max_raw <- max(c(
    grid_df$ADF_delta_plot,
    orig_curve_df$y_plot,
    if (nrow(collapse_curve_df) > 0) collapse_curve_df$y_plot else -Inf
  ), na.rm = TRUE)
  
  if (!is.finite(left_min_raw) ||
      !is.finite(left_max_raw) ||
      left_max_raw == left_min_raw) {
    left_min_raw <- min(grid_df$ADF_delta_plot, na.rm = TRUE)
    left_max_raw <- max(grid_df$ADF_delta_plot, na.rm = TRUE) + 1
  }
  
  left_range <- left_max_raw - left_min_raw
  left_min_plot <- left_min_raw - 0.10 * left_range
  left_max_plot <- left_max_raw + 0.20 * left_range
  
  right_min <- min(price_df$Price, na.rm = TRUE)
  right_max <- max(price_df$Price, na.rm = TRUE)
  
  if (right_max == right_min) {
    right_max <- right_min + 1
  }
  
  # Keep price curve in the lower portion of the left-axis panel.
  price_band_low <- left_min_plot + 0.01 * (left_max_plot - left_min_plot)
  price_band_high <- left_min_plot + 0.38 * (left_max_plot - left_min_plot)
  
  price_df$Price_scaled <- price_band_low +
    (price_df$Price - right_min) *
    (price_band_high - price_band_low) /
    (right_max - right_min)
  
  price_inverse <- function(z) {
    right_min +
      (z - price_band_low) *
      (right_max - right_min) /
      (price_band_high - price_band_low)
  }
  
  date_span_days <- as.integer(max(price_df$Date) - min(price_df$Date))
  label_offset_left <- max(10, round(0.035 * date_span_days))
  label_offset_right <- max(10, round(0.035 * date_span_days))
  
  left_breaks <- pretty(c(left_min_plot, left_max_plot), n = 6)
  right_breaks_vals <- pretty(c(right_min, right_max), n = 5)
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = grid_df,
      ggplot2::aes(x = Date, y = ADF_delta_plot),
      color = adf_col,
      linewidth = 1.15
    ) +
    ggplot2::geom_line(
      data = orig_curve_df,
      ggplot2::aes(x = Date, y = y_plot),
      color = boundary_col,
      linetype = "solid",
      linewidth = 1.15,
      alpha = 0.95
    ) +
    {
      if (nrow(collapse_curve_df) > 0) {
        ggplot2::geom_line(
          data = collapse_curve_df,
          ggplot2::aes(x = Date, y = y_plot),
          color = boundary2_col,
          linetype = "solid",
          linewidth = 1.15,
          alpha = 0.95
        )
      }
    } +
    ggplot2::geom_line(
      data = price_df,
      ggplot2::aes(x = Date, y = Price_scaled),
      color = price_col,
      linewidth = 1.10,
      alpha = 0.95
    ) +
    ggplot2::scale_x_date(
      date_breaks = "3 months",
      date_labels = "%Y %b"
    ) +
    ggplot2::scale_y_continuous(
      name = NULL,
      limits = c(left_min_plot, left_max_plot),
      breaks = left_breaks,
      labels = scales::number_format(accuracy = 0.1),
      sec.axis = ggplot2::sec_axis(
        ~ price_inverse(.),
        breaks = right_breaks_vals,
        labels = right_breaks_vals,
        name = NULL
      )
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      title = company_name
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(
        size = 26,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 10)
      ),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.title.y.right = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 14, angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 14, color = adf_col),
      axis.text.y.right = ggplot2::element_text(size = 14, color = price_col),
      axis.ticks.y = ggplot2::element_line(color = adf_col),
      axis.ticks.y.right = ggplot2::element_line(color = price_col),
      axis.line.y.left = ggplot2::element_line(color = adf_col),
      axis.line.y.right = ggplot2::element_line(color = price_col),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(12, 55, 12, 55)
    )
  
  if (!is.na(orig_idx)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = as.numeric(grid_df$Date[orig_idx]),
        linetype = "solid",
        linewidth = 1.1,
        color = orig_col
      ) +
      ggplot2::annotate(
        "label",
        x = grid_df$Date[orig_idx] - label_offset_left,
        y = left_max_plot - 0.055 * (left_max_plot - left_min_plot),
        label = paste0("Origination — ", format(grid_df$Date[orig_idx], "%b %Y")),
        hjust = 1,
        vjust = 1,
        size = 4.8,
        fontface = "bold",
        color = orig_col,
        fill = "white",
        label.size = 0.25,
        label.r = grid::unit(0.18, "lines"),
        label.padding = grid::unit(0.20, "lines")
      )
  }
  
  if (!is.na(collapse_idx)) {
    collapse_label_x <- collapse_date + label_offset_right
    collapse_label_y <- left_max_plot - 0.135 * (left_max_plot - left_min_plot)
    
    p <- p +
      ggplot2::geom_vline(
        xintercept = as.numeric(collapse_date),
        linetype = "solid",
        linewidth = 1.1,
        color = collapse_col
      ) +
      ggplot2::annotate(
        "label",
        x = collapse_label_x,
        y = collapse_label_y,
        label = paste0("Collapse — ", collapse_label),
        hjust = 0,
        vjust = 1,
        size = 4.6,
        fontface = "bold",
        color = collapse_col,
        fill = "white",
        label.size = 0.25,
        label.r = grid::unit(0.18, "lines"),
        label.padding = grid::unit(0.20, "lines")
      )
  }
  
  episodes_df <- detect_bubble_episodes(
    grid_df = grid_df,
    M = M,
    R = R,
    bridge_days = bridge_days
  )
  
  print_episode_summary(
    company_name = company_name,
    episodes_df = episodes_df,
    M = M,
    R = R,
    bridge_days = bridge_days
  )
  
  invisible(list(
    grid_df = grid_df,
    episodes_df = episodes_df,
    plot = p,
    bridge_info = list(
      r_hat_e = if (!is.na(orig_idx)) grid_df$s[orig_idx] else NA_real_,
      bridge_days = bridge_days,
      R = R,
      bridge_start_date = if (!is.na(orig_idx)) grid_df$Date[orig_idx] else as.Date(NA),
      bridge_end_date = if (!is.na(orig_idx)) grid_df$Date[orig_idx] + bridge_days else as.Date(NA),
      collapse_within_bridge = collapse_within_bridge,
      collapse_bridge_date = collapse_bridge_date,
      collapse_bridge_label = collapse_bridge_label,
      collapse_idx = collapse_idx,
      collapse_date = collapse_date,
      collapse_label = collapse_label,
      collapse_regime = collapse_regime
    )
  ))
}

# Run the plotting function while suppressing individual plot rendering.
get_bubble_result <- function(...) {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp, width = 10, height = 8)
  
  out <- NULL
  
  tryCatch(
    {
      invisible(utils::capture.output({
        out <- plot_bubble_episodes_horizontal(...)
      }))
    },
    finally = {
      grDevices::dev.off()
      
      if (file.exists(tmp)) {
        unlink(tmp)
      }
    }
  )
  
  out
}

# Optional helper to manually nudge origination labels in selected panels.
# This affects display only and has no effect on date-stamping.
nudge_origination_label <- function(p, dy = 0.04) {
  for (i in seq_along(p$layers)) {
    layer_data <- p$layers[[i]]$data
    
    if (is.null(layer_data) || !is.data.frame(layer_data)) {
      next
    }
    
    if (!("label" %in% names(layer_data)) || !("y" %in% names(layer_data))) {
      next
    }
    
    idx <- grepl("origination", layer_data$label, ignore.case = TRUE)
    
    if (!any(idx)) {
      next
    }
    
    y_range <- diff(range(layer_data$y, na.rm = TRUE))
    
    if (!is.finite(y_range) || y_range == 0) {
      y_range <- 1
    }
    
    layer_data$y[idx] <- layer_data$y[idx] + dy * y_range
    p$layers[[i]]$data <- layer_data
  }
  
  p
}

make_additional_crypto_plot <- function(save_figure = SAVE_FIGURES) {
  
  # Cardano.
  out_ada <- get_bubble_result(
    ticker = "ADA-USD",
    company_name = "Cardano",
    start_date = as.Date("2020-01-01"),
    end_date = as.Date("2023-01-01"),
    L = 60,
    M = 60,
    R = 30,
    bridge_days = 60
  )
  
  # Solana.
  out_sol <- get_bubble_result(
    ticker = "SOL-USD",
    company_name = "Solana",
    start_date = as.Date("2020-01-01"),
    end_date = as.Date("2023-01-01"),
    L = 60,
    M = 60,
    R = 30,
    bridge_days = 60
  )
  
  # Binance Coin.
  out_bnb <- get_bubble_result(
    ticker = "BNB-USD",
    company_name = "Binance Coin",
    start_date = as.Date("2020-04-01"),
    end_date = as.Date("2023-01-01"),
    L = 10,
    M = 30,
    R = 10,
    bridge_days = 50
  )
  
  # Dogecoin.
  out_doge <- get_bubble_result(
    ticker = "DOGE-USD",
    company_name = "Dogecoin",
    start_date = as.Date("2020-01-01"),
    end_date = as.Date("2023-01-01"),
    L = 60,
    M = 30,
    R = 30,
    bridge_days = 60
  )
  
  p_ada <- out_ada$plot +
    scale_x_date(date_breaks = "6 months", date_labels = "%Y %b") +
    ggtitle("Cardano") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  
  p_sol <- out_sol$plot +
    scale_x_date(date_breaks = "6 months", date_labels = "%Y %b") +
    ggtitle("Solana") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  
  p_bnb <- out_bnb$plot +
    scale_x_date(date_breaks = "6 months", date_labels = "%Y %b") +
    ggtitle("Binance Coin") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  
  p_doge <- out_doge$plot +
    scale_x_date(date_breaks = "6 months", date_labels = "%Y %b") +
    ggtitle("Dogecoin") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))
  
  # Display-only label adjustments.
  p_ada <- nudge_origination_label(p_ada, dy = 10)
  p_bnb <- nudge_origination_label(p_bnb, dy = -0.4)
  
  # Layout:
  # left column: Cardano, Dogecoin
  # right column: Solana, Binance Coin
  middle_panel <- arrangeGrob(
    grobs = list(p_ada, p_doge, p_sol, p_bnb),
    layout_matrix = rbind(
      c(1, 3),
      c(2, 4)
    )
  )
  
  left_label <- textGrob(
    "Recursive SV-ADF Statistic",
    rot = 90,
    gp = gpar(col = "#1f77b4", fontsize = 25, fontface = "bold")
  )
  
  right_label <- textGrob(
    "Crypto-coin Price",
    rot = 270,
    gp = gpar(col = "#d62728", fontsize = 25, fontface = "bold")
  )
  
  figure_a2 <- arrangeGrob(
    grobs = list(left_label, middle_panel, right_label),
    ncol = 3,
    widths = c(0.04, 0.92, 0.04)
  )
  
  grid.newpage()
  grid.draw(figure_a2)
  
  if (isTRUE(save_figure)) {
    grDevices::pdf(FIGURE_A2_PDF, width = 18, height = 16)
    grid.newpage()
    grid.draw(figure_a2)
    grDevices::dev.off()
    
    if (isTRUE(SAVE_PNG_PREVIEW)) {
      grDevices::png(FIGURE_A2_PNG, width = 18, height = 16, units = "in", res = 300)
      grid.newpage()
      grid.draw(figure_a2)
      grDevices::dev.off()
    }
  }
  
  invisible(list(
    Cardano = out_ada,
    Solana = out_sol,
    Binance_Coin = out_bnb,
    Dogecoin = out_doge,
    plot = figure_a2
  ))
}

################################################################################
# Run all Online Appendix outputs
################################################################################

figure_a1_output <- make_log_gap_comparison_plot(
  seed = 123,
  n_grid = seq(50, 500, by = 50),
  B = 1000,
  X0 = 1,
  r_e = 0.4,
  r_f = 0.6,
  c_par = 1,
  alpha = 0.5,
  sigma = 1,
  d_par = 1,
  eta_par = 0.5,
  sigma2_0 = 1,
  save_figure = SAVE_FIGURES
)

figure_a2_output <- make_additional_crypto_plot(
  save_figure = SAVE_FIGURES
)

################################################################################
# End of script
################################################################################
