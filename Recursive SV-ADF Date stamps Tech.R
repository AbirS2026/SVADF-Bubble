################################################################################
# Recursive SV-ADF Date Stamps: Large Technology Firms
#
# Replication code for Figure 7 of:
# "Is There an AI Bubble? Robust Date-Stamping for Periods of Exuberance"
#
# Figure 7: Recursive SV-ADF statistics and stock prices for six large
# technology firms:
#   Alphabet, Meta, Tesla, Apple, Microsoft, and Amazon.
#
# Data source:
# Yahoo Finance daily adjusted closing prices.
#
# Notes:
# - Package installation is intentionally omitted; see README.md.
# - Set SAVE_FIGURES <- TRUE to save the combined figure.
# - The empirical date windows are fixed for reproducibility.
# - The algorithmic date-stamping rule is kept separate from display-only
#   rescaling used to make the panels visually readable.
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

OUTPUT_PDF <- file.path(OUTPUT_DIR, "six_stock_bubbles_tech.pdf")
OUTPUT_PNG <- file.path(OUTPUT_DIR, "six_stock_bubbles_tech.png")

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

# Convert a date into the month label used for collapse reporting.
# If the collapse occurs after the middle of the month, the label is rounded
# forward to the next month.
collapse_month_label <- function(date_value) {
  day_value <- as.integer(format(date_value, "%d"))
  
  if (day_value < 15) {
    format(date_value, "%b-%Y")
  } else {
    format(seq(date_value, by = "month", length.out = 2)[2], "%b-%Y")
  }
}

# Fetch adjusted closing prices from Yahoo Finance.
# The function stops if Yahoo does not return usable data.
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
  
  # SV-ADF thresholds used in the paper.
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

# Detect the first bubble origination date.
# Origination is identified when the recursive statistic exceeds log(ns)/10
# for M consecutive recursive windows.
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
# A candidate collapse must first pass a screen, ADF_delta < log(ns), and then
# remain below log(ns)/2 for R consecutive recursive windows.
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

# Apply a display-only adjustment to pre-origination spikes.
#
# Important:
# This does NOT affect bubble detection. It only rescales short-lived
# pre-origination spikes in the plotted ADF statistic so that the final figure
# remains readable when a panel contains very large transient excursions.
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
    
    # Collapse within the bridge period: ADF < log(ns)/10.
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
    
    # Collapse after the bridge period: log(ns) screen plus R-run below log(ns)/2.
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

################################################################################
# Main plotting function
################################################################################

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
    
    # Threshold shown before origination: log(tau)/10.
    tau_plot <- 1:grid_df$tau[orig_idx]
    orig_curve_df <- data.frame(
      Date = as.Date(dates[tau_plot + 1]),
      y = log(tau_plot) / 10
    )
    
    # Collapse threshold is displayed after the bridge plus R-period confirmation.
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
    
    # Collapse inside the bridge: first crossing below log(ns)/10.
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
    
    # Collapse after the bridge: log(ns) screen plus R-run below log(ns)/2.
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
  
  # --------------------------------------------------------------------------
  # Display-only transformations
  # --------------------------------------------------------------------------
  
  grid_df <- make_display_adf_series(
    grid_df = grid_df,
    orig_idx = orig_idx,
    compression_factor = 6
  )
  
  orig_curve_df$y_plot <- orig_curve_df$y
  collapse_curve_df$y_plot <- collapse_curve_df$y
  
  # --------------------------------------------------------------------------
  # Price series for right-hand axis
  # --------------------------------------------------------------------------
  
  price_df <- data.frame(
    Date = as.Date(dates),
    Price = X,
    stringsAsFactors = FALSE
  )
  
  # Colors match the main paper convention.
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
  
  # Keep the stock-price curve in the lower part of the left-axis panel.
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
  
  # --------------------------------------------------------------------------
  # Plot
  # --------------------------------------------------------------------------
  
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
  
  # Mark origination.
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
  
  # Mark collapse.
  if (!is.na(collapse_idx)) {
    collapse_label_x <- collapse_date + label_offset_right
    collapse_label_y <- left_max_plot - 0.135 * (left_max_plot - left_min_plot)
    
    # Meta-specific label adjustment to avoid overlap in the final figure.
    if (identical(company_name, "Meta")) {
      collapse_label_x <- collapse_date + 2.2 * label_offset_right
      collapse_label_y <- left_max_plot - 0.08 * (left_max_plot - left_min_plot)
    }
    
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
# This is useful because Figure 7 is a combined six-panel figure.
get_bubble_result <- function(...) {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp, width = 10, height = 6)
  
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

################################################################################
# Generate individual company plots
################################################################################

common_settings <- list(
  L = 60,
  M = 60,
  R = 30,
  bridge_days = 90
)

out_goog <- get_bubble_result(
  ticker = "GOOGL",
  company_name = "Alphabet",
  start_date = as.Date("2022-01-01"),
  end_date = as.Date("2026-02-01"),
  L = common_settings$L,
  M = common_settings$M,
  R = common_settings$R,
  bridge_days = common_settings$bridge_days
)

out_meta <- get_bubble_result(
  ticker = "META",
  company_name = "Meta",
  start_date = as.Date("2022-03-01"),
  end_date = as.Date("2026-01-01"),
  L = common_settings$L,
  M = common_settings$M,
  R = common_settings$R,
  bridge_days = common_settings$bridge_days
)

out_tsla <- get_bubble_result(
  ticker = "TSLA",
  company_name = "Tesla",
  start_date = as.Date("2018-01-01"),
  end_date = as.Date("2022-01-01"),
  L = common_settings$L,
  M = common_settings$M,
  R = common_settings$R,
  bridge_days = common_settings$bridge_days
)

out_aapl <- get_bubble_result(
  ticker = "AAPL",
  company_name = "Apple",
  start_date = as.Date("2022-01-01"),
  end_date = as.Date("2026-02-01"),
  L = common_settings$L,
  M = common_settings$M,
  R = common_settings$R,
  bridge_days = common_settings$bridge_days
)

out_msft <- get_bubble_result(
  ticker = "MSFT",
  company_name = "Microsoft",
  start_date = as.Date("2022-01-01"),
  end_date = as.Date("2026-03-01"),
  L = common_settings$L,
  M = common_settings$M,
  R = common_settings$R,
  bridge_days = common_settings$bridge_days
)

out_amzn <- get_bubble_result(
  ticker = "AMZN",
  company_name = "Amazon",
  start_date = as.Date("2022-01-01"),
  end_date = as.Date("2026-02-01"),
  L = common_settings$L,
  M = common_settings$M,
  R = common_settings$R,
  bridge_days = common_settings$bridge_days
)

################################################################################
# Construct Figure 7 layout
################################################################################

p_goog <- out_goog$plot +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y %b") +
  ggtitle("Alphabet") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))

p_meta <- out_meta$plot +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y %b") +
  ggtitle("Meta") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))

p_tsla <- out_tsla$plot +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y %b") +
  ggtitle("Tesla") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))

p_aapl <- out_aapl$plot +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y %b") +
  ggtitle("Apple") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))

p_msft <- out_msft$plot +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y %b") +
  ggtitle("Microsoft") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))

p_amzn <- out_amzn$plot +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y %b") +
  ggtitle("Amazon") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 18))

# Left column: Alphabet, Meta, Tesla.
left_column <- arrangeGrob(
  grobs = list(p_goog, p_meta, p_tsla),
  ncol = 1
)

# Right column: Apple, Microsoft, Amazon.
right_column <- arrangeGrob(
  grobs = list(p_aapl, p_msft, p_amzn),
  ncol = 1
)

middle_panel <- arrangeGrob(
  grobs = list(left_column, right_column),
  ncol = 2
)

left_axis_label <- textGrob(
  "Recursive SV-ADF Statistic",
  rot = 90,
  gp = gpar(col = "#1f77b4", fontsize = 25, fontface = "bold")
)

right_axis_label <- textGrob(
  "Stock Price",
  rot = 270,
  gp = gpar(col = "#d62728", fontsize = 25, fontface = "bold")
)

figure7 <- arrangeGrob(
  grobs = list(left_axis_label, middle_panel, right_axis_label),
  ncol = 3,
  widths = c(0.04, 0.92, 0.04)
)

################################################################################
# Display and save Figure 7
################################################################################

grid.newpage()
grid.draw(figure7)

if (isTRUE(SAVE_FIGURES)) {
  grDevices::pdf(OUTPUT_PDF, width = 16, height = 20)
  grid.newpage()
  grid.draw(figure7)
  grDevices::dev.off()
}

if (isTRUE(SAVE_FIGURES) && isTRUE(SAVE_PNG_PREVIEW)) {
  grDevices::png(OUTPUT_PNG, width = 16, height = 20, units = "in", res = 300)
  grid.newpage()
  grid.draw(figure7)
  grDevices::dev.off()
}

################################################################################
# End of script
################################################################################
