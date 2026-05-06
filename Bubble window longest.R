plot_bubble_episodes_longest <- function(
    ticker = "TSLA",
    company_name = "TESLA",
    search_start = as.Date("2018-01-01"),
    search_end   = as.Date("2026-01-01"),
    M = 60,                         # consecutive points above log(ns)/10 for origination
    R = 30,                         # consecutive points below log(ns)/2 to confirm post-bridge collapse
    bridge_days = 90,               # calendar days to delay switching after origination
    horizons_years = c(1, 2, 3, 4), # window lengths to scan (years)
    step_months = 1,                # step size for rolling windows (months)
    plot_best = TRUE                # plot only the single global best (longest bubble) window
) {
  # -----------------------------
  # 1) Download full data once
  # -----------------------------
  xts_obj <- quantmod::getSymbols(
    ticker, src = "yahoo", from = search_start, to = search_end, auto.assign = FALSE
  )
  price_xts_full <- quantmod::Ad(xts_obj)
  price_xts_full <- na.omit(price_xts_full)
  if (NROW(price_xts_full) < 50) stop("Not enough observations in the chosen search range.")
  
  # --------------------------------------------------------
  # 2) Analyze one window using log(ns)/10 and log(ns)/2 only
  #    - Origination: M-run above log(ns)/10
  #    - Collapse within bridge: ADF < log(ns)/10 between orig and (orig+bridge_days)
  #    - Post-bridge collapse: first j >= bridge_end with an R-run below log(ns)/2
  #    - Plot thresholds as:
  #         log(ns)/10 up to origination (inclusive)
  #         log(ns)/2 only after bridge_end + 30 days (no connecting line)
  # --------------------------------------------------------
  analyze_window <- function(window_xts) {
    window_xts <- na.omit(window_xts)
    X <- as.numeric(window_xts)
    dates <- zoo::index(window_xts)
    n <- length(X) - 1
    if (n < 20) return(NULL)
    
    # Grid
    s_grid <- seq(0.1, 1, by = 1 / n)
    tau_grid <- floor(n * s_grid)
    grid_df <- unique(data.frame(s = s_grid, tau = tau_grid))
    grid_df <- grid_df[grid_df$tau >= 2, ]
    if (nrow(grid_df) == 0) return(NULL)
    
    # Recursive DF
    delta_hat_vec <- sapply(grid_df$tau, function(tau) {
      y_tau <- X[2:(tau + 1)]
      xlag_tau <- X[1:tau]
      y_tilde_tau <- y_tau - mean(y_tau)
      xlag_tilde_tau <- xlag_tau - mean(xlag_tau)
      denom <- sum(xlag_tilde_tau^2)
      if (denom == 0) return(NA_real_)
      sum(xlag_tilde_tau * y_tilde_tau) / denom
    })
    grid_df$delta_hat <- delta_hat_vec
    grid_df$ADF_delta <- grid_df$tau * (grid_df$delta_hat - 1)
    grid_df$Date <- as.Date(dates[grid_df$tau + 1])
    
    # Thresholds
    ns_vec <- n * grid_df$s
    grid_df$boundary_orig <- log(ns_vec) / 10
    grid_df$boundary_collapse <- log(ns_vec) / 2
    
    keep <- is.finite(grid_df$ADF_delta) &
      is.finite(grid_df$boundary_orig) &
      is.finite(grid_df$boundary_collapse)
    grid_df <- grid_df[keep, ]
    if (nrow(grid_df) == 0) return(NULL)
    N <- nrow(grid_df)
    
    # Origination (M-run above log(ns)/10)
    above_orig <- grid_df$ADF_delta > grid_df$boundary_orig
    orig_idx <- NA_integer_
    if (N >= M) {
      for (j in 1:(N - M + 1)) {
        if (all(above_orig[j:(j + M - 1)])) {
          orig_idx <- j
          break
        }
      }
    }
    
    if (is.na(orig_idx)) {
      return(list(
        grid_df = grid_df,
        unified_parts = list(
          orig_curve_df = data.frame(),
          collapse_curve_df = data.frame()
        ),
        episode = NULL,
        X = X,
        dates = dates
      ))
    }
    
    # Dates around bridge
    bridge_start_date <- grid_df$Date[orig_idx]
    bridge_end_date <- bridge_start_date + bridge_days
    idx_bridge_end <- which(grid_df$Date >= bridge_end_date)[1]
    if (is.na(idx_bridge_end)) idx_bridge_end <- N
    
    # Collapse within bridge: ADF < log(ns)/10 between orig+1 and bridge_end-1
    collapse_idx <- NA_integer_
    collapse_date <- as.Date(NA)
    collapse_regime <- NA_character_
    collapse_label <- NA_character_
    
    collapse_month_label <- function(d) {
      dday <- as.integer(format(d, "%d"))
      if (dday < 15) {
        format(d, "%b-%Y")
      } else {
        format(seq(d, by = "month", length.out = 2)[2], "%b-%Y")
      }
    }
    
    j_start <- orig_idx + 1
    j_end <- idx_bridge_end - 1
    if (j_start <= j_end) {
      for (jj in j_start:j_end) {
        if (grid_df$ADF_delta[jj] < grid_df$boundary_orig[jj]) {
          collapse_idx <- jj
          collapse_date <- grid_df$Date[jj]
          collapse_regime <- "bridge"
          collapse_label <- collapse_month_label(collapse_date)
          break
        }
      }
    }
    
    # If no bridge collapse, post-bridge collapse with log(ns)/2 R-run
    if (is.na(collapse_idx) && idx_bridge_end <= N) {
      last_possible_start <- N - R + 1
      if (last_possible_start >= idx_bridge_end) {
        for (j in idx_bridge_end:last_possible_start) {
          run_idx <- j:(j + R - 1)
          if (all(grid_df$ADF_delta[run_idx] < grid_df$boundary_collapse[run_idx])) {
            collapse_idx <- j
            collapse_date <- grid_df$Date[j]
            collapse_regime <- "post_bridge"
            collapse_label <- collapse_month_label(collapse_date)
            break
          }
        }
      }
    }
    
    # Bubble length
    orig_date <- grid_df$Date[orig_idx]
    if (is.na(collapse_idx)) {
      bubble_length_days <- as.integer(grid_df$Date[N] - orig_date)
      bubble_length_points <- N - orig_idx
    } else {
      bubble_length_days <- as.integer(collapse_date - orig_date)
      bubble_length_points <- collapse_idx - orig_idx
    }
    
    # Parts for plotting thresholds
    orig_curve_df <- data.frame(
      Date = grid_df$Date[1:orig_idx],
      y = grid_df$boundary_orig[1:orig_idx]
    )
    
    # Brown threshold line starts at origination + bridge_days + 30 days
    collapse_curve_start_date <- bridge_start_date + bridge_days + R
    idx_collapse_curve_start <- which(grid_df$Date >= collapse_curve_start_date)[1]
    
    collapse_curve_df <- if (!is.na(idx_collapse_curve_start) && idx_collapse_curve_start <= N) {
      data.frame(
        Date = grid_df$Date[idx_collapse_curve_start:N],
        y = grid_df$boundary_collapse[idx_collapse_curve_start:N]
      )
    } else {
      data.frame(Date = as.Date(character()), y = numeric())
    }
    
    ep <- data.frame(
      orig_idx = orig_idx,
      orig_date = orig_date,
      orig_label = format(orig_date, "%b-%Y"),
      collapse_idx = if (!is.na(collapse_idx)) collapse_idx else NA_integer_,
      collapse_date = if (!is.na(collapse_idx)) collapse_date else as.Date(NA),
      collapse_label = if (!is.na(collapse_idx)) collapse_label else NA_character_,
      collapse_regime = if (!is.na(collapse_idx)) collapse_regime else NA_character_,
      bubble_length_days = bubble_length_days,
      bubble_length_points = bubble_length_points,
      stringsAsFactors = FALSE
    )
    
    list(
      grid_df = grid_df,
      unified_parts = list(
        orig_curve_df = orig_curve_df,
        collapse_curve_df = collapse_curve_df,
        orig_idx = orig_idx,
        idx_bridge_end = idx_bridge_end,
        idx_collapse_curve_start = idx_collapse_curve_start
      ),
      episode = ep,
      X = X,
      dates = dates
    )
  }
  
  # -----------------------------
  # 3) Plot best window
  # -----------------------------
  plot_best_window <- function(window_xts, analysis_out, title_text, company_name) {
    grid_df <- analysis_out$grid_df
    parts <- analysis_out$unified_parts
    ep <- analysis_out$episode
    
    price_df <- data.frame(
      Date = as.Date(zoo::index(window_xts)),
      Price = as.numeric(window_xts)
    )
    
    left_min_raw <- min(c(
      grid_df$ADF_delta,
      if (nrow(parts$orig_curve_df) > 0) parts$orig_curve_df$y else Inf,
      if (nrow(parts$collapse_curve_df) > 0) parts$collapse_curve_df$y else Inf
    ), na.rm = TRUE)
    
    left_max_raw <- max(c(
      grid_df$ADF_delta,
      if (nrow(parts$orig_curve_df) > 0) parts$orig_curve_df$y else -Inf,
      if (nrow(parts$collapse_curve_df) > 0) parts$collapse_curve_df$y else -Inf
    ), na.rm = TRUE)
    
    if (!is.finite(left_min_raw) || !is.finite(left_max_raw) || left_max_raw == left_min_raw) {
      left_min_raw <- min(grid_df$ADF_delta, na.rm = TRUE)
      left_max_raw <- max(grid_df$ADF_delta, na.rm = TRUE) + 1
    }
    
    left_range <- left_max_raw - left_min_raw
    left_min_plot <- left_min_raw - 0.10 * left_range
    left_max_plot <- left_max_raw + 0.18 * left_range
    
    right_min <- min(price_df$Price, na.rm = TRUE)
    right_max <- max(price_df$Price, na.rm = TRUE)
    if (right_max == right_min) right_max <- right_min + 1
    
    price_band_low  <- left_min_plot + 0.02 * (left_max_plot - left_min_plot)
    price_band_high <- left_min_plot + 0.55 * (left_max_plot - left_min_plot)
    
    price_df$Price_scaled <- price_band_low +
      (price_df$Price - right_min) * (price_band_high - price_band_low) / (right_max - right_min)
    
    price_inverse <- function(z) {
      right_min + (z - price_band_low) * (right_max - right_min) / (price_band_high - price_band_low)
    }
    
    cols <- list(
      adf = "#1f77b4",
      b_orig = "#ff7f0e",
      b_collapse = "#8c564b",
      price = "#d62728",
      orig = "#2ca02c",
      coll = "#9467bd"
    )
    
    date_span_days <- as.integer(max(price_df$Date) - min(price_df$Date))
    label_offset_left  <- max(8, round(0.03 * date_span_days))
    label_offset_right <- max(8, round(0.03 * date_span_days))
    
    p <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = grid_df,
        ggplot2::aes(x = Date, y = ADF_delta),
        color = cols$adf, linewidth = 1
      ) +
      { if (nrow(parts$orig_curve_df) > 0)
        ggplot2::geom_line(
          data = parts$orig_curve_df,
          ggplot2::aes(x = Date, y = y),
          color = cols$b_orig,
          linetype = "dotted",
          linewidth = 0.9
        )
      } +
      { if (nrow(parts$collapse_curve_df) > 0)
        ggplot2::geom_line(
          data = parts$collapse_curve_df,
          ggplot2::aes(x = Date, y = y),
          color = cols$b_collapse,
          linetype = "dotted",
          linewidth = 0.9
        )
      } +
      ggplot2::geom_line(
        data = price_df,
        ggplot2::aes(x = Date, y = Price_scaled),
        color = cols$price, linewidth = 0.95, alpha = 0.95
      ) +
      ggplot2::scale_x_date(date_breaks = "3 months", date_labels = "%Y %b") +
      ggplot2::scale_y_continuous(
        name = NULL,
        limits = c(left_min_plot, left_max_plot),
        labels = NULL,
        sec.axis = ggplot2::sec_axis(~ price_inverse(.), name = NULL, labels = NULL)
      ) +
      ggplot2::labs(title = title_text, x = NULL, y = NULL) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "none",
        axis.title.x = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_blank(),
        axis.title.y.right = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
        axis.text.y.right = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.ticks.y.right = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        plot.margin = ggplot2::margin(10, 45, 10, 45)
      )
    
    # Origination marker and label
    if (!is.null(ep) && nrow(ep) == 1 && !is.na(ep$orig_date[1])) {
      orig_date_k <- as.Date(ep$orig_date[1])
      p <- p +
        ggplot2::geom_vline(
          xintercept = as.numeric(orig_date_k),
          linetype = "solid",
          linewidth = 1,
          color = cols$orig
        ) +
        ggplot2::annotate(
          "label",
          x = orig_date_k - label_offset_left,
          y = left_max_plot - 0.06 * (left_max_plot - left_min_plot),
          label = paste0("origination: ", ep$orig_label[1]),
          hjust = 1, vjust = 1, size = 3.2,
          color = cols$orig, fill = "white", label.size = 0.2
        )
    }
    
    # Collapse marker and label
    if (!is.null(ep) && nrow(ep) == 1 && !is.na(ep$collapse_date[1])) {
      collapse_date_k <- as.Date(ep$collapse_date[1])
      p <- p +
        ggplot2::geom_vline(
          xintercept = as.numeric(collapse_date_k),
          linetype = "solid",
          linewidth = 1,
          color = cols$coll
        ) +
        ggplot2::annotate(
          "label",
          x = collapse_date_k + label_offset_right,
          y = left_max_plot - 0.12 * (left_max_plot - left_min_plot),
          label = paste0("collapse: ", ep$collapse_label[1]),
          hjust = 0, vjust = 1, size = 3.2,
          color = cols$coll, fill = "white", label.size = 0.2
        )
    }
    
    print(p)
    invisible(p)
  }
  
  # -----------------------------
  # 4) Scan all rolling windows, pick longest bubble
  # -----------------------------
  all_windows <- list()
  rec_id <- 1
  
  for (h in horizons_years) {
    candidate_starts <- seq(from = search_start, to = search_end, by = paste(step_months, "months"))
    candidate_starts <- candidate_starts[candidate_starts <= (search_end - 365 * h)]
    if (length(candidate_starts) == 0) next
    
    for (s0 in candidate_starts) {
      e0 <- s0 + 365 * h
      window_xts <- price_xts_full[paste0(s0, "/", e0)]
      window_xts <- na.omit(window_xts)
      if (NROW(window_xts) < 30) next
      
      out <- analyze_window(window_xts)
      if (is.null(out) || is.null(out$episode)) next
      ep <- out$episode
      if (nrow(ep) != 1 || is.na(ep$orig_idx[1])) next
      
      all_windows[[rec_id]] <- data.frame(
        horizon_years = h,
        window_start = as.Date(min(zoo::index(window_xts))),
        window_end = as.Date(max(zoo::index(window_xts))),
        origination = as.Date(ep$orig_date[1]),
        collapse = if (!is.na(ep$collapse_date[1])) as.Date(ep$collapse_date[1]) else as.Date(NA),
        orig_label = ep$orig_label[1],
        collapse_label = if (!is.na(ep$collapse_label[1])) ep$collapse_label[1] else NA_character_,
        bubble_length_days = ep$bubble_length_days[1],
        bubble_length_points = ep$bubble_length_points[1],
        collapse_regime = if (!is.na(ep$collapse_regime[1])) ep$collapse_regime[1] else NA_character_,
        stringsAsFactors = FALSE
      )
      attr(all_windows[[rec_id]], "analysis") <- out
      attr(all_windows[[rec_id]], "window_xts") <- window_xts
      rec_id <- rec_id + 1
    }
  }
  
  if (length(all_windows) == 0) {
    cat("\nNo windows with any origination/collapse found under log(ns)/10 and log(ns)/2 rules.\n")
    return(invisible(NULL))
  }
  
  all_windows_df <- do.call(rbind, all_windows)
  ord <- order(-all_windows_df$bubble_length_days, all_windows_df$origination)
  best_idx <- ord[1]
  best_row <- all_windows_df[best_idx, ]
  best_analysis <- attr(all_windows[[best_idx]], "analysis")
  best_window_xts <- attr(all_windows[[best_idx]], "window_xts")
  
  cat("\nLongest bubble window for", company_name, "\n")
  print(best_row)
  
  if (isTRUE(plot_best)) {
    title_text <- paste0(
      company_name, " - Longest bubble window (",
      best_row$window_start, " to ", best_row$window_end, ")"
    )
    plot_best_window(
      window_xts = best_window_xts,
      analysis_out = best_analysis,
      title_text = title_text,
      company_name = company_name
    )
  }
  
  invisible(list(
    best_window = best_row,
    all_windows = all_windows_df,
    best_analysis = best_analysis
  ))
}