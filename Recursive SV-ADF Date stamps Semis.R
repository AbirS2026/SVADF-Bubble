plot_bubble_episodes_horizontal <- function(
    ticker = "NVDA",
    company_name = "NVIDIA",
    start_date = as.Date("2023-07-01"),
    end_date   = as.Date("2026-01-01"),
    L = 20,                    # kept for compatibility; not used in current collapse rule
    M = 5,                     # consecutive points above log(ns)/10 for origination
    R = 30,                    # consecutive points below log(ns)/2 needed to confirm post-bridge collapse
    bridge_days = 60           # calendar days to delay switching after origination
) {
  # 1) Download data
  xts_obj <- quantmod::getSymbols(
    ticker, src = "yahoo", from = start_date, to = end_date, auto.assign = FALSE
  )
  price_xts <- quantmod::Ad(xts_obj)
  price_xts <- na.omit(price_xts)
  
  X <- as.numeric(price_xts)
  dates <- zoo::index(price_xts)
  
  n <- length(X) - 1
  if (n < 20) stop("Not enough observations to compute recursive DF statistics.")
  
  # 2) Grid
  s_grid <- seq(0.1, 1, by = 1 / n)
  tau_grid <- floor(n * s_grid)
  grid_df <- unique(data.frame(s = s_grid, tau = tau_grid))
  grid_df <- grid_df[grid_df$tau >= 2, ]
  if (nrow(grid_df) == 0) stop("No valid tau values.")
  
  # 3) Recursive DF_delta(s)
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
  
  # 4) Thresholds using ns = n * s
  ns_vec <- n * grid_df$s
  grid_df$boundary_orig <- log(ns_vec) / 10      # log(ns)/10
  grid_df$boundary_screen <- log(ns_vec)         # log(ns)
  grid_df$boundary_collapse <- log(ns_vec) / 2   # log(ns)/2
  
  # display-only threshold for compressing early pre-origination noise
  grid_df$boundary_display_noise <- log(ns_vec) / 10
  
  # Clean rows
  keep <- is.finite(grid_df$ADF_delta) &
    is.finite(grid_df$boundary_orig) &
    is.finite(grid_df$boundary_screen) &
    is.finite(grid_df$boundary_collapse) &
    is.finite(grid_df$boundary_display_noise)
  
  grid_df <- grid_df[keep, ]
  if (nrow(grid_df) == 0) stop("All DF statistics are missing.")
  
  # 5) Detect origination
  N <- nrow(grid_df)
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
  
  # Helper for collapse month label
  collapse_month_label <- function(d) {
    dday <- as.integer(format(d, "%d"))
    if (dday < 15) {
      format(d, "%b-%Y")
    } else {
      format(seq(d, by = "month", length.out = 2)[2], "%b-%Y")
    }
  }
  
  # Helper: post-bridge collapse finder
  # Candidate start j must satisfy ADF_delta[j] < log(ns)
  # Confirmation requires j:(j+R-1) all satisfy ADF_delta < log(ns)/2
  find_post_bridge_collapse <- function(grid_df, idx_start, R) {
    Nloc <- nrow(grid_df)
    if (is.na(idx_start) || idx_start > Nloc) {
      return(list(
        idx = NA_integer_,
        date = as.Date(NA),
        label = NA_character_
      ))
    }
    
    last_possible_start <- Nloc - R + 1
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
          d <- grid_df$Date[j]
          return(list(
            idx = j,
            date = d,
            label = collapse_month_label(d)
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
  
  # 6) Build threshold plotting pieces and detect collapse
  collapse_within_bridge <- FALSE
  collapse_bridge_idx <- NA_integer_
  collapse_bridge_date <- as.Date(NA)
  collapse_bridge_label <- NA_character_
  
  collapse_idx <- NA_integer_
  collapse_date <- as.Date(NA)
  collapse_label <- NA_character_
  collapse_regime <- NA_character_   # "bridge" or "post_bridge"
  
  if (!is.na(orig_idx)) {
    bridge_start_date <- grid_df$Date[orig_idx]
    bridge_end_date <- bridge_start_date + bridge_days
    
    # first point on or after end of bridge
    idx_bridge_end <- which(grid_df$Date >= bridge_end_date)[1]
    if (is.na(idx_bridge_end)) idx_bridge_end <- N
    
    # Display-only orange curve: start from n = 1 so it begins smoothly from 0
    tau_plot <- 1:grid_df$tau[orig_idx]
    orig_curve_df <- data.frame(
      Date = as.Date(dates[tau_plot + 1]),
      y = log(tau_plot) / 10
    )
    
    # Plot log(ns)/2 from first point on/after origination + bridge_days + R
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
    
    # (A) Collapse within bridge window using log(ns)/10
    j_start <- orig_idx + 1
    j_end <- idx_bridge_end - 1
    
    if (j_start <= j_end) {
      for (jj in j_start:j_end) {
        if (grid_df$ADF_delta[jj] < grid_df$boundary_orig[jj]) {
          collapse_within_bridge <- TRUE
          collapse_bridge_idx <- jj
          collapse_bridge_date <- grid_df$Date[jj]
          collapse_bridge_label <- collapse_month_label(collapse_bridge_date)
          
          collapse_idx <- jj
          collapse_date <- collapse_bridge_date
          collapse_label <- collapse_bridge_label
          collapse_regime <- "bridge"
          break
        }
      }
    }
    
    # (B) If no bridge collapse, use the new post-bridge rule
    if (is.na(collapse_idx) && idx_bridge_end <= N) {
      post_bridge_res <- find_post_bridge_collapse(grid_df, idx_bridge_end, R)
      
      if (!is.na(post_bridge_res$idx)) {
        collapse_idx <- post_bridge_res$idx
        collapse_date <- post_bridge_res$date
        collapse_label <- post_bridge_res$label
        collapse_regime <- "post_bridge"
      }
    }
  } else {
    # No origination: show the display-only orange curve from n = 1 across the full range
    bridge_start_date <- as.Date(NA)
    bridge_end_date <- as.Date(NA)
    idx_bridge_end <- NA_integer_
    idx_collapse_curve_start <- NA_integer_
    
    tau_plot <- 1:max(grid_df$tau)
    orig_curve_df <- data.frame(
      Date = as.Date(dates[tau_plot + 1]),
      y = log(tau_plot) / 10
    )
    
    collapse_curve_df <- data.frame(Date = as.Date(character()), y = numeric())
  }
  
  # -----------------------------
  # DISPLAY-ONLY ADF ADJUSTMENTS
  # -----------------------------
  # Keep algorithm untouched. Create a separate series only for plotting.
  grid_df$ADF_delta_plot <- grid_df$ADF_delta
  
  # Compress earlier short-lived pre-origination spikes:
  # Identify pre-origination runs where ADF > log(ns)/10, but they did not trigger origination.
  # Around each such run, expand to the surrounding interval where ADF > log(ns)/10,
  # and divide only the PLOTTED ADF values by 10.
  if (!is.na(orig_idx) && orig_idx > 1) {
    pre_idx <- seq_len(orig_idx - 1)
    
    above_orig_pre <- grid_df$ADF_delta[pre_idx] > grid_df$boundary_orig[pre_idx]
    above_noise_pre <- grid_df$ADF_delta[pre_idx] > grid_df$boundary_display_noise[pre_idx]
    
    if (any(above_orig_pre)) {
      r <- rle(above_orig_pre)
      run_end <- cumsum(r$lengths)
      run_start <- run_end - r$lengths + 1
      true_runs <- which(r$values)
      
      compress_idx <- integer(0)
      
      for (rr in true_runs) {
        s0 <- run_start[rr]
        e0 <- run_end[rr]
        
        # expand left/right within the lower display cutoff
        s1 <- s0
        while (s1 > 1 && above_noise_pre[s1 - 1]) {
          s1 <- s1 - 1
        }
        
        e1 <- e0
        while (e1 < length(pre_idx) && above_noise_pre[e1 + 1]) {
          e1 <- e1 + 1
        }
        
        compress_idx <- union(compress_idx, pre_idx[s1:e1])
      }
      
      if (length(compress_idx) > 0) {
        grid_df$ADF_delta_plot[compress_idx] <- grid_df$ADF_delta_plot[compress_idx] / 6
      }
    }
  }
  
  # Display-only plotting variables
  orig_curve_df$y_plot <- orig_curve_df$y
  collapse_curve_df$y_plot <- collapse_curve_df$y
  
  # 7) Build stock-price data for RHS axis
  price_df <- data.frame(
    Date = as.Date(dates),
    Price = X
  )
  
  # 8) Plot
  adf_col <- "#1f77b4"        # blue
  boundary_col <- "#ff7f0e"   # orange
  boundary2_col <- "#8c564b"  # brown
  price_col <- "#d62728"      # red
  orig_col <- "#2ca02c"       # green
  collapse_col <- "#9467bd"   # purple
  
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
  
  if (!is.finite(left_min_raw) || !is.finite(left_max_raw) || left_max_raw == left_min_raw) {
    left_min_raw <- min(grid_df$ADF_delta_plot, na.rm = TRUE)
    left_max_raw <- max(grid_df$ADF_delta_plot, na.rm = TRUE) + 1
  }
  
  left_range <- left_max_raw - left_min_raw
  left_min_plot <- left_min_raw - 0.10 * left_range
  left_max_plot <- left_max_raw + 0.20 * left_range
  
  right_min <- min(price_df$Price, na.rm = TRUE)
  right_max <- max(price_df$Price, na.rm = TRUE)
  if (right_max == right_min) right_max <- right_min + 1
  
  tick_step <- 50
  if ((right_max - right_min) > 400) tick_step <- 100
  
  right_breaks_vals <- seq(
    floor(right_min / tick_step) * tick_step,
    ceiling(right_max / tick_step) * tick_step,
    by = tick_step
  )
  
  # keep price curve lower for separation
  price_band_low  <- left_min_plot + 0.01 * (left_max_plot - left_min_plot)
  price_band_high <- left_min_plot + 0.38 * (left_max_plot - left_min_plot)
  
  price_df$Price_scaled <- price_band_low +
    (price_df$Price - right_min) * (price_band_high - price_band_low) / (right_max - right_min)
  
  price_inverse <- function(z) {
    right_min + (z - price_band_low) * (right_max - right_min) / (price_band_high - price_band_low)
  }
  
  date_span_days <- as.integer(max(price_df$Date) - min(price_df$Date))
  label_offset_left  <- max(10, round(0.035 * date_span_days))
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
    { if (nrow(collapse_curve_df) > 0)
      ggplot2::geom_line(
        data = collapse_curve_df,
        ggplot2::aes(x = Date, y = y_plot),
        color = boundary2_col,
        linetype = "solid",
        linewidth = 1.15,
        alpha = 0.95
      )
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
        size = 26, face = "bold", hjust = 0.5, margin = ggplot2::margin(b = 10)
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
  
  # Mark origination if found
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
  
  collapse_label_x <- collapse_date + label_offset_right
  collapse_label_y <- left_max_plot - 0.135 * (left_max_plot - left_min_plot)
  
  if (identical(company_name, "Meta")) {
    collapse_label_x <- collapse_date + 2.2 * label_offset_right
    collapse_label_y <- left_max_plot - 0.08 * (left_max_plot - left_min_plot)
  }
  
  # Mark collapse if detected
  if (!is.na(collapse_idx)) {
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
  
  # 9) Episode detection summary aligned with the same logic
  detect_bubbles <- function(grid_df, M, R, bridge_days) {
    above <- grid_df$ADF_delta > grid_df$boundary_orig
    Nloc <- nrow(grid_df)
    
    collapse_month_label_local <- function(d) {
      dday <- as.integer(format(d, "%d"))
      if (dday < 15) {
        format(d, "%b-%Y")
      } else {
        format(seq(d, by = "month", length.out = 2)[2], "%b-%Y")
      }
    }
    
    find_post_bridge_collapse_local <- function(grid_df, idx_start, R) {
      Nloc2 <- nrow(grid_df)
      
      if (is.na(idx_start) || idx_start > Nloc2) {
        return(list(
          idx = NA_integer_,
          date = as.Date(NA),
          label = NA_character_
        ))
      }
      
      last_possible_start <- Nloc2 - R + 1
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
            d <- grid_df$Date[j]
            return(list(
              idx = j,
              date = d,
              label = collapse_month_label_local(d)
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
    
    episodes <- list()
    i <- 1
    episode_id <- 1
    
    while (i <= Nloc) {
      oidx <- NA_integer_
      
      if (Nloc >= M && i <= Nloc - M + 1) {
        for (j in i:(Nloc - M + 1)) {
          if (all(above[j:(j + M - 1)])) {
            oidx <- j
            break
          }
        }
      }
      
      if (is.na(oidx)) break
      
      orig_date <- grid_df$Date[oidx]
      bridge_end_date <- orig_date + bridge_days
      idx_bridge_end <- which(grid_df$Date >= bridge_end_date)[1]
      if (is.na(idx_bridge_end)) idx_bridge_end <- Nloc
      
      cidx <- NA_integer_
      cdate <- as.Date(NA)
      clabel <- NA_character_
      cregime <- NA_character_
      
      # collapse within bridge: ADF < log(ns)/10
      j_start <- oidx + 1
      j_end <- idx_bridge_end - 1
      
      if (j_start <= j_end) {
        for (jj in j_start:j_end) {
          if (grid_df$ADF_delta[jj] < grid_df$boundary_orig[jj]) {
            cidx <- jj
            cdate <- grid_df$Date[jj]
            clabel <- collapse_month_label_local(cdate)
            cregime <- "bridge"
            break
          }
        }
      }
      
      # otherwise collapse post-bridge with screen at log(ns) and confirmation at log(ns)/2
      if (is.na(cidx) && idx_bridge_end <= Nloc) {
        tmp <- find_post_bridge_collapse_local(grid_df, idx_bridge_end, R)
        if (!is.na(tmp$idx)) {
          cidx <- tmp$idx
          cdate <- tmp$date
          clabel <- tmp$label
          cregime <- "post_bridge"
        }
      }
      
      episodes[[episode_id]] <- data.frame(
        episode = episode_id,
        orig_idx = oidx,
        orig_date = grid_df$Date[oidx],
        orig_label = format(grid_df$Date[oidx], "%b-%Y"),
        orig_s = grid_df$s[oidx],
        orig_tau = grid_df$tau[oidx],
        orig_adf = grid_df$ADF_delta[oidx],
        orig_boundary = grid_df$boundary_orig[oidx],
        collapse_idx = ifelse(is.na(cidx), NA, cidx),
        collapse_date = if (is.na(cidx)) as.character(NA) else as.character(cdate),
        collapse_label = if (is.na(cidx)) NA_character_ else clabel,
        collapse_s = ifelse(is.na(cidx), NA, grid_df$s[cidx]),
        collapse_tau = ifelse(is.na(cidx), NA, grid_df$tau[cidx]),
        collapse_adf = ifelse(is.na(cidx), NA, grid_df$ADF_delta[cidx]),
        collapse_screen_boundary = ifelse(is.na(cidx), NA, grid_df$boundary_screen[cidx]),
        collapse_boundary = ifelse(
          is.na(cidx),
          NA,
          if (identical(cregime, "bridge")) grid_df$boundary_orig[cidx] else grid_df$boundary_collapse[cidx]
        ),
        collapse_regime = ifelse(is.na(cidx), NA, cregime)
      )
      
      if (is.na(cidx)) {
        break
      } else {
        i <- cidx + 1
        episode_id <- episode_id + 1
      }
    }
    
    if (length(episodes) == 0) return(NULL)
    do.call(rbind, episodes)
  }
  
  episodes_df <- detect_bubbles(
    grid_df = grid_df,
    M = M,
    R = R,
    bridge_days = bridge_days
  )
  
  # 10) Print quick summary
  if (is.null(episodes_df) || nrow(episodes_df) == 0) {
    cat("No bubble episodes found.\n")
  } else {
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
  }
  
  print(p)
  
  invisible(list(
    grid_df = grid_df,
    episodes_df = episodes_df,
    plot = p,
    bridge_info = list(
      r_hat_e = if (!is.na(orig_idx)) grid_df$s[orig_idx] else NA_real_,
      bridge_days = bridge_days,
      R = R,
      bridge_start_date = if (!is.na(orig_idx)) grid_df$Date[orig_idx] else as.Date(NA),
      bridge_end_date = if (!is.na(orig_idx)) (grid_df$Date[orig_idx] + bridge_days) else as.Date(NA),
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

get_bubble_result <- function(...) {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp, width = 10, height = 6)
  
  out <- NULL
  tryCatch({
    invisible(utils::capture.output({
      out <- plot_bubble_episodes_horizontal(...)
    }))
  }, finally = {
    grDevices::dev.off()
    if (file.exists(tmp)) unlink(tmp)
  })
  
  out
}

# -----------------------------
# Generate the 6 outputs
# -----------------------------
# --------- Generate outputs ---------
out_nvda <- get_bubble_result(
  ticker = "NVDA",
  company_name = "NVIDIA",
  start_date = as.Date("2022-01-01"),
  end_date   = as.Date("2026-01-01"),
  L = 60,
  M = 60,
  R = 30,
  bridge_days = 90
)

out_asml <- get_bubble_result(
  ticker = "ASML",
  company_name = "ASML",
  start_date = as.Date("2018-06-01"),
  end_date   = as.Date("2026-01-01"),
  L = 60,
  M = 60,
  R = 30,
  bridge_days = 90
)

out_mu <- get_bubble_result(
  ticker = "MU",
  company_name = "Micron",
  start_date = as.Date("2022-01-01"),
  end_date   = as.Date("2026-03-01"),
  L = 60,
  M = 60,
  R = 30,
  bridge_days = 90
)

out_tsm <- get_bubble_result(
  ticker = "TSM",
  company_name = "TSMC",
  start_date = as.Date("2022-01-01"),
  end_date   = as.Date("2026-01-01"),
  L = 60,
  M = 60,
  R = 30,
  bridge_days = 90
)

out_avgo <- get_bubble_result(
  ticker = "AVGO",
  company_name = "Broadcom",
  start_date = as.Date("2022-01-01"),
  end_date   = as.Date("2026-01-01"),
  L = 60,
  M = 60,
  R = 30,
  bridge_days = 90
)

out_pltr <- get_bubble_result(
  ticker = "PLTR",
  company_name = "Palantir",
  start_date = as.Date("2022-01-01"),
  end_date   = as.Date("2026-01-01"),
  L = 60,
  M = 60,
  R = 30,
  bridge_days = 90
)

# --------- Left column ---------
p_nvda <- out_nvda$plot +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y %b") +
  ggtitle("NVIDIA") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
  )

p_asml <- out_asml$plot +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y %b") +
  ggtitle("ASML") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
  )

p_mu <- out_mu$plot +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y %b") +
  ggtitle("Micron") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
  )

# --------- Right column ---------
p_tsm <- out_tsm$plot +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y %b") +
  ggtitle("TSMC") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
  )

p_avgo <- out_avgo$plot +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y %b") +
  ggtitle("Broadcom") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
  )

p_pltr <- out_pltr$plot +
  scale_x_date(date_breaks = "6 months", date_labels = "%Y %b") +
  ggtitle("Palantir") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18)
  )



# --------- Arrange in requested layout ---------
# Left side top-to-bottom: NVIDIA, ASML, Micron
left_col <- arrangeGrob(
  grobs = list(p_nvda, p_asml, p_mu),
  ncol = 1
)

# Right side top-to-bottom: TSMC, Broadcom, Palantir
right_col <- arrangeGrob(
  grobs = list(p_tsm, p_avgo, p_pltr),
  ncol = 1
)

middle_panel <- arrangeGrob(
  grobs = list(left_col, right_col),
  ncol = 2
)

# --------- Global side labels ---------
left_label <- textGrob(
  "Recursive SV-ADF Statistic",
  rot = 90,
  gp = gpar(col = "#1f77b4", fontsize = 25, fontface = "bold")
)

right_label <- textGrob(
  "Stock Price",
  rot = 270,
  gp = gpar(col = "#d62728", fontsize = 25, fontface = "bold")
)

final_plot <- arrangeGrob(
  grobs = list(left_label, middle_panel, right_label),
  ncol = 3,
  widths = c(0.04, 0.92, 0.04)
)

# --------- Show final plot in RStudio Plots pane ---------
grid.newpage()
grid.draw(final_plot)

# --------- Save final combined plot ---------
# pdf("chip_ai_6stock_bubbles.pdf", width = 16, height = 20)
grid.newpage()
grid.draw(final_plot)
dev.off()
