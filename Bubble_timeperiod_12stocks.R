################################################################################
# Bubble_timeperiod_12stocks.R
#
# Replication code for Figure 11 and Table 2 of:
# "Is There an AI Bubble? Robust Date-Stamping for Periods of Exuberance"
#
# Figure 11:
#   Estimated AR coefficient delta_n and explosive rate gamma_n, with confidence
#   intervals, for the twelve AI-exposed technology and semiconductor firms.
#
# Table 2:
#   Bubble/no-bubble classifications for the twelve firms across the two main
#   empirical periods:
#     - 2018--2020
#     - 2022--2026
#
# Data source:
# Yahoo Finance daily adjusted closing prices.
#
# Notes:
# - Package installation is intentionally omitted; see README.md.
# - Set SAVE_FIGURES <- TRUE to save Figure 11.
# - Set SAVE_TABLES <- TRUE to save the Table 2 decision outputs.
################################################################################

library(quantmod)
library(ggplot2)
library(dplyr)

SAVE_FIGURES <- TRUE
SAVE_TABLES <- TRUE
SAVE_PNG_PREVIEW <- TRUE

OUTPUT_DIR <- "figures"
TABLE_DIR <- "tables"

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

if (!dir.exists(TABLE_DIR)) {
  dir.create(TABLE_DIR, recursive = TRUE)
}

FIGURE11_PDF <- file.path(OUTPUT_DIR, "delta_gamma_dual_axis_plot.pdf")
FIGURE11_PNG <- file.path(OUTPUT_DIR, "delta_gamma_dual_axis_plot.png")

TABLE2_CSV <- file.path(TABLE_DIR, "table2_bubble_decisions_two_periods.csv")
TABLE2_LATEX_ROWS <- file.path(TABLE_DIR, "table2_latex_rows.txt")

ALL_TICKERS <- c(
  "AAPL", "MSFT", "NVDA", "GOOGL", "AMZN", "META",
  "TSLA", "TSM", "AVGO", "PLTR", "ASML", "MU"
)

COMPANY_LABELS <- c(
  AAPL  = "Apple",
  MSFT  = "Microsoft",
  NVDA  = "Nvidia",
  GOOGL = "Alphabet",
  AMZN  = "Amazon",
  META  = "Meta",
  TSLA  = "Tesla",
  TSM   = "TSMC",
  AVGO  = "Broadcom",
  PLTR  = "Palantir",
  ASML  = "ASML",
  MU    = "Micron"
)

fetch_adjusted_prices <- function(ticker, start_date, end_date) {
  xts_obj <- tryCatch(
    quantmod::getSymbols(
      Symbols = ticker,
      src = "yahoo",
      from = as.Date(start_date),
      to = as.Date(end_date),
      auto.assign = FALSE,
      warnings = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(xts_obj)) {
    return(NULL)
  }
  
  price_xts <- tryCatch(quantmod::Ad(xts_obj), error = function(e) NULL)
  
  if (is.null(price_xts)) {
    price_xts <- tryCatch(quantmod::Cl(xts_obj), error = function(e) NULL)
  }
  
  if (is.null(price_xts)) {
    return(NULL)
  }
  
  na.omit(price_xts)
}

resolve_ticker <- function(x) {
  name_to_ticker <- c(
    "apple" = "AAPL",
    "microsoft" = "MSFT",
    "nvidia" = "NVDA",
    "alphabet" = "GOOGL",
    "google" = "GOOGL",
    "amazon" = "AMZN",
    "meta" = "META",
    "tesla" = "TSLA",
    "tsmc" = "TSM",
    "taiwan semiconductor" = "TSM",
    "broadcom" = "AVGO",
    "palantir" = "PLTR",
    "asml" = "ASML",
    "micron" = "MU"
  )
  
  x_trim <- trimws(x)
  x_upper <- toupper(x_trim)
  x_lower <- tolower(x_trim)
  
  if (x_upper %in% names(COMPANY_LABELS)) {
    return(x_upper)
  }
  
  if (x_lower %in% names(name_to_ticker)) {
    return(unname(name_to_ticker[x_lower]))
  }
  
  stop(sprintf("Unknown stock name/ticker: '%s'", x))
}

decision_symbol <- function(x) {
  ifelse(x == "Bubble", "\\bubble", "\\nobubble")
}

compute_delta_gamma_stats <- function(price_xts, lambda = 0.05) {
  X <- as.numeric(na.omit(price_xts))
  
  if (length(X) < 3) {
    return(c(
      T = NA_real_,
      delta_hat = NA_real_,
      gamma_hat = NA_real_,
      delta_lower = NA_real_,
      delta_upper = NA_real_,
      gamma_lower = NA_real_,
      gamma_upper = NA_real_,
      regime_code = NA_real_
    ))
  }
  
  y <- X[-1]
  xlag <- X[-length(X)]
  Tn <- length(y)
  
  delta_hat <- sum(y * xlag) / sum(xlag^2)
  
  diff_from_one <- abs(delta_hat - 1)
  
  gamma_hat <- if (diff_from_one <= .Machine$double.eps) {
    NA_real_
  } else {
    -log(diff_from_one) / log(Tn)
  }
  
  delta_lower <- delta_upper <- NA_real_
  gamma_lower <- gamma_upper <- NA_real_
  regime_code <- NA_real_
  
  cv_normal <- qnorm(1 - lambda / 2)
  cv_cauchy <- qcauchy(1 - lambda / 2)
  
  if (!is.na(delta_hat) && !is.na(gamma_hat)) {
    
    if (delta_hat < 1) {
      regime_code <- 0
      
      delta_half_width <- cv_normal * sqrt(2) / (Tn^((1 + gamma_hat) / 2))
      delta_lower <- delta_hat - delta_half_width
      delta_upper <- delta_hat + delta_half_width
      
      gamma_half_width <- cv_normal * sqrt(2) /
        ((Tn^((1 - gamma_hat) / 2)) * log(Tn))
      gamma_lower <- gamma_hat - gamma_half_width
      gamma_upper <- gamma_hat + gamma_half_width
      
    } else if (delta_hat > 1) {
      regime_code <- 1
      
      delta_half_width <- cv_cauchy * 2 /
        ((Tn^gamma_hat) * (delta_hat^Tn))
      delta_lower <- delta_hat - delta_half_width
      delta_upper <- delta_hat + delta_half_width
      
      gamma_half_width <- cv_cauchy * 2 /
        (((1 + 1 / (Tn^gamma_hat))^Tn) * log(Tn))
      gamma_lower <- gamma_hat - gamma_half_width
      gamma_upper <- gamma_hat + gamma_half_width
    }
  }
  
  c(
    T = Tn,
    delta_hat = delta_hat,
    gamma_hat = gamma_hat,
    delta_lower = delta_lower,
    delta_upper = delta_upper,
    gamma_lower = gamma_lower,
    gamma_upper = gamma_upper,
    regime_code = regime_code
  )
}

scale_to_band <- function(x, xmin, xmax, low, high) {
  if (!is.finite(xmin) || !is.finite(xmax) || xmax == xmin) {
    return(rep((low + high) / 2, length(x)))
  }
  
  low + (x - xmin) * (high - low) / (xmax - xmin)
}

make_delta_gamma_dual_plot <- function(
    start_date = "2022-01-01",
    end_date = "2026-03-01",
    lambda = 0.05,
    save_figure = SAVE_FIGURES
) {
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  
  price_list <- lapply(ALL_TICKERS, function(ticker) {
    px <- fetch_adjusted_prices(ticker, start_date, end_date)
    
    if (is.null(px) || NROW(px) < 3) {
      warning(paste("Could not load sufficient data for", ticker))
      return(NULL)
    }
    
    px
  })
  
  names(price_list) <- ALL_TICKERS
  
  est_mat <- t(sapply(ALL_TICKERS, function(ticker) {
    px <- price_list[[ticker]]
    
    if (is.null(px)) {
      return(c(
        T = NA_real_,
        delta_hat = NA_real_,
        gamma_hat = NA_real_,
        delta_lower = NA_real_,
        delta_upper = NA_real_,
        gamma_lower = NA_real_,
        gamma_upper = NA_real_,
        regime_code = NA_real_
      ))
    }
    
    compute_delta_gamma_stats(px, lambda = lambda)
  }))
  
  results_df <- data.frame(
    Ticker = ALL_TICKERS,
    Company = unname(COMPANY_LABELS[ALL_TICKERS]),
    T = est_mat[, "T"],
    Delta_Hat = est_mat[, "delta_hat"],
    Gamma_Hat = est_mat[, "gamma_hat"],
    Delta_CI_Lower = est_mat[, "delta_lower"],
    Delta_CI_Upper = est_mat[, "delta_upper"],
    Gamma_CI_Lower = est_mat[, "gamma_lower"],
    Gamma_CI_Upper = est_mat[, "gamma_upper"],
    Regime = ifelse(
      est_mat[, "regime_code"] == 0,
      "MIR (<1)",
      ifelse(est_mat[, "regime_code"] == 1, "MER (>1)", "Boundary (=1)")
    ),
    row.names = NULL
  )
  
  results_df$Company <- factor(
    results_df$Company,
    levels = unname(COMPANY_LABELS[ALL_TICKERS])
  )
  
  delta_band_low <- 0.08
  delta_band_high <- 0.42
  gamma_band_low <- 0.58
  gamma_band_high <- 0.92
  
  dmin <- min(results_df$Delta_CI_Lower, results_df$Delta_Hat, na.rm = TRUE)
  dmax <- max(results_df$Delta_CI_Upper, results_df$Delta_Hat, na.rm = TRUE)
  
  gmin <- min(results_df$Gamma_CI_Lower, results_df$Gamma_Hat, na.rm = TRUE)
  gmax <- max(results_df$Gamma_CI_Upper, results_df$Gamma_Hat, na.rm = TRUE)
  
  results_df$Delta_Hat_Plot <- scale_to_band(
    results_df$Delta_Hat, dmin, dmax, delta_band_low, delta_band_high
  )
  results_df$Delta_CI_Lower_Plot <- scale_to_band(
    results_df$Delta_CI_Lower, dmin, dmax, delta_band_low, delta_band_high
  )
  results_df$Delta_CI_Upper_Plot <- scale_to_band(
    results_df$Delta_CI_Upper, dmin, dmax, delta_band_low, delta_band_high
  )
  
  results_df$Gamma_Hat_Plot <- scale_to_band(
    results_df$Gamma_Hat, gmin, gmax, gamma_band_low, gamma_band_high
  )
  results_df$Gamma_CI_Lower_Plot <- scale_to_band(
    results_df$Gamma_CI_Lower, gmin, gmax, gamma_band_low, gamma_band_high
  )
  results_df$Gamma_CI_Upper_Plot <- scale_to_band(
    results_df$Gamma_CI_Upper, gmin, gmax, gamma_band_low, gamma_band_high
  )
  
  delta_breaks_vals <- pretty(c(dmin, dmax), n = 5)
  gamma_breaks_vals <- pretty(c(gmin, gmax), n = 5)
  
  delta_breaks_pos <- scale_to_band(
    delta_breaks_vals, dmin, dmax, delta_band_low, delta_band_high
  )
  gamma_breaks_pos <- scale_to_band(
    gamma_breaks_vals, gmin, gmax, gamma_band_low, gamma_band_high
  )
  
  delta_one_pos <- scale_to_band(1, dmin, dmax, delta_band_low, delta_band_high)
  gamma_one_pos <- scale_to_band(1, gmin, gmax, gamma_band_low, gamma_band_high)
  
  delta_bar_color <- "#0047AB"
  delta_dot_color <- "#d62728"
  gamma_bar_color <- "black"
  gamma_dot_color <- "#2E8B57"
  
  p <- ggplot(results_df, aes(x = Company)) +
    geom_errorbar(
      aes(ymin = Delta_CI_Lower_Plot, ymax = Delta_CI_Upper_Plot),
      color = delta_bar_color,
      width = 0.18,
      linewidth = 0.8
    ) +
    geom_point(
      aes(y = Delta_Hat_Plot),
      color = delta_dot_color,
      size = 2.5
    ) +
    geom_hline(
      yintercept = delta_one_pos,
      linetype = "dashed",
      color = delta_bar_color
    ) +
    geom_errorbar(
      aes(ymin = Gamma_CI_Lower_Plot, ymax = Gamma_CI_Upper_Plot),
      color = gamma_bar_color,
      width = 0.18,
      linewidth = 0.8
    ) +
    geom_point(
      aes(y = Gamma_Hat_Plot),
      color = gamma_dot_color,
      size = 2.5
    ) +
    geom_hline(
      yintercept = gamma_one_pos,
      linetype = "dashed",
      color = gamma_bar_color
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      name = expression("Estimated AR coefficient and C.I. for " * delta[n]),
      breaks = delta_breaks_pos,
      labels = round(delta_breaks_vals, 3),
      sec.axis = sec_axis(
        ~ .,
        name = expression("Estimated explosive growth rate and C.I. for " * gamma[n]),
        breaks = gamma_breaks_pos,
        labels = round(gamma_breaks_vals, 3)
      )
    ) +
    labs(x = NULL, title = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y.left = element_text(size = 20, color = delta_bar_color),
      axis.title.y.right = element_text(size = 20, color = gamma_bar_color),
      axis.text.x = element_text(size = 18, angle = 45, hjust = 1, face = "bold"),
      axis.text.y.left = element_text(size = 20, color = delta_bar_color),
      axis.text.y.right = element_text(size = 20, color = gamma_bar_color),
      axis.ticks.y.left = element_line(color = delta_bar_color),
      axis.ticks.y.right = element_line(color = gamma_bar_color),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  
  print(p)
  
  if (isTRUE(save_figure)) {
    ggsave(FIGURE11_PDF, plot = p, width = 12, height = 10)
    
    if (isTRUE(SAVE_PNG_PREVIEW)) {
      ggsave(FIGURE11_PNG, plot = p, width = 12, height = 10, dpi = 300)
    }
  }
  
  invisible(list(
    estimates = results_df,
    plot = p
  ))
}

compute_df_delta <- function(X) {
  X <- as.numeric(X)
  X <- X[is.finite(X)]
  
  if (length(X) < 3) {
    return(list(
      n_eff = NA_real_,
      delta_hat = NA_real_,
      df_delta = NA_real_
    ))
  }
  
  Y <- X[-1]
  Xlag <- X[-length(X)]
  n_eff <- length(Y)
  
  denom <- sum(Xlag^2)
  
  if (!is.finite(denom) || denom <= 0) {
    return(list(
      n_eff = n_eff,
      delta_hat = NA_real_,
      df_delta = NA_real_
    ))
  }
  
  delta_hat <- sum(Y * Xlag) / denom
  df_delta <- n_eff * (delta_hat - 1)
  
  list(
    n_eff = n_eff,
    delta_hat = delta_hat,
    df_delta = df_delta
  )
}

simulate_null_unit_root_sv <- function(
    n_eff,
    d_par = 1,
    eta_par = 0.1,
    x0 = 1,
    sigma2_0 = 1
) {
  if (!is.finite(n_eff) || n_eff < 20) {
    stop("n_eff must be at least 20 for the log(log n) volatility specification.")
  }
  
  loglog_n <- log(log(n_eff))
  
  if (!is.finite(loglog_n) || loglog_n == 0) {
    stop("log(log(n_eff)) is not finite; increase sample size.")
  }
  
  phi_n <- 1 - d_par / loglog_n
  
  X <- numeric(n_eff + 1)
  X[1] <- x0
  
  log_sigma2_prev <- log(sigma2_0)
  
  for (t in seq_len(n_eff)) {
    eta_t <- rnorm(1, mean = 0, sd = eta_par)
    log_sigma2_t <- phi_n * log_sigma2_prev + eta_t
    sigma_t <- sqrt(exp(log_sigma2_t))
    u_t <- sigma_t * rnorm(1)
    
    X[t + 1] <- X[t] + u_t
    log_sigma2_prev <- log_sigma2_t
  }
  
  X
}

bubble_decisions_ar1 <- function(
    stock_names,
    start_date,
    end_date,
    B = 1000,
    d_par = 1,
    eta_par = 0.1,
    x0_boot = 1,
    sigma2_0_boot = 1,
    normalize_start = TRUE,
    seed = NULL,
    verbose = TRUE
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  tickers <- vapply(stock_names, resolve_ticker, character(1))
  
  boot_cache <- new.env(parent = emptyenv())
  
  bootstrap_cv_90 <- function(n_eff) {
    key <- paste(n_eff, B, d_par, eta_par, x0_boot, sigma2_0_boot, sep = "|")
    
    if (exists(key, envir = boot_cache, inherits = FALSE)) {
      return(get(key, envir = boot_cache, inherits = FALSE))
    }
    
    bootstrap_values <- numeric(B)
    
    for (b in seq_len(B)) {
      X_star <- simulate_null_unit_root_sv(
        n_eff = n_eff,
        d_par = d_par,
        eta_par = eta_par,
        x0 = x0_boot,
        sigma2_0 = sigma2_0_boot
      )
      
      bootstrap_values[b] <- compute_df_delta(X_star)$df_delta
    }
    
    bootstrap_values <- bootstrap_values[is.finite(bootstrap_values)]
    
    out <- if (length(bootstrap_values) == 0) {
      list(cv_90 = NA_real_, boot_vals = numeric(0))
    } else {
      list(
        cv_90 = unname(quantile(bootstrap_values, probs = 0.90, na.rm = TRUE, type = 7)),
        boot_vals = bootstrap_values
      )
    }
    
    assign(key, out, envir = boot_cache)
    out
  }
  
  out_list <- vector("list", length(tickers))
  
  for (i in seq_along(tickers)) {
    ticker <- tickers[i]
    
    if (isTRUE(verbose)) {
      message(sprintf("Processing %s (%d/%d) ...", ticker, i, length(tickers)))
    }
    
    price_xts <- fetch_adjusted_prices(ticker, start_date, end_date)
    
    if (is.null(price_xts)) {
      out_list[[i]] <- data.frame(
        Ticker = ticker,
        Company = unname(COMPANY_LABELS[ticker]),
        Start_Date = as.character(as.Date(start_date)),
        End_Date = as.character(as.Date(end_date)),
        N_Obs = NA_integer_,
        N_Eff = NA_integer_,
        Delta_Hat = NA_real_,
        DF_Delta = NA_real_,
        Threshold_LogN_10 = NA_real_,
        Threshold_LogLogN_100 = NA_real_,
        Bootstrap_CV_90 = NA_real_,
        Decision_LogN_10 = NA_character_,
        Decision_LogLogN_100 = NA_character_,
        Decision_Bootstrap_90 = NA_character_,
        Bootstrap_P_Value = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    
    X <- as.numeric(price_xts)
    X <- X[is.finite(X)]
    
    if (normalize_start && length(X) > 0) {
      X <- X / X[1]
    }
    
    obs <- compute_df_delta(X)
    n_eff <- obs$n_eff
    
    threshold_log_n_10 <- log(n_eff) / 10
    threshold_loglog_n_100 <- log(log(n_eff)) / 100
    
    boot_obj <- bootstrap_cv_90(n_eff)
    cv90 <- boot_obj$cv_90
    boot_vals <- boot_obj$boot_vals
    
    bootstrap_p_value <- if (length(boot_vals) > 0 && is.finite(obs$df_delta)) {
      mean(boot_vals >= obs$df_delta)
    } else {
      NA_real_
    }
    
    decision_log_n_10 <- if (
      is.finite(obs$df_delta) &&
      obs$df_delta > threshold_log_n_10
    ) {
      "Bubble"
    } else {
      "No bubble"
    }
    
    decision_loglog_n_100 <- if (
      is.finite(obs$df_delta) &&
      obs$df_delta > threshold_loglog_n_100
    ) {
      "Bubble"
    } else {
      "No bubble"
    }
    
    decision_bootstrap <- if (
      is.finite(obs$df_delta) &&
      is.finite(cv90) &&
      obs$df_delta > cv90
    ) {
      "Bubble"
    } else {
      "No bubble"
    }
    
    out_list[[i]] <- data.frame(
      Ticker = ticker,
      Company = unname(COMPANY_LABELS[ticker]),
      Start_Date = as.character(as.Date(start_date)),
      End_Date = as.character(as.Date(end_date)),
      N_Obs = length(X),
      N_Eff = n_eff,
      Delta_Hat = obs$delta_hat,
      DF_Delta = obs$df_delta,
      Threshold_LogN_10 = threshold_log_n_10,
      Threshold_LogLogN_100 = threshold_loglog_n_100,
      Bootstrap_CV_90 = cv90,
      Decision_LogN_10 = decision_log_n_10,
      Decision_LogLogN_100 = decision_loglog_n_100,
      Decision_Bootstrap_90 = decision_bootstrap,
      Bootstrap_P_Value = bootstrap_p_value,
      stringsAsFactors = FALSE
    )
  }
  
  results <- do.call(rbind, out_list)
  rownames(results) <- NULL
  results
}

make_table2_decisions <- function(pre2020_results, post2022_results) {
  pre_tbl <- pre2020_results |>
    select(
      Ticker,
      Company,
      SV_ADF_2018_2020 = Decision_LogN_10,
      PWY_2018_2020 = Decision_LogLogN_100
    )
  
  post_tbl <- post2022_results |>
    select(
      Ticker,
      SV_ADF_2022_2026 = Decision_LogN_10,
      PWY_2022_2026 = Decision_LogLogN_100
    )
  
  table2 <- pre_tbl |>
    left_join(post_tbl, by = "Ticker") |>
    mutate(
      Company = factor(Company, levels = unname(COMPANY_LABELS[ALL_TICKERS]))
    ) |>
    arrange(Company)
  
  table2
}

write_table2_latex_rows <- function(table2, file = TABLE2_LATEX_ROWS) {
  latex_rows <- apply(table2, 1, function(row) {
    paste0(
      row[["Company"]], " & ",
      decision_symbol(row[["SV_ADF_2018_2020"]]), " & ",
      decision_symbol(row[["PWY_2018_2020"]]), " & ",
      decision_symbol(row[["SV_ADF_2022_2026"]]), " & ",
      decision_symbol(row[["PWY_2022_2026"]]), " \\\\"
    )
  })
  
  writeLines(latex_rows, con = file)
  invisible(latex_rows)
}

figure11_output <- make_delta_gamma_dual_plot(
  start_date = "2022-01-01",
  end_date = "2026-03-01",
  lambda = 0.05,
  save_figure = SAVE_FIGURES
)

figure11_estimates <- figure11_output$estimates

decisions_2022_2026 <- bubble_decisions_ar1(
  stock_names = ALL_TICKERS,
  start_date = "2022-01-01",
  end_date = "2026-01-01",
  B = 1000,
  d_par = 1,
  eta_par = 0.1,
  seed = 123,
  verbose = TRUE
)

decisions_2018_2020 <- bubble_decisions_ar1(
  stock_names = ALL_TICKERS,
  start_date = "2018-06-01",
  end_date = "2020-06-01",
  B = 1000,
  d_par = 1,
  eta_par = 0.1,
  seed = 123,
  verbose = TRUE
)

decisions_2024_2026 <- bubble_decisions_ar1(
  stock_names = ALL_TICKERS,
  start_date = "2024-01-01",
  end_date = "2026-01-01",
  B = 1000,
  d_par = 1,
  eta_par = 0.1,
  seed = 123,
  verbose = TRUE
)

table2_decisions <- make_table2_decisions(
  pre2020_results = decisions_2018_2020,
  post2022_results = decisions_2022_2026
)

print(table2_decisions)

if (isTRUE(SAVE_TABLES)) {
  write.csv(table2_decisions, TABLE2_CSV, row.names = FALSE)
  write_table2_latex_rows(table2_decisions, TABLE2_LATEX_ROWS)
  
  write.csv(
    decisions_2022_2026,
    file.path(TABLE_DIR, "bubble_decisions_2022_2026_full.csv"),
    row.names = FALSE
  )
  
  write.csv(
    decisions_2018_2020,
    file.path(TABLE_DIR, "bubble_decisions_2018_2020_full.csv"),
    row.names = FALSE
  )
  
  write.csv(
    decisions_2024_2026,
    file.path(TABLE_DIR, "bubble_decisions_2024_2026_robustness_full.csv"),
    row.names = FALSE
  )
}
