
make_delta_gamma_dual_plot <- function(start_date = "2022-01-01",
                                       end_date   = "2026-03-01",
                                       lambda     = 0.05) {
  # -------------------------
  # Tickers and labels
  # -------------------------
  tickers <- c("AAPL","MSFT","NVDA","GOOGL","AMZN","META",
               "TSLA","TSM","AVGO","PLTR","ASML","MU")
  
  labels <- c(
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
  
  start_date <- as.Date(start_date)
  end_date   <- as.Date(end_date)
  
  # -------------------------
  # Critical values
  # -------------------------
  cv_lambda  <- qnorm(1 - lambda / 2)
  cvC_lambda <- qcauchy(1 - lambda / 2)
  
  # -------------------------
  # Download data
  # -------------------------
  suppressWarnings(
    getSymbols(tickers, src = "yahoo", from = start_date, to = end_date, auto.assign = TRUE)
  )
  
  # -------------------------
  # Compute estimates and CIs
  # -------------------------
  compute_stats_and_cis <- function(price_xts, lambda = 0.05) {
    X <- as.numeric(na.omit(Ad(price_xts)))
    
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
    
    y    <- X[-1]
    xlag <- X[-length(X)]
    Tn   <- length(y)
    
    delta_hat <- sum(y * xlag) / sum(xlag^2)
    
    diff_from_one <- abs(delta_hat - 1)
    if (diff_from_one <= .Machine$double.eps) {
      gamma_hat <- NA_real_
    } else {
      gamma_hat <- -log(diff_from_one) / log(Tn)
    }
    
    delta_lower <- delta_upper <- NA_real_
    gamma_lower <- gamma_upper <- NA_real_
    regime_code <- NA_real_
    
    if (!is.na(delta_hat) && !is.na(gamma_hat)) {
      if (delta_hat < 1) {
        regime_code <- 0
        
        delta_half_width <- cv_lambda * 2 / (Tn^((1 + gamma_hat) / 2))
        delta_lower <- delta_hat - delta_half_width
        delta_upper <- delta_hat + delta_half_width
        
        gamma_half_width <- cv_lambda * sqrt(2) / ((Tn^((1 - gamma_hat) / 2)) * log(Tn))
        gamma_lower <- gamma_hat - gamma_half_width
        gamma_upper <- gamma_hat + gamma_half_width
        
      } else if (delta_hat > 1) {
        regime_code <- 1
        
        delta_half_width <- cvC_lambda * 2 / ((Tn^gamma_hat) * (delta_hat^Tn))
        delta_lower <- delta_hat - delta_half_width
        delta_upper <- delta_hat + delta_half_width
        
        gamma_half_width <- cvC_lambda * 2 / (((1 + 1 / (Tn^gamma_hat))^Tn) * log(Tn))
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
  
  est_mat <- t(sapply(tickers, function(tkr) compute_stats_and_cis(get(tkr), lambda = lambda)))
  
  results_df <- data.frame(
    Ticker = tickers,
    Company = unname(labels[tickers]),
    T = est_mat[, "T"],
    Delta_Hat = est_mat[, "delta_hat"],
    Gamma_Hat = est_mat[, "gamma_hat"],
    Delta_CI_Lower = est_mat[, "delta_lower"],
    Delta_CI_Upper = est_mat[, "delta_upper"],
    Gamma_CI_Lower = est_mat[, "gamma_lower"],
    Gamma_CI_Upper = est_mat[, "gamma_upper"],
    Regime = ifelse(est_mat[, "regime_code"] == 0, "MIR (<1)",
                    ifelse(est_mat[, "regime_code"] == 1, "MER (>1)", "Boundary (=1)")),
    row.names = NULL
  )
  
  results_df$Company <- factor(results_df$Company, levels = unname(labels[tickers]))
  
  # -------------------------
  # Colors
  # -------------------------
  delta_bar_col  <- "#0047AB"  # blue
  delta_dot_col  <- "#d62728"  # red
  gamma_bar_col  <- "black"
  gamma_dot_col  <- "#2E8B57"  # yellow
  
  # -------------------------
  # Separate vertical bands so the two CI sets never intersect
  # -------------------------
  dmin <- min(results_df$Delta_CI_Lower, results_df$Delta_Hat, na.rm = TRUE)
  dmax <- max(results_df$Delta_CI_Upper, results_df$Delta_Hat, na.rm = TRUE)
  
  gmin <- min(results_df$Gamma_CI_Lower, results_df$Gamma_Hat, na.rm = TRUE)
  gmax <- max(results_df$Gamma_CI_Upper, results_df$Gamma_Hat, na.rm = TRUE)
  
  # lower band for delta, upper band for gamma
  delta_band_low   <- 0.08
  delta_band_high  <- 0.42
  gamma_band_low   <- 0.58
  gamma_band_high  <- 0.92
  
  scale_to_band <- function(x, xmin, xmax, low, high) {
    if (!is.finite(xmin) || !is.finite(xmax) || xmax == xmin) {
      return(rep((low + high) / 2, length(x)))
    }
    low + (x - xmin) * (high - low) / (xmax - xmin)
  }
  
  results_df$Delta_Hat_plot      <- scale_to_band(results_df$Delta_Hat,      dmin, dmax, delta_band_low, delta_band_high)
  results_df$Delta_CI_Lower_plot <- scale_to_band(results_df$Delta_CI_Lower, dmin, dmax, delta_band_low, delta_band_high)
  results_df$Delta_CI_Upper_plot <- scale_to_band(results_df$Delta_CI_Upper, dmin, dmax, delta_band_low, delta_band_high)
  
  results_df$Gamma_Hat_plot      <- scale_to_band(results_df$Gamma_Hat,      gmin, gmax, gamma_band_low, gamma_band_high)
  results_df$Gamma_CI_Lower_plot <- scale_to_band(results_df$Gamma_CI_Lower, gmin, gmax, gamma_band_low, gamma_band_high)
  results_df$Gamma_CI_Upper_plot <- scale_to_band(results_df$Gamma_CI_Upper, gmin, gmax, gamma_band_low, gamma_band_high)
  
  # pretty breaks for each axis
  delta_breaks_vals <- pretty(c(dmin, dmax), n = 5)
  gamma_breaks_vals <- pretty(c(gmin, gmax), n = 5)
  
  delta_breaks_pos <- scale_to_band(delta_breaks_vals, dmin, dmax, delta_band_low, delta_band_high)
  gamma_breaks_pos <- scale_to_band(gamma_breaks_vals, gmin, gmax, gamma_band_low, gamma_band_high)
  delta_one_pos <- scale_to_band(1, dmin, dmax, delta_band_low, delta_band_high)
  gamma_one_pos <- scale_to_band(1, gmin, gmax, gamma_band_low, gamma_band_high)
  
  # -------------------------
  # Plot
  # -------------------------
  p <- ggplot(results_df, aes(x = Company)) +
    # delta: blue CI bars + red dots
    geom_errorbar(
      aes(ymin = Delta_CI_Lower_plot, ymax = Delta_CI_Upper_plot),
      color = delta_bar_col,
      width = 0.18,
      linewidth = 0.8
    ) +
    geom_point(
      aes(y = Delta_Hat_plot),
      color = delta_dot_col,
      size = 2.5
    ) +
    geom_hline(yintercept = delta_one_pos, linetype = "dashed", color = delta_bar_col) +
    geom_hline(yintercept = gamma_one_pos, linetype = "dashed", color = gamma_bar_col) +
    
    # gamma: black CI bars + yellow dots
    geom_errorbar(
      aes(ymin = Gamma_CI_Lower_plot, ymax = Gamma_CI_Upper_plot),
      color = gamma_bar_col,
      width = 0.18,
      linewidth = 0.8
    ) +
    geom_point(
      aes(y = Gamma_Hat_plot),
      color = gamma_dot_col,
      size = 2.5
    ) +
    
    scale_y_continuous(
      limits = c(0, 1),
      name =expression("Estimated AR coefficient and C.I. for " * delta[n]),
      breaks = delta_breaks_pos,
      labels = round(delta_breaks_vals, 3),
      sec.axis = sec_axis(
        ~ .,
        name = expression("Estimated explosive growth rate and C.I. for " * gamma[n]),
        breaks = gamma_breaks_pos,
        labels = round(gamma_breaks_vals, 3)
      )
    ) +
    labs(
      x = "",
      title = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y.left  = element_text(size = 20, color = delta_bar_col),
      axis.title.y.right = element_text(size = 20, color = gamma_bar_col),
      axis.text.x = element_text(size = 18, angle = 45, hjust = 1, face = "bold"),
      axis.text.y.left  = element_text(size = 20, color = delta_bar_col),
      axis.text.y.right = element_text(size = 20, color = gamma_bar_col),
      axis.ticks.y.left  = element_line(color = delta_bar_col),
      axis.ticks.y.right = element_line(color = gamma_bar_col),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  
  print(p)
  ggsave("delta_gamma_dual_axis_plot.pdf", plot = p, width = 12, height = 10)
  
  return(results_df)
}

make_delta_gamma_dual_plot(
  start_date = "2022-01-01",
  end_date   = "2026-03-01",
  lambda     = 0.05
)

########################
# 2 time perios bubbles vs no bubble



bubble_decisions_ar1 <- function(stock_names,
                                 start_date,
                                 end_date,
                                 B = 1000,
                                 d_par = 1,
                                 eta_par = 0.5,
                                 x0_boot = 1,
                                 sigma2_0_boot = 1,
                                 normalize_start = TRUE,
                                 seed = NULL,
                                 verbose = TRUE) {
  # ----------------------------
  # Packages
  # ----------------------------
  if (!requireNamespace("quantmod", quietly = TRUE)) {
    stop("Please install the 'quantmod' package.")
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  # ----------------------------
  # Default ticker set + labels
  # ----------------------------
  default_tickers <- c("AAPL","MSFT","NVDA","GOOGL","AMZN","META",
                       "TSLA","TSM","AVGO","PLTR","ASML","MU")
  
  labels <- c(
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
  
  # Allow either ticker or company name
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
  
  resolve_ticker <- function(x) {
    x_trim <- trimws(x)
    x_up   <- toupper(x_trim)
    x_low  <- tolower(x_trim)
    
    if (x_up %in% names(labels)) return(x_up)
    if (x_low %in% names(name_to_ticker)) return(unname(name_to_ticker[x_low]))
    
    stop(sprintf("Unknown stock name/ticker: '%s'", x))
  }
  
  tickers <- vapply(stock_names, resolve_ticker, character(1))
  
  # ----------------------------
  # Helper: compute DF_delta
  # Model: X_t = delta X_{t-1} + u_t
  # OLS without intercept
  # DF_delta = n * (delta_hat - 1)
  # ----------------------------
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
    
    Y    <- X[-1]
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
    df_delta  <- n_eff * (delta_hat - 1)
    
    list(
      n_eff = n_eff,
      delta_hat = delta_hat,
      df_delta = df_delta
    )
  }
  
  # ----------------------------
  # Helper: simulate under H0: delta = 1
  # with stochastic volatility:
  # log sigma_t^2 = phi_n log sigma_{t-1}^2 + eta_t
  # ----------------------------
  simulate_null_unit_root_sv <- function(n_eff,
                                         d_par = 1,
                                         eta_par = 0.1,
                                         x0 = 1,
                                         sigma2_0 = 1) {
    if (!is.finite(n_eff) || n_eff < 20) {
      stop("n_eff must be at least 20 for the log(log n) volatility specification.")
    }
    
    ll_n <- log(log(n_eff))
    if (!is.finite(ll_n) || ll_n == 0) {
      stop("log(log(n_eff)) is not finite; increase sample size.")
    }
    
    phi_n <- 1 - d_par / ll_n
    
    X <- numeric(n_eff + 1)
    X[1] <- x0
    
    logsig2_prev <- log(sigma2_0)
    
    for (t in seq_len(n_eff)) {
      eta_t     <- rnorm(1, mean = 0, sd = eta_par)
      logsig2_t <- phi_n * logsig2_prev + eta_t
      sigma_t   <- sqrt(exp(logsig2_t))
      u_t       <- sigma_t * rnorm(1)
      
      # H0: delta = 1 always
      X[t + 1] <- X[t] + u_t
      logsig2_prev <- logsig2_t
    }
    
    X
  }
  
  # ----------------------------
  # Helper: bootstrap 95% CV for given n_eff
  # Caches by n_eff so if several firms have same sample size
  # we do not recompute
  # ----------------------------
  boot_cache <- new.env(parent = emptyenv())
  
  bootstrap_cv_95 <- function(n_eff,
                              B = 1000,
                              d_par = 1,
                              eta_par = 0.1,
                              x0 = 1,
                              sigma2_0 = 1) {
    key <- paste(n_eff, B, d_par, eta_par, x0, sigma2_0, sep = "|")
    
    if (exists(key, envir = boot_cache, inherits = FALSE)) {
      return(get(key, envir = boot_cache, inherits = FALSE))
    }
    
    vals <- numeric(B)
    
    for (b in seq_len(B)) {
      X_star <- simulate_null_unit_root_sv(
        n_eff = n_eff,
        d_par = d_par,
        eta_par = eta_par,
        x0 = x0,
        sigma2_0 = sigma2_0
      )
      
      vals[b] <- compute_df_delta(X_star)$df_delta
    }
    
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0) {
      out <- list(cv_95 = NA_real_, boot_vals = numeric(0))
    } else {
      out <- list(
        cv_95 = unname(stats::quantile(vals, probs = 0.90, na.rm = TRUE, type = 7)),
        boot_vals = vals
      )
    }
    
    assign(key, out, envir = boot_cache)
    out
  }
  
  # ----------------------------
  # Main loop over stocks
  # ----------------------------
  out_list <- vector("list", length(tickers))
  
  for (i in seq_along(tickers)) {
    tkr <- tickers[i]
    
    if (verbose) {
      message(sprintf("Processing %s (%d/%d) ...", tkr, i, length(tickers)))
    }
    
    px_xts <- tryCatch(
      quantmod::getSymbols(
        Symbols = tkr,
        src = "yahoo",
        from = as.Date(start_date),
        to   = as.Date(end_date),
        auto.assign = FALSE,
        warnings = FALSE
      ),
      error = function(e) NULL
    )
    
    if (is.null(px_xts)) {
      out_list[[i]] <- data.frame(
        Ticker = tkr,
        Company = unname(labels[tkr]),
        Start_Date = as.character(as.Date(start_date)),
        End_Date = as.character(as.Date(end_date)),
        N_Obs = NA_integer_,
        N_Eff = NA_integer_,
        Delta_Hat = NA_real_,
        DF_Delta = NA_real_,
        Threshold_LogN_10 = NA_real_,
        Threshold_LogLogN_100 = NA_real_,
        Bootstrap_CV_95 = NA_real_,
        Decision_LogN_10 = NA_character_,
        Decision_LogLogN_100 = NA_character_,
        Decision_Bootstrap_95 = NA_character_,
        Bootstrap_P_Value = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    
    X <- as.numeric(quantmod::Ad(px_xts))
    X <- X[is.finite(X)]
    
    if (length(X) < 3) {
      out_list[[i]] <- data.frame(
        Ticker = tkr,
        Company = unname(labels[tkr]),
        Start_Date = as.character(as.Date(start_date)),
        End_Date = as.character(as.Date(end_date)),
        N_Obs = length(X),
        N_Eff = NA_integer_,
        Delta_Hat = NA_real_,
        DF_Delta = NA_real_,
        Threshold_LogN_10 = NA_real_,
        Threshold_LogLogN_100 = NA_real_,
        Bootstrap_CV_95 = NA_real_,
        Decision_LogN_10 = NA_character_,
        Decision_LogLogN_100 = NA_character_,
        Decision_Bootstrap_95 = NA_character_,
        Bootstrap_P_Value = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }
    
    # Normalize to start at 1 so the observed series and bootstrap null
    # are on a comparable scale, matching your simulation setup
    if (normalize_start) {
      X <- X / X[1]
    }
    
    obs <- compute_df_delta(X)
    n_eff <- obs$n_eff
    
    thr1 <- log(n_eff) / 10
    thr2 <- log(log(n_eff)) / 100
    
    boot_obj <- bootstrap_cv_95(
      n_eff = n_eff,
      B = B,
      d_par = d_par,
      eta_par = eta_par,
      x0 = x0_boot,
      sigma2_0 = sigma2_0_boot
    )
    
    cv95 <- boot_obj$cv_95
    boot_vals <- boot_obj$boot_vals
    
    p_boot <- if (length(boot_vals) > 0 && is.finite(obs$df_delta)) {
      mean(boot_vals >= obs$df_delta)
    } else {
      NA_real_
    }
    
    dec1 <- if (is.finite(obs$df_delta) && obs$df_delta > thr1) "Bubble" else "No bubble"
    dec2 <- if (is.finite(obs$df_delta) && obs$df_delta > thr2) "Bubble" else "No bubble"
    dec3 <- if (is.finite(obs$df_delta) && is.finite(cv95) && obs$df_delta > cv95) "Bubble" else "No bubble"
    
    out_list[[i]] <- data.frame(
      Ticker = tkr,
      Company = unname(labels[tkr]),
      Start_Date = as.character(as.Date(start_date)),
      End_Date = as.character(as.Date(end_date)),
      N_Obs = length(X),
      N_Eff = n_eff,
      Delta_Hat = obs$delta_hat,
      DF_Delta = obs$df_delta,
      Threshold_LogN_10 = thr1,
      Threshold_LogLogN_100 = thr2,
      Bootstrap_CV_95 = cv95,
      Decision_LogN_10 = dec1,
      Decision_LogLogN_100 = dec2,
      Decision_Bootstrap_95 = dec3,
      Bootstrap_P_Value = p_boot,
      stringsAsFactors = FALSE
    )
  }
  
  results <- do.call(rbind, out_list)
  rownames(results) <- NULL
  results
}

all_12 <- c("AAPL","MSFT","NVDA","GOOGL","AMZN","META",
            "TSLA","TSM","AVGO","PLTR","ASML","MU")

#Post 2022 - Nvidia, Alphabet, TSMC, Micron
res_all <- bubble_decisions_ar1(
  stock_names = all_12,
  start_date  = "2022-01-01",
  end_date    = "2026-01-01",
  B           = 1000,
  d_par       = 1,
  eta_par     = 0.1,
  seed        = 123
)

print(res_all[, c("Ticker", "Company",
                  "DF_Delta",
                  "Decision_LogN_10",
                  "Decision_LogLogN_100",
                  "Decision_Bootstrap_95")])

#Post 2024 - Nvidia, Alphabet, TSMC, Micron
res_all <- bubble_decisions_ar1(
  stock_names = all_12,
  start_date  = "2024-01-01",
  end_date    = "2026-01-01",
  B           = 1000,
  d_par       = 1,
  eta_par     = 0.1,
  seed        = 123
)
print(res_all[, c("Ticker", "Company",
                  "DF_Delta",
                  "Decision_LogN_10",
                  "Decision_LogLogN_100",
                  "Decision_Bootstrap_95")])

#Pre 2020 - Tesla only significant!
res_all <- bubble_decisions_ar1(
  stock_names = all_12,
  start_date  = "2018-06-01",
  end_date    = "2020-06-01",
  B           = 1000,
  d_par       = 1,
  eta_par     = 0.1,
  seed        = 123
)
print(res_all[, c("Ticker", "Company",
                  "DF_Delta",
                  "Decision_LogN_10",
                  "Decision_LogLogN_100",
                  "Decision_Bootstrap_95")])
