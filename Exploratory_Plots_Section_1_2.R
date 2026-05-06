################################################################################
# Exploratory_Plots_Section_1_2.R
#
# Replication code for Figures 1--6 of:
# "Is There an AI Bubble? Robust Date-Stamping for Periods of Exuberance"
#
# Data source:
# Yahoo Finance daily adjusted closing prices.
#
# Notes:
# - Package installation is intentionally omitted; see README.md.
# - Set SAVE_FIGURES <- TRUE to write all PDF files to disk.
# - The end date is fixed for reproducibility.
################################################################################

# ------------------------------------------------------------------------------
# Libraries
# ------------------------------------------------------------------------------

library(quantmod)
library(ggplot2)
library(dplyr)
library(zoo)
library(grid)
library(ggh4x)
library(RColorBrewer)

# ------------------------------------------------------------------------------
# Global settings
# ------------------------------------------------------------------------------

SAVE_FIGURES <- TRUE

PAPER_END_DATE <- as.Date("2026-04-10")
OUTPUT_DIR <- "."

if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

# Fetch adjusted prices from Yahoo Finance.
# If Adjusted Close is unavailable, fall back to Close.
fetch_price_series <- function(symbol, from, to) {
  xt <- tryCatch(
    suppressWarnings(
      getSymbols(
        symbol,
        src = "yahoo",
        from = from,
        to = to,
        auto.assign = FALSE
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(xt)) {
    return(NULL)
  }
  
  px <- tryCatch(Ad(xt), error = function(e) NULL)
  
  if (is.null(px)) {
    px <- tryCatch(Cl(xt), error = function(e) NULL)
  }
  
  if (is.null(px)) {
    return(NULL)
  }
  
  na.omit(px)
}

# Convert an xts price series to a clean data frame.
price_to_df <- function(price_xts, name) {
  data.frame(
    Date = as.Date(index(price_xts)),
    Price = as.numeric(price_xts),
    Name = name,
    stringsAsFactors = FALSE
  ) |>
    filter(!is.na(Price))
}

# Rolling standard deviation of price levels.
compute_rolling_sd <- function(price_xts, k_lags = 50) {
  rollapplyr(
    data = price_xts,
    width = k_lags + 1,
    FUN = function(x) sd(as.numeric(x), na.rm = TRUE),
    fill = NA,
    align = "right"
  )
}

# Find the nearest observed trading-day price to a target calendar date.
nearest_price_point <- function(df, target_date) {
  idx <- which.min(abs(as.numeric(df$Date - target_date)))
  df[idx, c("Date", "Price")]
}

# Save a plot only when SAVE_FIGURES is TRUE.
save_figure <- function(filename, plot, width, height, dpi = 300) {
  if (isTRUE(SAVE_FIGURES)) {
    ggsave(
      filename = file.path(OUTPUT_DIR, filename),
      plot = plot,
      width = width,
      height = height,
      dpi = dpi,
      units = "in"
    )
  }
}

# A common clean theme used across the exploratory figures.
theme_paper <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.title.x = element_blank()
    )
}

################################################################################
# Figure 1
# Nasdaq, Nvidia, Broadcom, and Tesla over the past 10 years
################################################################################

fig1_start_date <- PAPER_END_DATE - 10 * 365
fig1_end_date <- PAPER_END_DATE

# Nasdaq proxy priority:
# 1. Nasdaq Composite (^IXIC)
# 2. Nasdaq-100 (^NDX)
# 3. Nasdaq-100 ETF (QQQ)
index_candidates <- c("^IXIC", "^NDX", "QQQ")
index_series <- NULL
index_label <- NULL

for (symbol in index_candidates) {
  candidate_series <- fetch_price_series(symbol, fig1_start_date, fig1_end_date)
  
  if (!is.null(candidate_series) && NROW(candidate_series) > 0) {
    index_series <- candidate_series
    index_label <- switch(
      symbol,
      "^IXIC" = "Nasdaq",
      "^NDX"  = "Nasdaq-100",
      "QQQ"   = "Nasdaq-100 (QQQ)"
    )
    break
  }
}

if (is.null(index_series)) {
  stop("Could not load ^IXIC, ^NDX, or QQQ from Yahoo Finance.")
}

fig1_tickers <- c(
  Nvidia = "NVDA",
  Broadcom = "AVGO",
  Tesla = "TSLA"
)

fig1_stock_list <- lapply(names(fig1_tickers), function(company) {
  symbol <- fig1_tickers[[company]]
  xt <- fetch_price_series(symbol, fig1_start_date, fig1_end_date)
  
  if (is.null(xt) || NROW(xt) == 0) {
    stop(paste("Failed to load", symbol, "from Yahoo Finance."))
  }
  
  price_to_df(xt, company)
})

fig1_df <- bind_rows(
  price_to_df(index_series, index_label),
  bind_rows(fig1_stock_list)
)

fig1_df$Name <- factor(
  fig1_df$Name,
  levels = c(index_label, "Nvidia", "Broadcom", "Tesla")
)

fig1_colors <- setNames(
  c("#0047AB", "#B30000", "#D4AA00", "#006400"),
  c(index_label, "Nvidia", "Broadcom", "Tesla")
)

fig1_vertical_lines <- as.Date(c("2020-01-01", "2025-01-01"))

fig1 <- ggplot(fig1_df, aes(x = Date, y = Price, color = Name)) +
  geom_line(linewidth = 0.9) +
  geom_vline(
    xintercept = fig1_vertical_lines,
    linetype = "solid",
    color = "black",
    linewidth = 0.9,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Name, nrow = 2, ncol = 2, scales = "free_y") +
  scale_color_manual(values = fig1_colors, guide = "none") +
  scale_x_date(
    date_breaks = "2 years",
    date_labels = "%Y"
  ) +
  labs(
    x = NULL,
    y = "Index / Stock Price (USD)"
  ) +
  theme_paper(base_size = 12) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13)
  )

print(fig1)

save_figure(
  filename = "nasdaq_nvda_avgo_tsla_2x2.pdf",
  plot = fig1,
  width = 12,
  height = 8
)

################################################################################
# Figure 2
# Broadcom and Nvidia price levels with rolling price volatility
################################################################################

fig2_tickers <- c(
  Broadcom = "AVGO",
  Nvidia = "NVDA"
)

fig2_start_date <- PAPER_END_DATE - 4 * 365
fig2_end_date <- PAPER_END_DATE
fig2_k_lags <- 50

fig2_list <- lapply(names(fig2_tickers), function(company) {
  symbol <- fig2_tickers[[company]]
  px_xts <- fetch_price_series(symbol, fig2_start_date, fig2_end_date)
  
  if (is.null(px_xts) || NROW(px_xts) == 0) {
    stop(paste("Failed to load", symbol, "from Yahoo Finance."))
  }
  
  sd_xts <- compute_rolling_sd(px_xts, k_lags = fig2_k_lags)
  colnames(sd_xts) <- "RollingSD"
  
  price_df <- data.frame(
    Date = as.Date(index(px_xts)),
    Price = as.numeric(px_xts),
    Company = company,
    stringsAsFactors = FALSE
  )
  
  sd_df <- data.frame(
    Date = as.Date(index(sd_xts)),
    RollingSD = as.numeric(sd_xts),
    stringsAsFactors = FALSE
  )
  
  left_join(price_df, sd_df, by = "Date")
})

fig2_df <- bind_rows(fig2_list) |>
  filter(!is.na(Price))

fig2_df$Company <- factor(fig2_df$Company, levels = c("Broadcom", "Nvidia"))

# Scale rolling volatility onto the left price axis for visual comparison.
target_fraction <- 0.50

price_top <- quantile(fig2_df$Price, 0.98, na.rm = TRUE)
sd_top <- quantile(fig2_df$RollingSD, 0.98, na.rm = TRUE)

scale_factor <- if (is.finite(sd_top) && sd_top > 0) {
  target_fraction * price_top / sd_top
} else {
  1
}

stock_color <- "#1f77b4"
volatility_color <- "#d62728"

fig2 <- ggplot(fig2_df, aes(x = Date)) +
  geom_line(
    aes(y = Price, color = "Stock Price (USD)"),
    linewidth = 0.9
  ) +
  geom_line(
    aes(y = RollingSD * scale_factor, color = "Volatility of Stock Price"),
    linewidth = 0.9,
    alpha = 0.85,
    na.rm = TRUE
  ) +
  facet_wrap(~ Company, nrow = 1, ncol = 2, scales = "free_y") +
  scale_y_continuous(
    name = "Stock Price (USD)",
    sec.axis = sec_axis(
      ~ . / scale_factor,
      name = "Volatility of Stock Price"
    )
  ) +
  scale_x_date(
    date_breaks = "6 months",
    date_labels = "%Y %b"
  ) +
  scale_color_manual(
    values = c(
      "Stock Price (USD)" = stock_color,
      "Volatility of Stock Price" = volatility_color
    ),
    guide = "none"
  ) +
  labs(x = NULL, color = NULL) +
  theme_paper(base_size = 12) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.title.y.left = element_text(size = 15, color = stock_color),
    axis.text.y.left = element_text(size = 12, color = stock_color),
    axis.title.y.right = element_text(size = 15, color = volatility_color),
    axis.text.y.right = element_text(size = 12, color = volatility_color),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1)
  )

print(fig2)

save_figure(
  filename = "broadcom_nvidia_price_rolling_sd_51pt_facet.pdf",
  plot = fig2,
  width = 12,
  height = 5.5
)

################################################################################
# Figure 3
# Nvidia and Tesla bubble episodes: stock-price panels with event annotations
################################################################################

fig3_panel_info <- data.frame(
  Name = c("Nvidia", "Tesla"),
  Symbol = c("NVDA", "TSLA"),
  Start = as.Date(c("2021-10-01", "2018-07-01")),
  End = as.Date(c("2024-09-30", "2021-06-30")),
  stringsAsFactors = FALSE
)

fig3_price_list <- lapply(seq_len(nrow(fig3_panel_info)), function(i) {
  symbol <- fig3_panel_info$Symbol[i]
  company <- fig3_panel_info$Name[i]
  start <- fig3_panel_info$Start[i]
  end <- fig3_panel_info$End[i]
  
  px_xts <- fetch_price_series(symbol, start, end)
  
  if (is.null(px_xts) || NROW(px_xts) == 0) {
    stop(paste("Failed to load", symbol, "from Yahoo Finance."))
  }
  
  price_to_df(px_xts, company)
})

fig3_df <- bind_rows(fig3_price_list)
fig3_df$Name <- factor(fig3_df$Name, levels = c("Nvidia", "Tesla"))

fig3_events <- data.frame(
  Name = c("Nvidia", "Nvidia", "Tesla", "Tesla"),
  EventDate = as.Date(c(
    "2024-01-01", "2024-08-31",
    "2020-06-15", "2021-03-10"
  )),
  EventType = c("origination", "collapse", "origination", "collapse"),
  Label = c(
    "bubble origination", "collapse",
    "bubble origination", "collapse"
  ),
  x_offset = c(-115, 95, -115, 95),
  y_mult = c(0.20, -0.18, 0.20, -0.16),
  hjust = c(1, 0, 1, 0),
  stringsAsFactors = FALSE
)

fig3_events$Name <- factor(fig3_events$Name, levels = c("Nvidia", "Tesla"))

fig3_stock_colors <- c(
  "Nvidia" = "#0047AB",
  "Tesla" = "#D4AA00"
)

fig3_event_colors <- c(
  "origination" = "#B30000",
  "collapse" = "#006400"
)

fig3_annotation_list <- lapply(seq_len(nrow(fig3_events)), function(i) {
  event_i <- fig3_events[i, ]
  df_i <- fig3_df |>
    filter(Name == event_i$Name)
  
  nearest_pt <- nearest_price_point(df_i, event_i$EventDate)
  
  y_range <- range(df_i$Price, na.rm = TRUE)
  y_span <- diff(y_range)
  
  data.frame(
    Name = event_i$Name,
    EventDate = event_i$EventDate,
    EventType = event_i$EventType,
    Label = event_i$Label,
    x_point = nearest_pt$Date,
    y_point = nearest_pt$Price,
    x_text = nearest_pt$Date + event_i$x_offset,
    y_text = nearest_pt$Price + event_i$y_mult * y_span,
    hjust = event_i$hjust,
    stringsAsFactors = FALSE
  )
})

fig3_ann <- bind_rows(fig3_annotation_list)
fig3_ann$Name <- factor(fig3_ann$Name, levels = c("Nvidia", "Tesla"))

fig3 <- ggplot(fig3_df, aes(x = Date, y = Price)) +
  geom_line(aes(color = Name), linewidth = 1.1, show.legend = FALSE) +
  geom_vline(
    data = fig3_events,
    aes(xintercept = EventDate, color = EventType),
    linetype = "dashed",
    linewidth = 1.0,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = fig3_ann,
    aes(
      x = x_text,
      y = y_text,
      xend = x_point,
      yend = y_point,
      color = EventType
    ),
    arrow = arrow(length = unit(0.16, "inches"), type = "closed"),
    linewidth = 0.8,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = fig3_ann,
    aes(x = x_point, y = y_point, color = EventType),
    size = 2.8,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = fig3_ann,
    aes(
      x = x_text,
      y = y_text,
      label = Label,
      color = EventType,
      hjust = hjust
    ),
    size = 4.5,
    fontface = "bold",
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Name, nrow = 1, ncol = 2, scales = "free") +
  scale_color_manual(values = c(fig3_stock_colors, fig3_event_colors)) +
  scale_x_date(
    date_breaks = "6 months",
    date_labels = "%Y %b",
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  scale_y_continuous(
    breaks = function(x) sort(unique(c(0, pretty(x, n = 5)))),
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  labs(
    x = NULL,
    y = "Stock Price (USD)"
  ) +
  coord_cartesian(clip = "off") +
  theme_paper(base_size = 13) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 15),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    panel.grid.major.x = element_line(color = "grey85"),
    panel.grid.major.y = element_line(color = "grey88"),
    panel.spacing = unit(1.2, "lines"),
    plot.margin = margin(10, 55, 10, 10)
  )

print(fig3)

save_figure(
  filename = "nvda_tsla_bubbles_updated.pdf",
  plot = fig3,
  width = 12,
  height = 5.5
)

################################################################################
# Figure 4
# Bitcoin and Nvidia: PWY versus SV-ADF bubble date-stamps
################################################################################

fig4_panel_info <- data.frame(
  Name = c("Bitcoin", "Nvidia"),
  Symbol = c("BTC-USD", "NVDA"),
  Start = as.Date(c("2020-01-01", "2022-08-01")),
  End = as.Date(c("2021-09-30", "2025-01-31")),
  stringsAsFactors = FALSE
)

fig4_price_list <- lapply(seq_len(nrow(fig4_panel_info)), function(i) {
  symbol <- fig4_panel_info$Symbol[i]
  company <- fig4_panel_info$Name[i]
  start <- fig4_panel_info$Start[i]
  end <- fig4_panel_info$End[i]
  
  px_xts <- fetch_price_series(symbol, start, end)
  
  if (is.null(px_xts) || NROW(px_xts) == 0) {
    stop(paste("Failed to load", symbol, "from Yahoo Finance."))
  }
  
  price_to_df(px_xts, company)
})

fig4_df <- bind_rows(fig4_price_list)
fig4_df$Name <- factor(fig4_df$Name, levels = c("Bitcoin", "Nvidia"))

fig4_price_colors <- c(
  "Bitcoin" = "#0047AB",
  "Nvidia" = "#D4AA00"
)

fig4_event_colors <- c(
  "origination" = "#B30000",
  "collapse" = "#006400"
)

fig4_events <- data.frame(
  Name = c(
    "Bitcoin", "Bitcoin", "Bitcoin", "Bitcoin",
    "Nvidia", "Nvidia", "Nvidia", "Nvidia"
  ),
  Method = c(
    "PWY", "PWY", "SV-ADF", "SV-ADF",
    "PWY", "PWY", "SV-ADF", "SV-ADF"
  ),
  EventDate = as.Date(c(
    "2020-11-01", "2020-12-20", "2021-01-01", "2021-05-16",
    "2023-05-01", "2023-11-01", "2024-01-01", "2024-08-01"
  )),
  EventType = c(
    "origination", "collapse", "origination", "collapse",
    "origination", "collapse", "origination", "collapse"
  ),
  Label = c(
    "PWY origination", "PWY collapse",
    "SV-ADF origination", "SV-ADF collapse",
    "PWY origination", "PWY collapse",
    "SV-ADF origination", "SV-ADF collapse"
  ),
  x_offset = c(
    -75, -140, 70, 85,
    -70, -40, 55, 65
  ),
  y_mult = c(
    0.15, 0.22, -0.14, -0.22,
    -0.03, -0.18, 0.08, -0.22
  ),
  hjust = c(
    1, 1, 0, 0,
    1, 1, 0, 0
  ),
  stringsAsFactors = FALSE
)

fig4_events$Name <- factor(fig4_events$Name, levels = c("Bitcoin", "Nvidia"))

fig4_annotation_list <- lapply(seq_len(nrow(fig4_events)), function(i) {
  event_i <- fig4_events[i, ]
  df_i <- fig4_df |>
    filter(Name == event_i$Name)
  
  nearest_pt <- nearest_price_point(df_i, event_i$EventDate)
  
  y_range <- range(df_i$Price, na.rm = TRUE)
  y_span <- diff(y_range)
  
  data.frame(
    Name = event_i$Name,
    Method = event_i$Method,
    EventDate = event_i$EventDate,
    EventType = event_i$EventType,
    Label = event_i$Label,
    x_point = nearest_pt$Date,
    y_point = nearest_pt$Price,
    x_text = nearest_pt$Date + event_i$x_offset,
    y_text = nearest_pt$Price + event_i$y_mult * y_span,
    hjust = event_i$hjust,
    stringsAsFactors = FALSE
  )
})

fig4_ann <- bind_rows(fig4_annotation_list)
fig4_ann$Name <- factor(fig4_ann$Name, levels = c("Bitcoin", "Nvidia"))

fig4 <- ggplot(fig4_df, aes(x = Date, y = Price)) +
  geom_line(aes(color = Name), linewidth = 1.1, show.legend = FALSE) +
  geom_vline(
    data = fig4_events,
    aes(xintercept = EventDate, color = EventType),
    linewidth = 1.0,
    linetype = "dashed",
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = fig4_ann,
    aes(
      x = x_text,
      y = y_text,
      xend = x_point,
      yend = y_point,
      color = EventType
    ),
    arrow = arrow(length = unit(0.16, "inches"), type = "closed"),
    linewidth = 0.8,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = fig4_ann,
    aes(x = x_point, y = y_point, color = EventType),
    size = 2.8,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = fig4_ann,
    aes(
      x = x_text,
      y = y_text,
      label = Label,
      color = EventType,
      hjust = hjust
    ),
    size = 4.5,
    fontface = "bold",
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Name, nrow = 1, ncol = 2, scales = "free") +
  facetted_pos_scales(
    x = list(
      scale_x_date(
        limits = as.Date(c("2020-01-01", "2021-09-30")),
        date_breaks = "3 months",
        date_labels = "%Y %b",
        expand = expansion(mult = c(0.02, 0.02))
      ),
      scale_x_date(
        limits = as.Date(c("2023-01-01", "2025-01-31")),
        breaks = seq.Date(
          from = as.Date("2023-02-01"),
          to = as.Date("2025-01-31"),
          by = "3 months"
        ),
        date_labels = "%Y %b",
        expand = expansion(mult = c(0.02, 0.02))
      )
    )
  ) +
  scale_color_manual(values = c(fig4_price_colors, fig4_event_colors)) +
  scale_y_continuous(
    breaks = function(x) sort(unique(c(0, pretty(x, n = 5)))),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(
    x = NULL,
    y = "Price / Index (USD)"
  ) +
  coord_cartesian(clip = "off") +
  theme_paper(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 15),
    panel.grid.major.x = element_line(color = "grey85"),
    panel.grid.major.y = element_line(color = "grey88"),
    plot.margin = margin(8, 65, 8, 8)
  )

print(fig4)

save_figure(
  filename = "bitcoin_nvidia_bubble_pwy_svadf_combined.pdf",
  plot = fig4,
  width = 12,
  height = 5.5
)

################################################################################
# Figure 5
# Twelve AI-exposed technology and semiconductor stocks over the past 10 years
################################################################################

fig5_start_date <- PAPER_END_DATE - 10 * 365
fig5_end_date <- PAPER_END_DATE

fig5_tickers <- c(
  "AAPL", "GOOGL", "MSFT", "AMZN", "META", "TSLA",
  "NVDA", "TSM", "AVGO", "ASML", "MU", "PLTR"
)

fig5_labels <- c(
  AAPL = "Apple",
  GOOGL = "Alphabet",
  MSFT = "Microsoft",
  AMZN = "Amazon",
  META = "Meta",
  TSLA = "Tesla",
  NVDA = "Nvidia",
  TSM = "TSMC",
  AVGO = "Broadcom",
  ASML = "ASML",
  MU = "Micron",
  PLTR = "Palantir"
)

fig5_price_list <- lapply(fig5_tickers, function(symbol) {
  px_xts <- fetch_price_series(symbol, fig5_start_date, fig5_end_date)
  
  if (is.null(px_xts) || NROW(px_xts) == 0) {
    warning(paste("Failed to load", symbol, "from Yahoo Finance."))
    return(NULL)
  }
  
  price_to_df(px_xts, unname(fig5_labels[[symbol]]))
})

fig5_df <- Filter(Negate(is.null), fig5_price_list) |>
  bind_rows() |>
  filter(!is.na(Price))

fig5_facet_levels <- unname(fig5_labels[fig5_tickers])
fig5_df$Name <- factor(fig5_df$Name, levels = fig5_facet_levels)

fig5_colors <- brewer.pal(12, "Set3")
names(fig5_colors) <- fig5_facet_levels

darken_color <- function(color, factor = 0.75) {
  rgb_values <- col2rgb(color)
  rgb(t(pmax(0, rgb_values * factor)), maxColorValue = 255)
}

fig5_colors_deep <- sapply(fig5_colors, darken_color, factor = 0.75)
names(fig5_colors_deep) <- names(fig5_colors)

fig5 <- ggplot(fig5_df, aes(x = Date, y = Price, color = Name)) +
  geom_line(linewidth = 1.1) +
  geom_vline(
    xintercept = as.Date("2020-01-01"),
    linetype = "solid",
    color = "#111111",
    linewidth = 1.0
  ) +
  facet_wrap(~ Name, nrow = 4, ncol = 3, scales = "free_y") +
  scale_color_manual(values = fig5_colors_deep, guide = "none") +
  scale_x_date(
    date_breaks = "2 years",
    date_labels = "%Y"
  ) +
  labs(
    x = NULL,
    y = "Stock Price (USD)"
  ) +
  theme_paper(base_size = 12) +
  theme(
    strip.text = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 16),
    panel.grid.major = element_line(color = "grey88"),
    panel.grid.minor = element_line(color = "grey93"),
    panel.spacing = unit(1.0, "lines")
  )

print(fig5)

save_figure(
  filename = "twelve_stocks_4x3_10y_solid_vline.pdf",
  plot = fig5,
  width = 12,
  height = 16
)

################################################################################
# Figure 6
# Nvidia volatility with PWY and SV-ADF bubble windows
################################################################################

fig6_ticker <- "NVDA"
fig6_company <- "Nvidia"
fig6_k_lags <- 70

fig6_start_date <- as.Date("2022-01-01")
fig6_end_date <- as.Date("2024-12-31")

fig6_px_xts <- fetch_price_series(fig6_ticker, fig6_start_date, fig6_end_date)

if (is.null(fig6_px_xts) || NROW(fig6_px_xts) == 0) {
  stop("Failed to load NVDA from Yahoo Finance.")
}

fig6_sd_xts <- compute_rolling_sd(fig6_px_xts, k_lags = fig6_k_lags)
colnames(fig6_sd_xts) <- "RollingSD"

fig6_df <- data.frame(
  Date = as.Date(index(fig6_sd_xts)),
  RollingSD = as.numeric(fig6_sd_xts),
  stringsAsFactors = FALSE
) |>
  filter(!is.na(RollingSD))

fig6_windows <- data.frame(
  xmin = as.Date(c("2023-05-15", "2024-01-01")),
  xmax = as.Date(c("2023-11-01", "2024-08-31")),
  fill = c("grey80", "grey55"),
  label = c("PWY bubble", "SV-ADF bubble"),
  xlab = as.Date(c("2023-08-10", "2024-05-01")),
  stringsAsFactors = FALSE
)

fig6_y_range <- range(fig6_df$RollingSD, na.rm = TRUE)
fig6_y_span <- diff(fig6_y_range)
fig6_label_y <- max(fig6_df$RollingSD, na.rm = TRUE) - 0.06 * fig6_y_span

fig6 <- ggplot(fig6_df, aes(x = Date, y = RollingSD)) +
  annotate(
    "rect",
    xmin = fig6_windows$xmin[1],
    xmax = fig6_windows$xmax[1],
    ymin = -Inf,
    ymax = Inf,
    fill = fig6_windows$fill[1],
    alpha = 0.35
  ) +
  annotate(
    "rect",
    xmin = fig6_windows$xmin[2],
    xmax = fig6_windows$xmax[2],
    ymin = -Inf,
    ymax = Inf,
    fill = fig6_windows$fill[2],
    alpha = 0.35
  ) +
  geom_line(color = "#1f77b4", linewidth = 1.0) +
  annotate(
    "text",
    x = fig6_windows$xlab[1],
    y = fig6_label_y,
    label = fig6_windows$label[1],
    fontface = "bold",
    size = 5
  ) +
  annotate(
    "text",
    x = fig6_windows$xlab[2],
    y = fig6_label_y,
    label = fig6_windows$label[2],
    fontface = "bold",
    size = 5
  ) +
  scale_x_date(
    date_breaks = "6 months",
    date_labels = "%Y %b"
  ) +
  labs(
    x = NULL,
    y = "Volatility of Stock Price",
    title = fig6_company
  ) +
  theme_paper(base_size = 12) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.y = element_text(size = 15, color = "#1f77b4"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, color = "#1f77b4"),
    legend.position = "none"
  )

print(fig6)

save_figure(
  filename = "nvidia_volatility_bubble_windows_only_vol.pdf",
  plot = fig6,
  width = 12,
  height = 5.5
)

################################################################################
# End of script
################################################################################
