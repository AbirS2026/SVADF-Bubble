
# install.packages(c("quantmod","ggplot2","dplyr"))
library(quantmod)
library(ggplot2)
library(dplyr)

# Date range: last 10 years ############.  Figure 1 ###########################
end_date <- as.Date("2026-04-10")
start_date <- end_date - 10 * 365

# Helper: fetch Adjusted prices safely; fall back to Close if needed
fetch_series <- function(sym, from, to) {
  xt <- tryCatch(
    suppressWarnings(getSymbols(sym, src = "yahoo", from = from, to = to, auto.assign = FALSE)),
    error = function(e) NULL
  )
  if (is.null(xt)) return(NULL)
  adj <- tryCatch(Ad(xt), error = function(e) NULL)
  if (is.null(adj)) adj <- tryCatch(Cl(xt), error = function(e) NULL)
  adj
}

# Choose Nasdaq proxy in order: ^IXIC -> ^NDX -> QQQ
index_candidates <- c("^IXIC", "^NDX", "QQQ")
index_series <- NULL
index_label <- NULL

for (sym in index_candidates) {
  xt <- fetch_series(sym, start_date, end_date)
  if (!is.null(xt) && NROW(xt) > 0) {
    index_series <- xt
    index_label <- switch(sym,
                          "^IXIC" = "Nasdaq",
                          "^NDX"  = "Nasdaq-100",
                          "QQQ"   = "Nasdaq-100 (QQQ)")
    break
  }
}

if (is.null(index_series)) stop("Could not load ^IXIC, ^NDX, or QQQ from Yahoo.")

# Fetch individual stocks
nvda_xt <- fetch_series("NVDA", start_date, end_date)
avgo_xt <- fetch_series("AVGO", start_date, end_date)
tsla_xt <- fetch_series("TSLA", start_date, end_date)

if (is.null(nvda_xt) || is.null(avgo_xt) || is.null(tsla_xt)) {
  stop("Failed to load NVDA, AVGO, or TSLA.")
}

# Build data frames
df_index <- data.frame(
  Date = index(index_series),
  Price = as.numeric(index_series),
  Name = index_label
)

df_nvda <- data.frame(
  Date = index(nvda_xt),
  Price = as.numeric(nvda_xt),
  Name = "Nvidia"
)

df_avgo <- data.frame(
  Date = index(avgo_xt),
  Price = as.numeric(avgo_xt),
  Name = "Broadcom"
)

df_tsla <- data.frame(
  Date = index(tsla_xt),
  Price = as.numeric(tsla_xt),
  Name = "Tesla"
)

plot_df <- bind_rows(df_index, df_nvda, df_avgo, df_tsla) |>
  filter(!is.na(Price))

# Ensure 2x2 facet order
plot_df$Name <- factor(plot_df$Name, levels = c(index_label, "Nvidia", "Broadcom", "Tesla"))

# Distinct colors
cols <- setNames(
  c("#0047AB", "#B30000", "#D4AA00", "#006400"),
  c(index_label, "Nvidia", "Broadcom", "Tesla")
)

# Vertical lines
vlines <- as.Date(c("2020-01-01", "2025-01-01"))

# Plot
p <- ggplot(plot_df, aes(x = Date, y = Price, color = Name)) +
  geom_line(linewidth = 0.9) +
  geom_vline(
    xintercept = vlines,
    linetype = "solid",
    color = "black",
    linewidth = 0.9,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Name, nrow = 2, ncol = 2, scales = "free_y") +
  scale_color_manual(values = cols, guide = "none") +
  scale_x_date(
    date_breaks = "2 years",
    date_labels = "%Y"
  ) +
  labs(
    x = "",
    y = "Index / Stock Price (USD)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    panel.grid.minor = element_blank()
  )

print(p)

# ggsave("nasdaq_nvda_avgo_tsla_2x2.pdf", p, width = 12, height = 8)

##################################################################################

######################################################################################
##Figure 2 ## Rolling Volatility
tickers <- c("AVGO", "NVDA")

# Named display mapping
company_names <- c(
  AVGO = "Broadcom",
  NVDA = "Nvidia"
)

# Rolling window: k_lags + 1 points
k_lags <- 50

# Date range: last 4 years
years_back <- 4
end_date <- as.Date("2026-04-10")
start_date <- end_date - years_back * 365

# ---------------------------
# Data download
# ---------------------------
suppressWarnings(
  getSymbols(tickers, src = "yahoo", from = start_date, to = end_date, auto.assign = TRUE)
)

# ---------------------------
# Helper: rolling standard deviation of price levels
# ---------------------------
compute_sd_price <- function(px_xts, k_lags = 50) {
  rollapplyr(
    data = px_xts,
    width = k_lags + 1,
    FUN = function(x) sd(as.numeric(x), na.rm = TRUE),
    fill = NA,
    align = "right"
  )
}

# ---------------------------
# Build combined data with prices and rolling SD
# ---------------------------
stocks_df <- lapply(tickers, function(tkr) {
  px_xts <- Ad(get(tkr))
  sd50 <- compute_sd_price(px_xts, k_lags = k_lags)
  colnames(sd50) <- "SD50"
  
  df_p <- data.frame(
    Date = index(px_xts),
    Price = as.numeric(px_xts),
    Ticker = tkr,
    stringsAsFactors = FALSE
  )
  
  df_sd <- data.frame(
    Date = index(sd50),
    SD50 = as.numeric(sd50)
  )
  
  df <- left_join(df_p, df_sd, by = "Date")
  df$Company <- unname(company_names[df$Ticker])
  df
}) |>
  bind_rows() |>
  filter(!is.na(Price))

stocks_df$Company[is.na(stocks_df$Company)] <- stocks_df$Ticker[is.na(stocks_df$Company)]
stocks_df$Company <- factor(stocks_df$Company, levels = c("Broadcom", "Nvidia"))

# ---------------------------
# Scale factor to overlay SD on price axis
# ---------------------------
target_frac <- 0.5   # try 0.25 to 0.40 if you want even smaller/larger

price_top <- quantile(stocks_df$Price, 0.98, na.rm = TRUE)
sd_top    <- quantile(stocks_df$SD50, 0.98, na.rm = TRUE)

sf <- if (is.finite(sd_top) && sd_top > 0) {
  target_frac * price_top / sd_top
} else {
  1
}
# ---------------------------
# Colors
# ---------------------------
stock_col <- "#1f77b4"
vol_col   <- "#d62728"

# ---------------------------
# Plot
# ---------------------------
p <- ggplot(stocks_df, aes(x = Date)) +
  geom_line(aes(y = Price, color = "Stock Price (USD)"), linewidth = 0.9) +
  geom_line(
    aes(y = SD50 * sf, color = "Volatility of Stock Price"),
    linewidth = 0.9,
    alpha = 0.85,
    na.rm = TRUE
  ) +
  facet_wrap(~ Company, nrow = 1, ncol = 2, scales = "free_y") +
  scale_y_continuous(
    name = "Stock Price (USD)",
    sec.axis = sec_axis(~ . / sf, name = "Volatility of Stock Price")
  ) +
  scale_x_date(
    breaks = seq.Date(
      from = as.Date(format(min(stocks_df$Date, na.rm = TRUE), "%Y-02-01")),
      to   = as.Date(format(max(stocks_df$Date, na.rm = TRUE), "%Y-12-31")),
      by   = "6 months"
    ),
    date_labels = "%Y %b"
  ) +
  scale_color_manual(
    values = c(
      "Stock Price (USD)" = stock_col,
      "Volatility of Stock Price" = vol_col
    )
  ) +
  labs(
    x = "",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.title.y.left  = element_text(size = 15, color = stock_col),
    axis.text.y.left   = element_text(size = 12, color = stock_col),
    axis.title.y.right = element_text(size = 15, color = vol_col),
    axis.text.y.right  = element_text(size = 12, color = vol_col),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    legend.position = "none",
    legend.text = element_text(size = 13),
    legend.key.width = unit(1.5, "cm"),
    panel.grid.minor = element_blank()
  )

print(p)

# ggsave("broadcom_nvidia_price_rolling_sd_51pt_facet.pdf", p, width = 12, height = 5.5, dpi = 300)
###############################################################################################

#Figure 4  Bitcoin Nvidia PWY Comparison



library(grid)
library(ggh4x)   # install.packages("ggh4x") if needed

panel_info <- data.frame(
  Name   = c("Bitcoin", "Nvidia"),
  Symbol = c("BTC-USD", "NVDA"),
  Start  = as.Date(c("2020-01-01", "2022-08-01")),
  End    = as.Date(c("2021-09-30", "2025-01-31")),
  stringsAsFactors = FALSE
)

# -----------------------------
# Download data
# -----------------------------
plot_list <- lapply(seq_len(nrow(panel_info)), function(i) {
  sym   <- panel_info$Symbol[i]
  nm    <- panel_info$Name[i]
  start <- panel_info$Start[i]
  end   <- panel_info$End[i]
  
  xt <- fetch_series(sym, start, end)
  if (is.null(xt)) stop(paste("Failed to load", sym, "from Yahoo."))
  
  data.frame(
    Date  = index(xt),
    Price = as.numeric(xt),
    Name  = nm,
    stringsAsFactors = FALSE
  )
})

plot_df <- bind_rows(plot_list) %>%
  filter(!is.na(Price))

plot_df$Name <- factor(plot_df$Name, levels = c("Bitcoin", "Nvidia"))

# -----------------------------
# Colors
# -----------------------------
price_cols <- c(
  "Bitcoin" = "#0047AB",
  "Nvidia"  = "#D4AA00"
)

event_cols <- c(
  "origination" = "#B30000",
  "collapse"    = "#006400"
)

# -----------------------------
# Event dates + label positions
# -----------------------------
events <- data.frame(
  Name      = c("Bitcoin", "Bitcoin", "Bitcoin", "Bitcoin",
                "Nvidia", "Nvidia", "Nvidia", "Nvidia"),
  Method    = c("PWY", "PWY", "SV-ADF", "SV-ADF",
                "PWY", "PWY", "SV-ADF", "SV-ADF"),
  EventDate = as.Date(c("2020-11-01", "2020-12-20", "2021-01-01", "2021-05-16",
                        "2023-05-01", "2023-11-01", "2024-01-01", "2024-08-01")),
  EventType = c("origination", "collapse", "origination", "collapse",
                "origination", "collapse", "origination", "collapse"),
  Label     = c("PWY origination", "PWY collapse",
                "SV-ADF origination", "SV-ADF collapse",
                "PWY origination", "PWY collapse",
                "SV-ADF origination", "SV-ADF collapse"),
  x_offset  = c(-75, -140,  70,  85,
                -70,  -40,  55,  65),
  
  y_mult    = c( 0.15,  0.22, -0.14, -0.22,
                 -0.03, -0.18,  0.08, -0.22),
  
  hjust     = c(1, 1, 0, 0,
                1, 1, 0, 0),
  stringsAsFactors = FALSE
)

events$Name <- factor(events$Name, levels = c("Bitcoin", "Nvidia"))

# -----------------------------
# Find nearest trading-day point
# -----------------------------
nearest_point <- function(df, d) {
  idx <- which.min(abs(as.numeric(df$Date - d)))
  df[idx, c("Date", "Price")]
}

# -----------------------------
# Build annotation data
# -----------------------------
ann_list <- lapply(seq_len(nrow(events)), function(i) {
  ev <- events[i, ]
  df_sub <- plot_df %>% filter(Name == ev$Name)
  
  pt <- nearest_point(df_sub, ev$EventDate)
  
  yrng  <- range(df_sub$Price, na.rm = TRUE)
  yspan <- diff(yrng)
  
  data.frame(
    Name      = ev$Name,
    Method    = ev$Method,
    EventDate = ev$EventDate,
    EventType = ev$EventType,
    Label     = ev$Label,
    x_point   = pt$Date,
    y_point   = pt$Price,
    x_text    = pt$Date + ev$x_offset,
    y_text    = pt$Price + ev$y_mult * yspan,
    hjust     = ev$hjust,
    stringsAsFactors = FALSE
  )
})

ann_df <- bind_rows(ann_list)
ann_df$Name <- factor(ann_df$Name, levels = c("Bitcoin", "Nvidia"))

# -----------------------------
# Plot
# -----------------------------
p <- ggplot(plot_df, aes(x = Date, y = Price)) +
  geom_line(aes(color = Name), linewidth = 1.1, show.legend = FALSE) +
  geom_vline(
    data = events,
    aes(xintercept = EventDate, color = EventType),
    linewidth = 1.0,
    linetype = "dashed",
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = ann_df,
    aes(
      x = x_text, y = y_text,
      xend = x_point, yend = y_point,
      color = EventType
    ),
    arrow = arrow(length = unit(0.16, "inches"), type = "closed"),
    linewidth = 0.8,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = ann_df,
    aes(x = x_point, y = y_point, color = EventType),
    size = 2.8,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = ann_df,
    aes(
      x = x_text, y = y_text,
      label = Label, color = EventType, hjust = hjust
    ),
    size = 4.5,
    fontface = "bold",
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  facet_wrap(~ Name, nrow = 1, ncol = 2, scales = "free") +
  facetted_pos_scales(
    x = list(
      # Bitcoin panel
      scale_x_date(
        limits = as.Date(c("2020-01-01", "2021-09-30")),
        date_breaks = "3 months",
        date_labels = "%Y %b",
        expand = expansion(mult = c(0.02, 0.02))
      ),
      # Nvidia panel
      scale_x_date(
        limits = as.Date(c("2023-01-01", "2025-01-31")),
        breaks = seq.Date(
          from = as.Date("2023-02-01"),
          to   = as.Date("2025-01-31"),
          by   = "3 months"
        ),
        date_labels = "%Y %b",
        expand = expansion(mult = c(0.02, 0.02))
      )
    )
  ) +
  scale_color_manual(values = c(price_cols, event_cols)) +
  scale_y_continuous(
    breaks = function(x) sort(unique(c(0, pretty(x, n = 5)))),
    expand = expansion(mult = c(0, 0.06))
  ) +
  labs(x = "", y = "Price / Index (USD)") +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 15),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey85"),
    panel.grid.major.y = element_line(color = "grey88"),
    plot.margin = margin(8, 65, 8, 8)
  )

print(p)  


# ggsave("bitcoin_nvidia_bubble_pwy_svadf_combined.pdf", plot = p, width = 12, height = 5.5, dpi = 300)
###################################################################################
#Figure 3 Nvidia Tesla SDADF Results 



# ----------------------------
# Helper: fetch adjusted prices safely
# ----------------------------
fetch_series <- function(sym, from, to) {
  xt <- tryCatch(
    suppressWarnings(
      getSymbols(sym, src = "yahoo", from = from, to = to, auto.assign = FALSE)
    ),
    error = function(e) NULL
  )
  
  if (is.null(xt)) stop(paste("Could not download", sym))
  
  px <- tryCatch(Ad(xt), error = function(e) NULL)
  if (is.null(px)) px <- Cl(xt)
  na.omit(px)
}

# ----------------------------
# Match original plot ranges
# ----------------------------
nvda_start <- as.Date("2021-10-01")
nvda_end   <- as.Date("2024-09-30")

tsla_start <- as.Date("2018-07-01")
tsla_end   <- as.Date("2021-06-30")

# ----------------------------
# Download data
# ----------------------------
nvda_xt <- fetch_series("NVDA", from = nvda_start, to = nvda_end)
tsla_xt <- fetch_series("TSLA", from = tsla_start, to = tsla_end)

# ----------------------------
# Build plotting data
# ----------------------------
df_nvda <- data.frame(
  Date = index(nvda_xt),
  Price = as.numeric(nvda_xt),
  Name = "Nvidia"
)

df_tsla <- data.frame(
  Date = index(tsla_xt),
  Price = as.numeric(tsla_xt),
  Name = "Tesla"
)

plot_df <- bind_rows(df_nvda, df_tsla)
plot_df$Name <- factor(plot_df$Name, levels = c("Nvidia", "Tesla"))

# ----------------------------
# Bubble event dates + label placement
# ----------------------------
events_df <- data.frame(
  Name = c("Nvidia", "Nvidia", "Tesla", "Tesla"),
  EventDate = as.Date(c(
    "2024-01-01", "2024-08-31",
    "2020-06-15", "2021-03-10"
  )),
  EventType = c("origination", "collapse", "origination", "collapse"),
  Label = c("bubble origination", "collapse",
            "bubble origination", "collapse"),
  
  # Controls horizontal label placement.
  # Negative = text to the left of event point; positive = text to the right.
  x_offset = c(-115, 95, -115, 95),
  
  # Controls vertical label placement relative to the price range.
  y_mult = c(0.20, -0.18, 0.20, -0.16),
  
  # Text alignment.
  hjust = c(1, 0, 1, 0),
  
  stringsAsFactors = FALSE
)

events_df$Name <- factor(events_df$Name, levels = c("Nvidia", "Tesla"))

# ----------------------------
# Colors
# ----------------------------
stock_cols <- c(
  "Nvidia" = "#0047AB",
  "Tesla"  = "#D4AA00"
)

event_cols <- c(
  "origination" = "#B30000",
  "collapse"    = "#006400"
)

# ----------------------------
# Find nearest trading-day point
# ----------------------------
nearest_point <- function(df, d) {
  idx <- which.min(abs(as.numeric(df$Date - d)))
  df[idx, c("Date", "Price")]
}

# ----------------------------
# Build annotation data
# ----------------------------
ann_list <- lapply(seq_len(nrow(events_df)), function(i) {
  ev <- events_df[i, ]
  df_sub <- plot_df %>% filter(Name == ev$Name)
  
  pt <- nearest_point(df_sub, ev$EventDate)
  
  yrng  <- range(df_sub$Price, na.rm = TRUE)
  yspan <- diff(yrng)
  
  data.frame(
    Name      = ev$Name,
    EventDate = ev$EventDate,
    EventType = ev$EventType,
    Label     = ev$Label,
    x_point   = pt$Date,
    y_point   = pt$Price,
    x_text    = pt$Date + ev$x_offset,
    y_text    = pt$Price + ev$y_mult * yspan,
    hjust     = ev$hjust,
    stringsAsFactors = FALSE
  )
})

ann_df <- bind_rows(ann_list)
ann_df$Name <- factor(ann_df$Name, levels = c("Nvidia", "Tesla"))

# ----------------------------
# Plot
# ----------------------------
p <- ggplot(plot_df, aes(x = Date, y = Price)) +
  
  # Stock price line
  geom_line(aes(color = Name), linewidth = 1.1, show.legend = FALSE) +
  
  # Dashed vertical event lines
  geom_vline(
    data = events_df,
    aes(xintercept = EventDate, color = EventType),
    linetype = "dashed",
    linewidth = 1.0,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  
  # Slanted arrows from labels to nearest stock-price point
  geom_segment(
    data = ann_df,
    aes(
      x = x_text, y = y_text,
      xend = x_point, yend = y_point,
      color = EventType
    ),
    arrow = arrow(length = unit(0.16, "inches"), type = "closed"),
    linewidth = 0.8,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  
  # Event points
  geom_point(
    data = ann_df,
    aes(x = x_point, y = y_point, color = EventType),
    size = 2.8,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  
  # Bold labels
  geom_text(
    data = ann_df,
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
  
  scale_color_manual(
    values = c(stock_cols, event_cols)
  ) +
  
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
  
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 15),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey85"),
    panel.grid.major.y = element_line(color = "grey88"),
    panel.spacing = unit(1.2, "lines"),
    plot.margin = margin(10, 55, 10, 10)
  )

print(p)

# ggsave(
#   "nvda_tsla_bubbles_updated.pdf",
#   plot = p,
#   width = 12,
#   height = 5.5,
#   units = "in"
# )


####################################################################################
#Figure 5

library(RColorBrewer)

# Date range: last 10 years
end_date <- Sys.Date()
start_date <- end_date - 10 * 365

# Helper: fetch Adjusted prices safely; fall back to Close if needed
fetch_series <- function(sym, from, to) {
  xt <- tryCatch(
    suppressWarnings(getSymbols(sym, src = "yahoo", from = from, to = to, auto.assign = FALSE)),
    error = function(e) NULL
  )
  if (is.null(xt)) return(NULL)
  adj <- tryCatch(Ad(xt), error = function(e) NULL)
  if (is.null(adj)) adj <- tryCatch(Cl(xt), error = function(e) NULL)
  adj
}

# Tickers and display labels
tickers <- c(
  "AAPL","GOOGL","MSFT","AMZN","META","TSLA",
  "NVDA","TSM","AVGO","ASML","MU","PLTR"
)

labels <- c(
  AAPL  = "Apple",
  GOOGL = "Alphabet",
  MSFT  = "Microsoft",
  AMZN  = "Amazon",
  META  = "Meta",
  TSLA  = "Tesla",
  NVDA  = "Nvidia",
  TSM   = "TSMC",
  AVGO  = "Broadcom",
  ASML  = "ASML",
  MU    = "Micron",
  PLTR  = "Palantir"
)

# Build data frame for all tickers
stock_dfs <- lapply(tickers, function(sym) {
  xt <- fetch_series(sym, start_date, end_date)
  if (is.null(xt) || NROW(xt) == 0) {
    warning(paste("Failed to load", sym))
    return(NULL)
  }
  data.frame(
    Date = index(xt),
    Price = as.numeric(xt),
    Name = unname(labels[[sym]]),
    stringsAsFactors = FALSE
  )
})

stock_dfs <- Filter(Negate(is.null), stock_dfs)
plot_df <- bind_rows(stock_dfs) |>
  filter(!is.na(Price))

# Facet order matches valuation order:
# first 6 companies, then 6 semiconductors
facet_levels <- unname(labels[tickers])
plot_df$Name <- factor(plot_df$Name, levels = facet_levels)

# Distinct colors for 12 panels
cols <- brewer.pal(12, "Set3")
names(cols) <- facet_levels

darken <- function(col, factor = 0.75) {
  rgb_vals <- col2rgb(col)
  rgb(t(pmax(0, rgb_vals * factor)), maxColorValue = 255)
}

cols_deep <- sapply(cols, darken, factor = 0.75)
names(cols_deep) <- names(cols)

# Solid vertical line at 2020-01-01
vlines <- as.Date(c("2020-01-01"))

p <- ggplot(plot_df, aes(x = Date, y = Price, color = Name)) +
  geom_line(linewidth = 1.1) +
  geom_vline(
    xintercept = vlines,
    linetype = "solid",
    color = "#111111",
    linewidth = 1.0
  ) +
  facet_wrap(~ Name, nrow = 4, ncol = 3, scales = "free_y") +
  scale_color_manual(values = cols_deep, guide = "none") +
  scale_x_date(
    date_breaks = "2 years",
    date_labels = "%Y"
  ) +
  labs(
    x = NULL,
    y = "Stock Price (USD)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text = element_text(size = 20, face = "bold"),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 16, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 16),
    panel.grid.major = element_line(color = "grey88"),
    panel.grid.minor = element_line(color = "grey93"),
    panel.spacing = unit(1.0, "lines")
  )

print(p)
# # Save to PDF
# ggsave(
#   "twelve_stocks_4x3_10y_solid_vline.pdf",
#   plot = p,
#   width = 12,
#   height = 16,
#   units = "in"
# )
#################################################################################

#####ONLY VOL NO PRICE 
# Figure  6
# ---------------------------
# Helper: rolling standard deviation of price levels
# ---------------------------
compute_sd_price <- function(px_xts, k_lags = 50) {
  rollapplyr(
    data = px_xts,
    width = k_lags + 1,
    FUN = function(x) sd(as.numeric(x), na.rm = TRUE),
    fill = NA,
    align = "right"
  )
}

# ---------------------------
# Nvidia data
# ---------------------------
ticker <- "NVDA"
company_name <- "Nvidia"
k_lags <- 70

start_date <- as.Date("2022-01-01")
end_date   <- as.Date("2024-12-31")

suppressWarnings(
  getSymbols(ticker, src = "yahoo", from = start_date, to = end_date, auto.assign = TRUE)
)

px_xts <- Ad(get(ticker))
sd50   <- compute_sd_price(px_xts, k_lags = k_lags)
colnames(sd50) <- "SD50"

df <- data.frame(
  Date = index(sd50),
  SD50 = as.numeric(sd50),
  stringsAsFactors = FALSE
) |>
  filter(!is.na(SD50))

# ---------------------------
# Bubble windows
# ---------------------------
bubble_windows <- data.frame(
  xmin  = as.Date(c("2023-05-15", "2024-01-01")),
  xmax  = as.Date(c("2023-11-01", "2024-08-31")),
  fill  = c("grey80", "grey55"),
  label = c("PWY bubble", "SV-ADF bubble"),
  xlab  = as.Date(c("2023-08-10", "2024-05-01")),
  stringsAsFactors = FALSE
)

# label heights
yrng  <- range(df$SD50, na.rm = TRUE)
yspan <- diff(yrng)
y_lab <- max(df$SD50, na.rm = TRUE) - 0.06 * yspan

# ---------------------------
# Plot
# ---------------------------
p <- ggplot(df, aes(x = Date, y = SD50)) +
  annotate(
    "rect",
    xmin = bubble_windows$xmin[1], xmax = bubble_windows$xmax[1],
    ymin = -Inf, ymax = Inf,
    fill = bubble_windows$fill[1], alpha = 0.35
  ) +
  annotate(
    "rect",
    xmin = bubble_windows$xmin[2], xmax = bubble_windows$xmax[2],
    ymin = -Inf, ymax = Inf,
    fill = bubble_windows$fill[2], alpha = 0.35
  ) +
  geom_line(color = "#1f77b4", linewidth = 1.0) +
  annotate(
    "text",
    x = bubble_windows$xlab[1], y = y_lab,
    label = bubble_windows$label[1],
    fontface = "bold", size = 5
  ) +
  annotate(
    "text",
    x = bubble_windows$xlab[2], y = y_lab,
    label = bubble_windows$label[2],
    fontface = "bold", size = 5
  ) +
  scale_x_date(
    date_breaks = "6 months",
    date_labels = "%Y %b"
  ) +
  labs(
    x = "",
    y = "Volatility of Stock Price",
    title = company_name
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 15, color = "#1f77b4"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, color = "#1f77b4"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

print(p)

# ggsave("nvidia_volatility_bubble_windows_only_vol.pdf", p, width = 12, height = 5.5, dpi = 300)

