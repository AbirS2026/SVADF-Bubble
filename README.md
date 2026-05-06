# SV-ADF Bubble Detection Code Repository

This repository contains the replication code for the paper:

**Is There an AI Bubble? Robust Date-Stamping for Periods of Exuberance**

Paper: [arXiv:2604.12062](https://arxiv.org/abs/2604.12062)

The code implements the stochastic-volatility-robust ADF (SV-ADF) bubble detection and date-stamping procedures developed in the paper, together with comparison exercises against the PWY procedure. The scripts generate the figures and tables used in the main paper and the Online Appendix.

---

## Repository Structure

### 1. Exploratory plots for the main empirical narrative

File: [`Exploratory_Plots_Section_1_2.R`](https://github.com/AbirS2026/SVADF-Bubble/blob/main/Exploratory_Plots_Section_1_2.R)

This script generates the exploratory stock-price and volatility plots used in the early empirical sections of the paper. These figures summarize the main empirical motivation: the sharp run-up in AI-exposed technology and semiconductor stocks, together with pronounced time-varying volatility.

Relevant outputs:

- Figures 1--6 of the main paper
- Exploratory plots discussed in Sections 1--3 of the paper

These figures provide the main story of the paper in a nutshell.

---

### 2. Recursive SV-ADF date-stamping for technology stocks

File: [`Recursive SV-ADF Date stamps Tech.R`](https://github.com/AbirS2026/SVADF-Bubble/blob/main/Recursive%20SV-ADF%20Date%20stamps%20Tech.R)

This script applies the recursive SV-ADF date-stamping procedure to large technology firms. It identifies bubble origination and collapse dates for the technology-stock sample.

Relevant paper output:

- Figure 7 of the main paper

---

### 3. Recursive SV-ADF date-stamping for semiconductor and AI-infrastructure stocks

File: [`Recursive SV-ADF Date stamps Semis.R`](https://github.com/AbirS2026/SVADF-Bubble/blob/main/Recursive%20SV-ADF%20Date%20stamps%20Semis.R)

This script applies the recursive SV-ADF date-stamping procedure to semiconductor and AI-infrastructure firms, including Nvidia, TSMC, Broadcom, ASML, Micron, and Palantir.

Relevant paper output:

- Figure 8 of the main paper

---

### 4. Nasdaq, Bitcoin, and Ethereum bubble analysis

File: [`Nasdaq_crypto_bubbles.R`](https://github.com/AbirS2026/SVADF-Bubble/blob/main/Nasdaq_crypto_bubbles.R)

This script applies the recursive SV-ADF and PWY comparison framework to Nasdaq, Bitcoin, and Ethereum. It illustrates how accounting for stochastic volatility changes the dating of bubble origination and collapse relative to homoskedastic PWY procedures.

Relevant paper outputs:

- Figure 9 of the main paper
- Figure 10 of the main paper

---

### 5. Threshold-selection simulations

File: [`Threshold selections.R`](https://github.com/AbirS2026/SVADF-Bubble/blob/main/Threshold%20selections.R)

This script generates the simulation evidence used to motivate and calibrate the origination and collapse thresholds used by the SV-ADF procedure.

Relevant outputs:

- Figure 12 of the main paper
- Figure A.3 of the Online Appendix
- Table A.1 of the Online Appendix

---

### 6. Optimal bubble-window selection

File: [`Bubble window longest.R`](https://github.com/AbirS2026/SVADF-Bubble/blob/main/Bubble%20window%20longest.R)

This script implements the rolling-window search procedure used to select the optimal bubble window for a given asset over a specified time frame. The criterion selects the window with the most persistent detected bubble episode, measured by bubble duration or intensity across candidate rolling windows.

This file is useful for applying the proposed procedure to any individual stock or asset over a user-specified period.

---

### 7. Bubble classification across time periods for the 12-stock sample

File: [`Bubble_timeperiod_12stocks.R`](https://github.com/AbirS2026/SVADF-Bubble/blob/main/Bubble_timeperiod_12stocks.R)

This script compares bubble classifications across the two main empirical subperiods considered in the paper. It applies the SV-ADF and PWY procedures to the twelve AI-exposed firms and produces the period-level bubble/no-bubble classifications.

Relevant paper outputs:

- Figure 11 of the main paper
- Table 2 of the main paper

---

### 8. Additional crypto and PWY-assumption counterexamples

File: [`Appendix Crypto.R`](https://github.com/AbirS2026/SVADF-Bubble/blob/main/Appendix%20Crypto.R)

This script generates additional Online Appendix results. It includes simulation evidence illustrating why the PWY maintained assumption on the reset magnitude can be problematic, and it also provides additional cryptocurrency applications beyond Bitcoin and Ethereum.

Relevant Online Appendix outputs:

- Figure A.1: simulation evidence against the PWY reset assumption
- Figure A.2: additional cryptocurrency bubble-detection results

---

### 9. Table-generating code

File: [`Tables_generating_code.R`](https://github.com/AbirS2026/SVADF-Bubble/blob/main/Tables_generating_code.R)

This script generates the simulation tables reported in the main paper. It contains the Monte Carlo experiments comparing the SV-ADF procedure with the PWY benchmark in terms of identification power, collapse detection, and estimation accuracy.

Relevant paper outputs:

- Table 3 of the main paper
- Table 4 of the main paper
- Table 5 of the main paper

---

## How to Use the Code

Each script can be run independently in R. The empirical scripts download daily adjusted closing prices from Yahoo Finance, usually through the `quantmod` package.

A typical workflow is:

```r
source("Recursive SV-ADF Date stamps Tech.R")
```

For the semiconductor and AI-infrastructure date-stamping results:

```r
source("Recursive SV-ADF Date stamps Semis.R")
```

For the Nasdaq, Bitcoin, and Ethereum results:

```r
source("Nasdaq_crypto_bubbles.R")
```

For the table simulations:

```r
source("Tables_generating_code.R")
```

Some scripts generate PDF figures directly using `ggsave()` or `pdf()`. Output files are saved in the working directory unless otherwise specified.

---

## Required R Packages

A comprehensive package set for running the scripts is:

```r
install.packages(c(
  "quantmod",
  "ggplot2",
  "dplyr",
  "zoo",
  "xts",
  "RColorBrewer",
  "gridExtra",
  "ggh4x",
  "scales",
  "lubridate",
  "purrr",
  "stringr",
  "tidyr",
  "readr",
  "httr",
  "jsonlite"
))
```

Some scripts may only require a subset of these packages. The `grid` package is also used in some plotting routines, but it is included with base R and usually does not need to be installed separately.

---

## Data Source

The empirical analysis uses daily adjusted closing prices from Yahoo Finance:

[https://finance.yahoo.com](https://finance.yahoo.com)

Market-cap information reported in the paper is also based on Yahoo Finance.

Because Yahoo Finance data may be updated or revised over time, exact numerical outputs may differ slightly depending on the date on which the scripts are run.

---

## Citation

If you use this code, please cite the paper:

Sarkar, A. and Wells, M. T.  
**Is There an AI Bubble? Robust Date-Stamping for Periods of Exuberance.**  
arXiv:2604.12062.  
[https://arxiv.org/abs/2604.12062](https://arxiv.org/abs/2604.12062)

BibTeX entry for the paper:

```bibtex
@misc{sarkarwells2026aibubble,
  author       = {Sarkar, Abir and Wells, Martin T.},
  title        = {Is There an AI Bubble? Robust Date-Stamping for Periods of Exuberance},
  year         = {2026},
  eprint       = {2604.12062},
  archivePrefix = {arXiv},
  primaryClass = {econ.EM},
  url          = {https://arxiv.org/abs/2604.12062}
}
```

The empirical data source is Yahoo Finance:

Yahoo! Finance.  
**Yahoo historical data.**  
[https://finance.yahoo.com](https://finance.yahoo.com)

BibTeX entry for the Yahoo Finance data source:

```bibtex
@online{yahooNVDA,
  author  = {{Yahoo! Finance}},
  title   = {Yahoo historical data},
  year    = {2026},
  url     = {https://finance.yahoo.com},
  urldate = {2026-03-15}
}
```

---

## Notes

The code is organized to make each figure and table reproducible from a corresponding R script. File names in this repository are chosen to match the empirical and simulation components of the paper as closely as possible.

Yahoo Finance data are downloaded dynamically by the empirical scripts. Since Yahoo Finance may update historical prices, adjusted prices, and market-cap fields over time, small differences from the paper’s reported numbers may occur if the scripts are run at a later date.
