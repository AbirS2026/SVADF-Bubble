# SV-ADF Bubble Detection Code Repository

This repository contains the replication code for the paper:

**Is There an AI Bubble? Robust Date-Stamping for Periods of Exuberance**

Paper: [arXiv:2604.12062](https://arxiv.org/abs/2604.12062)

The code implements the stochastic-volatility-robust ADF (SV-ADF) bubble detection and date-stamping procedures developed in the paper, together with comparison exercises against the PWY procedure. The scripts generate the figures and tables used in the main paper and Online Appendix.

---

## Repository Structure

### 1. Exploratory plots for the main empirical narrative

File: [`Exploratory_Plots_Section_1_2.R`](https://github.com/AbirS2026/SVADF-Bubble/blob/main/Exploratory_Plots_Section_1_2.R)

This script generates the exploratory stock-price and volatility plots used in the early empirical sections of the paper. These figures summarize the main empirical motivation: the sharp run-up in AI-exposed technology and semiconductor stocks, together with pronounced time-varying volatility.

Relevant outputs include the plots discussed in Sections 1--3 of the paper.

---

### 2. Recursive SV-ADF date-stamping for technology stocks

File: [`Recursive SV-ADF Date stamps Tech.R`](https://github.com/AbirS2026/SVADF-Bubble/blob/main/Recursive%20SV-ADF%20Date%20stamps%20Tech.R)

This script applies the recursive SV-ADF date-stamping procedure to large technology firms. It identifies bubble origination and collapse dates for the technology-stock sample.

This file is used for the technology-stock panels in the main empirical date-stamping results.

Relevant paper output:

- Figure 7 of the paper

---

### 3. Recursive SV-ADF date-stamping for semiconductor and AI-infrastructure stocks

File: [`Recursive SV-ADF Date stamps Semis.R`](https://github.com/AbirS2026/SVADF-Bubble/blob/main/Recursive%20SV-ADF%20Date%20stamps%20Semis.R)

This script applies the recursive SV-ADF date-stamping procedure to semiconductor and AI-infrastructure firms, including Nvidia, TSMC, Broadcom, ASML, Micron, and Palantir.

Relevant paper output:

- Figure 8 of the paper

---

### 4. Nasdaq, Bitcoin, and Ethereum bubble analysis

File: [`Nasdaq_crypto_bubbles.R`](https://github.com/AbirS2026/SVADF-Bubble/blob/main/Nasdaq_crypto_bubbles.R)

This script applies the recursive SV-ADF and PWY comparison framework to Nasdaq, Bitcoin, and Ethereum. It illustrates how accounting for stochastic volatility changes the dating of bubble origination and collapse relative to homoskedastic PWY procedures.

Relevant paper outputs:

- Figure 9 of the paper
- Figure 10 of the paper

---

### 5. Threshold-selection simulations

File: [`Threshold selections.R`](https://github.com/AbirS2026/SVADF-Bubble/blob/main/Threshold%20selections.R)

This script generates the simulation evidence used to motivate and calibrate the origination and collapse thresholds used by the SV-ADF procedure.

Relevant outputs:

- Figure 12 of the paper
- Figure A.3 of the Online Appendix
- Table A.1 of the Online Appendix

---

### 6. Optimal bubble-window selection

File: [`Bubble window longest.R`](https://github.com/AbirS2026/SVADF-Bubble/blob/main/Bubble%20window%20longest.R)

This script implements the rolling-window search procedure used to select the optimal bubble window for a given asset over a specified time frame. The criterion selects the window with the most persistent detected bubble episode, measured by bubble duration or intensity over candidate rolling windows.

This file is useful for applying the proposed procedure to any individual stock or asset over a user-specified period.

---

### 7. Bubble classification across time periods for the 12-stock sample

File: [`Bubble_timeperiod_12stocks.R`](https://github.com/AbirS2026/SVADF-Bubble/blob/main/Bubble_timeperiod_12stocks.R)

This script compares bubble classifications across the two main empirical subperiods considered in the paper. It applies the SV-ADF and PWY procedures to the twelve AI-exposed firms and produces the period-level bubble/no-bubble classifications.

Relevant paper outputs:

- Figure 11 of the paper
- Table 2 of the paper

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

- Table 3
- Table 4
- Table 5

---

## How to Use the Code

Each script can be run independently in R. The empirical scripts download price data from Yahoo Finance, usually through the `quantmod` package.

A typical workflow is:

```r
source("Recursive SV-ADF Date stamps Tech.R")
