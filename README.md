# Combined_window

# Demographic Forecasting (Mortality & ASFR): Rolling, Expanding, and Combined Windows

This repository provides **point** and **interval** forecasting functions for age–specific **mortality** and **fertility** (ASFR), together with cleaned example datasets for 23 countries. The code is organised so you can drop in any dataset of class `demogdata` and re–use the same forecasting workflow.

---

## Repository layout

- **`Point_forecast_function.R`** — Functions for **point forecasts**. The script ends with a runnable example using **AUS Female** mortality from `demog_list`.
- **`Interval forecast via SD method_function.R`** — Functions for **interval forecasts** via a standard deviation–based method. The script ends with a runnable example using **AUS Female** mortality from `demog_list`.
- **`demog_list.rds`** — Mortality dataset: for each of **23 countries**, male and female age schedules from **0 to 95+**, with both **raw** and **smoothed** versions. Each series is a `demogdata` object that has already been cleaned using the **demograph** workflow.
- **`demog_list_asfr.rds`** — Age–specific fertility rates (ASFR): for each country, female ages **15–49+**, also with **raw** and **smoothed** versions. Each series is a `demogdata` object.

> **Country keys** follow the HMD/HFD codes (e.g., `AUS`, `CAN`, `JPN`).

---

## Data description

- **Mortality (`demog_list`)**: a named list where each entry is a country code (e.g., `AUS`). Within each country you will find male and female `demogdata` objects, each available in **raw** and **smoothed** form. Ages span **0–95+**.
- **Fertility (`demog_list_asfr`)**: a named list with `demogdata` objects for **ASFR** (female ages **15–49+**), again provided in **raw** and **smoothed** form.
- The code is **dataset–agnostic**: you may substitute your own `demogdata` objects (from HMD/HFD or elsewhere). As long as the input is of class `demogdata`, the functions will operate in the same way.

---

## Quick start

1) **Install/Load dependencies** (adjust to your environment):
```r
# Example: install if missing
# install.packages(c("ftsa", "demography"))  # or 'demograph' if your setup uses that package name

library(ftsa)
# library(demography)  # or library(demograph)
```

2) **Load the example datasets**:
```r
demog_list        <- readRDS("demog_list.rds")
demog_list_asfr   <- readRDS("demog_list_asfr.rds")
```

3) **Source the forecasting functions** (both scripts include runnable examples at the bottom):
```r
source("Point_forecast_function.R")
source("Interval forecast via SD method_function.R")
```

4) **Run the built–in example (AUS Female mortality)**  
Open each script and execute the **example section at the end**. As a reminder, AUS Female mortality objects typically appear like:
```r
# Access pattern (illustrative; exact structure depends on your list wrapper)
# Smoothed female rates for Australia:
# demog_list$AUS$Female_smooth$rate[[1]]
```

---

## Accuracy metrics

- **Point forecasts:** MAFE/RMSFE (or your project’s chosen metrics).
- **Interval forecasts:** empirical coverage probability (ECP), conditional coverage deviation (CPD), and the **interval score**. These are reported across horizons $h=1{:}H$.
  
> The example scripts reproduce the AUS Female mortality experiments and show how accuracy is computed.

---

## Reproducibility

- All random seeds (when present) are set within the scripts.
- You may change the window length $W$ and forecast horizon $H$ in the script configuration blocks.

---

## Using your own data

As long as your inputs are `demogdata` objects with the expected components (e.g., `rate`, `pop`, `age`, `year`), you can replace `demog_list.rds`/`demog_list_asfr.rds` with your own files and re–run the same scripts without changing function internals.

---


## Contact / Citation

If you use this code or datasets, please cite the repository and acknowledge the source mortality/fertility databases (HMD/HFD). For questions or issues, open a GitHub issue in this repository.
