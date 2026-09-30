# evalbehtransmission
Evaluation of nowcasting and forecasting using behavioural data including CoMix and Google Mobility

## Data Sources

- **Google COVID-19 Community Mobility Reports**: Google LLC. https://www.google.com/covid19/mobility/
- **CoMix social contact data (UK)**: https://zenodo.org/records/13684044
- **inc2prev** (incidence, R(t) and modelled positivity for England, fitted to the ONS COVID-19 Infection Survey): https://github.com/epiforecasts/inc2prev
- **OxCGRT** (policy stringency, used for period definitions only): https://github.com/OxCGRT/covid-policy-dataset

## Repository state (September 2026)

`main` holds the Chapter 1 pipeline. Work in progress sits on branches:

| Branch | PR | What it holds |
|---|---|---|
| `51-inc2prev-generation-interval` | #56 (draft) | Switch between EpiEstim and inc2prev generation intervals. `gi_type` argument still to add (#57) |
| `54-comix-population-weights` | #60 (draft) | CoMix mean contacts weighted to the ONS mid-2020 UK population. Review before rerunning the pipeline |
| `appendix-scripts` | #59 (draft) | Scripts for the upgrading report's Appendix A.1 figures |
| `55-synthetic-outbreak` | — | Synthetic outbreak plan (`synthetic/ch1_synthetic_experiment.qmd`, #55) |

## Chapter 1

Fits a renewal-equation GAM to England infection incidence and tests whether CoMix contacts or Google Mobility improve 1-4 week forecasts. Four models are compared: baseline, contacts, mobility, combined. Negative binomial used with April 2020 to January 2021 study window (pre-mass-vaccination).

Run script in this order, or use `analysis/run_pipeline.R`

| Script | What it does |
|---|---|
| `process_mobility.R` | Google Mobility, UK national rows |
| `process_comix.R` | CoMix contact matrix eigenvalue per survey round, and daily mean contacts (14-day trailing, age-standardised) |
| `ch1_data.R` | Creates a single dataset with date, incidence, contacts, and mobility |
| `ch1_covariates.R` | Processes inputs for modelling: Mobility composite stream with trailing 7-day mean; both covariates z-scored |
| `ch1_periods.R` | Defines periods based on lockdowns, tiered restrictions, relaxation etc. from OxCGRT. Descriptive, not used as covariate |
| `ch1_diagnostics.R` | Testing whether Poisson or NegBin is more suitable based on residuals, autocorrelation, Rt validation against inc2prev |
| `ch1_gam.R` | Functions for running GAM models. Calculates generation interval, creates data frame for each fit, defines renewal formula and fits models |
| `ch1_rolling.R` | Runs the forecasts. 27 weekly origins (2020-07-01 to 2020-12-30), 8-week trailing training window, 4-week horizon, 200 trajectories each |
| `ch1_scoring.R` | Calculates CRPS on log and natural scales, decomposition, bias, interval coverage, by horizon and period |
| `ch1_window_plots.R` | Plot model coefficients and deviance explained across different windows, showing strength and uncertainty in relationships |
| `ch1_forecast_plots.R` | Produces forecast fan plots by period and rolling origin, with a zoomed in version where there's high variation in performance |
| `ch1_descriptive.R` | Produces plots of covariates (z-scored) against Rt |

`R/` contains the remote path for inc2prev and shared choices for which mobility streams or CoMix covariate to retain, passed through subsequent analysis and plot names.

Data goes to `data-processed/`, tables and figures to `outputs/ch1/`

### Other scripts

Not part of the pipeline:

| Script | Status |
|---|---|
| `eda.R` | Early exploration, superseded by `ch1_descriptive.R` |
| `model_rtglm.R`, `model_rtgam.R` | Pre-pipeline GLM/GAM forecasts that `ch1_gam.R` and `ch1_rolling.R` were built from (#48) |
| `process_ons.R` | Reads a local copy of inc2prev estimates, superseded by `ch1_data.R` reading inc2prev remotely |
| `age_incidence.R` | Age-stratified prevalence and incidence plots from inc2prev |
| `publication_tables.R` | Renders score tables from `outputs/ch1/` as Word tables for the report |
| `ch1_ccf.R` | Stub for cross-correlation by period (#12) |

`stashed/` holds earlier experiments (mvgam, forecast comparisons), kept for reference.

## Running

```r
renv::restore()
```

```sh
Rscript analysis/run_pipeline.R
```

Everything is fetched remotely, so a fresh clone needs no local data. Each processing script has a `use_remote` toggle for working from local copies. `process_comix.R` downloads ~150MB and takes a few minutes potentially. After that, everything runs within a couple of minutes.
