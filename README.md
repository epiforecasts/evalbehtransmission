# evalbehtransmission
Evaluation of nowcasting and forecasting using behavioural data including CoMix and Google Mobility

## Data Sources

- **Google COVID-19 Community Mobility Reports**: Google LLC. https://www.google.com/covid19/mobility/
- **CoMix social contact data (UK)**: https://zenodo.org/records/13684044

## Chapter 1

Fits a renewal-equation GAM to England infection incidence and tests whether CoMix contacts or Google Mobility improve 1-4 week forecasts. Four models are compared: baseline, contacts, mobility, combined. Negative binomial used with April 2020 to January 2021 study window (pre-vaccine).

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

## Running

renv::restore()
Rscript analysis/run_pipeline.R

Everything is fetched remotely, so a fresh close needs no local data. Each processing script has a `use_remote` toggle for working from local copies. `process_comix.R` downloads ~150MB and takes a few minutes potentially. After that, everything runs within a couple minutes.
