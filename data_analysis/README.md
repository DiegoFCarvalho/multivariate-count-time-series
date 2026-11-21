# Data analysis: forecast comparison on the Kim et al. (2025) data

Compares five-day-ahead forecasts of four models on daily ratings (1-5) of one student (30 items, 90 days, five groups C1-C5 of six items): S-vines with exact ecdf margins (`svines_exact`), S-vines with scaled ecdf margins (`svines_scaled`), MDACP (`mdacp1`, rank 2) and the latent Gaussian model of Kim et al. (`kim`, r = 1). Each model is fitted on days 1-85; 500 free draws of days 86-90 are scored by the mean absolute error.

## Run

```
cd data_analysis
Rscript run_data_analysis.R          # all stages; optional argument K=<draws>
```

`run_data_analysis.R` runs the stages below in order and skips a stage whose output files already exist (`KIM_FORCE=1` recomputes). It works from any working directory.

| Stage | Code | Output |
|---|---|---|
| `data` | `prep_data.R` | `data/` |
| `svines` | `run_models.R`, `fit_and_draw.R`, `models/model_svines.R` | `results/kim_<group>/raw/` |
| `mdacp` | `run_models.R`, `fit_and_draw.R`, `models/model_mdacp.R` | `results/kim_<group>/raw/` |
| `kim` | `kim_forecast.R` | `results/kim_<group>/raw/` |
| `metrics` | `metrics.R` | `results/abs_err_*.csv` |
| `plots` | `plot_forecast_errors.R`, `plot_data_overview.R`, `plot_results_table.R` | `plots/` |
| `diagnostics` | `diagnostics.R` | `results/draw_support*.csv`, `results/fit_diagnostics*.csv`, two PNG tables in `plots/` |

To reproduce the stored outputs from scratch, run with empty `KIM_RESULTS_DIR`, `KIM_PLOTS_DIR`, `KIM_DATA_DIR` (or delete the old outputs). Repeated runs gave identical outputs apart from fit-time fields.

## Files

| File | Role |
|---|---|
| `run_data_analysis.R` | Driver: checks package versions, sources the other files, runs the stages |
| `config.R` | All settings, paths and environment variables |
| `prep_data.R` | Writes the group data files from the Kim et al. data |
| `run_models.R` | Runs one S-vine or MDACP model on every group and saves draws and fit records |
| `fit_and_draw.R` | Fit and draw engine: seeds, parallel workers, absolute-error helpers |
| `kim_forecast.R` | Fit and draws of the Kim et al. model |
| `metrics.R` | Mean absolute errors, running means, results table (CSV, TEX) |
| `plot_forecast_errors.R` | Convergence and item-day figures; shared plot helpers |
| `plot_data_overview.R` | Data overview figure |
| `plot_results_table.R` | Reference image of the values of Table 7.2 (PNG) |
| `diagnostics.R` | Draw support and fit diagnostics (CSV, PNG tables) |
| `models/model_interface.R` | Model contract and timeout helper |
| `models/model_registry.R` | List of the three fitted models (`svines_exact`, `svines_scaled`, `mdacp1`) |
| `models/model_svines.R` | S-vine model with exact or scaled ecdf margins |
| `models/model_mdacp.R` | MDACP model |
| `kim_repo/Code/*.R` | Code of Kim et al. (third party, three patches below) |

## Data

The data are `kim_repo/Code/Data.rda` of the repository of Kim et al. (https://github.com/yk748/latent-Gaussian-TSM, commit 06b887c, own `LICENSE`; no download needed): daily self-ratings of 30 personality items (Borkenau & Ostendorf, 1998) of one student over 90 days. The original 0-6 scale is merged to 1-5 (0 into 1, 6 into 5); the file holds only the integers 1 to 5 and no missing values (checked by `prep_data()`). Days 1-85 are used for estimation and days 86-90 for evaluating the forecasts.

The 30 items form five groups of six items (`config.R:22-28`). All groups come from the same student, so they are not independent samples.

| Group | Items |
|---|---|
| C1 | lazy, industrious, persistent, reckless, changeable, responsible |
| C2 | dynamic, sociable, shy, silent, lively, reserved |
| C3 | selfish, goodnatured, domineering, helpful, obstinate, considerate |
| C4 | unimaginative, witty, knowledgeable, prudent, fanciless, uninformed |
| C5 | irritable, emotionallystable, calm, badtempered, resistant, vulnerable |

`data/` holds `kim_C1` to `kim_C5` as `.rds` (6 x 90 integer matrix, item names as row names) and as `.csv` (values only, no header). The code reads only the `.rds` files; the `.csv` files are copies for use outside R. To regenerate `data/`, run `run_data_analysis.R` with `KIM_STAGES=data KIM_FORCE=1`.

## Output

The thesis tables are typed by hand in LaTeX.

| Output | Thesis |
|---|---|
| `plots/data_overview.png` | Figure 7.1 |
| `plots/abs_err_convergence_C1_K500.png` ... `C5` | Figures 7.2, 7.4, 7.6, 7.8 and 7.10 |
| `plots/abs_err_item_day_C1.png` ... `C5` | Figures 7.3, 7.5, 7.7, 7.9 and 7.11 |
| `results/abs_err_results_table.csv` | Values of Table 7.2 |
| `plots/abs_err_results_table.png` | Reference image of the values of Table 7.2 (the thesis table is typed by hand) |
| `results/draw_support_table.csv` | Values of Table 7.1 (the thesis shows only the MDACP and scaled-ecdf columns; "Total" = "Pooled") |
| `results/draw_support.csv` | Values quoted in the text (pooled shares, maximum draw) |
| `results/fit_diagnostics_kim.csv` | Values of the two negative variances quoted in the text |

All other files in `results/` and `plots/` are intermediate or alternative outputs not used in the thesis. `results/abs_err_results_table.tex` has a different column order than the thesis table.

Raw draws and fit records are in `results/kim_<group>/raw/` (`raw_<model>__oos_free_K500.rds`, `fit_<model>_K500.rds`). `results/sessionInfo_both.txt` is written by the driver.

Stored results are from the full run (K = 500, all groups). `results/run_both.log` and `run_plots.log` are logs of that run; they show the old folder name `data_analysis_merged` in absolute paths (stored result, not edited).

## Runtime

Full run (K = 500, five groups, 4 workers; stage times from `results/run_both.log`): `svines` 1148 s, `mdacp` 1447 s, `kim` 100 s, `metrics` 1 s, `plots` 27 s, `diagnostics` 5 s; about 46 minutes. The `mdacp` time varies between runs.

## Environment

- R 4.6.0, Windows 11.
- The driver stops unless svines 0.3.0, rvinecopulib 1.0.0.1.0 and kde1d 1.1.1 are loaded (`config.R:13`).
- The default R library of the machine has svines 0.2.7 and rvinecopulib 0.7.3.1.0, so the required versions are installed into a separate folder. The copies in `C:/Work/r_lib030` come from CRAN (`Repository: CRAN` in their DESCRIPTION).
- The install commands below were not tested on a clean machine.

```r
lib <- "C:/Work/r_lib030"; dir.create(lib)       # any folder
.libPaths(c(lib, .libPaths()))
# kde1d 1.1.1 is an archive version (CRAN now has 1.2.2): source install, needs Rtools on Windows
install.packages("https://cran.r-project.org/src/contrib/Archive/kde1d/kde1d_1.1.1.tar.gz",
                 repos = NULL, type = "source", lib = lib)
# svines 0.3.0 and rvinecopulib 1.0.0.1.0 are the current CRAN versions (Windows binaries exist)
install.packages(c("rvinecopulib", "svines"), lib = lib)
```

- If CRAN has newer versions, take the exact versions from `https://cran.r-project.org/src/contrib/Archive/<package>/` with `type = "source"`.
- Source installs need Rtools matching the R version (Windows) or a C++17 compiler (Linux/macOS).
- Check with `packageVersion("svines", lib.loc = lib)`.
- Point the driver at the folder with `KIM_LIB_DIR=<folder>` (default `C:/Work/r_lib030`, `config.R:11`). The folder is put first in the library paths and exported as `R_LIBS` to the workers.
- The MDACP and Kim code do not use svines or rvinecopulib.

| Package | Version |
|---|---|
| svines | 0.3.0 |
| rvinecopulib | 1.0.0.1.0 |
| kde1d | 1.1.1 |
| univariateML (dependency of svines) | 1.5.0 |
| gamlss.dist (MDACP) | 6.1-1 |
| mvtnorm | 1.3-7 |
| Matrix | 1.7-5 |
| polynom | 1.4-1 |
| TruncatedNormal | 2.3 |
| callr | 3.7.6 |
| future, future.apply | 1.75.0, 1.20.2 |
| withr | 3.0.2 |
| ggplot2, patchwork | 4.0.3, 1.3.2 |

`polynom` and `TruncatedNormal` are used only by the Kim code (`polynom`: Hermite polynomials in `Latent_Gauss_Link.R`; `TruncatedNormal`: `rtmvnorm` and `pmvnorm` in `Latent_DFM_Forecast.R`). Nothing installs them automatically, so install them with `install.packages()`.

The simulation study ran on a VM with R 4.4.1, the same svines and rvinecopulib versions, and other versions of the remaining packages (see `simulation_study/README.md`).

## Environment variables and arguments

| Name | Default | Effect |
|---|---|---|
| `KIM_ROOT_DIR` | directory of `config.R` | Code directory |
| `KIM_OUT_DIR` | `KIM_ROOT_DIR` | Base for `results/` and `plots/` |
| `KIM_DATA_DIR` | `<root>/data` | Data folder |
| `KIM_REPO_DIR` | `<root>/kim_repo` | Kim et al. copy |
| `KIM_RESULTS_DIR` | `<out>/results` | Results folder |
| `KIM_PLOTS_DIR` | `<out>/plots` | Plots folder |
| `KIM_LIB_DIR` | `C:/Work/r_lib030` | Library with svines 0.3.0 and rvinecopulib 1.0.0.1.0; empty string = none |
| `KIM_GROUPS` | all five groups | Subset, e.g. `C1,C2` |
| `KIM_WORKERS` | 4 | Worker processes |
| `KIM_STAGES` | all seven | Stages to run, e.g. `mdacp,kim` |
| `KIM_FORCE` | unset | `1` recomputes stages whose outputs exist |
| `K=<n>` (argument) | 500 | Number of draws |

## Hardcoded settings

| File:line | What | Impact |
|---|---|---|
| `config.R:11` | Library path `C:/Work/r_lib030` | Must exist or be overridden with `KIM_LIB_DIR`, else the version check stops the run |
| `config.R:13` | Required package versions | Driver stops on any other version |
| `config.R:15-20` | 90 days, 85 training days, 6 items, file prefix `kim`, 5 categories | Data split and dimensions |
| `config.R:22-28` | Item names of groups C1-C5 | Group composition |
| `config.R:31` | K = 500 draws | Monte Carlo error; `K=` argument overrides |
| `config.R:33` | `base_seed` = 1 | All random streams |
| `config.R:34` | Fit timeout 600 s | A longer fit is recorded as `timed_out` |
| `config.R:39-40` | Kim `r = 1`, VAR order 1 | Kim model |
| `config.R:41-42` | Stream per model (S-vines share stream 1, MDACP 2, Kim 3), `p_max_seed` = 2 | Seed layout; changing it changes all draws |
| `diagnostics.R:69` | Tolerance 1e-8 for the recomputed Kim fit | Run stops if exceeded |
| `diagnostics.R:174` | Column widths and height of the draw-support table PNG | Appearance only |
| `diagnostics.R:181` | Column widths and height of the fit-diagnostics table PNG | Appearance only |
| `kim_forecast.R:55` | `identy_opt = 2` | Kim identification |
| `kim_repo/Code/Latent_DFM_Forecast.R:22-23` | 5-day window, 100 particles | Kim forecast |
| `metrics.R:14-18` | The five plotted series (model, variant, label) | Table and figure columns |
| `metrics.R:21` | Clipping rule `pmin(pmax(round(x), 1), 5)` | Clipped errors |
| `metrics.R:98` | Tolerance 1e-10 for the item-day check | Run stops if exceeded |
| `models/model_interface.R:5` | Default timeout 600 s (used if the option is unset) | Same as `config.R:34` |
| `models/model_mdacp.R:118` | `maxit = 1000` | MDACP estimate |
| `models/model_mdacp.R:186` | Lower bound 1e-3 for phi | MDACP estimate |
| `models/model_mdacp.R:192-199` | Start values | MDACP estimate |
| `models/model_registry.R:6-8` | S-vines `p = 2`, MDACP `r = 2` | Model orders |
| `plot_data_overview.R:36` | Figure size 12 x 8 in | Appearance only |
| `plot_results_table.R:34` | Figure size 10.5 x 3.1 in | Appearance only |
| `plot_forecast_errors.R:8-12` | Colours, line types, point shapes of the five series | Appearance only |
| `plot_forecast_errors.R:47-48` | 300 dpi | Appearance only |
| `plot_forecast_errors.R:83` | Convergence figure size 10 x 5.2 in | Appearance only |
| `plot_forecast_errors.R:125` | Item-day figure size 12 x 4.6 in | Appearance only |
| `run_data_analysis.R:49` | `results/sessionInfo_both.txt` | Output name |
| `fit_and_draw.R:108` | One BLAS thread per worker | Speed only |

## Models directory

`models/` is separate from `simulation_study/models/` on purpose:

| Aspect | `data_analysis/models` | `simulation_study/models` |
|---|---|---|
| Models | `svines_exact`, `svines_scaled`, `mdacp1` | Also MINAR and the count autoregression |
| Timeout | Option `est_timeout_seconds`, default 600 s | `est_timeout_seconds` from `get_config()` |
| `true_moments` slot | None | Present |
| `new_model()` checks | `stopifnot` on `sim`, `p`, `est` | Same checks with error messages, plus `true_moments` |
| Model class | `forecast_model` | `sim_study_model` |

## Model implementation

- S-vines (Markov order 2, all variables discrete). `svines_exact`: `svinecop` + `svine_dist` + `svine_sim`; margins `cf = cumsum(counts)/85` of days 1-85 with `cf[5] = 1`, quantile `1 + findInterval(u, cf, left.open = TRUE)` (as in `Latent_Gauss_Link.R`). `svines_scaled`: `svines::svine(..., margin_families = "empirical")` (cdf times n/(n+1), interpolating quantile); its draws are rounded when scored.
- MDACP: rank 2, double Poisson margins, Gaussian copula. The `optim` fit (L-BFGS-B) runs in a `callr` child; random draws are made in the parent.
- Kim: r = 1 factor, VAR order 1, 5 categories, 100 particles, filter over the last 5 days.
- Draws: every draw has its own random-number stream (`base_seed` 1), so results do not depend on K, the worker count or the other models.
- Clipping (post hoc): no-op for Kim and `svines_exact`, rounding for `svines_scaled`, rounding and clipping to 1..5 for MDACP.

## Patches to the Kim code

- `Latent_DFM_Estim.R`: the dimension d is taken from `Cov_Z`.
- `Latent_DFM_Forecast.R`: array dimensions are kept for r = 1 (`drop = FALSE`, `matrix()`).
- `Latent_DFM_Forecast.R`: covariance matrices are symmetrised before `rtmvnorm`.
