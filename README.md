# Multivariate Count Time Series: thesis code

R code for the figures and numerical results of the thesis. Each folder has its
own README with run instructions, output, packages and hardcoded settings;
start there. The folders are independent of each other.

| Folder | Content |
|---|---|
| `chapter_2/` | Figures of the count distributions, Chapter 2 |
| `jittering_test/` | Simulation comparing the exact discrete likelihood with jittering, Section 6.1 |
| `simulation_study/` | Five data-generating processes, four estimation methods, Section 6.2 |
| `data_analysis/` | Forecast comparison of S-vines, MDACP and the Kim et al. model on personality rating data, Chapter 7 |

`master_thesis.Rproj` is the RStudio project file.

## Reproducibility

All outputs can be reproduced by running the scripts in each folder. The one
practical exception is the estimation of the simulation study. It was run on a
Linux VM (20 workers) because it needs more compute than the local machine
has. The scripts can also be run locally, but the stored results in
`simulation_study/results/` come from the VM runs and were copied over.

Each VM run had a time limit, so the T = 500 case was run twice: first DGP 1,
2, 3 and 5, then DGP 4 (`svine1`) alone, each run with all four methods
(details in `simulation_study/README.md`).

The tables and true moments of the simulation study are recomputed locally
from the stored results.

## Software

Local runs used Windows 11 and R 4.6.0, the VM R 4.4.1.

| Part | svines | rvinecopulib |
|---|---|---|
| Simulation study, estimation (VM) | 0.3.0 | 1.0.0.1.0 |
| Simulation study, `make_tables.R` | any | any |
| Simulation study, `true_moments.R` | 0.3.0 | 1.0.0.1.0 |
| Data analysis | 0.3.0 | 1.0.0.1.0 |
| Jittering test | not used | 0.7.3.1.0 |
| Chapter 2 | not used | not used |

The package versions used to produce the stored output of `true_moments.R` are
not recorded. The data analysis also needs kde1d 1.1.1 and loads these
packages from a separate library (`KIM_LIB_DIR`, see `data_analysis/README.md`).
The default R library of the local machine has svines 0.2.7 and rvinecopulib
0.7.3.1.0. The other VM package versions are listed in
`simulation_study/README.md`.

## Thesis figures and tables

Figures are the PNG files below. Tables are typed by hand in the thesis; the
listed files contain their values.

| Thesis | Script | File |
|---|---|---|
| Figures 2.1 to 2.5 | `chapter_2/distribution_plots.R` | `chapter_2/chapter2_plots/{binomial,poisson,geom,negbin1,negbin2}_pmf_cdf.png` |
| Figure 6.1 | `jittering_test/jittering_sim.R` | `jittering_test/mse_convergence_lambda_combined.png` |
| Table 6.1 | `jittering_test/jittering_sim.R` | `jittering_test/mse_summary_table.csv` |
| Tables 6.2 to 6.11 | `simulation_study/make_tables.R` | `simulation_study/results/tables/tab_dgp{1..5}_{num,gof}.csv` |
| Figure 7.1 | `data_analysis/plot_data_overview.R` | `data_analysis/plots/data_overview.png` |
| Table 7.1 | `data_analysis/diagnostics.R` | `data_analysis/results/draw_support_table.csv` |
| Table 7.2 | `data_analysis/metrics.R` | `data_analysis/results/abs_err_results_table.csv` |
| Figures 7.2, 7.4, 7.6, 7.8, 7.10 | `data_analysis/plot_forecast_errors.R` | `data_analysis/plots/abs_err_convergence_C{1..5}_K500.png` |
| Figures 7.3, 7.5, 7.7, 7.9, 7.11 | `data_analysis/plot_forecast_errors.R` | `data_analysis/plots/abs_err_item_day_C{1..5}.png` |

## Known deviations

- The true moments of DGP 5 are Monte Carlo values that depend on the random
  number state of the VM run, which cannot be recreated. Recomputed locally,
  9 of 274 compared numbers differ from the thesis, by at most 0.027
  (details in `simulation_study/README.md`).
- `data_analysis/models` and `simulation_study/models` are separate copies on
  purpose. The differences are listed in `data_analysis/README.md`.

## Third-party code

`data_analysis/kim_repo` contains the code and data of Kim et al. (source
commit 06b887c) with its own `LICENSE`. Three small patches to two files are
listed in `data_analysis/README.md`.