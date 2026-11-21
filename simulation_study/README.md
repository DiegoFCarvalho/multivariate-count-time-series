# Simulation study

Compares four estimation methods (`minar1_poi`, `minar1_negbin`, `mdacp1`, `svine1`) for multivariate count time series on five data-generating processes (DGPs). For each DGP, replication and method it fits the model, simulates forecast paths and records moment-based goodness-of-fit (GoF) statistics. There are no input data; all parameters are in the scripts.

## Where the results come from

The estimation (`run_simulation.R`) was run on a Linux VM (20 workers, R 4.4.1) because it needs more compute than the local machine has. It can also be run locally, but the stored results in `results/` come from the VM and were copied over. A local run writes `results/<dgp>.rds`, `results/<dgp>_gof_table.png`, `results/<dgp>_diagnostics_table.png` and `results/run.log` directly into `results/` (`run_simulation.R:54-55, 389, 418, 567`). The stored files are in `results/T100/`, `results/T500/` and `results/tables/`, so none of them is overwritten.

The VM runs had a time limit, so the T = 500 case was run twice. The first run (`run_log_T500.txt`, `run_T500.log`) covered DGP 1, 2, 3 and 5 (`minar1_poi`, `minar1_negbin`, `mdacp1`, `mcar1`); its log ends where DGP 4 (`svine1`) starts. The second run (the files with `svine1_rerun` in their names) covered DGP 4 (`svine1`) only. Both runs fitted all four methods, so the method `svine1` is also in the first run. `results/T500/svine1_T500.rds` comes from the second run.

## Run

All scripts expect `simulation_study/` as working directory.

```
cd simulation_study
Rscript install_packages.R   # installs and checks packages, used to set up the VM
Rscript run_simulation.R     # estimation, writes results/<dgp>.rds
Rscript true_moments.R       # true moments, parameter and stationarity tables
Rscript make_tables.R        # tables from the stored results
```

`true_moments.R` takes a few minutes, `make_tables.R` about 10 seconds. `run_simulation.R` wrote `results/<dgp>.rds`; the files were renamed by hand to `results/T100/<dgp>_T100.rds` and `results/T500/<dgp>_T500.rds` (T = 500 runs: `T = 500` in `config.R`, all else unchanged).

| Script | Role |
|---|---|
| `config.R` | Parameters of a run |
| `gof.R` | Sample mean, autocovariance and rank-correlation statistics of one path |
| `models/model_interface.R` | Model contract (`sim`, `est`, `true_moments`) |
| `models/model_minar.R` | MINAR model (DGPs 1-2, methods `minar1_poi`, `minar1_negbin`) |
| `models/model_mdacp.R` | MDACP with Gaussian copula (DGP 3, method `mdacp1`) |
| `models/model_svines.R` | S-vine model (DGP 4, method `svine1`) |
| `models/model_multiv_count_ar.R` | Count autoregression with copula (DGP 5, simulation only) |
| `models/model_registry.R` | Sources the model files, builds the model list |
| `run_simulation.R` | Driver; also defines the aggregation functions and true parameters |
| `true_moments.R` | True means and autocovariances, parameter table, stationarity checks |
| `make_tables.R` | Tables from the stored results (reads `run_simulation.R` by parsing, does not run it) |
| `install_packages.R` | Installs and checks packages on the VM |

## DGPs

| DGP | Name | Model | Thesis tables |
|---|---|---|---|
| DGP 1 | `minar1_poi` | MINAR(1), Poisson family | Tables 6.2 and 6.3 |
| DGP 2 | `minar1_negbin` | MINAR(1), negative binomial family | Tables 6.4 and 6.5 |
| DGP 3 | `mdacp1` | MDACP(1), rank 2 | Tables 6.6 and 6.7 |
| DGP 4 | `svine1` | S-vine, Markov order 2 | Tables 6.8 and 6.9 |
| DGP 5 | `mcar1` | Count autoregression with copula (simulation only) | Tables 6.10 and 6.11 |

Every DGP except `mcar1` is also fitted as a method.

## Output

The thesis tables are typed by hand in LaTeX.

| Output | Script | Thesis |
|---|---|---|
| `results/T<T>/<dgp>_T<T>.rds` | `run_simulation.R` | Source of all tables below. Data frame with columns `i`, `est_met_name`, `statistic`, `lag`, `value`, `fit_status`, `stationary`, `fit_time`; the configuration is not stored in it |
| `results/T<T>/<dgp>_gof_table_T<T>.png`, `<dgp>_diagnostics_table_T<T>.png`, logs | `run_simulation.R` | Not in the thesis |
| `results/tables/tab_dgp{1..5}_num.csv` | `make_tables.R` | Values of Tables 6.2, 6.4, 6.6, 6.8 and 6.10 |
| `results/tables/tab_dgp{1..5}_gof.csv` | `make_tables.R` | Values of Tables 6.3, 6.5, 6.7, 6.9 and 6.11 |
| `results/tables/tab_dgp{1..5}_{num,gof}.png` | `make_tables.R` | Images of the same tables, not in the thesis |
| `results/tables/true_moments.csv`, `.png`, `dgp_parameters.csv`, `stationarity_checks.csv`, `.png` | `true_moments.R` | No thesis figure; reference values of the DGPs |

Log files: `results/T100/run_log_T100.txt`, `run_T100.log`; `results/T500/run_log_T500.txt`, `run_T500.log`; `results/T500/run_log_T500_svine1_rerun.txt`, `run_T500_svine1_rerun.log` (second T = 500 run).

## Environment

All three runs used R 4.4.1, svines 0.3.0, rvinecopulib 1.0.0.1.0, gridExtra 2.3.1, ggplot2 3.5.1, callr 3.7.6, mvtnorm 1.2.5, Matrix 1.7.0, gamlss.dist 6.1.11, univariateML 1.5.0 and copula 1.1.7. They differ only in future and future.apply: the T = 100 run and the first T = 500 run used future 1.76.0 and future.apply 1.20.2, the second (later) T = 500 run used future 1.80.0 and future.apply 1.21.0.

`install_packages.R` only enforces minimum versions (R >= 4.4.1, rvinecopulib >= 1.0.0.1.0, svines >= 0.3.0) and installs the other packages at their newest CRAN version, so the versions above are taken from the run logs.

Local scripts and vine package versions:
- `make_tables.R`: uses only base R, ggplot2, grid, gridExtra, grDevices. Vine package versions do not matter: run with svines 0.2.7 / rvinecopulib 0.7.3.1.0 and with svines 0.3.0 / rvinecopulib 1.0.0.1.0, it produced byte-identical CSV files, identical to the stored ones.
- `true_moments.R`: use svines 0.3.0 and rvinecopulib 1.0.0.1.0 (the versions of the VM run). The package versions used to produce the stored true-moment output are not recorded.

## Hardcoded settings

| File:line | What | Impact |
|---|---|---|
| `config.R:5-12`, `:14` | `n` = 500 replications, `d` = 3, `T` = 100 (500 for the T = 500 runs), `T_burnin` = 200, `K` = 100, lags 0:3, `seed` = 1, `n_workers` = 20, `output_dir` = `results/` | Whole run. Set `T` by hand for the T = 500 runs; rename the output files by hand |
| `config.R:13` | `est_timeout_seconds` = 480 (matches the thesis; the data analysis uses 600) | A longer fit is recorded as `timed_out` |
| `config.R:16-17` | DGP and method names | Which models run |
| `run_simulation.R:428` | Workers from `config$n_workers` | Parallelism |
| `run_simulation.R:430` | `set.seed(config$seed)` | Random-number state of everything after it |
| `run_simulation.R:435-537` | True parameters of the five DGPs | DGP definitions |
| `run_simulation.R:491-495` | Three `rpois(20, .)` calls for throwaway margin data | Consume random numbers; shift the stream of everything after them (see Known deviations) |
| `true_moments.R:80-86` | 1000 paths, `T` = 500, burn-in 200, true means (5, 7, 9) for the S-vine DGP | Monte Carlo ground truth for DGP 4 |
| `true_moments.R:93` | `set.seed(2026L)` | Seed of the DGP 4 Monte Carlo estimate |
| `make_tables.R:5-7` | DGP names, horizons 100 and 500, method order | Which result files are read |
| `make_tables.R:64` | Paths `results/T<T>/<dgp>_T<T>.rds` | Result file names |
| `install_packages.R` (top) | Minimum versions R 4.4.1, rvinecopulib 1.0.0.1.0, svines 0.3.0 | VM provisioning only |

`make_tables.R` and `true_moments.R` read `run_simulation.R` by parsing it. Keep the seven aggregation helpers as top-level `name <- function` and the `est_mets <-` assignment in place when editing it.

## Known deviations from the thesis

The off-diagonal true moments of DGP 5 (`mcar1`) are Monte Carlo estimates without their own seed; they depend on the state of the random-number stream when `true_moments()` is called. That state differs between `true_moments.R` and the VM run, where the 500 paths of the DGPs processed before `mcar1` are simulated first, and the stream state of the VM run cannot be recovered. Regenerated values differ from the thesis in 9 of 274 compared numbers, by up to 0.027. In a diagnostic run over 20 test seeds the entries varied by 0.002 to 0.024.
