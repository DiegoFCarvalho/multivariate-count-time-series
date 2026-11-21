# Jittering test

Compares the exact discrete likelihood with jittered data when fitting a vine copula to Poisson-margin count data, using the mean squared error (MSE) of the estimated copula parameter as error measure. The script needs no input data and takes about 2 minutes 49 seconds (4 workers).

## Run

    cd jittering_test
    Rscript jittering_sim.R

The script finds its own directory and writes next to itself, so it can be started from any working directory.

## Output

| File | Thesis |
|---|---|
| `mse_convergence_lambda_combined.png` | Figure 6.1 |
| `mse_summary_table.csv` | Values of Table 6.1 |
| `mse_summary_table.png` | Values of Table 6.1 |
| `tables/mse_per_replicate.csv`, `tables/true_parameters.csv` | Not in the thesis |

With the packages below, all output files are reproduced exactly (byte for byte).

## Packages

R 4.6.0, rvinecopulib 0.7.3.1.0, ggplot2 4.0.3, patchwork 1.3.2, gridExtra 2.3. (rvinecopulib 1.0.0.1.0, as used by the data analysis, was not tested here.)

## Hardcoded settings

| File:line | What | Impact |
|---|---|---|
| `jittering_sim.R:754, 767, 780` | Poisson lambda = 5, 50, 250 | One experiment per lambda |
| `jittering_sim.R:748, 761, 774` | n = 500 | Sample size |
| `jittering_sim.R:749, 762, 775` | Dimension 4 | Vine dimension |
| `jittering_sim.R:750, 763, 776` | D-vine | Vine structure |
| `jittering_sim.R:751-752, 764-765, 777-778` | Clayton copula, parameter 2 | True dependence |
| `jittering_sim.R:755, 768, 781` | 200 replicates | Monte Carlo error |
| `jittering_sim.R:756, 769, 782` | 4 workers | Speed only |
| `jittering_sim.R:757, 770, 783` | `base_seed` = 2028 | Replicate k uses seed `base_seed + k` for the data (line 477) and `base_seed + 100000 + k` for the jitter (line 492) |
| `jittering_sim.R:820`, `:878` | Output sizes and 300 dpi of the table and figure | Appearance |
