# Chapter 2 figures

PMF and CDF panels of the count distributions in Chapter 2 of the thesis. The
script needs no input data and takes about 6 seconds.

## Run

    cd chapter_2
    Rscript distribution_plots.R

Run from `chapter_2/`. The script writes to `chapter2_plots/` relative to the
working directory and creates the folder if it is missing.

## Output

| File in `chapter2_plots/` | Thesis |
|---|---|
| `binomial_pmf_cdf.png` | Figure 2.1 |
| `poisson_pmf_cdf.png` | Figure 2.2 |
| `geom_pmf_cdf.png` | Figure 2.3 |
| `negbin1_pmf_cdf.png` | Figure 2.4 |
| `negbin2_pmf_cdf.png` | Figure 2.5 |

With the packages below, the five PNG files are reproduced exactly.

## Packages

R 4.6.0, ggplot2 4.0.3, patchwork 1.3.2.

## Hardcoded settings

| File:line | What | Impact |
|---|---|---|
| `distribution_plots.R:6-10` | Three colours | Appearance |
| `distribution_plots.R:151-152` | Binomial: 40 trials, p = 0.25, 0.5, 0.75 | Figure 2.1 |
| `distribution_plots.R:172` | Poisson: lambda = 5, 15, 25 | Figure 2.2 |
| `distribution_plots.R:193` | Geometric: p = 0.3, 0.5, 0.7 | Figure 2.3 |
| `distribution_plots.R:214-215` | Negative binomial: r = 10, p = 0.3, 0.5, 0.7 | Figure 2.4 |
| `distribution_plots.R:236-237` | Negative binomial: p = 0.5, r = 10, 20, 30 | Figure 2.5 |
| `distribution_plots.R:175, 196, 218, 240` | Support truncated at the 0.999 quantile (`qpois`, `qgeom`, `qnbinom`) | Plotted range of infinite supports |
| `distribution_plots.R:79, 86, 104, 132` | Zoom window mean +/- 4 SD | Plotted range |
| `distribution_plots.R:141-143` | Output folder `chapter2_plots`, 12 x 5.5 in, 300 dpi | Output location and size |
