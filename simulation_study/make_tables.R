# Simulation tables of the thesis (numerical and GoF statistics per DGP and T).
# Run from simulation_study/. Reads results/T<T>/<dgp>_T<T>.rds.
# Writes results/tables/tab_<dgp>_{num,gof}.csv and .png.

dgps <- c(dgp1 = "minar1_poi", dgp2 = "minar1_negbin", dgp3 = "mdacp1", dgp4 = "svine1", dgp5 = "mcar1")
horizons <- c(100L, 500L)
est_met_order <- c("minar1_poi", "minar1_negbin", "mdacp1", "svine1")

# Evaluates only the aggregation functions of run_simulation.R; sourcing it would run the driver.
helper_names <- c("finite_mean", "compute_gof_table", "compute_convergence_rate", "compute_nonconvergence_rate",
                  "compute_timeout_rate", "compute_stationarity_rate", "compute_fit_time_table")
agg <- new.env()
for (e in parse(file = "run_simulation.R")) {
  if (is.call(e) && identical(e[[1]], as.name("<-")) && is.name(e[[2]]) && as.character(e[[2]]) %in% helper_names) {
    eval(e, agg)
  }
}
stopifnot(all(helper_names %in% ls(agg)))

tab_dir <- file.path("results", "tables")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

numeric_table <- function(df, horizon) {
  mets <- est_met_order
  fit_time <- agg$compute_fit_time_table(df)
  data.frame(
    T = horizon, est_met = mets,
    convergence_rate    = unname(agg$compute_convergence_rate(df)[mets]),
    nonconvergence_rate = unname(agg$compute_nonconvergence_rate(df)[mets]),
    timeout_rate        = unname(agg$compute_timeout_rate(df)[mets]),
    stationarity_rate   = unname(agg$compute_stationarity_rate(df)[mets]),
    fit_time_mean       = unname(fit_time[mets, "mean"]),
    fit_time_sd         = unname(fit_time[mets, "sd"]),
    row.names = NULL)
}

gof_table <- function(df, horizon) {
  tab <- agg$compute_gof_table(df)[est_met_order, , drop = FALSE]
  data.frame(T = horizon, est_met = est_met_order, as.data.frame(tab), row.names = NULL)
}

write_csv <- function(df, file) {
  for (j in which(vapply(df, is.double, logical(1)))) {
    df[[j]] <- ifelse(is.finite(df[[j]]), sprintf("%.17g", df[[j]]), "NA")
  }
  write.csv(df, file.path(tab_dir, file), row.names = FALSE, quote = FALSE)
}

# "-" marks a statistic that no converged fit could provide, as in the thesis tables.
fmt <- function(x, digits) ifelse(is.finite(x), formatC(round(x, digits), format = "f", digits = digits), "-")

write_png <- function(df, file) {
  grob <- gridExtra::tableGrob(df, rows = NULL)
  width <- sum(grid::convertWidth(grob$widths, "in", valueOnly = TRUE)) + 0.3
  height <- sum(grid::convertHeight(grob$heights, "in", valueOnly = TRUE)) + 0.3
  ggplot2::ggsave(file.path(tab_dir, file), plot = grob, width = width, height = height, units = "in", dpi = 300)
}

# Null device for measuring the table grobs, so no Rplots.pdf is written.
grDevices::pdf(NULL)

for (dgp in names(dgps)) {
  runs <- lapply(horizons, function(h) {
    readRDS(file.path("results", paste0("T", h), sprintf("%s_T%d.rds", dgps[[dgp]], h)))
  })

  num <- do.call(rbind, Map(numeric_table, runs, horizons))
  write_csv(num, sprintf("tab_%s_num.csv", dgp))
  write_png(data.frame(
    T = num$T, `Est.Met.` = num$est_met,
    `Convergence rate` = fmt(num$convergence_rate, 3), `Non-convergence rate` = fmt(num$nonconvergence_rate, 3),
    `Time-out rate` = fmt(num$timeout_rate, 3), `Stationarity rate` = fmt(num$stationarity_rate, 3),
    `Mean fit time (s)` = fmt(num$fit_time_mean, 2), `Fit time sd (s)` = fmt(num$fit_time_sd, 2),
    check.names = FALSE), sprintf("tab_%s_num.png", dgp))

  gof <- do.call(rbind, Map(gof_table, runs, horizons))
  write_csv(gof, sprintf("tab_%s_gof.csv", dgp))
  gof_png <- data.frame(T = gof$T, `Est.Met.` = gof$est_met, check.names = FALSE)
  for (col in setdiff(names(gof), c("T", "est_met"))) gof_png[[col]] <- fmt(gof[[col]], 4)
  write_png(gof_png, sprintf("tab_%s_gof.png", dgp))
}

invisible(grDevices::dev.off())
