# Driver: forecast comparison of svines, MDACP and the Kim et al. (2025) model.
# Run `Rscript run_data_analysis.R [K=<draws>]` from any directory.

.this_dir <- function() {
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a), winslash = "/")) else getwd()
}
.config_dir <- .this_dir()
source(file.path(.config_dir, "config.R"))

cl <- commandArgs(trailingOnly = TRUE)
if (any(grepl("^K=", cl))) K <- as.integer(sub("^K=", "", grep("^K=", cl, value = TRUE)))

for (p in names(required_versions))
  if (as.character(packageVersion(p)) != required_versions[[p]])
    stop(sprintf("%s %s loaded, %s required (see KIM_LIB_DIR)", p, packageVersion(p), required_versions[[p]]))
message("S-vine models: ", paste(svine_main_models, collapse = ", "))
message("root_dir:    ", root_dir)
message("results_dir: ", results_dir)
message("plots_dir:   ", plots_dir)
message(sprintf("K = %d, n_train = %d, skip_if_exists = %s", K, n_train, skip_if_exists))

for (f in c("fit_and_draw.R", "run_models.R", "prep_data.R", "kim_forecast.R", "metrics.R",
            "plot_forecast_errors.R", "plot_data_overview.R", "plot_results_table.R", "diagnostics.R"))
  source(file.path(root_dir, f))

run_if_needed <- function(name, done, fun) {
  if (!run_stage[[name]]) return(invisible(NULL))
  if (skip_if_exists && done()) { message("[", name, "] outputs exist, skipped"); return(invisible(NULL)) }
  t0 <- Sys.time(); message("[", name, "] start")
  fun()
  message(sprintf("[%s] done in %.0fs", name, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
all_groups <- names(groups)
exists_all <- function(model) all(file.exists(sapply(all_groups, raw_file, model = model, K = K)))

run_if_needed("data", function() all(file.exists(file.path(data_dir, c(sprintf("%s_%s.rds", data_prefix, all_groups),
                                                                       sprintf("%s_%s.csv", data_prefix, all_groups))))),
              prep_data)
run_if_needed("svines", function() all(sapply(svine_main_models, exists_all)),
              function() for (m in svine_main_models) run_model_stage(m, all_groups, K))
run_if_needed("mdacp", function() exists_all(mdacp_model), function() run_model_stage(mdacp_model, all_groups, K))
run_if_needed("kim", function() exists_all("kim"), function() run_kim_stage(all_groups, K))
run_if_needed("metrics", function() metrics_done(K), function() run_metrics(K))
run_if_needed("plots", function() plots_done(K), function() run_plots(K))
run_if_needed("diagnostics", diagnostics_done, function() run_diagnostics(K))
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
invisible(lapply(names(required_versions), loadNamespace))   # so that sessionInfo() lists them
writeLines(capture.output(sessionInfo()), file.path(results_dir, "sessionInfo_both.txt"))
