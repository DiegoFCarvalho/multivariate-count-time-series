# Forecast draws of one registry model (svines or MDACP) on every group. Defines functions only.
# Output per group: <results_dir>/<data_prefix>_<group>/raw/{raw,fit}_<model>_K<K>.rds

group_dataset_id <- function(g) sprintf("%s_%s", data_prefix, g)
group_raw_dir    <- function(g) file.path(results_dir, group_dataset_id(g), "raw")
raw_file         <- function(model, g, K) file.path(group_raw_dir(g), sprintf("raw_%s__oos_free_K%d.rds", model, K))
fit_file         <- function(model, g, K) file.path(group_raw_dir(g), sprintf("fit_%s_K%d.rds", model, K))

read_group <- function(g) read_data(file.path(data_dir, sprintf("%s_%s.rds", data_prefix, g)))

run_model_stage <- function(model, groups_sel = names(groups), K) {
  for (g in groups_sel) {
    if (skip_if_exists && file.exists(raw_file(model, g, K))) next
    x <- read_group(g)
    stopifnot(nrow(x) == dimension, ncol(x) == n_days)
    t0  <- Sys.time()
    res <- run_evaluation(x, root_dir, model, n_train, K, n_workers, base_seed, est_timeout,
                          stream_ids = stream_ids[[model]], p_max = p_max_seed)
    f <- res$fits[[model]]
    dir.create(group_raw_dir(g), showWarnings = FALSE, recursive = TRUE)
    saveRDS(res$results[[model]], raw_file(model, g, K))
    saveRDS(list(model = model, status = f$status, stationary = f$stationary, fit_time = f$fit_time,
                 theta_hat = f$theta, K = K, n_train = n_train, base_seed = base_seed,
                 stream_id = stream_ids[[model]], p_max = p_max_seed,
                 versions = c(R = as.character(getRversion()),
                              vapply(c("svines", "rvinecopulib", "kde1d"),
                                     function(p) as.character(packageVersion(p)), character(1)))),
            fit_file(model, g, K))
    message(sprintf("[%s] %s done in %.0fs (fit status %s)", model, g,
                    as.numeric(difftime(Sys.time(), t0, units = "secs")), f$status))
  }
  invisible(NULL)
}
