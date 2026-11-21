# Settings for the forecast comparison; paths relative to root_dir, overridable by KIM_* variables.

root_dir     <- normalizePath(Sys.getenv("KIM_ROOT_DIR", unset = if (exists(".config_dir")) .config_dir else getwd()),
                              winslash = "/", mustWork = TRUE)
out_dir      <- Sys.getenv("KIM_OUT_DIR", unset = root_dir)
data_dir     <- Sys.getenv("KIM_DATA_DIR", unset = file.path(root_dir, "data"))
kim_repo_dir <- Sys.getenv("KIM_REPO_DIR", unset = file.path(root_dir, "kim_repo"))
results_dir  <- Sys.getenv("KIM_RESULTS_DIR", unset = file.path(out_dir, "results"))
plots_dir    <- Sys.getenv("KIM_PLOTS_DIR", unset = file.path(out_dir, "plots"))
# Library with svines 0.3.0 and rvinecopulib 1.0.0.1.0; exported as R_LIBS for the workers.
lib_extra    <- Sys.getenv("KIM_LIB_DIR", unset = "C:/Work/r_lib030")
if (nzchar(lib_extra)) { .libPaths(c(lib_extra, .libPaths())); Sys.setenv(R_LIBS = lib_extra) }
required_versions <- c(svines = "0.3.0", rvinecopulib = "1.0.0.1.0", kde1d = "1.1.1")

n_days      <- 90L      # days per series
n_train     <- 85L      # days used for estimation
horizon     <- n_days - n_train
dimension   <- 6L       # items per group
data_prefix <- "kim"    # data files are <data_prefix>_<group>.rds / .csv
n_cat       <- 5L       # rating levels 1..5

groups <- list(
  C1 = c("lazy", "industrious", "persistent", "reckless", "changeable", "responsible"),
  C2 = c("dynamic", "sociable", "shy", "silent", "lively", "reserved"),
  C3 = c("selfish", "goodnatured", "domineering", "helpful", "obstinate", "considerate"),
  C4 = c("unimaginative", "witty", "knowledgeable", "prudent", "fanciless", "uninformed"),
  C5 = c("irritable", "emotionallystable", "calm", "badtempered", "resistant", "vulnerable")
)
if (nzchar(Sys.getenv("KIM_GROUPS"))) groups <- groups[strsplit(Sys.getenv("KIM_GROUPS"), ",")[[1]]]

K           <- 500L                                     # forecast draws per model and group
n_workers   <- as.integer(Sys.getenv("KIM_WORKERS", unset = "4"))
base_seed   <- 1L
est_timeout <- 600                                      # seconds per model fit
# S-vines (Markov order 2): svines_exact = exact ecdf margins, svines_scaled = svines' built-in empirical margins
# (cdf times n/(n+1), interpolating quantile).
svine_main_models <- c("svines_exact", "svines_scaled")
mdacp_model <- "mdacp1"                                 # MDACP, rank 2
kim_r       <- 1L                                       # Kim model: factors
kim_p       <- 1L                                       # Kim model: VAR order
stream_ids  <- c(svines_exact = 1L, svines_scaled = 1L, mdacp1 = 2L, kim = 3L)    # random-number stream per model (both S-vine models share one)
p_max_seed  <- 2L                                       # history length of the seed layout

# A stage whose output files exist is skipped; KIM_FORCE=1 recomputes.
# KIM_STAGES (e.g. "mdacp,kim") selects stages.
run_stage <- c(data = TRUE, svines = TRUE, mdacp = TRUE, kim = TRUE, metrics = TRUE, plots = TRUE, diagnostics = TRUE)
if (nzchar(Sys.getenv("KIM_STAGES")))
  run_stage[] <- names(run_stage) %in% strsplit(Sys.getenv("KIM_STAGES"), ",")[[1]]
skip_if_exists <- !identical(Sys.getenv("KIM_FORCE"), "1")
