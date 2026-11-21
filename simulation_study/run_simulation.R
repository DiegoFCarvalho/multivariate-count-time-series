# Simulation driver: GoF computation per DGP (compute_gofs), then a loop over the DGPs.
# Run from simulation_study/.

source("config.R")
source("models/model_registry.R")
source("gof.R")

library(gridExtra)
library(ggplot2)
library(future)
library(future.apply)


# Multisession workers are fresh R sessions. Static global detection cannot find model helpers
# reached only through a stored closure, so each worker sources the setup files once, from the
# main session's project root.
project_root <- getwd()
worker_setup <- function() {
  if (!exists("all_models", envir = .GlobalEnv, inherits = FALSE)) {
    setwd(project_root)
    source("config.R")
    source("models/model_registry.R")
    source("gof.R")
  }
}


# Elementwise mean over K replicates of a list of named lists of d x d matrices (one per lag).
# NA entries are ignored, so a degenerate replicate does not poison the pooled matrix; an entry
# that is NA in every replicate stays NA.
average_gamma <- function(gamma_list) {
  lag_names <- names(gamma_list[[1]])
  out <- lapply(lag_names, function(h) {
    mats   <- lapply(gamma_list, `[[`, h)
    sums   <- Reduce(`+`, lapply(mats, function(m) replace(m, is.na(m), 0)))
    counts <- Reduce(`+`, lapply(mats, function(m) !is.na(m)))
    sums / replace(counts, counts == 0, NA_real_)
  })
  names(out) <- lag_names
  out
}


# Simulates T + T_burnin observations and drops the first T_burnin, so the recorded path is
# not dominated by the transient from the artificial x_init.
simulate_with_burnin <- function(model, theta, x_init, T, T_burnin) {
  full <- model$sim(theta, x_init, T = T + T_burnin)
  full[, (T_burnin + 1):(T + T_burnin), drop = FALSE]
}


# Appends one timestamped line to run.log in the output directory. Safe to call from workers.
log_line <- function(config, text) {
  cat(sprintf("%s %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text),
      file = file.path(config$output_dir, "run.log"), append = TRUE)
}


# All work for one Est.Met against one DGP path. Returns long-format rows: one per statistic,
# with GoF2 and GoF3 expanded to one row per lag. Called with plain lapply() inside
# compute_gofs()'s parallel loop, so parallelism is never nested.
evaluate_est_met <- function(est_met, est_met_name, x_data,
                             true_moments, true_rank_corr, config, i, dgp_name) {
  # An uncaught error in any est() is recorded as a failed fit and logged, so one bad fit
  # cannot abort the run.
  fit_start  <- Sys.time()
  fit_result <- tryCatch(
    est_met$est(x_data),
    error = function(e) {
      log_line(config, sprintf(
        "est() error | dgp=%s model=%s T=%s replication=%d base_seed=%s | %s",
        dgp_name, est_met_name, config$T, i, config$seed, conditionMessage(e)))
      list(theta = NULL, status = "failed", stationary = NA,
           fit_time = as.numeric(difftime(Sys.time(), fit_start, units = "secs")))
    }
  )
  status <- fit_result$status
  stationary <- fit_result$stationary
  fit_time <- fit_result$fit_time

  # Single definition of a results row, shared by the early returns and the normal path.
  row <- function(statistic, lag, value) {
    data.frame(i = i, est_met_name = est_met_name, statistic = statistic,
               lag = lag, value = value, fit_status = status,
               stationary = stationary, fit_time = fit_time, stringsAsFactors = FALSE)
  }

  # Only converged fits are scored. "hit_max_iter" may still return a theta, so checking for a
  # NULL theta alone would let an unconverged fit contribute GoF values. Unscored fits give NA
  # rows with the same shape as normal rows, which the tapply-based aggregation relies on.
  if (is.null(fit_result$theta) || status != "converged") {
    return(rbind(
      row("GoF1", NA_integer_, NA_real_),
      row("GoF2", as.integer(config$lags), NA_real_),
      row("GoF3", as.integer(config$lags), NA_real_)
    ))
  }

  # Non-stationary fits are still simulated and scored; `stationary` is only a diagnostic column.

  # Forecasts start from the last est_met$p columns of the observed series. No burn-in here:
  # it would discard that conditioning on the real history.
  hist_w <- est_met$p
  x_init_est <- x_data[, (ncol(x_data) - hist_w + 1):ncol(x_data), drop = FALSE]

  # A forecast that diverges mid-path signals "mdacp_forecast_breakdown"; only that condition is
  # absorbed (the path is dropped and a warning issued) and other errors propagate. Pooling below
  # divides by the number of surviving paths.
  sim_results <- lapply(seq_len(config$K), function(k) {
    tryCatch(
      est_met$sim(fit_result$theta, x_init = x_init_est, T = config$T),
      mdacp_forecast_breakdown = function(e) {
        warning(sprintf(
          "evaluate_est_met(): '%s' forecast breakdown on replicate i=%d, k=%d/%d -- %s",
          est_met_name, i, k, config$K, conditionMessage(e)
        ), call. = FALSE)
        NULL
      }
    )
  })
  sim_paths <- Filter(Negate(is.null), sim_results)

  # If every forecast broke down, fit_status stays as reported by est(): it describes the fit,
  # not the forecast.
  if (length(sim_paths) == 0L) {
    return(rbind(
      row("GoF1", NA_integer_, NA_real_),
      row("GoF2", as.integer(config$lags), NA_real_),
      row("GoF3", as.integer(config$lags), NA_real_)
    ))
  }

  # Per-path moments, averaged over the K paths.
  moms       <- lapply(sim_paths, sample_moments, lags = config$lags)
  mu_list    <- lapply(moms, `[[`, "mu")
  gamma_list <- lapply(moms, `[[`, "gamma")
  pooled_mu    <- Reduce(`+`, mu_list) / length(mu_list)
  pooled_gamma <- average_gamma(gamma_list)

  rank_corrs        <- lapply(sim_paths, sample_rank_corr, lags = config$lags)
  pooled_rank_corr  <- average_gamma(rank_corrs)

  # GoF1: mean bias, GoF2: autocovariance bias, GoF3: rank-correlation bias.
  gof1 <- mean_bias(true_moments$mu, pooled_mu)
  gof2 <- gamma_bias(true_moments$gamma, pooled_gamma)
  gof3 <- rank_corr_bias(true_rank_corr, pooled_rank_corr)

  rbind(
    row("GoF1", NA_integer_, gof1),
    row("GoF2", as.integer(names(gof2)), as.numeric(gof2)),
    row("GoF3", as.integer(names(gof3)), as.numeric(gof3))
  )
}


# Full GoF computation for one DGP; returns the long-format data frame with columns
# i, est_met_name, statistic, lag, value, fit_status, stationary, fit_time.
#   dgp: the DGP model; theta_true: its true parameters; est_mets: named list of Est.Met models.
# Ground truth without a closed form (rank correlation always; gamma for the S-vine) is the
# average over the n DGP paths, so the function runs in two passes:
#   1. Simulate and cache all n paths with their moments, serially, because the ground truth
#      needs every path before pass 2 starts.
#   2. Evaluate the Est.Mets, parallelized across replications, reusing the cached paths.
compute_gofs <- function(dgp, theta_true, est_mets, config) {
  T_total <- config$T

  d <- config$d

  # NULL x_init gives a stationary start. A zero start is far from stationary for the S-vine
  # DGP and its transient decays too slowly for the burn-in to remove.
  x_init_dgp      <- NULL
  x_full_list    <- vector("list", config$n)
  moments_cache   <- vector("list", config$n)
  rank_corr_cache <- vector("list", config$n)
  for (i in seq_len(config$n)) {
    x_full <- simulate_with_burnin(dgp, theta_true, x_init = x_init_dgp,
                                   T = T_total, T_burnin = config$T_burnin)
    x_full_list[[i]]     <- x_full
    moments_cache[[i]]   <- sample_moments(x_full, lags = config$lags)
    rank_corr_cache[[i]] <- sample_rank_corr(x_full, lags = config$lags)
  }

  true_moments <- dgp$true_moments(theta_true)
  if (is.null(true_moments$mu)) {
    mu_list <- lapply(moments_cache, `[[`, "mu")
    true_moments$mu <- Reduce(`+`, mu_list) / config$n
  }
  if (is.null(true_moments$gamma)) {
    gamma_list <- lapply(moments_cache, `[[`, "gamma")
    true_moments$gamma <- average_gamma(gamma_list)
  }
  true_rank_corr <- average_gamma(rank_corr_cache)

  # Degenerate DGP paths: paths with a constant row, and per lag the paths whose rank-correlation
  # matrix has an NA entry (ignored entrywise in the reference above).
  n_constant <- sum(vapply(x_full_list, function(x) {
    any(apply(x, 1, function(r) all(r == r[1])))
  }, logical(1)))
  n_dropped  <- vapply(names(rank_corr_cache[[1]]), function(h) {
    sum(vapply(rank_corr_cache, function(rc) anyNA(rc[[h]]), logical(1)))
  }, integer(1))
  log_line(config, sprintf("DGP %s: paths with constant row = %d; paths dropped from rank-corr reference per lag: %s",
                           dgp$name, n_constant,
                           paste0("lag", names(n_dropped), "=", n_dropped, collapse = ", ")))

  # future.seed = TRUE gives each replication its own L'Ecuyer stream, assigned by position.
  results_list <- future.apply::future_lapply(seq_len(config$n), function(i) {
    worker_setup()

    x_full <- x_full_list[[i]]

    do.call(rbind, lapply(names(est_mets), function(nm) {
      evaluate_est_met(est_mets[[nm]], nm, x_full,
                       true_moments, true_rank_corr, config, i, dgp$name)
    }))
  }, future.seed = TRUE)

  do.call(rbind, results_list)
}


# Each takes the long-format results_df and returns one summary. Console output and the PNG
# tables both use them, so the two cannot drift apart.

# Mean over finite entries only; NA if none are finite.
finite_mean <- function(x) {
  fx <- x[is.finite(x)]
  if (length(fx) == 0L) NA_real_ else mean(fx)
}

# Rows = est_met, columns = GoF1, GoF2_lag*, GoF3_lag*; cells = finite_mean() over replications.
compute_gof_table <- function(results_df) {
  df <- results_df
  df$stat_label <- ifelse(df$statistic %in% c("GoF2", "GoF3"),
                          paste0(df$statistic, "_lag", df$lag),
                          df$statistic)

  # tapply (not aggregate) keeps all-NA columns instead of dropping them.
  tab <- tapply(df$value, list(df$est_met_name, df$stat_label), FUN = finite_mean)

  cols <- colnames(tab)
  gof2 <- grep("^GoF2_lag", cols, value = TRUE)
  gof2 <- gof2[order(as.integer(sub("GoF2_lag", "", gof2)))]
  gof3 <- grep("^GoF3_lag", cols, value = TRUE)
  gof3 <- gof3[order(as.integer(sub("GoF3_lag", "", gof3)))]
  ordered <- c(intersect("GoF1", cols), gof2, gof3)
  tab[, ordered, drop = FALSE]
}

# Same shape as compute_gof_table(); cells count the non-finite values excluded from each mean.
compute_nonfinite_counts <- function(results_df) {
  df <- results_df
  df$stat_label <- ifelse(df$statistic %in% c("GoF2", "GoF3"),
                          paste0(df$statistic, "_lag", df$lag),
                          df$statistic)
  tapply(df$value, list(df$est_met_name, df$stat_label),
        FUN = function(x) sum(!is.finite(x)))
}

# Fit-outcome rates use all replications as denominator and sum to 1 per est_met across
# converged, timed out and other. Stationarity and fit time are computed over converged fits only.

# Proportion of replications with fit_status "converged".
compute_convergence_rate <- function(results_df) {
  diag_df <- unique(results_df[, c("i", "est_met_name", "fit_status")])
  tapply(diag_df$fit_status, diag_df$est_met_name,
        FUN = function(s) mean(s == "converged"))
}

# Proportion with fit_status neither "converged" nor "timed_out".
compute_nonconvergence_rate <- function(results_df) {
  diag_df <- unique(results_df[, c("i", "est_met_name", "fit_status")])
  tapply(diag_df$fit_status, diag_df$est_met_name,
        FUN = function(s) mean(!(s %in% c("converged", "timed_out"))))
}

# Proportion with fit_status "timed_out".
compute_timeout_rate <- function(results_df) {
  diag_df <- unique(results_df[, c("i", "est_met_name", "fit_status")])
  tapply(diag_df$fit_status, diag_df$est_met_name,
        FUN = function(s) mean(s == "timed_out"))
}

# Proportion of converged fits with stationary == TRUE; NA for an est_met with no converged fit.
# Pre-filled with NA because tapply omits empty groups.
compute_stationarity_rate <- function(results_df) {
  diag_df <- unique(results_df[, c("i", "est_met_name", "fit_status", "stationary")])
  diag_df <- diag_df[diag_df$fit_status == "converged", ]
  est_met_names <- unique(results_df$est_met_name)
  out <- setNames(rep(NA_real_, length(est_met_names)), est_met_names)
  if (nrow(diag_df) > 0) {
    r <- tapply(diag_df$stationary, diag_df$est_met_name, FUN = mean)
    out[names(r)] <- r
  }
  out
}

# Matrix (mean, sd, min, max) of fit time per est_met over converged fits only; NA where none
# converged (sd is also NA for a single fit).
compute_fit_time_table <- function(results_df) {
  time_df <- unique(results_df[, c("i", "est_met_name", "fit_status", "fit_time")])
  time_df <- time_df[time_df$fit_status == "converged", ]

  est_met_names <- unique(results_df$est_met_name)
  na_col <- setNames(rep(NA_real_, length(est_met_names)), est_met_names)
  fit_time_mean <- na_col; fit_time_sd <- na_col; fit_time_min <- na_col; fit_time_max <- na_col

  if (nrow(time_df) > 0) {
    m  <- tapply(time_df$fit_time, time_df$est_met_name, FUN = mean)
    s  <- tapply(time_df$fit_time, time_df$est_met_name, FUN = sd)
    mn <- tapply(time_df$fit_time, time_df$est_met_name, FUN = min)
    mx <- tapply(time_df$fit_time, time_df$est_met_name, FUN = max)
    fit_time_mean[names(m)]  <- m
    fit_time_sd[names(s)]    <- s
    fit_time_min[names(mn)]  <- mn
    fit_time_max[names(mx)]  <- mx
  }
  cbind(mean = fit_time_mean, sd = fit_time_sd, min = fit_time_min, max = fit_time_max)
}


# Console summary for one DGP's results; formatting only.
print_output <- function(results_df) {
  tab <- compute_gof_table(results_df)
  cat("GoF summary (mean over i):\n")
  print(round(tab, 4))

  # Shown only when non-finite values were excluded from the means above.
  nonfinite_counts <- compute_nonfinite_counts(results_df)
  if (any(nonfinite_counts > 0, na.rm = TRUE)) {
    cat("\nWARNING: non-finite (NA/NaN/Inf) values excluded from the means above\n")
    cat("(count per est_met/statistic):\n")
    print(nonfinite_counts)
  }

  stationarity_rate <- compute_stationarity_rate(results_df)
  cat("\nStationarity rate (proportion of CONVERGED replicates with stationary == TRUE; NA if none converged):\n")
  print(round(stationarity_rate, 4))

  convergence_rate    <- compute_convergence_rate(results_df)
  nonconvergence_rate <- compute_nonconvergence_rate(results_df)
  timeout_rate        <- compute_timeout_rate(results_df)
  cat("\nFit outcome rates (proportion of ALL replicates; Convergence + Non-Convergence + Timeout = 1):\n")
  print(round(cbind(`Convergence Rate` = convergence_rate,
                    `Non-Convergence Rate` = nonconvergence_rate,
                    `Timeout Rate` = timeout_rate), 4))

  fit_time_tab <- compute_fit_time_table(results_df)
  cat("\nFit time diagnostic, seconds (per est_met, over CONVERGED replicates only; NA if none converged):\n")
  print(round(fit_time_tab[, c("mean", "sd"), drop = FALSE], 4))

  invisible(tab)
}


# Both build a data frame from the shared aggregation helpers, render it with tableGrob() and
# save it with ggsave() into output_dir. "Est.Met." is an ordinary first column and the title
# names the DGP, so each image is readable on its own.

# Table width in inches from each column's character content; needed because the diagnostics
# headers vary a lot in length.
estimate_table_width_in <- function(tab_df, in_per_char = 0.11, padding_in = 0.35) {
  col_widths <- vapply(seq_along(tab_df), function(j) {
    max(nchar(colnames(tab_df)[j]), nchar(format(tab_df[[j]])))
  }, numeric(1))
  sum(col_widths) * in_per_char + padding_in * ncol(tab_df)
}

# Shared builder: table grob with a title above it, saved as PNG.
save_table_png <- function(tab_df, title, out_path, width, height) {
  table_grob <- gridExtra::tableGrob(tab_df, rows = NULL)
  title_grob <- grid::textGrob(title, gp = grid::gpar(fontsize = 14, fontface = "bold"))
  combined <- gridExtra::arrangeGrob(table_grob, top = title_grob)

  ggsave(out_path, plot = combined, width = width, height = height,
        units = "in", dpi = 300)
}

# Table 1: GoF summary. Columns follow compute_gof_table(), so a wider config$lags needs no change.
save_gof_table_png <- function(results_df, dgp_name, output_dir) {
  tab <- round(compute_gof_table(results_df), 4)

  tab_df <- data.frame(`Est.Met.` = rownames(tab), as.data.frame(tab),
                       check.names = FALSE, row.names = NULL)

  save_table_png(
    tab_df,
    title = sprintf("%s — GoF Summary", dgp_name),
    out_path = file.path(output_dir, paste0(dgp_name, "_gof_table.png")),
    width = 2 + 1.1 * ncol(tab), height = 1.6 + 0.4 * nrow(tab)
  )
}

# Table 2: stationarity, outcome rates and fit times.
save_diagnostics_table_png <- function(results_df, dgp_name, output_dir) {
  stationarity_rate    <- compute_stationarity_rate(results_df)
  convergence_rate     <- compute_convergence_rate(results_df)
  nonconvergence_rate  <- compute_nonconvergence_rate(results_df)
  timeout_rate         <- compute_timeout_rate(results_df)
  fit_time_tab         <- compute_fit_time_table(results_df)

  # Index everything by the same names so rows line up.
  est_met_names <- rownames(fit_time_tab)
  tab_df <- data.frame(
    `Est.Met.`             = est_met_names,
    `Stationarity Rate`    = round(stationarity_rate[est_met_names], 4),
    `Convergence Rate`     = round(convergence_rate[est_met_names], 4),
    `Non-Convergence Rate` = round(nonconvergence_rate[est_met_names], 4),
    `Timeout Rate`         = round(timeout_rate[est_met_names], 4),
    `Mean Fit Time (s)`    = round(fit_time_tab[est_met_names, "mean"], 4),
    `Fit Time SD (s)`      = round(fit_time_tab[est_met_names, "sd"], 4),
    check.names = FALSE, row.names = NULL
  )

  save_table_png(
    tab_df,
    title = sprintf("%s — Diagnostics", dgp_name),
    out_path = file.path(output_dir, paste0(dgp_name, "_diagnostics_table.png")),
    width = estimate_table_width_in(tab_df), height = 1.6 + 0.4 * nrow(tab_df)
  )
}


config <- get_config()

# multisession works on every OS (unlike forking). The plan is reset to sequential on exit to
# stop the workers, including after an error.
future::plan(future::multisession, workers = config$n_workers)

set.seed(config$seed)



# minar1_poi and minar1_negbin share A so the innovation family is the only difference.
A_minar1 <- matrix(c(0.30, 0.15, 0.09,
                     0.18, 0.27, 0.23,
                     0.06, 0.10, 0.25), nrow = 3, byrow = TRUE)

lambda_minar1 <- c(1, 1, 5, 3)
theta_true_minar1_poi <- list(
  A = list(A_minar1), innovation_params = list(lambda = lambda_minar1)
)

# NegBin innovation means match the Poisson ones, so only dispersion differs.
r_minar1 <- 5
p_minar1 <- c(0.15, 0.3, 0.25)
theta_true_minar1_negbin <- list(
  A = list(A_minar1), innovation_params = list(r = r_minar1, p = p_minar1)
)


# mdacp1: the top 2 x 2 block of delta is the identity; only row 3 is free.
gamma_mdacp1 <- matrix(c(0.12, 0.10,
                         0.09, 0.15,
                         0.10, 0.08), nrow = 3, ncol = 2, byrow = TRUE)
delta_mdacp1 <- matrix(c(1,    0,
                         0,    1,
                         0.6,  0.3), nrow = 3, ncol = 2, byrow = TRUE)
alpha_mdacp1 <- c(0.12, 0.25, 0.18)
beta_mdacp1  <- c(0.25, 0.35, 0.40)
phi_mdacp1   <- c(2.25, 2.15, 1.85)
omega_mdacp1 <- c(2.25, 1.75, 3.00)
Sigma_mdacp1 <- matrix(c(1.00,  0.25, -0.20,
                         0.25,  1.00,  0.15,
                        -0.20,  0.15,  1.00), nrow = 3)
theta_true_mdacp1 <- list(
  omega = omega_mdacp1, gamma = gamma_mdacp1, delta = delta_mdacp1,
  alpha = alpha_mdacp1, beta = beta_mdacp1, phi = phi_mdacp1, Sigma = Sigma_mdacp1
)


# svine1: unlike the other models, theta is a svine_dist object, which svine() would only
# produce by fitting, so it is built directly from rvinecopulib's and svines' constructors.
# The extended vine has (p + 1) * d0 variables. For tree t, svinecop_dist() uses the first
# min(d0, d - t) pair copulas and tiles them across the remaining time-shifted copies, which is
# what makes the vine stationary.
# Input slot i of tree t lands at output edge ((i - 1 + t) %% d0) + 1, so the slots below are
# ordered to hit these edges:
#   tree 1, slot 1 -> V3-V2 at the same time
#   tree 1, slot 2 -> V2-V1 at the same time
#   tree 1, slot 3 -> V1 lag-1 serial edge
#   tree 4, slot 3 -> V1 lag-2 edge, conditioned on the full previous cross-section
# All four are Clayton copulas.
# All other slots are independence copulas, so the DGP has exactly these four non-trivial
# dependencies (an S-vine cannot be truncated).
# Margins are Poisson. Each is a univariateML object fitted to throwaway
# data and then overwritten with the chosen parameters via `[<-`, so it keeps the attributes
# svines needs.
d0_svine1 <- 3L; p_svine1 <- 2L; d_svine1 <- d0_svine1 * (p_svine1 + 1L)

margin_v1_svine1 <- univariateML::mlpois(rpois(20, 5))
margin_v1_svine1[] <- 5
margin_v2_svine1 <- univariateML::mlpois(rpois(20, 7))
margin_v2_svine1[] <- 7
margin_v3_svine1 <- univariateML::mlpois(rpois(20, 9))
margin_v3_svine1[] <- 9

indep_dd_svine1 <- rvinecopulib::bicop_dist("indep", var_types = c("d", "d"))
tree_lengths_svine1 <- vapply(seq_len(d_svine1 - 1L),
                              function(t) min(d0_svine1, d_svine1 - t), integer(1))
pair_copulas_svine1 <- lapply(tree_lengths_svine1,
                              function(n) rep(list(indep_dd_svine1), n))
pair_copulas_svine1[[1]][[1]] <- rvinecopulib::bicop_dist(
  "clayton", 0, 2,   var_types = c("d", "d"))
pair_copulas_svine1[[1]][[2]] <- rvinecopulib::bicop_dist(
  "clayton", 0, 2.5, var_types = c("d", "d"))
pair_copulas_svine1[[1]][[3]] <- rvinecopulib::bicop_dist(
  "clayton", 0, 3,   var_types = c("d", "d"))

pair_copulas_svine1[[4]][[3]] <- rvinecopulib::bicop_dist(
  "clayton", 0, 2,   var_types = c("d", "d"))

cs_structure_svine1 <- rvinecopulib::dvine_structure(1:3)
copula_svine1 <- svines::svinecop_dist(
  pair_copulas_svine1, cs_structure_svine1, p = p_svine1,
  out_vertices = 1:3, in_vertices = 1:3, var_types = rep("d", d0_svine1)
)
theta_true_svine1 <- svines::svine_dist(
  margins = list(margin_v1_svine1, margin_v2_svine1, margin_v3_svine1),
  copula  = copula_svine1
)


# mcar1: DGP only. A and B are asymmetric, and the Gaussian copula has distinct negative
# correlations per pair. All entries of omega, A, B must be positive and the stationarity
# condition is the matrix 1-norm sum documented in model_multiv_count_ar.R.
A_mcar1 <- matrix(c(0.25, 0.17, 0.15,
                    0.12, 0.20, 0.19,
                    0.14, 0.13, 0.15), nrow = 3, byrow = TRUE)
B_mcar1 <- matrix(c(0.20, 0.08, 0.09,
                    0.05, 0.15, 0.03,
                    0.07, 0.04, 0.30), nrow = 3, byrow = TRUE)
omega_mcar1 <- c(0.5, 1.0, 2.75)
theta_true_mcar1 <- list(
  omega = omega_mcar1, A = A_mcar1, B = B_mcar1,
  copula_param = c(-0.30, -0.45, -0.20)
)


# DGPs and Est.Mets, selected by name from config and resolved through the registry
theta_lookup <- list(minar1_poi    = theta_true_minar1_poi,
                     minar1_negbin = theta_true_minar1_negbin,
                     mdacp1        = theta_true_mdacp1,
                     mcar1         = theta_true_mcar1,
                     svine1        = theta_true_svine1)

DGPs  <- setNames(lapply(config$dgp_names, function(nm) all_models[[nm]]), config$dgp_names)
Theta <- theta_lookup[config$dgp_names]

est_mets <- setNames(lapply(config$est_met_names, function(nm) all_models[[nm]]),
                     config$est_met_names)

dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)

tryCatch({
  for (i in seq_along(DGPs)) {
    dgp_name <- names(DGPs)[i]
    cat(sprintf("\n=== DGP: %s ===\n", dgp_name))

    result <- compute_gofs(DGPs[[i]], Theta[[i]], est_mets, get_config())
    print_output(result)

    save_gof_table_png(result, dgp_name, config$output_dir)
    save_diagnostics_table_png(result, dgp_name, config$output_dir)

    # Saved per DGP so a later failure cannot lose earlier results.
    out_path <- file.path(config$output_dir, paste0(dgp_name, ".rds"))
    saveRDS(result, out_path)
    cat(sprintf("saved: %s\n", out_path))
  }
}, finally = {
  future::plan(future::sequential)
})
