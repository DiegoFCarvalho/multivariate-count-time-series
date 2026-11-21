# Forecast draws of the Kim et al. (2025) latent Gaussian dynamic factor model, for every group.
# Defines functions only. Writes <results_dir>/<data_prefix>_<group>/raw/{raw,fit}_kim_K<K>.rds.

load_kim_code <- function() {
  code_dir <- file.path(kim_repo_dir, "Code")
  suppressMessages({
    library(polynom); library(mvtnorm); library(Matrix)
  })
  load(file.path(code_dir, "Polys.rda"), envir = .GlobalEnv)
  for (f in c("Latent_Gauss_Cov.R", "Latent_Gauss_Link.R", "Latent_Gauss_invLink.R",
              "Latent_DFM_Model.R", "Latent_DFM_Estim.R", "Latent_DFM_Forecast.R"))
    source(file.path(code_dir, f), chdir = TRUE)
}

symmetrize_matrix <- function(S) (S + t(S)) / 2

set_stream <- function(state) assign(".Random.seed", state, envir = .GlobalEnv)

# n forecast paths of H days: sample a particle by its weight, draw the factor state from the
# filtered Gaussian, run the factor VAR forward, add idiosyncratic noise, and map the latent value
# to a rating with the estimated thresholds. Returns an integer array [n, d, H].
kim_draws <- function(parts, model, r, n, H) {
  d   <- nrow(model$Lambda)
  Psi <- matrix(model$Psi[, , 1], r, r)
  j   <- sample.int(length(parts$w_T), n, replace = TRUE, prob = parts$w_T)
  Y   <- t(parts$Y_tilde_h[, j, drop = FALSE]) +
         mvtnorm::rmvnorm(n, rep(0, r), symmetrize_matrix(matrix(parts$Q_tilde_h, r, r)), method = "eigen")
  x <- array(NA_integer_, c(n, d, H))
  for (h in seq_len(H)) {
    Y <- Y %*% t(Psi) + mvtnorm::rmvnorm(n, rep(0, r), symmetrize_matrix(matrix(model$Cov_eta, r, r)), method = "eigen")
    Z <- Y %*% t(model$Lambda) + mvtnorm::rmvnorm(n, rep(0, d), symmetrize_matrix(model$Cov_eps), method = "eigen")
    for (i in seq_len(d)) x[, i, h] <- 1L + findInterval(Z[, i], parts$threshold[[i]], left.open = TRUE)
  }
  x
}

# fit and K draws for one group; returns list(raw_result, fit_record)
kim_group <- function(g, K, seeds_kim) {
  x     <- read_group(g)
  d     <- nrow(x)
  stopifnot(d == dimension, ncol(x) == n_days)
  X_t   <- x[, seq_len(n_train), drop = FALSE]
  n_is  <- n_train - p_max_seed
  B     <- 2L + n_is + horizon

  setTimeLimit(elapsed = est_timeout, transient = TRUE)
  on.exit(setTimeLimit(), add = TRUE)

  # estimation and filter use the first stream state
  set_stream(seeds_kim[[1L]])
  t0    <- Sys.time()
  Link  <- Latent_Gauss_Link(d, n_train, "multinom", X_t, n_cat)
  Cov_X <- Latent_Gauss_Cov(d, n_train, kim_p, X_t)
  Cov_Z <- Latent_Gauss_invLink(d, kim_p, Cov_X, Link$Link)
  Est   <- Latent_DFM_Estim(r = kim_r, kim_p, Cov_Z, identy_opt = 2, shift = TRUE)
  model <- list(Lambda = Est$Lambda, Cov_eps = diag(abs(Est$Cov_eps), nrow = d),
                Psi = Est$Psi, Cov_eta = Est$Cov_eta, Param = Link$Param)
  t_fit <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  t1 <- Sys.time()
  F  <- Latent_DFM_Forecast(model, kim_p, kim_r, d, n_train, X_t, horizon, "multinom", n_cat)
  t_fc <- as.numeric(difftime(Sys.time(), t1, units = "secs"))

  win   <- dim(F$Y_tilde)[2]
  parts <- list(w_T = F$w_t[win, ], Y_tilde_h = matrix(F$Y_tilde[, win, ], nrow = kim_p * kim_r),
                Q_tilde_h = matrix(F$Q_tilde[, , win], nrow = kim_p * kim_r), threshold = F$threshold)

  # one stream state per draw
  raw <- array(NA_real_, c(K, d, horizon))
  for (k in seq_len(K)) {
    set_stream(seeds_kim[[(k - 1L) * B + n_is + 1L + 2L]])
    raw[k, , ] <- kim_draws(parts, model, kim_r, 1L, horizon)[1, , ]
  }

  raw_result <- list(model = "kim", mode = "oos_free", K = K, scored_cols = n_train + seq_len(horizon),
                     series_names = rownames(x), truth = x[, n_train + seq_len(horizon), drop = FALSE],
                     raw = raw, fit_status = "converged")
  fit_record <- list(model = "kim", r = kim_r, p = kim_p, n_train = n_train, K = K, base_seed = base_seed,
                     stream_id = stream_ids[["kim"]], p_max = p_max_seed, ESS = F$ESS,
                     Lambda = model$Lambda, Psi = model$Psi, Cov_eta = model$Cov_eta,
                     Cov_eps = diag(model$Cov_eps), Param = Link$Param, threshold = F$threshold,
                     seconds = c(fit = t_fit, forecast = t_fc))
  list(raw_result = raw_result, fit_record = fit_record)
}

run_kim_stage <- function(groups_sel = names(groups), K) {
  load_kim_code()
  B <- 2L + (n_train - p_max_seed) + horizon
  seeds_kim <- make_stream_seeds(base_seed, max(stream_ids), K, B)[[stream_ids[["kim"]]]]
  for (g in groups_sel) {
    if (skip_if_exists && file.exists(raw_file("kim", g, K))) next
    t0  <- Sys.time()
    res <- tryCatch(kim_group(g, K, seeds_kim), error = function(e) e)
    if (inherits(res, "error")) { message(sprintf("[kim] %s FAILED: %s", g, conditionMessage(res))); next }
    dir.create(group_raw_dir(g), showWarnings = FALSE, recursive = TRUE)
    saveRDS(res$raw_result, raw_file("kim", g, K))
    saveRDS(res$fit_record, fit_file("kim", g, K))
    message(sprintf("[kim] %s done in %.0fs (minimum ESS %.0f)", g,
                    as.numeric(difftime(Sys.time(), t0, units = "secs")), min(res$fit_record$ESS)))
  }
  invisible(NULL)
}
