# Out-of-sample forecast draws of the registry models on one d x T matrix; defines functions only.
# Each model is fitted on the first n_train columns; K free paths over the remaining columns are
# then drawn from the fitted model, conditional on the last p training columns.
# Every draw has its own L'Ecuyer-CMRG stream derived from base_seed (one stream per model slot), so
# results do not depend on K, the worker count or the other models.

# mean absolute difference per column, averaged over the d series (rating points)
abs_err <- function(a, b) colMeans(abs(a - b))

mean_or_na <- function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)

# K x n grid of abs_err per draw and forecast day for raw [K, d, n] and truth [d, n]
abs_err_grid <- function(raw, truth) {
  K <- dim(raw)[1]; d <- dim(raw)[2]; n <- dim(raw)[3]
  t(vapply(seq_len(K), function(k) abs_err(matrix(raw[k, , ], d, n), truth), numeric(n)))
}


# sources the model definitions; est() and the MDACP child process need cwd = root_dir
load_project <- function(root_dir, est_timeout) {
  withr::local_dir(root_dir)
  options(est_timeout_seconds = est_timeout)
  source("models/model_registry.R")
}

# workers are fresh processes
worker_setup <- function(root_dir, est_timeout) {
  if (!exists("all_models", envir = .GlobalEnv, inherits = FALSE)) {
    setwd(root_dir)
    load_project(root_dir, est_timeout)
  }
}

get_model <- function(name) get("all_models", envir = .GlobalEnv)[[name]]

read_data <- function(path) {
  x <- if (grepl("\\.rds$", path, ignore.case = TRUE)) readRDS(path) else as.matrix(read.csv(path, header = FALSE))
  x <- as.matrix(x)
  storage.mode(x) <- "integer"
  x
}

# d x Tn matrix; column s is the MDACP conditional mean mu_s that scores x_s (from the data only)
mdacp_mu_state <- function(theta, x) {
  AB <- .mdacp_build_AB(theta)
  d  <- length(theta$omega)
  mu <- as.vector(solve(diag(d) - AB$A - AB$B, theta$omega))
  out <- matrix(NA_real_, nrow = d, ncol = ncol(x))
  out[, 1L] <- mu
  for (s in seq_len(ncol(x))[-1L]) {
    mu <- as.vector(theta$omega + AB$A %*% x[, s - 1L] + AB$B %*% mu)
    out[, s] <- mu
  }
  out
}


# per model slot, 1 + K*B states: position 1 = fit, then one sub-stream per unit
make_stream_seeds <- function(base_seed, n_models, K, B) {
  set.seed(base_seed, kind = "L'Ecuyer-CMRG")
  stream <- .Random.seed
  lapply(seq_len(n_models), function(m) {
    stream <<- parallel::nextRNGStream(stream)
    states <- vector("list", 1L + K * B)
    s <- stream
    states[[1L]] <- s
    for (j in seq_len(K * B)) {
      s <- parallel::nextRNGSubStream(s)
      states[[j + 1L]] <- s
    }
    states
  })
}


# one free path over the test columns from the last p training columns; returns d x (Tn - N)
draw_unit <- function(model, mc, x, N, Tn) {
  p  <- model$p
  th <- mc$theta
  if (mc$stateful) th$mu_init <- mc$mu_state[, N]
  out <- model$sim(th, x[, (N - p + 1L):N, drop = FALSE], Tn - N)
  if (any(!is.finite(out))) stop("non-finite simulated values")
  matrix(out, nrow = nrow(x), ncol = Tn - N)
}

fit_unit <- function(m, shared) {
  worker_setup(shared$root_dir, shared$est_timeout)
  get_model(shared$model_names[m])$est(shared$x_train)
}

# NULL if the draw fails
run_unit <- function(u, shared) {
  worker_setup(shared$root_dir, shared$est_timeout)
  m <- shared$u_model[u]
  tryCatch(draw_unit(get_model(shared$model_names[m]), shared$mcs[[m]], shared$x, shared$N, shared$Tn),
           error = function(e) NULL)
}


# Fits `model_names` on the first n_train columns of x and draws K paths over the remaining columns.
# stream_ids gives each model its random-number stream; p_max fixes the seed layout's history.
# Returns list(results, fits); results[[model]] is in the layout of the raw draw files.
run_evaluation <- function(x, root_dir, model_names, n_train, K, n_workers, base_seed, est_timeout,
                           stream_ids, p_max) {
  library(future)
  library(future.apply)
  withr::local_dir(root_dir)
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
  if (!exists("all_models", envir = .GlobalEnv, inherits = FALSE)) load_project(root_dir, est_timeout)
  options(est_timeout_seconds = est_timeout)

  d  <- nrow(x); Tn <- ncol(x); N <- as.integer(n_train); K <- as.integer(K)
  models <- lapply(model_names, get_model)
  stopifnot(p_max >= max(vapply(models, function(m) m$p, integer(1))))
  n_is  <- N - p_max
  n_oos <- Tn - N
  B     <- 2L + n_is + n_oos
  test_cols <- N + seq_len(n_oos)
  seeds <- make_stream_seeds(as.integer(base_seed), max(stream_ids), K, B)[stream_ids]

  if (n_workers > 1L) future::plan(future::multisession, workers = n_workers) else future::plan(future::sequential)
  on.exit(future::plan(future::sequential), add = TRUE)

  # stage 1: fit
  fits <- future_lapply(seq_along(models), fit_unit,
                        shared = list(root_dir = root_dir, est_timeout = est_timeout, model_names = model_names,
                                      x_train = x[, seq_len(N), drop = FALSE]),
                        future.seed = lapply(seeds, `[[`, 1L))
  names(fits) <- model_names
  fit_ok <- vapply(fits, function(f) !is.null(f$theta), logical(1))
  mcs <- lapply(seq_along(models), function(m) {
    if (!fit_ok[m]) return(NULL)
    th <- fits[[m]]$theta
    stateful <- !is.null(th$mu_init)
    list(theta = th, stateful = stateful, mu_state = if (stateful) mdacp_mu_state(th, x) else NULL)
  })

  # stage 2: one unit per model and draw
  units <- do.call(rbind, c(list(data.frame(m = integer(0), k = integer(0))),
                            lapply(which(fit_ok), function(m) data.frame(m = m, k = seq_len(K)))))
  unit_seeds <- lapply(seq_len(nrow(units)), function(u) seeds[[units$m[u]]][[(units$k[u] - 1L) * B + n_is + 1L + 2L]])
  draws <- future_lapply(seq_len(nrow(units)), run_unit,
                         shared = list(root_dir = root_dir, est_timeout = est_timeout, model_names = model_names,
                                       mcs = mcs, x = x, N = N, Tn = Tn, u_model = units$m),
                         future.seed = unit_seeds, future.scheduling = 4)

  raws <- lapply(models, function(m) array(NA_real_, c(K, d, n_oos)))
  for (u in seq_along(draws)) if (!is.null(draws[[u]])) raws[[units$m[u]]][units$k[u], , ] <- draws[[u]]

  results <- lapply(seq_along(models), function(m)
    list(model = model_names[m], mode = "oos_free", K = K, scored_cols = test_cols, series_names = rownames(x),
         truth = x[, test_cols, drop = FALSE], raw = raws[[m]], fit_status = fits[[m]]$status))
  names(results) <- model_names
  list(results = results, fits = fits)
}
