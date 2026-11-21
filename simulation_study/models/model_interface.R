# Contract for every DGP / Est.Met. model (slots sim, est, true_moments).
# The simulation loop never branches on a model's name.
# Notation: d = number of components, T = number of time points, p = lags of history.
# Every model declares its exact history length p; history matrices (x_init) are d x p,
# with past observations in the columns.
#
#   sim(theta, x_init, T)
#     Returns a d x T count matrix. x_init has exactly p columns; its meaning is
#     model-specific. x_init = NULL means a model-specific stationary start. Required.
#
#   est(x_data)
#     Fits the model to a d x T count matrix. Returns a list with:
#       theta      : fitted parameters (model-specific), NULL if status is "failed"/"timed_out".
#       status     : "converged", "failed", "hit_max_iter" or "timed_out".
#       stationary : TRUE/FALSE, NA if undetermined or no theta.
#       fit_time   : seconds elapsed (the timeout value for "timed_out"); never omitted.
#
#   true_moments(theta)
#     Returns list(mu, gamma): d-vector of means and a list of d x d autocovariance
#     matrices per lag. Either may be NULL; values may be approximate.


# Runs `expr` under a base-R elapsed-time limit; returns list(timed_out, value).
# Ordinary errors are re-thrown. A timeout is detected by elapsed time, not by message
# text, since the message is locale-dependent. setTimeLimit() cannot interrupt a single
# long compiled call, so models whose fit is compiled code need another mechanism.
with_est_timeout <- function(expr, timeout_seconds) {
  setTimeLimit(elapsed = timeout_seconds, transient = TRUE)
  on.exit(setTimeLimit())
  t_start <- Sys.time()

  tryCatch(
    list(timed_out = FALSE, value = expr),
    error = function(e) {
      elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
      if (elapsed >= timeout_seconds - 0.5) {
        list(timed_out = TRUE, value = NULL)
      } else {
        stop(e)
      }
    }
  )
}


# Builds a model object and checks that `sim` is a function, `p` a positive integer,
# and the optional slots functions or NULL.
new_model <- function(name,
                      sim,
                      p,
                      est = NULL,
                      true_moments = NULL) {

  if (missing(sim) || !is.function(sim)) {
    stop("new_model(): `sim` must be provided and must be a function.")
  }

  if (missing(p) || length(p) != 1L || is.na(p) ||
      p != as.integer(p) || p <= 0L) {
    stop("new_model(): `p` must be a single positive integer.")
  }
  p <- as.integer(p)

  if (!is.null(est) && !is.function(est)) {
    stop("new_model(): `est` must be a function or NULL.")
  }
  if (!is.null(true_moments) && !is.function(true_moments)) {
    stop("new_model(): `true_moments` must be a function or NULL.")
  }

  model <- list(
    name         = name,
    p            = p,
    sim          = sim,
    est          = est,
    true_moments = true_moments
  )
  class(model) <- "sim_study_model"
  model
}
