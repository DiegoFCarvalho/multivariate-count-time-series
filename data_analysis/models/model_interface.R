# Model object: list with slots sim(theta, x_init, T), est(x_data) and the history length p.
# est() returns list(theta, status, stationary, fit_time).

# Seconds allowed per est() call, set with options(est_timeout_seconds = ...).
est_timeout_seconds <- function() getOption("est_timeout_seconds", 600)

# Runs `expr` under an elapsed-time limit; returns list(timed_out, value).
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

new_model <- function(name, sim, p, est = NULL) {
  stopifnot(is.function(sim), length(p) == 1L, p == as.integer(p), p > 0L)
  if (!is.null(est)) stopifnot(is.function(est))
  model <- list(name = name, p = as.integer(p), sim = sim, est = est)
  class(model) <- "forecast_model"
  model
}
