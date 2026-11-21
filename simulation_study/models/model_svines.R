# S-vine model of Markov order p (interface in model_interface.R).
# svines uses T x d data, so d x T input is transposed at the boundary.

# Quantiles of one fitted margin, mirroring svines' own dispatch with public functions only.
.svines_margin_quantile <- function(margin, probs) {
  if (inherits(margin, "univariateML")) {
    univariateML::qml(probs, margin)
  } else {
    margin$q(probs)
  }
}

# Marginal mean as the integral of the quantile function over a midpoint grid.
.svines_marginal_mean <- function(margin, n_grid = 1e5L) {
  probs <- (seq_len(n_grid) - 0.5) / n_grid
  mean(.svines_margin_quantile(margin, probs))
}


# S-vines cannot be truncated, so the Markov order p is the only structural argument.
make_svines_model <- function(name, p) {

  # theta is a svine_dist object (fitted by est() or built by hand)
  sim <- function(theta, x_init, T) {
    past <- if (is.null(x_init)) NULL else t(x_init)

    # svine_sim() returns an n x d matrix for rep = 1
    out <- svines::svine_sim(n = T, rep = 1, model = theta, past = past)

    if (length(dim(out)) == 3L) out <- out[, , 1L]
    t(out)
  }

  est <- function(x_data) {
    start_time <- Sys.time()

    d <- nrow(x_data)

    # The fit call is not pre-wrapped in its own tryCatch(), so that with_est_timeout()
    # can tell a timeout from an ordinary error.
    timeout_seconds <- get_config()$est_timeout_seconds
    to <- tryCatch(
      with_est_timeout(
        svines::svine(
          data      = t(x_data),
          p         = p,
          var_types = rep("d", d)
        ),
        timeout_seconds
      ),
      error = function(e) e
    )

    # No convergence flag is exposed, so a completed fit counts as converged.
    if (inherits(to, "error")) {
      fit_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      return(list(theta = NULL, status = "failed", stationary = NA, fit_time = fit_time))
    }
    if (to$timed_out) {
      return(list(theta = NULL, status = "timed_out", stationary = NA, fit_time = timeout_seconds))
    }
    fit <- to$value
    # Stationary by construction.
    fit_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    list(theta = fit, status = "converged", stationary = TRUE, fit_time = fit_time)
  }

  # gamma is NULL: autocovariances have no closed form here.
  true_moments <- function(theta) {
    margins <- theta$margins
    mu <- vapply(margins, .svines_marginal_mean, numeric(1))
    list(mu = mu, gamma = NULL)
  }

  new_model(
    name         = name,
    sim          = sim,
    p            = p,
    est          = est,
    true_moments = true_moments
  )
}
