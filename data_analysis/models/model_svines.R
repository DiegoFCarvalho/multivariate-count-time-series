# S-vine model of Markov order p. Data are d x T; svines uses T x d, so data are transposed.
# margins = "exact_ecdf": relative-frequency margins on 1..5 as in Latent_Gauss_Link.R, built from exported
#   svines functions (svinecop, svine_dist).
# margins = "scaled_ecdf": svines' built-in empirical margins (margin_families = "empirical").

# Discrete margin of one item from its training values; F(5) = 1, no n/(n+1) factor.
kim_margin <- function(x) {
  stopifnot(all(x %in% 1:5))
  cf <- cumsum(tabulate(x, 5L)) / length(x); cf[5] <- 1
  m <- list(p     = function(v) cf[v],
            p_sub = function(v) c(0, cf)[v],
            q     = function(u) 1L + findInterval(u, cf, left.open = TRUE))
  structure(m, type = "discrete", class = c("list", "svine_margin"), cf = cf, model = "empirical")
}

# pseudo-observations [F(x), F(x-)] in the column layout svinecop() expects for discrete variables
kim_pseudo_obs <- function(xt, margins) {
  cbind(vapply(seq_along(margins), function(j) margins[[j]]$p(xt[, j]), numeric(nrow(xt))),
        vapply(seq_along(margins), function(j) margins[[j]]$p_sub(xt[, j]), numeric(nrow(xt))))
}

fit_kim_svine <- function(xt, p) {
  margins <- lapply(seq_len(ncol(xt)), function(j) kim_margin(xt[, j]))
  names(margins) <- colnames(xt)
  copula <- svines::svinecop(kim_pseudo_obs(xt, margins), p = p, var_types = rep("d", ncol(xt)))
  svines::svine_dist(margins, copula)
}

make_svines_model <- function(name, p, margins = c("exact_ecdf", "scaled_ecdf")) {
  margins <- match.arg(margins)

  # theta is a svine_dist object; x_init is d x p
  sim <- function(theta, x_init, T) {
    past <- if (is.null(x_init)) NULL else t(x_init)
    out <- svines::svine_sim(n = T, rep = 1, model = theta, past = past)
    if (length(dim(out)) == 3L) out <- out[, , 1L]
    t(out)
  }

  fit_svine <- function(x_data) {
    d <- nrow(x_data)
    switch(margins,
      exact_ecdf  = fit_kim_svine(t(x_data), p),
      scaled_ecdf = svines::svine(data = t(x_data), p = p, var_types = rep("d", d), margin_families = "empirical"))
  }

  est <- function(x_data) {
    start_time <- Sys.time()

    timeout_seconds <- est_timeout_seconds()
    to <- tryCatch(
      with_est_timeout(fit_svine(x_data), timeout_seconds),
      error = function(e) e
    )

    if (inherits(to, "error")) {
      fit_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      return(list(theta = NULL, status = "failed", stationary = NA, fit_time = fit_time))
    }
    if (to$timed_out) {
      return(list(theta = NULL, status = "timed_out", stationary = NA, fit_time = timeout_seconds))
    }
    fit_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    list(theta = to$value, status = "converged", stationary = TRUE, fit_time = fit_time)
  }

  new_model(name = name, sim = sim, p = p, est = est)
}
