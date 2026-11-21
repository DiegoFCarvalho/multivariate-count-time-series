# Multivariate INAR(p) model (MINAR).
# sim() and true_moments() work for any p; est() (conditional ML) only for p = 1.
# est() needs sample_moments() (gof.R) and get_config() to be sourced by the caller.
# theta = list(A = list of p d x d thinning matrices,
#              innovation_params = family-specific):
#   Poisson: lambda, length d + 1, lambda[1] is the shared factor.
#   NegBin:  list(r, p) with length-d p; p_0 = 1 - sum(p) is derived.
# The innovation family is fixed at construction, not stored in theta.
# Series are d x T; x_init is d x p with the most recent observation in the last column.


# Closed-form mean m and covariance Sigma_E of one innovation vector.
# Poisson: X_i = Y_0 + Y_i with a shared Poisson Y_0.
# NegBin: Gamma-Poisson mixture with a shared Gamma variate.
.minar_innovation_moments <- function(family, params, d) {
  if (family == "poisson") {
    lambda   <- params$lambda
    lambda_0 <- lambda[1]
    lambda_i <- lambda[-1]
    m       <- lambda_0 + lambda_i
    Sigma_E <- lambda_0 * matrix(1, nrow = d, ncol = d) + diag(lambda_i, nrow = d)
  } else if (family == "negbin") {
    r    <- params$r
    pvec <- params$p
    p0   <- 1 - sum(pvec)
    m       <- (r / p0) * pvec
    Sigma_E <- (r / p0^2) * (pvec %*% t(pvec) + p0 * diag(pvec, nrow = d))
  } else {
    stop("Unsupported innovation family: ", family)
  }
  list(m = m, Sigma_E = Sigma_E)
}


# Draws T iid innovation vectors (d x T) with the construction used in the moments above.
.minar_draw_innovations <- function(T, d, family, params) {
  if (family == "poisson") {
    lambda <- params$lambda
    Y0 <- rpois(T, lambda[1])
    E <- matrix(0L, nrow = d, ncol = T)
    for (i in seq_len(d)) {
      E[i, ] <- rpois(T, lambda[i + 1]) + Y0
    }
    E
  } else if (family == "negbin") {
    r    <- params$r
    pvec <- params$p
    p0   <- 1 - sum(pvec)
    Lambda <- rgamma(T, shape = r, rate = 1)
    E <- matrix(0L, nrow = d, ncol = T)
    for (i in seq_len(d)) {
      E[i, ] <- rpois(T, lambda = Lambda * pvec[i] / p0)
    }
    E
  } else {
    stop("Unsupported innovation family: ", family)
  }
}


# Companion (block) thinning matrix of size dp x dp: top block row [A_1, ..., A_p],
# identity blocks on the sub-diagonal.
.minar_companion_A <- function(A_list, d, p) {
  dp <- d * p
  A_star <- matrix(0, nrow = dp, ncol = dp)

  for (j in seq_len(p)) {
    cols <- ((j - 1) * d + 1):(j * d)
    A_star[1:d, cols] <- A_list[[j]]
  }

  if (p >= 2L) {
    for (b in seq_len(p - 1L)) {
      rows <- (b * d + 1):((b + 1) * d)
      cols <- ((b - 1) * d + 1):(b * d)
      A_star[rows, cols] <- diag(d)
    }
  }

  A_star
}


# Integer matrix power by repeated multiplication; M^0 is the identity.
.minar_mat_power <- function(M, h) {
  n <- nrow(M)
  if (h == 0L) return(diag(n))
  P <- M
  if (h > 1L) {
    for (k in seq_len(h - 1L)) P <- P %*% M
  }
  P
}



# Convolution of two pmf vectors indexed from 0.
.minar_conv2 <- function(p, q) {
  n <- length(p); m <- length(q)
  out <- numeric(n + m - 1L)
  for (i in seq_along(p)) {
    if (p[i] != 0) out[i:(i + m - 1L)] <- out[i:(i + m - 1L)] + p[i] * q
  }
  out
}


# pmf of sum_j Binom(x_prev[j], A_row[j]), indexed from 0. Its support is [0, sum(x_prev)],
# which bounds the likelihood's convolution sum.
.minar_thinning_pmf <- function(x_prev, A_row) {
  d <- length(x_prev)
  pmf <- dbinom(0:x_prev[1], x_prev[1], A_row[1])
  if (d >= 2L) {
    for (j in 2:d) {
      pmf <- .minar_conv2(pmf, dbinom(0:x_prev[j], x_prev[j], A_row[j]))
    }
  }
  pmf
}


# Joint innovation pmf evaluated for every row of the nrow x d grid EG at once.
# Poisson sums over the shared component; NegBin is closed form, evaluated in log space.
.minar_fE_grid <- function(EG, family, params) {
  n <- nrow(EG)
  d <- ncol(EG)

  if (family == "poisson") {
    lambda_0 <- params$lambda[1]
    lambda_i <- params$lambda[-1]

    row_min <- do.call(pmin, lapply(seq_len(d), function(i) EG[, i]))
    kmax <- max(row_min)

    total <- numeric(n)
    if (kmax >= 0) {
      for (k in 0:kmax) {
        valid <- row_min >= k
        if (!any(valid)) next
        diffs <- EG[valid, , drop = FALSE] - k
        prod_i <- Reduce(`*`, lapply(seq_len(d), function(i) dpois(diffs[, i], lambda_i[i])))
        total[valid] <- total[valid] + dpois(k, lambda_0) * prod_i
      }
    }
    total
  } else if (family == "negbin") {
    r    <- params$r
    pvec <- params$p
    p0   <- 1 - sum(pvec)
    se   <- rowSums(EG)
    lfact_sum <- rowSums(lfactorial(EG))
    log_pmf <- lgamma(r + se) - lgamma(r) - lfact_sum +
      r * log(p0) + as.vector(EG %*% log(pvec))
    exp(log_pmf)
  } else {
    stop("Unsupported innovation family: ", family)
  }
}


# Cartesian product of d integer ranges as an nrow x d matrix. Built with rep() because
# the likelihood calls this too often for expand.grid()'s data.frame overhead.
.minar_cartesian_grid <- function(ranges) {
  d    <- length(ranges)
  lens <- vapply(ranges, length, integer(1))
  n    <- prod(lens)

  EG <- matrix(0, nrow = n, ncol = d)
  rep_each <- 1L
  for (i in seq_len(d)) {
    EG[, i] <- rep(ranges[[i]], each = rep_each, length.out = n)
    rep_each <- rep_each * lens[i]
  }
  EG
}


# Per-time exact innovation ranges for each component. They depend only on x_data, so
# they are computed once per est() call. Only ranges are cached: the joint grid depends
# on theta through the mass filter in .minar_loglik().
# Entry [[t-1]] is list(e_ranges) for t in 2:Tn.
.minar_precompute_grids <- function(x_data) {
  d  <- nrow(x_data)
  Tn <- ncol(x_data)

  lapply(2:Tn, function(t) {
    x_prev <- x_data[, t - 1]
    x_t    <- x_data[, t]
    Stot   <- sum(x_prev)

    e_ranges <- lapply(seq_len(d), function(i) {
      lo <- max(0, x_t[i] - Stot)
      hi <- x_t[i]
      lo:hi
    })

    list(e_ranges = e_ranges)
  })
}


# Per-component probability mass left out when truncating the likelihood grid. This makes
# the likelihood an approximation of the exact convolution sum.
.MINAR_ALPHA_MASS_THRESHOLD <- 1e-4


# Conditional log-likelihood sum_t log f(x_t | x_{t-1}; theta), with mass-based grid truncation.
# Each component's range is cut to the smallest window around the pmf mode holding at least
# 1 - alpha_mass of the mass; the window grows until it covers the full range, so mass needed
# to explain the observation is never dropped. Truncating per component before forming the
# product keeps the joint grid small.
.minar_loglik <- function(A, family, params, x_data, grids, alpha_mass = .MINAR_ALPHA_MASS_THRESHOLD) {
  d  <- nrow(x_data)
  Tn <- ncol(x_data)
  ll <- 0
  target <- 1 - alpha_mass

  for (t in 2:Tn) {
    x_prev <- x_data[, t - 1]
    x_t    <- x_data[, t]
    e_ranges_exact <- grids[[t - 1L]]$e_ranges

    f_list <- lapply(seq_len(d), function(i) .minar_thinning_pmf(x_prev, A[i, ]))

    # Grow a window from the mode toward the neighbour with more mass.
    e_ranges_filtered <- lapply(seq_len(d), function(i) {
      rng <- e_ranges_exact[[i]]
      n   <- length(rng)
      if (n == 1L) return(rng)

      k     <- x_t[i] - rng
      probs <- f_list[[i]][k + 1L]

      mode_idx <- which.max(probs)
      lo_idx <- mode_idx; hi_idx <- mode_idx
      cum <- probs[mode_idx]
      while (cum < target && (lo_idx > 1L || hi_idx < n)) {
        left_val  <- if (lo_idx > 1L) probs[lo_idx - 1L] else -Inf
        right_val <- if (hi_idx < n) probs[hi_idx + 1L] else -Inf
        if (left_val >= right_val) {
          lo_idx <- lo_idx - 1L
          cum <- cum + left_val
        } else {
          hi_idx <- hi_idx + 1L
          cum <- cum + right_val
        }
      }
      rng[lo_idx:hi_idx]
    })

    EG <- .minar_cartesian_grid(e_ranges_filtered)

    K <- matrix(x_t, nrow = nrow(EG), ncol = d, byrow = TRUE) - EG
    prod_f <- Reduce(`*`, lapply(seq_len(d), function(i) f_list[[i]][K[, i] + 1L]))

    fE_vals <- .minar_fE_grid(EG, family, params)

    total <- sum(prod_f * fE_vals)
    # Clamp so log() stays finite for the optimizer.
    if (!is.finite(total) || total <= 0) total <- .Machine$double.xmin
    ll <- ll + log(total)
  }

  ll
}


# Parameter-vector layout for optim():
#   [1 .. d^2]: A, row-major
#   then Poisson: lambda_0, ..., lambda_d; NegBin: r, p_1, ..., p_d
.minar_flatten <- function(A, family, params) {
  vecA  <- as.vector(t(A))
  extra <- if (family == "poisson") params$lambda else c(params$r, params$p)
  c(vecA, extra)
}

.minar_unflatten <- function(theta_vec, d, family) {
  A    <- matrix(theta_vec[seq_len(d^2)], nrow = d, ncol = d, byrow = TRUE)
  rest <- theta_vec[(d^2 + 1L):length(theta_vec)]
  if (family == "poisson") {
    list(A = A, innovation_params = list(lambda = rest))
  } else {
    list(A = A, innovation_params = list(r = rest[1], p = rest[-1]))
  }
}


# Box constraints for L-BFGS-B over the flattened layout; eps is the interior margin.
# NegBin's sum(p_i) <= 1 - eps is not a box constraint and is enforced as a penalty in nll().
.minar_build_bounds <- function(d, family, eps) {
  lower <- c(rep(eps, d^2), rep(eps, d + 1L))
  upper <- if (family == "poisson") {
    c(rep(1 - eps, d^2), rep(Inf, d + 1L))
  } else {
    c(rep(1 - eps, d^2), Inf, rep(1 - eps, d))
  }
  list(lower = lower, upper = upper)
}


# Method-of-moments starting values from the lag-0 and lag-1 sample moments, clipped into the
# optimizer's box. Unclipped values can be invalid at small T.
.minar_mom_start <- function(x_data, family, eps) {
  d <- nrow(x_data)

  sm        <- sample_moments(x_data, lags = c(0, 1))
  mu_bar    <- sm$mu
  Sigma_bar <- sm$gamma[["0"]]
  Gamma1    <- sm$gamma[["1"]]

  A_hat <- Gamma1 %*% solve(Sigma_bar)
  m_hat <- as.vector((diag(d) - A_hat) %*% mu_bar)

  D_mat <- matrix(0, nrow = d, ncol = d)
  for (i in seq_len(d)) {
    D_mat[i, i] <- sum(mu_bar * A_hat[i, ] * (1 - A_hat[i, ]))
  }
  Sigma_E_hat <- Sigma_bar - A_hat %*% Sigma_bar %*% t(A_hat) - D_mat

  off_diag <- function(M) M[row(M) != col(M)]

  # Match m_hat and Sigma_E_hat to .minar_innovation_moments().
  if (family == "poisson") {
    lambda_0_hat <- mean(off_diag(Sigma_E_hat))
    lambda_i_hat <- m_hat - lambda_0_hat
    params_hat <- list(lambda = c(lambda_0_hat, lambda_i_hat))
  } else {
    R_mat <- matrix(0, nrow = d, ncol = d)
    for (i in seq_len(d)) for (j in seq_len(d)) {
      if (i != j) R_mat[i, j] <- (m_hat[i] * m_hat[j]) / Sigma_E_hat[i, j]
    }
    r_hat <- mean(off_diag(R_mat))
    denom <- r_hat + sum(m_hat)
    p_hat <- m_hat / denom
    params_hat <- list(r = r_hat, p = p_hat)
  }

  A_clipped <- pmin(pmax(A_hat, eps), 1 - eps)

  if (family == "poisson") {
    params_clipped <- list(lambda = pmax(params_hat$lambda, eps))
  } else {
    r_clipped <- max(params_hat$r, eps)
    p_clipped <- pmin(pmax(params_hat$p, eps), 1 - eps)
    if (sum(p_clipped) > 1 - eps) {
      # Rescaling may push a single p_i slightly below eps; acceptable for a starting value.
      p_clipped <- p_clipped * (1 - eps) / sum(p_clipped)
    }
    params_clipped <- list(r = r_clipped, p = p_clipped)
  }

  list(A = A_clipped, innovation_params = params_clipped)
}


# Constructor: slots are closures over the Markov order p and the innovation family.
make_minar_model <- function(name, p, family) {

  if (missing(name) || missing(p) || missing(family)) {
    stop("make_minar_model(): `name`, `p`, and `family` are all required.")
  }
  if (!(family %in% c("poisson", "negbin"))) {
    stop('make_minar_model(): `family` must be "poisson" or "negbin".')
  }

  # theta$A must hold exactly p matrices. Returns the T newly simulated columns.
  sim <- function(theta, x_init, T) {
    A  <- theta$A
    pp <- length(A)
    d  <- nrow(A[[1]])

    E <- .minar_draw_innovations(T, d, family, theta$innovation_params)

    # Buffer: pp history columns from x_init, then T new columns. NULL x_init starts from zeros.
    X <- matrix(0L, nrow = d, ncol = T + pp)
    if (!is.null(x_init)) X[, seq_len(pp)] <- x_init

    for (t in seq_len(T)) {
      col <- pp + t
      contrib <- numeric(d)

      for (j in seq_len(pp)) {
        lag_col <- X[, col - j]
        Aj <- A[[j]]

        # Binomial thinning of each component of the lagged vector, summed per row.
        for (i in seq_len(d)) {
          contrib[i] <- contrib[i] + sum(rbinom(d, size = lag_col, prob = Aj[i, ]))
        }
      }

      X[, col] <- contrib + E[, t]
    }

    X[, (pp + 1):(pp + T), drop = FALSE]
  }

  # Exact moments: write the model as a dp-dimensional MINAR(1) in companion form,
  # solve for the mean and the Lyapunov equation, and read off the top-left blocks.
  true_moments <- function(theta) {
    A  <- theta$A
    pp <- length(A)
    d  <- nrow(A[[1]])
    dp <- d * pp

    A_star <- .minar_companion_A(A, d, pp)

    # Innovation moments padded to the companion level (only the first block is nonzero).
    im <- .minar_innovation_moments(family, theta$innovation_params, d)
    m_star <- c(im$m, rep(0, d * (pp - 1L)))
    Sigma_E_star <- matrix(0, nrow = dp, ncol = dp)
    Sigma_E_star[1:d, 1:d] <- im$Sigma_E

    mu_star <- solve(diag(dp) - A_star, m_star)

    V_star <- Sigma_E_star
    for (i in seq_len(dp)) {
      V_star[i, i] <- V_star[i, i] +
        sum(mu_star * A_star[i, ] * (1 - A_star[i, ]))
    }

    # vec(Sigma) = (I - A %x% A)^{-1} vec(V)
    dp2 <- dp^2
    vecV <- as.vector(V_star)
    vecSigma <- solve(diag(dp2) - kronecker(A_star, A_star), vecV)
    Sigma_star <- matrix(vecSigma, nrow = dp, ncol = dp)

    # Gamma_star(h) = A_star^h Sigma_star
    lags <- 0:3
    gamma <- lapply(lags, function(h) {
      Gh_star <- .minar_mat_power(A_star, h) %*% Sigma_star
      Gh_star[1:d, 1:d, drop = FALSE]
    })
    names(gamma) <- as.character(lags)

    list(mu = mu_star[1:d], gamma = gamma)
  }

  # est() is only attached for p = 1.
  est_slot <- NULL
  if (p == 1L) {
    est_slot <- function(x_data) {
      start_time <- Sys.time()

      d   <- nrow(x_data)
      eps <- 1e-3   # interior margin shared by the starting-value clipping and the bounds

      start       <- .minar_mom_start(x_data, family, eps)
      theta_start <- .minar_flatten(start$A, family, start$innovation_params)
      bounds      <- .minar_build_bounds(d, family, eps)
      theta_start <- pmin(pmax(theta_start, bounds$lower), bounds$upper)

      # Computed once here because the exact ranges depend only on x_data.
      grids <- .minar_precompute_grids(x_data)

      nll <- function(theta_vec) {
        th <- .minar_unflatten(theta_vec, d, family)

        if (family == "negbin") {
          # Penalty for the sum(p_i) constraint, which is not a box constraint.
          p0 <- 1 - sum(th$innovation_params$p)
          if (p0 < eps) return(1e10)
        }

        -.minar_loglik(th$A, family, th$innovation_params, x_data, grids)
      }

      # L-BFGS-B supports box constraints and finite-difference gradients; constrOptim() needs
      # an analytic gradient. optim() must not be wrapped in its own tryCatch(): the inner
      # handler would turn a timeout into a plain error before with_est_timeout() sees it.
      # The outer tryCatch() handles what with_est_timeout() re-throws.
      timeout_seconds <- get_config()$est_timeout_seconds
      to <- tryCatch(
        with_est_timeout(
          optim(par = theta_start, fn = nll, gr = NULL,
                method = "L-BFGS-B",
                lower = bounds$lower, upper = bounds$upper,
                control = list(maxit = 1000, factr = 1e9, pgtol = 1e-5)),
          timeout_seconds
        ),
        error = function(e) e
      )

      if (inherits(to, "error")) {
        fit_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        return(list(theta = NULL, status = "failed", stationary = NA, fit_time = fit_time))
      }
      if (to$timed_out) {
        return(list(theta = NULL, status = "timed_out", stationary = NA, fit_time = timeout_seconds))
      }
      fit <- to$value

      # optim()$convergence: 0 converged, 1 maxit reached, anything else failed.
      status <- if (!is.null(fit$convergence) && fit$convergence == 0L) {
        "converged"
      } else if (!is.null(fit$convergence) && fit$convergence == 1L) {
        "hit_max_iter"
      } else {
        "failed"
      }

      if (status == "failed") {
        fit_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        return(list(theta = NULL, status = status, stationary = NA, fit_time = fit_time))
      }

      fitted       <- .minar_unflatten(fit$par, d, family)
      A_hat_fitted <- fitted$A

      # Stationary if the spectral radius of A is below 1.
      stationary <- max(Mod(eigen(A_hat_fitted, only.values = TRUE)$values)) < 1

      fit_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      list(
        theta = list(A = list(A_hat_fitted), innovation_params = fitted$innovation_params),
        status = status,
        stationary = stationary,
        fit_time = fit_time
      )
    }
  }

  new_model(
    name         = name,
    sim          = sim,
    p            = p,
    est          = est_slot,
    true_moments = true_moments
  )
}
