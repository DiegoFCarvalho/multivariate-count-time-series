# MDACP model: double Poisson margins with a Gaussian copula.
# Order p = 1, mean recursion mu_t = omega + A x_{t-1} + B mu_{t-1},
# A = gamma delta' + diag(alpha) (rank r), B = diag(beta).
# Needs mvtnorm, gamlss.dist, callr and Matrix.
# theta = list(omega, gamma, delta, alpha, beta, phi, Sigma, mu_init); series d x T, x_init d x 1.

# A and B from the reduced-rank factorization.
.mdacp_build_AB <- function(theta) {
  d <- length(theta$omega)
  A <- theta$gamma %*% t(theta$delta) + diag(theta$alpha, nrow = d)
  B <- diag(theta$beta, nrow = d)
  list(A = A, B = B)
}


# optim() layout: omega, gamma (row-major), free rows of delta (row-major), alpha, beta, phi.
# The top r rows of delta are fixed to the identity (identification).
.mdacp_flatten <- function(omega, gamma, delta, alpha, beta, phi) {
  d <- nrow(gamma)
  r <- ncol(gamma)
  vec_gamma <- as.vector(t(gamma))
  rows_free <- if (r < d) (r + 1L):d else integer(0)
  vec_delta_free <- as.vector(t(delta[rows_free, , drop = FALSE]))
  c(omega, vec_gamma, vec_delta_free, alpha, beta, phi)
}

.mdacp_unflatten <- function(theta_vec, d, r) {
  idx <- 0L
  take <- function(n) {
    out <- theta_vec[(idx + 1L):(idx + n)]
    idx <<- idx + n
    out
  }

  omega      <- take(d)
  gamma      <- matrix(take(d * r), nrow = d, ncol = r, byrow = TRUE)
  rows_free  <- if (r < d) (r + 1L):d else integer(0)
  delta_free <- take(length(rows_free) * r)
  delta      <- matrix(0, nrow = d, ncol = r)
  n_id       <- min(r, d)
  delta[seq_len(n_id), seq_len(n_id)] <- diag(n_id)
  if (length(rows_free) > 0L) {
    delta[rows_free, ] <- matrix(delta_free, nrow = length(rows_free), ncol = r, byrow = TRUE)
  }
  alpha <- take(d)
  beta  <- take(d)
  phi   <- take(d)

  list(omega = omega, gamma = gamma, delta = delta, alpha = alpha, beta = beta, phi = phi)
}


# L-BFGS-B bounds: non-negative parameters keep the conditional mean non-negative.
.mdacp_build_bounds <- function(d, r, eps) {
  n_free_delta <- max(d - r, 0L) * r
  lower <- c(rep(0, d), rep(0, d * r), rep(0, n_free_delta), rep(0, d), rep(0, d), rep(eps, d))
  upper <- c(rep(Inf, d), rep(Inf, d * r), rep(Inf, n_free_delta), rep(Inf, d), rep(Inf, d), rep(Inf, d))
  list(lower = lower, upper = upper)
}


# Conditional-mean recursion over an observed series; the first column conditions and is not scored.
# Returns mu0, mu_mat (d x (Tn-1)) and mu_T (last mean).
.mdacp_mu_sequence <- function(omega, A, B, x_data) {
  d  <- length(omega)
  Tn <- ncol(x_data)

  mu0 <- as.vector(solve(diag(d) - A - B, omega))

  mu_prev <- mu0
  x_prev  <- x_data[, 1]
  mu_mat  <- matrix(NA_real_, nrow = d, ncol = Tn - 1L)

  for (t in 2:Tn) {
    mu_t <- as.vector(omega + A %*% x_prev + B %*% mu_prev)
    mu_mat[, t - 1L] <- mu_t
    if (any(!is.finite(mu_t))) break
    x_prev  <- x_data[, t]
    mu_prev <- mu_t
  }

  list(mu0 = mu0, mu_mat = mu_mat, mu_T = mu_prev)
}


# Part A's optim() fit, run in a fresh child process so that a timeout can kill it.
.mdacp_optim_subprocess <- function(project_root, x_data, d, r, theta_start, bounds) {
  setwd(project_root)
  source("models/model_mdacp.R")

  Tn <- ncol(x_data)

  # Negative log-likelihood of the marginals; invalid trials return a large penalty.
  nll <- function(theta_vec) {
    th <- .mdacp_unflatten(theta_vec, d, r)
    AB <- .mdacp_build_AB(th)

    mu_seq <- tryCatch(
      .mdacp_mu_sequence(th$omega, AB$A, AB$B, x_data),
      error = function(e) NULL
    )
    if (is.null(mu_seq)) return(1e10)
    if (any(!is.finite(mu_seq$mu0)) || any(mu_seq$mu0 <= 0)) return(1e10)
    if (any(!is.finite(mu_seq$mu_mat)) || any(mu_seq$mu_mat <= 0)) return(1e10)

    x_scored  <- x_data[, 2:Tn, drop = FALSE]
    sigma_rep <- rep(1 / th$phi, times = Tn - 1L)

    ll <- sum(gamlss.dist::dDPO(as.vector(x_scored), mu = as.vector(mu_seq$mu_mat),
                                sigma = sigma_rep, log = TRUE))
    if (!is.finite(ll)) return(1e10)
    -ll
  }

  optim(par = theta_start, fn = nll, gr = NULL,
        method = "L-BFGS-B",
        lower = bounds$lower, upper = bounds$upper,
        control = list(maxit = 1000))
}


# Constructor; r is the rank of the gamma / delta factorization.
make_mdacp_model <- function(name, r = 1) {

  # Per step: update mu_t, draw correlated normal scores, invert the double Poisson cdf.
  sim <- function(theta, x_init, T) {
    omega <- theta$omega
    AB    <- .mdacp_build_AB(theta)
    A     <- AB$A
    B     <- AB$B
    phi   <- theta$phi
    Sigma <- theta$Sigma
    d     <- length(omega)

    mu_prev <- if (!is.null(theta$mu_init)) {
      as.vector(theta$mu_init)
    } else {
      as.vector(solve(diag(d) - A - B, omega))
    }

    x_prev <- if (is.null(x_init)) rep(0, d) else as.vector(x_init)

    X <- matrix(0L, nrow = d, ncol = T)

    for (t in seq_len(T)) {
      mu_t <- as.vector(omega + A %*% x_prev + B %*% mu_prev)

      # A diverging recursion is signalled with a dedicated condition class.
      if (any(!is.finite(mu_t)) || any(mu_t <= 0)) {
        stop(structure(
          class = c("mdacp_forecast_breakdown", "error", "condition"),
          list(
            message = sprintf(
              "model_mdacp.R sim() ('%s'): mu_t non-finite/non-positive at t=%d (mu_t = %s) -- recursion diverged mid-forecast, aborting this sim() call",
              name, t, paste(round(mu_t, 4), collapse = ", ")
            ),
            call = sys.call(-1L),
            t    = t,
            mu_t = mu_t
          )
        ))
      }

      Z_t <- as.vector(mvtnorm::rmvnorm(1, mean = rep(0, d), sigma = Sigma))

      U_tilde_t <- pnorm(Z_t)

      x_t <- vapply(seq_len(d), function(i) {
        gamlss.dist::qDPO(U_tilde_t[i], mu = mu_t[i], sigma = 1 / phi[i])
      }, numeric(1))

      X[, t] <- x_t

      x_prev  <- x_t
      mu_prev <- mu_t
    }

    X
  }

  est <- function(x_data) {
    start_time <- Sys.time()

    d   <- nrow(x_data)
    Tn  <- ncol(x_data)
    eps <- 1e-3   # lower bound margin for phi

    # Part A: marginal parameters (start values from sample moments)
    sample_mean <- rowMeans(x_data)
    sample_var  <- apply(x_data, 1, stats::var)

    gamma_start <- matrix(0.1, nrow = d, ncol = r)
    delta_start <- matrix(0, nrow = d, ncol = r)
    n_id_start  <- min(r, d)
    delta_start[seq_len(n_id_start), seq_len(n_id_start)] <- diag(n_id_start)
    rows_free   <- if (r < d) (r + 1L):d else integer(0)
    if (length(rows_free) > 0L) delta_start[rows_free, ] <- 0.1
    alpha_start <- rep(0, d)
    beta_start  <- rep(0.3, d)

    AB_start <- .mdacp_build_AB(list(
      omega = rep(0, d), gamma = gamma_start, delta = delta_start,
      alpha = alpha_start, beta = beta_start
    ))
    omega_start <- as.vector((diag(d) - AB_start$A - AB_start$B) %*% sample_mean)
    omega_start <- pmax(omega_start, eps)

    phi_start <- pmax(sample_mean / sample_var, eps)

    theta_start <- .mdacp_flatten(omega_start, gamma_start, delta_start, alpha_start, beta_start, phi_start)
    bounds      <- .mdacp_build_bounds(d, r, eps)

    # The fit runs in a callr subprocess; project_root is the working directory.
    timeout_seconds <- est_timeout_seconds()
    project_root    <- getwd()

    to <- tryCatch(
      list(timed_out = FALSE, value = callr::r(
        func    = .mdacp_optim_subprocess,
        args    = list(project_root = project_root, x_data = x_data, d = d, r = r,
                       theta_start = theta_start, bounds = bounds),
        timeout = timeout_seconds
      )),
      callr_timeout_error = function(e) list(timed_out = TRUE, value = NULL),
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

    status <- if (!is.null(fit$convergence) && fit$convergence == 0L) {
      "converged"
    } else {
      "hit_max_iter"
    }

    fitted    <- .mdacp_unflatten(fit$par, d, r)
    AB_hat    <- .mdacp_build_AB(fitted)
    A_hat     <- AB_hat$A
    B_hat     <- AB_hat$B
    omega_hat <- fitted$omega
    gamma_hat <- fitted$gamma
    delta_hat <- fitted$delta
    alpha_hat <- fitted$alpha
    beta_hat  <- fitted$beta
    phi_hat   <- fitted$phi

    # Part B: Gaussian copula correlation
    mu_seq_hat <- tryCatch(
      .mdacp_mu_sequence(omega_hat, A_hat, B_hat, x_data),
      error = function(e) NULL
    )
    if (is.null(mu_seq_hat) ||
        any(!is.finite(mu_seq_hat$mu_mat)) || any(mu_seq_hat$mu_mat <= 0)) {
      fit_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      return(list(theta = NULL, status = "failed", stationary = NA, fit_time = fit_time))
    }

    mu_mat_hat    <- mu_seq_hat$mu_mat
    x_scored      <- x_data[, 2:Tn, drop = FALSE]
    sigma_hat_rep <- rep(1 / phi_hat, times = Tn - 1L)

    # Jittered uniforms u = F(x-1) + f(x) U.
    U       <- matrix(runif(d * (Tn - 1L)), nrow = d, ncol = Tn - 1L)
    F_lower <- gamlss.dist::pDPO(as.vector(x_scored) - 1, mu = as.vector(mu_mat_hat),
                                 sigma = sigma_hat_rep)
    f_at_x  <- gamlss.dist::dDPO(as.vector(x_scored), mu = as.vector(mu_mat_hat),
                                 sigma = sigma_hat_rep)
    u_tilde <- F_lower + f_at_x * as.vector(U)

    # Normal scores, second-moment matrix, unit diagonal.
    z          <- qnorm(u_tilde)
    z_mat      <- matrix(z, nrow = d, ncol = Tn - 1L)
    Sigma_hat  <- (z_mat %*% t(z_mat)) / (Tn - 1L)

    Sigma_hat_norm <- Sigma_hat / sqrt(outer(diag(Sigma_hat), diag(Sigma_hat)))

    # Nearest valid correlation matrix.
    Sigma_hat_norm <- Matrix::nearPD(Sigma_hat_norm, corr = TRUE, base.matrix = TRUE)$mat

    # Stationary if the spectral radius of A + B is below 1.
    stationary <- max(Mod(eigen(A_hat + B_hat, only.values = TRUE)$values)) < 1

    fit_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

    list(
      theta = list(
        omega   = omega_hat,
        gamma   = gamma_hat,
        delta   = delta_hat,
        alpha   = alpha_hat,
        beta    = beta_hat,
        phi     = phi_hat,
        Sigma   = Sigma_hat_norm,
        # last fitted mean, seed for forecasts
        mu_init = mu_seq_hat$mu_T
      ),
      status     = status,
      stationary = stationary,
      fit_time   = fit_time
    )
  }

  new_model(name = name, sim = sim, p = 1L, est = est)
}
