# MDACP model (multivariate autoregressive conditional double Poisson), Gaussian copula.
# est() uses inference functions for margins: Part A fits the marginal parameters,
# Part B estimates the copula correlation from Part A's jittered residuals.
# Fixed structure: p = 1, and A = gamma %*% t(delta) + Diag(alpha), B = Diag(beta).
# gamlss.dist's DPO functions take sigma = 1/phi, so 1/phi is passed wherever phi appears.
# Requires mvtnorm, gamlss.dist, callr and Matrix.
# theta = list(omega (d), gamma (d x r), delta (d x r), alpha (d), beta (d), phi (d),
#              Sigma (d x d correlation), mu_init (d, optional)).
# The top r x r block of delta must be the identity; .mdacp_build_AB() does not enforce it.
# mu_init is the conditional-mean state preceding x_init; if absent, sim() starts from the
# unconditional mean.
# Series are d x T; x_init is d x 1.


# Integer matrix power by repeated multiplication; M^0 is the identity.
.mdacp_mat_power <- function(M, h) {
  n <- nrow(M)
  if (h == 0L) return(diag(n))
  P <- M
  if (h > 1L) {
    for (k in seq_len(h - 1L)) P <- P %*% M
  }
  P
}


# Rebuilds the d x d matrices A and B from the reduced-rank factorization.
# diag(x, nrow = d) avoids R's scalar-as-size ambiguity when d = 1.
.mdacp_build_AB <- function(theta) {
  d <- length(theta$omega)
  A <- theta$gamma %*% t(theta$delta) + diag(theta$alpha, nrow = d)
  B <- diag(theta$beta, nrow = d)
  list(A = A, B = B)
}



# Parameter-vector layout for optim():
#   omega (d), gamma (row-major, d x r), delta[(r+1):d, ] (row-major), alpha (d), beta (d), phi (d)
# gamma %*% t(delta) is unchanged under gamma -> gamma M, delta -> delta solve(t(M)) for any
# invertible r x r M. Fixing the top r rows of delta to the identity removes exactly that
# redundancy. If r >= d, delta[1:d, 1:d] is fixed to the identity and no delta entries are free.
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


# Box constraints for L-BFGS-B over the flattened layout. Non-negative parameters give
# A, B >= 0 entrywise, which keeps the conditional mean non-negative. This rules out negative
# entries in A; negative cross-sectional dependence remains possible through the copula.
.mdacp_build_bounds <- function(d, r, eps) {
  n_free_delta <- max(d - r, 0L) * r
  lower <- c(rep(0, d), rep(0, d * r), rep(0, n_free_delta), rep(0, d), rep(0, d), rep(eps, d))
  upper <- c(rep(Inf, d), rep(Inf, d * r), rep(Inf, n_free_delta), rep(Inf, d), rep(Inf, d), rep(Inf, d))
  list(lower = lower, upper = upper)
}


# Conditional-mean recursion over an observed series. The first column of x_data is a fixed
# conditioning observation and is not scored; the recursion starts at the unconditional mean
# for the current (omega, A, B).
# Returns mu0 (starting mean), mu_mat (d x (Tn-1), aligned with x_data[, 2:Tn]) and mu_T
# (last mean, the seed for forecasting from the end of the data).
# solve() may error for singular (I - A - B); a non-finite mu_t stops the loop and leaves the
# remaining columns NA. Callers check for both.
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


# Part A's optim() fit, executed as the callr::r() function in a fresh child process. The
# parent kills the child on timeout, which works even while the child is inside a single
# long compiled call (dDPO), where setTimeLimit() cannot interrupt. The child shares no
# state with the parent, so it changes to the project root and re-sources this file.
# The fit is deterministic; the random jittering of Part B stays in the parent.
.mdacp_optim_subprocess <- function(project_root, x_data, d, r, theta_start, bounds) {
  setwd(project_root)
  source("models/model_mdacp.R")

  Tn <- ncol(x_data)

  # Negative log-likelihood of the marginal model. Invalid trials return a large penalty
  # instead of erroring inside dDPO().
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


# Constructor. p is always 1; r is the rank of the gamma/delta factorization and is only used
# by est(); sim() and true_moments() read it from theta$gamma.
make_mdacp_model <- function(name, r = 1) {

  if (missing(name)) {
    stop("make_mdacp_model(): `name` is required.")
  }

  # Recursion per step: update mu_t, draw correlated normal scores, map to uniforms, and
  # invert the double Poisson cdf (gamlss.dist's qDPO is the generalized inverse).
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

    # NULL x_init starts from a zero observation (mu_prev above is the mean-based state).
    x_prev <- if (is.null(x_init)) rep(0, d) else as.vector(x_init)

    X <- matrix(0L, nrow = d, ncol = T)

    for (t in seq_len(T)) {
      mu_t <- as.vector(omega + A %*% x_prev + B %*% mu_prev)

      # The recursion can diverge mid-forecast for a valid fitted theta. Signal it with a
      # dedicated condition class so the caller can absorb this case only, instead of letting
      # qDPO() fail with an opaque error.
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

  # Approximate closed-form moments. Two approximations: the double Poisson mean equals mu
  # only approximately (exact at phi = 1), and the copula-induced conditional covariance has
  # no closed form, so Omega below is an asymptotic approximation.
  true_moments <- function(theta) {
    omega <- theta$omega
    AB    <- .mdacp_build_AB(theta)
    A     <- AB$A
    B     <- AB$B
    phi   <- theta$phi
    Sigma <- theta$Sigma
    d     <- length(omega)

    # mu = (I - A - B)^{-1} omega
    mu <- as.vector(solve(diag(d) - A - B, omega))

    # marginal variances mu / phi
    sigma2 <- mu / phi

    V_DP <- diag(sigma2, nrow = d)

    # Omega ~= V^{1/2} Sigma V^{1/2}
    V_DP_sqrt <- diag(sqrt(sigma2), nrow = d)
    Omega <- V_DP_sqrt %*% Sigma %*% V_DP_sqrt

    # Lyapunov equation Sigma_mu = (A+B) Sigma_mu (A+B)' + A Omega A', solved exactly by
    # vectorization.
    AB2 <- A + B
    d2  <- d^2
    AOA <- A %*% Omega %*% t(A)
    vec_Sigma_mu <- solve(diag(d2) - kronecker(AB2, AB2), as.vector(AOA))
    Sigma_mu <- matrix(vec_Sigma_mu, nrow = d, ncol = d)

    # lag-0 autocovariance
    Sigma_X <- Omega + Sigma_mu

    # Gamma(h) = (A+B)^h Sigma_X - (A+B)^(h-1) B Omega for h >= 1
    lags <- 0:3
    gamma_out <- lapply(lags, function(h) {
      if (h == 0L) return(Sigma_X)
      .mdacp_mat_power(AB2, h) %*% Sigma_X -
        .mdacp_mat_power(AB2, h - 1L) %*% B %*% Omega
    })
    names(gamma_out) <- as.character(lags)

    list(mu = mu, gamma = gamma_out)
  }

  est <- function(x_data) {
    start_time <- Sys.time()

    d   <- nrow(x_data)
    Tn  <- ncol(x_data)
    eps <- 1e-3   # lower bound margin for phi

    # Part A: marginal parameters

    # omega_start makes the implied stationary mean equal the sample mean;
    # phi_start comes from Var ~= mu / phi.
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

    # The fit runs in a callr subprocess so that a timeout can kill it at the OS level (see
    # .mdacp_optim_subprocess()). Do not copy this to a model whose fit setTimeLimit() can
    # already interrupt. A callr timeout has its own condition class; callr_timeout_error is
    # listed first because it also inherits from "error". project_root is the current working
    # directory, which is the project root in every calling context.
    timeout_seconds <- get_config()$est_timeout_seconds
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

    # optim()$convergence: 0 is converged, anything else is hit_max_iter.
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
    # Recompute the mean path at the fitted parameters.
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

    # Jittered uniforms u = F(x-1) + f(x) * U with U ~ Unif(0,1).
    U       <- matrix(runif(d * (Tn - 1L)), nrow = d, ncol = Tn - 1L)
    F_lower <- gamlss.dist::pDPO(as.vector(x_scored) - 1, mu = as.vector(mu_mat_hat),
                                 sigma = sigma_hat_rep)
    f_at_x  <- gamlss.dist::dDPO(as.vector(x_scored), mu = as.vector(mu_mat_hat),
                                 sigma = sigma_hat_rep)
    u_tilde <- F_lower + f_at_x * as.vector(U)

    # Normal scores, their second-moment matrix, then normalization to unit diagonal.
    z          <- qnorm(u_tilde)
    z_mat      <- matrix(z, nrow = d, ncol = Tn - 1L)
    Sigma_hat  <- (z_mat %*% t(z_mat)) / (Tn - 1L)

    Sigma_hat_norm <- Sigma_hat / sqrt(outer(diag(Sigma_hat), diag(Sigma_hat)))

    # Project to the nearest valid correlation matrix in case of non-PSD sampling noise.
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
        # Last fitted mean: the seed for forecasting from the end of the data.
        mu_init = mu_seq_hat$mu_T
      ),
      status     = status,
      stationary = stationary,
      fit_time   = fit_time
    )
  }

  new_model(
    name         = name,
    sim          = sim,
    p            = 1L,
    est          = est,
    true_moments = true_moments
  )
}
