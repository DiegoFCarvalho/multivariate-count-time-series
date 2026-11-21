# Multivariate count autoregression with a copula for cross-sectional dependence.
# DGP only: sim() and true_moments(), no est().
# Fixed structure: p = 1, linear conditional mean, copula family fixed at construction.
# Requires the copula package.
# theta = list(omega (d, > 0), A (d x d, > 0), B (d x d, > 0), copula_param):
#   (X_{i,t} | F_{t-1}) ~ Poisson(lambda_{i,t}), lambda_t = omega + A lambda_{t-1} + B X_{t-1}.
# Positive omega, A, B keep the Poisson means positive and are part of what a valid theta is.
# copula_param by family:
#   gaussian: (rho_12, rho_13, rho_23), filling the lower triangle column-major, as
#             copula:::p2P() does; must give a positive definite correlation matrix.
#   clayton:  scalar alpha > 0.
#   gumbel:   scalar alpha >= 1.
# Series are d x T; x_init is d x 1 and is the initial intensity lambda_0, not a count.
# Stationarity requires max(colSums(abs(A))) + max(colSums(abs(B))) < 1. This is stronger than
# the spectral-radius checks of the other models and is not called in this file.


# Integer matrix power by repeated multiplication; M^0 is the identity.
.mcar_mat_power <- function(M, h) {
  n <- nrow(M)
  if (h == 0L) return(diag(n))
  P <- M
  if (h > 1L) {
    for (k in seq_len(h - 1L)) P <- P %*% M
  }
  P
}


# Builds the copula object. dispstr = "un" gives a free correlation per pair for the Gaussian
# family. Neither normalCopula() nor rCopula() rejects an invalid correlation matrix (rCopula()
# only warns), so positive definiteness is checked here via the eigenvalues.
.mcar_build_copula <- function(copula_family, copula_param, d) {
  if (copula_family == "gaussian") {
    P <- diag(0, nrow = d)
    P[lower.tri(P)] <- copula_param
    P <- P + t(P)
    diag(P) <- 1

    eigenvalues <- eigen(P, symmetric = TRUE, only.values = TRUE)$values
    if (any(eigenvalues <= 1e-8)) {
      stop("mcar copula_param does not give a positive semi-definite ",
           "correlation matrix (smallest eigenvalue = ", min(eigenvalues), ")")
    }
  }

  switch(
    copula_family,
    gaussian = copula::normalCopula(param = copula_param, dim = d, dispstr = "un"),
    clayton  = copula::claytonCopula(param = copula_param, dim = d),
    gumbel   = copula::gumbelCopula(param = copula_param, dim = d),
    stop("Unsupported copula_family: ", copula_family,
        ' (supported: "gaussian", "clayton", "gumbel")')
  )
}


# Draws n_draws independent count vectors given the same lambda (n_draws x d integer matrix).
# Waiting-time construction: with copula rows U, Y = -log(U) / lambda and
# X_i = max{n : sum of the first n waiting times <= 1}, 0 if the first already exceeds 1.
# If a replication's cumulative sum does not cross 1 within the batch, the batch size is
# doubled and fresh draws are made for the unresolved replications only; this is valid
# because the draws are iid, and avoids truncating large counts.
# All replications of a round share one rCopula() call and the cumsum/crossing step runs as two
# apply() calls over arrays, which removes per-replication call overhead.
.mcar_draw_counts_batch <- function(lambda, cop, d, n_draws = 1L) {
  K_per_rep <- max(10L, ceiling(5 * max(lambda)))

  X <- matrix(NA_integer_, nrow = n_draws, ncol = d)
  unresolved <- seq_len(n_draws)

  while (length(unresolved) > 0L) {
    n_unresolved <- length(unresolved)
    U <- copula::rCopula(n_unresolved * K_per_rep, cop)
    Y <- -log(U) / matrix(lambda, nrow = nrow(U), ncol = d, byrow = TRUE)

    # Rows of Y are consecutive K_per_rep-row blocks, one per replication, which matches
    # array()'s column-major fill.
    Y_arr <- array(Y, dim = c(K_per_rep, n_unresolved, d))
    S_arr <- apply(Y_arr, c(2, 3), cumsum)
    hits  <- apply(S_arr > 1, c(2, 3), function(col) which(col)[1])

    resolved <- rowSums(is.na(hits)) == 0L
    if (any(resolved)) {
      X[unresolved[resolved], ] <- hits[resolved, , drop = FALSE] - 1L
    }

    unresolved <- which(is.na(X[, 1]))
    if (length(unresolved) > 0L) K_per_rep <- K_per_rep * 2L
  }

  X
}


# Same recursion as sim(), but returns the d x T matrix of lambda values used to generate each
# count. sim()'s return type cannot carry this, and true_moments() needs it.
.mcar_lambda_path <- function(theta, cop, lambda0, T) {
  omega <- theta$omega
  A     <- theta$A
  B     <- theta$B
  d     <- length(omega)

  lambda_prev <- lambda0
  Lambda_used <- matrix(0, nrow = d, ncol = T)

  for (t in seq_len(T)) {
    Lambda_used[, t] <- lambda_prev
    x_t <- as.vector(.mcar_draw_counts_batch(lambda_prev, cop, d, n_draws = 1L))
    lambda_prev <- as.vector(omega + A %*% lambda_prev + B %*% x_t)
  }

  Lambda_used
}


# Constructor. p is always 1; copula_family is "gaussian", "clayton" or "gumbel".
make_mcar_model <- function(name, copula_family = "gaussian") {

  if (missing(name)) {
    stop("make_mcar_model(): `name` is required.")
  }

  # Each iteration draws the count from the current lambda first and then updates lambda.
  # This is the reverse of model_mdacp.R, which updates mu before drawing.
  sim <- function(theta, x_init, T) {
    omega <- theta$omega
    A     <- theta$A
    B     <- theta$B
    d     <- length(omega)

    cop <- .mcar_build_copula(copula_family, theta$copula_param, d)

    # NULL x_init starts from the stationary mean intensity.
    lambda_prev <- if (is.null(x_init)) {
      as.vector(solve(diag(d) - A - B, omega))
    } else {
      as.vector(x_init)
    }

    X <- matrix(0L, nrow = d, ncol = T)

    for (t in seq_len(T)) {
      x_t <- as.vector(.mcar_draw_counts_batch(lambda_prev, cop, d, n_draws = 1L))
      X[, t] <- x_t

      lambda_prev <- as.vector(omega + A %*% lambda_prev + B %*% x_t)
    }

    X
  }

  # Mean and diagonal of Omega are exact; the off-diagonal of Omega is Monte Carlo estimated.
  true_moments <- function(theta) {
    omega <- theta$omega
    A     <- theta$A
    B     <- theta$B
    d     <- length(omega)

    # mu = (I - A - B)^{-1} omega. The stationarity condition of this model implies the
    # spectral-radius condition this formula needs.
    mu <- as.vector(solve(diag(d) - A - B, omega))

    # Omega = E[Cov(X_t | F_{t-1})]. Its diagonal equals mu for Poisson margins. The
    # off-diagonal depends on the copula and has no closed form: along a burnt-in lambda path,
    # draw K_omega fresh independent count vectors at each kept lambda_t, take their sample
    # covariance, and average over t. Fresh draws are needed because the recursion produces only
    # one count per t. The burn-in is implemented locally since model files must not depend on
    # run_simulation.R.
    cop <- .mcar_build_copula(copula_family, theta$copula_param, d)

    T_gt         <- 2000L
    T_burnin_gt  <- 200L
    K_omega      <- 200L

    lambda_path <- .mcar_lambda_path(theta, cop, lambda0 = mu, T = T_gt + T_burnin_gt)
    lambda_path <- lambda_path[, (T_burnin_gt + 1L):(T_gt + T_burnin_gt), drop = FALSE]

    cond_cov_sum <- matrix(0, nrow = d, ncol = d)
    for (t in seq_len(T_gt)) {
      X_sample <- .mcar_draw_counts_batch(lambda_path[, t], cop, d, n_draws = K_omega)
      cond_cov_sum <- cond_cov_sum + stats::cov(X_sample)
    }
    Omega_hat <- cond_cov_sum / T_gt
    diag(Omega_hat) <- mu

    # Lyapunov equation solved exactly by vectorization:
    #   Sigma_lambda = (A+B) Sigma_lambda (A+B)' + B Omega B',  Sigma_X = Omega + Sigma_lambda.
    # B (not A) multiplies Omega here, the reverse of model_mdacp.R.
    AB <- A + B
    d2 <- d^2
    BOB <- B %*% Omega_hat %*% t(B)
    vec_Sigma_lambda <- solve(diag(d2) - kronecker(AB, AB), as.vector(BOB))
    Sigma_lambda <- matrix(vec_Sigma_lambda, nrow = d, ncol = d)
    Sigma_X <- Omega_hat + Sigma_lambda

    # Gamma(h) = (A+B)^h Sigma_X - (A+B)^(h-1) A Omega for h >= 1. A (not B) appears in the
    # subtracted term; this differs from the Sigma_X formula above on purpose.
    lags <- 0:3
    gamma_out <- lapply(lags, function(h) {
      if (h == 0L) return(Sigma_X)
      .mcar_mat_power(AB, h) %*% Sigma_X -
        .mcar_mat_power(AB, h - 1L) %*% A %*% Omega_hat
    })
    names(gamma_out) <- as.character(lags)

    list(mu = mu, gamma = gamma_out)
  }

  new_model(
    name         = name,
    sim          = sim,
    p            = 1L,
    true_moments = true_moments
  )
}
