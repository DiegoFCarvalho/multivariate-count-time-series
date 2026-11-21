# Moment-based goodness-of-fit statistics for one path; the caller pools across replicates.

# Sample mean and autocovariances of one d x T path.
# gamma[[h]][i, j] = (1/T) * sum_t (x[i, t+h] - mu_i) (x[j, t] - mu_j), via stats::acf().
sample_moments <- function(x, lags = 0:3) {
  d  <- nrow(x)
  mu <- rowMeans(x)

  ac <- stats::acf(t(x), lag.max = max(lags), type = "covariance",
                   demean = TRUE, plot = FALSE)

  # matrix() keeps a d x d matrix when d = 1
  gamma <- lapply(lags, function(h) matrix(ac$acf[h + 1, , ], nrow = d, ncol = d))
  names(gamma) <- as.character(lags)

  list(mu = mu, gamma = gamma)
}


# GoF1: squared Euclidean distance between mean vectors; NA if mu_true is NULL.
mean_bias <- function(mu_true, mu_hat) {
  if (is.null(mu_true)) return(NA_real_)
  diff <- mu_true - mu_hat
  sum(diff * diff)
}


# GoF2: squared Frobenius distance per lag between two lists of d x d matrices.
# Lags with a missing true matrix, or all lags if gamma_true is NULL, give NA.
gamma_bias <- function(gamma_true, gamma_hat) {
  lag_names <- names(gamma_hat)

  if (is.null(gamma_true)) {
    return(stats::setNames(rep(NA_real_, length(lag_names)), lag_names))
  }

  vapply(lag_names, function(h) {
    G_true <- gamma_true[[h]]
    G_hat  <- gamma_hat[[h]]
    if (is.null(G_true) || is.null(G_hat)) return(NA_real_)
    D <- G_true - G_hat
    sum(D * D)
  }, numeric(1))
}


# Sample rank-correlation matrices of one d x T path, one per lag.
# Entry [i, j] is cor(x_i,t+h, x_j,t), with the same orientation as sample_moments().
sample_rank_corr <- function(x, lags = 0:3, method = "spearman") {
  d   <- nrow(x)
  n_t <- ncol(x)

  rc <- lapply(lags, function(h) {
    later   <- x[, (h + 1):n_t, drop = FALSE]
    earlier <- x[, 1:(n_t - h), drop = FALSE]
    matrix(stats::cor(t(later), t(earlier), method = method), nrow = d, ncol = d)
  })
  names(rc) <- as.character(lags)

  rc
}


# GoF3: same computation as gamma_bias(), applied to rank-correlation matrices.
rank_corr_bias <- function(true_corr, hat_corr) {
  gamma_bias(true_corr, hat_corr)
}
