# Prints the true moments (mu, Gamma(h)) of every DGP; run from simulation_study/.
# Reads config.R, models/ and run_simulation.R; writes CSV and PNG tables to results/tables/.

source("config.R")
source("models/model_registry.R")

# Reuse run_simulation.R's theta_true_* objects: evaluate its top-level expressions only up to the
# `est_mets <- ...` assignment, so the driver loop never runs.
.load_true_theta <- function() {
  exprs <- parse(file = "run_simulation.R")

  is_est_mets_assign <- function(e) {
    is.call(e) && identical(as.character(e[[1]]), "<-") &&
      identical(as.character(e[[2]]), "est_mets")
  }
  cutoff <- Position(is_est_mets_assign, exprs, right = TRUE)
  if (is.na(cutoff)) {
    stop("True_moments.R: run_simulation.R's `est_mets <- ...` setup ",
         "assignment was not found -- its structure may have changed; ",
         "update this script's cutoff logic to match.")
  }

  env <- new.env(parent = globalenv())
  for (i in seq_len(cutoff)) eval(exprs[[i]], envir = env)

  # run_simulation.R's setup starts a `future` worker pool; stop it, this script runs sequentially.
  future::plan(future::sequential)

  env
}

run_sim_env <- .load_true_theta()
DGPs  <- run_sim_env$DGPs
Theta <- run_sim_env$Theta

# Keep each DGP's moments for the output files written at the end.
moments_record <- list()
record_moments <- function(f, nm) {
  force(f); force(nm)
  function(theta) {
    m <- f(theta)
    moments_record[[nm]] <<- m
    m
  }
}
for (nm in names(DGPs)) DGPs[[nm]]$true_moments <- record_moments(DGPs[[nm]]$true_moments, nm)

# gamma's lags come from the returned list's names; each model's true_moments() sets them.
for (dgp_name in names(DGPs)) {
  cat("\n==================== ", dgp_name, " ====================\n", sep = "")

  moments <- DGPs[[dgp_name]]$true_moments(Theta[[dgp_name]])

  cat("mu:\n")
  if (is.null(moments$mu)) {
    cat("  NULL (no closed-form mean for this model)\n")
  } else {
    print(moments$mu)
  }

  cat("\ngamma:\n")
  if (is.null(moments$gamma)) {
    cat("  no closed-form autocovariance for this model\n")
  } else {
    for (lag_name in names(moments$gamma)) {
      cat("  Gamma(", lag_name, "):\n", sep = "")
      print(moments$gamma[[lag_name]])
    }
  }
}

# svine1 has no closed-form autocovariance: GoF2's truth is the mean sample autocovariance over
# simulated paths (x_init = NULL; burn-in only guards start-up effects).
# n_rep = 1000 paths, T_path = 500, burnin = 200; margins Poisson(5), Poisson(7), Poisson(9).
# Centered on the known margin means (exact by Sklar's theorem): the sample mean would bias the
# variance down under positive autocorrelation.
svine1_model <- DGPs[["svine1"]]
svine1_theta <- Theta[["svine1"]]

n_rep_svine1  <- 1000L
T_path_svine1 <- 500L
burnin_svine1 <- 200L
T_total_svine1 <- T_path_svine1 + burnin_svine1
lags_svine1   <- 0:3

true_mu_svine1 <- c(5, 7, 9)   # Poisson(5), Poisson(7), Poisson(9)

d_svine1 <- length(svine1_theta$margins)   # one margin per cross-sectional variable

gamma_sum_svine1 <- lapply(lags_svine1, function(h) matrix(0, nrow = d_svine1, ncol = d_svine1))
names(gamma_sum_svine1) <- as.character(lags_svine1)

set.seed(2026L)   # fixed seed -- this estimate is meant to be read off once, reproducibly
mean_sum_svine1 <- rep(0, d_svine1)
for (r in seq_len(n_rep_svine1)) {
  x_full <- svine1_model$sim(svine1_theta, NULL, T = T_total_svine1)
  x  <- x_full[, (burnin_svine1 + 1L):T_total_svine1, drop = FALSE]
  mean_sum_svine1 <- mean_sum_svine1 + rowMeans(x)
  xc <- x - true_mu_svine1   # d x T_path, centered on the true (not sample) mean

  for (h in lags_svine1) {
    # Gamma(h)[i,j] = (1/T_path) * sum_{t=1}^{T_path-h} xc[i,t+h] * xc[j,t]
    g_h <- if (h == 0L) {
      (xc %*% t(xc)) / T_path_svine1
    } else {
      (xc[, (h + 1L):T_path_svine1, drop = FALSE] %*%
         t(xc[, 1:(T_path_svine1 - h), drop = FALSE])) / T_path_svine1
    }
    gamma_sum_svine1[[as.character(h)]] <- gamma_sum_svine1[[as.character(h)]] + g_h
  }
}
print(mean_sum_svine1 / n_rep_svine1)   # should be about 5, 7, 9
gamma_est_svine1 <- lapply(gamma_sum_svine1, `/`, n_rep_svine1)

cat("\n==================== svine1 (estimated, GoF2 ground truth) ====================\n")
cat("gamma (Monte-Carlo estimate, averaged over ", n_rep_svine1,
    " replicate paths, T = ", T_path_svine1, ", burnin = ", burnin_svine1, "):\n", sep = "")
for (lag_name in names(gamma_est_svine1)) {
  cat("  Gamma(", lag_name, "):\n", sep = "")
  print(gamma_est_svine1[[lag_name]])
}

# Output files: CSVs with full precision, PNG tables rounded as in the thesis.
out_dir <- file.path("results", "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_write_csv <- function(df, file) {
  df$value <- sprintf("%.17g", df$value)
  write.csv(df, file.path(out_dir, file), row.names = FALSE, quote = FALSE)
}

out_png <- function(df, file, width, height) {
  grob <- gridExtra::tableGrob(df, rows = NULL)
  ggplot2::ggsave(file.path(out_dir, file), plot = grob, width = width, height = height,
                  units = "in", dpi = 300)
}

out_fmt <- function(x, digits) formatC(round(x, digits), format = "f", digits = digits)

out_matrix_long <- function(m, dgp, quantity, method, lag) {
  m <- as.matrix(m)
  data.frame(dgp = dgp, quantity = quantity, method = method, lag = lag,
             i = rep(seq_len(nrow(m)), ncol(m)), j = rep(seq_len(ncol(m)), each = nrow(m)),
             value = as.vector(m))
}

out_gamma_method <- c(minar1_poi = "analytic", minar1_negbin = "analytic", mdacp1 = "analytic_approximation",
                      mcar1 = "analytic_with_monte_carlo_omega", svine1 = "monte_carlo")

# True moments of every DGP; svine1 autocovariances are the Monte Carlo estimate above.
out_moments <- list()
for (out_name in names(moments_record)) {
  out_m <- moments_record[[out_name]]
  out_moments[[length(out_moments) + 1L]] <-
    out_matrix_long(matrix(out_m$mu, nrow = 1L), out_name, "mu", "analytic", NA_integer_)
  out_g <- if (out_name == "svine1") gamma_est_svine1 else out_m$gamma
  for (out_lag in names(out_g)) {
    out_moments[[length(out_moments) + 1L]] <-
      out_matrix_long(out_g[[out_lag]], out_name, "gamma", out_gamma_method[[out_name]], as.integer(out_lag))
  }
}
out_moments[[length(out_moments) + 1L]] <-
  out_matrix_long(matrix(mean_sum_svine1 / n_rep_svine1, nrow = 1L), "svine1", "mu_monte_carlo", "monte_carlo", NA_integer_)
out_moments <- do.call(rbind, out_moments)
out_moments$i[out_moments$quantity %in% c("mu", "mu_monte_carlo")] <- NA_integer_
out_write_csv(out_moments, "true_moments.csv")

out_wide_row <- function(dgp, quantity, lag, i) {
  s <- out_moments[out_moments$dgp == dgp & out_moments$quantity == quantity &
                     (is.na(lag) | out_moments$lag == lag) & (is.na(i) | out_moments$i == i), ]
  s <- s[order(s$j), ]
  data.frame(DGP = dgp, Quantity = quantity, Lag = ifelse(is.na(lag), "", as.character(lag)),
             i = ifelse(is.na(i), "", as.character(i)), `j = 1` = out_fmt(s$value[1], 2),
             `j = 2` = out_fmt(s$value[2], 2), `j = 3` = out_fmt(s$value[3], 2), check.names = FALSE)
}
out_png_wide <- list()
for (out_name in names(DGPs)) {
  out_png_wide[[length(out_png_wide) + 1L]] <- out_wide_row(out_name, "mu", NA_integer_, NA_integer_)
  for (out_lag in 0:3) {
    for (out_i in 1:3) {
      out_png_wide[[length(out_png_wide) + 1L]] <- out_wide_row(out_name, "gamma", out_lag, out_i)
    }
  }
}
out_png_wide <- do.call(rbind, out_png_wide)
out_png(out_png_wide, "true_moments.png", width = 6, height = 0.23 * nrow(out_png_wide) + 0.6)

# True parameters of every DGP.
out_param <- list()
out_add <- function(dgp, parameter, m, label = "") {
  m <- as.matrix(m)
  out_param[[length(out_param) + 1L]] <<- data.frame(
    dgp = dgp, parameter = parameter, i = rep(seq_len(nrow(m)), ncol(m)),
    j = rep(seq_len(ncol(m)), each = nrow(m)), value = as.vector(m), label = label)
}
out_add("minar1_poi", "A", Theta$minar1_poi$A[[1]])
out_add("minar1_poi", "lambda", matrix(Theta$minar1_poi$innovation_params$lambda, ncol = 1L))
out_add("minar1_negbin", "A", Theta$minar1_negbin$A[[1]])
out_add("minar1_negbin", "r", matrix(Theta$minar1_negbin$innovation_params$r, ncol = 1L))
out_add("minar1_negbin", "p", matrix(Theta$minar1_negbin$innovation_params$p, ncol = 1L))
for (out_par in c("omega", "gamma", "delta", "alpha", "beta", "phi", "Sigma")) {
  out_add("mdacp1", out_par, if (is.matrix(Theta$mdacp1[[out_par]])) Theta$mdacp1[[out_par]] else
    matrix(Theta$mdacp1[[out_par]], ncol = 1L))
}
for (out_par in c("omega", "A", "B", "copula_param")) {
  out_add("mcar1", out_par, if (is.matrix(Theta$mcar1[[out_par]])) Theta$mcar1[[out_par]] else
    matrix(Theta$mcar1[[out_par]], ncol = 1L))
}
for (out_k in seq_along(Theta$svine1$margins)) {
  out_add("svine1", "margin_lambda", matrix(unname(coef(Theta$svine1$margins[[out_k]])), ncol = 1L),
          label = "poisson")
  out_param[[length(out_param)]]$i <- out_k
}
out_pc <- Theta$svine1$copula$pair_copulas
for (out_tree in seq_along(out_pc)) {
  for (out_edge in seq_along(out_pc[[out_tree]])) {
    out_bc <- out_pc[[out_tree]][[out_edge]]
    out_param[[length(out_param) + 1L]] <- data.frame(
      dgp = "svine1", parameter = "pair_copula", i = out_tree, j = out_edge,
      value = if (length(out_bc$parameters)) out_bc$parameters[1] else NA_real_,
      label = paste(out_bc$family, out_bc$rotation))
  }
}
out_param <- do.call(rbind, out_param)
out_write_csv(out_param, "dgp_parameters.csv")

# Stationarity: MINAR and MDACP need spectral radius < 1, count autoregression 1-norm sum < 1.
out_rho <- function(m) max(Mod(eigen(m, only.values = TRUE)$values))
out_mdacp_ab <- .mdacp_build_AB(Theta$mdacp1)
out_norm1 <- function(m) max(colSums(abs(m)))
out_stat <- data.frame(
  dgp = c("minar1_poi", "minar1_negbin", "mdacp1", "mcar1", "mcar1", "mcar1"),
  check = c("spectral_radius_A", "spectral_radius_A", "spectral_radius_A_plus_B",
            "one_norm_A", "one_norm_B", "one_norm_A_plus_one_norm_B"),
  value = c(out_rho(Theta$minar1_poi$A[[1]]), out_rho(Theta$minar1_negbin$A[[1]]),
            out_rho(out_mdacp_ab$A + out_mdacp_ab$B),
            out_norm1(Theta$mcar1$A), out_norm1(Theta$mcar1$B),
            out_norm1(Theta$mcar1$A) + out_norm1(Theta$mcar1$B)),
  bound = 1)
out_stat$satisfied <- c(out_stat$value[1:3] < out_stat$bound[1:3], NA, NA, out_stat$value[6] < out_stat$bound[6])
out_stat_csv <- out_stat
out_stat_csv$satisfied <- as.character(out_stat_csv$satisfied)
out_stat_csv$value <- sprintf("%.17g", out_stat_csv$value)
write.csv(out_stat_csv, file.path(out_dir, "stationarity_checks.csv"), row.names = FALSE, quote = FALSE)

out_stat_png <- data.frame(
  DGP = out_stat$dgp, Check = out_stat$check,
  Value = ifelse(grepl("^one_norm", out_stat$check), out_fmt(out_stat$value, 2), out_fmt(out_stat$value, 4)),
  Bound = "1", Satisfied = ifelse(is.na(out_stat$satisfied), "", ifelse(out_stat$satisfied, "yes", "no")))
out_png(out_stat_png, "stationarity_checks.png", width = 7.5, height = 3)
