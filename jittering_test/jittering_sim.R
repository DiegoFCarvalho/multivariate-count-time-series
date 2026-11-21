# Validation of discrete vine-copula estimation: exact discrete likelihood vs. jittered data (MSE).
# Writes mse_summary_table.{png,csv}, mse_convergence_lambda_combined.png and tables/*.csv here.

library(rvinecopulib)
library(ggplot2)


# Generate discrete multivariate data

# Canonical family name per alias (e.g. "pois" = "poisson"), so later switches need one case.
standardize_discrete_family <- function(family) {
  family <- tolower(gsub("[- ]", "_", family))
  if (family %in% c("bernoulli", "binary")) return("bernoulli")
  if (family %in% c("binomial", "binom")) return("binomial")
  if (family %in% c("poisson", "pois")) return("poisson")
  if (family %in% c("negative_binomial", "nbinomial", "nbinom")) {
    return("negative_binomial")
  }
  if (family %in% c("geometric", "geom")) return("geometric")

  stop(
    "Unsupported marginal distribution: '", family, "'. ",
    "Use bernoulli, binomial, poisson, negative_binomial, or geometric."
  )
}


# Default parameters for one family; used only when marginal_parameters = NULL.
default_discrete_parameters <- function(family) {
  switch(
    family,
    bernoulli = list(prob = 0.5),
    binomial = list(size = 5, prob = 0.5),
    poisson = list(lambda = 2),
    negative_binomial = list(size = 5, prob = 0.5),
    geometric = list(prob = 0.5)
  )
}


# Inverse CDF of one discrete margin (bernoulli = binomial with size 1).
discrete_margin_quantile <- function(u, family, parameters) {
  arguments <- c(list(p = u), parameters)
  switch(
    family,
    bernoulli = do.call(qbinom, c(arguments, list(size = 1))),
    binomial = do.call(qbinom, arguments),
    poisson = do.call(qpois, arguments),
    negative_binomial = do.call(qnbinom, arguments),
    geometric = do.call(qgeom, arguments)
  )
}


# CDF of one discrete margin; F(x - 1) is the left limit F(x-) rvinecopulib needs for discrete data.
discrete_margin_cdf <- function(x, family, parameters) {
  arguments <- c(list(q = x), parameters)
  switch(
    family,
    bernoulli = do.call(pbinom, c(arguments, list(size = 1))),
    binomial = do.call(pbinom, arguments),
    poisson = do.call(ppois, arguments),
    negative_binomial = do.call(pnbinom, arguments),
    geometric = do.call(pgeom, arguments)
  )
}


# Simulate discrete data from a known vine copula with the same pair copula on every edge.
# copula_parameter: scalar, or vector for multi-parameter families (t: c(rho, df)).
# marginal_parameters: NULL (defaults), one named list for all margins, or d named lists.
# Returns the data, the 2d-column copula data, the latent uniforms and the true model.
generate_discrete_copula_data <- function(
  n,
  dimension,
  vine_type = c("dvine", "cvine"),
  copula_family = "clayton",
  copula_parameter = 2,
  marginal_distributions = rep("poisson", dimension),
  marginal_parameters = NULL,
  rotation = 0,
  seed = NULL
) {
  if (length(n) != 1L || n < 3 || n != as.integer(n)) {
    stop("'n' must be an integer of at least 3.")
  }
  if (
    length(dimension) != 1L ||
    dimension < 2 ||
    dimension != as.integer(dimension)
  ) {
    stop("'dimension' must be an integer of at least 2.")
  }

  n <- as.integer(n)
  dimension <- as.integer(dimension)
  vine_type <- match.arg(tolower(vine_type[1L]), c("dvine", "cvine"))

  # rvinecopulib names these "gaussian" and "t"; the long names are accepted as aliases.
  copula_family <- tolower(copula_family[1L])
  if (copula_family == "normal") copula_family <- "gaussian"
  if (copula_family == "student") copula_family <- "t"

  if (length(marginal_distributions) == 1L) {
    marginal_distributions <- rep(
      marginal_distributions,
      dimension
    )
  }
  if (length(marginal_distributions) != dimension) {
    stop(
      "'marginal_distributions' must contain one family per variable."
    )
  }
  marginal_distributions <- unname(vapply(
    marginal_distributions,
    standardize_discrete_family,
    character(1L)
  ))

  if (is.null(marginal_parameters)) {
    marginal_parameters <- lapply(
      marginal_distributions,
      default_discrete_parameters
    )
  } else if (!is.list(marginal_parameters)) {
    stop("'marginal_parameters' must be a list.")
  } else if (!is.null(names(marginal_parameters))) {
    marginal_parameters <- rep(
      list(marginal_parameters),
      dimension
    )
  }

  if (
    length(marginal_parameters) != dimension ||
    !all(vapply(marginal_parameters, is.list, logical(1L)))
  ) {
    stop(
      "'marginal_parameters' must contain one parameter list per variable."
    )
  }

  # Pair copula with the known true parameter, marked discrete; the same one on every edge.
  common_pair_copula <- bicop_dist(
    family = copula_family,
    rotation = rotation,
    parameters = copula_parameter,
    var_types = c("d", "d")
  )

  # d - 1 trees with d - 1, ..., 1 edges; the common pair copula is replicated over this layout.
  number_edges_per_tree <- dimension - seq_len(dimension - 1L)
  pair_copulas <- lapply(
    number_edges_per_tree,
    function(number_edges) {
      rep(list(common_pair_copula), number_edges)
    }
  )

  # Reversing the order passed to cvine_structure() gives the intended C-vine root order.
  vine_structure <- if (vine_type == "dvine") {
    dvine_structure(seq_len(dimension))
  } else {
    cvine_structure(rev(seq_len(dimension)))
  }

  true_model <- vinecop_dist(
    pair_copulas = pair_copulas,
    structure = vine_structure,
    var_types = rep("d", dimension)
  )

  if (!is.null(seed)) set.seed(seed)

  latent_uniform <- rvinecop(
    n = n,
    vinecop = true_model
  )

  discrete_data <- vapply(
    seq_len(dimension),
    function(variable_index) {
      discrete_margin_quantile(
        u = latent_uniform[, variable_index],
        family = marginal_distributions[variable_index],
        parameters = marginal_parameters[[variable_index]]
      )
    },
    numeric(n)
  )
  colnames(discrete_data) <- paste0("x", seq_len(dimension))
  discrete_data <- as.data.frame(discrete_data)

  # A discrete observation is an interval on the copula scale, so rvinecopulib needs F_j(x_j) and
  # F_j(x_j-) = F_j(x_j - 1): copula_data = cbind(upper, lower) with 2d columns.
  cdf_upper <- vapply(
    seq_len(dimension),
    function(variable_index) {
      discrete_margin_cdf(
        x = discrete_data[[variable_index]],
        family = marginal_distributions[variable_index],
        parameters = marginal_parameters[[variable_index]]
      )
    },
    numeric(n)
  )

  cdf_lower <- vapply(
    seq_len(dimension),
    function(variable_index) {
      discrete_margin_cdf(
        x = discrete_data[[variable_index]] - 1,
        family = marginal_distributions[variable_index],
        parameters = marginal_parameters[[variable_index]]
      )
    },
    numeric(n)
  )

  colnames(cdf_upper) <- names(discrete_data)
  colnames(cdf_lower) <- paste0(names(discrete_data), "_left")

  list(
    data = discrete_data,
    copula_data = cbind(cdf_upper, cdf_lower),
    latent_uniform = latent_uniform,
    true_model = true_model,
    true_parameters = as.data.frame(summary(true_model)),
    marginal_distributions = marginal_distributions,
    marginal_parameters = marginal_parameters
  )
}


# Estimate the known model and calculate parameter MSE

# Fit the pair-copula parameters to discrete copula data (2d columns: F_j(x_j), F_j(x_j-)).
# The true model is the fixed template.
estimate_discrete_copula <- function(
  copula_data,
  true_model,
  cores = 1L
) {
  # The structure order has one element per variable, which gives d for the 2d-column check.
  dimension <- length(true_model$structure$order)
  copula_data <- as.matrix(copula_data)

  if (nrow(copula_data) < 3L) {
    stop("'copula_data' must contain at least three observations.")
  }
  if (ncol(copula_data) != 2L * dimension) {
    stop(
      "'copula_data' must contain d CDF columns followed by ",
      "d left-limit CDF columns."
    )
  }
  if (anyNA(copula_data) || any(!is.finite(copula_data))) {
    stop("'copula_data' must not contain missing or infinite values.")
  }

  # Structure, families and rotations stay fixed, so tree t, edge e is the same pair copula in the
  # true and fitted models (needed for an edge-by-edge MSE).
  # vinecop(vinecop_object = ...) needs class "vinecop"; the generating "vinecop_dist" gets the
  # fitting class so its parameters can be updated.
  model_template <- true_model
  class(model_template) <- c("vinecop", "vinecop_dist")

  fitted_model <- vinecop(
    data = copula_data,
    vinecop_object = model_template,
    par_method = "mle",
    cores = max(1L, as.integer(cores))
  )

  fitted_model
}


# Compare true and estimated parameters per edge; MSE = mean[(estimated - true)^2] over parameters.
# Families differ in parameter count: one row per (tree, edge, parameter_index) after expansion.
calculate_parameter_mse <- function(true_model, fitted_model) {
  extract_parameters <- function(model, value_name) {
    model_summary <- as.data.frame(summary(model))

    parameter_rows <- lapply(
      seq_len(nrow(model_summary)),
      function(row_index) {
        parameter_values <- as.numeric(
          model_summary$parameters[[row_index]]
        )

        data.frame(
          tree = model_summary$tree[row_index],
          edge = model_summary$edge[row_index],
          # parameter_index separates the parameters of one pair copula (e.g. rho and df for t).
          parameter_index = seq_along(parameter_values),
          value = parameter_values
        )
      }
    )

    result <- do.call(rbind, parameter_rows)
    names(result)[names(result) == "value"] <- value_name
    result
  }

  comparison <- merge(
    extract_parameters(true_model, "true_parameter"),
    extract_parameters(fitted_model, "estimated_parameter"),
    by = c("tree", "edge", "parameter_index"),
    all = FALSE,
    sort = TRUE
  )

  comparison$squared_error <- with(
    comparison,
    (estimated_parameter - true_parameter)^2
  )

  if (
    nrow(comparison) == 0L ||
    any(!is.finite(comparison$squared_error))
  ) {
    stop("The true and estimated parameters could not be compared.")
  }

  list(
    mse = mean(comparison$squared_error),
    comparison = comparison
  )
}


# Continuous extension via jittering

# Independent Uniform(0,1) jitter per column: x_tilde = x + (U - 1); `seed` affects only the jitter.
jitter_discrete_data <- function(discrete_data, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  n <- nrow(discrete_data)
  dimension <- ncol(discrete_data)

  # Independent noise per column adds no dependence beyond the discrete data.
  jittered_data <- vapply(
    seq_len(dimension),
    function(variable_index) {
      x <- discrete_data[[variable_index]]
      u <- runif(n)
      x + (u - 1)
    },
    numeric(n)
  )
  colnames(jittered_data) <- paste0(names(discrete_data), "_tilde")

  as.data.frame(jittered_data)
}


# pmf of one discrete margin; used by jittered_margin_cdf().
discrete_margin_pmf <- function(x, family, parameters) {
  arguments <- c(list(x = x), parameters)
  switch(
    family,
    bernoulli = do.call(dbinom, c(arguments, list(size = 1))),
    binomial = do.call(dbinom, arguments),
    poisson = do.call(dpois, arguments),
    negative_binomial = do.call(dnbinom, arguments),
    geometric = do.call(dgeom, arguments)
  )
}


# Exact cdf of X_tilde = X + (U - 1), U ~ Uniform(0, 1) independent of X:
#   F_tilde(x) = F(floor(x)) + f(floor(x) + 1) * (x - floor(x)),
# used in estimate_jittered_copula() instead of empirical pseudo-observations.
jittered_margin_cdf <- function(x, family, parameters) {
  vapply(
    x,
    function(xi) {
      floor_x <- floor(xi)
      discrete_margin_cdf(floor_x, family, parameters) +
        discrete_margin_pmf(floor_x + 1, family, parameters) *
          (xi - floor_x)
    },
    numeric(1L)
  )
}


# Fit pair-copula parameters to jittered data with the exact jittered cdf, not pseudo-observations.
estimate_jittered_copula <- function(
  jittered_data,
  marginal_distributions,
  marginal_parameters,
  true_model,
  cores = 1L
) {
  dimension <- length(true_model$structure$order)
  jittered_data <- as.data.frame(jittered_data)
  n <- nrow(jittered_data)

  transformed_data <- vapply(
    seq_len(dimension),
    function(variable_index) {
      jittered_margin_cdf(
        x = jittered_data[[variable_index]],
        family = marginal_distributions[variable_index],
        parameters = marginal_parameters[[variable_index]]
      )
    },
    numeric(n)
  )
  colnames(transformed_data) <- names(jittered_data)

  # Same template as true_model with variable types switched from "d" to "c" (var_types of the
  # vinecop_dist and of each nested bicop_dist).
  continuous_template <- true_model
  continuous_template$var_types <- rep("c", dimension)
  continuous_template$pair_copulas <- lapply(
    true_model$pair_copulas,
    function(tree_copulas) {
      lapply(tree_copulas, function(pair_copula) {
        pair_copula$var_types <- rep("c", length(pair_copula$var_types))
        pair_copula
      })
    }
  )
  class(continuous_template) <- c("vinecop", "vinecop_dist")

  # The jittered scale has no left-limit ambiguity, so d columns suffice (not 2d).
  fitted_model <- vinecop(
    data = transformed_data,
    vinecop_object = continuous_template,
    par_method = "mle",
    cores = max(1L, as.integer(cores))
  )

  fitted_model
}


# Monte Carlo comparison: count estimation vs. jittering

# Exact vs. jittered estimation for n_replicates on PSOCK workers (mclapply() is serial on Windows).
# Workers use cores = 1; per-replicate seeds: base_seed + k (data), base_seed + 100000 + k (jitter).
# Returns the long-format results, a per-method MSE summary and the true parameters.
compare_count_vs_jittered <- function(
  n = 1000,
  dimension = 4,
  vine_type = "dvine",
  copula_family = "clayton",
  copula_parameter = 2,
  marginal_distributions = rep("poisson", dimension),
  marginal_parameters = NULL,
  rotation = 0,
  cores = 1L,
  n_replicates = 200,
  base_seed = 2028
) {
  n_replicates <- as.integer(n_replicates)

  # Inputs are passed explicitly (no closure), so workers do not rely on environment serialization.
  run_one_replicate <- function(
    k, n, dimension, vine_type, copula_family, copula_parameter,
    marginal_distributions, marginal_parameters, rotation, base_seed
  ) {
    simulation <- generate_discrete_copula_data(
      n = n,
      dimension = dimension,
      vine_type = vine_type,
      copula_family = copula_family,
      copula_parameter = copula_parameter,
      marginal_distributions = marginal_distributions,
      marginal_parameters = marginal_parameters,
      rotation = rotation,
      seed = base_seed + k
    )

    count_fit <- estimate_discrete_copula(
      copula_data = simulation$copula_data,
      true_model = simulation$true_model,
      cores = 1L
    )
    count_comparison <- calculate_parameter_mse(
      true_model = simulation$true_model,
      fitted_model = count_fit
    )$comparison

    jittered_data <- jitter_discrete_data(
      simulation$data,
      seed = base_seed + 100000L + k
    )
    jittered_fit <- estimate_jittered_copula(
      jittered_data = jittered_data,
      marginal_distributions = simulation$marginal_distributions,
      marginal_parameters = simulation$marginal_parameters,
      true_model = simulation$true_model,
      cores = 1L
    )
    jittered_comparison <- calculate_parameter_mse(
      true_model = simulation$true_model,
      fitted_model = jittered_fit
    )$comparison

    rbind(
      data.frame(replicate = k, method = "count", count_comparison),
      data.frame(replicate = k, method = "jittered", jittered_comparison)
    )
  }

  n_workers <- max(1L, as.integer(cores))
  cl <- parallel::makePSOCKcluster(n_workers)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, library(rvinecopulib))
  parallel::clusterExport(
    cl,
    varlist = c(
      "generate_discrete_copula_data", "estimate_discrete_copula",
      "estimate_jittered_copula", "calculate_parameter_mse",
      "jitter_discrete_data", "jittered_margin_cdf", "discrete_margin_cdf",
      "discrete_margin_pmf", "discrete_margin_quantile",
      "standardize_discrete_family", "default_discrete_parameters"
    ),
    envir = globalenv()
  )

  # parLapply() has no progress callback, so progress is reported only before and after the batch.
  cat(sprintf(
    "compare_count_vs_jittered: running %d replicates across %d workers...\n",
    n_replicates, n_workers
  ))
  replicate_rows <- parallel::parLapply(
    cl, seq_len(n_replicates), run_one_replicate,
    n = n, dimension = dimension, vine_type = vine_type,
    copula_family = copula_family, copula_parameter = copula_parameter,
    marginal_distributions = marginal_distributions,
    marginal_parameters = marginal_parameters, rotation = rotation,
    base_seed = base_seed
  )
  cat("compare_count_vs_jittered: done\n")

  results <- do.call(rbind, replicate_rows)
  rownames(results) <- NULL

  # The true model is identical across replicates, so the reference values come from results.
  true_parameters <- unique(
    results[, c("tree", "edge", "parameter_index", "true_parameter")]
  )
  true_parameters <- true_parameters[
    order(
      true_parameters$tree,
      true_parameters$edge,
      true_parameters$parameter_index
    ),
  ]
  rownames(true_parameters) <- NULL

  # One MSE value per replicate and method, as in calculate_parameter_mse().
  per_replicate_mse <- aggregate(
    squared_error ~ replicate + method,
    data = results,
    FUN = mean
  )
  names(per_replicate_mse)[names(per_replicate_mse) == "squared_error"] <-
    "mse"

  summary_table <- do.call(rbind, lapply(
    split(per_replicate_mse, per_replicate_mse$method),
    function(method_rows) {
      data.frame(
        method = method_rows$method[1L],
        mean = mean(method_rows$mse),
        median = median(method_rows$mse),
        sd = sd(method_rows$mse),
        mc_se = sd(method_rows$mse) / sqrt(n_replicates)
      )
    }
  ))
  rownames(summary_table) <- NULL

  list(
    results = results,
    summary = summary_table,
    true_parameters = true_parameters
  )
}


# Cumulative-average MSE per method (log scale, running Monte Carlo SE ribbon) over replicates.
plot_mse_convergence <- function(mc_results, y_log = TRUE) {
  # One MSE value per replicate and method.
  per_replicate_mse <- aggregate(
    squared_error ~ replicate + method,
    data = mc_results,
    FUN = mean
  )
  names(per_replicate_mse)[names(per_replicate_mse) == "squared_error"] <-
    "mse"

  # Running mean and SE (sd / sqrt(k)) over the first k replicates; SE is 0 at k = 1 (sd undefined).
  convergence_data <- do.call(rbind, lapply(
    split(per_replicate_mse, per_replicate_mse$method),
    function(method_rows) {
      method_rows <- method_rows[order(method_rows$replicate), ]
      k <- seq_len(nrow(method_rows))
      cumulative_mean <- cumsum(method_rows$mse) / k

      running_se <- vapply(
        k,
        function(i) {
          if (i < 2L) return(0)
          sd(method_rows$mse[seq_len(i)]) / sqrt(i)
        },
        numeric(1L)
      )

      data.frame(
        replicate_count = k,
        method = method_rows$method[1L],
        cumulative_mean = cumulative_mean,
        running_se = running_se
      )
    }
  ))
  rownames(convergence_data) <- NULL

  # Relabel the method codes; the labels also key the manual palette.
  convergence_data$method <- factor(
    convergence_data$method,
    levels = c("count", "jittered"),
    labels = c("Count (exact discrete)", "Jittered (continuized)")
  )

  convergence_data$ribbon_min <-
    convergence_data$cumulative_mean - convergence_data$running_se
  convergence_data$ribbon_max <-
    convergence_data$cumulative_mean + convergence_data$running_se
  if (y_log) {
    # Floor the ribbon's lower bound above zero so scale_y_log10() does not drop segments.
    convergence_data$ribbon_min <- pmax(convergence_data$ribbon_min, 1e-8)
  }

  y_axis_label <- "MSE"
  if (y_log) y_axis_label <- paste0(y_axis_label, " (log scale)")

  method_colors <- c(
    "Count (exact discrete)" = "#0033FF",
    "Jittered (continuized)" = "#CC0000"
  )

  convergence_plot <- ggplot(
    convergence_data,
    aes(
      x = replicate_count,
      y = cumulative_mean,
      color = method,
      fill = method,
      group = method
    )
  ) +
    geom_ribbon(
      aes(ymin = ribbon_min, ymax = ribbon_max),
      alpha = 0.15,
      color = NA
    ) +
    # Both lines are solid; methods differ by color alone and are identified in the figure caption.
    geom_line(linewidth = 1, linetype = "solid") +
    labs(
      x = "Monte Carlo replicates",
      y = y_axis_label,
      color = "Estimation approach",
      fill = "Estimation approach"
    )

  if (y_log) convergence_plot <- convergence_plot + scale_y_log10()

  # Use the thesis-wide common_theme if defined, else a default; the manual palette is applied last.
  convergence_plot <- convergence_plot +
    (if (exists("common_theme")) common_theme else theme_minimal()) +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    theme(legend.position = "none")

  convergence_plot
}


# Jittering vs. exact discrete, lambda sensitivity (thesis Section \ref{sec:jittering_simulation})

# Runnable record of the thesis study: Poisson margins, 4-dimensional Clayton D-vine, theta = 2,
# n = 500, 200 replicates per lambda, base_seed = 2028;
# regenerates mse_convergence_lambda_combined.png.

library(patchwork)
library(gridExtra)

# Script directory (machine-independent output_dir): Rscript --file, then the call-stack 'ofile',
# then rstudioapi. Stops instead of using getwd(), which could silently misplace the figures.
get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }

  for (frame in sys.frames()) {
    if (!is.null(frame$ofile)) {
      return(dirname(normalizePath(frame$ofile)))
    }
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    active_path <- tryCatch(
      rstudioapi::getActiveDocumentContext()$path,
      error = function(e) ""
    )
    if (nzchar(active_path)) {
      return(dirname(normalizePath(active_path)))
    }
  }

  stop(
    "Could not determine the script's own directory (not run via Rscript, ",
    "not sourced, and no active RStudio document found). Run this file with ",
    "Rscript, source() it (or click Source in RStudio), or set output_dir ",
    "manually before this line."
  )
}

output_dir <- get_script_dir()

# Theme copied from chapter_2/distribution_plots.R (sourcing that file fails here: its
# ggsave(path = "chapter2_plots") needs a relative folder).
# common_theme is the alias plot_mse_convergence() looks for.
common_theme <- theme_minimal(base_size = 13) +
  theme(
    axis.title       = element_text(size = 13),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    legend.position  = "none"
  )

# Three lambda settings, all else fixed; each call re-simulates 200 data sets (base_seed + 1..200)
# and fits both estimators.
result_lambda5 <- compare_count_vs_jittered(
  n = 500,
  dimension = 4,
  vine_type = "dvine",
  copula_family = "clayton",
  copula_parameter = 2,
  marginal_distributions = "poisson",
  marginal_parameters = list(lambda = 5),
  n_replicates = 200,
  cores = 4,
  base_seed = 2028
)

result_lambda50 <- compare_count_vs_jittered(
  n = 500,
  dimension = 4,
  vine_type = "dvine",
  copula_family = "clayton",
  copula_parameter = 2,
  marginal_distributions = "poisson",
  marginal_parameters = list(lambda = 50),
  n_replicates = 200,
  cores = 4,
  base_seed = 2028
)

result_lambda250 <- compare_count_vs_jittered(
  n = 500,
  dimension = 4,
  vine_type = "dvine",
  copula_family = "clayton",
  copula_parameter = 2,
  marginal_distributions = "poisson",
  marginal_parameters = list(lambda = 250),
  n_replicates = 200,
  cores = 4,
  base_seed = 2028
)

# Thesis Table \ref{tab:jittering_mse}: mean and median per-replicate squared error per lambda.
extract_summary_rows <- function(result, lambda) {
  summary_table <- result$summary
  data.frame(
    lambda = lambda,
    Method = c("Discrete", "Jittered"),
    MSE = round(
      c(
        summary_table$mean[summary_table$method == "count"],
        summary_table$mean[summary_table$method == "jittered"]
      ),
      3
    ),
    `Median SE` = round(
      c(
        summary_table$median[summary_table$method == "count"],
        summary_table$median[summary_table$method == "jittered"]
      ),
      3
    ),
    check.names = FALSE
  )
}

mse_summary_table <- rbind(
  extract_summary_rows(result_lambda5, 5),
  extract_summary_rows(result_lambda50, 50),
  extract_summary_rows(result_lambda250, 250)
)

mse_summary_grob <- gridExtra::tableGrob(mse_summary_table, rows = NULL)

ggsave(
  file.path(output_dir, "mse_summary_table.png"),
  plot = mse_summary_grob, width = 6, height = 3, units = "in", dpi = 300
)

# Full-precision CSVs of the summary table and of the per-replicate values behind it.
tables_dir <- file.path(output_dir, "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

lambda_results <- list(`5` = result_lambda5, `50` = result_lambda50, `250` = result_lambda250)

write_full_precision <- function(df, path) {
  for (j in which(vapply(df, is.double, logical(1)))) df[[j]] <- sprintf("%.17g", df[[j]])
  write.csv(df, path, row.names = FALSE, quote = FALSE)
}

write_full_precision(
  do.call(rbind, lapply(names(lambda_results), function(lambda) {
    result <- lambda_results[[lambda]]
    summary_table <- result$summary[match(c("count", "jittered"), result$summary$method), ]
    data.frame(
      lambda = as.numeric(lambda),
      Method = c("Discrete", "Jittered"),
      MSE = summary_table$mean,
      `Median SE` = summary_table$median,
      sd = summary_table$sd,
      mc_se = summary_table$mc_se,
      n_replicates = length(unique(result$results$replicate)),
      check.names = FALSE
    )
  })),
  file.path(output_dir, "mse_summary_table.csv")
)

write_full_precision(
  do.call(rbind, lapply(names(lambda_results), function(lambda) {
    per_replicate <- aggregate(squared_error ~ replicate + method, data = lambda_results[[lambda]]$results, FUN = mean)
    data.frame(lambda = as.numeric(lambda), per_replicate[order(per_replicate$method, per_replicate$replicate), ],
               row.names = NULL)
  })),
  file.path(tables_dir, "mse_per_replicate.csv")
)

write_full_precision(
  do.call(rbind, lapply(names(lambda_results), function(lambda) {
    data.frame(lambda = as.numeric(lambda), lambda_results[[lambda]]$true_parameters)
  })),
  file.path(tables_dir, "true_parameters.csv")
)

plot_lambda5 <- plot_mse_convergence(result_lambda5$results, y_log = TRUE)
plot_lambda50 <- plot_mse_convergence(result_lambda50$results, y_log = TRUE)
plot_lambda250 <- plot_mse_convergence(result_lambda250$results, y_log = TRUE)

plot_lambda_combined <- plot_lambda5 + plot_lambda50 + plot_lambda250 +
  plot_layout(ncol = 3)

# Filename and size as used for the thesis figure.
ggsave(
  file.path(output_dir, "mse_convergence_lambda_combined.png"),
  plot = plot_lambda_combined, width = 14, height = 4.5, units = "in",
  dpi = 300
)
