# PMF and CDF plots of the count distributions in Chapter 2 of the thesis.
library(ggplot2)
library(patchwork)

# Color per parameter value
colors <- c(
  "#CC0000",
  "#009933",
  "#0033FF"
)

common_theme <- list(
  theme_minimal(base_size = 13) +
    theme(
      axis.title       = element_text(size = 13),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      legend.position  = "none"
    ),
  scale_color_manual(values = colors),
  scale_fill_manual(values = colors)
)


# PMF/CDF data frame for a discrete distribution.
# Truncate infinite supports yourself, e.g. at qpois(0.9999, lambda).
build_dist_data <- function(support, dfun, pfun, ..., group_label) {
  data.frame(
    x     = support,
    pmf   = dfun(support, ...),
    cdf   = pfun(support, ...),
    group = factor(group_label)
  )
}

# Horizontal offset per group so overlaid lollipops do not overlap.
add_dodge_offset <- function(data, dodge_width = 0.3, group_col = "group") {
  groups  <- levels(data[[group_col]])
  n       <- length(groups)
  offsets <- (seq_len(n) - (n + 1) / 2) * dodge_width
  offset_map <- setNames(offsets, groups)
  data$x_dodge <- data$x + offset_map[as.character(data[[group_col]])]
  data
}

# Plateau segments and closed/open jump points for the CDF step plot.
build_cdf_step_geometry <- function(data, group_col = "group") {
  groups <- split(data, data[[group_col]])
  segs   <- list()
  closed <- list()
  open   <- list()

  for (g in names(groups)) {
    df <- groups[[g]]
    df <- df[order(df$x), ]
    n  <- nrow(df)

    segs[[g]] <- data.frame(
      x = df$x[-n], xend = df$x[-1],
      y = df$cdf[-n], yend = df$cdf[-n],
      group = df[[group_col]][-n]
    )
    closed[[g]] <- data.frame(x = df$x, y = df$cdf, group = df[[group_col]])
    open[[g]] <- data.frame(
      x = df$x[-1], y = df$cdf[-n],
      group = df[[group_col]][-n]
    )
  }

  list(
    segments = do.call(rbind, segs),
    closed   = do.call(rbind, closed),
    open     = do.call(rbind, open)
  )
}

# Shared zoom window: mean +/- zoom_sd*SD over groups, clamped to the support.
compute_zoom_xlim <- function(data, mean_data, zoom_sd = 4) {
  lo <- min(mean_data$mean - zoom_sd * mean_data$sd)
  hi <- max(mean_data$mean + zoom_sd * mean_data$sd)
  c(max(lo, min(data$x)), min(hi, max(data$x)))
}

# PMF lollipops, zoomed to the CDF window so the pair stays comparable.
plot_pmf_lollipop <- function(data, mean_data, zoom_sd = 4,
                               dodge_width = 0.25, point_size = 2) {
  data <- add_dodge_offset(data, dodge_width)
  xlim <- compute_zoom_xlim(data, mean_data, zoom_sd)

  ggplot(data, aes(x = x_dodge, color = group)) +
    geom_segment(aes(xend = x_dodge, y = 0, yend = pmf), linewidth = 0.7) +
    geom_point(aes(y = pmf), size = point_size) +
    geom_vline(
      data = mean_data, aes(xintercept = mean, color = group),
      linetype = "dotted", linewidth = 0.8, show.legend = FALSE
    ) +
    coord_cartesian(xlim = xlim) +
    labs(x = "k", y = "P(X = k)", color = NULL) +
    common_theme
}

# CDF step function zoomed to mean +/- zoom_sd*SD, where the open/closed circles are visible.
plot_cdf_step_zoom <- function(data, mean_data, zoom_sd = 4, point_size = 2.2) {
  geo  <- build_cdf_step_geometry(data)
  xlim <- compute_zoom_xlim(data, mean_data, zoom_sd)

  ggplot() +
    geom_segment(
      data = geo$segments,
      aes(x = x, xend = xend, y = y, yend = yend, color = group),
      linewidth = 0.9
    ) +
    geom_point(
      data = geo$closed, aes(x = x, y = y, color = group),
      shape = 16, size = point_size          # closed = value attained
    ) +
    geom_point(
      data = geo$open, aes(x = x, y = y, color = group),
      shape = 1, size = point_size, stroke = 1.1  # open = left-limit
    ) +
    geom_vline(
      data = mean_data, aes(xintercept = mean, color = group),
      linetype = "dotted", linewidth = 0.8, show.legend = FALSE
    ) +
    coord_cartesian(xlim = xlim) +
    labs(x = "k", y = "P(X <= k)", color = NULL) +
    common_theme
}

# Combined PMF | CDF plot; mean_data (columns `mean`, `sd`) sets the shared zoom window.
combined_pmf_cdf_plot <- function(data, mean_data, zoom_sd = 4, dodge_width = 0.25) {
  pmf_plot <- plot_pmf_lollipop(data, mean_data, zoom_sd, dodge_width)
  cdf_plot <- plot_cdf_step_zoom(data, mean_data, zoom_sd)

  pmf_plot + cdf_plot +
    plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "none")
}

save_combined_plot <- function(plot, filename, folder = "chapter2_plots",
                                width = 12, height = 5.5, dpi = 300) {
  dir.create("chapter2_plots", showWarnings = FALSE)
  ggsave(filename = filename, path = folder, plot = plot,
         width = width, height = height, dpi = dpi)
}


# Binomial(n, p), varying p

n_trials <- 40
p_values <- c(0.25, 0.5, 0.75)

plot_data_binom <- do.call(rbind, lapply(p_values, function(p) {
  build_dist_data(0:n_trials, dbinom, pbinom, size = n_trials, prob = p,
                   group_label = paste0("p = ", p))
}))

mean_data_binom <- data.frame(
  group = factor(paste0("p = ", p_values)),
  mean  = n_trials * p_values,
  sd    = sqrt(n_trials * p_values * (1 - p_values))
)

combined_plot_binom <- combined_pmf_cdf_plot(plot_data_binom, mean_data_binom)
print(combined_plot_binom)
save_combined_plot(combined_plot_binom, "binomial_pmf_cdf.png")


# Poisson(lambda), varying mean

lambda_values <- c(5, 15, 25)

plot_data_pois <- do.call(rbind, lapply(lambda_values, function(lam) {
  x_max <- qpois(0.999, lam)  # infinite support -> truncate at a high quantile
  build_dist_data(0:x_max, dpois, ppois, lambda = lam,
                   group_label = paste0("lambda = ", lam))
}))

mean_data_pois <- data.frame(
  group = factor(paste0("lambda = ", lambda_values)),
  mean  = lambda_values,
  sd    = sqrt(lambda_values)
)

combined_plot_pois <- combined_pmf_cdf_plot(plot_data_pois, mean_data_pois)
print(combined_plot_pois)
save_combined_plot(combined_plot_pois, "poisson_pmf_cdf.png")


# Geometric(p), varying p

p_values <- c(0.3, 0.5, 0.7)

plot_data_geom <- do.call(rbind, lapply(p_values, function(p) {
  x_max <- qgeom(0.999, p)  # infinite support -> truncate at a high quantile
  build_dist_data(0:x_max, dgeom, pgeom, prob = p,
                   group_label = paste0("p = ", p))
}))

mean_data_geom <- data.frame(
  group = factor(paste0("p = ", p_values)),
  mean  = (1 - p_values) / p_values,
  sd    = sqrt(1 - p_values) / p_values
)

combined_plot_geom <- combined_pmf_cdf_plot(plot_data_geom, mean_data_geom)
print(combined_plot_geom)
save_combined_plot(combined_plot_geom, "geom_pmf_cdf.png")


# NegBin(r, p), fixed r, varying p

r_fixed  <- 10
p_values <- c(0.3, 0.5, 0.7)

plot_data_negbin1 <- do.call(rbind, lapply(p_values, function(p) {
  x_max <- qnbinom(0.999, size = r_fixed, prob = p)
  build_dist_data(0:x_max, dnbinom, pnbinom, size = r_fixed, prob = p,
                   group_label = paste0("p = ", p))
}))

mean_data_negbin1 <- data.frame(
  group = factor(paste0("p = ", p_values)),
  mean  = r_fixed * (1 - p_values) / p_values,
  sd    = sqrt(r_fixed * (1 - p_values)) / p_values
)

combined_plot_negbin1 <- combined_pmf_cdf_plot(plot_data_negbin1, mean_data_negbin1)
print(combined_plot_negbin1)
save_combined_plot(combined_plot_negbin1, "negbin1_pmf_cdf.png")


# NegBin(r, p), fixed p, varying r

p_fixed  <- 0.5
r_values <- c(10, 20, 30)

plot_data_negbin2 <- do.call(rbind, lapply(r_values, function(r) {
  x_max <- qnbinom(0.999, size = r, prob = p_fixed)
  build_dist_data(0:x_max, dnbinom, pnbinom, size = r, prob = p_fixed,
                   group_label = paste0("r = ", r))
}))

mean_data_negbin2 <- data.frame(
  group = factor(paste0("r = ", r_values)),
  mean  = r_values * (1 - p_fixed) / p_fixed,
  sd    = sqrt(r_values * (1 - p_fixed)) / p_fixed
)

combined_plot_negbin2 <- combined_pmf_cdf_plot(plot_data_negbin2, mean_data_negbin2)
print(combined_plot_negbin2)
save_combined_plot(combined_plot_negbin2, "negbin2_pmf_cdf.png")
