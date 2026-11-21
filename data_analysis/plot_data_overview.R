# Data overview of the 30 items as daily step profiles (black: estimation, gray: forecast days).
# Defines functions only.

check_overview_data <- function() {
  for (g in names(groups)) {
    x <- readRDS(file.path(data_dir, sprintf("%s_%s.rds", data_prefix, g)))
    stopifnot(identical(rownames(x), groups[[g]]), all(!is.na(x)), all(x >= 1 & x <= n_cat & x == round(x)),
              ncol(x) == n_days, all(rowSums(!is.na(x)) == n_days))
  }
}

plot_data_overview <- function() {
  library(ggplot2); library(patchwork)
  check_overview_data()
  panels <- list()
  for (gi in seq_along(groups)) {
    g <- names(groups)[gi]
    x <- readRDS(file.path(data_dir, sprintf("%s_%s.rds", data_prefix, g)))
    for (i in seq_len(nrow(x))) {
      tr <- seq_len(n_train); te <- n_train + seq_len(horizon)
      p <- ggplot() +
        annotate("rect", xmin = n_train + 0.5, xmax = n_days + 0.5, ymin = 0, ymax = n_cat, fill = "gray92") +
        geom_polygon(data = step_polygon(tr, x[i, tr]), aes(x, y), fill = "#000000", color = NA) +
        geom_polygon(data = step_polygon(te, x[i, te]), aes(x, y), fill = "#7F7F7F", color = NA) +
        scale_x_continuous(limits = c(0.5, n_days + 0.5), breaks = c(30, 60, 90), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, n_cat), breaks = c(1, 3, 5), expand = c(0, 0)) +
        labs(title = rownames(x)[i], x = if (gi == length(groups)) "t" else NULL, y = if (i == 1L) g else NULL) +
        theme_minimal(base_size = 8) +
        theme(legend.position = "none", axis.text = element_text(size = 6), axis.title = element_text(size = 8),
              plot.title = element_text(size = 8, hjust = 0.5), panel.grid.minor = element_blank(),
              panel.grid.major = element_line(color = "grey92"), panel.border = element_rect(fill = NA, color = "grey70"),
              plot.margin = margin(2, 3, 2, 3))
      panels[[length(panels) + 1L]] <- p
    }
  }
  save_plot(wrap_plots(panels, ncol = dimension), "data_overview", 12, 8)
}
