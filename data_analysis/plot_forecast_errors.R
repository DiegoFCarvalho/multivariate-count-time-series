# Figures from the saved results and data files (PNG, 300 dpi); defines functions only.
#   abs_err_convergence_<group>_K<K>.png  MAE of the draws vs. Monte Carlo replicates:
#                                         running mean with a +/- 1 SE band, five series of series_spec
#   abs_err_item_day_<group>.png          per item and day; rows: data, absolute error (five series)

# Style of the five series of series_spec. Green, red and blue are the colours of the original figures; the scaled-ecdf
# S-vines are purple. Only MDACP clipped is dashed; point shapes also separate the series.
series_style <- data.frame(
  color    = c("#009900", "#CC0000", "#CC0000", "#0033FF", "#9B1FFF"),
  linetype = c("solid", "solid", "dashed", "solid", "solid"),
  shape    = c(16, 17, 2, 15, 5),
  stringsAsFactors = FALSE)
rownames(series_style) <- series_spec$label

# data frame of the five series from a long table with columns model and variant
select_series <- function(df) {
  out <- lapply(seq_len(nrow(series_spec)), function(k) {
    z <- df[df$model == series_spec$model[k] & df$variant == series_spec$variant[k], ]
    if (nrow(z)) { z$series <- series_spec$label[k]; z }
  })
  out <- do.call(rbind, out)
  out$series <- factor(out$series, levels = series_spec$label)
  out
}

# manual scales for the aesthetics a plot maps (a scale for an unmapped aesthetic triggers a warning)
series_scales <- function(aesthetics = c("color", "linetype")) {
  lab <- series_spec$label
  all <- list(color    = ggplot2::scale_color_manual(values = setNames(series_style$color, lab), breaks = lab),
              fill     = ggplot2::scale_fill_manual(values = setNames(series_style$color, lab), breaks = lab),
              linetype = ggplot2::scale_linetype_manual(values = setNames(series_style$linetype, lab), breaks = lab),
              shape    = ggplot2::scale_shape_manual(values = setNames(series_style$shape, lab), breaks = lab))
  unname(all[aesthetics])
}

plot_theme <- function(base_size = 14) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(plot.title       = ggplot2::element_text(size = base_size, hjust = 0.5),
                   axis.title       = ggplot2::element_text(size = base_size),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major = ggplot2::element_line(color = "grey90"),
                   legend.position  = "top")
}

save_plot <- function(plot, name, width, height) {
  dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(file.path(plots_dir, paste0(name, ".png")), plot = plot, width = width, height = height,
                  units = "in", dpi = 300)
}

# vertices of the stepped area under the ratings y of the days t (bars of width 1, no gaps)
step_polygon <- function(t, y) {
  n <- length(t)
  data.frame(x = c(t[1] - 0.5, as.vector(rbind(t - 0.5, t + 0.5)), t[n] + 0.5),
             y = c(0, as.vector(rbind(y, y)), 0))
}

plot_files <- function(K) {
  g <- names(groups)
  file.path(plots_dir, c(sprintf("abs_err_convergence_%s_K%d.png", g, K), sprintf("abs_err_item_day_%s.png", g),
                         "data_overview.png", "abs_err_results_table.png"))
}
plots_done <- function(K) all(file.exists(plot_files(K)))

conv_breaks <- function(K) {
  b <- pretty(c(1, K))
  unique(c(1, b[b > 1 & b <= K]))
}

plot_abs_err_convergence <- function(K) {
  library(ggplot2)
  run <- read.csv(file.path(results_dir, "abs_err_running.csv"))
  for (g in names(groups)) {
    df <- select_series(run[run$group == g, ])
    p <- ggplot(df, aes(x = k, y = running_mean, color = series, fill = series, linetype = series, group = series)) +
      geom_ribbon(aes(ymin = running_mean - running_se, ymax = running_mean + running_se), alpha = 0.10, color = NA, na.rm = TRUE) +
      geom_line(linewidth = 0.9, na.rm = TRUE) +
      scale_x_continuous(limits = c(1, K), breaks = conv_breaks(K)) +
      series_scales(c("color", "fill", "linetype")) +
      labs(x = "Monte Carlo replicates", y = "Mean absolute error", color = NULL, fill = NULL, linetype = NULL) +
      guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1), fill = "none") +
      plot_theme() + theme(legend.key.width = grid::unit(1.6, "lines"), legend.text = element_text(size = 10))
    save_plot(p, sprintf("abs_err_convergence_%s_K%d", g, K), 10, 5.2)
  }
}

plot_abs_err_item_day <- function() {
  library(ggplot2); library(patchwork)
  md    <- select_series(read.csv(file.path(results_dir, "abs_err_item_day.csv")))
  ymax  <- ceiling(max(md$abs_err_item_day, na.rm = TRUE))   # shared error axis
  days  <- n_train + 0:horizon
  fdays <- n_train + seq_len(horizon)
  small <- theme(axis.text = element_text(size = 6), axis.title = element_text(size = 8),
                 plot.title = element_text(size = 9, hjust = 0.5),
                 panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey92"),
                 panel.border = element_rect(fill = NA, color = "grey70"), plot.margin = margin(2, 3, 2, 3))
  for (g in names(groups)) {
    x <- readRDS(file.path(data_dir, sprintf("%s_%s.rds", data_prefix, g)))
    panels <- list()
    for (r in 1:2) for (i in seq_along(groups[[g]])) {
      item <- groups[[g]][i]
      if (r == 1L) {                                   # last observed day in black, forecast days in gray
        p <- ggplot() +
          geom_polygon(data = step_polygon(n_train, x[item, n_train]), aes(x, y), fill = "#000000", color = NA) +
          geom_polygon(data = step_polygon(fdays, x[item, fdays]), aes(x, y), fill = "#7F7F7F", color = NA) +
          scale_y_continuous(limits = c(0, n_cat), breaks = c(1, 3, 5), expand = c(0, 0))
        ylab <- "Data"
      } else {
        df <- md[md$group == g & md$item == item, ]
        df <- df[order(df$series, df$day), ]
        p <- ggplot(df, aes(day, abs_err_item_day, color = series, linetype = series, shape = series, group = series)) +
          geom_line(linewidth = 0.6) + geom_point(size = 1.4, stroke = 0.6) +
          series_scales(c("color", "linetype", "shape")) + labs(color = NULL, linetype = NULL, shape = NULL) +
          scale_y_continuous(limits = c(0, ymax), breaks = 0:ymax, expand = expansion(mult = c(0, 0.02)))
        ylab <- "Mean absolute error"
      }
      p <- p + scale_x_continuous(breaks = days, limits = range(days) + c(-0.5, 0.5), expand = c(0, 0)) +
        labs(title = if (r == 1L) item else NULL, x = if (r == 2L) "t" else NULL, y = if (i == 1L) ylab else NULL) +
        theme_minimal(base_size = 8) + small
      panels[[length(panels) + 1L]] <- p
    }
    fig <- wrap_plots(panels, ncol = length(groups[[g]]), guides = "collect") &
      theme(legend.position = "bottom", legend.text = element_text(size = 8), legend.key.width = grid::unit(1.8, "lines")) &
      guides(color = guide_legend(nrow = 1), linetype = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
    save_plot(fig, sprintf("abs_err_item_day_%s", g), 12, 4.6)
  }
  invisible(ymax)
}

run_plots <- function(K) {
  plot_abs_err_convergence(K)
  plot_abs_err_item_day()
  plot_data_overview()
  plot_results_table()
  invisible(NULL)
}
