# Results table as a PNG (ggplot2 only), from results/abs_err_results_table.csv.
# Columns are the five series of series_spec. Defines plot_results_table().

plot_results_table <- function() {
  library(ggplot2)
  long <- read.csv(file.path(results_dir, "abs_err_results_table.csv"))
  rows <- c(names(groups), "Mean over groups")
  nc   <- nrow(series_spec)
  cw   <- 2.6; gw <- 2.4
  xc   <- c(gw / 2, gw + cw * (seq_len(nc) - 0.5))          # column centres: group, then the value columns
  xr   <- gw + cw * nc
  yr   <- rev(seq_along(rows)) - 0.5
  cells <- do.call(rbind, lapply(seq_along(rows), function(i) {
    s <- do.call(rbind, lapply(seq_len(nc), function(k)
      long[long$group == rows[i] & long$model == series_spec$model[k] & long$variant == series_spec$variant[k], ]))
    data.frame(x = xc[-1], y = yr[i], label = sprintf("%.3f (%.3f)", s$mean, s$se))
  }))
  top <- length(rows) + 2
  txt <- rbind(
    data.frame(x = xc[1], y = top - 1, label = "Group"),
    data.frame(x = xc[-1], y = top - 1, label = gsub(" \\(", "\n(", series_spec$label)),
    data.frame(x = xc[1], y = yr, label = rows),
    cells)
  hl <- function(y, x0 = 0, x1 = xr) data.frame(x = x0, xend = x1, y = y, yend = y)
  lines <- rbind(hl(top), hl(top - 2), hl(1), hl(0))
  vx <- c(0, gw, xr)
  vl <- data.frame(x = vx, xend = vx, y = 0, yend = top)
  p <- ggplot() +
    geom_segment(data = lines, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.4) +
    geom_segment(data = vl, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.4) +
    geom_text(data = txt, aes(x = x, y = y, label = label), size = 3.6, lineheight = 0.9) +
    coord_cartesian(xlim = c(-0.15, xr + 0.15), ylim = c(-0.15, top + 0.2), expand = FALSE) +
    theme_void() + theme(plot.background = element_rect(fill = "white", color = NA))
  save_plot(p, "abs_err_results_table", 10.5, 3.1)
}
