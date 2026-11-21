# Mean absolute error of the draws (abs_err) from the saved forecast draws. Defines functions only.
# For draw k and day t, abs_err is the mean over items of |draw - observed|; the group value is the
# mean over days, then over draws (rating points). Outputs in results_dir:
#   abs_err_table.csv     group value, Monte Carlo SE, clipped fraction per group, model, variant
#   abs_err_item_day.csv  mean over the draws of |draw - observed| per item and forecast day
#   abs_err_running.csv   running mean and SE over the first k draws
#   abs_err_results_table.csv           long table (mean and SE per group, model and variant)
#   abs_err_results_table.tex           the five series of `series_spec` as columns (mean and SE per group)
# Models: svines_exact, svines_scaled, mdacp, kim. Variants: raw draws, and draws clipped to {1,...,5}
# (rounded for svines_scaled; Kim's draws are in 1..5 and svines_exact's draws are integers in 1..5, so clipping
# changes neither).

# The five series of the results table and the figures: label, model and draw variant.
series_spec <- data.frame(
  label   = c("Latent Gaussian", "MDACP unclipped", "MDACP clipped", "S-vines (exact ecdf)", "S-vines (scaled ecdf)"),
  model   = c("kim", "mdacp", "mdacp", "svines_exact", "svines_scaled"),
  variant = c("raw", "raw", "clipped", "clipped", "clipped"),
  stringsAsFactors = FALSE)

# post-hoc clipping of draws to the rating range
clip_draws <- function(x) pmin(pmax(round(x), 1), 5)

# running mean and SE over the valid draws among the first k
running_stats <- function(v) {
  valid   <- !is.na(v)
  n_valid <- cumsum(valid)
  v0      <- replace(v, !valid, 0)
  mean_k  <- ifelse(n_valid >= 1L, cumsum(v0) / n_valid, NA_real_)
  var_k   <- ifelse(n_valid >= 2L, (cumsum(v0^2) - n_valid * mean_k^2) / (n_valid - 1L), NA_real_)
  data.frame(k = seq_along(v), running_mean = mean_k, running_se = sqrt(pmax(var_k, 0)) / sqrt(n_valid))
}

metric_files <- function() file.path(results_dir, c("abs_err_table.csv", "abs_err_item_day.csv", "abs_err_running.csv",
                                                     "abs_err_results_table.csv", "abs_err_results_table.tex"))
metrics_done <- function(K) all(file.exists(metric_files()))

# long table (group, variant, model, mean, se) with the mean over groups appended (mean of the
# group means; its SE is the mean of the group SEs)
results_table_long <- function(tab) {
  long <- data.frame(group = tab$group, variant = tab$variant, model = tab$model, mean = tab$abs_err_mean, se = tab$abs_err_se)
  agg  <- aggregate(cbind(mean, se) ~ variant + model, long, mean)
  rbind(long, data.frame(group = "Mean over groups", agg[c("variant", "model", "mean", "se")]))
}

# LaTeX tabular with the series of series_spec as columns
write_results_table_tex <- function(long, file) {
  rows <- c(names(groups), "Mean over groups")
  val  <- function(g, k) long[long$group == g & long$model == series_spec$model[k] & long$variant == series_spec$variant[k], ]
  lines <- c(sprintf("\\begin{tabular}{|c|%s}", paste(rep("c|", nrow(series_spec)), collapse = "")), "\\hline",
             paste("Group &", paste(series_spec$label, collapse = " & "), "\\\\ \\hline"))
  for (g in rows) {
    v <- lapply(seq_len(nrow(series_spec)), function(k) val(g, k))
    cells <- sapply(v, function(s) sprintf("%.3f (%.3f)", s$mean, s$se))
    lines <- c(lines, paste(if (g == "Mean over groups") "\\hline Mean over groups" else g, "&",
                            paste(cells, collapse = " & "), "\\\\"))
  }
  writeLines(c(lines, "\\hline", "\\end{tabular}"), file)
}

run_metrics <- function(K) {
  tab <- list(); runn <- list(); itemday <- list()
  for (g in names(groups)) {
    x     <- read_group(g)
    truth <- x[, n_train + seq_len(horizon), drop = FALSE]
    items <- rownames(x)
    draws <- list(svines_exact  = readRDS(raw_file("svines_exact", g, K))$raw,
                  svines_scaled = readRDS(raw_file("svines_scaled", g, K))$raw,
                  mdacp         = readRDS(raw_file(mdacp_model, g, K))$raw,
                  kim           = readRDS(raw_file("kim", g, K))$raw)
    variants <- list(svines_exact = c("raw", "clipped"), svines_scaled = c("raw", "clipped"),
                     mdacp = c("raw", "clipped"), kim = "raw")

    for (m in names(draws)) for (v in variants[[m]]) {
      raw  <- if (v == "clipped") clip_draws(draws[[m]]) else draws[[m]]
      draw <- apply(abs_err_grid(raw, truth), 1L, mean_or_na)
      row  <- data.frame(group = g, model = m, variant = v, K = K,
                         abs_err_mean = mean_or_na(draw), abs_err_se = sd(draw, na.rm = TRUE) / sqrt(sum(!is.na(draw))),
                         clipped_fraction = if (v == "clipped") mean(draws[[m]] != raw, na.rm = TRUE) else 0,
                         n_failed_draws = sum(is.na(draw)))
      tab[[length(tab) + 1L]] <- row
      if (m == "kim") { row$variant <- "clipped"; tab[[length(tab) + 1L]] <- row }   # Kim: clipped = raw
      runn[[length(runn) + 1L]] <- cbind(data.frame(group = g, model = m, variant = v), running_stats(draw))

      ae <- abs(raw - aperm(array(truth, c(dim(truth), dim(raw)[1])), c(3L, 1L, 2L)))
      me <- apply(ae, c(2L, 3L), mean_or_na)
      for (vv in if (m == "kim") c("raw", "clipped") else v)
        itemday[[length(itemday) + 1L]] <- data.frame(group = g, item = rep(items, times = horizon),
            day = rep(n_train + seq_len(horizon), each = length(items)), horizon = rep(seq_len(horizon), each = length(items)),
            model = m, variant = vv, abs_err_item_day = as.vector(me))
    }
  }
  tab <- do.call(rbind, tab); runn <- do.call(rbind, runn); itemday <- do.call(rbind, itemday)

  # the item-day values must reproduce the group values
  chk <- aggregate(abs_err_item_day ~ group + model + variant, itemday, mean)
  chk <- merge(chk, tab[c("group", "model", "variant", "abs_err_mean")], by = c("group", "model", "variant"))
  max_diff <- max(abs(chk$abs_err_item_day - chk$abs_err_mean))
  if (max_diff > 1e-10) stop("abs_err_item_day does not reproduce abs_err_table")
  message(sprintf("check: max |item-day mean - group value| = %.3g", max_diff))

  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  wr <- function(df, f) write.csv(df, file.path(results_dir, f), row.names = FALSE)
  wr(tab, "abs_err_table.csv"); wr(itemday, "abs_err_item_day.csv"); wr(runn, "abs_err_running.csv")
  long <- results_table_long(tab)
  wr(long, "abs_err_results_table.csv")
  write_results_table_tex(long, file.path(results_dir, "abs_err_results_table.tex"))
  invisible(NULL)
}
