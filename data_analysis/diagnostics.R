# Diagnostics of the saved draws and fits. Defines functions only.
# Reads saved draws and fit records; Kim variances before abs() are recomputed (estimation only).
# Outputs in results_dir:
#   draw_support.csv            draw counts and percentages by value range, per model and group
#   draw_support_items.csv      per item: percentage of draws equal to 0 and above 5
#   fit_diagnostics_kim.csv     per group and item: Cov_eps before abs(), categories, Psi radius
#   fit_diagnostics_svines.csv  per S-vine model and group: margin per item, margin and pair-copula family counts, copula parameters
#   draw_support_table.csv, fit_diagnostics_table.csv   the displayed tables
# Tables as PNG in plots_dir: draw_support_table.png, fit_diagnostics_table.png

diagnostics_models <- function() c(svine_main_models, mdacp_model, "kim")

diagnostics_files <- function() c(file.path(results_dir, c("draw_support.csv", "draw_support_items.csv",
  "fit_diagnostics_kim.csv", "fit_diagnostics_svines.csv", "draw_support_table.csv", "fit_diagnostics_table.csv")),
  file.path(plots_dir, c("draw_support_table.png", "fit_diagnostics_table.png")))
diagnostics_done <- function() all(file.exists(diagnostics_files()))


support_row <- function(model, group, x) {
  n <- length(x)
  r <- data.frame(model = model, group = group, n_values = n,
                  n_zero = sum(x == 0), n_in_1_5 = sum(x %in% 1:5), n_above_5 = sum(x > 5), n_below_0 = sum(x < 0),
                  n_non_integer = sum(x != round(x)), min = min(x), max = max(x))
  r$pct_zero <- 100 * r$n_zero / n; r$pct_in_1_5 <- 100 * r$n_in_1_5 / n; r$pct_above_5 <- 100 * r$n_above_5 / n
  r$pct_below_0 <- 100 * r$n_below_0 / n; r$pct_non_integer <- 100 * r$n_non_integer / n
  r$pct_outside_1_5 <- 100 * (n - r$n_in_1_5) / n
  r
}

diag_draw_support <- function(K) {
  rows <- list(); items <- list()
  for (m in diagnostics_models()) {
    pooled <- NULL
    for (g in names(groups)) {
      r <- readRDS(raw_file(m, g, K)); x <- r$raw
      rows[[length(rows) + 1L]] <- support_row(m, g, x)
      pooled <- c(pooled, as.vector(x))
      items[[length(items) + 1L]] <- data.frame(model = m, group = g, item = r$series_names,
        pct_zero = sapply(seq_len(dim(x)[2]), function(j) 100 * mean(x[, j, ] == 0)),
        pct_above_5 = sapply(seq_len(dim(x)[2]), function(j) 100 * mean(x[, j, ] > 5)))
    }
    rows[[length(rows) + 1L]] <- support_row(m, "pooled", pooled)
  }
  list(support = do.call(rbind, rows), items = do.call(rbind, items))
}


# same calls as kim_group(); estimation only
kim_estimate <- function(g) {
  x     <- read_group(g)
  d     <- nrow(x)
  X_t   <- x[, seq_len(n_train), drop = FALSE]
  Link  <- Latent_Gauss_Link(d, n_train, "multinom", X_t, n_cat)
  Cov_X <- Latent_Gauss_Cov(d, n_train, kim_p, X_t)
  Cov_Z <- Latent_Gauss_invLink(d, kim_p, Cov_X, Link$Link)
  list(x = X_t, Est = Latent_DFM_Estim(r = kim_r, kim_p, Cov_Z, identy_opt = 2, shift = TRUE))
}

diag_kim_fit <- function(K) {
  load_kim_code()
  out <- list()
  for (g in names(groups)) {
    setTimeLimit(elapsed = est_timeout, transient = TRUE)
    e <- kim_estimate(g); setTimeLimit()
    Est <- e$Est; X_t <- e$x
    saved <- readRDS(fit_file("kim", g, K))
    max_diff <- max(abs(abs(Est$Cov_eps) - saved$Cov_eps), abs(Est$Lambda - saved$Lambda),
                    abs(Est$Psi - saved$Psi), abs(Est$Cov_eta - saved$Cov_eta))
    if (max_diff > 1e-8) stop(sprintf("recomputed Kim fit of %s differs from the saved record (%.3g)", g, max_diff))
    counts   <- t(apply(X_t, 1L, tabulate, nbins = n_cat))
    observed <- rowSums(counts > 0)
    out[[g]] <- data.frame(group = g, item = rownames(X_t), cov_eps_before_abs = Est$Cov_eps,
      negative = Est$Cov_eps < 0, n_negative_group = sum(Est$Cov_eps < 0),
      spectral_radius_psi = max(Mod(eigen(matrix(Est$Psi[, , 1], kim_r, kim_r))$values)),
      n_categories_observed = observed, unobserved_category = observed < n_cat,
      missing_categories = apply(counts, 1L, function(v) paste(which(v == 0), collapse = " ")),
      n_items_unobserved_group = sum(observed < n_cat), max_abs_diff_saved_fit = max_diff)
  }
  do.call(rbind, out)
}


# empirical margins; the exact ecdf margins carry their cumulative frequencies in attr "cf"
margin_label  <- function(m) "empirical"
margin_params <- function(m) {
  if (!is.null(attr(m, "cf"))) paste(sprintf("%.3f", diff(c(0, attr(m, "cf")))), collapse = "; ")
  else ""
}

diag_svines_fit <- function(K) {
  out <- list()
  for (m in svine_main_models) for (g in names(groups)) {
    th  <- readRDS(fit_file(m, g, K))$theta_hat
    nms <- rownames(read_group(g))
    fam <- vapply(th$margins, margin_label, "")
    par <- vapply(th$margins, margin_params, "")
    pc  <- unlist(lapply(th$copula$pair_copulas, function(tr) vapply(tr, function(b) b$family, "")))
    cnt <- function(v, kind) { t <- table(v); data.frame(model = m, group = g, kind = kind, item = NA_character_,
                               family = names(t), params = NA_character_, count = as.integer(t)) }
    out[[paste(m, g)]] <- rbind(
      data.frame(model = m, group = g, kind = "margin", item = nms, family = fam, params = par, count = NA_integer_),
      cnt(fam, "margin_family_count"), cnt(pc, "pair_copula_count"),
      data.frame(model = m, group = g, kind = "copula_parameters", item = NA_character_, family = NA_character_,
                 params = NA_character_, count = as.integer(th$copula$npars)))
  }
  do.call(rbind, out)
}


fmt_pct <- function(v) sprintf("%.1f", v)

draw_support_table <- function(sup) {
  rows <- c(names(groups), "pooled")
  pick <- function(m, col) sapply(rows, function(g) sup[[col]][sup$model == m & sup$group == g])
  data.frame(Group = ifelse(rows == "pooled", "Pooled", rows),
             exact_eq0 = fmt_pct(pick("svines_exact", "pct_zero")), exact_gt5 = fmt_pct(pick("svines_exact", "pct_above_5")),
             exact_nonint = fmt_pct(pick("svines_exact", "pct_non_integer")),
             scaled_eq0 = fmt_pct(pick("svines_scaled", "pct_zero")), scaled_gt5 = fmt_pct(pick("svines_scaled", "pct_above_5")),
             scaled_nonint = fmt_pct(pick("svines_scaled", "pct_non_integer")),
             mdacp_eq0 = fmt_pct(pick(mdacp_model, "pct_zero")), mdacp_gt5 = fmt_pct(pick(mdacp_model, "pct_above_5")),
             mdacp_nonint = fmt_pct(pick(mdacp_model, "pct_non_integer")),
             kim_outside = fmt_pct(pick("kim", "pct_outside_1_5")), row.names = NULL)
}

fit_diagnostics_table <- function(kimfit, svfit) {
  do.call(rbind, lapply(names(groups), function(g) {
    k <- kimfit[kimfit$group == g, ]
    neg <- k[k$negative, ]
    neg_txt <- if (nrow(neg) == 0L) "0" else
      sprintf("%d (%s)", nrow(neg), paste(sprintf("%s: %.3f", neg$item, neg$cov_eps_before_abs), collapse = ", "))
    np <- function(m) svfit$count[svfit$model == m & svfit$group == g & svfit$kind == "copula_parameters"]
    data.frame(Group = g, kim_negative_variances = neg_txt,
               kim_spectral_radius_psi = sprintf("%.3f", k$spectral_radius_psi[1]),
               kim_items_unobserved_category = as.character(k$n_items_unobserved_group[1]),
               svines_copula_parameters = sprintf("%d / %d", np("svines_exact"), np("svines_scaled")))
  }))
}

# PNG table; blocks = list(label, columns), cw = column widths, sep_before = rows with a rule above
table_png <- function(name, tab, blocks, cw, sep_before = NULL, notes, height) {
  library(ggplot2)
  n  <- nrow(tab); nc <- ncol(tab) - 1L; gw <- 2.4
  xr <- gw + sum(cw); xc <- gw + cumsum(cw) - cw / 2
  xe <- gw + cumsum(cw)                              # right edge of each value column
  yr <- rev(seq_len(n)) - 0.5; top <- n + 2
  b_end <- cumsum(sapply(blocks, function(b) length(b$cols)))
  b_start <- c(1L, head(b_end, -1L) + 1L)
  b_x <- mapply(function(s, e) mean(c(xe[s] - cw[s], xe[e])), b_start, b_end)
  cells <- do.call(rbind, lapply(seq_len(n), function(i) data.frame(x = xc, y = yr[i], label = unlist(tab[i, -1]))))
  txt <- rbind(data.frame(x = b_x, y = top - 0.5, label = sapply(blocks, `[[`, "label")),
               data.frame(x = xc, y = top - 1.5, label = unlist(lapply(blocks, `[[`, "cols"))),
               data.frame(x = gw / 2, y = top - 1.5, label = "Group"),
               data.frame(x = gw / 2, y = yr, label = tab[[1]]), cells)
  hl <- function(y, x0 = 0, x1 = xr) data.frame(x = x0, xend = x1, y = y, yend = y)
  lines <- rbind(hl(top), hl(top - 1, gw + 0.1), hl(top - 2), hl(0),
                 if (!is.null(sep_before)) hl(n - sep_before + 1))
  vx <- c(0, gw, xe[b_end]); vl <- data.frame(x = vx, xend = vx, y = 0, yend = top)
  p <- ggplot() +
    geom_segment(data = lines, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.4) +
    geom_segment(data = vl, aes(x = x, xend = xend, y = y, yend = yend), linewidth = 0.4) +
    geom_text(data = txt, aes(x = x, y = y, label = label), size = 3.6, lineheight = 0.9) +
    (if (length(notes)) annotate("text", x = 0, y = -0.6 - 0.5 * (seq_along(notes) - 1), label = notes, hjust = 0, size = 2.6)) +
    coord_cartesian(xlim = c(-0.15, xr + 0.15), ylim = c(if (length(notes)) -0.6 - 0.5 * length(notes) else -0.15, top + 0.2), expand = FALSE) +
    theme_void() + theme(plot.background = element_rect(fill = "white", color = NA))
  save_plot(p, name, xr * 0.74, height)
}

plot_draw_support_table <- function(sup, t1 = draw_support_table(sup)) {
  table_png("draw_support_table", t1,
            list(list(label = "S-vines (exact ecdf)", cols = c("= 0 (%)", "> 5 (%)", "non-integer (%)")),
                 list(label = "S-vines (scaled ecdf)", cols = c("= 0 (%)", "> 5 (%)", "non-integer (%)")),
                 list(label = "MDACP", cols = c("= 0 (%)", "> 5 (%)", "non-integer (%)")),
                 list(label = "Latent Gaussian", cols = "outside {1,...,5} (%)")),
            cw = c(rep(c(1.9, 1.9, 2.8), 3), 3.2), sep_before = nrow(t1), notes = character(0), height = 2.8)
}

# The S-vines copula parameter counts stay in fit_diagnostics_table.csv; the PNG shows the Latent Gaussian columns only.
plot_fit_diagnostics_table <- function(t2) {
  table_png("fit_diagnostics_table", t2[, names(t2) != "svines_copula_parameters"],
            list(list(label = "Latent Gaussian", cols = c("Negative variances\n(item: value)", "Spectral radius\nof Psi", "Items with unobserved\ncategory"))),
            cw = c(4.6, 2.4, 3.0), sep_before = NULL, notes = character(0), height = 2.8)
}

run_diagnostics <- function(K) {
  dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
  wr <- function(df, f) write.csv(df, file.path(results_dir, f), row.names = FALSE)
  ds <- diag_draw_support(K)
  kimfit <- diag_kim_fit(K); svfit <- diag_svines_fit(K)
  wr(ds$support, "draw_support.csv"); wr(ds$items, "draw_support_items.csv")
  wr(kimfit, "fit_diagnostics_kim.csv"); wr(svfit, "fit_diagnostics_svines.csv")

  t1 <- draw_support_table(ds$support); t2 <- fit_diagnostics_table(kimfit, svfit)
  wr(t1, "draw_support_table.csv"); wr(t2, "fit_diagnostics_table.csv")

  plot_draw_support_table(ds$support)
  plot_fit_diagnostics_table(t2)
  invisible(NULL)
}
