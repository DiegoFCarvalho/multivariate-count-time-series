# Writes one d x T data file (.rds and .csv) per item group from the Kim et al. data.
# Defines prep_data().

prep_data <- function() {
  e <- new.env()
  load(file.path(kim_repo_dir, "Code", "Data.rda"), envir = e)
  x <- e$data
  stopifnot(is.matrix(x), nrow(x) == 30, ncol(x) == n_days, !anyNA(x), all(x == round(x)), all(x >= 1 & x <= n_cat))

  dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
  for (g in names(groups)) {
    stopifnot(length(groups[[g]]) == dimension, all(groups[[g]] %in% rownames(x)))
    xg <- x[groups[[g]], , drop = FALSE]
    storage.mode(xg) <- "integer"
    saveRDS(xg, file.path(data_dir, sprintf("%s_%s.rds", data_prefix, g)))
    # csv: values only, one row per item, no header
    write.table(xg, file.path(data_dir, sprintf("%s_%s.csv", data_prefix, g)),
                sep = ",", row.names = FALSE, col.names = FALSE)
  }
  invisible(NULL)
}
