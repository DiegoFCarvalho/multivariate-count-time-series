# Installs the missing CRAN packages the simulation needs; run once before the simulation.

if (getRversion() < "4.4.1") {
  stop("R >= 4.4.1 is required. Found R ", getRversion(), ".")
}

repos <- getOption("repos")
if (is.null(repos) || identical(unname(repos["CRAN"]), "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

# Package -> minimum version. svines needs rvinecopulib >= 1.0: its C++ code is compiled
# against rvinecopulib's headers, and builds against older headers compute a wrong
# discrete log-likelihood without any error. rvinecopulib is therefore installed first.
cran_pkgs <- c(
  gridExtra      = "0",
  ggplot2        = "0",
  future         = "0",
  future.apply   = "0",
  callr          = "0",
  mvtnorm        = "0",
  Matrix         = "0",
  "gamlss.dist"  = "0",
  univariateML   = "0",
  copula         = "0"
)
rvinecopulib_min <- "1.0.0.1.0"
svines_min       <- "0.3.0"

# Does not load the package: a loaded DLL cannot be overwritten on Windows.
is_ok <- function(pkg, min_version) {
  pkg %in% rownames(installed.packages()) && packageVersion(pkg) >= min_version
}

# `source` defaults to the package name; for an archive tarball pass its URL, repos = NULL
# and type = "source". Warnings are muffled so that they do not abort the installation.
install_or_stop <- function(pkg, hint = "", source = pkg, ...) {
  tryCatch(
    withCallingHandlers(install.packages(source, ...),
                        warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) NULL
  )
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Could not install '", pkg, "' from CRAN. ", hint,
         "See the installation messages above.", call. = FALSE)
  }
}

compiler_hint <- paste(
  "Building from source needs a compiler: Rtools on Windows,",
  "Xcode command line tools on macOS, build-essential (and libgsl-dev for 'gsl') on Linux."
)

# Linux only: check the system tools that source builds need.
# Prints what is missing and stops; it never installs system packages itself.
# A source build is needed on Linux unless a binary CRAN mirror is configured.
check_linux_build_tools <- function(need_gsl) {
  is_linux <- identical(Sys.info()[["sysname"]], "Linux")
  if (!is_linux || !identical(getOption("pkgType"), "source")) return(invisible(TRUE))
  missing <- character()
  if (!nzchar(Sys.which("g++")))  missing <- c(missing, "C++17 compiler (g++): Debian/Ubuntu: build-essential")
  if (!nzchar(Sys.which("make"))) missing <- c(missing, "make: Debian/Ubuntu: build-essential")
  if (need_gsl) {
    gsl_version <- if (nzchar(Sys.which("gsl-config"))) system2("gsl-config", "--version", stdout = TRUE) else NA
    if (is.na(gsl_version) || package_version(gsl_version) < "2.5") {
      missing <- c(missing, "GNU Scientific Library >= 2.5 with gsl-config: Debian/Ubuntu: libgsl-dev, RHEL/Fedora: gsl-devel")
    }
  }
  if (!nzchar(Sys.which("cmake"))) {
    cat("Note: 'cmake' was not found; the package 'nloptr' needs it if no system libnlopt is installed ",
        "(Debian/Ubuntu: cmake).\n", sep = "")
  }
  if (length(missing) > 0L) {
    stop("Missing system requirements for building R packages from source:\n  - ",
         paste(missing, collapse = "\n  - "),
         "\nAsk the system administrator to install them, then rerun this script.", call. = FALSE)
  }
  invisible(TRUE)
}

# gsl (needed by copula): the current CRAN release requires R >= 4.5.0, so on older R the
# newest release that supports it, 2.1-8 (R >= 4.0.0), is installed from the CRAN archive.
# It is installed before copula so that copula does not try to pull the incompatible release.
# Where a binary repository provides gsl for this R (Windows, macOS), the plain install works.
gsl_archive_version <- "2.1-8"
check_linux_build_tools(need_gsl = !is_ok("gsl", "0"))
if (!is_ok("gsl", "0")) {
  cat("Installing gsl\n")
  tryCatch(
    withCallingHandlers(install.packages("gsl"),
                        warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) NULL
  )
  if (!is_ok("gsl", "0") && getRversion() < "4.5.0") {
    install_or_stop(
      "gsl",
      hint = paste0("Archived gsl ", gsl_archive_version,
                    " needs the GNU Scientific Library >= 2.5 and a compiler. "),
      source = sprintf("https://cran.r-project.org/src/contrib/Archive/gsl/gsl_%s.tar.gz",
                       gsl_archive_version),
      repos = NULL, type = "source")
  }
  if (!is_ok("gsl", "0")) stop("Could not install 'gsl'. ", compiler_hint, call. = FALSE)
}

for (pkg in names(cran_pkgs)) {
  if (!is_ok(pkg, cran_pkgs[[pkg]])) {
    cat("Installing", pkg, "\n")
    install_or_stop(pkg, hint = compiler_hint)
  }
}

# rvinecopulib first; a newly installed or upgraded one forces svines to be rebuilt.
rebuild_svines <- FALSE
if (!is_ok("rvinecopulib", rvinecopulib_min)) {
  cat("Installing rvinecopulib\n")
  install_or_stop("rvinecopulib", hint = compiler_hint)
  if (!is_ok("rvinecopulib", rvinecopulib_min)) {
    stop("rvinecopulib >= ", rvinecopulib_min, " is required but CRAN provided ",
         as.character(packageVersion("rvinecopulib")), " for this R version.", call. = FALSE)
  }
  rebuild_svines <- TRUE
}

if (rebuild_svines || !is_ok("svines", svines_min)) {
  cat("Installing svines\n")
  install_or_stop("svines", hint = compiler_hint)
}


# Compares svines' discrete log-likelihood with a brute-force value: Markov order 1,
# one Poisson margin, Gaussian copula, where P(x_{t-1}, x_t) is a rectangle probability
# of the copula.
svines_loglik_matches_brute_force <- function(tol = 1e-6) {
  lambda <- 5
  bc_disc <- rvinecopulib::bicop_dist("gaussian", 0, 0.6, var_types = c("d", "d"))
  bc_cont <- rvinecopulib::bicop_dist("gaussian", 0, 0.6)
  cop <- svines::svinecop_dist(list(list(bc_disc)), rvinecopulib::cvine_structure(1), p = 1,
                               out_vertices = 1, in_vertices = 1, var_types = "d")
  margin <- univariateML::mlpois(rpois(20, lambda))
  margin[] <- lambda
  attr(margin, "type") <- "discrete"
  class(margin) <- c(class(margin), "svine_margin")
  model <- svines::svine_dist(margins = list(margin), copula = cop)

  old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
  set.seed(11)
  x <- as.integer(svines::svine_sim(n = 400, rep = 1, model = model))
  if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = globalenv())

  F <- function(k) ppois(k, lambda)
  C <- function(u, v) rvinecopulib::pbicop(cbind(u, v), bc_cont)
  a <- x[-length(x)]
  b <- x[-1]
  P <- C(F(a), F(b)) - C(F(a - 1), F(b)) - C(F(a), F(b - 1)) + C(F(a - 1), F(b - 1))
  brute <- sum(log(P)) - sum(dpois(x[2:(length(x) - 1)], lambda, log = TRUE))

  pkg_value <- svines::svine_loglik(matrix(x, ncol = 1), model)
  abs(pkg_value - brute) < tol
}

if (!svines_loglik_matches_brute_force()) {
  cat("svines log-likelihood check failed; rebuilding svines from source.\n")
  detach_ok <- tryCatch({ unloadNamespace("svines"); TRUE }, error = function(e) FALSE)
  install_or_stop("svines", hint = compiler_hint, type = "source")
  if (!detach_ok || !svines_loglik_matches_brute_force()) {
    stop("svines still computes a wrong discrete log-likelihood after a source rebuild. ",
         "Restart R and rerun this script; if it persists, check that rvinecopulib >= ",
         rvinecopulib_min, " is the version visible first in .libPaths().", call. = FALSE)
  }
}


all_pkgs <- c(names(cran_pkgs), "rvinecopulib", "svines")
cat("\nR", as.character(getRversion()), "\n")
for (pkg in all_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Package '", pkg, "' is not available.", call. = FALSE)
  cat(sprintf("  %-14s %s\n", pkg, as.character(packageVersion(pkg))))
}
cat("\nAll packages installed; svines log-likelihood check passed.\n")
