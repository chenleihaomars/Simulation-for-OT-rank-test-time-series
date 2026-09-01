options(repos = c(CRAN = "https://cloud.r-project.org"), error = NULL)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

script_path <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(arg)) {
    return(normalizePath(sub("^--file=", "", arg[1L]), winslash = "/"))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    path <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(path)) return(normalizePath(path, winslash = "/"))
  }

  path <- vapply(
    sys.frames(),
    function(frame) if (is.null(frame$ofile)) "" else as.character(frame$ofile[1L]),
    character(1L)
  )
  path <- path[nzchar(path)]
  if (length(path)) {
    return(normalizePath(tail(path, 1L), winslash = "/"))
  }

  normalizePath(file.choose(), winslash = "/")
}

ensure_package <- function(package, optional = FALSE) {
  if (requireNamespace(package, quietly = TRUE)) return(invisible(TRUE))

  type <- if (.Platform$OS.type == "windows") "binary" else getOption("pkgType")
  installed <- tryCatch({
    install.packages(package, dependencies = NA, type = type)
    requireNamespace(package, quietly = TRUE)
  }, error = function(e) FALSE)

  if (!installed && !optional) {
    stop("Could not install package '", package, "'.")
  }
  invisible(installed)
}

root <- dirname(script_path())
core_file <- file.path(root, "simulation_code.R")
if (!file.exists(core_file)) {
  stop("simulation_code.R must be in the same folder as run_simulation.R.")
}

ensure_package("clue")
ensure_package("writexl")
ensure_package("iGraphMatch", optional = TRUE)
source(core_file, local = .GlobalEnv)

ncores <- parallel::detectCores(logical = TRUE)
if (is.na(ncores)) ncores <- 32L

cfg <- list(
  reps = 1000L,
  n = c(250L, 500L, 1000L, 2000L, 4000L),
  p = c(2L, 3L, 5L),
  M = c(1L, 3L, 5L),
  cases = seq_len(16L),
  alpha = 0.05,
  burnin = 1000L,
  null_draws = 100000L,
  smooth_exp_tau = 0.25,
  ncores = max(1L, min(32L, ncores - 1L)),
  wave_per_worker = 2L,
  solver = "auto",              # auto, iGraphMatch, TreeDist, or clue
  allow_slow_solver = FALSE,
  resume = TRUE,
  outdir = root
)

results <- run_experiment(cfg, core_file)
selected <- results[results$M == 3L, ]

cat("\nMain size results (M = 3):\n")
print(selected[selected$role == "size_main", ], row.names = FALSE)
cat("\nMain power results (M = 3):\n")
print(selected[selected$role == "power_main", ], row.names = FALSE)
cat("\nRobustness results (M = 3):\n")
print(selected[grepl("_robustness$", selected$role), ], row.names = FALSE)
cat("\nSpreadsheet:\n", file.path(root, "simulation_results.xlsx"), "\n", sep = "")
