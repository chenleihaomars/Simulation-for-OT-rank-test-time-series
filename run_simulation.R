
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

script <- unlist(lapply(sys.frames(), function(frame) frame$ofile))
sourced <- length(script) > 0L
if (!sourced) {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  script <- gsub("~+~", " ", script, fixed = TRUE)
}
root <- if (length(script)) dirname(tail(script, 1L)) else getwd()
core_file <- file.path(root, "simulation_code.R")
if (!file.exists(core_file) && file.exists("simulation_code.R"))
  core_file <- "simulation_code.R"
if (!file.exists(core_file)) {
  if (!interactive())
    stop("Cannot find simulation_code.R. Keep both R files in the extracted folder.",
         call. = FALSE)
  message("Select simulation_code.R from the unzipped simulation folder.")
  core_file <- file.choose()
}
core_file <- normalizePath(core_file, mustWork = TRUE)
root <- dirname(core_file)

packages <- c("clue", "writexl")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
source(core_file)

cfg <- list(
  reps = 1000L,
  p = c(2L, 3L, 5L),
  M = c(1L, 3L, 5L),
  alpha = 0.05,
  burnin = 1000L,
  smooth_exp_tau = 0.25,
  ncores = max(1L, min(32L, parallel::detectCores() - 1L)),
  batch_size = 128L
)

stages <- if (sourced) character() else commandArgs(TRUE)
if (!length(stages)) stages <- c("observed", "fitted")
stages <- match.arg(stages, c("observed", "fitted"), several.ok = TRUE)
results <- file.path(root, "results")
for (setting in stages) {
  cfg$setting <- setting
  cfg$n <- if (setting == "fitted") c(250L, 500L, 1000L, 2000L) else
    c(125L, 250L, 500L, 1000L)
  cfg$cases <- if (setting == "fitted") 1:14 else 1:7
  cfg$outdir <- file.path(results, setting)
  run_experiment(cfg, core_file)
}
