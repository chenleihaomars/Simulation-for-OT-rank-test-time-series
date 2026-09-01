# Simulation code

This repository contains the R implementation of the simulation study for the split-sample residual serial-independence test based on empirical optimal-transport ranks and tensor-cosine features.

## Files

- **`simulation_code.R`** contains the data-generating processes, model fitting, OT assignment, test statistics, Gaussian-chaos calibration, parallel execution, checkpointing, and output routines.
- **`run_simulation.R`** is the launcher. It locates `simulation_code.R`, checks the required packages, specifies the simulation configuration, runs the experiment, and prints selected summaries.

Keep both files in the same directory.

## Running an experiment

From a terminal, run

```bash
Rscript run_simulation.R
```

The script can also be sourced from RStudio. Simulation settings are collected in the `cfg` list in `run_simulation.R`. In particular, edit `reps`, `n`, `p`, `M`, `cases`, `ncores`, and `outdir` as needed. To run all implemented designs, use

```r
cases = CASE_INFO$case
```

The runner requires `clue` and `writexl`; it attempts to install `iGraphMatch` as a fast assignment solver. `TreeDist` is also supported when installed. With `solver = "auto"`, available fast solvers are checked against `clue` and the fastest valid solver is selected.

## Output

The output directory contains

- `simulation_results.xlsx`: complete summaries, design settings, timings, and solver information;
- `summary_results.csv`: rejection rates and Monte Carlo standard errors;
- `raw_results.rds`: replication-level results;
- hashed checkpoint files used when `resume = TRUE`.

A failed replication is not discarded: the run stops and writes a diagnostic CSV.

## Implementation notes

- Model parameters are estimated only from the training block; the fitted recursion is then continued through the test block to construct the test residuals.
- Deterministic substream seeds make results invariant to worker scheduling. For fixed `(p, n, replication)`, all cases share one independent OT reference grid, and its tensor-cosine features are evaluated once and reused.
- VAR coefficients are fitted by pivoted-QR least squares (`tol = 1e-10`). The structured VARMA model uses bounded conditional sum of squares, and the GARCH model uses variance-targeted Gaussian QMLE with an unconstrained parameterization enforcing `alpha + beta < 0.995`; fallback optimizers are used after primary-optimizer failure.
- Gaussian-chaos critical values are simulated once for each distinct `(p, K)`. Repeated product weights are grouped after rounding to 14 significant digits, and empirical quantiles use R's type-8 convention.

The master seed is fixed in `simulation_code.R` by `SIMULATION_SEED`.
