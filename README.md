# OT-rank time-series simulations

## Run

Extract the ZIP, keep the files together, and run in R:

```r
source(file.choose())
```

Select `run_simulation.R`. The observed study runs first, followed by the fitted study. Missing packages (`clue` and `writexl`) are installed automatically.

Alternatively, run from a terminal in the extracted folder:

```sh
Rscript run_simulation.R          # both studies
Rscript run_simulation.R observed # observed only
Rscript run_simulation.R fitted   # fitted only
```

## Settings

| Setting | Designs | Sample sizes | Tests |
|---|---|---|---|
| Observed | O1–O7 | 125, 250, 500, 1000 | Pairwise aggregate and block maximum |
| Fitted | F1–F14 | 250, 500, 1000, 2000 | Learned projection |

Both studies use **1,000 replications per design cell**, dimensions **2, 3, 5**, horizons **1, 3, 5**, and level **0.05**. All tests use analytical calibration. For fitted models, `n` is the total sample size before splitting.

Settings are specified in `run_simulation.R`; the methods are implemented in `simulation_code.R`. Change `ncores` to choose the number of parallel workers. The default uses up to 32 workers. Seeds are fixed, so results do not depend on worker scheduling.

## Output and resuming

Results are saved in `results/observed` and `results/fitted`. Each folder contains CSV and Excel summaries of rejection frequencies and simulation standard errors, replication-level statistics, and records of the settings, software versions, and run times.

Checkpoints are saved after each batch. Run the same launcher again to resume with the same code and settings.
