#' scComBat-Seq: Master Reproducibility Script
#'
#' Reproduces all paper figures and result tables.
#'
#' Usage:
#'   Rscript run_all.R            # full run (all cells)
#'   Rscript run_all.R --fast     # dev run (10% subset per dataset)
#'
#' Outputs go into results/<dataset_id>/ and results/simulation/.
#' Figures go into results/<dataset_id>/figures/.

library(here)

args     <- commandArgs(trailingOnly = TRUE)
FAST     <- "--fast" %in% args
SUBSET   <- if (FAST) 0.1 else 1.0

message(if (FAST) "[dev mode] subset = 0.1" else "[full run] subset = 1.0")

# ── Load all function modules ─────────────────────────────────────────────────
source(here("R", "core", "combat_nb_math.R"))
source(here("R", "core", "preprocessing.R"))
source(here("R", "core", "qc_checks.R"))
source(here("R", "core", "evaluation.R"))
source(here("R", "core", "combat_variants.R"))
source(here("R", "core", "benchmark_methods.R"))
source(here("R", "visualization", "plot_umap.R"))
source(here("R", "simulation", "simulate_splatter.R"))
source(here("R", "simulation", "simulate_scdesign3.R"))
source(here("analysis", "run_combat_variants.R"))

# ── Section 1 (Methods): Simulation experiments ───────────────────────────────
message("\n===== SECTION 1: SIMULATION =====")
source(here("analysis", "run_simulation.R"))

# ── Section 2 (Results): Real datasets ───────────────────────────────────────
message("\n===== SECTION 2: REAL DATA =====")

DATASETS <- c("dataset1", "dataset2", "dataset4", "dataset5",
              "dataset6", "dataset7", "dataset8", "dataset9", "dataset10")

for (ds in DATASETS) {
  tryCatch(
    run_dataset_analysis(ds, subset = SUBSET),
    error = function(e) message("[WARN] ", ds, " failed: ", conditionMessage(e))
  )
}

message("\n===== DONE =====")
message("Results in: ", here("results"))
