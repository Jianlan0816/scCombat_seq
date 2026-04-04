#' Simulation Experiments for the scComBat-Seq Paper
#'
#' Generates all simulation scenarios described in the Methods section,
#' runs ComBat variants and competitor benchmarks on each, and saves metrics.
#'
#' Scenarios:
#'   S1 — Balanced, 5k cells, no covariate
#'   S2 — Imbalanced cell type composition, 5k cells
#'   S3 — Balanced with categorical covariate ("Young"/"Old"), 5k cells
#'   S4 — Large balanced, 40k cells (use subset for dev)
#'   SD — scDesign3-based realistic simulation (requires dataset1_sce.rds)

library(here)

source(here("R", "core", "evaluation.R"))
source(here("R", "core", "combat_variants.R"))
source(here("R", "core", "benchmark_methods.R"))
source(here("R", "core", "qc_checks.R"))
source(here("R", "simulation", "simulate_splatter.R"))
source(here("R", "simulation", "simulate_scdesign3.R"))
source(here("R", "visualization", "plot_umap.R"))

SIM_DIR    <- here("results", "simulation")
SUBSET_DEV <- 0.2   # change to 1.0 for final paper run

run_simulation_scenario <- function(id, seu, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if ("Group" %in% colnames(seu@meta.data) && !"celltype" %in% colnames(seu@meta.data)) {
    seu$celltype <- seu$Group
  }
  check_data_characteristics(seu)

  combat_res <- ComBat_combo(
    seu,
    subset       = SUBSET_DEV,
    print_raw    = TRUE,
    combat_seq   = TRUE,
    combat_scseq = TRUE,
    combat_pcseq = TRUE,
    combat_ind   = TRUE,
    save_path    = file.path(out_dir, "combat_umaps")
  )

  bench_res <- benchmark_batch_methods(
    seu,
    subset    = SUBSET_DEV,
    save_path = file.path(out_dir, "benchmark_umaps")
  )

  saveRDS(combat_res, file.path(out_dir, "combat_results.rds"))
  saveRDS(bench_res,  file.path(out_dir, "benchmark_results.rds"))

  all_metrics <- rbind(
    cbind(source = "combat",    summarise_metrics(combat_res)),
    cbind(source = "benchmark", summarise_metrics(bench_res))
  )
  write.csv(all_metrics, file.path(out_dir, "all_metrics.csv"), row.names = FALSE)
  message(sprintf("[%s] Metrics saved.", id))
  print(all_metrics)

  assemble_figure(c(combat_res, bench_res),
                  save_path = file.path(out_dir, "figures"))

  invisible(list(combat = combat_res, benchmark = bench_res))
}

# ── Splatter scenarios ────────────────────────────────────────────────────────

message("=== S1: Balanced, 5k cells ===")
seu_S1 <- simulate_condition("S1", n_cells = 5000, output_dir = SIM_DIR)
run_simulation_scenario("S1", seu_S1, file.path(SIM_DIR, "S1"))

message("=== S2: Imbalanced composition, 5k cells ===")
seu_S2 <- simulate_condition("S2", n_cells = 5000, batch_imbalance = TRUE,
                              output_dir = SIM_DIR)
run_simulation_scenario("S2", seu_S2, file.path(SIM_DIR, "S2"))

message("=== S3: Balanced + categorical covariate, 5k cells ===")
seu_S3 <- simulate_condition("S3", n_cells = 5000,
                              include_covariate_categorical = TRUE,
                              output_dir = SIM_DIR)
run_simulation_scenario("S3", seu_S3, file.path(SIM_DIR, "S3"))

message("=== S4: Large balanced, 40k cells ===")
seu_S4 <- simulate_condition("S4", n_cells = 40000, output_dir = SIM_DIR)
run_simulation_scenario("S4", seu_S4, file.path(SIM_DIR, "S4"))

# ── scDesign3 scenario ────────────────────────────────────────────────────────
sce_path <- here("Data", "dataset1", "dataset1_sce.rds")
if (file.exists(sce_path)) {
  message("=== SD: scDesign3 realistic simulation ===")
  library(scDesign3)
  library(SingleCellExperiment)
  example_sce <- readRDS(sce_path)

  seu_SD <- simulate_scdesign3(example_sce,
                               percent_genes_modify    = 0.3,
                               batch_effect_mean       = 2,
                               celltype_effect_strength = 0.8)
  saveRDS(seu_SD, file.path(SIM_DIR, "simulation_SD_seurat.rds"))
  run_simulation_scenario("SD", seu_SD, file.path(SIM_DIR, "SD"))
} else {
  message("Skipping scDesign3 scenario: ", sce_path, " not found")
}

message("=== All simulation scenarios complete ===")
