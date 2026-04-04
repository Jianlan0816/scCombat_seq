library(here)

source(here("R", "core", "combat_nb_math.R"))
source(here("R", "core", "preprocessing.R"))
source(here("R", "core", "qc_checks.R"))
source(here("R", "core", "evaluation.R"))
source(here("R", "core", "combat_variants.R"))
source(here("R", "visualization", "plot_umap.R"))
source(here("R", "simulation", "simulate_splatter.R"))

# ── TEST 1: Simulation ────────────────────────────────────────────────────────
message("\n========== TEST 1: Simulation (S1) ==========")
seu_sim <- simulate_condition("S1", n_cells = 2000, n_genes = 500,
                              output_dir = here("results", "test_simulation"))

message("\n--- QC check ---")
check_data_characteristics(seu_sim)

message("\n--- Running ComBat_combo (pcseq only) ---")
res_sim <- ComBat_combo(seu_sim,
                        combat_pcseq = TRUE,
                        save_path    = here("results", "test_simulation", "umaps"))

message("\n--- Simulation metrics ---")
print(summarise_metrics(res_sim))

# ── TEST 2: Dataset 1 ─────────────────────────────────────────────────────────
message("\n========== TEST 2: Dataset 1 (DC, UC3) ==========")
source(here("analysis", "run_combat_variants.R"))

res_ds1 <- run_dataset_analysis(
  "dataset1",
  variants      = c("combat_seq", "combat_scseq", "combat_pcseq", "combat_ind"),
  run_benchmark = FALSE,
  subset        = 0.5,
  output_dir    = here("results")
)

message("\n--- Dataset 1 metrics ---")
print(summarise_metrics(res_ds1$combat))
