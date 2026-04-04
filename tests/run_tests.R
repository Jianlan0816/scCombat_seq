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

message("\n--- Dataset 1 metrics (standard) ---")
print(summarise_metrics(res_ds1$combat))

# ── TEST 3: Adaptive covariate on Dataset 1 ───────────────────────────────────
message("\n========== TEST 3: Adaptive covariate (Option A) on Dataset 1 ==========")

# Load dataset1 directly (same as run_dataset_analysis does internally)
cfg    <- yaml::read_yaml(here("analysis","config","datasets.yaml"))$datasets[[1]]
counts <- as.matrix(read.table(here(cfg$counts_file), header=TRUE, row.names=1,
                               sep="\t", check.names=FALSE))
meta   <- read.table(here(cfg$sample_file), header=TRUE, row.names=1,
                     sep="\t", stringsAsFactors=FALSE)
filt   <- filter_data_mtx(counts, is_filter_cells=TRUE, min_genes=300,
                           is_filter_genes=TRUE, min_cells=10)
common <- intersect(colnames(filt$counts), rownames(meta))
seu_ds1 <- Seurat::CreateSeuratObject(
  counts    = Matrix::Matrix(filt$counts[, common], sparse=TRUE),
  meta.data = meta[common, , drop=FALSE]
)
seu_ds1 <- Seurat::NormalizeData(seu_ds1, verbose=FALSE) |>
  Seurat::FindVariableFeatures(verbose=FALSE) |>
  Seurat::ScaleData(verbose=FALSE) |>
  Seurat::RunPCA(npcs=20, verbose=FALSE) |>
  Seurat::RunUMAP(dims=1:20, verbose=FALSE)

message("\n--- combat_scseq ORIGINAL (all celltypes in covariate, ignores confounding) ---")
res_orig  <- ComBat_combo(seu_ds1, combat_scseq = TRUE,
                          use_adaptive_covariate = FALSE, save_path = NULL)

message("\n--- combat_scseq ADAPTIVE (drops CD141 + CD1C from covariate) ---")
res_adapt <- ComBat_combo(seu_ds1, combat_scseq = TRUE,
                          use_adaptive_covariate = TRUE, save_path = NULL)

message("\n--- combat_seq (batch only, for reference) ---")
res_seq   <- ComBat_combo(seu_ds1, combat_seq = TRUE, save_path = NULL)

message("\n--- Comparison ---")
comp <- rbind(
  cbind(variant = "combat_seq      ", summarise_metrics(res_seq)),
  cbind(variant = "scseq_original  ", summarise_metrics(res_orig)),
  cbind(variant = "scseq_adaptive  ", summarise_metrics(res_adapt))
)
print(comp)
