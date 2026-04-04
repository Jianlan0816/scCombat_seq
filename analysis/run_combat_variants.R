#' Run scComBat-Seq variants on a single dataset
#'
#' Reads dataset configuration from datasets.yaml, loads and preprocesses the
#' data, applies all four ComBat variants plus competitor benchmarking, and
#' saves results to output_dir/dataset_id/.
#'
#' Usage (interactive):
#'   source("analysis/run_combat_variants.R")
#'   run_dataset_analysis("dataset1")
#'
#' Usage (command line):
#'   Rscript analysis/run_combat_variants.R --dataset dataset1 --subset 0.1

library(here)
library(yaml)
library(Seurat)
library(Matrix)

source(here("R", "core", "preprocessing.R"))
source(here("R", "core", "qc_checks.R"))
source(here("R", "core", "evaluation.R"))
source(here("R", "core", "combat_variants.R"))
source(here("R", "core", "benchmark_methods.R"))
source(here("R", "visualization", "plot_umap.R"))

# ── Data loaders ─────────────────────────────────────────────────────────────

.load_single_matrix <- function(cfg, root) {
  counts <- as.matrix(read.table(file.path(root, cfg$counts_file),
                                 header = TRUE, row.names = 1, sep = "\t",
                                 check.names = FALSE))
  meta   <- read.table(file.path(root, cfg$sample_file),
                       header = TRUE, row.names = 1, sep = "\t",
                       stringsAsFactors = FALSE)
  # Standardise column names
  if (!is.null(cfg$celltype_col) && cfg$celltype_col != "celltype" &&
      cfg$celltype_col %in% colnames(meta)) {
    meta$celltype <- meta[[cfg$celltype_col]]
  }
  if (!is.null(cfg$batch_col) && cfg$batch_col != "batch" &&
      cfg$batch_col %in% colnames(meta)) {
    meta$batch <- meta[[cfg$batch_col]]
  }
  list(counts = counts, meta = meta)
}

.load_split_batches <- function(cfg, root) {
  b1_counts <- as.matrix(read.table(file.path(root, cfg$batch1_counts),
                                    header = TRUE, row.names = 1, sep = "\t",
                                    check.names = FALSE))
  b2_counts <- as.matrix(read.table(file.path(root, cfg$batch2_counts),
                                    header = TRUE, row.names = 1, sep = "\t",
                                    check.names = FALSE))
  b1_meta   <- read.table(file.path(root, cfg$batch1_meta),
                           header = TRUE, row.names = 1, sep = "\t",
                           stringsAsFactors = FALSE)
  b2_meta   <- read.table(file.path(root, cfg$batch2_meta),
                           header = TRUE, row.names = 1, sep = "\t",
                           stringsAsFactors = FALSE)

  # Standardise
  for (m in list(b1_meta, b2_meta)) {
    if (!is.null(cfg$celltype_col) && cfg$celltype_col != "celltype" &&
        cfg$celltype_col %in% colnames(m)) m$celltype <- m[[cfg$celltype_col]]
    if (!"batch" %in% colnames(m)) m$batch <- NA
  }
  b1_meta$batch <- "Batch1"
  b2_meta$batch <- "Batch2"

  # Intersect genes
  common_genes <- intersect(rownames(b1_counts), rownames(b2_counts))
  counts <- cbind(b1_counts[common_genes, ], b2_counts[common_genes, ])
  meta   <- rbind(b1_meta, b2_meta)
  list(counts = counts, meta = meta)
}

.load_rds <- function(cfg, root) {
  counts <- readRDS(file.path(root, cfg$counts_file))
  if (!is.matrix(counts)) counts <- as.matrix(counts)
  meta   <- read.table(file.path(root, cfg$sample_file),
                       header = TRUE, row.names = 1, sep = "\t",
                       stringsAsFactors = FALSE)
  if (!is.null(cfg$celltype_col) && cfg$celltype_col != "celltype" &&
      cfg$celltype_col %in% colnames(meta)) {
    meta$celltype <- meta[[cfg$celltype_col]]
  }
  if (!is.null(cfg$batch_col) && cfg$batch_col != "batch" &&
      cfg$batch_col %in% colnames(meta)) {
    meta$batch <- meta[[cfg$batch_col]]
  }
  list(counts = counts, meta = meta)
}

# ── Main runner ───────────────────────────────────────────────────────────────

#' Run all scComBat-Seq variants + benchmarks for one dataset
#'
#' @param dataset_id   ID string matching an entry in datasets.yaml (e.g. "dataset1")
#' @param config_path  Path to datasets.yaml (default: analysis/config/datasets.yaml)
#' @param variants     Character vector of variants to run. Options:
#'                     "combat_seq", "combat_scseq", "combat_pcseq", "combat_ind"
#' @param run_benchmark Logical; also run Harmony/LIGER/Seurat comparison
#' @param subset       Fraction of cells to use (1.0 = all; use 0.1 for quick dev runs)
#' @param output_dir   Root output directory; results go into output_dir/dataset_id/
#' @return Invisibly returns list with $combat and $benchmark results
run_dataset_analysis <- function(dataset_id,
                                 config_path   = here("analysis", "config", "datasets.yaml"),
                                 variants      = c("combat_seq", "combat_scseq",
                                                   "combat_pcseq", "combat_ind"),
                                 run_benchmark = TRUE,
                                 subset        = 1.0,
                                 output_dir    = here("results")) {

  config   <- yaml::read_yaml(config_path)
  ds_list  <- config$datasets
  cfg      <- ds_list[[which(sapply(ds_list, `[[`, "id") == dataset_id)]]
  if (is.null(cfg)) stop("Dataset '", dataset_id, "' not found in ", config_path)

  message("\n========== ", cfg$name, " (", dataset_id, ") ==========")

  root <- here()
  raw  <- switch(cfg$data_format,
    single_matrix = .load_single_matrix(cfg, root),
    split_batches = .load_split_batches(cfg, root),
    rds           = .load_rds(cfg, root),
    stop("Unknown data_format: ", cfg$data_format)
  )

  # Filter
  filtered <- filter_data_mtx(
    raw$counts,
    is_filter_cells = isTRUE(cfg$filter_cells),
    min_genes       = cfg$min_genes  %||% 300,
    is_filter_genes = isTRUE(cfg$filter_genes),
    min_cells       = cfg$min_cells  %||% 10
  )

  # Build Seurat object
  common_cells <- intersect(colnames(filtered$counts), rownames(raw$meta))
  seu <- CreateSeuratObject(
    counts    = Matrix::Matrix(filtered$counts[, common_cells], sparse = TRUE),
    meta.data = raw$meta[common_cells, , drop = FALSE]
  )

  # QC check
  check_data_characteristics(seu)

  # Pre-compute PCA (required for combat_pcseq)
  seu <- NormalizeData(seu, verbose = FALSE) |>
    FindVariableFeatures(verbose = FALSE) |>
    ScaleData(verbose = FALSE) |>
    RunPCA(npcs = cfg$pca_dims %||% 20, verbose = FALSE) |>
    RunUMAP(dims = seq_len(cfg$pca_dims %||% 20), verbose = FALSE)

  out_dir <- file.path(output_dir, dataset_id)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # For log-normalised datasets (8, 9): only pcseq and raw ComBat apply
  if (cfg$input_type == "log_normalized") {
    message("[i] log_normalized input: only combat_pcseq is applicable")
    variants <- intersect(variants, "combat_pcseq")
  }

  # Run ComBat variants
  combat_res <- ComBat_combo(
    seu,
    subset        = subset,
    print_raw     = TRUE,
    combat_seq    = "combat_seq"   %in% variants,
    combat_scseq  = "combat_scseq" %in% variants,
    combat_pcseq  = "combat_pcseq" %in% variants,
    combat_ind    = "combat_ind"   %in% variants,
    save_path     = file.path(out_dir, "combat_umaps")
  )
  saveRDS(combat_res, file.path(out_dir, "combat_results.rds"))

  # Save metrics table
  metrics <- summarise_metrics(combat_res)
  write.csv(metrics, file.path(out_dir, "combat_metrics.csv"), row.names = FALSE)
  message("ComBat metrics:")
  print(metrics)

  # Run competitor benchmarks
  bench_res <- NULL
  if (run_benchmark) {
    bench_res <- benchmark_batch_methods(
      seu,
      subset      = subset,
      save_path   = file.path(out_dir, "benchmark_umaps")
    )
    saveRDS(bench_res, file.path(out_dir, "benchmark_results.rds"))
    bench_metrics <- summarise_metrics(bench_res)
    write.csv(bench_metrics, file.path(out_dir, "benchmark_metrics.csv"), row.names = FALSE)
    message("Benchmark metrics:")
    print(bench_metrics)
  }

  # Assemble figure
  all_results <- c(combat_res, bench_res)
  assemble_figure(all_results, save_path = file.path(out_dir, "figures"))

  invisible(list(combat = combat_res, benchmark = bench_res))
}

# ── NULL-coalescing operator ───────────────────────────────────────────────
`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Command-line interface ────────────────────────────────────────────────────
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  ds   <- grep("^--dataset",  args, value = TRUE)
  sub  <- grep("^--subset",   args, value = TRUE)

  dataset_id <- if (length(ds))  sub("^--dataset=?\\s*", "", ds[1])  else stop("--dataset required")
  subset_val <- if (length(sub)) as.numeric(sub("^--subset=?\\s*", "", sub[1])) else 1.0

  run_dataset_analysis(dataset_id, subset = subset_val)
}
