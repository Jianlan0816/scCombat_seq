# Batch Correction Benchmarking for Single-Cell RNA-seq

This repository contains R scripts to simulate single-cell datasets with batch effects and covariates, apply various batch correction methods (including the ComBat family and integration-based methods like Harmony, LIGER, and Seurat), and evaluate performance using UMAP, kBET, and LISI metrics.

## 📁 Repository Structure

| File | Description |
|------|-------------|
| `simulate_batch_conditions.R` | Simulates scRNA-seq data with configurable batch structure, covariates, and expression effects. |
| `combat_combo_seu.R` | Applies ComBat-based batch correction strategies directly to raw counts or embeddings in Seurat objects. |
| `Benchmark_batch_methods.R` | Applies Harmony, LIGER, and Seurat v5 integration workflows and evaluates using UMAP, kBET, and LISI. |
| `pre_check.R` | Performs quality control on Seurat objects and checks for key characteristics such as raw counts, metadata, and library size across batches. |

## 🚀 Usage

### 1. Simulate Datasets

```r
source("simulate_batch_conditions.R")
simulate_condition(
  id = "S1",
  n_cells = 5000,
  batch_imbalance = FALSE,
  include_covariate_categorical = TRUE
)
