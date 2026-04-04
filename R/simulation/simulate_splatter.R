#' Splatter-based Simulation of Single-cell RNA-seq with Batch Effects
#'
#' Generates synthetic scRNA-seq data with controlled batch effects, optional
#' cell type imbalance, and optional categorical covariates.

library(splatter)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)

#' Simulate a single batch-effect scenario
#'
#' @param id                         Scenario identifier string (e.g. "S1"). Used for output filenames.
#' @param n_cells                    Total cells across both batches (default 5000)
#' @param n_genes                    Number of genes (default 1000)
#' @param n_celltypes                Number of cell types (default 3)
#' @param batch_imbalance            Logical; use group_probs_batch1/2 instead of equal proportions
#' @param group_probs_batch1         Cell type proportions in Batch1 (length must equal n_celltypes)
#' @param group_probs_batch2         Cell type proportions in Batch2
#' @param include_covariate_categorical Logical; add a "Young"/"Old" categorical covariate
#' @param output_dir                 Directory for saved .rds and UMAP .png files
#' @return Seurat object with batch, Group (→ celltype), and optional Covariate_cat in metadata
simulate_condition <- function(id,
                               n_cells         = 5000,
                               n_genes         = 1000,
                               n_celltypes     = 3,
                               batch_imbalance = FALSE,
                               group_probs_batch1          = c(0.7, 0.2, 0.1),
                               group_probs_batch2          = c(0.1, 0.1, 0.8),
                               include_covariate_categorical = FALSE,
                               output_dir      = ".") {

  stopifnot(length(group_probs_batch1) == n_celltypes)
  stopifnot(length(group_probs_batch2) == n_celltypes)

  set.seed(as.integer(substr(gsub("[^0-9]", "", id), 1, 5)) + 100)

  probs1 <- if (batch_imbalance) group_probs_batch1 else rep(1 / n_celltypes, n_celltypes)
  probs2 <- if (batch_imbalance) group_probs_batch2 else rep(1 / n_celltypes, n_celltypes)

  params <- newSplatParams()
  params <- setParam(params, "nGenes", n_genes)

  half <- n_cells %/% 2
  sim1 <- splatSimulateGroups(params, batchCells = half, group.prob = probs1, verbose = FALSE)
  sim2 <- splatSimulateGroups(params, batchCells = half, group.prob = probs2, verbose = FALSE)
  sim1$batch <- "Batch1"
  sim2$batch <- "Batch2"
  colnames(sim1) <- paste0("Batch1_", colnames(sim1))
  colnames(sim2) <- paste0("Batch2_", colnames(sim2))

  # Optional: categorical covariate
  if (include_covariate_categorical) {
    for (g in unique(sim1$Group)) {
      g1 <- which(sim1$Group == g)
      g2 <- which(sim2$Group == g)
      sim1$Covariate_cat[g1] <- sample(c("Young", "Old"), length(g1), replace = TRUE)
      sim2$Covariate_cat[g2] <- sample(c("Young", "Old"), length(g2), replace = TRUE)
    }
    sim1$Covariate_cat <- factor(sim1$Covariate_cat)
    sim2$Covariate_cat <- factor(sim2$Covariate_cat)
  }

  all_exprs <- cbind(counts(sim1), counts(sim2))
  all_meta  <- as.data.frame(rbind(colData(sim1), colData(sim2)))
  rownames(all_meta) <- colnames(all_exprs)

  # Inject DEGs for categorical covariate
  ground_truth_degs <- list()
  if (include_covariate_categorical) {
    for (g in unique(all_meta$Group)) {
      cells_g    <- rownames(all_meta)[all_meta$Group == g]
      young      <- cells_g[all_meta[cells_g, "Covariate_cat"] == "Young"]
      old        <- cells_g[all_meta[cells_g, "Covariate_cat"] == "Old"]
      if (length(young) > 0 && length(old) > 0) {
        de_genes <- sample(rownames(all_exprs), size = floor(0.05 * nrow(all_exprs)))
        ground_truth_degs[[g]]     <- de_genes
        all_exprs[de_genes, young] <- all_exprs[de_genes, young] * 1.2
        all_exprs[de_genes, old]   <- all_exprs[de_genes, old]   * 0.85
      }
    }
  }

  # Technical batch effect: 20% of genes in Batch2 get 1.2× uplift
  batch2_cells      <- rownames(all_meta)[all_meta$batch == "Batch2"]
  batch_effect_genes <- sample(rownames(all_exprs), size = floor(0.2 * nrow(all_exprs)))
  all_exprs[batch_effect_genes, batch2_cells] <-
    all_exprs[batch_effect_genes, batch2_cells] * 1.2

  all_exprs <- round(all_exprs)
  all_exprs <- Matrix::Matrix(all_exprs, sparse = TRUE)

  # Rename Group → celltype in metadata
  all_meta$celltype <- all_meta$Group

  seu <- CreateSeuratObject(counts = all_exprs, meta.data = all_meta)
  seu <- NormalizeData(seu, verbose = FALSE) |>
    FindVariableFeatures(verbose = FALSE) |>
    ScaleData(verbose = FALSE) |>
    RunPCA(npcs = 20, verbose = FALSE) |>
    RunUMAP(dims = 1:20, verbose = FALSE)

  # Batch effect strength: % variance in PC1 explained by batch
  pca_scores <- Embeddings(seu, "pca")[, 1, drop = FALSE]
  lm_res     <- summary(lm(pca_scores ~ as.factor(seu$batch)))
  message(sprintf("[%s] Batch effect (PC1 R²): %.1f%% (p = %.3g)",
                  id, 100 * lm_res$adj.r.squared,
                  lm_res$coefficients[2, 4]))

  # Plots
  plot_cols <- list(DimPlot(seu, group.by = "celltype") + ggtitle("Cell Type"),
                    DimPlot(seu, group.by = "batch")    + ggtitle("Batch"))
  if (include_covariate_categorical && "Covariate_cat" %in% colnames(seu@meta.data)) {
    plot_cols <- c(plot_cols,
                   list(DimPlot(seu, group.by = "Covariate_cat") + ggtitle("Covariate")))
  }
  g <- plot_grid(plotlist = plot_cols, ncol = length(plot_cols))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(output_dir, paste0("UMAP_", id, ".png")),
         g, width = 5 * length(plot_cols), height = 5)

  # Save outputs
  saveRDS(seu, file.path(output_dir, paste0("simulation_", id, "_seurat.rds")))
  saveRDS(batch_effect_genes, file.path(output_dir, paste0("simulation_", id, "_batch_genes.rds")))
  if (length(ground_truth_degs) > 0) {
    saveRDS(ground_truth_degs, file.path(output_dir, paste0("simulation_", id, "_degs.rds")))
  }
  message("Saved: simulation_", id, "_seurat.rds")

  invisible(seu)
}
