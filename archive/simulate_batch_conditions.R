library(splatter)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(patchwork)
library(scater)

simulate_condition <- function(id,
                               n_cells = 5000,
                               n_genes = 1000,
                               n_celltypes = 3,
                               batch_imbalance = FALSE,
                               raw = TRUE,
                               include_covariate_categorical = FALSE,
                               group_probs_batch1 = c(0.7, 0.2, 0.1),
                               group_probs_batch2 = c(0.1, 0.1, 0.8)) {
  
  set.seed(as.integer(substr(gsub("[^0-9]", "", id), 1, 5)) + 100)
  
  stopifnot(length(group_probs_batch1) == n_celltypes)
  stopifnot(length(group_probs_batch2) == n_celltypes)
  
  batch_cells <- rep(n_cells / 2, 2)
  params <- newSplatParams()
  params <- setParam(params, "nGenes", n_genes)
  
  probs1 <- if (batch_imbalance) group_probs_batch1 else rep(1 / n_celltypes, n_celltypes)
  probs2 <- if (batch_imbalance) group_probs_batch2 else rep(1 / n_celltypes, n_celltypes)
  
  sim1 <- splatSimulateGroups(params, batchCells = batch_cells[1], group.prob = probs1, verbose = FALSE)
  sim2 <- splatSimulateGroups(params, batchCells = batch_cells[2], group.prob = probs2, verbose = FALSE)
  sim1$batch <- "Batch1"
  sim2$batch <- "Batch2"
  
  colnames(sim1) <- paste0("Batch1_", colnames(sim1))
  colnames(sim2) <- paste0("Batch2_", colnames(sim2))
  
  # Optional: categorical covariate
  if (include_covariate_categorical) {
    for (g in unique(sim1$Group)) {
      g_idx1 <- which(sim1$Group == g)
      g_idx2 <- which(sim2$Group == g)
      
      sim1$Covariate_cat[g_idx1] <- sample(c("Young", "Old"), length(g_idx1), replace = TRUE)
      sim2$Covariate_cat[g_idx2] <- sample(c("Young", "Old"), length(g_idx2), replace = TRUE)
    }
    
    sim1$Covariate_cat <- factor(sim1$Covariate_cat)
    sim2$Covariate_cat <- factor(sim2$Covariate_cat)
  }
  
  exprs1 <- counts(sim1)
  exprs2 <- counts(sim2)
  all_exprs <- cbind(exprs1, exprs2)
  all_meta <- rbind(colData(sim1), colData(sim2))
  rownames(all_meta) <- colnames(all_exprs)
  
  ground_truth_degs <- list()
  
  # Optional: inject biological DEGs (if Covariate_cat exists)
  if (include_covariate_categorical) {
    for (g in unique(all_meta$Group)) {
      cells_g <- rownames(all_meta)[all_meta$Group == g]
      young_cells <- cells_g[all_meta[cells_g, "Covariate_cat"] == "Young"]
      old_cells <- cells_g[all_meta[cells_g, "Covariate_cat"] == "Old"]
      
      if (length(young_cells) > 0 && length(old_cells) > 0) {
        de_genes <- sample(rownames(all_exprs), size = floor(0.05 * nrow(all_exprs)))
        ground_truth_degs[[g]] <- de_genes
        all_exprs[de_genes, young_cells] <- all_exprs[de_genes, young_cells] * 1.2
        all_exprs[de_genes, old_cells] <- all_exprs[de_genes, old_cells] * 0.85
      }
    }
  }
  
  # Always apply batch effect (technical)
  batch2_cells <- rownames(all_meta)[all_meta$batch == "Batch2"]
  batch_effect_genes <- sample(rownames(all_exprs), size = floor(0.2 * nrow(all_exprs)))
  all_exprs[batch_effect_genes, batch2_cells] <- all_exprs[batch_effect_genes, batch2_cells] * 1.2
  
  all_exprs <- round(all_exprs)
  all_exprs <- Matrix::Matrix(all_exprs, sparse = TRUE)
  
  sim_combined <- SingleCellExperiment(assays = list(counts = all_exprs), colData = all_meta)
  
  expr_matrix <- counts(sim_combined)
  seu <- CreateSeuratObject(counts = expr_matrix, meta.data = as.data.frame(colData(sim_combined)))
  
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, npcs = 20)
  seu <- RunUMAP(seu, dims = 1:20)
  
  plots <- list()
  
  if ("pca" %in% names(seu@reductions) || "umap" %in% names(seu@reductions)) {
    plots[[1]] <- tryCatch(DimPlot(seu, group.by = "Group", label = TRUE) + ggtitle("UMAP by Cell Type"), error = function(e) NULL)
    plots[[2]] <- tryCatch(DimPlot(seu, group.by = "batch") + ggtitle("UMAP by Batch"), error = function(e) NULL)
    if (include_covariate_categorical) {
      plots[[3]] <- tryCatch(DimPlot(seu, group.by = "Covariate_cat") + ggtitle("UMAP by Covariate_cat"), error = function(e) NULL)
    }
    
    g <- cowplot::plot_grid(plotlist = plots, ncol = length(plots))
    ggsave(paste0("UMAP_", id, ".png"), g, width = 5 * length(plots), height = 5)
  } else {
    message("⚠️ Skipping DimPlot: no dimensional reduction available (raw = TRUE)")
  }
  
  g <- cowplot::plot_grid(plotlist = plots, ncol = length(plots))
  ggsave(paste0("UMAP_", id, ".png"), g, width = 5 * length(plots), height = 5)
  
  saveRDS(seu, file = paste0("simulation_", id, "_seurat.rds"))
  saveRDS(batch_effect_genes, file = paste0("simulation_", id, "_batch_effect_genes.rds"))
  if (length(ground_truth_degs) > 0) {
    saveRDS(ground_truth_degs, file = paste0("simulation_", id, "_ground_truth_degs.rds"))
  }
  
  message("✅ Saved: simulation_", id, "_seurat.rds and UMAP_", id, ".png")
  
  if ("pca" %in% names(seu@reductions)) {
    pca_scores <- Embeddings(seu, reduction = "pca")[, 1, drop = FALSE]
    lm_result <- summary(lm(pca_scores ~ as.factor(seu$batch)))
    batch_r2 <- lm_result$adj.r.squared
    batch_pval <- lm_result$coefficients[2, 4]
    message(sprintf("📊 Batch effect strength (PC1 variance explained): %.2f%% (p = %.3g)",
                    100 * batch_r2, batch_pval))
  }
}



# Run 10 scenarios
# Balanced scenarios
simulate_condition("S1", 
                   raw = TRUE,  
                   n_cells = 5000,  
                   batch_imbalance = FALSE)

simulate_condition("S3", 
                   raw = TRUE,  
                   n_cells = 40000, 
                   batch_imbalance = FALSE, 
                   include_celltype = TRUE)
#simulate_condition("S4", raw = TRUE,  n_cells = 5000,  batch_imbalance = FALSE, include_celltype = FALSE, include_covariate = FALSE)

simulate_condition("S4", raw = FALSE, return_normalized = TRUE)
#simulate_condition("S6", raw = FALSE, return_normalized = TRUE, include_celltype = FALSE, include_covariate = FALSE)
#simulate_condition("S5",
#                   raw = TRUE,
#                   n_cells = 5000,
#                   batch_imbalance = FALSE,
#                   include_celltype = TRUE,
#                   include_covariate_continuous = TRUE)
#simulate_condition("S7", raw = TRUE,  n_cells = 40000, batch_imbalance = FALSE, include_celltype = FALSE, include_covariate = FALSE)
#simulate_condition("S5",raw = FALSE, n_cells = 40000, batch_imbalance = FALSE, include_celltype = TRUE,  include_covariate = TRUE, return_normalized = TRUE)
simulate_condition("S3",
                   raw = TRUE,
                   n_cells = 5000,
                   batch_imbalance = FALSE,
                   include_covariate_categorical = TRUE)
# Imbalanced scenarios (S2, S8) using special function
simulate_condition("S2", n_cells = 5000, batch_imbalance = TRUE)
simulate_condition("S6", n_cells = 40000, batch_imbalance = TRUE)
