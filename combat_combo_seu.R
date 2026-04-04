#' Adjust for batch effects in single-cell RNA-seq data using multiple ComBat strategies
#'
#' ComBat_combo provides four correction options: ComBat_seq (batch only), ComBat_seq with cell type, 
#' ComBat on PCA embeddings, and per-cell type ComBat_seq.
#'
#' @param seu A Seurat object containing batches and cell type metadata. Assumes PCA and UMAP have been run.
#' @param subset Fraction of cells to retain for faster computation (default = 1.0)
#' @param print_raw Logical; whether to plot UMAP of uncorrected data (default = FALSE)
#' @param ComBat_seq Apply ComBat_seq using batch only (default = FALSE)
#' @param Combat_scseq Apply ComBat_seq using batch + cell type covariates (default = FALSE)
#' @param Combat_pcseq Apply ComBat on PCA embeddings using batch + cell type (default = FALSE)
#' @param Combat_ind Apply ComBat_seq separately for each cell type across batches (default = FALSE)

#' @param save_path Directory to save UMAP plots. Set NULL to skip saving (default = NULL)
#' @param score_clusters Whether to compute ARI/NMI comparing to raw cell types (default = TRUE)
#'
#' @return A named list with Seurat objects and ARI/NMI scores per method
#' 
library(dplyr)
library(Seurat)
library(sva)
library(umap)
library(uwot)
library(ggplot2)
library(cowplot)
library(aricode)
library(mclust)
library(patchwork) 
library(kBET)
library(clue)
library(future)
library(lisi)

ComBat_combo <- function(seu, subset = 1.0, print_raw = FALSE, 
                         combat_seq = FALSE, combat_scseq = FALSE, 
                         combat_pcseq = FALSE, combat_ind = FALSE,
                         save_path = NULL, score_clusters = TRUE) {
  
  set.seed(123)
  
  cells_to_use <- sample(Cells(seu), size = floor(subset * length(Cells(seu))))
  seu <- subset(seu, cells = cells_to_use)
  counts <- GetAssayData(seu, layer = "counts")
  metadata <- seu@meta.data
  metadata$celltype <- as.factor(metadata$celltype)
  results <- list()
  raw_celltype <- seu$celltype
  
  
  plot_and_score <- function(seu_obj, method, raw_celltype = NULL) {
    # Preprocessing
    seu_obj <- NormalizeData(seu_obj) %>% 
      FindVariableFeatures() %>% 
      ScaleData() %>% 
      RunPCA(npcs = 20, reduction.name = "pca", overwrite = TRUE) %>% 
      RunUMAP(dims = 1:20, umap.method = "uwot", metric = "cosine")
    
    # Plot UMAP
    p1 <- DimPlot(seu_obj, group.by = "celltype")
    p2 <- DimPlot(seu_obj, group.by = "batch")
    title <- ggdraw() + draw_label(paste0(method, " UMAP"), fontface = 'bold')
    plot_combined <- cowplot::plot_grid(title, cowplot::plot_grid(p1, p2, ncol = 2), rel_heights = c(0.1, 1), nrow = 2)
    
    if (!is.null(save_path)) {
      save_file <- file.path(save_path, paste0(method, "_UMAP.png"))
      print(paste("Saving plot to:", save_file))
      ggsave(
        filename = save_file,
        plot = plot_combined,
        width = 10,
        height = 5
      )
    }
    
    # Batch mixing metrics (kBET, LISI)
    meta <- seu_obj@meta.data
    emb <- Embeddings(seu_obj, reduction = "pca")
    
    # kBET
    # Safe k0 range
    k0 <- min(max(round(0.05 * nrow(emb)), 10), 50)
    
    # Ensure proper batch labels
    batch_vector <- as.factor(meta$batch)
    
    # Run kBET
    kbet_res <- kBET(emb, batch = batch_vector, k0 = k0, plot = FALSE)
    
    # Handle potential NA
    if (!is.null(kbet_res$results) && "kBET.pvalue.test" %in% colnames(kbet_res$results)) {
      kbet_acceptance <- mean(kbet_res$results$kBET.pvalue.test > 0.05, na.rm = TRUE)
    } else {
      kbet_acceptance <- NA
    }
    
    # LISI
    lisi_scores <- compute_lisi(emb, meta, c("batch","celltype"))
    batch_lisi <- mean(lisi_scores$batch)
    celltype_lisi <- mean(lisi_scores$celltype)
    
    return(list(
      seurat = seu_obj,
      kBET_acceptance = kbet_acceptance,
      batch_LISI = batch_lisi,
      celltype_LISI = celltype_lisi
    ))
  }
  
  if (print_raw) {
    results$raw = plot_and_score(seu, "Raw")
  }
  
  if (combat_seq) {
    combat_counts <- sva::ComBat_seq(as.matrix(counts), batch = as.factor(metadata$batch))
    seu_cb <- Seurat::CreateSeuratObject(counts = combat_counts, meta.data = metadata)
    results$combat_seq <- plot_and_score(seu_cb, "combat_seq")
  }
  
  if (combat_scseq) {
    if ("Covariate_cont" %in% colnames(sim@meta.data)){
      mod <- model.matrix(~ celltype + Covariate_cont, data = metadata)}
    else{
      mod <- model.matrix(~ celltype, data = metadata)}
    combat_counts <- sva::ComBat_seq(as.matrix(counts), batch = as.factor(metadata$batch), covar_mod = mod)
    seu_cbsc <- CreateSeuratObject(counts = combat_counts, meta.data = metadata)
    results$combat_scseq <- plot_and_score(seu_cbsc, "combat_scseq")
  }
  
  if (combat_pcseq) {
    pcs <- seu@reductions$pca@cell.embeddings
    if ("Covariate_cont" %in% colnames(sim@meta.data)){
      mod <- model.matrix(~ celltype + Covariate_cont, data = metadata)}
    else{
      mod <- model.matrix(~ celltype, data = metadata)}
    corrected_pcs <- sva::ComBat(dat = t(pcs), batch = metadata$batch, mod = mod)
    corrected_pcs <- t(corrected_pcs)
    
    umap_coords <- as.data.frame(umap::umap(corrected_pcs, a = 1, b = 1)$layout)
    umap_coords$celltype <- seu$celltype  # ✅ Matches corrected_pcs
    umap_coords$batch <- seu$batch
    
    # plots
    p1 <- ggplot(umap_coords, aes(x = V1, y = V2, color = celltype)) +
      geom_point(size = 1) + theme_classic() + ggtitle("Cell Type")
    p2 <- ggplot(umap_coords, aes(x = V1, y = V2, color = batch)) +
      geom_point(size = 1) + theme_classic() + ggtitle("Batch")
    title <- ggdraw() + draw_label(paste0('combat_pcseq', " UMAP"), fontface = 'bold')
    plot_combined <- cowplot::plot_grid(title, cowplot::plot_grid(p1, p2, ncol = 2), rel_heights = c(0.1, 1), nrow = 2)
    
    if (!is.null(save_path)) {
      save_file <- file.path(save_path, paste0("combat_pcseq", "_UMAP.png"))
      print(paste("Saving plot to:", save_file))
      
      ggsave(
        filename = save_file,
        plot = plot_combined,   # Explicit plot object
        width = 10,
        height = 5
      )
    }
  }
  
  if (combat_ind) {
    all_celltypes <- unique(metadata$celltype)
    corrected_data <- data.frame(row.names = rownames(counts))
    combined_metadata <- NULL
    
    for (ct in all_celltypes) {
      idx <- metadata$celltype == ct
      data_ct <- counts[, idx]
      meta_ct <- metadata[idx, ]
      if (length(unique(meta_ct$batch)) > 1) {
        adjusted <- ComBat_seq(as.matrix(data_ct), batch = as.factor(meta_ct$batch))
      } else {
        adjusted <- data_ct
      }
      corrected_data <- cbind(corrected_data, data.frame(adjusted))
      combined_metadata <- rbind(combined_metadata, meta_ct)
    }
    
    colnames(corrected_data) <- rownames(combined_metadata)
    seu_cb_ind <- CreateSeuratObject(counts = corrected_data, meta.data = combined_metadata)
    results$combat_ind <- plot_and_score(seu_cb_ind, "combat_ind")
  }
  return(results)
}

