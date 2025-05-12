#' Benchmark Multiple Batch Correction Methods for Single-cell RNA-seq
#'
#' Applies Harmony, LIGER, Seurat v3 Integration, and scMerge to a Seurat object.
#' Returns corrected Seurat objects, along with ARI and NMI scores for each method.
#'
#' @param seu A Seurat object with raw counts, and metadata containing \"batch\" and \"celltype\".
#' @param subset A float between 0 and 1 to subsample cells for benchmarking. Default is 1 (all cells).
#' @param save_path Optional path to save UMAP plots. Default is NULL.
#' @param score_clusters Logical; if TRUE, computes ARI/NMI using k-means on PCA. Default is TRUE.
#'
#' @return A named list where each element is a list with:
#'   \item{seurat}{Corrected Seurat object}
#'   \item{ARI}{Adjusted Rand Index against celltype}
#'   \item{NMI}{Normalized Mutual Information against celltype}
#'
#' @import Seurat
#' @import harmony
#' @import rliger
#' @import aricode
#' @import mclust
#' @import scMerge
#' @import ggplot2
#' @import cowplot
#' @export
#' 
library(Seurat)
library(harmony)
library(rliger)
library(aricode)
library(mclust)
#library(scMerge)
library(cowplot)
library(ggplot2)
library(uwot)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(Seurat)
library(cowplot)
library(lisi)
library(kBET)
library(ggplot2)

benchmark_batch_methods <- function(seu, subset = 1.0, 
                                    save_path = NULL, score_clusters = TRUE) {
  set.seed(123)
  cells_to_use <- sample(Cells(seu), size = floor(subset * length(Cells(seu))))
  seu <- subset(seu, cells = cells_to_use)
  seu$batch <- as.factor(seu$batch)  # ensure factor
  #raw_celltype <- seu$celltype
  # 👇 Preserve cell names for clean metadata restoration later
  raw_celltype <- setNames(seu$celltype, Cells(seu))
  raw_batch <- setNames(seu$batch, Cells(seu))
  
  results <- list()
  
  plot_and_score <- function(seu_obj, method) {
    assay_now <- DefaultAssay(seu_obj)
    
    # Preprocessing depending on method
    if (method %in% c("seurat5")) {
      if (!"pca" %in% names(seu_obj@reductions)) {
        message("Running NormalizeData, FindVariableFeatures, and ScaleData for ", method)
        seu_obj <- NormalizeData(seu_obj)
        seu_obj <- FindVariableFeatures(seu_obj)
        seu_obj <- ScaleData(seu_obj)
        seu_obj <- RunPCA(seu_obj)
        seu_obj <- RunUMAP(seu_obj, dims = 1:20)
      }
      reduction_name <- "umap"
    } else if (method == "liger") {
      message("Skipping NormalizeData and FindVariableFeatures for LIGER")
      seu_obj <- NormalizeData(seu_obj)
      seu_obj <- ScaleData(seu_obj)
      seu_obj <- RunPCA(seu_obj)
      seu_obj <- RunUMAP(seu_obj, dims = 1:20)
      reduction_name <- "umap"
    } else if (method == "harmony") {
      if (!"harmony" %in% names(seu_obj@reductions)) {
        stop("Harmony reduction not found in seu_obj. Please run Harmony integration first.")
      }
      reduction_name <- "harmony"
      # Run UMAP on harmony
      if (!"umap" %in% names(seu_obj@reductions)) {
        seu_obj <- RunUMAP(seu_obj, reduction = "harmony", dims = 1:20)
      }
    } else {
      reduction_name <- if ("umap" %in% names(seu_obj@reductions)) "umap" else "pca"
    }
    
    # UMAP plots
    p1 <- DimPlot(seu_obj, reduction = reduction_name, group.by = "celltype")
    p2 <- DimPlot(seu_obj, reduction = reduction_name, group.by = "batch")
    title <- ggdraw() + draw_label(paste0(method, " UMAP"), fontface = 'bold')
    plot_combined <- plot_grid(title, plot_grid(p1, p2, ncol = 2), rel_heights = c(0.1, 1), nrow = 2)
    
    if (!is.null(save_path)) {
      save_file <- file.path(save_path, paste0(method, "_UMAP.png"))
      message("Saving plot to: ", save_file)
      ggsave(filename = save_file, plot = plot_combined, width = 10, height = 5)
    }
    
    # Get embeddings and metadata
    embedding <- Embeddings(seu_obj, reduction = reduction_name)
    metadata <- seu_obj@meta.data
    
    # --- ARI & NMI ---
    # Create clusters if not already
    if (!"seurat_clusters" %in% colnames(seu_obj@meta.data)) {
      # Determine max available PCs
      n_pcs <- ncol(Embeddings(seu_obj, reduction = "pca"))
      dims_to_use <- 1:min(20, n_pcs)
      
      if (length(dims_to_use) < 1) {
        stop("PCA reduction exists but no usable dimensions found.")
      }
      
      seu_obj <- FindNeighbors(seu_obj, reduction = "pca", dims = dims_to_use)
      seu_obj <- FindClusters(seu_obj)
    }
    
    ARI <- adjustedRandIndex(seu_obj$seurat_clusters, seu_obj$celltype)
    NMI_score <- NMI(seu_obj$seurat_clusters, seu_obj$celltype)
    # Return everything
    return(list(
      seu_obj = seu_obj,
      ARI = ARI,
      NMI = NMI_score
    ))
  }
  
  
  # ── 1. Harmony ──
  print("-------------Harmony----------------")
  seu_harmony <- NormalizeData(seu) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(npcs = 20)
  seu_harmony <- RunHarmony(seu_harmony, group.by.vars = "batch")
  seu_harmony <- RunUMAP(seu_harmony, reduction = "harmony", dims = 1:20)
  results$harmony <- plot_and_score(seu_harmony, "harmony")
  
  # ── 2. LIGER ──
  # Split and prepare
  print("--------------Liger------------")
  
  seu_list <- SplitObject(seu, split.by = "batch")
  rawList <- lapply(seu_list, function(obj) {
    GetAssayData(obj, assay = "RNA", slot = "counts")
  })
  # 1. Create and normalize LIGER object
  pbmcLiger <- createLiger(rawList, organism = "human") %>%
    rliger::normalize() %>%
    rliger::selectGenes() %>%
    rliger::scaleNotCenter()
  
  # 2. Run integration and alignment
  pbmcLiger <- runIntegration(pbmcLiger, k = 50)
  pbmcLiger <- alignFactors(pbmcLiger, method = "centroidAlign")
  pbmcLiger <- runUMAP(pbmcLiger)
  
  # Convert directly
  seu_liger <- ligerToSeurat(pbmcLiger)
  print(unique(seu_liger$dataset))
  # --- Restore benchmark metadata ---
  # Find overlapping cell names
  # Extract dataset (batch) assignments used as prefixes
  batch_prefixes <- unique(seu_liger$dataset)
  # Escape ( and ) in batch names
  batch_prefixes_fixed <- gsub("([\\(\\)])", "\\\\\\1", batch_prefixes)
  # Now build the regex pattern
  pattern <- paste0("^(", paste0(batch_prefixes_fixed, collapse = "|"), ")_")
  cat("fixed pattern =", pattern, "\n")
  # Apply the corrected pattern to clean cell names
  colnames(seu_liger) <- sub(pattern, "", colnames(seu_liger))
  # Check if they match now
  print(colnames(seu_liger)[1:5])
  print(names(raw_celltype)[1:5])
  # Restore clean metadata using stripped barcodes
  common_cells <- intersect(colnames(seu_liger), names(raw_celltype))
  seu_liger <- subset(seu_liger, cells = common_cells)
  seu_liger$celltype <- raw_celltype[common_cells]
  seu_liger$batch <- raw_batch[common_cells]
  print(seu_liger)
  # Then run your scoring function
  results$liger <- plot_and_score(seu_liger, "liger")
  
  # ── 3. Seurat v5 Integration ──
  # Split object by batch
  seu_list <- SplitObject(seu, split.by = "batch")
  
  # Standard preprocessing for each object
  seu_list <- lapply(seu_list, function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x)
    x
  })
  
  # Identify integration features
  features <- SelectIntegrationFeatures(object.list = seu_list)
  
  # Prepare for integration using RPCA
  seu_list <- lapply(seu_list, function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(
      x, 
      features = features, 
      npcs = min(20, ncol(x) - 1, length(features) - 1), 
      verbose = FALSE,
      reduction.name = "pca"  # ✅ Required for RPCA
    )
    x
  })
  
  anchors <- FindIntegrationAnchors(
    object.list = seu_list,
    anchor.features = features,
    reduction = "rpca",
    dims = 1:10
  )
  
  seu_integrated <- IntegrateData(
    anchorset = anchors,
    dims = 1:10,
    k.weight = 10
  )
  
  # 🔁 Work on integrated assay only — DO NOT re-normalize
  DefaultAssay(seu_integrated) <- "integrated"
  # Fix for layer assignment in Seurat v5
  seu_integrated <- SetAssayData(
    object = seu_integrated,
    assay = "integrated",
    layer = "counts",
    new.data = GetAssayData(seu_integrated, layer = "data")
  )
  
  seu_integrated <- ScaleData(seu_integrated, verbose = FALSE)
  seu_integrated <- RunPCA(seu_integrated, npcs = 10)
  seu_integrated <- RunUMAP(seu_integrated, dims = 1:10)
  
  results$seurat5 <- plot_and_score(seu_integrated, "seurat5")
  
  return(results)
}
