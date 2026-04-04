#' Benchmark Multiple Batch Correction Methods Against scComBat-Seq
#'
#' Applies Harmony, LIGER, and Seurat v5 RPCA integration to a Seurat object.
#' Returns corrected Seurat objects with kBET and LISI metrics for each method.
#'
#' @param seu           Seurat object with 'batch' and 'celltype' in metadata
#' @param subset        Fraction of cells to retain for benchmarking (default 1.0)
#' @param returnHarmony Run Harmony integration (default TRUE)
#' @param returnLiger   Run LIGER integration (default TRUE)
#' @param returnSeurat  Run Seurat v5 RPCA integration (default TRUE)
#' @param save_path     Directory for UMAP PNGs, or NULL to skip
#' @return Named list; each element has: seurat, kBET_acceptance, batch_LISI, celltype_LISI

library(Seurat)
library(harmony)
library(rliger)
library(ggplot2)
library(cowplot)

source(file.path(dirname(sys.frame(1)$ofile), "evaluation.R"))

benchmark_batch_methods <- function(seu,
                                    subset        = 1.0,
                                    returnHarmony = TRUE,
                                    returnLiger   = TRUE,
                                    returnSeurat  = TRUE,
                                    save_path     = NULL) {
  set.seed(123)
  if (subset < 1.0) {
    cells_to_use <- sample(Cells(seu), size = floor(subset * ncol(seu)))
    seu <- subset(seu, cells = cells_to_use)
  }

  seu$batch     <- as.factor(seu$batch)
  raw_celltype  <- setNames(seu$celltype, Cells(seu))
  raw_batch     <- setNames(seu$batch,    Cells(seu))
  has_cat       <- "Covariate_cat" %in% colnames(seu@meta.data)
  if (has_cat) raw_cat <- setNames(seu$Covariate_cat, Cells(seu))

  results <- list()

  # ── 1. Harmony ────────────────────────────────────────────────────────────
  if (returnHarmony) {
    message("--- Harmony ---")
    seu_harmony <- NormalizeData(seu, verbose = FALSE) |>
      FindVariableFeatures(verbose = FALSE) |>
      ScaleData(verbose = FALSE) |>
      RunPCA(npcs = 20, verbose = FALSE)
    seu_harmony <- RunHarmony(seu_harmony, group.by.vars = "batch", verbose = FALSE)
    seu_harmony <- RunUMAP(seu_harmony, reduction = "harmony", dims = 1:20, verbose = FALSE)
    results$harmony <- score_and_plot(seu_harmony, "harmony",
                                       preprocess = FALSE,
                                       reduction  = "harmony",
                                       save_path  = save_path)
  }

  # ── 2. LIGER ─────────────────────────────────────────────────────────────
  if (returnLiger) {
    message("--- LIGER ---")
    seu_list <- SplitObject(seu, split.by = "batch")
    rawList  <- lapply(seu_list, function(obj) {
      mat <- GetAssayData(obj, assay = "RNA", layer = "counts")
      if (inherits(mat, "DelayedMatrix")) as.matrix(mat) else mat
    })

    pbmcLiger <- rliger::createLiger(rawList) |>
      rliger::normalize() |>
      rliger::selectGenes() |>
      rliger::scaleNotCenter() |>
      rliger::runIntegration(k = 20) |>
      rliger::alignFactors(method = "centroidAlign")

    seu_liger <- rliger::ligerToSeurat(pbmcLiger)
    seu_liger <- RunUMAP(seu_liger, reduction = "inmf", dims = 1:20, verbose = FALSE)

    common_cells      <- intersect(colnames(seu_liger), names(raw_celltype))
    seu_liger         <- subset(seu_liger, cells = common_cells)
    seu_liger$celltype <- raw_celltype[common_cells]
    seu_liger$batch    <- raw_batch[common_cells]
    if (has_cat) seu_liger$Covariate_cat <- raw_cat[common_cells]

    results$liger <- score_and_plot(seu_liger, "liger",
                                     preprocess = FALSE,
                                     reduction  = "umap",
                                     save_path  = save_path)
  }

  # ── 3. Seurat v5 RPCA Integration ────────────────────────────────────────
  if (returnSeurat) {
    message("--- Seurat v5 RPCA ---")
    seu_list <- SplitObject(seu, split.by = "batch")
    seu_list <- lapply(seu_list, function(x) {
      x <- NormalizeData(x, verbose = FALSE)
      x <- FindVariableFeatures(x, verbose = FALSE)
      x
    })

    features <- SelectIntegrationFeatures(object.list = seu_list)
    seu_list <- lapply(seu_list, function(x) {
      x <- ScaleData(x, features = features, verbose = FALSE)
      x <- RunPCA(x, features = features,
                  npcs = min(20, ncol(x) - 1, length(features) - 1),
                  reduction.name = "pca", verbose = FALSE)
      x
    })

    anchors <- FindIntegrationAnchors(object.list = seu_list,
                                      anchor.features = features,
                                      reduction = "rpca", dims = 1:10)
    seu_integrated <- IntegrateData(anchorset = anchors, dims = 1:10, k.weight = 10)

    DefaultAssay(seu_integrated) <- "integrated"
    seu_integrated <- SetAssayData(seu_integrated, assay = "integrated", layer = "counts",
                                   new.data = GetAssayData(seu_integrated, layer = "data"))
    seu_integrated <- ScaleData(seu_integrated, verbose = FALSE)
    seu_integrated <- RunPCA(seu_integrated, npcs = 10, verbose = FALSE)
    seu_integrated <- RunUMAP(seu_integrated, dims = 1:10, verbose = FALSE)

    results$seurat5 <- score_and_plot(seu_integrated, "seurat5",
                                       preprocess = FALSE,
                                       reduction  = "umap",
                                       save_path  = save_path)
  }

  return(results)
}
