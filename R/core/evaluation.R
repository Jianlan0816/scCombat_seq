#' Shared Evaluation Utilities for scComBat-Seq
#'
#' Common scoring and plotting logic used by both combat_variants.R and benchmark_methods.R.
#' Source this file before sourcing either of those.

library(Seurat)
library(ggplot2)
library(cowplot)
library(kBET)
library(lisi)

#' Preprocess a Seurat object, plot UMAP, and compute batch mixing metrics
#'
#' Runs NormalizeData → FindVariableFeatures → ScaleData → RunPCA → RunUMAP
#' (unless \code{preprocess = FALSE}), then computes kBET and LISI scores.
#'
#' @param seu_obj Seurat object with 'batch' and 'celltype' in metadata
#' @param method Character label used in plot title and output filenames
#' @param preprocess Logical; if FALSE skip preprocessing (e.g. for Harmony/LIGER output)
#' @param reduction Name of reduction to use for plotting and metrics (default "pca")
#' @param save_path Directory to save UMAP PNG, or NULL to skip
#' @param k0_frac Fraction of cells used as neighbourhood size for kBET (default 0.05)
#' @return Named list: seurat (processed object), kBET_acceptance, batch_LISI, celltype_LISI
score_and_plot <- function(seu_obj, method,
                           preprocess = TRUE,
                           reduction  = "pca",
                           save_path  = NULL,
                           k0_frac    = 0.05) {
  if (preprocess) {
    seu_obj <- NormalizeData(seu_obj, verbose = FALSE) |>
      FindVariableFeatures(verbose = FALSE) |>
      ScaleData(verbose = FALSE) |>
      RunPCA(npcs = 20, reduction.name = "pca", verbose = FALSE) |>
      RunUMAP(dims = 1:20, umap.method = "uwot", metric = "cosine", verbose = FALSE)
    reduction <- "umap"
  }

  p1 <- DimPlot(seu_obj, reduction = reduction, group.by = "celltype") + NoLegend()
  p2 <- DimPlot(seu_obj, reduction = reduction, group.by = "batch") + NoLegend()
  title_row <- ggdraw() + draw_label(paste0(method, " — UMAP"), fontface = "bold")
  combined  <- plot_grid(title_row, plot_grid(p1, p2, ncol = 2), nrow = 2, rel_heights = c(0.12, 1))

  if (!is.null(save_path)) {
    dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
    ggsave(file.path(save_path, paste0(method, "_UMAP.png")),
           plot = combined, width = 10, height = 5, dpi = 150)
  }

  meta <- seu_obj@meta.data
  emb  <- Embeddings(seu_obj, reduction = reduction)

  # kBET
  k0         <- min(max(round(k0_frac * nrow(emb)), 10), 50)
  batch_vec  <- as.factor(meta$batch)
  kbet_res   <- kBET(emb, batch = batch_vec, k0 = k0, plot = FALSE)
  if (!is.null(kbet_res$results) && "kBET.pvalue.test" %in% colnames(kbet_res$results)) {
    kbet_acc <- mean(kbet_res$results$kBET.pvalue.test > 0.05, na.rm = TRUE)
  } else {
    kbet_acc <- NA
  }

  # LISI
  lisi_scores  <- compute_lisi(emb, meta, c("batch", "celltype"))
  batch_lisi   <- mean(lisi_scores$batch,    na.rm = TRUE)
  ct_lisi      <- mean(lisi_scores$celltype, na.rm = TRUE)

  list(
    seurat          = seu_obj,
    kBET_acceptance = kbet_acc,
    batch_LISI      = batch_lisi,
    celltype_LISI   = ct_lisi
  )
}

#' Summarise a results list as a data.frame of metrics
#'
#' @param results_list Named list returned by ComBat_combo or benchmark_batch_methods
#' @return data.frame with columns: method, kBET_acceptance, batch_LISI, celltype_LISI
summarise_metrics <- function(results_list) {
  rows <- lapply(names(results_list), function(nm) {
    r <- results_list[[nm]]
    data.frame(
      method          = nm,
      kBET_acceptance = r$kBET_acceptance,
      batch_LISI      = r$batch_LISI,
      celltype_LISI   = r$celltype_LISI,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
