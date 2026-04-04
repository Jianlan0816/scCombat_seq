#' scComBat-Seq: Four ComBat Correction Variants for Single-Cell RNA-seq
#'
#' Applies up to four correction strategies to a Seurat object and returns
#' corrected objects with batch mixing metrics.
#'
#' @param seu         Seurat object. Must have 'batch' and 'celltype' in metadata
#'                    and a pre-computed PCA reduction (used by combat_pcseq).
#' @param subset      Fraction of cells to retain (default 1.0; use < 1 for faster dev runs)
#' @param print_raw   Logical; include uncorrected baseline in results (default FALSE)
#' @param combat_seq  Apply ComBat-Seq using batch only
#' @param combat_scseq Apply ComBat-Seq with batch + celltype (+ optional Covariate_cont)
#' @param combat_pcseq Apply ComBat to PCA embeddings with batch + celltype
#' @param combat_ind  Apply ComBat-Seq per celltype then merge
#' @param save_path   Directory for UMAP PNGs, or NULL to skip
#' @param score_clusters Logical; compute kBET and LISI (default TRUE)
#' @return Named list; each element has: seurat, kBET_acceptance, batch_LISI, celltype_LISI

library(Seurat)
library(sva)
library(umap)
library(dplyr)
library(ggplot2)
library(cowplot)

source(file.path(dirname(sys.frame(1)$ofile), "evaluation.R"))

ComBat_combo <- function(seu,
                         subset        = 1.0,
                         print_raw     = FALSE,
                         combat_seq    = FALSE,
                         combat_scseq  = FALSE,
                         combat_pcseq  = FALSE,
                         combat_ind    = FALSE,
                         save_path     = NULL,
                         score_clusters = TRUE) {

  set.seed(123)
  if (subset < 1.0) {
    cells_to_use <- sample(Cells(seu), size = floor(subset * ncol(seu)))
    seu <- subset(seu, cells = cells_to_use)
  }

  counts   <- GetAssayData(seu, layer = "counts")
  metadata <- seu@meta.data
  metadata$celltype <- as.factor(metadata$celltype)
  results  <- list()

  # ── Helper: build covariate model matrix ─────────────────────────────────
  make_mod <- function(meta) {
    if ("Covariate_cont" %in% colnames(meta)) {
      model.matrix(~ celltype + Covariate_cont, data = meta)
    } else {
      model.matrix(~ celltype, data = meta)
    }
  }

  # ── Raw baseline ──────────────────────────────────────────────────────────
  if (print_raw) {
    results$raw <- score_and_plot(seu, "raw", preprocess = TRUE, save_path = save_path)
  }

  # ── Variant 1: ComBat-Seq, batch only ────────────────────────────────────
  if (combat_seq) {
    corrected <- sva::ComBat_seq(as.matrix(counts),
                                 batch = as.factor(metadata$batch))
    seu_cb <- CreateSeuratObject(counts = corrected, meta.data = metadata)
    results$combat_seq <- score_and_plot(seu_cb, "combat_seq",
                                         preprocess = TRUE, save_path = save_path)
  }

  # ── Variant 2: ComBat-Seq, batch + celltype (+ optional covariate) ───────
  if (combat_scseq) {
    mod      <- make_mod(metadata)
    corrected <- sva::ComBat_seq(as.matrix(counts),
                                 batch     = as.factor(metadata$batch),
                                 covar_mod = mod)
    seu_cbsc <- CreateSeuratObject(counts = corrected, meta.data = metadata)
    results$combat_scseq <- score_and_plot(seu_cbsc, "combat_scseq",
                                            preprocess = TRUE, save_path = save_path)
  }

  # ── Variant 3: ComBat on PCA embeddings ──────────────────────────────────
  if (combat_pcseq) {
    if (!"pca" %in% names(seu@reductions)) {
      stop("combat_pcseq requires a pre-computed PCA. Run RunPCA() on seu first.")
    }
    pcs      <- seu@reductions$pca@cell.embeddings
    mod      <- make_mod(metadata)
    corr_pcs <- sva::ComBat(dat = t(pcs),
                             batch = metadata$batch,
                             mod   = mod)
    corr_pcs <- t(corr_pcs)

    umap_coords      <- as.data.frame(umap::umap(corr_pcs, a = 1, b = 1)$layout)
    colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
    umap_coords$celltype  <- metadata$celltype
    umap_coords$batch     <- metadata$batch

    p1 <- ggplot(umap_coords, aes(UMAP_1, UMAP_2, color = celltype)) +
      geom_point(size = 0.5) + theme_classic() + ggtitle("Cell Type") + NoLegend()
    p2 <- ggplot(umap_coords, aes(UMAP_1, UMAP_2, color = batch)) +
      geom_point(size = 0.5) + theme_classic() + ggtitle("Batch") + NoLegend()
    title_row <- ggdraw() + draw_label("combat_pcseq — UMAP", fontface = "bold")
    combined  <- plot_grid(title_row, plot_grid(p1, p2, ncol = 2), nrow = 2, rel_heights = c(0.12, 1))

    if (!is.null(save_path)) {
      dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
      ggsave(file.path(save_path, "combat_pcseq_UMAP.png"),
             plot = combined, width = 10, height = 5, dpi = 150)
    }

    # Metrics computed on corrected PCs (not UMAP)
    k0       <- min(max(round(0.05 * nrow(corr_pcs)), 10), 50)
    kbet_res <- kBET::kBET(corr_pcs, batch = as.factor(metadata$batch), k0 = k0, plot = FALSE)
    if (!is.null(kbet_res$results) && "kBET.pvalue.test" %in% colnames(kbet_res$results)) {
      kbet_acc <- mean(kbet_res$results$kBET.pvalue.test > 0.05, na.rm = TRUE)
    } else {
      kbet_acc <- NA
    }
    lisi_scores <- lisi::compute_lisi(corr_pcs, metadata, c("batch", "celltype"))

    results$combat_pcseq <- list(
      seurat          = NULL,   # no Seurat object for PC-space variant
      kBET_acceptance = kbet_acc,
      batch_LISI      = mean(lisi_scores$batch,    na.rm = TRUE),
      celltype_LISI   = mean(lisi_scores$celltype, na.rm = TRUE),
      corrected_pcs   = corr_pcs,
      umap_coords     = umap_coords
    )
  }

  # ── Variant 4: ComBat-Seq applied per celltype, then merged ──────────────
  if (combat_ind) {
    all_celltypes    <- unique(metadata$celltype)
    corrected_list   <- list()
    combined_metadata <- NULL

    for (ct in all_celltypes) {
      idx     <- metadata$celltype == ct
      data_ct <- counts[, idx, drop = FALSE]
      meta_ct <- metadata[idx, , drop = FALSE]

      if (length(unique(meta_ct$batch)) > 1) {
        adjusted <- sva::ComBat_seq(as.matrix(data_ct),
                                    batch = as.factor(meta_ct$batch))
      } else {
        adjusted <- as.matrix(data_ct)
      }
      corrected_list[[ct]]  <- adjusted
      combined_metadata     <- rbind(combined_metadata, meta_ct)
    }

    corrected_data <- do.call(cbind, corrected_list)
    # Restore original cell order
    corrected_data <- corrected_data[, rownames(combined_metadata)]

    seu_cb_ind <- CreateSeuratObject(counts = corrected_data,
                                     meta.data = combined_metadata)
    results$combat_ind <- score_and_plot(seu_cb_ind, "combat_ind",
                                          preprocess = TRUE, save_path = save_path)
  }

  return(results)
}
