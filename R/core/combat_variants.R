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
#' @param min_batch_coverage Minimum number of batches a celltype must appear in
#'   (with >= min_cells_per_ct_batch cells) to be included in the covariate model.
#'   Default = all batches (full coverage required). Set lower to allow partial coverage.
#' @param min_cells_per_ct_batch Minimum cells a celltype needs in a batch to count
#'   as "represented" in that batch (default = 5).
#' @return Named list; each element has: seurat, kBET_acceptance, batch_LISI, celltype_LISI

library(Seurat)
library(sva)
library(umap)
library(dplyr)
library(ggplot2)
library(cowplot)

# evaluation.R must be sourced before this file (via run_all.R or run_combat_variants.R)

ComBat_combo <- function(seu,
                         subset                = 1.0,
                         print_raw             = FALSE,
                         combat_seq            = FALSE,
                         combat_scseq          = FALSE,
                         combat_pcseq          = FALSE,
                         combat_ind            = FALSE,
                         save_path             = NULL,
                         score_clusters        = TRUE,
                         min_batch_coverage     = NULL,   # default: all batches required
                         min_cells_per_ct_batch = 5,
                         use_adaptive_covariate = TRUE) {

  set.seed(123)
  if (subset < 1.0) {
    cells_to_use <- sample(Cells(seu), size = floor(subset * ncol(seu)))
    seu <- subset(seu, cells = cells_to_use)
  }

  counts   <- GetAssayData(seu, layer = "counts")
  metadata <- seu@meta.data
  metadata$celltype <- as.factor(metadata$celltype)
  results  <- list()

  n_batches <- length(unique(metadata$batch))
  if (is.null(min_batch_coverage)) min_batch_coverage <- n_batches

  # When adaptive covariate is disabled, always include all celltypes
  if (!use_adaptive_covariate) {
    get_usable_celltypes <- function(meta) {
      list(usable = levels(meta$celltype), dropped = character(0))
    }
  }

  # ── Helper: identify which celltypes are safe to use as covariates ────────
  get_usable_celltypes <- function(meta) {
    ct_batch_counts <- table(meta$celltype, meta$batch)
    # A celltype is "covered" in a batch if it has >= min_cells_per_ct_batch cells
    ct_covered_batches <- rowSums(ct_batch_counts >= min_cells_per_ct_batch)
    usable <- names(ct_covered_batches[ct_covered_batches >= min_batch_coverage])
    dropped <- setdiff(levels(meta$celltype), usable)
    if (length(dropped) > 0) {
      message(sprintf(
        "[!] Adaptive covariate: dropping %d confounded celltype(s) from model: %s",
        length(dropped), paste(dropped, collapse = ", ")
      ))
      message(sprintf(
        "    (present in < %d/%d batches with >= %d cells)",
        min_batch_coverage, n_batches, min_cells_per_ct_batch
      ))
    } else {
      message("[OK] All celltypes sufficiently represented across batches.")
    }
    list(usable = usable, dropped = dropped)
  }

  # ── Helper: build adaptive covariate model matrix ─────────────────────────
  # For balanced celltypes: full celltype model (biology is protected).
  # For confounded celltypes: intercept-only rows (batch effect estimated from
  # balanced types is applied, but no celltype term — we can't separate
  # batch from biology for these types anyway).
  make_mod <- function(meta) {
    # If adaptive covariate is disabled, use standard full model immediately
    if (!use_adaptive_covariate) {
      if ("Covariate_cont" %in% colnames(meta)) {
        return(model.matrix(~ celltype + Covariate_cont, data = meta))
      } else {
        return(model.matrix(~ celltype, data = meta))
      }
    }

    ct_info <- get_usable_celltypes(meta)

    if (length(ct_info$dropped) == 0) {
      # All celltypes balanced — use standard full model
      if ("Covariate_cont" %in% colnames(meta)) {
        return(model.matrix(~ celltype + Covariate_cont, data = meta))
      } else {
        return(model.matrix(~ celltype, data = meta))
      }
    }

    # Build model matrix using only the balanced celltypes
    meta_usable <- meta[meta$celltype %in% ct_info$usable, , drop = FALSE]
    meta_usable$celltype <- droplevels(as.factor(meta_usable$celltype))
    if ("Covariate_cont" %in% colnames(meta_usable)) {
      mod_usable <- model.matrix(~ celltype + Covariate_cont, data = meta_usable)
    } else {
      mod_usable <- model.matrix(~ celltype, data = meta_usable)
    }

    # Full model matrix for all cells: intercept-only for confounded types
    mod_full <- matrix(0, nrow = nrow(meta), ncol = ncol(mod_usable))
    rownames(mod_full) <- rownames(meta)
    colnames(mod_full) <- colnames(mod_usable)
    mod_full[, "(Intercept)"] <- 1                          # intercept for everyone
    mod_full[rownames(meta_usable), ] <- mod_usable         # full model for balanced types

    return(mod_full)
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
    if (!use_adaptive_covariate) {
      # Standard: one ComBat_seq call with all celltypes in covariate
      mod       <- make_mod(metadata)
      corrected <- sva::ComBat_seq(as.matrix(counts),
                                   batch     = as.factor(metadata$batch),
                                   covar_mod = mod)
    } else {
      # Adaptive (Option A): two-stage correction
      #   Balanced celltypes → ComBat_seq with celltype covariate
      #   Confounded celltypes → ComBat_seq with batch only (can't protect biology)
      ct_info <- get_usable_celltypes(metadata)

      if (length(ct_info$dropped) == 0) {
        # All balanced — identical to standard
        mod       <- model.matrix(~ celltype, data = metadata)
        corrected <- sva::ComBat_seq(as.matrix(counts),
                                     batch     = as.factor(metadata$batch),
                                     covar_mod = mod)
      } else {
        # Stage 1: correct balanced celltypes with celltype covariate
        idx_bal   <- metadata$celltype %in% ct_info$usable
        meta_bal  <- metadata[idx_bal, , drop = FALSE]
        meta_bal$celltype <- droplevels(as.factor(meta_bal$celltype))
        mod_bal   <- if ("Covariate_cont" %in% colnames(meta_bal))
                       model.matrix(~ celltype + Covariate_cont, data = meta_bal)
                     else
                       model.matrix(~ celltype, data = meta_bal)
        corr_bal  <- sva::ComBat_seq(as.matrix(counts[, idx_bal]),
                                     batch     = as.factor(meta_bal$batch),
                                     covar_mod = mod_bal)

        # Stage 2: correct confounded celltypes with batch only
        idx_conf  <- !idx_bal
        meta_conf <- metadata[idx_conf, , drop = FALSE]
        if (length(unique(meta_conf$batch)) > 1) {
          corr_conf <- sva::ComBat_seq(as.matrix(counts[, idx_conf]),
                                       batch = as.factor(meta_conf$batch))
        } else {
          message("    Confounded celltypes are single-batch; leaving uncorrected.")
          corr_conf <- as.matrix(counts[, idx_conf])
        }

        # Merge back into original cell order
        corrected <- matrix(0, nrow = nrow(counts), ncol = ncol(counts),
                            dimnames = dimnames(counts))
        corrected[, idx_bal]  <- corr_bal
        corrected[, idx_conf] <- corr_conf
      }
    }

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
