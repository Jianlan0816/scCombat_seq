#' Data Quality Checks for scComBat-Seq
#'
#' Pre-flight validation of a Seurat object before applying batch correction.

library(Seurat)
library(Matrix)

#' Check Seurat object characteristics and recommend a correction method
#'
#' Validates raw count availability, metadata completeness, batch×celltype balance,
#' library size distributions, and prints a method recommendation.
#'
#' @param seu_obj      A Seurat object
#' @param celltype_col Column name for cell type labels (default "celltype")
#' @param batch_col    Column name for batch labels (default "batch")
#' @param covariate_col Column name for optional categorical covariate (default "Covariate_cat")
#' @return Invisibly returns a named logical vector of check results
check_data_characteristics <- function(seu_obj,
                                       celltype_col  = "celltype",
                                       batch_col     = "batch",
                                       covariate_col = "Covariate_cat") {
  message("=== scComBat-Seq: Data Characteristics ===")

  DefaultAssay(seu_obj) <- "RNA"

  counts_mat <- tryCatch(
    GetAssayData(seu_obj, layer = "counts"),
    error = function(e) NULL
  )

  has_raw <- !is.null(counts_mat) &&
    (is(counts_mat, "dgCMatrix") || is.matrix(counts_mat)) &&
    all(counts_mat == floor(counts_mat))

  if (has_raw) {
    message("[OK] Raw integer counts detected.")
  } else {
    message("[!]  Raw counts NOT available (likely TPM/CPM). ComBat-Seq is NOT applicable.")
  }

  meta         <- seu_obj@meta.data
  has_celltype <- celltype_col %in% colnames(meta)
  has_batch    <- batch_col    %in% colnames(meta)
  has_cov      <- covariate_col %in% colnames(meta)

  if (has_celltype) message("[OK] Cell type column '", celltype_col, "' found.")
  else              message("[!]  Cell type column '", celltype_col, "' NOT found.")
  if (has_batch)    message("[OK] Batch column '", batch_col, "' found.")
  else              message("[!]  Batch column '", batch_col, "' NOT found.")
  if (has_cov)      message("[OK] Covariate column '", covariate_col, "' found.")

  if (has_celltype && has_batch) {
    ct_batch  <- table(meta[[celltype_col]], meta[[batch_col]])
    n_batches <- ncol(ct_batch)
    message("\n[i] Cell type × batch distribution:")
    print(ct_batch)

    # Cramér's V: 0 = no confounding, 1 = perfect confounding
    chi2      <- suppressWarnings(chisq.test(ct_batch)$statistic)
    n_total   <- sum(ct_batch)
    k         <- min(nrow(ct_batch), n_batches)
    cramers_v <- as.numeric(sqrt(chi2 / (n_total * (k - 1))))
    message(sprintf("[i]  Cramér's V (batch×celltype confounding): %.3f  %s",
                    cramers_v,
                    ifelse(cramers_v > 0.5,
                           "<-- HIGH: celltype covariate likely to hurt",
                    ifelse(cramers_v > 0.2,
                           "<-- MODERATE: use with caution",
                           "<-- LOW: safe to use celltype covariate"))))

    # Per-celltype batch coverage (>= 5 cells to count as represented)
    min_cells_thresh <- 5
    ct_coverage  <- rowSums(ct_batch >= min_cells_thresh)
    confounded   <- names(ct_coverage[ct_coverage < n_batches])
    if (length(confounded) > 0) {
      message("[!]  Confounded celltypes (absent/sparse in >=1 batch): ",
              paste(confounded, collapse = ", "))
      message("     --> Recommended: use combat_scseq with adaptive covariate",
              " (min_batch_coverage param) or combat_seq / combat_pcseq.")
    } else {
      message("[OK] All cell types sufficiently represented in every batch.")
    }
  }

  n_cells <- ncol(seu_obj)
  message(sprintf("[i] Total cells: %d", n_cells))
  if (n_cells > 30000) {
    message("[!]  Large dataset. Consider subset < 1.0 or combat_pcseq for speed.")
  }

  if (has_raw && has_batch) {
    lib_size     <- Matrix::colSums(counts_mat)
    meta$lib_size <- lib_size
    message("\n[i] Library size by batch:")
    print(tapply(meta$lib_size, meta[[batch_col]], summary))
    if (has_cov) {
      grp <- interaction(meta[[batch_col]], meta[[covariate_col]], drop = TRUE)
      message("[i] Library size by batch × covariate:")
      print(tapply(meta$lib_size, grp, summary))
    }
  }

  message("\n[→] Method recommendation:")
  if (!has_raw) {
    message("    Use combat_pcseq, Harmony, or Seurat RPCA (log-normalised input)")
  } else if (has_celltype && has_batch) {
    message("    combat_scseq (preferred) or combat_ind for imbalanced designs")
  } else if (has_batch) {
    message("    combat_seq (no celltype prior) or Harmony")
  } else {
    message("    Insufficient metadata — cannot recommend ComBat variant")
  }

  invisible(c(has_raw = has_raw, has_celltype = has_celltype, has_batch = has_batch))
}
