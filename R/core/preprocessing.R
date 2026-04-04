#' Preprocessing Utilities for scComBat-Seq
#'
#' Functions for filtering, normalizing, and converting raw count matrices.

#' Convert a large sparse matrix to dense in chunks to avoid memory spikes
#'
#' @param sparse_matrix A sparse matrix (dgCMatrix or similar)
#' @param by Number of columns per chunk (default 100,000)
#' @return Dense matrix
preprocess_big <- function(sparse_matrix, by = 100000) {
  n_col <- ncol(sparse_matrix)
  n     <- floor(n_col / by)

  if (n < 1) {
    return(as.matrix(sparse_matrix))
  }

  res <- NULL
  for (i in 1:n) {
    mat <- as.matrix(sparse_matrix[, ((i - 1) * by + 1):(i * by)])
    res <- if (is.null(res)) mat else cbind(res, mat)
    rm(mat)
  }
  if (n_col > n * by) {
    res <- cbind(res, as.matrix(sparse_matrix[, (n * by + 1):n_col]))
  }
  rm(sparse_matrix)
  return(res)
}

#' Filter cells and genes from a raw count matrix, then median-normalise and log-transform
#'
#' @param myData Matrix of raw counts (genes x cells)
#' @param is_filter_cells Logical; filter cells with fewer than min_genes detected genes
#' @param min_genes Minimum number of genes per cell (used when is_filter_cells = TRUE)
#' @param is_filter_genes Logical; filter genes detected in fewer than min_cells cells
#' @param min_cells Minimum number of cells per gene (used when is_filter_genes = TRUE)
#' @return List with two elements:
#'   \item{counts}{Filtered raw count matrix (for ComBat-Seq input)}
#'   \item{log_norm}{Median-normalised, log-transformed matrix (for dimensionality reduction)}
filter_data_mtx <- function(myData,
                             is_filter_cells = FALSE, min_genes = 300,
                             is_filter_genes = FALSE, min_cells = 10) {
  if (is_filter_cells) {
    num_genes  <- colSums(myData > 0)
    cells_use  <- names(num_genes[num_genes > min_genes])
    myData     <- myData[, cells_use]
  }
  if (is_filter_genes) {
    num_cells  <- rowSums(myData > 0)
    genes_use  <- names(num_cells[num_cells > min_cells])
    myData     <- myData[genes_use, ]
  }

  message("Median-normalising counts and log-transforming")
  col_sums   <- colSums(myData)
  med_trans  <- median(col_sums)
  norm_counts <- med_trans * scale(myData, center = FALSE, scale = col_sums)
  log_norm   <- as.data.frame(log(norm_counts + 1))

  list(counts = myData, log_norm = log_norm)
}

#' Create and preprocess a Seurat v5 object from a raw count matrix
#'
#' @param myData Raw count matrix (genes x cells)
#' @param mySample Sample/cell metadata data.frame with at minimum 'batch' and 'celltype' columns
#' @param min_cells Minimum cells per gene for Seurat filtering
#' @param min_genes Minimum genes per cell for Seurat filtering
#' @param regress_umi Logical; whether to regress out nCount_RNA during ScaleData
#' @return Preprocessed Seurat object (normalised, scaled, variable genes identified)
filter_data_seurat <- function(myData, mySample,
                               min_cells = 10, min_genes = 300,
                               regress_umi = FALSE) {
  library(Seurat)

  seu <- CreateSeuratObject(
    counts    = myData,
    project   = "scComBat_seq",
    min.cells = min_cells,
    min.features = min_genes
  )

  cells_use <- Cells(seu)
  mySample  <- mySample[cells_use, , drop = FALSE]

  for (col in intersect(c("batch", "batchlb", "celltype"), colnames(mySample))) {
    seu[[col]] <- mySample[[col]]
  }

  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e4)

  if (regress_umi) {
    seu <- ScaleData(seu, vars.to.regress = "nCount_RNA")
  } else {
    seu <- ScaleData(seu)
  }
  seu <- FindVariableFeatures(seu)

  return(seu)
}

#' Track and export runtime to a TSV file
#'
#' @param t1 Start time (from Sys.time())
#' @param t2 End time
#' @param label Label for this timing record
#' @param output_path Output file path (without extension)
runtime_export <- function(t1, t2, label, output_path) {
  secs <- as.numeric(difftime(t2, t1, units = "secs"))
  mins <- as.numeric(difftime(t2, t1, units = "mins"))
  message(sprintf("[%s] %.1f secs / %.2f mins", label, secs, mins))
  df <- data.frame(label = label, secs = secs, mins = mins)
  write.table(df, file = paste0(output_path, "_runtime.txt"),
              row.names = FALSE, col.names = TRUE, sep = "\t")
}
