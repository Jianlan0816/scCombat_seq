#' scDesign3-based Realistic Single-cell RNA-seq Simulation with Batch Effects
#'
#' Fits marginal NB models to a reference SCE object, optionally modifies batch
#' and celltype coefficients to inject realistic technical variation, then
#' generates a new count matrix and returns it as a Seurat object.

library(scDesign3)
library(SingleCellExperiment)
library(Seurat)

#' Inject subtle celltype-specific expression shifts into a Seurat object
#'
#' Multiplies gene counts by a per-gene random factor drawn from N(1, strength)
#' independently for each cell type.
#'
#' @param seu_obj   Seurat object with 'celltype' in metadata
#' @param strength  SD of the Normal multiplier distribution (default 0.1)
#' @param seed      Random seed
#' @return Seurat object with modified counts and recomputed PCA/UMAP
inject_subtle_celltype_effects <- function(seu_obj, strength = 0.1, seed = 123) {
  set.seed(seed)
  counts_mat <- GetAssayData(seu_obj, layer = "counts")
  if (!"celltype" %in% colnames(seu_obj@meta.data)) {
    stop("'celltype' column not found in metadata")
  }
  meta_df   <- seu_obj@meta.data
  celltypes <- unique(meta_df$celltype)
  message("Injecting subtle celltype-specific shifts (strength = ", strength, ")")
  for (ct in celltypes) {
    idx    <- which(meta_df$celltype == ct)
    shift  <- rnorm(nrow(counts_mat), mean = 1, sd = strength)
    shift[shift < 0] <- 0
    counts_mat[, idx] <- round(counts_mat[, idx] * shift)
  }
  counts_mat[counts_mat < 0] <- 0
  seu_obj <- SetAssayData(seu_obj, layer = "counts", new.data = counts_mat)
  seu_obj <- NormalizeData(seu_obj, verbose = FALSE) |>
    FindVariableFeatures(verbose = FALSE) |>
    ScaleData(verbose = FALSE) |>
    RunPCA(npcs = 20, verbose = FALSE) |>
    RunUMAP(dims = 1:20, verbose = FALSE)
  seu_obj
}

#' Simulate scRNA-seq data with realistic batch effects using scDesign3
#'
#' Workflow:
#' 1. Construct data and fit marginal NB models from a reference SCE.
#' 2. Modify batch coefficients for a random subset of genes (partial batch effect).
#' 3. Optionally shrink celltype coefficients to make biology more subtle.
#' 4. Generate new counts via copula-based simulation.
#' 5. Return a preprocessed Seurat object.
#'
#' @param example_sce          Reference SingleCellExperiment object with 'celltype' and 'batch' in colData
#' @param assay_use            Assay to use for fitting (default "counts")
#' @param percent_genes_modify Fraction of genes to receive modified batch effects (default 0.3)
#' @param batch_effect_mean    Mean batch shift (log scale) for modified genes (default 2)
#' @param batch_effect_sd      SD of batch shift (default 0.5)
#' @param celltype_shrink      Factor to shrink celltype coefficients (0 = remove, 1 = keep; default 0.2)
#' @param batch_amplify        Factor to amplify batch coefficients (default 2)
#' @param copula               Copula type: "gaussian" or "vine" (vine is more realistic, slower)
#' @param n_cores              Cores for parallel fitting (default 2)
#' @param celltype_batch_balanced Logical; if FALSE, force unbalanced batch×celltype design
#' @param celltype_effect_strength Subtle celltype shift strength passed to inject_subtle_celltype_effects
#' @param seed                 Random seed
#' @return Preprocessed Seurat object with 'batch' and 'celltype' in metadata
simulate_scdesign3 <- function(example_sce,
                               assay_use                 = "counts",
                               percent_genes_modify      = 0.3,
                               batch_effect_mean         = 2,
                               batch_effect_sd           = 0.5,
                               celltype_shrink           = 0.2,
                               batch_amplify             = 2,
                               copula                    = "gaussian",
                               n_cores                   = 2,
                               celltype_batch_balanced   = TRUE,
                               celltype_effect_strength  = 0.8,
                               seed                      = 123) {

  set.seed(seed)

  # Step 1: Construct data and fit marginal models (batch only in mu_formula)
  BATCH_data <- construct_data(
    sce            = example_sce,
    assay_use      = assay_use,
    celltype       = "celltype",
    pseudotime     = NULL,
    spatial        = NULL,
    other_covariates = "batch",
    corr_by        = "1"
  )

  BATCH_marginal <- fit_marginal(
    data       = BATCH_data,
    predictor  = "gene",
    mu_formula = "celltype + batch",
    sigma_formula = "1",
    family_use = "nb",
    n_cores    = n_cores,
    usebam     = FALSE
  )

  # Step 2: Modify batch + celltype coefficients for selected genes
  all_genes    <- seq_along(BATCH_marginal)
  modify_genes <- sample(all_genes, size = round(percent_genes_modify * length(all_genes)))

  BATCH_marginal_mod <- BATCH_marginal
  for (i in modify_genes) {
    coefs    <- BATCH_marginal_mod[[i]]$fit$coefficients
    n_coef   <- length(coefs)
    # Amplify batch shift
    BATCH_marginal_mod[[i]]$fit$coefficients[n_coef] <-
      rnorm(1, mean = batch_effect_mean, sd = batch_effect_sd)
    # Increase dispersion
    if (!is.null(BATCH_marginal_mod[[i]]$fit$theta)) {
      BATCH_marginal_mod[[i]]$fit$theta <-
        BATCH_marginal_mod[[i]]$fit$theta * runif(1, 0.5, 0.8)
    }
  }

  # Shrink celltype coefficients globally (makes biology more subtle)
  for (i in seq_along(BATCH_marginal_mod)) {
    model <- BATCH_marginal_mod[[i]]
    if (is.null(model) || !is.list(model) || is.null(model$fit)) next
    coefs       <- model$fit$coefficients
    idx_ct      <- grep("^celltype", names(coefs))
    idx_batch   <- grep("^batch",    names(coefs))
    coefs[idx_ct]    <- coefs[idx_ct]    * celltype_shrink
    coefs[idx_batch] <- coefs[idx_batch] * batch_amplify
    BATCH_marginal_mod[[i]]$fit$coefficients <- coefs
  }

  # Step 3: Extract parameters and fit copula
  BATCH_para <- extract_para(
    sce          = example_sce,
    marginal_list = BATCH_marginal_mod,
    n_cores      = n_cores,
    family_use   = "nb",
    new_covariate = BATCH_data$newCovariate,
    data         = BATCH_data$dat
  )

  BATCH_copula <- fit_copula(
    sce          = example_sce,
    assay_use    = assay_use,
    marginal_list = BATCH_marginal_mod,
    family_use   = "nb",
    copula       = copula,
    n_cores      = n_cores,
    input_data   = BATCH_data$dat
  )

  # Step 4: Generate new counts
  set.seed(seed)
  new_counts <- simu_new(
    sce            = example_sce,
    mean_mat       = BATCH_para$mean_mat,
    sigma_mat      = BATCH_para$sigma_mat,
    zero_mat       = BATCH_para$zero_mat,
    quantile_mat   = NULL,
    copula_list    = BATCH_copula$copula_list,
    n_cores        = n_cores,
    family_use     = "nb",
    input_data     = BATCH_data$dat,
    new_covariate  = BATCH_data$newCovariate,
    important_feature = BATCH_copula$important_feature,
    filtered_gene  = BATCH_data$filtered_gene
  )

  # Step 5: Build metadata
  meta_df <- as.data.frame(BATCH_data$newCovariate)
  if ("cell_type" %in% colnames(meta_df) && !"celltype" %in% colnames(meta_df)) {
    colnames(meta_df)[colnames(meta_df) == "cell_type"] <- "celltype"
  }

  if (!celltype_batch_balanced) {
    major_ct      <- names(sort(table(meta_df$celltype), decreasing = TRUE))[1]
    idx_major     <- which(meta_df$celltype == major_ct)
    idx_minor     <- which(meta_df$celltype != major_ct)
    meta_df$batch[idx_major] <- "batch1"
    meta_df$batch[idx_minor] <- "batch2"
  }
  meta_df$batch <- as.factor(meta_df$batch)
  stopifnot(all(rownames(meta_df) == colnames(new_counts)))

  # Step 6: Create Seurat object
  seu <- CreateSeuratObject(counts = new_counts, meta.data = meta_df)
  seu <- NormalizeData(seu, verbose = FALSE) |>
    FindVariableFeatures(verbose = FALSE) |>
    ScaleData(verbose = FALSE) |>
    RunPCA(verbose = FALSE) |>
    RunUMAP(dims = 1:20, verbose = FALSE)

  # Inject subtle celltype effects
  if (celltype_effect_strength > 0) {
    seu <- inject_subtle_celltype_effects(seu, strength = celltype_effect_strength, seed = seed)
  }

  # Report batch effect strength
  pcs     <- Embeddings(seu, "pca")
  lm_res  <- summary(lm(pcs[, 1] ~ meta_df$batch))
  message(sprintf("Batch effect (PC1 R²): %.1f%%", 100 * lm_res$adj.r.squared))

  invisible(seu)
}
