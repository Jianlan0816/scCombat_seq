# make sure seurat has batch, celltype, counts
#Instead of modifying all genes, randomly alter batch coefficients for 20-30% of genes only.
#Use smaller batch effect sizes: mean ~1-2 instead of 5.
#Ensure cell types are balanced across batches.
#Shift mean + dispersion	Realistic batch effects on counts and variability	Tests robust batch correction - theta
#Subtle cell type differences	Small biological differences	Forces batch effects to be visible - fit_marginal
#Vine copula	Stronger gene correlation structure	Matches real scRNA-seq dependencies - fit_copula
#inject small but real cell type effects - inject_subtle_celltype_effects

library(SingleCellExperiment)
library(Seurat)
library(BiocParallel)
library(scDesign3)
# -------------
# INPUT PARAMETERS
# -------------
percent_genes_to_modify <- 0.3   # 30% genes 30% - 70%
batch_effect_mean <- 2           # batch effect size mean 2 - 4
batch_effect_sd <- 0.5           # batch effect size sd
celltype_batch_balanced <- TRUE  # TRUE: balanced, FALSE: unbalanced

# -------------
# 1. Load example SCE
# -------------
example_sce <- readRDS(url("https://figshare.com/ndownloader/files/40581965"))
# Load the SCE object
example_sce <- readRDS("/Users/Jianlan/Desktop/Lab/TB-rutgers/combat/Jinmiao/batch_effect/Data/dataset1/dataset1_sce.rds")
print(example_sce)
# ------------- simulate a new data with batch effect information ----------------
set.seed(123)
simu_res <- scdesign3(sce = example_sce, 
                      assay_use = "TPM", 
                      celltype = "celltype", 
                      pseudotime = NULL, 
                      spatial = NULL, 
                      other_covariates = c("batch"), 
                      mu_formula = "celltype + batch", 
                      sigma_formula = "1", 
                      family_use = "nb", 
                      n_cores = 2, 
                      usebam = FALSE, 
                      corr_formula = "1", 
                      copula = "gaussian", 
                      DT = TRUE, 
                      pseudo_obs = FALSE, 
                      return_model = FALSE)
# Create Seurat object
# -------------
simu_sce <- SingleCellExperiment(list(counts = simu_res$new_count), colData = BATCH_data$newCovariate)
logcounts(simu_sce) <- log1p(counts(simu_sce))

simu_seurat <- CreateSeuratObject(counts = counts(simu_sce), meta.data = as.data.frame(BATCH_data$newCovariate))
check_data_characteristics(simu_seurat)

simu_seurat <- NormalizeData(simu_seurat)
simu_seurat <- FindVariableFeatures(simu_seurat)
simu_seurat <- ScaleData(simu_seurat)
simu_seurat <- RunPCA(simu_seurat)
simu_seurat <- RunUMAP(simu_seurat, dims = 1:20)
# -------------
# 2. Construct data and fit marginal models
# -------------
# After constructing BATCH_data
BATCH_data <- construct_data(
  sce = example_sce,
  assay_use = "TPM",
  celltype = "celltype",
  pseudotime = NULL,
  spatial = NULL,
  other_covariates = c("batch"),
  corr_by = "1"
)

# Now fit marginal model safely
BATCH_marginal <- fit_marginal(
  data = BATCH_data,
  predictor = "gene",
  mu_formula = "celltype + batch",   # Only batch now
  sigma_formula = "1",
  family_use = "nb",
  n_cores = 2,
  usebam = FALSE
)

# -------------
# 3. Modify batch effect coefficients for selected genes
# -------------
set.seed(123)
all_genes <- seq_along(BATCH_marginal)
modify_genes <- sample(all_genes, size = round(percent_genes_to_modify * length(all_genes)))

BATCH_marginal_alter <- BATCH_marginal
for (i in modify_genes) {
  n_coef <- length(BATCH_marginal_alter[[i]]$fit$coefficients)
  
  # ✅ Shift mean (same as before)
  BATCH_marginal_alter[[i]]$fit$coefficients[n_coef] <- 
    rnorm(1, mean = batch_effect_mean, sd = batch_effect_sd)
  
  # ✅ Now also shift dispersion (sigma)
  if (!is.null(BATCH_marginal_alter[[i]]$fit$theta)) {
    # theta = 1/dispersion, lower theta → higher dispersion
    BATCH_marginal_alter[[i]]$fit$theta <- 
      BATCH_marginal_alter[[i]]$fit$theta * runif(1, 0.5, 0.8)  # Increase dispersion 20–50%
  }
}

BATCH_marginal_adjusted <- BATCH_marginal

for (i in seq_along(BATCH_marginal_adjusted)) {
  model <- BATCH_marginal_adjusted[[i]]
  
  if (is.null(model) || !is.list(model) || is.null(model$fit)) next
  
  coefs <- model$fit$coefficients
  names_coefs <- names(coefs)
  
  # Scale down celltype-related coefficients
  idx_celltype <- grep("^celltype", names_coefs)
  coefs[idx_celltype] <- coefs[idx_celltype] * 0.2  # shrink to 20%
  
  # Optionally boost batch effect
  idx_batch <- grep("^batch", names_coefs)
  coefs[idx_batch] <- coefs[idx_batch] * 2  # amplify batch effect
  
  model$fit$coefficients <- coefs
  BATCH_marginal_adjusted[[i]] <- model
}

# -------------
# 4. Simulate New Data
# -------------
BATCH_para_alter <- extract_para(
  sce = example_sce,
  marginal_list = BATCH_marginal_adjusted,
  n_cores = 2,
  family_use = "nb",
  new_covariate = BATCH_data$newCovariate,
  data = BATCH_data$dat
)

BATCH_copula <- fit_copula(
  sce = example_sce,
  assay_use = "counts",
  marginal_list = BATCH_marginal_adjusted,
  family_use = "nb",
  copula = "gaussian",  # ✅ Use "vine" instead of "gaussian"for realistic gene-gene dependencies. 
  n_cores = 2,
  input_data = BATCH_data$dat
)

set.seed(123)
BATCH_newcount_alter <- simu_new(
  sce = example_sce,
  mean_mat = BATCH_para_alter$mean_mat,
  sigma_mat = BATCH_para_alter$sigma_mat,
  zero_mat = BATCH_para_alter$zero_mat,
  quantile_mat = NULL,
  copula_list = BATCH_copula$copula_list,
  n_cores = 2,
  family_use = "nb",
  input_data = BATCH_data$dat,
  new_covariate = BATCH_data$newCovariate,
  important_feature = BATCH_copula$important_feature,
  filtered_gene = BATCH_data$filtered_gene
)

# -------------
# 5. Adjust cell type - batch balance if needed
# -------------
meta_df <- as.data.frame(BATCH_data$newCovariate)
colnames(meta_df)[colnames(meta_df) == "cell_type"] <- "celltype"

if (!celltype_batch_balanced) {
  # Force unbalanced: assign most of one celltype into one batch
  major_celltype <- names(sort(table(meta_df$cell_type), decreasing = TRUE))[1]
  idx_major <- which(meta_df$cell_type == major_celltype)
  idx_minor <- which(meta_df$cell_type != major_celltype)
  
  meta_df$batch[idx_major] <- sample(c("batch1"), length(idx_major), replace = TRUE)
  meta_df$batch[idx_minor] <- sample(c("batch2"), length(idx_minor), replace = TRUE)
}

meta_df$batch <- as.factor(meta_df$batch)
# Ensure rownames match columns
stopifnot(all(rownames(meta_df) == colnames(BATCH_newcount_alter)))

# -------------
# 6. Create Seurat object
# -------------
simu_sce <- SingleCellExperiment(list(counts = BATCH_newcount_alter), colData = meta_df)
logcounts(simu_sce) <- log1p(counts(simu_sce))

simu_seurat <- CreateSeuratObject(counts = counts(simu_sce), meta.data = meta_df)
simu_seurat <- NormalizeData(simu_seurat)
simu_seurat <- FindVariableFeatures(simu_seurat)
simu_seurat <- ScaleData(simu_seurat)
simu_seurat <- RunPCA(simu_seurat)
simu_seurat <- RunUMAP(simu_seurat, dims = 1:20)

inject_subtle_celltype_effects <- function(seu_obj, strength = 0.1, seed = 123) {
  set.seed(seed)
  
  # Get counts
  counts_mat <- GetAssayData(seu_obj, slot = "counts")
  
  # Check celltype metadata exists
  if (!"celltype" %in% colnames(seu_obj@meta.data)) {
    stop("No 'celltype' column found in metadata!")
  }
  
  meta_df <- seu_obj@meta.data
  celltypes <- unique(meta_df$celltype)
  
  message("Injecting subtle celltype-specific shifts...")
  
  for (celltype in celltypes) {
    idx_cells <- which(meta_df$celltype == celltype)
    
    # For each gene, generate a random multiplier ~ N(1, strength)
    random_shift <- rnorm(nrow(counts_mat), mean = 1, sd = strength)
    random_shift[random_shift < 0] <- 0  # avoid negative scaling
    
    # Scale counts for these cells
    counts_mat[, idx_cells] <- round(counts_mat[, idx_cells] * random_shift)
  }
  
  # Ensure counts are non-negative integers
  counts_mat[counts_mat < 0] <- 0
  
  # Update Seurat object
  seu_obj <- SetAssayData(seu_obj, slot = "counts", new.data = counts_mat)
  
  # Optional: re-run normalization
  seu_obj <- NormalizeData(seu_obj)
  seu_obj <- FindVariableFeatures(seu_obj)
  seu_obj <- ScaleData(seu_obj)
  seu_obj <- RunPCA(seu_obj, npcs = 20)
  seu_obj <- RunUMAP(seu_obj, dims = 1:20)
  
  return(seu_obj)
}
simu_seurat <- inject_subtle_celltype_effects(simu_seurat, strength = 0.8)

# -------------
# 7. Output
# -------------
print(simu_seurat)
DimPlot(simu_seurat, group.by = "batch")
DimPlot(simu_seurat, group.by = "celltype")

# Optional: save the simulated dataset
# saveRDS(simu_seurat, file = "simulated_batch_dataset.rds")
pcs <- Embeddings(simu_seurat, "pca")
summary(lm(pcs[,1] ~ meta_df$batch))  # % variance in PC1 due to batch

res_combat <- ComBat_combo(simu_seurat, print_raw = FALSE, combat_seq = FALSE, combat_scseq = TRUE,
                           combat_pcseq = TRUE, combat_ind = TRUE,
                           save_path = "/Users/Jianlan/Desktop/Lab/TB-rutgers/combat/Jinmiao/batch_effect/results/combat_combo")

res_benchmark <- benchmark_batch_methods(simu_seurat,
                                         save_path = "/Users/Jianlan/Desktop/Lab/TB-rutgers/combat/Jinmiao/batch_effect/results/benchmark_methods")
