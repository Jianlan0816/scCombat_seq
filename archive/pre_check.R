check_data_characteristics <- function(seu_obj, celltype_col = "celltype", batch_col = "batch", covariate_col = "Covariate_cat") {
  message("=== Checking Seurat Object Characteristics ===")
  
  # Set default assay
  DefaultAssay(seu_obj) <- "RNA"
  
  # Get raw counts with Seurat v5-compatible syntax
  counts_mat <- tryCatch({
    GetAssayData(seu_obj, layer = "counts")
  }, error = function(e) {
    NULL
  })
  
  # 1. Raw count check
  has_raw_counts <- !is.null(counts_mat) && is(counts_mat, "dgCMatrix") && all(counts_mat == floor(counts_mat))
  if (has_raw_counts) {
    message("[✓] Raw counts are available.")
  } else {
    message("[!] Raw counts NOT available — likely TPM/CPM. ComBat-Seq is NOT applicable.")
  }
  
  # 2. Metadata check
  meta <- seu_obj@meta.data
  has_celltype <- celltype_col %in% colnames(meta)
  has_batch <- batch_col %in% colnames(meta)
  has_covariate <- covariate_col %in% colnames(meta)
  
  if (has_celltype) message("[✓] Cell type labels found in metadata.") else message("[!] Cell type labels NOT found.")
  if (has_batch) message("[✓] Batch labels found in metadata.") else message("[!] Batch labels NOT found.")
  if (has_covariate) message(paste0("[✓] Categorical covariate '", covariate_col, "' found in metadata."))
  
  # 3. Cell type × batch distribution
  if (has_celltype && has_batch) {
    ct_batch_table <- table(meta[[celltype_col]], meta[[batch_col]])
    print(ct_batch_table)
    if (any(ct_batch_table == 0)) {
      message("[!] Some cell types missing in certain batches → imbalance detected.")
    } else {
      message("[✓] All cell types are represented in each batch.")
    }
  }
  
  # 4. Cell number
  n_cells <- ncol(seu_obj)
  message(paste("[i] Total number of cells:", n_cells))
  if (n_cells > 30000) message("[!] Large dataset detected. Consider downsampling or using PCA-based correction.")
  
  # 5. Library size check
  if (has_raw_counts && has_batch) {
    lib_size <- Matrix::colSums(counts_mat)
    meta$lib_size <- lib_size
    
    message("[i] Checking library size across batches...")
    print(tapply(meta$lib_size, meta[[batch_col]], summary))
    
    if (has_covariate) {
      group_combo <- interaction(meta[[batch_col]], meta[[covariate_col]], drop = TRUE)
      message("[i] Checking library size across batch + ", covariate_col, " combinations...")
      print(tapply(meta$lib_size, group_combo, summary))
    }
  }
  
  # 6. Method recommendation
  if (!has_raw_counts) {
    message("→ Suggested method: ComBat-PCA, Harmony, or Seurat integration (on log-normalized data)")
  } else if (has_celltype && has_batch) {
    message("→ Suggested method: ComBat-ScSeq with celltype as covariate")
  } else if (has_batch && !has_celltype) {
    message("→ Suggested method: ComBat-ScSeq or ComBat-SVASeq (with surrogate variable estimation)")
  } else {
    message("→ Metadata insufficient for any supervised ComBat variant.")
  }
}



library(splatter)
library(SingleCellExperiment)
library(Seurat)
library(ggplot2)
library(patchwork)

expr_matrix <- as.matrix(myFilteredData)
expr_matrix <- Matrix(expr_matrix, sparse = TRUE)
seu <- CreateSeuratObject(counts = expr_matrix,
                          meta.data = mySample)
seu <- SetAssayData(
  object = seu,
  assay = "RNA",
  layer = "data",
  new.data = expr_matrix
)
# 3. Set default assay
DefaultAssay(seu) <- "RNA"

# 4. (Optional) log-transform if TPM is raw
#seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 1e6)

# 5. Continue with standard pipeline
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 20)
seu <- RunUMAP(seu, dims = 1:20)

# 6. Plot by batch or any meta column
DimPlot(seu, group.by = "batch") + theme_classic()

sim = readRDS("~/Desktop/Lab/TB-rutgers/combat/Jinmiao/batch_effect/simulation/simulation_S3_seurat.rds")
if ("Group" %in% colnames(sim@meta.data)) {
  colnames(sim@meta.data)[colnames(sim@meta.data) == "Group"] <- "celltype"
}
colnames(sim@meta.data)
check_data_characteristics(sim)

sim1_result = ComBat_combo(sim, print_raw = TRUE, combat_seq = TRUE, combat_scseq = TRUE,
                           combat_pcseq = TRUE, combat_ind = TRUE,
                           save_path = "/Users/Jianlan/Desktop/Lab/TB-rutgers/combat/Jinmiao/batch_effect/results/combat_combo")
sim1_result_benchmark = benchmark_batch_methods(sim, returnHarmony = TRUE, returnLiger = TRUE, returnSeurat = TRUE,
                                                save_path = "/Users/Jianlan/Desktop/Lab/TB-rutgers/combat/Jinmiao/batch_effect/results/benchmark_methods")
