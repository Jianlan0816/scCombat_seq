library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(patchwork)

# Assuming your Seurat object is named `seu`
# and has a metadata column named `batch`
library(scuttle)  # Bioconductor package
library(Seurat)

# Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(simu_seurat)
counts_mat <- assay(sce, "counts")
batch_vector <- sce$batch  # adjust if your batch column has a different name

# Downsample to match library size across batches
counts_ds <- downsampleBatches(counts_mat, batch = batch_vector)

# Put downsampled counts back into SCE
assay(sce, "counts") <- counts_ds

# Convert back to Seurat
simu_seurat <- as.Seurat(sce, counts = "counts")

# 1. Compute library size per cell (total UMI counts)
simu_seurat$library_size <- Matrix::colSums(GetAssayData(simu_seurat, slot = "counts", assay = "RNA"))

# 2. Density plot of library sizes by batch
p_density <- ggplot(simu_seurat@meta.data, aes(x = library_size, fill = batch)) +
  geom_density(alpha = 0.4) +
  scale_x_log10() +  # log scale for better visualization
  theme_minimal() +
  labs(title = "Library Size Density by Batch", x = "Library Size (log10)", y = "Density")

# 3. Violin plot for per-batch comparison
p_violin <- ggplot(simu_seurat@meta.data, aes(x = batch, y = library_size, fill = batch)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_y_log10() +
  theme_minimal() +
  labs(title = "Library Size Distribution per Batch", y = "Library Size (log10)")

# 4. UMAP colored by library size
simu_seurat <- NormalizeData(simu_seurat)  # Necessary before running PCA/UMAP
simu_seurat <- FindVariableFeatures(simu_seurat)
simu_seurat <- ScaleData(simu_seurat)
simu_seurat <- RunPCA(simu_seurat)
simu_seurat <- RunUMAP(simu_seurat, dims = 1:10)

p_umap_batch <- DimPlot(simu_seurat, group.by = "batch", reduction = "umap") + 
  ggtitle("UMAP by Batch")

p_umap_libsize <- FeaturePlot(simu_seurat, features = "library_size", reduction = "umap") +
  ggtitle("UMAP by Library Size")

# 5. Combine and visualize
(p_density | p_violin) / (p_umap_batch | p_umap_libsize)

simu_seurat@meta.data <- simu_seurat@meta.data %>%
  dplyr::rename(celltype = cell_type)
