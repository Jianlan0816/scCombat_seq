library(Seurat)
library(ggplot2)
library(patchwork)
library(scales)
library(RColorBrewer)

# Custom function: create transparent UMAPs for a Seurat object
generate_umap_pair <- function(seu, method_name, alpha_val = 0.6, pt_size = 0.3) {
  # Determine which reduction embedding to use
  reduction_to_use <- if (method_name == "harmony") {
    "harmony"
  } else if (method_name == "liger") {
    "umap"  # The reduction name in Seurat converted from LIGER is usually "inmf"
  } else {
    "umap"
  }
  emb_df <- as.data.frame(Embeddings(seu, reduction = reduction_to_use))
  colnames(emb_df)[1:2] <- c("umap_1", "umap_2")  # ensure standard naming
  
  # Add metadata
  emb_df$celltype <- seu$celltype
  emb_df$batch <- seu$batch
  if ("Covariate_cat" %in% colnames(seu@meta.data)){
    emb_df$cat <- seu$Covariate_cat}
  
  # Colors
  celltypes <- unique(emb_df$celltype)
  batches <- unique(emb_df$batch)
  if ("cat" %in% colnames(emb_df)){
  cats <- unique(emb_df$cat)
  cat_cols <- setNames(brewer.pal(max(6, length(cats)), "Set1"), cats)
  p_cat <- ggplot(emb_df, aes(umap_1, umap_2, color = cat)) +
    geom_point(alpha = alpha_val, size = pt_size) +
    scale_color_manual(values = cat_cols) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = margin(1,1,1,1)
    )
  }
  
  celltype_cols <- setNames(brewer.pal(max(3, length(celltypes)), "Set2"), celltypes)
  batch_cols <- setNames(brewer.pal(max(3, length(batches)), "Dark2"), batches)
  
  # Celltype UMAP
  p_celltype <- ggplot(emb_df, aes(umap_1, umap_2, color = celltype)) +
    geom_point(alpha = alpha_val, size = pt_size) +
    scale_color_manual(values = celltype_cols) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = margin(1,1,1,1)
    )
  p_celltype
  # Batch UMAP
  p_batch <- ggplot(emb_df, aes(umap_1, umap_2, color = batch)) +
    geom_point(alpha = alpha_val, size = pt_size) +
    scale_color_manual(values = batch_cols) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_blank(),
      plot.margin = margin(1,1,1,1)
    )
  p_batch
  # Pair with title
  title <- ggdraw() + draw_label(method_name, fontface = 'bold', hjust = 0.5, size = 12)
  if ("cat" %in% colnames(emb_df)){
  pair_plot <- cowplot::plot_grid(p_celltype, p_batch, p_cat, ncol = 3, rel_widths = c(0.5,0.5,0.5))}
  else{pair_plot <- cowplot::plot_grid(p_celltype, p_batch, ncol = 2, rel_widths = c(0.5,0.5))}
  full_panel <- cowplot::plot_grid(title, pair_plot, ncol = 1, rel_heights = c(0.2, 0.9))
  return(full_panel)
}

results_list <- c(sim1_result, sim1_result_benchmark)
for (i in seq_along(results_list)) {
  method <- method_names[i]
  print(method)
  print(results_list[[i]][["kBet_acceptance"]])
  print(results_list[[i]][["batch_LISI"]])
  print(results_list[[i]][["celltype_LISI"]])
}
# Generate all paired panels
paired_panels <- list()
method_names <- names(results_list)

for (i in seq_along(results_list)) {
  method <- method_names[i]
  print(method)
  seu <- results_list[[i]]$seurat
  print(colnames(seu@meta.data))
  paired_panels[[i]] <- generate_umap_pair(seu, method, alpha_val = 0.08, pt_size = 0.2)
}

# Arrange all in 2 rows
all_pairs <- cowplot::plot_grid(plotlist = paired_panels, nrow = 1)

# Shared legends (use one example Seurat object)
seu_example <- results_list[[1]]$seurat
umap_df <- as.data.frame(Embeddings(seu_example, "umap"))
umap_df$celltype <- seu_example$celltype
umap_df$batch <- seu_example$batch

if ("Covariate_cat" %in% colnames(seu_example@meta.data)){
umap_df$cat <- seu_example$Covariate_cat
cat_levels <- unique(umap_df$cat)
cat_cols <- setNames(brewer.pal(max(6, length(cat_levels)), "Set1"), cat_levels)
legend_cat <- ggplot(umap_df, aes(umap_1, umap_2, color = cat)) +
  geom_point(size = 3) +
  scale_color_manual(values = cat_cols) +
  theme_void() +
  theme(legend.position = "right")
}

# Define legend plots
celltype_levels <- unique(umap_df$celltype)
batch_levels <- unique(umap_df$batch)

celltype_cols <- setNames(brewer.pal(max(3, length(celltype_levels)), "Set2"), celltype_levels)
batch_cols <- setNames(brewer.pal(max(3, length(batch_levels)), "Dark2"), batch_levels)

legend_celltype <- ggplot(umap_df, aes(umap_1, umap_2, color = celltype)) +
  geom_point(size = 3) +
  scale_color_manual(values = celltype_cols) +
  theme_void() +
  theme(legend.position = "right")

legend_batch <- ggplot(umap_df, aes(umap_1, umap_2, color = batch)) +
  geom_point(size = 3) +
  scale_color_manual(values = batch_cols) +
  theme_void() +
  theme(legend.position = "right")

if (exists("legend_cat")){
  legend_combined <- cowplot::plot_grid(get_legend(legend_celltype),
                                        get_legend(legend_batch),
                                        get_legend(legend_cat),
                                        ncol = 1)
  }else{legend_combined <- cowplot::plot_grid(get_legend(legend_celltype),
                             get_legend(legend_batch),
                             ncol = 1)}

# Final layout with legends
final_plot <- cowplot::plot_grid(all_pairs, legend_combined, rel_widths = c(0.9, 0.1))
# 🔹 Save to PNG
png("/Users/Jianlan/Desktop/Lab/TB-rutgers/combat/Jinmiao/batch_effect/simulation/sim3.png", 
    width = 6888, height = 500, res = 300)
print(all_pairs)
dev.off()

png("/Users/Jianlan/Desktop/Lab/TB-rutgers/combat/Jinmiao/batch_effect/simulation/legend3.png", 
    width = 500, height = 800, res = 300)
print(legend_combined)
dev.off()
