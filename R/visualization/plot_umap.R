#' Publication-quality UMAP Panel Generator for scComBat-Seq
#'
#' Produces paired celltype/batch UMAP panels (and optional covariate panel)
#' for a single correction method.

library(Seurat)
library(ggplot2)
library(cowplot)
library(RColorBrewer)

#' Generate a multi-panel UMAP plot for one method
#'
#' @param seu         Seurat object with 'celltype' and 'batch' in metadata
#' @param method_name Character label shown as panel title
#' @param reduction   Name of the reduction to use; "harmony" uses harmony embeddings,
#'                    otherwise "umap" is used
#' @param alpha_val   Point transparency (default 0.6)
#' @param pt_size     Point size (default 0.3)
#' @return A ggplot/cowplot object (2 or 3 columns: celltype, batch, optional covariate)
generate_umap_pair <- function(seu, method_name,
                               reduction = NULL,
                               alpha_val = 0.6,
                               pt_size   = 0.3) {
  if (is.null(reduction)) {
    reduction <- if (method_name == "harmony") "harmony" else "umap"
  }

  emb_df           <- as.data.frame(Embeddings(seu, reduction = reduction))
  colnames(emb_df)[1:2] <- c("umap_1", "umap_2")
  emb_df$celltype  <- seu$celltype
  emb_df$batch     <- seu$batch
  has_cat          <- "Covariate_cat" %in% colnames(seu@meta.data)
  if (has_cat) emb_df$cat <- seu$Covariate_cat

  celltypes     <- unique(emb_df$celltype)
  batches       <- unique(emb_df$batch)
  n_ct          <- max(3, length(celltypes))
  n_batch       <- max(3, length(batches))

  celltype_cols <- setNames(RColorBrewer::brewer.pal(min(n_ct, 8),   "Set2"), celltypes)
  batch_cols    <- setNames(RColorBrewer::brewer.pal(min(n_batch, 8), "Dark2"), batches)

  base_theme <- theme_classic() + theme(
    legend.position = "none",
    axis.title       = element_blank(),
    axis.ticks       = element_blank(),
    axis.text        = element_blank(),
    plot.margin      = margin(1, 1, 1, 1)
  )

  p_ct <- ggplot(emb_df, aes(umap_1, umap_2, color = celltype)) +
    geom_point(alpha = alpha_val, size = pt_size) +
    scale_color_manual(values = celltype_cols) + base_theme

  p_bt <- ggplot(emb_df, aes(umap_1, umap_2, color = batch)) +
    geom_point(alpha = alpha_val, size = pt_size) +
    scale_color_manual(values = batch_cols) + base_theme

  panels <- list(p_ct, p_bt)

  if (has_cat) {
    cats     <- unique(emb_df$cat)
    n_cat    <- max(3, length(cats))
    cat_cols <- setNames(RColorBrewer::brewer.pal(min(n_cat, 8), "Set1"), cats)
    p_cat    <- ggplot(emb_df, aes(umap_1, umap_2, color = cat)) +
      geom_point(alpha = alpha_val, size = pt_size) +
      scale_color_manual(values = cat_cols) + base_theme
    panels <- c(panels, list(p_cat))
  }

  title_row  <- ggdraw() + draw_label(method_name, fontface = "bold", hjust = 0.5, size = 12)
  pair_panel <- plot_grid(plotlist = panels, ncol = length(panels))
  plot_grid(title_row, pair_panel, ncol = 1, rel_heights = c(0.15, 1))
}

#' Build a shared-legend figure from a results list
#'
#' @param results_list Named list; each element must contain a $seurat Seurat object.
#'                     Elements with NULL $seurat (e.g. combat_pcseq) are silently skipped.
#' @param alpha_val    Point transparency passed to generate_umap_pair
#' @param pt_size      Point size passed to generate_umap_pair
#' @param save_path    If non-NULL, saves the main panels + legend PNGs here
#' @return Invisible NULL (side-effect: plots/saves figures)
assemble_figure <- function(results_list,
                            alpha_val  = 0.08,
                            pt_size    = 0.2,
                            save_path  = NULL) {
  valid <- Filter(function(r) !is.null(r$seurat), results_list)
  if (length(valid) == 0) {
    message("No Seurat objects found in results_list.")
    return(invisible(NULL))
  }

  paired_panels <- lapply(names(valid), function(nm) {
    generate_umap_pair(valid[[nm]]$seurat, nm,
                       alpha_val = alpha_val, pt_size = pt_size)
  })

  all_pairs <- plot_grid(plotlist = paired_panels, nrow = 1)

  # Build shared legends from first valid object
  seu_ex   <- valid[[1]]$seurat
  umap_df  <- as.data.frame(Embeddings(seu_ex, "umap"))
  colnames(umap_df)[1:2] <- c("umap_1", "umap_2")
  umap_df$celltype <- seu_ex$celltype
  umap_df$batch    <- seu_ex$batch

  ct_levels    <- unique(umap_df$celltype)
  batch_levels <- unique(umap_df$batch)
  ct_cols      <- setNames(brewer.pal(min(max(3, length(ct_levels)), 8),    "Set2"),  ct_levels)
  batch_cols   <- setNames(brewer.pal(min(max(3, length(batch_levels)), 8), "Dark2"), batch_levels)

  leg_ct    <- get_legend(ggplot(umap_df, aes(umap_1, umap_2, color = celltype)) +
                            geom_point() + scale_color_manual(values = ct_cols) +
                            theme_void() + theme(legend.position = "right"))
  leg_batch <- get_legend(ggplot(umap_df, aes(umap_1, umap_2, color = batch)) +
                            geom_point() + scale_color_manual(values = batch_cols) +
                            theme_void() + theme(legend.position = "right"))

  leg_panels <- list(leg_ct, leg_batch)
  has_cat    <- "Covariate_cat" %in% colnames(seu_ex@meta.data)
  if (has_cat) {
    umap_df$cat <- seu_ex$Covariate_cat
    cat_levels  <- unique(umap_df$cat)
    cat_cols    <- setNames(brewer.pal(min(max(3, length(cat_levels)), 8), "Set1"), cat_levels)
    leg_cat     <- get_legend(ggplot(umap_df, aes(umap_1, umap_2, color = cat)) +
                                geom_point() + scale_color_manual(values = cat_cols) +
                                theme_void() + theme(legend.position = "right"))
    leg_panels  <- c(leg_panels, list(leg_cat))
  }

  legend_combined <- plot_grid(plotlist = leg_panels, ncol = 1)
  final_plot      <- plot_grid(all_pairs, legend_combined, rel_widths = c(0.9, 0.1))

  if (!is.null(save_path)) {
    dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
    n_methods <- length(valid)
    png(file.path(save_path, "figure_all_methods.png"),
        width = 1500 * n_methods, height = 500, res = 300)
    print(all_pairs)
    dev.off()
    png(file.path(save_path, "legend.png"), width = 500, height = 800, res = 300)
    print(legend_combined)
    dev.off()
    message("Saved figure to ", save_path)
  }

  invisible(final_plot)
}
