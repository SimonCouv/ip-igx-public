make_heatmap_list <- function(corr_like, feat_hc_TOM, color_overview, modality_df, row_names_fontsize=5, show_module_legend = FALSE, split=NULL, heatmap_body_legend_title="spearman corr"){
  # makes heatmap list in the sense of list(heatmap, heatmap_anno)
  # set blue-red color mapping for use in heatmaps (where values in heatmap body are corr_likes)
  corcol <-define_colors(as.vector(corr_like), breaks = c(-1,0,1), quant_colors = c("blue", "white", "red"))
  
  
  dend_TOM <- as.dendrogram(feat_hc_TOM)
  
  
  # TODO: deal with grey module separately
  module_hm <- HeatmapAnnotation(df=color_overview,
                                 which = "row",
                                 show_legend = show_module_legend, 
                                 width = unit(3,"cm"),
                                 col =  map(as_tibble(color_overview), ~define_colors(.x)),
                                 gp = gpar(col = "white", lwd = 0.5))
  
  # browser()
  
  modality_hm <- HeatmapAnnotation(df=modality_df,
                                   which = "row",
                                   show_legend = TRUE,
                                   col=map(as_tibble(modality_df), ~define_colors(.x)),
                                   gp = gpar(col = "white", lwd = 0.5))
  
  core_hm_modality <- ComplexHeatmap::Heatmap(matrix=corr_like,
                                              name=heatmap_body_legend_title,
                                              col= corcol,
                                              cluster_rows = dend_TOM,
                                              cluster_columns = rev(dend_TOM),
                                              row_names_gp = gpar(fontsize = row_names_fontsize),
                                              column_dend_side = "top",
                                              show_column_names = FALSE,
                                              right_annotation = modality_hm,
                                              split=split)
  # browser()
  hm_list <- module_hm + core_hm_modality
  hm_list
}