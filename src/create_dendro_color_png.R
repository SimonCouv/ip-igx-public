create_dendro_color_png <- function(dendrocolors, feat_hc, fig_path, plt_t, feat_modality=NA){
  

  png(fig_path)
  if (ncol(dendrocolors)<2){
    p <- plotDendroAndColors(feat_hc, 
                        dendrocolors, 
                        # groupLabels = c("module", "modality"),
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        autoColorHeight = FALSE, colorHeight = 0.1,
                        main = plt_t)
  } else {
    p <- plotDendroAndColors(feat_hc, 
                        dendrocolors, 
                        # groupLabels = c("module", "modality"),
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        autoColorHeight = FALSE, colorHeight = 0.1,
                        rowText = feat_modality,
                        textPositions = 2,
                        rowTextAlignment = "center",
                        main = plt_t)
  }
  dev.off()
  
}