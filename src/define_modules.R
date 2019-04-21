define_modules <- function(feat_hc_TOM, dissTOM, minModuleSize = 20, deep=1){
  
  tol21 <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")
  
  # Module identification using dynamic tree cut:
  dynamic_mods_TOM <- dynamicTreeCut::cutreeHybrid(dendro = feat_hc_TOM, distM = dissTOM,
                                              deepSplit = deep, pamRespectsDendro = TRUE,
                                              minClusterSize = minModuleSize)
  
  
  palette_generator <- colorRampPalette(tol21)
  
  n_labels <- length(unique(dynamic_mods_TOM$labels))
  label_colors <- palette_generator(n_labels)
  
  # define dendrocolors

  dendrocolors <- data.frame(
    module = labels2colors(
      labels = dynamic_mods_TOM$labels,
      colorSeq = label_colors)
  )
  
  


    
  
  
  return(list(dendrocolors=dendrocolors, dynamic_mods_TOM=dynamic_mods_TOM, feat_hc_TOM=feat_hc_TOM))
}