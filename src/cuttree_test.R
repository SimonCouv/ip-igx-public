# aim: test whether cuttreeHybrid returns labels in original order, or in order of feat_hc
# conclusion: module labels maintain the order in the adj/TOM/distance matrices

toy <- tibble(a=runif(50),
              b=5*a+rnorm(50),
              c=5*a+rnorm(50),
              d=runif(50),
              e=runif(50),
              f=5*a+rnorm(50),
              g=5*a+rnorm(50),
              h=runif(50),
              i=runif(50),
              )

adj <- adjacency(datExpr = toy,
                 type="signed hybrid",
                 power = 6,
                 corFnc = "cor")
diss_adj <- dist(adj)

ComplexHeatmap::Heatmap(adj,
                        clustering_method_columns = "average",
                        clustering_method_rows = "average")

TOM <- TOMsimilarity(adj)
dissTOM <-1 - TOM

feat_hc_TOM <- hclust(as.dist(dissTOM), method = "average")

# Module identification using dynamic tree cut:
dynamic_mods_TOM <- dynamicTreeCut::cutreeHybrid(dendro = feat_hc_TOM, distM = dissTOM,
                                                 deepSplit = 4, pamRespectsDendro = TRUE,
                                                 minClusterSize = 2)
dynamic_mods_TOM$labels
