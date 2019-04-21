build_TOM <- function(datExpr, networktype, softpower, corfnc){
  adj <- adjacency(datExpr = datExpr,
                   type = networktype,
                   power=softpower,
                   corFnc = corfnc)
  TOM <- TOMsimilarity(adj)
  TOM
}