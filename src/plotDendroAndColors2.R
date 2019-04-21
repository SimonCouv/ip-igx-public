plotDendroAndColors2 <- function (dendro, colors, groupLabels = NULL, rowText = NULL, 
          rowTextAlignment = c("left", "center", "right"), rowTextIgnore = NULL, 
          textPositions = NULL, setLayout = TRUE, autoColorHeight = TRUE, 
          colorHeight = 0.2, colorHeightBase = 0.2, colorHeightMax = 0.6, 
          rowWidths = NULL, dendroLabels = NULL, addGuide = FALSE, 
          guideAll = FALSE, guideCount = 50, guideHang = 0.2, addTextGuide = FALSE, 
          cex.colorLabels = 0.8, cex.dendroLabels = 0.9, cex.rowText = 0.8, 
          marAll = c(1, 5, 3, 1), saveMar = TRUE, abHeight = NULL, 
          abCol = "red", ...) 
{
  oldMar = par("mar")
  if (!is.null(dim(colors))) {
    nRows = dim(colors)[2]
  }
  else nRows = 1
  if (!is.null(rowText)) 
    nRows = nRows + if (is.null(textPositions)) 
      nRows
  else length(textPositions)
  if (autoColorHeight) 
    colorHeight = colorHeightBase + (colorHeightMax - colorHeightBase) * 
      (1 - exp(-(nRows - 1)/6))
  if (setLayout) 
    layout(matrix(c(1:2), 2, 1), heights = c(1 - colorHeight, 
                                             colorHeight))
  par(mar = c(0, marAll[2], marAll[3], marAll[4]))
  plot(dendro, labels = dendroLabels, cex = cex.dendroLabels, 
       ...)
  if (addGuide) 
    addGuideLines(dendro, count = if (guideAll) 
      length(dendro$height) + 1
      else guideCount, hang = guideHang)
  if (!is.null(abHeight)) 
    abline(h = abHeight, col = abCol)
  par(mar = c(marAll[1], marAll[2], 0, marAll[4]))
  plotColorUnderTree(dendro, colors, groupLabels, cex.rowLabels = cex.colorLabels, 
                     rowText = rowText, rowTextAlignment = rowTextAlignment, 
                     rowTextIgnore = rowTextIgnore, textPositions = textPositions, 
                     cex.rowText = cex.rowText, rowWidths = rowWidths, addTextGuide = addTextGuide)
  if (saveMar) 
    par(mar = oldMar)
  
  # browser()
  # p1 <- recordPlot()
  # p2 <- ggplot(iris, aes(Sepal.Length, Sepal.Width))+geom_point()
  # cowplot::plot_grid(p1,p2)
  
  
}
