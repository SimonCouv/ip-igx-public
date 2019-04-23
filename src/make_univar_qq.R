make_univar_qq <- function(v_pval, fig_path){
  png(fig_path, res=600,
      width = 6, height = 6, units = "in")
  qqman::qq(v_pval)
  dev.off()
  
}