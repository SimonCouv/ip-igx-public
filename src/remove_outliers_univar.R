remove_outliers_univar <- function(mat, fold=3){
  # set observations to NA if they are > fold*sd from mean 
  # dat: features in cols, samples in rows
  
  # set default for fold
  # (1-pnorm(3))*200
  
  mat[scale(mat, center=TRUE, scale=TRUE)>fold] <- NA
  mat
}