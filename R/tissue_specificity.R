tissue_specificity <- function(data_matrix, method = "stouffer") {
  
  z_cols <- scale(data_matrix)
  z_rows <- t(scale(t(data_matrix))) 
  
  Z <- (z_cols + z_rows) / sqrt(2)
  
  return(Z)
}
