#' Gene specificity
#' 
#' Calculating specificity of a given gene to a given sample using Stouffer's method. For a given data matrix, it calculates the z-transform along the rows (gene space) and the columns (sample space) independently which are then converted to a summarized Stouffer z-score which is given by:
#' 
#' Z_stouffer = (Z_rows + Z_columns) / sqrt(2)
#' 
#' @param data_matrix Gene expression data
#' @param method Stouffer's method (default)
#' @return Gene specificity matrix (as a sparse Z-score matrix)
#' 
#' @examples 
#' D <- replicate(expr = rnorm(n = 10, mean = 5, sd = 0.5), n = 10, simplify = TRUE)
#' z <- gene_specificity(D)
gene_specificity <- function(data_matrix, method = "stouffer", ...) {
  
  if (method != "stouffer") stop("Unknown method.")
  if (missing(method)) method <- "stouffer"
  
  z_cols <- scale(data_matrix)
  z_rows <- t(scale(t(data_matrix))) 
  
  Z <- (z_cols + z_rows) / sqrt(2)
  
  Z[Z < 1] <- 0
  Z <- Matrix::Matrix(Z, sparse = TRUE)
  
  return(Z)
} 