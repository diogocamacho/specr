rnaseq_filtering <- function(data_matrix, min_counts, min_samples, ...) {
  
  if(missing(min_counts)) min_counts <- 5
  if(missing(min_samples)) min_samples <- 0.75
  
  scount <- ceiling(ncol(data_matrix) * min_samples)
  
  xx <- apply(data_matrix, 1, function(y) length(which(y > min_counts)))
  
  gids <- which(xx >= scount)
  
  return(gids)
}
