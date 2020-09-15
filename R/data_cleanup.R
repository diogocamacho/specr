#' Data pre-processing
#' 
#' Cleaning up data for usage in `specr`. This will perform the following:
#' 
#' - for count data (default), it will eliminate genes with zero counts
#' - for count data (default), it will eliminate genes with less than a given number of reads (default = 5) in a defined number of samples (default = 1 sample)
#' - for all data, it will restrict analyses to the highest variant genes
#' 
#' @param data_matrix Gene expression data matrix
#' @param data_type Gene expression data type, as "rma" or "count" (default).
#' @param entrez_ids EntrezIDs for genes in the gene expression matrix, in the same order as the rows of the data matrix
#' @return Processed data matrix and cleaned up EntrezIDs.
#' 
data_cleanup <- function(data_matrix, data_type, entrez_ids, min_counts, min_samples) {
  
  if (missing(data_matrix)) stop("Need gene expression data.")
  
  if(class(data_matrix) != "matrix") data_matrix <- as.matrix(data_matrix)
  
  if (missing(data_type)) {
    message("No data type defined. Assuming count data")
    data_type <- "count" 
  }
  
  if (missing(min_counts)) { 
    message("No minimal number of counts defined. Setting to default")
    min_counts <- 5 
  }
  
  if (missing(min_samples)) {
    message("No minimal number of samples defined. Setting to default")
    min_samples <- 0.05
    sample_count <- ceiling(ncol(data_matrix) * min_samples)
  } else if (min_samples == 0) {
    sample_count <- 1
  }
  
  
  if (data_type == "count") {
    
    message("Identifying transcripts with zero counts...")
    sum_counts <- rowSums(data_matrix)
    zero_counts <- which(sum_counts == 0)
    
    message("Identifying transcript representation across samples...")
    count_mins <- integer(length = nrow(data_matrix)) 
    for (i in seq(1, nrow(data_matrix))) {
      count_mins[i] <- length(which(data_matrix[i, ] > min_counts))
    }
    # count_mins <- apply(data_matrix, 1, function(y) length(which(y > min_counts)))
    keep_genes <- which(count_mins >= sample_count)
    keep_genes <- setdiff(keep_genes, zero_counts)
    
    message("Discarding transcripts with no EntrezIDs...")
    genes_entrez <- which(!is.na(entrez_ids))
    keep_genes <- intersect(keep_genes, genes_entrez)
    
    data_matrix <- data_matrix[keep_genes, ]
    entrez_ids <- entrez_ids[keep_genes]
    
    message("Defining unique EntrezIDs set...")
    uentrez <- unique(entrez_ids)
    len <- sapply(uentrez, function(y) length(which(entrez_ids == y)))
    mult_ent <- which(len != 1)
    
    highest_var <- integer(length = length(mult_ent))
    nix_ids <- vector(mode = "list", length = length(mult_ent))
    for (i in seq(1, length(mult_ent))) {
      a1 <- which(entrez_ids == uentrez[mult_ent[i]])
      a2 <- which.max(diag(var(t(data_matrix[a1, ]))))
      highest_var[i] <- a1[a2]
      nix_ids[[i]] <- a1
    }
    nix_ids <- unlist(nix_ids)
    
    keep_ids <- setdiff(seq(1, length(entrez_ids)), nix_ids)
    keep_ids <- sort(union(keep_ids, highest_var))
    
    clean_data<- data_matrix[keep_ids, ]
    clean_genes <- entrez_ids[keep_ids]
    
    E <- list(expression_data = clean_data, entrez_ids = clean_genes)
    
  } else {
    message("Discarding transcripts with no EntrezIDs...")
    nix <- which(!is.na(entrez_ids))
    if (length(nix) != 0) {
      keep_genes <- setdiff(seq(1, length(entrez_ids)), nix)
    }
    
    data_matrix <- data_matrix[keep_genes, ]
    entrez_ids <- entrez_ids[keep_genes]
    
    message("Defining unique EntrezIDs set...")
    uentrez <- unique(entrez_ids)
    len <- sapply(uentrez, function(y) length(which(entrez_ids == y)))
    mult_ent <- which(len != 1)
    
    highest_var <- integer(length = length(mult_ent))
    nix_ids <- vector(mode = "list", length = length(mult_ent))
    for (i in seq(1, length(mult_ent))) {
      a1 <- which(entrez_ids == uentrez[mult_ent[i]])
      a2 <- which.max(diag(var(t(data_matrix[a1, ]))))
      highest_var[i] <- a1[a2]
      nix_ids[[i]] <- a1
    }
    nix_ids <- unlist(nix_ids)
    
    keep_ids <- setdiff(seq(1, length(entrez_ids)), nix_ids)
    keep_ids <- sort(union(keep_ids, highest_var))
    
    clean_data <- data_matrix[keep_ids, ]
    clean_genes <- entrez_ids[keep_ids]
    
    E <- list(expression_data = clean_data, entrez_ids = clean_genes)
  }
  
  return(E)
}
