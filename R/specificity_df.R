#' Specificity data frame
#' 
#' Generates a tibble from the specificity matrix (see \code{\link{gene_specificity}}) and gene annotations (see \code{\link{annotate_genes}}.)
#' 
#' @param specificity_matrix Matrix of gene-sample specificities determined by \code{\link{gene_specificity}}.
#' @param gene_annotations Annotations of genes as defined in \code{\link{annotate_genes}}
#' @return A tibble with EntrezIDs, gene symbol, gene description, sample (as defined in the column names of the data matrix), and z-score.
specificity_df <- function(specificity_matrix, gene_annotations, ...) {
  
  df <- tibble::tibble(entrez_ids = rep(gene_annotations$entrez_id, ncol(specificity_matrix)),
                       gene_symbol = rep(gene_annotations$gene_symbol, ncol(specificity_matrix)),
                       gene_description = rep(gene_annotations$gene_description, ncol(specificity_matrix)),
                       sample = as.vector(sapply(colnames(specificity_matrix), rep, nrow(specificity_matrix))),
                       z_score = as.vector(specificity_matrix)) %>% 
    dplyr::filter(., z_score > 1)
  
  return(df)
}