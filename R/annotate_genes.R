#' Annotate genes
#' 
#' Given a list of EntrezIDs, this function will map out gene symbols and gene descriptions using the `AnnotationDbi` package. For now it assumes gene IDs from human.
#' 
#' @param entrez_ids A vector with EntrezIDs
#' @return A tibble with EntrezIDs, gene symbol, and gene name
#' 
#' @examples 
#' ann <- annotate_genes("100287102")
annotate_genes <- function(entrez_ids, ...) {
  
  if(missing(entrez_ids)) stop("Need EntrezIDs")
  
  X <- AnnotationDbi::select(x = org.Hs.eg.db, 
                             keys = entrez_ids, 
                             column = c("SYMBOL", "GENENAME"), 
                             keytype = "ENTREZID")
  
  X <- tibble::tibble(entrez_id = X$ENTREZID,
                      gene_symbol = X$SYMBOL,
                      gene_description = X$GENENAME)
  
  return(X)
  
}