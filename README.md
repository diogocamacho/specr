# Gene-Sample Specificity

The `specifyR` package is aimed at defining the specificity of expression of a given gene in an data set. In order to define specificty of a gene to a given tissue, the algorithm applies a dual z-transform to the expression data (one on the rows, and onw on the columns), and using Stouffer's method, we combine both z-scores for any gij to define a corrected z-score value (see ()[] for an exploration of Stouffer's method.) 

## Installation
The `specifyR` package can be easily installed using the `remotes` package. In R, do:

```{r}
library(remotes)
remotes::install_github("diogocamacho/specifyR")
library(specifyR)
```




```{r}
nix <- c(which(samples$sub_tissue == "Cells - Transformed fibroblasts"),
  which(samples$sub_tissue == "Cells - EBV-transformed lymphocytes"),
  which(samples$sub_tissue == "Whole Blood"))


GTEX$expression_data <- GTEX$expression_data[, -nix]
samples <- samples[-nix, ]

Z <- gene_specificity(data_matrix = GTEX$expression_data)

gtex_specificity <- tibble::tibble(entrez_ids = rep(gene_annotations$entrez_id, ncol(Z)),
                     gene_symbol = rep(gene_annotations$gene_symbol, ncol(Z)),
                     gene_description = rep(gene_annotations$gene_description, ncol(Z)),
                     tissue = as.vector(sapply(samples$tissue, rep, nrow(Z))),
                     sub_tissue = as.vector(sapply(samples$sub_tissue, rep, nrow(Z))),
                     gender = as.vector(sapply(samples$gender, rep, nrow(Z))),
                     age = as.vector(sapply(samples$age, rep, nrow(Z))),
                     z_score = as.vector(Z)) %>% 
  dplyr::filter(., z_score > 1)
```