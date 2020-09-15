#' Gene specificity query on GTEx data
#' 
#' 
#' 
gtex_gene_query <- function(gene, query = "symbol", min_specificity, min_frequency, ...) {
  
  if (missing(min_specificity)) min_specificity <- 5
  if (missing(min_frequency)) min_frequency <- 0.1
  
  # this filters the gene of interest
  a1 <- gtex_specificity %>% dplyr::filter(., gene_symbol == gene, z_score >= min_specificity) 
  
  # this counts how many samples per sub-tissue exist
  a2 <- samples %>% dplyr::group_by(., sub_tissue, gender, age) %>% dplyr::count()
  
  # combine both
  a3 <- a1 %>% dplyr::mutate(., sample_count = a2$n[match(paste(a1$sub_tissue, a1$gender, a1$age), paste(a2$sub_tissue, a2$gender, a2$age))])
  
  # reduce to single mappings
  a4 <- a3 %>% 
    dplyr::group_by(., sub_tissue, gender, age, sample_count) %>% 
    dplyr::count() #%>% 
    # dplyr::mutate(., proportion = n / sample_count) %>% 
    # dplyr::filter(., proportion >= min_frequency)
    
  # plot distribution of specificity
  a5 <- a4 %>% group_by(., sub_tissue, gender) %>% dplyr::add_tally(n) %>% dplyr::distinct(., sub_tissue, n) 
  
  tmp1 <- tibble::tibble(tissue = unique(gtex_specificity$sub_tissue),
                         count = 0)
  
  tmp2 <- tibble::tibble(tissue = unique(gtex_specificity$sub_tissue),
                         count = 0)
  
  tmp1$count[which(tmp1$tissue %in% a5$sub_tissue[a5$gender == "male"])] <- a5$n[a5$gender == "male"]
  tmp1$normalized <- (tmp1$count - min(tmp1$count)) / (max(tmp1$count) - min(tmp1$count))
  tmp2$count[which(tmp2$tissue %in% a5$sub_tissue[a5$gender == "female"])] <- a5$n[a5$gender == "female"]
  tmp2$normalized <- (tmp2$count - min(tmp2$count)) / (max(tmp2$count) - min(tmp2$count))
  tmp1$gender <- "male"
  tmp2$gender <- "female"
  
  plot_df <- bind_rows(tmp1, tmp2)
  
  
  p <- plot_df %>%
    ggplot() +
    geom_bar(stat = "identity", aes(x = tissue, y = normalized, fill = tissue), color = "black") + 
    scale_fill_viridis_d() +
    facet_grid(gender ~ ., scales = "free_y") +
    labs(title = gene, x = NULL, y = "Relative specificity") +
    theme_bw() + 
    theme(axis.text.x = element_text(size = 12, color = "black", angle = 90, hjust = 1, vjust = 0.5),
          axis.text.y = element_text(size = 12, color = "black"),
          plot.title = element_text(size = 24, face = "bold", color = "black"),
          axis.ticks.length = unit(0.2, "cm"),
          axis.title.y = element_text(size = 12, color = "black"),
          panel.grid = element_blank(),
          legend.position = "none", 
          strip.text.y = element_text(size = 12), 
          strip.background = element_blank())
  
  print(p)
  
  return(a4)
  
}