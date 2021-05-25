#' Add scraps output to existing seurat object as new assay
#' 
#' @param file scraps output table containing 
#' @param object seurat object
#' @param n_min  minimum number of observations for a site to be kept
#' @param alt_only if TRUE, only keep genes with alternative polyA sites
#' 
scraps_to_seurat <- function(file, object, 
                             assay_name = "Asite",
                             n_min = 1,
                             alt_only = TRUE) {
  # read file
  dat <- read_tsv(file)
  
  # filter if total counts too low
  dat <- dat %>% group_by(gene) %>% 
    filter(sum(count) >= n_min) %>% 
    ungroup()
  
  # if TRUE, toss observations where gene only has 1 polyA site
  if (alt_only) {
    dat2 <- dat %>% separate(gene, 
                             sep = "_", 
                             into = c("gene_name", "gene_info"), 
                             remove = FALSE,
                             extra = "merge") %>% 
      distinct(gene_name, gene_info, gene) 
    
    keeps2 <- dat2 %>% 
      group_by(gene_name) %>% 
      summarize(n = n()) %>% 
      filter(n > 1) %>% 
      pull(gene_name)
    
    keeps <- dat2 %>% filter(gene_name %in% keeps2) %>% 
      pull(gene)
    
    dat <- dat %>% filter(gene %in% keeps)
  }
  
  # transform to matrix
  mat <- dat %>% pivot_wider(names_from = cell, values_from = count, values_fill = 0) %>% 
    column_to_rownames("gene")
  
  # keep consistent cell ids
  emptys <- setdiff(Cells(object), colnames(mat))
  
  if (length(emptys) > 0) {
    empty_mat <- matrix(data = 0, nrow = nrow(mat), ncol = length(emptys))
    colnames(empty_mat) <- emptys
    mat <- cbind(mat, empty_mat)
  }
  
  mat <- mat[, Cells(object)]
  object[[assay_name]] <- CreateAssayObject(counts = mat)
  
  return(object)
}

s_small3_A <- scraps_to_seurat("sample_R2_counts.tsv.gz",
                               clustifyr::s_small3,
                               alt_only = FALSE)
