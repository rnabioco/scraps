#' Run DEXSeq differential expression on polyA sites
#' 
#' @param mat PA site count matrix
#' @param cell_ids1 cell barcodes for group 1
#' @param cell_ids2 cell barcodes for group 2
#' @param fc_cut min fold change cutoff
#' @param padj_cut max adjusted p value cutoff
#' @param pseudo_n number of pseudobulk profiles to make for each group
#' @param types filter to only these types of PA sites, set to NULL to use all
#' @return DEXSeq result table
#' 
PA_DEXSeq <- function(mat,
                      cell_ids1,
                      cell_ids2, 
                      fc_cut = 0.25, 
                      padj_cut = 0.05,
                      pseudo_n = 6,
                      seed = 1,
                      type = c("3'UTR")) {
  # filter by types
  if (!is.null(type)) {
    mat <- mat[str_detect(rownames(mat), type), ]
  }
  
  # only genes with 2 + PA sites
  dat2 <- data.frame(gene = rownames(mat)) %>% 
    separate(gene, 
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
  mat <- mat[keeps, ]
  
  # pseudobulk
  set.seed(seed)
  cells_1 <- split(sample(cell_ids1), 1:length(cell_ids1)%%pseudo_n)
  set.seed(seed)
  cells_2 <- split(sample(cell_ids2), 1:length(cell_ids2)%%pseudo_n)
  
  profile_1 <- matrix(nrow = nrow(mat), ncol = length(cells_1))
  for (i in 1:length(cells_1)) {
    set_temp <- cells_1[[i]]
    mat_temp <- mat[, set_temp]
    if (length(set_temp) > 1) {
      profile_1[, i] <- Matrix::rowSums(mat_temp)
    } else {
      warning("not few cells, not pseudobulk")
      profile_1[, i] <- mat_temp
    }
  }
  rownames(profile_1) <- rownames(mat)
  colnames(profile_1) <- str_c("Population1_", 1:length(cells_1))
  
  profile_2 <- matrix(nrow = nrow(mat), ncol = length(cells_2))
  for (i in 1:length(cells_2)) {
    set_temp <- cells_2[[i]]
    mat_temp <- mat[, set_temp]
    if (length(set_temp) > 1) {
      profile_2[, i] <- Matrix::rowSums(mat_temp)
    } else {
      warning("not few cells, not pseudobulk")
      profile_2[, i] <- mat_temp
    }
  }
  rownames(profile_2) <- rownames(mat)
  colnames(profile_2) <- str_c("Population2_", 1:length(cells_2))

  profile_full <- cbind(profile_1, profile_2)
  
  # prepare DEXSeq
  st <- data.frame(row.names = c(colnames(profile_1), colnames(profile_2)),
                   condition = c(rep("target", ncol(profile_1)),
                                 rep("comparison", ncol(profile_2))))
  
  ft <- data.frame(gene = rownames(profile_full))  %>% 
    separate(gene, 
             sep = "_", 
             into = c("gene_name", "gene_info"), 
             remove = FALSE,
             extra = "merge")
  rownames(ft) <- rownames(profile_full)
  
  dxd <- DEXSeq::DEXSeqDataSet(profile_full, 
                               sampleData = st, groupID = ft$gene_name,
                               featureID = ft$gene_info, design= ~sample+exon+condition:exon)
  # run
  dxd <- DEXSeq::estimateSizeFactors(dxd, locfunc = genefilter::shorth)
  dxd <- DEXSeq::estimateDispersions(dxd)
  dxd <- DEXSeq::testForDEU(dxd)
  dxd <- DEXSeq::estimateExonFoldChanges(dxd)
  dxr1 <- DEXSeq::DEXSeqResults(dxd)

  mat1 <- mat[, cell_ids1, drop = FALSE]
  mat1[mat1 > 0] <- 1
  mat2 <- mat[, cell_ids2, drop = FALSE]
  mat2[mat2 > 0] <- 1
  
  pct <- data.frame(
    full_site = rownames(mat),
    pct_1 = Matrix::rowMeans(mat1),
    pct_2 = Matrix::rowMeans(mat2))
  
  sig <- as.data.frame(dxr1) %>% 
    select(gene = groupID, padj, log2FC = log2fold_target_comparison) %>% 
    arrange(padj) %>% 
    rownames_to_column("full_site") %>% 
    mutate(full_site = str_replace(full_site, ":", "_")) %>% 
    left_join(pct, by = "full_site") %>% 
    filter(padj <= padj_cut, abs(log2FC) >= fc_cut) %>% 
    filter(gene != "NA")
}