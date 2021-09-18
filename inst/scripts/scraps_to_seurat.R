#' Read scraps output to matrix
#' 
#' @param file scraps output table containing 
#' @param n_min  minimum number of observations for a site to be kept
#' @param alt_only if TRUE, only keep genes with alternative polyA sites
#' @param cell_ids if given, use only these cell barcodes, and fill in empty ones
#' @param types filter to only these types of PA sites, set to NULL to use all
#' @param types2 another filter to only these types of PA sites, set to NULL to use all
#' @param pf data.frame with SAF first column and position factor by gene, calculated by `parse_saf_pf`
#' @return count matrix
#' @examples 
#' s_small <- scraps_to_matrix("sample_R2_counts.tsv.gz",
#'                             SeuratObject::pbmc_small,
#'                             alt_only = FALSE)
#' 
scraps_to_matrix <- function(file,
                             n_min = 5,
                             alt_only = FALSE,
                             cell_ids = NULL,
                             types = "3'UTR",
                             types2 = NULL,
                             pf = NULL) {
  # read file
  message("read file")
  dat <- read_tsv(file)
  if (!is.null(cell_ids)) {
    dat <- dat %>% 
      filter(cell %in% cell_ids)
  }
  
  # filter if total counts too low
  dat <- dat %>% group_by(gene) %>% 
    filter(sum(count) >= n_min) %>% 
    ungroup()
  
  if (!is.null(types)) {
    message("filter for ", str_c(types, collapse = ". "))
    dat <- dat %>% filter(str_detect(gene, types))
  }
  
  if (!is.null(types2)) {
    message("filter for ", str_c(types2, collapse = ". "))
    dat <- dat %>% filter(str_detect(gene, types2))
  }
  
  # if TRUE, toss observations where gene only has 1 polyA site
  if (alt_only) {
    message("filter for genes with alternative sites")
    dat2 <- dat %>% separate(gene, 
                             sep = "_", 
                             into = c("gene_name", "gene_info"), 
                             remove = FALSE,
                             extra = "merge") %>% 
      distinct(gene, gene_name, gene_info)
    
    keeps2 <- dat2$gene_name %>% table() %>%
      as.data.frame() %>%
      setNames(c("gene_name", "n")) %>%
      filter(n > 1) %>%
      pull(gene_name) %>% 
      as.vector()
    
    keeps <- dat2 %>% filter(gene_name %in% keeps2) %>% 
      pull(gene)
    
    dat <- dat %>% filter(gene %in% keeps)
  }
  
  # use psi values if pf is used
  if (!is.null(pf)) {
    message("calculate normalized PSI value")
    dat <- dat %>% inner_join(pf, by = c("gene" = "GeneID")) %>% 
      mutate(gene_name = str_remove(gene, "_.+")) %>% 
      mutate(count2 = pf * count) %>% 
      group_by(gene_name, cell) %>% 
      summarize(count = sum(count2)/sum(count)) %>% 
      ungroup() %>% 
      rename(gene = "gene_name")
  }
  
  # transform to matrix
  message("convert to matrix")
  mat <- dat %>% pivot_wider(names_from = cell, values_from = count, values_fill = 0) %>% 
    column_to_rownames("gene")

  # keep consistent cell ids
  if (!is.null(cell_ids)) {
    emptys <- setdiff(cell_ids, colnames(mat))
    
    if (length(emptys) > 0) {
      message("filling empty 0s")
      empty_mat <- matrix(data = 0, nrow = nrow(mat), ncol = length(emptys))
      colnames(empty_mat) <- emptys
      mat <- cbind(mat, empty_mat)
    }
    
    mat <- mat[, cell_ids]
  }
  
  return(mat)
}

#' Add scraps output to existing seurat object as new assay
#' 
#' @param file scraps output table containing 
#' @param object seurat object
#' @param n_min  minimum number of observations for a site to be kept
#' @param alt_only if TRUE, only keep genes with alternative polyA sites
#' @param types filter to only these types of PA sites, set to NULL to use all
#' @param pf data.frame with SAF first column and position factor by gene, calculated by `parse_saf_pf`
#' @return Seurat object with inserted assay
#' @examples 
#' s_small <- scraps_to_seurat("sample_R2_counts.tsv.gz",
#'                             SeuratObject::pbmc_small,
#'                             alt_only = FALSE)
#' 
scraps_to_seurat <- function(file, object, 
                             assay_name = "Asite",
                             n_min = 5,
                             alt_only = TRUE,
                             types = "3'UTR",
                             pf = NULL) {
  cell_ids <- Cells(object)
  
  mat <- scraps_to_matrix(file,
                          object, 
                          assay_name = assay_name,
                          n_min = n_min,
                          alt_only = alt_only,
                          cell_ids = cell_ids,
                          types = types,
                          pf = pf)

  object[[assay_name]] <- CreateAssayObject(counts = mat)
  
  return(object)
}

#' Parse SAF annotation file for position factor (pf) values
#' 
#' @param file SAF file used with Scraps
#' @param alt_only if TRUE, only keep genes with alternative polyA sites
#' @return data.frame with SAF first column and position factor by gene
#' @examples 
#' saf <- parse_saf_pf("ref/polyadb32.hg38.saf.gz")
#' 
parse_saf_pf <- function(file,
                         alt_only = TRUE) {
  bed1 <- parse_saf(file)
  
  # remove symbols that come from different chromosomes
  diffchrom <- bed1 %>% distinct(gene, chrom) %>% 
    group_by(gene) %>% 
    summarise(n = n()) %>% 
    filter(n > 1) %>%
    pull(gene)
  
  # position factor is 0-1
  bed2_plus <- bed1 %>% 
    filter(strand == "+") %>% 
    filter(!(gene %in% diffchrom)) %>% 
    group_by(gene) %>% 
    arrange(.by_group = TRUE, start) %>% 
    mutate(pf = row_number()) %>% 
    mutate(pf = (pf - 1) / max(pf - 1)) %>% 
    ungroup()
  
  bed2_minus <- bed1 %>% 
    filter(strand == "-") %>% 
    filter(!(gene %in% diffchrom)) %>% 
    group_by(gene) %>% 
    arrange(.by_group = TRUE, desc(start)) %>% 
    mutate(pf = row_number()) %>% 
    mutate(pf = (pf - 1) / max(pf - 1)) %>% 
    ungroup()
  
  bed3 <- bind_rows(bed2_plus, bed2_minus)
  
  if (alt_only) {
    bed3 <- bed3 %>% group_by(gene) %>% 
      filter(n() > 1) %>% 
      ungroup()
  }
  
  bed3 %>% select(GeneID, pf)
}

parse_saf <- function(file) {
  # read SAF and split for gene symbol
  saf <- read_tsv(file) %>% 
    mutate(GeneID2 = str_remove(GeneID, Chr)) %>% 
    separate(GeneID2, 
             sep = "_", 
             into = c("gene", 
                      "genbank", 
                      "id", 
                      "chrom", 
                      "pos",
                      "strand",
                      "class"),
             remove = FALSE) %>% 
    mutate(pos = as.numeric(pos)) %>% 
    filter(gene != "NA") %>% 
    filter(!is.na(pos)) %>% 
    mutate(Start = Start - 1)
  
  # convert to bed format
  bed1 <- saf %>% select(chrom = Chr, start = Start, end = End, strand, gene, GeneID)
  bed1
}

parse_gtf <- function(file) {
  gtf <- rtracklayer::import(file)
  tbl_gtf <- as.tbl_intervalcustom(gtf)
}

as.tbl_intervalcustom <- function(x) {
  res <- tibble(
    chrom = as.character(x@seqnames),
    start = x@ranges@start - 1,
    end = x@ranges@start - 1 + x@ranges@width,
    name = x@elementMetadata$gene_name,
    transcript = x@elementMetadata$transcript_id,
    transcript_type =  x@elementMetadata$transcript_type,
    score = rep(".", length(x)),
    strand = as.character(x@strand),
    type = x@elementMetadata$type,
    exon_number = as.numeric(x@elementMetadata$exon_number)
  )
  res <- mutate(res, strand = ifelse(strand == "*", ".", strand))
  res
}

