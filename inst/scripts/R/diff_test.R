# GeneID parsing lives in geneid.R (parse_geneid); source it alongside this file
local({
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!exists("parse_geneid", mode = "function") && !is.null(f)) {
    source(file.path(dirname(f), "geneid.R"))
  }
})

#' Split site ids into gene (DEXSeq groupID) and the remainder (featureID)
#'
#' @param ids character vector of scraps GeneIDs (matrix row names)
#' @return data.frame with gene, gene_name, gene_info, delim; gene_info is the
#'   id with the leading gene symbol and delimiter removed, so
#'   `str_c(gene_name, delim, gene_info)` restores the id
#' @keywords internal
.split_site_ids <- function(ids) {
  if (!exists("parse_geneid", mode = "function")) {
    stop("parse_geneid() not found: source inst/scripts/R/geneid.R first")
  }
  p <- parse_geneid(ids, multi = "first")
  delim <- ifelse(is.na(p$delim), ";", p$delim)
  gene_name <- ifelse(is.na(p$gene), ids, p$gene)
  prefix <- str_c(gene_name, delim)
  gene_info <- ifelse(startsWith(ids, prefix), str_sub(ids, nchar(prefix) + 1), ids)
  data.frame(gene = ids, gene_name = gene_name, gene_info = gene_info,
             delim = delim, stringsAsFactors = FALSE)
}

#' Run DEXSeq differential expression on polyA sites
#'
#' @param mat PA site count matrix
#' @param cell_ids1 cell barcodes for group 1
#' @param cell_ids2 cell barcodes for group 2
#' @param fc_cut min fold change cutoff
#' @param padj_cut max adjusted p value cutoff
#' @param pseudo_n number of pseudobulk profiles to make for each group
#' @param seed random seed for the pseudobulk split
#' @param type filter to only these types of PA sites, set to NULL to use all
#' @param indep_fil independent filtering for DEXSeq
#' @param sep deprecated; the GeneID delimiter is detected per row
#' @param groups1 vector of grouping for cell_ids1
#' @param groups2 vector of grouping for cell_ids2
#' @import dplyr tidyr stringr
#' @return DEXSeq result table
#' @export
PA_DEXSeq <- function(mat,
                      cell_ids1,
                      cell_ids2,
                      fc_cut = 0.25,
                      padj_cut = 0.1,
                      pseudo_n = 4,
                      seed = 1,
                      type = c("3'UTR"),
                      indep_fil = TRUE,
                      sep = NULL,
                      groups1 = NA,
                      groups2 = NA) {
  if (!is.null(type)) {
    mat <- mat[str_detect(rownames(mat), type), ]
  }
  dat2 <- .split_site_ids(rownames(mat))
  keeps2 <- dat2$gene_name %>% table() %>% as.data.frame() %>%
    setNames(c("gene_name", "n")) %>%
    dplyr::filter(n > 1) %>%
    pull(gene_name) %>% 
    unique() %>% 
    as.vector()
  message("number of genes with 2+ PA sites: ", length(keeps2))
  keeps <- dat2 %>% dplyr::filter(gene_name %in% keeps2) %>%
    pull(gene) %>% 
    unique()
  message("number of peaks tested: ", length(keeps))
  mat <- mat[keeps, ]
  
  if (!all(is.na(groups1)) && !all(is.na(groups2))) {
    ngroup <- unique(groups1) %>% length()
    message("using ", ngroup, " custom groups for pseudo-bulks")
    cells_1 <- split(cell_ids1, groups1)
    cells_2 <- split(cell_ids2, groups2)
  } else {
    message("prepare ", pseudo_n, " pseudo-bulks per sample")
    set.seed(seed)
    cells_1 <- split(sample(cell_ids1), 1:length(cell_ids1)%%pseudo_n)
    set.seed(seed)
    cells_2 <- split(sample(cell_ids2), 1:length(cell_ids2)%%pseudo_n)
  }
  
  profile_1 <- matrix(nrow = nrow(mat), ncol = length(cells_1))
  for (i in 1:length(cells_1)) {
    set_temp <- cells_1[[i]]
    mat_temp <- mat[, set_temp]
    if (length(set_temp) > 1) {
      profile_1[, i] <- Matrix::rowSums(mat_temp)
    }
    else {
      warning("not pseudobulk")
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
    }
    else {
      warning("not pseudobulk")
      profile_2[, i] <- mat_temp
    }
  }
  rownames(profile_2) <- rownames(mat)
  colnames(profile_2) <- str_c("Population2_", 1:length(cells_2))
  profile_full <- cbind(profile_1, profile_2)
  message("run DEXSeq prep...")
  st <- data.frame(row.names = c(colnames(profile_1), colnames(profile_2)), 
                   condition = c(rep("target", ncol(profile_1)), rep("comparison", 
                                                                     ncol(profile_2))))
  ft <- .split_site_ids(rownames(profile_full))
  rownames(ft) <- rownames(profile_full)
  # DEXSeq names rows groupID:featureID; map those back to the original ids
  id_map <- setNames(ft$gene, str_c(ft$gene_name, ":", ft$gene_info))
  dxd <- suppressWarnings(DEXSeq::DEXSeqDataSet(profile_full, sampleData = st, 
                                                groupID = ft$gene_name, featureID = ft$gene_info, design = ~sample + 
                                                  exon + condition:exon))
  dxd <- DEXSeq::estimateSizeFactors(dxd, locfunc = genefilter::shorth)
  dxd <- DEXSeq::estimateDispersions(dxd, BPPARAM = BiocParallel::MulticoreParam())
  message("run DEXSeq test...")
  dxd <- DEXSeq::testForDEU(dxd, BPPARAM = BiocParallel::MulticoreParam())
  dxd <- DEXSeq::estimateExonFoldChanges(dxd, BPPARAM = BiocParallel::MulticoreParam())
  dxr1 <- DEXSeq::DEXSeqResults(dxd, independentFiltering = indep_fil)
  message("get output...")
  mat1 <- mat[, cell_ids1, drop = FALSE]
  mat1[mat1 > 0] <- 1
  mat2 <- mat[, cell_ids2, drop = FALSE]
  mat2[mat2 > 0] <- 1
  pct <- data.frame(full_site = rownames(mat), pct_1 = Matrix::rowMeans(mat1), 
                    pct_2 = Matrix::rowMeans(mat2))
  sig <- as.data.frame(dxr1) %>% dplyr::select(gene = groupID, padj, 
                                               log2FC = log2fold_target_comparison) %>% arrange(padj) %>% 
    rownames_to_column("full_site") %>%
    mutate(full_site = unname(id_map[full_site])) %>%
    left_join(pct, by = "full_site") %>%
    dplyr::filter(padj <= padj_cut, abs(log2FC) >= fc_cut) %>%
    dplyr::filter(gene != "NA")
  sig
}