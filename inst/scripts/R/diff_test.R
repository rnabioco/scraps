#' Run DEXSeq differential expression on polyA sites
#' 
#' @param mat PA site count matrix
#' @param cell_ids1 cell barcodes for group 1
#' @param cell_ids2 cell barcodes for group 2
#' @param fc_cut min fold change cutoff
#' @param padj_cut max adjusted p value cutoff
#' @param pseudo_n number of pseudobulk profiles to make for each group
#' @param types filter to only these types of PA sites, set to NULL to use all
#' @param indep_fil independent filtering for DEXSeq
#' @return DEXSeq result table
#' 
PA_DEXSeq <- function(mat,
                      cell_ids1,
                      cell_ids2, 
                      fc_cut = 0.5, 
                      padj_cut = 0.1,
                      pseudo_n = 6,
                      seed = 1,
                      type = c("3'UTR"),
                      indep_fil = TRUE) {
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
             extra = "merge")
  
  keeps2 <- dat2$gene_name %>% table() %>%
    as.data.frame() %>%
    setNames(c("gene_name", "n")) %>%
    filter(n > 1) %>%
    pull(gene_name) %>% 
    unique() %>% 
    as.vector()
  message("number of genes with 2+ PA sites: ", length(keeps2))
  
  keeps <- dat2 %>% filter(gene_name %in% keeps2) %>% 
    pull(gene) %>%
    unique()
  message("number of peaks tested: ", length(keeps))
  
  mat <- mat[keeps, ]
  
  # pseudobulk
  message("prepare ", pseudo_n, " pseudo-bulks per sample")
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
  message("run DEXSeq")
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
  dxr1 <- DEXSeq::DEXSeqResults(dxd, independentFiltering = indep_fil)

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

#' prepare motif index from oRNAment files
#' 
#' @param rna cDNA_oRNAment.csv.gz file
#' @param id1 string_to_int_ID_conversion.csv.gz file
#' @param id2 RBP_id_encoding.csv.gz file
#' @param min_score filtering score threshold, default at 0.8
#' @param out saved file name
#' @return motif idx table

# rna <- "/Users/rf/Downloads/Homo_sapiens_cDNA_oRNAment.csv.gz"
# id1 <- "/Users/rf/Downloads/Homo_sapiens_string_to_int_ID_conversion.csv.gz"
# id2 <- "/Users/rf/Downloads/RBP_id_encoding.csv.gz"
# motifs_prep (rna, id1, id2)
motifs_prep <- function(rna,
                      id1,
                      id2,
                      min_score = 0.5,
                      out = "hs_motifs_05") {
  cdna <- read_csv(rna,
                   col_names = c("ensembl_gene_id", "ensembl_transcript_id", "gene_biotype", "transcript_biotype", "transcript_position", "RBPid", "score", "unpaired_probability", "chrom", "region", "start", "end")) %>% 
    filter(region == "3;3")
  conv <- read_csv(id1,
                   col_names = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "ensembl_transcript_id_INT", "ensembl_gene_id,INT"))
  rbps <- read_csv(id2,
                   col_names = c("RBPid", "RBP"))
  
  cdna <- cdna %>% filter(score >= 0.8) 
  
  cdna2 <- cdna %>% 
    inner_join(conv, by = c("ensembl_transcript_id" = "ensembl_transcript_id_INT")) %>% 
    inner_join(rbps %>% mutate(RBPid = as.character(RBPid))) %>% 
    select(chrom, start, end, RBP) %>% 
    mutate(RBP = str_remove(RBP, "\\(.+"))
  
  write_csv(cdna2, paste0(out, "_oRNAament_rbp.csv.gz"))
}

#' find motif differences from 1 pair polyA sites
#' 
#' @param motifs motif file
#' @param site1 
#' @param site2
#' @return motif vector
#' 
# m1 <- read_csv("/Users/rf/scraps_merge/hs_motifs_oRNAament_rbp.csv.gz")
# m2 <- PA_motifs(motifs = m1)
PA_motifs <- function(site1 = "chr14_77512017",
                      site2 = "chr14_77505997", motifs) {
  chr1 <- str_remove(site1, "_.+")
  chr2 <- str_remove(site2, "_.+")
  if (chr1 != chr2) {
    warning("not same chromosome")
    return(NA)
  }
  start <- str_remove(site1, ".+_") %>% as.numeric()
  end <- str_remove(site2, ".+_") %>% as.numeric()
  if (start > end) {
    temp <- end
    end <- start
    start <- temp
  }

  temp <- tibble(chrom = chr1,
                 start = start,
                 end = end)
  valr::bed_intersect(temp, motifs) %>% pull(`RBP.y`) %>% unique() %>% sort() %>% str_c(collapse = ",")
}

#' find motif differences from all differential expression on polyA sites
#' 
#' @param motifs motif file
#' @param res DEXseq results data.frame
#' @return motif table
#' 
# m1 <- read_csv("/Users/rf/scraps_merge/hs_motifs_oRNAament_rbp.csv.gz")
# res <- read_tsv("/Users/rf/scraps_merge/scraps_dex.tsv")
# res2 <- DE_PA_motifs(m1, res)
DE_PA_motifs <- function(motifs,
                         res) {
  res2 <- res %>% separate(full_site, remove = FALSE, sep = "_", into = c("gene", 
                                                                  "genbank", 
                                                                  "id", 
                                                                  "chrom", 
                                                                  "pos",
                                                                  "strand",
                                                                  "class")) %>% 
    mutate(site = str_c(chrom, pos, sep = "_")) %>% 
    select(gene, site)
  
  res3 <- res2 %>%
    group_by(gene) %>%
    filter(n() >= 2) %>% 
    arrange(site) %>% 
    do(as_data_frame(t(combn(.$site, m = 2)))) %>%
    ungroup() %>%
    mutate(diff_motifs = mapply(PA_motifs, V1, V2, MoreArgs = list(motifs = motifs)), comp = str_c(V1, "_vs_", V2))
  
  res3
}
