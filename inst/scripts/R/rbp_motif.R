# GeneID parsing lives in geneid.R (parse_geneid); source it alongside this file
local({
  f <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!exists("parse_geneid", mode = "function") && !is.null(f)) {
    source(file.path(dirname(f), "geneid.R"))
  }
})

#' prepare motif index from oRNAment files
#' 
#' @param rna cDNA_oRNAment.csv.gz file
#' @param id1 string_to_int_ID_conversion.csv.gz file
#' @param id2 RBP_id_encoding.csv.gz file
#' @param min_score filtering score threshold, default at 0.8
#' @param out saved file name
#' @return motif idx table
#' @import dplyr tidyr stringr readr
#' @export
motifs_prep <- function(rna = "http://rnabiology.ircm.qc.ca/BIF/oRNAment/static/Homo_sapiens_cDNA_oRNAment.csv.gz",
                        id1 = "http://rnabiology.ircm.qc.ca/BIF/oRNAment/static/Homo_sapiens_string_to_int_ID_conversion.csv.gz",
                        id2 = "http://rnabiology.ircm.qc.ca/BIF/oRNAment/static/RBP_id_encoding.csv.gz",
                        min_score = 0.5,
                        out = "hs_motifs_05") {
  cdna <- read_csv(rna,
                   col_names = c("ensembl_gene_id", "ensembl_transcript_id", "gene_biotype", "transcript_biotype", "transcript_position", "RBPid", "score", "unpaired_probability", "chrom", "region", "start", "end")) %>% 
    filter(region == "3;3")
  conv <- read_csv(id1,
                   col_names = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "ensembl_transcript_id_INT", "ensembl_gene_id,INT"))
  rbps <- read_csv(id2,
                   col_names = c("RBPid", "RBP"))
  
  cdna <- cdna %>% filter(score >= min_score)
  
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
#' @param site1 site key `chrom_pos` for site 1 (chrom may itself contain `_`)
#' @param site2 site key `chrom_pos` for site 2
#' @return motif vector
#' @import dplyr tidyr stringr
#' @export
PA_motifs <- function(site1 = "chr14_77512017",
                      site2 = "chr14_77505997", motifs) {
  chr1 <- str_remove(site1, "_[^_]+$")
  chr2 <- str_remove(site2, "_[^_]+$")
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
#' @import dplyr tidyr stringr
#' @export
DE_PA_motifs <- function(motifs,
                         res) {
  if (!exists("parse_geneid", mode = "function")) {
    stop("parse_geneid() not found: source inst/scripts/R/geneid.R first")
  }
  p <- parse_geneid(res$full_site, multi = "first")
  res2 <- res %>%
    mutate(gene = p$gene, site = str_c(p$chrom, p$pos, sep = "_")) %>%
    select(gene, site)

  res3 <- res2 %>%
    group_by(gene) %>%
    filter(n() >= 2) %>%
    arrange(site) %>%
    do(as_tibble(t(combn(.$site, m = 2)), .name_repair = "minimal")) %>%
    ungroup() %>%
    mutate(diff_motifs = mapply(PA_motifs, V1, V2, MoreArgs = list(motifs = motifs)), comp = str_c(V1, "_vs_", V2))
  
  res3
}