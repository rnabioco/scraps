#' Expand a polyA-site SAF with gene-region features (5'UTR, CDS, intron, 3'UTR, non-coding)
#'
#' Region features that do not overlap a polyA site are added with GeneIDs of
#' the form `gene{D}NA{D}NA{D}chrom{D}NA{D}strand{D}<nc|cds|intron|utr3|utr5>`,
#' using the same delimiter `{D}` as the input SAF (';' human / '_' mouse) so
#' the output is consistently delimited. Requires `parse_saf` / `parse_gtf`
#' (scraps_to_seurat.R), `parse_geneid` (geneid.R) and valr.
#'
#' @param saf_file polyA-site SAF (scraps 7-field GeneID encoding)
#' @param gtf_file GTF annotation matching the SAF's genome build
#' @param save_file output SAF path
#' @param slop bp added on both sides of each polyA-site window before subtraction
#' @import readr dplyr stringr valr
#' @return invisibly, the written table
#' @export
expand_region_saf <- function(saf_file, gtf_file, save_file, slop = 0) {
  saf_names <- read_tsv(saf_file, n_max = 2) %>% colnames()
  print(saf_names)
  saf <- parse_saf(saf_file)
  delim <- geneid_delim(saf$GeneID)
  region_id <- function(name, chrom, strand, type) {
    str_c(name, "NA", "NA", chrom, "NA", strand, type, sep = delim)
  }
  gtf <- parse_gtf(gtf_file)
  
  if (slop != 0) {
    saf <- saf %>% mutate(start = start - slop, end = end + slop)
  }
  
  # non-coding
  ncgenes <- setdiff(gtf %>% filter(type == "gene") %>% pull(name), 
                gtf %>% filter(type == "CDS") %>% pull(name))
  nc <- gtf %>% filter(name %in% ncgenes, type == "exon")
  nc_notpolya <- bed_subtract(nc, saf)
  nc_notpolya2 <- nc_notpolya %>% 
    mutate(GeneID = region_id(name, chrom, strand, "nc")) %>% 
    rename(gene = "name") %>% 
    dplyr::select(colnames(saf))
  
  # cds
  cds <- gtf %>% filter(type == "CDS")
  cds_notpolya <- bed_subtract(cds, saf)
  cds_notpolya2 <- cds_notpolya %>% 
    mutate(GeneID = region_id(name, chrom, strand, "cds")) %>% 
    rename(gene = "name") %>% 
    dplyr::select(colnames(saf))
  
  # intron
  intron <- bed_subtract(gtf %>% filter(type == "gene"),
                         gtf %>% filter(type == "exon"))
  intron_notpolya <- bed_subtract(intron, saf)
  intron_notpolya2 <- intron_notpolya %>% 
    mutate(GeneID = region_id(name, chrom, strand, "intron")) %>% 
    rename(gene = "name") %>% 
    dplyr::select(colnames(saf))
  
  # utr3
  utr <- gtf %>% filter(type == "UTR")
  stops <- gtf %>% filter(type == "stop_codon")
  utr3 <- utr %>% left_join(stops, by = "transcript") %>%
    filter((strand.x == "+" & start.x >= start.y) | (strand.x == "-" & end.x <= end.y)) %>% 
    dplyr::select(name = name.x, chrom = chrom.x, strand = strand.x, start = start.x, end = end.x)
  utr3_notpolya <- bed_subtract(utr3, saf)
  utr3_notpolya2 <- utr3_notpolya %>% 
    mutate(GeneID = region_id(name, chrom, strand, "utr3")) %>% 
    rename(gene = "name") %>% 
    dplyr::select(colnames(saf))
  
  # utr5
  utr5 <- utr %>% left_join(stops, by = "transcript") %>%
    filter((strand.x == "+" & start.x <= start.y) | (strand.x == "-" & end.x >= end.y)) %>% 
    dplyr::select(name = name.x, chrom = chrom.x, strand = strand.x, start = start.x, end = end.x)
  utr5_notpolya <- bed_subtract(utr5, saf)
  utr5_notpolya2 <- utr5_notpolya %>% 
    mutate(GeneID = region_id(name, chrom, strand, "utr5")) %>% 
    rename(gene = "name") %>% 
    dplyr::select(colnames(saf))
  
  # out
  full <- bind_rows(saf, nc_notpolya2, cds_notpolya2, utr3_notpolya2, utr5_notpolya2, intron_notpolya2) %>% 
    bed_sort() %>% 
    mutate(start = start + 1) %>% 
    dplyr::rename(Chr = chrom,
                  Start = start,
                  End = end,
                  Strand = strand) %>% 
    dplyr::select(saf_names)
  
  write_tsv(full, save_file)
  invisible(full)
}