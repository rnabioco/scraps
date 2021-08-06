expand_region_saf <- function(saf_file, gtf_file, save_file, slop = 0) {
  saf_names <- read_tsv(saf_file, n_max = 2) %>% colnames()
  print(saf_names)
  saf <- parse_saf(saf_file)
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
    mutate(GeneID = str_c(name, "NA", "NA", chrom, "NA", strand, "nc", sep = "_")) %>% 
    rename(gene = "name") %>% 
    dplyr::select(colnames(saf))
  
  # cds
  cds <- gtf %>% filter(type == "CDS")
  cds_notpolya <- bed_subtract(cds, saf)
  cds_notpolya2 <- cds_notpolya %>% 
    mutate(GeneID = str_c(name, "NA", "NA", chrom, "NA", strand, "cds", sep = "_")) %>% 
    rename(gene = "name") %>% 
    dplyr::select(colnames(saf))
  
  # intron
  intron <- bed_subtract(gtf %>% filter(type == "gene"),
                         gtf %>% filter(type == "exon"))
  intron_notpolya <- bed_subtract(intron, saf)
  intron_notpolya2 <- intron_notpolya %>% 
    mutate(GeneID = str_c(name, "NA", "NA", chrom, "NA", strand, "intron", sep = "_")) %>% 
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
    mutate(GeneID = str_c(name, "NA", "NA", chrom, "NA", strand, "utr3", sep = "_")) %>% 
    rename(gene = "name") %>% 
    dplyr::select(colnames(saf))
  
  # utr5
  utr5 <- utr %>% left_join(stops, by = "transcript") %>%
    filter((strand.x == "+" & start.x <= start.y) | (strand.x == "-" & end.x >= end.y)) %>% 
    dplyr::select(name = name.x, chrom = chrom.x, strand = strand.x, start = start.x, end = end.x)
  utr5_notpolya <- bed_subtract(utr5, saf)
  utr5_notpolya2 <- utr5_notpolya %>% 
    mutate(GeneID = str_c(name, "NA", "NA", chrom, "NA", strand, "utr5", sep = "_")) %>% 
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
}

expand_region_saf("/Users/rf/scraps/ref/polyadb32.hg38.saf.gz", 
                  "/Users/rf/gencode.v33.annotation.gtf", 
                  "polyadb32.hg38.expand2.saf.gz", 
                  slop = 2)
