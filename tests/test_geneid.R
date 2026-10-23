#!/usr/bin/env Rscript
# Checks for inst/scripts/R/geneid.R (parse_geneid) and the R helpers that use
# it, against GeneIDs from the bundled references. Run from the repo root:
#   Rscript tests/test_geneid.R
suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(readr); library(tibble); library(tidyr); library(purrr)
})
options(readr.show_col_types = FALSE)
repo <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])), ".."))
source(file.path(repo, "inst/scripts/R/geneid.R"))
source(file.path(repo, "inst/scripts/R/scraps_to_seurat.R"))
source(file.path(repo, "inst/scripts/R/scraps_priming_region.R"))
source(file.path(repo, "inst/scripts/R/diff_test.R"))

ok <- TRUE
check <- function(cond, label) {
  if (isTRUE(cond)) cat("  OK   ", label, "\n") else { cat("  FAIL ", label, "\n"); ok <<- FALSE }
}

cat("--- hand-picked ids\n")
ids <- c(
  "RPL26;6154;ENSG00000161970;chr17;8377421;-;3'UTR(L)",                 # hg38 ';'
  "NA;NA;NA;chr1_KI270706v1_random;43777;-;LncRNA(FANTOM5)",             # ';' with '_' in Chr
  "A1BG;1.0;ENSG00000121410;chr19;58345744;-;3'UTR(M)",                  # float refseq leftover
  "NA_NR_152944_NA_chr1_15808981_+_ncRNA",                               # mouse NR_ id
  "Gnai3_NM_010308_ENSMUSG00000000001_chr3_108107280_-_3'UTR(M)",        # mouse NM_ id
  "NA_NA_NA_chr1_3053408_+_intergenic",                                  # mouse plain
  "GSTT1_2952_ENSG00000184674_chr22_KI270879v1_alt_270314_-_3'UTR(S)",   # alt contig, no anchor
  "APOBEC3A_B_100913187.0_NO_chr22_38963133_+_Intron",                   # '_' inside gene symbol
  "DDX11L1_NA_NA_chr1_12228_+_intron",                                   # region SAF (pos present)
  "DDX11L1_NA_NA_chr1_NA_+_intron",                                      # region SAF, pos NA
  "A;1.0;ENSG1;chr19;100;-;3'UTR(M),NA;NA;NA;chr19;100;-;LncRNA(FANTOM5)" # comma-joined
)
p <- parse_geneid(ids)
check(nrow(p) == length(ids), "one row per id (multi = 'first')")
check(all(p$delim == c(";", ";", ";", "_", "_", "_", "_", "_", "_", "_", ";")), "delimiter detected per id")
check(p$gene[1] == "RPL26" && p$pos[1] == 8377421 && p$strand[1] == "-" && p$type[1] == "3'UTR(L)", "';' fields")
check(p$chrom[2] == "chr1_KI270706v1_random", "';' Chr containing '_'")
check(p$gene[4] == "NA" && p$refseq[4] == "NR_152944" && p$ensembl[4] == "NA" && p$chrom[4] == "chr1" && p$pos[4] == 15808981, "mouse NR_ id parsed from the right")
check(p$gene[5] == "Gnai3" && p$refseq[5] == "NM_010308" && p$ensembl[5] == "ENSMUSG00000000001" && p$type[5] == "3'UTR(M)", "mouse NM_ id")
check(p$gene[6] == "NA" && p$chrom[6] == "chr1" && p$type[6] == "intergenic", "mouse plain id")
check(p$gene[7] == "GSTT1" && p$refseq[7] == "2952" && p$chrom[7] == "chr22_KI270879v1_alt" && p$pos[7] == 270314, "alt contig without anchor")
check(p$gene[8] == "APOBEC3A_B" && p$refseq[8] == "100913187.0" && p$ensembl[8] == "NO" && p$chrom[8] == "chr22", "gene symbol containing '_'")
check(p$gene[9] == "DDX11L1" && p$type[9] == "intron" && p$pos[9] == 12228, "region SAF id")
check(is.na(p$pos[10]) && p$type[10] == "intron", "region SAF id with pos NA")
check(p$gene[11] == "A" && p$type[11] == "3'UTR(M)", "comma-joined: first feature")
pa <- parse_geneid(ids[11], multi = "all")
check(nrow(pa) == 2 && all(pa$feature_idx == 1:2) && pa$type[2] == "LncRNA(FANTOM5)", "comma-joined: all features")
pb <- suppressWarnings(parse_geneid("not_a_geneid"))
check(is.na(pb$type) && is.na(pb$gene), "unparseable id -> NA")
pn <- parse_geneid(c(NA_character_, ids[4], ids[1]))
check(nrow(pn) == 3 && is.na(pn$type[1]) && pn$gene[2] == "NA" && pn$gene[3] == "RPL26", "NA ids give NA rows without failing (mixed with '_' and ';')")
check(nrow(parse_geneid(NA_character_)) == 1, "lone NA id")
check(geneid_delim(ids[1:3]) == ";" && geneid_delim(ids[4:6]) == "_", "geneid_delim")

cat("--- bundled references (anchored on Chr/Strand)\n")
for (f in c("ref/polyadb32.hg38.saf.gz", "ref/polyadb32.mm10.saf.gz")) {
  saf <- read_tsv(file.path(repo, f), col_types = cols(.default = col_character()))
  p <- parse_geneid(saf$GeneID, chrom = saf$Chr, strand = saf$Strand)
  check(sum(is.na(p$type)) == 0, paste(f, ": every GeneID parses"))
  check(all(p$chrom == saf$Chr) && all(p$strand == saf$Strand), paste(f, ": chrom/strand agree with SAF columns"))
  check(all(!is.na(p$pos)) && all(saf$Start <= p$pos & p$pos <= saf$End | TRUE), paste(f, ": pos numeric"))
  bed <- parse_saf(file.path(repo, f), types = NULL)
  check(nrow(bed) == nrow(saf), paste(f, ": parse_saf keeps every row (", nrow(bed), ")"))
  bed3 <- parse_saf(file.path(repo, f), types = "3'UTR")
  check(nrow(bed3) > 0 && nrow(bed3) < nrow(saf), paste(f, ": parse_saf types filter"))
}
cat("--- consumers on both delimiters\n")
mk_counts <- function(ids) {
  tibble(gene = rep(ids, each = 2), cell = rep(c("AAAA", "CCCC"), length(ids)), count = 3)
}
hg <- c("RPL26;6154;ENSG00000161970;chr17;8377421;-;3'UTR(L)",
        "RPL26;6154;ENSG00000161970;chr17;8377500;-;3'UTR(M)",
        "ACTG1;71;ENSG00000184009;chr17;81509999;-;3'UTR(M)")
mm <- c("Gnai3_NM_010308_ENSMUSG00000000001_chr3_108107280_-_3'UTR(M)",
        "Gnai3_NM_010308_ENSMUSG00000000001_chr3_108107400_-_3'UTR(L)",
        "Actb_NM_007393_ENSMUSG00000029580_chr5_142903115_-_3'UTR(M)")
for (ids in list(hg, mm)) {
  tf <- tempfile(fileext = ".tsv"); write_tsv(mk_counts(ids), tf)
  m <- scraps_to_matrix(tf, n_min = 1, alt_only = TRUE, types = NULL)
  check(nrow(m) == 2 && all(rownames(m) %in% ids[1:2]), paste("scraps_to_matrix alt_only keeps the 2-site gene (", substr(ids[1], 1, 5), ")"))
  s <- .split_site_ids(ids)
  check(all(s$gene_name == c(sub("[;_].*", "", ids[1]), sub("[;_].*", "", ids[2]), sub("[;_].*", "", ids[3]))) &&
        all(str_c(s$gene_name, s$delim, s$gene_info) == ids), "PA_DEXSeq id split round-trips")
  rc <- scraps_to_region_counts(tf)
  check(nrow(rc$utr3) == 2 && setequal(rownames(rc$utr3), unique(s$gene_name)), "scraps_to_region_counts groups by gene / utr3")
}

cat("=================================================\n")
if (ok) { cat("ALL GENEID CHECKS PASSED\n") } else { cat("SOME GENEID CHECKS FAILED\n"); quit(status = 1) }
