#' Parse scraps GeneIDs (count-table `gene` or SAF `GeneID` values)
#'
#' A scraps GeneID encodes seven fields,
#' `gene{D}refseq{D}ensembl{D}chrom{D}pos{D}strand{D}type`, with delimiter `{D}`
#' = `;` (human references, the de novo discovery SAF) or `_` (mouse
#' references, the region SAF built by `expand_region_saf`). Count-table ids
#' may join several overlapping features with `,` (featureCounts `-O`). The
#' delimiter is detected per feature. `_`-delimited ids are parsed from the
#' right because `_` also occurs inside RefSeq ids (`NR_152944`), scaffold /
#' alt-contig names (`chr22_KI270879v1_alt`) and some gene symbols
#' (`APOBEC3A_B`); this mirrors `_split_geneid` in
#' `inst/scripts/annotate_sites.py`.
#'
#' @param x character vector of GeneIDs
#' @param multi `"first"` keeps the first comma-joined feature of each id (one
#'   row per input); `"all"` expands every feature (one row per feature, with
#'   `feature_idx`)
#' @param chrom optional SAF `Chr` column (same length as `x`) used as the
#'   chromosome anchor for `_` ids
#' @param strand optional SAF `Strand` column, validated against the parsed
#'   strand when given
#' @import dplyr stringr tibble
#' @return tibble with columns `GeneID` (input, unchanged), `feature_idx`,
#'   `delim`, `gene`, `refseq`, `ensembl`, `chrom`, `pos` (numeric; `NA` for
#'   `"NA"`), `strand`, `type`. Unparseable features give `NA` fields and one
#'   aggregated warning.
#' @export
parse_geneid <- function(x, multi = c("first", "all"), chrom = NULL, strand = NULL) {
  multi <- match.arg(multi)
  x <- as.character(x)
  feats <- stringr::str_split(x, ",")
  if (multi == "first") {
    feat <- vapply(feats, function(f) f[1], character(1))
    src <- seq_along(x)
    idx <- rep(1L, length(x))
  } else {
    lens <- lengths(feats)
    feat <- unlist(feats, use.names = FALSE)
    src <- rep(seq_along(x), lens)
    idx <- unlist(lapply(lens, seq_len), use.names = FALSE)
  }
  anchor_chrom <- if (is.null(chrom)) rep(NA_character_, length(x)) else as.character(chrom)
  anchor_strand <- if (is.null(strand)) rep(NA_character_, length(x)) else as.character(strand)

  n <- length(feat)
  out <- list(
    delim = rep(NA_character_, n), gene = rep(NA_character_, n),
    refseq = rep(NA_character_, n), ensembl = rep(NA_character_, n),
    chrom = rep(NA_character_, n), pos = rep(NA_character_, n),
    strand = rep(NA_character_, n), type = rep(NA_character_, n)
  )

  # ';' features: fixed 7-field split (vectorised). Missing ids stay NA rows.
  is_semi <- !is.na(feat) & stringr::str_detect(feat, stringr::fixed(";"))
  if (any(is_semi)) {
    ok <- is_semi & stringr::str_count(feat, stringr::fixed(";")) == 6L
    parts <- stringr::str_split_fixed(feat[ok], stringr::fixed(";"), 7)
    out$delim[ok] <- ";"
    out$gene[ok] <- parts[, 1]; out$refseq[ok] <- parts[, 2]; out$ensembl[ok] <- parts[, 3]
    out$chrom[ok] <- parts[, 4]; out$pos[ok] <- parts[, 5]; out$strand[ok] <- parts[, 6]
    out$type[ok] <- parts[, 7]
  }

  # '_' features: read from the right, per element
  is_us <- !is_semi & !is.na(feat) & feat != ""
  for (i in which(is_us)) {
    r <- .parse_geneid_underscore(feat[i], anchor_chrom[src[i]])
    if (is.null(r)) next
    out$delim[i] <- "_"
    out$gene[i] <- r$gene; out$refseq[i] <- r$refseq; out$ensembl[i] <- r$ensembl
    out$chrom[i] <- r$chrom; out$pos[i] <- r$pos; out$strand[i] <- r$strand
    out$type[i] <- r$type
  }

  bad <- is.na(out$type) & !is.na(feat)
  if (any(bad)) {
    warning(sum(bad), " GeneID feature(s) not parseable as 7 fields, e.g. ",
            paste(utils::head(unique(feat[bad]), 3), collapse = " | "))
  }
  mismatch <- !is.na(anchor_strand[src]) & !is.na(out$strand) & out$strand != anchor_strand[src]
  if (any(mismatch)) {
    warning(sum(mismatch), " GeneID(s) whose strand field disagrees with the SAF Strand column")
  }

  tibble::tibble(
    GeneID = x[src],
    feature_idx = idx,
    delim = out$delim,
    gene = out$gene,
    refseq = out$refseq,
    ensembl = out$ensembl,
    chrom = out$chrom,
    pos = suppressWarnings(as.numeric(ifelse(out$pos == "NA", NA_character_, out$pos))),
    strand = out$strand,
    type = out$type
  )
}

#' @keywords internal
.parse_geneid_underscore <- function(feat, anchor_chrom = NA_character_) {
  t <- strsplit(feat, "_", fixed = TRUE)[[1]]
  n <- length(t)
  if (n < 7) return(NULL)
  type <- t[n]; strand <- t[n - 1]; pos <- t[n - 2]
  if (!strand %in% c("+", "-", ".")) return(NULL)
  cs <- NA_integer_
  if (!is.na(anchor_chrom)) {
    ct <- strsplit(anchor_chrom, "_", fixed = TRUE)[[1]]
    k <- length(ct)
    cand <- n - 2 - k
    if (cand >= 4 && identical(t[cand:(n - 3)], ct)) cs <- cand
  }
  if (is.na(cs)) {
    hits <- which(startsWith(t, "chr"))
    hits <- hits[hits >= 4 & hits <= n - 3]
    cs <- if (length(hits)) hits[1] else n - 3L
    if (cs < 4) cs <- n - 3L
  }
  chrom <- paste(t[cs:(n - 3)], collapse = "_")
  L <- t[seq_len(cs - 1)]
  m <- length(L)
  if (m < 3) {
    gene <- paste(L, collapse = "_"); refseq <- NA_character_; ensembl <- NA_character_
  } else {
    ensembl <- L[m]
    if (m >= 4 && grepl("^[NX][MR]$", L[m - 2]) && grepl("^[0-9]+(\\.[0-9]+)?$", L[m - 1])) {
      refseq <- paste(L[m - 2], L[m - 1], sep = "_")
      gene <- paste(L[seq_len(m - 3)], collapse = "_")
    } else {
      refseq <- L[m - 1]
      gene <- paste(L[seq_len(m - 2)], collapse = "_")
    }
  }
  list(gene = gene, refseq = refseq, ensembl = ensembl, chrom = chrom,
       pos = pos, strand = strand, type = type)
}

#' Delimiter used by a set of GeneIDs
#'
#' @param x character vector of GeneIDs
#' @return `";"` if any id contains `;`, else `"_"`
#' @export
geneid_delim <- function(x) {
  if (any(stringr::str_detect(as.character(x), stringr::fixed(";")), na.rm = TRUE)) ";" else "_"
}
