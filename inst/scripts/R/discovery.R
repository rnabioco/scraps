#' Load de novo RT priming site annotations
#'
#' Reads the annotated sites table produced by the scraps discovery module
#' (`{results}/discovery/<unit>_sites.tsv.gz`, from `annotate_sites.py`).
#'
#' @param file path to a `*_sites.tsv.gz` discovery output
#' @param status optional character vector of statuses to keep
#'   (any of `known_pas`, `likely_internal_priming`, `potential_novel_pas`).
#'   Defaults to all.
#' @param min_umi minimum `umi_support` to retain a site (default 0)
#' @import readr dplyr
#' @return tibble of annotated RT priming sites
#' @export
read_discovery_sites <- function(file, status = NULL, min_umi = 0) {
  sites <- readr::read_tsv(
    file,
    col_types = readr::cols(
      chrom = readr::col_character(),
      summit = readr::col_integer(),
      strand = readr::col_character(),
      kde_score = readr::col_double(),
      umi_support = readr::col_double(),
      status = readr::col_character(),
      gene = readr::col_character(),
      id = readr::col_character(),
      class = readr::col_character(),
      polya_dist = readr::col_character(),
      pas_motif = readr::col_character(),
      pas_dist = readr::col_character(),
      a_content = readr::col_double(),
      a_run = readr::col_integer()
    )
  )
  sites <- dplyr::filter(sites, umi_support >= min_umi)
  if (!is.null(status)) {
    sites <- dplyr::filter(sites, status %in% !!status)
  }
  sites
}

#' Summarize discovered RT priming sites by status
#'
#' @param sites tibble from \code{read_discovery_sites}
#' @import dplyr
#' @return tibble with per-status counts and total UMI support
#' @export
summarize_discovery_sites <- function(sites) {
  sites %>%
    dplyr::group_by(status) %>%
    dplyr::summarize(
      n_sites = dplyr::n(),
      total_umi = sum(umi_support),
      median_kde = stats::median(kde_score),
      .groups = "drop"
    )
}
