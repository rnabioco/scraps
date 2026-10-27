#' Load de novo RT priming site annotations
#'
#' Reads the annotated sites table produced by the scraps discovery module
#' (`{results}/discovery/<unit>_sites.tsv.gz`, from `annotate_sites.py`).
#'
#' @param file path to a `*_sites.tsv.gz` discovery output
#' @param status optional character vector of statuses to keep (any of
#'   `known_pas`, `suspected_internal_priming`, `potential_novel_pas`,
#'   `potential_novel_no_pas`). Defaults to all. Tables written by older scraps
#'   versions use `likely_internal_priming`; it is read as
#'   `suspected_internal_priming` (re-run annotate_sites.py to split off the
#'   `potential_novel_no_pas` sites that the old label also covered).
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
  # legacy label (scraps < this release) -> current name
  sites <- dplyr::mutate(
    sites,
    status = dplyr::recode(status, likely_internal_priming = "suspected_internal_priming"),
    class = dplyr::recode(class, likely_internal_priming = "suspected_internal_priming")
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
