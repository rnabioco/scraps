#' Read scraps output to read counts by region type
#' 
#' @param files scraps output table files containing 
#' @param cell_ids if given, use only these cell barcodes, and fill in empty ones
#' @param batch batch prefix to use in order for list of files
#' @return count list for regions
#' @export
#' 
scraps_to_region_counts <- function(files,
                                    cell_ids = NULL,
                                    batch = NULL,
                                    batch_sep = "-") {
  ids <- cell_ids
  if (!is.null(batch)) {
    temp <- map(seq_along(files),
                function(x) read_tsv(files[x]) %>% mutate(cell = str_c(batch[x], cell, sep = batch_sep)))
  } else {
    temp <- map(seq_along(files),
                function(x) read_tsv(files[x]))
  }
  
  temp <- do.call(rbind, temp)
  
  if (!is.null(ids)) {
    temp <- temp %>% filter(cell %in% ids)
  }
  
  temp2 <- temp %>%
    separate(gene, 
             sep = "_", 
             into = c("gene", NA,
                      NA, NA, NA, NA,
                      "pos"))
  temp3 <- temp2 %>%
    separate(pos, 
             sep = ",",
             into = c("pos", NA))
  temp4 <- temp3 %>% mutate(pos2 = case_when(
    str_detect(pos, "3'UTR") ~ "utr3",
    pos == "utr3" ~ "utr3", 
    pos == "CDS" ~ "cds",
    pos == "cds" ~ "cds",
    pos == "intron" ~ "intron",
    pos == "Intron" ~ "intron",
    TRUE ~ "utr5")) %>% 
    group_by(gene, pos2, cell) %>% 
    summarise(count = sum(count))
  
  cds <- temp4 %>% ungroup() %>%
    filter(pos2 == "cds") %>%
    dplyr::select(-pos2) %>% pivot_wider(names_from = cell, values_from = count, values_fill = 0) %>% 
    column_to_rownames("gene")
  
  utr3 <- temp4 %>% ungroup() %>%
    filter(pos2 == "utr3") %>%
    dplyr::select(-pos2) %>% pivot_wider(names_from = cell, values_from = count, values_fill = 0) %>% 
    column_to_rownames("gene")
  
  utr5 <- temp4 %>% ungroup() %>%
    filter(pos2 == "utr5") %>%
    dplyr::select(-pos2) %>% pivot_wider(names_from = cell, values_from = count, values_fill = 0) %>% 
    column_to_rownames("gene")
  
  intron <- temp4 %>% ungroup() %>%
    filter(pos2 == "intron") %>%
    dplyr::select(-pos2) %>% pivot_wider(names_from = cell, values_from = count, values_fill = 0) %>% 
    column_to_rownames("gene")
  
  l <- list()
  l$utr5 <- utr5
  l$utr3 <- utr3
  l$cds <- cds
  l$intron <- intron
  
  l
}

#' Fill matrix with 0 rows or columns as needed
#' 
#' @param x target matrix
#' @param vec target row(or column) names to ensure exists and in right order
#' @param row if TRUE, target rows, if FALSE, columns
#' @return 0 filled and ordered matrix
#' @export
#' 
fill_mat <- function(x, vec, row = TRUE) {
  if (row) {
    emptys <- setdiff(vec, rownames(x))
    
    if (length(emptys) > 0) {
      empty_mat <- matrix(data = 0, ncol = ncol(x), nrow = length(emptys))
      rownames(empty_mat) <- emptys
      colnames(empty_mat) <- colnames(x)
      x <- rbind(x, empty_mat)
      x <- x[vec, ]
    } else {
      x <- x[vec, ]
    }
    x
  } else {
    emptys <- setdiff(vec, colnames(x))
    
    if (length(emptys) > 0) {
      empty_mat <- matrix(data = 0, nrow = nrow(x), ncol = length(emptys))
      colnames(empty_mat) <- emptys
      rownames(empty_mat) <- rownames(x)
      x <- cbind(x, empty_mat)
      x <- x[, vec]
    } else {
      x <- x[, vec]
    }
    x
  }
}

#' filter region count list and convert to matrices
#' 
#' @param countlist output from scraps_to_region_counts
#' @param n_min  minimum number of observations for a gene to be kept
#' @return count list for regions
#' @export
#' 
filter_countlist <- function(countlist,
                             n_min = 100) {
  shared_genes <- sapply(countlist, rownames, simplify = F) %>% unlist() %>% table()
  keep_gene <- shared_genes[shared_genes >= 2] %>% names()
  
  countlist$intron <- fill_mat(countlist$intron, keep_gene)
  countlist$cds <- fill_mat(countlist$cds, keep_gene)
  countlist$utr3 <- fill_mat(countlist$utr3, keep_gene)
  countlist$utr5 <- fill_mat(countlist$utr5, keep_gene)
  
  sums <- sapply(countlist,
                 function(x) colSums(x) %>% as.data.frame() %>% rownames_to_column("cell"), 
                 simplify = F)
  sums <- do.call(rbind, sums) %>% as.data.frame() %>% 
    group_by(cell) %>%
    summarize(total = sum(`.`))
  keep_ids <- sums %>% filter(total >= n_min) %>% pull(cell)
  
  countlist$intron <- fill_mat(countlist$intron, keep_ids, row = FALSE)
  countlist$cds <- fill_mat(countlist$cds, keep_ids, row = FALSE)
  countlist$utr3 <- fill_mat(countlist$utr3, keep_ids, row = FALSE)
  countlist$utr5 <- fill_mat(countlist$utr5, keep_ids, row = FALSE)
  
  countlist
}

#' filter region count list and convert to matrices
#' 
#' @param df one df from output from scraps_to_region_counts
#' @param dfgene ranked gene data.frame
#' @param vec vector of cluster identities
#' @param utr5 counts for region
#' @param utr3 counts for region
#' @param cds counts for region
#' @param intron counts for region
#' @param k number of bins
#' @return percentage binned averages for specific region
#' @export
#' 
dfnot_avg <- function(df, dfgene, vec, k,
                      utr5, utr3, cds, intron) {
  sapply(1:k,
         function(x) {
           message(x)
           genes <- dfgene %>% rownames_to_column("gene") %>% 
             filter(quantile == x) %>%
             pull(gene)
           
           a0 <- clustifyr::average_clusters(
             data.frame(as.list(colMeans(df[genes, ]))),
             vec,
             if_log = FALSE,
             output_log = FALSE) %>% 
             t() %>% 
             as.data.frame() 
           
           a1 <- clustifyr::average_clusters(
             data.frame(as.list(colMeans(intron[genes, ]))),
             vec,
             if_log = FALSE,
             output_log = FALSE) %>% 
             t() %>% 
             as.data.frame() 
           
           a2 <- clustifyr::average_clusters(
             data.frame(as.list(colMeans(utr3[genes, ]))),
             vec,
             if_log = FALSE,
             output_log = FALSE) %>% 
             t() %>% 
             as.data.frame() 
           
           a3 <- clustifyr::average_clusters(
             data.frame(as.list(colMeans(cds[genes, ]))),
             vec,
             if_log = FALSE,
             output_log = FALSE) %>% 
             t() %>% 
             as.data.frame() 
           
           a5 <- clustifyr::average_clusters(
             data.frame(as.list(colMeans(utr5[genes, ]))),
             vec,
             if_log = FALSE,
             output_log = FALSE) %>% 
             t() %>% 
             as.data.frame() 
           
           (a0)/(a1+a2+a3+a5)
         }, simplify = FALSE)
}

#' calculations from countlist, ready to plot
#' 
#' @param countlist output from scraps_to_region_counts and filter_countlist 
#' @param vec vector of cluster identities
#' @param genes if provided, only calculates based on these genes
#' @param k quantiles to use, binned by 2000 genes if not given
#' @return data.frame with count bins
#' @export
#' 
countlist_to_bins <- function(countlist,
                              vec,
                              genes = NA,
                              k = NA) {
  if (is.na(genes)) {
    if (is.na(k)) {
      k <- ceiling(nrow(countlist$utr3)/2000)
      message("binning genes by 2000")
    } else {
      message("using ", k, " bins")
    }
    utr3 <- countlist$utr3
    utr5 <- countlist$utr5
    intron <- countlist$intron
    cds <- countlist$cds
  } else {
    k <- 1
    genes <- intersect(rownames(countlist$utr3), genes)
    message("custom list of ", length(genes), " genes")
    utr3 <- countlist$utr3[genes, ]
    utr5 <- countlist$utr5[genes, ]
    intron <- countlist$intron[genes, ]
    cds <- countlist$cds[genes, ]
  }

  dfgene <- data.frame(list(rowSums(utr3) + rowSums(cds) + rowSums(intron) + rowSums(utr5))) %>%
    setNames("n") %>%
    mutate(quantile = ntile(n, k))
  
  dfnot2_cds <- do.call(cbind, 
                        dfnot_avg(cds, dfgene = dfgene, vec = vec, k = k,
                                  utr5 = utr5, utr3 = utr3, cds = cds, intron = intron)) %>% 
    rownames_to_column("id") %>% 
    setNames(c("id", 1:k)) %>% 
    pivot_longer(-id, names_to = "quantile", values_to = "frac") %>% 
    mutate(quantile = as.numeric(quantile)) %>% 
    mutate(loc = "cds")
  
  dfnot2_intron <- do.call(cbind, 
                           dfnot_avg(intron, dfgene = dfgene, vec = vec, k = k,
                                     utr5 = utr5, utr3 = utr3, cds = cds, intron = intron)) %>% 
    rownames_to_column("id") %>% 
    setNames(c("id", 1:k)) %>% 
    pivot_longer(-id, names_to = "quantile", values_to = "frac") %>% 
    mutate(quantile = as.numeric(quantile)) %>% 
    mutate(loc = "intron")
  
  dfnot2_utr3 <- do.call(cbind, 
                         dfnot_avg(utr3, dfgene = dfgene, vec = vec, k = k,
                                   utr5 = utr5, utr3 = utr3, cds = cds, intron = intron)) %>%
    rownames_to_column("id") %>%
    setNames(c("id", 1:k)) %>%
    pivot_longer(-id, names_to = "quantile", values_to = "frac") %>%
    mutate(quantile = as.numeric(quantile)) %>%
    mutate(loc = "utr3")
  
  dfnot2_utr5 <- do.call(cbind, 
                         dfnot_avg(utr5, dfgene = dfgene, vec = vec, k = k,
                                   utr5 = utr5, utr3 = utr3, cds = cds, intron = intron)) %>% 
    rownames_to_column("id") %>% 
    setNames(c("id", 1:k)) %>% 
    pivot_longer(-id, names_to = "quantile", values_to = "frac") %>% 
    mutate(quantile = as.numeric(quantile)) %>% 
    mutate(loc = "utr5")
  
  dfnot2 <- bind_rows(dfnot2_utr5, dfnot2_utr3, dfnot2_cds, dfnot2_intron) %>% 
    mutate(loc = factor(loc, levels = c("utr5", "cds", "intron", "utr3")))
}

#' calculation cohens_d effect size
#' 
#' @param x
#' @param y
#' @return cohens_d effect size
#' @export
#' 
cohens_d <- function(x, y) {
  lx <- length(x)- 1
  ly <- length(y)- 1
  md  <- abs(mean(x) - mean(y))        ## mean difference (numerator)
  csd <- lx * var(x) + ly * var(y)
  csd <- csd/(lx + ly)
  csd <- sqrt(csd)                     ## common sd computation
  if (csd == 0) {
    return(0)
  }
  cd  <- md/csd                        ## cohen's d
}

#' calculation mean diff and effect size of utr3 priming
#' 
#' @param df data.frame with count bins
#' @param target_loc string for target region, default utr3
#' @param col1 column designating broader grouping if needed
#' @param col2 column designating grouping, default to "id" from countlist_to_bins
#' @param ctrl target control set in col
#' @return data.frame with stats (mean diff and effecti size)
#' @export
#' 
bin_effectsize <- function(df,
                           target_loc = "utr3",
                           col1 = NA,
                           col2 = "id",
                           ctrl = NA) {
  temp <- df %>% mutate(quantile = as.numeric(quantile)) %>%
    filter(quantile != min(quantile), quantile != max(quantile)) %>% 
    filter(loc == target_loc) 
  if (is.na(col1)) {
    groupvar <- "quantile"
  } else {
    groupvar <- c(col1, "quantile")
  }
  if (is.na(col1)) {
    groupvar2 <- col2
  } else {
    groupvar2 <- c(col1, "quantile")
  }
  temp2 <- temp %>% group_by_at(all_of(vars(groupvar))) %>% 
    mutate(target := !!sym(col2)) %>% 
    mutate(score = frac[target == ctrl][1] - frac) %>% 
    group_by_at(all_of(vars(c(col1, col2) %>% na.omit()))) %>% 
    summarize(mean = mean(score),
              effectsize = cohens_d(score, rep(0, length(score))),
              p = t.test(score, rep(0, length(score)))$p.value) %>% 
    ungroup()
  temp2 %>% mutate(p = ifelse(is.nan(p), 1, p))
}

#' calculation mean diff and effect size of utr3 priming
#' 
#' @param countlist countlist
#' @param knn_matrix nearest neightbors in matrix format, can be calculated by enabling return.neighbor in Seurat::FindNeighbors
#' @param cell_ids cell barcodes to use
#' @return countlist after smoothing
#' @export
#' 
knn_smoothing <- function(countlist,
                          knn_matrix,
                          cell_ids) {
  countlist$utr5 <- countlist$utr5[,cell_ids]
  for (cell1 in cell_ids) {
    countlist$utr5[, cell1] <- countlist$utr5[, knn_matrix[cell1,]] %>% rowSums()
  }
  
  countlist$utr3 <- countlist$utr3[,cell_ids]
  for (cell1 in cell_ids) {
    countlist$utr3[, cell1] <- countlist$utr3[, knn_matrix[cell1,]] %>% rowSums()
  }
  
  message("intron...")
  countlist$intron <- countlist$intron[,cell_ids]
  for (cell1 in cell_ids) {
    countlist$intron[, cell1] <- countlist$intron[, knn_matrix[cell1,]] %>% rowSums()
  }
  
  countlist$cds <- countlist$cds[,cell_ids]
  for (cell1 in cell_ids) {
    countlist$cds[, cell1] <- countlist$cds[, knn_matrix[cell1,]] %>% rowSums()
  }
  
  return(countlist)
}