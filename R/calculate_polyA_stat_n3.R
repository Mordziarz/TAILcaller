#' Calculate per-entity group comparison p-values for polyA lengths
#'
#' \code{calculate_polyA_stat_n3} performs, for each unique entity in a 
#' polyA length table, a hypothesis test comparing distributions of 
#' \code{polyA_length} across groups. For each entity it:
#'   1. Checks if there are at least two groups and ≥3 total observations.
#'   2. Tests normality (Shapiro–Wilk or Lilliefors).
#'   3. Tests variance homogeneity (Levene’s test) if each group has ≥2 obs.
#'   4. Chooses one of:
#'      - One-way ANOVA (normal + homogeneous variances);
#'      - Welch ANOVA (normal + unequal variances);
#'      - Kruskal–Wallis (non-normal or insufficient obs).
#'   5. Records the p-value and test used.
#'   6. Adjusts all p-values for multiple testing.
#'
#' @param polyA_table A \code{data.frame} containing at least the columns 
#'   \code{"polyA_length"}, the grouping factor, and the entity identifier.
#'   Defaults to \code{get_gene_id_out}.
#' @param grouping_factor A string naming the column in \code{polyA_table} 
#'   that defines groups for comparison. Default \code{"group"}.
#' @param which_level A string naming the column in \code{polyA_table} 
#'   that defines the entity (e.g., gene, transcript). Default \code{"gene_id"}.
#' @param padj_method A string specifying the p-value adjustment method 
#'   passed to \code{\link[stats]{p.adjust}}. Default \code{"fdr"}.
#'
#' @author Mateusz Mazdziarz
#'
#' @importFrom stats shapiro.test kruskal.test aov oneway.test p.adjust
#' @importFrom nortest lillie.test
#' @importFrom car leveneTest
#' @import data.table

calculate_polyA_stat_n3 <- function(polyA_table = get_gene_id_out, grouping_factor = "group", which_level = "gene_id", padj_method = "fdr") {
  
  if (missing(polyA_table)) stop("'polyA_table' must be defined.")
  
  dt <- as.data.table(polyA_table)
  
  setnames(dt, c(grouping_factor, which_level), c("group_tmp", "unit_tmp"))
  dt[, group_tmp := as.factor(group_tmp)]
  
  message("Starting calculations...")
  start_time <- Sys.time()
  
  out_df <- dt[, {
    n_obs <- .N
    grp_counts <- .SD[, .N, by = group_tmp]$N
    n_groups <- length(grp_counts)
    
    p_val <- NA_real_
    chosen_test <- NA_character_
    
    if (n_groups >= 2) {
      if (n_obs < 3 || any(grp_counts < 2)) {
        p_val <- suppressWarnings(stats::kruskal.test(polyA_length ~ group_tmp, data = .SD)$p.value)
        chosen_test <- "Kruskal–Wallis"
      } else {
        norm_p <- .SD[, {
          x <- polyA_length
          if (length(x) >= 3) {
            if (length(unique(x)) == 1) 0
            else if (length(x) <= 5000) stats::shapiro.test(x)$p.value
            else nortest::lillie.test(x)$p.value
          } else 0
        }, by = group_tmp]$V1
        
        all_normal <- all(norm_p > 0.05, na.rm = TRUE)
        
        lev_p <- tryCatch({
          car::leveneTest(polyA_length ~ group_tmp, data = .SD)[["Pr(>F)"]][1]
        }, error = function(e) NA_real_)
        
        if (all_normal && !is.na(lev_p) && lev_p > 0.05) {
          aov_res <- stats::aov(polyA_length ~ group_tmp, data = .SD)
          p_val <- summary(aov_res)[[1]][["Pr(>F)"]][1]
          chosen_test <- "ANOVA"
        } else if (all_normal) {
          welch_res <- stats::oneway.test(polyA_length ~ group_tmp, data = .SD, var.equal = FALSE)
          p_val <- welch_res$p.value
          chosen_test <- "Welch ANOVA"
        } else {
          p_val <- suppressWarnings(stats::kruskal.test(polyA_length ~ group_tmp, data = .SD)$p.value)
          chosen_test <- "Kruskal–Wallis"
        }
      }
    }
    
    list(test_used = chosen_test, p_value = p_val)
    
  }, by = unit_tmp]
  
  setnames(out_df, "unit_tmp", which_level)
  out_df[, padj := stats::p.adjust(p_value, method = padj_method)]
  
  setnames(dt, c("group_tmp", "unit_tmp"), c(grouping_factor, which_level))
  
  end_time <- Sys.time()
  message(sprintf("Done in %.2f seconds.", as.numeric(difftime(end_time, start_time, units = "secs"))))
  
  return(as.data.frame(out_df))
}