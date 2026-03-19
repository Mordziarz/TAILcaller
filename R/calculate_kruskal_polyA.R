#' Calculate Kruskal–Wallis p-values per entity for polyA lengths
#'
#' \code{calculate_kruskal_polyA} computes, for each unique entity in a
#' polyA length table, the Kruskal–Wallis test p-value comparing
#' \code{polyA_length} distributions across groups. Only entities with
#' more than two groups are tested; others receive \code{NA}.
#'
#' @param polyA_table A \code{data.frame} containing at least the columns
#'   \code{"polyA_length"}, the grouping factor, and the entity identifier.
#'   Defaults to \code{get_gene_id_out}.
#' @param grouping_factor A string naming the column in \code{polyA_table}
#'   that defines the groups for comparison. Default \code{"group"}.
#' @param which_level A string naming the column in \code{polyA_table}
#'   that defines the entity (e.g., gene, transcript). Default
#'   \code{"gene_id"}.
#' @param padj_method A string specifying the p-value adjustment method
#'   passed to \code{\link[stats]{p.adjust}}. Default \code{"fdr"}.
#'
#' @author Mateusz Mazdziarz
#'
#' @importFrom stats kruskal.test p.adjust
#' @import data.table

calculate_kruskal_polyA <- function(polyA_table = get_gene_id_out, 
                                         grouping_factor = "group", 
                                         which_level = "gene_id", 
                                         padj_method = "fdr") {
  
  if (missing(polyA_table)) stop("'polyA_table' must be defined.")
  
  message("Starting to process the data...")
  start_time <- Sys.time()
  
  dt <- as.data.table(polyA_table)
  `.` <- list
  setnames(dt, c(grouping_factor, which_level), c("grp_tmp", "unit_tmp"))
  
  results <- dt[, {
    unique_grps <- uniqueN(grp_tmp)
    
    p_val <- NA_real_
    
    if (unique_grps > 2) {
      p_val <- suppressWarnings(
        stats::kruskal.test(polyA_length ~ grp_tmp, data = .SD)$p.value
      )
    }
    
    list(p_value = p_val)
    
  }, by = unit_tmp]
  
  setnames(results, "unit_tmp", which_level)
  results[, padj := stats::p.adjust(p_value, method = padj_method)]
  
  setnames(dt, c("grp_tmp", "unit_tmp"), c(grouping_factor, which_level))
  
  end_time <- Sys.time()
  message(sprintf("Processing complete. Time taken: %.2f seconds.", 
                  as.numeric(difftime(end_time, start_time, units = "secs"))))
  
  return(as.data.frame(results))
}