#' Calculate comparative statistics for polyA lengths between two groups
#'
#' \code{calculate_statistics} computes per-entity summary and statistical
#' comparisons of polyA tail lengths between a specified control and treated
#' group. For each unique entity (e.g., gene or transcript) it:
#'   1. Counts observations in control and treated groups.
#'   2. If both groups have ≥2 observations, performs a Wilcoxon rank-sum test.
#'   3. Calculates group means, fold change, mean difference, Cohen’s d, and
#'      categorizes effect size.
#'   4. Adjusts p-values across all entities using the specified method.
#'
#' @param polyA_table A \code{data.frame} produced by \code{get_polyA()}, which
#'   must contain at least the columns \code{"polyA_length"}, the grouping factor,
#'   and the entity identifier. Defaults to \code{get_gene_id_out}.
#' @param grouping_factor A string giving the column name in
#'   \code{polyA_table} that defines group membership (e.g., treatment vs control).
#'   Default \code{"group"}.
#' @param which_level A string giving the column name in \code{polyA_table} that
#'   defines the entity level for aggregation (e.g., \code{"gene_id"} or
#'   \code{"transcript_id"}). Default \code{"gene_id"}.
#' @param control_group A string naming the level of \code{grouping_factor} to
#'   be treated as the control group. Must be provided.
#' @param treated_group A string naming the level of \code{grouping_factor} to
#'   be treated as the treated group. Must be provided.
#' @param padj_method A string specifying the p-value adjustment method passed
#'   to \code{\link[stats]{p.adjust}} (e.g., \code{"fdr"}, \code{"BH"},
#'   \code{"bonferroni"}). Default \code{"fdr"}.
#'
#' @author Mateusz Mazdziarz
#'
#' @importFrom stats wilcox.test p.adjust
#' @import data.table

calculate_statistics <- function(polyA_table = get_gene_id_out, 
                                      grouping_factor = "group", 
                                      which_level = "gene_id", 
                                      control_group = NULL, 
                                      treated_group = NULL, 
                                      padj_method = "fdr") {
  
  if (missing(polyA_table)) stop("'polyA_table' must be defined.")
  if (is.null(control_group) | is.null(treated_group)) {
    stop("Both 'control_group' and 'treated_group' must be provided.")
  }
  
  message("Starting processing...")
  start_time <- Sys.time()
  
  dt <- as.data.table(polyA_table)
  `.` <- list
  dt_sub <- dt[get(grouping_factor) %in% c(control_group, treated_group)]
  
  setnames(dt_sub, c(grouping_factor, which_level), c("grp_tmp", "unit_tmp"))
  
  results <- dt_sub[, {
    ctr_data <- polyA_length[grp_tmp == control_group]
    trt_data <- polyA_length[grp_tmp == treated_group]
    
    n_ctr <- length(ctr_data)
    n_trt <- length(trt_data)

    res <- list(
      p_value = NA_real_, mean_group_ctr = NA_real_, mean_group_trt = NA_real_,
      diff_length = NA_real_, fold_change = NA_real_, cohen_d = NA_real_,
      cohen_effect = NA_character_
    )
    
    if (n_ctr >= 2 && n_trt >= 2) {
      m_ctr <- mean(ctr_data, na.rm = TRUE)
      m_trt <- mean(trt_data, na.rm = TRUE)
      v_ctr <- var(ctr_data, na.rm = TRUE)
      v_trt <- var(trt_data, na.rm = TRUE)
      
      res$p_value <- suppressWarnings(wilcox.test(ctr_data, trt_data)$p.value)
      
      res$mean_group_ctr <- m_ctr
      res$mean_group_trt <- m_trt
      res$diff_length    <- m_trt - m_ctr
      res$fold_change    <- if(m_ctr != 0) m_trt / m_ctr else NA_real_
      
      pooled_sd <- sqrt(((n_ctr - 1) * v_ctr + (n_trt - 1) * v_trt) / (n_ctr + n_trt - 2))
      if (!is.na(pooled_sd) && pooled_sd > 0) {
        cd <- abs((m_ctr - m_trt) / pooled_sd)
        res$cohen_d <- cd
        res$cohen_effect <- if(cd < 0.2) "small" else if(cd < 0.5) "medium" else "large"
      }
    }
    
    res
  }, by = unit_tmp]
  
  setnames(results, "unit_tmp", which_level)
  results[, padj := stats::p.adjust(p_value, method = padj_method)]
  results[, Log2FC := log2(fold_change)]
  
  message(sprintf("Processing completed in %.2f seconds.", 
                  as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
  
  return(as.data.frame(results))
}
