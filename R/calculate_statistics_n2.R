#' Calculate Statistical Differences Between Two Groups for Each Molecule
#'
#' For each unique molecule (defined by `which_level`), this function compares
#' poly(A) tail lengths between a `control_group` and a `treated_group` using
#' appropriate statistical tests:
#' \itemize{
#'   \item \strong{Student's t-test}: If data in both groups exhibit a normal distribution and have homogeneous variances.
#'   \item \strong{Welch's t-test}: If data in both groups exhibit a normal distribution but have unequal variances.
#'   \item \strong{Wilcoxon rank-sum test (Mann-Whitney U test)}: Otherwise (i.e., non-normality, or fewer than 2 observations in any group).
#' }
#' It also calculates mean tail lengths for each group, fold change,
#' difference in length, Cohen's d, and adjusts p-values for multiple comparisons.
#'
#' @param polyA_table A \code{data.frame} with per-read poly(A) tail lengths. Must contain:
#'   \describe{
#'     \item{\code{polyA_length} (numeric)}{Tail length of each read.}
#'     \item{\code{<grouping_factor>} (character or factor)}{Grouping variable for comparison, e.g., experimental condition.}
#'     \item{\code{<which_level>} (character or factor)}{Molecular identifier, e.g., gene_id, transcript_id.}
#'   }
#' @param grouping_factor Character; name of the column in \code{polyA_table} to use for defining groups. Defaults to \code{"group"}.
#' @param which_level Character; name of the column in \code{polyA_table} that uniquely identifies each molecule for which statistics should be calculated. Defaults to \code{"gene_id"}.
#' @param control_group Character; name of the control group level in `grouping_factor`.
#' @param treated_group Character; name of the treated group level in `grouping_factor`.
#' @param padj_method Character; method for p-value adjustment. See `p.adjust.methods` for options. Defaults to \code{"fdr"}.
#'
#' @author Mateusz Mazdziarz
#'
#' @importFrom stats wilcox.test t.test shapiro.test p.adjust sd mean
#' @importFrom nortest lillie.test
#' @importFrom car leveneTest
#' @importFrom dplyr %>% distinct
#' @importFrom rlang sym
#' @import data.table

calculate_statistics_n2 <- function(polyA_table = get_gene_id_out, 
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
  dt_filtered <- dt[get(grouping_factor) %in% c(control_group, treated_group)]
  
  setnames(dt_filtered, c(grouping_factor, which_level), c("grp_tmp", "unit_tmp"))
  dt_filtered[, grp_tmp := factor(grp_tmp, levels = c(control_group, treated_group))]
  
  results <- dt_filtered[, {
    ctr_vals <- polyA_length[grp_tmp == control_group]
    trt_vals <- polyA_length[grp_tmp == treated_group]
    
    n_ctr <- length(ctr_vals)
    n_trt <- length(trt_vals)
    
    res <- list(
      p_value = NA_real_, mean_group_ctr = NA_real_, mean_group_trt = NA_real_,
      diff_length = NA_real_, fold_change = NA_real_, cohen_d = NA_real_,
      cohen_effect = NA_character_, test_performed = "Insufficient data"
    )
    
    if (n_ctr >= 2 && n_trt >= 2) {
      m_ctr <- mean(ctr_vals)
      m_trt <- mean(trt_vals)
      v_ctr <- var(ctr_vals)
      v_trt <- var(trt_vals)
      
      res$mean_group_ctr <- m_ctr
      res$mean_group_trt <- m_trt
      res$diff_length    <- m_trt - m_ctr
      res$fold_change    <- if(m_ctr != 0) m_trt / m_ctr else NA_real_
      
      shapiro_p <- function(x) {
        if(length(x) < 3) return(0)
        if(uniqueN(x) == 1) return(0)
        if(length(x) <= 5000) shapiro.test(x)$p.value else nortest::lillie.test(x)$p.value
      }
      
      norm_ctr <- shapiro_p(ctr_vals)
      norm_trt <- shapiro_p(trt_vals)
      all_normal <- (norm_ctr > 0.05) && (norm_trt > 0.05)
      
      lev_p <- tryCatch({
        tmp_val <- c(ctr_vals, trt_vals)
        tmp_grp <- factor(c(rep("C", n_ctr), rep("T", n_trt)))
        car::leveneTest(tmp_val ~ tmp_grp, center = median)[["Pr(>F)"]][1]
      }, error = function(e) NA_real_)
      
      homo_var <- !is.na(lev_p) && lev_p > 0.05
      
      if (all_normal) {
        t_res <- t.test(ctr_vals, trt_vals, var.equal = homo_var)
        res$p_value <- t_res$p.value
        res$test_performed <- if(homo_var) "Student's t-test" else "Welch's t-test"
      } else {
        res$p_value <- wilcox.test(ctr_vals, trt_vals)$p.value
        res$test_performed <- "Wilcoxon rank-sum test"
      }
      
      pooled_sd <- sqrt(((n_ctr - 1) * v_ctr + (n_trt - 1) * v_trt) / (n_ctr + n_trt - 2))
      if (!is.na(pooled_sd) && pooled_sd > 0) {
        d <- abs(res$diff_length / pooled_sd)
        res$cohen_d <- d
        res$cohen_effect <- if(d < 0.2) "small" else if(d < 0.5) "medium" else "large"
      }
    }
    
    res
  }, by = unit_tmp]
  
  setnames(results, "unit_tmp", which_level)
  results[, padj := p.adjust(p_value, method = padj_method)]
  results[, Log2FC := log2(fold_change)]
  
  message(sprintf("Processing completed in %.2f seconds.", 
                  as.numeric(difftime(Sys.time(), start_time, units = "secs"))))
  
  return(as.data.frame(results))
}
