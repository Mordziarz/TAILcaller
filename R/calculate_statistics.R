#' Calculation of statistics and creation of the table
#'
#' @param polyA_table the table created using the get_polyA function.
#' @param grouping_factor the name of the column in the table that divides the experiment into groups.
#' @param which_level the name of the column by which the statistics should be grouped, either transcripts or genes.
#' @param control_group the name of the control group in "".
#' @param treated_group the name of the treated group in "".
#' @return a table object.
#' @export
#'

calculate_statistics <- function(
    polyA_table = get_gene_id_out,
    grouping_factor = "group",
    which_level = "gene_id",
    control_group = NULL,
    treated_group = NULL
) {
  if (missing(polyA_table)) stop("'polyA_table' must be defined.")
  if (missing(control_group) | missing(treated_group)) {
    stop("Both 'control_group' and 'treated_group' must be provided.")
  }
  
  message("Starting processing...")
  start_time <- Sys.time()
  
  how_molecules <- polyA_table[!duplicated(polyA_table[[which_level]]), ]
  how_molecules <- how_molecules[, which_level, drop = FALSE]
  
  results <- lapply(1:nrow(how_molecules), function(i) {
    molecule <- how_molecules[[which_level]][i]
    subset_data <- polyA_table[polyA_table[[which_level]] == molecule, ]
    
    n_ctr <- sum(subset_data[[grouping_factor]] == control_group)
    n_trt <- sum(subset_data[[grouping_factor]] == treated_group)
    
    if (n_ctr < 2 | n_trt < 2) return(list(
      p_value = NA,
      mean_ctr = NA,
      mean_trt = NA,
      fold_change = NA,
      cohen_d = NA,
      cohen_effect = NA,
      diff_length = NA
    ))
    
    ctr_data <- subset_data$polyA_length[subset_data[[grouping_factor]] == control_group]
    trt_data <- subset_data$polyA_length[subset_data[[grouping_factor]] == treated_group]
    
    p_val <- suppressWarnings(
      wilcox.test(ctr_data, trt_data)$p.value
    )
    
    mean_ctr <- mean(ctr_data, na.rm = TRUE)
    mean_trt <- mean(trt_data, na.rm = TRUE)
    
    fc <- ifelse(mean_ctr == 0, NA, mean_trt / mean_ctr)
    
    sd_ctr <- sd(ctr_data, na.rm = TRUE)
    sd_trt <- sd(trt_data, na.rm = TRUE)
    pooled_sd <- sqrt(((n_ctr-1)*sd_ctr^2 + (n_trt-1)*sd_trt^2)/(n_ctr + n_trt - 2))
    cd <- abs((mean_ctr - mean_trt)/pooled_sd)
    
    list(
      p_value = p_val,
      mean_ctr = mean_ctr,
      mean_trt = mean_trt,
      fold_change = fc,
      cohen_d = cd,
      cohen_effect = ifelse(cd < 0.2, "small", 
                            ifelse(cd < 0.5, "medium", "large")),
      diff_length = mean_trt - mean_ctr
    )
  })
  
  how_molecules$p_value <- sapply(results, `[[`, "p_value")
  how_molecules$mean_group_ctr <- sapply(results, `[[`, "mean_ctr")
  how_molecules$mean_group_trt <- sapply(results, `[[`, "mean_trt")
  how_molecules$diff_length <- sapply(results, `[[`, "diff_length")
  how_molecules$fold_change <- sapply(results, `[[`, "fold_change")
  how_molecules$cohen_d <- sapply(results, `[[`, "cohen_d")
  how_molecules$cohen_effect <- sapply(results, `[[`, "cohen_effect")
  how_molecules$padj <- p.adjust(how_molecules$p_value, method = "fdr")
  how_molecules$Log2FC <- log2(how_molecules$fold_change)
  
  message("Processing completed. Time: ", 
          round(difftime(Sys.time(), start_time, units = "mins"), 1), 
          " minutes")
  
  return(how_molecules)
}
