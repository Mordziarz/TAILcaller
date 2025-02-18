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

calculate_statistics <- function(polyA_table=get_gene_id_out, grouping_factor="group", which_level="gene_id",
                                 control_group=NULL, treated_group=NULL) {
  
  if (missing(polyA_table)) {
    stop("'polyA_table' must be defined.")
  }
  
  if (missing(control_group) | missing(treated_group)) {
    stop("Both 'control_group' and 'treated_group' must be provided.")
  }
  
  base::message("Starting to process the data and calculate statistics...")
  
  start_time <- base::Sys.time()
  
  how_molecules <- polyA_table[!duplicated(polyA_table[[which_level]]), ]
  how_molecules <- how_molecules[, which_level, drop=FALSE]
  
  p_values <- base::numeric(base::nrow(how_molecules))
  mean_group_ctr <- base::numeric(base::nrow(how_molecules))
  mean_group_trt <- base::numeric(base::nrow(how_molecules))
  fold_change <- base::numeric(base::nrow(how_molecules))
  cohen_d <- base::numeric(base::nrow(how_molecules))
  cohen_effect <- base::character(base::nrow(how_molecules))
  diff_length <- base::numeric(base::nrow(how_molecules))
  
  formula <- base::paste0("polyA_length ~ ", grouping_factor)
  
  for (i in 1:base::nrow(how_molecules)) {
    
    which_molecule <- how_molecules[[which_level]][i]
    
    polyA_table_unique <- polyA_table[polyA_table[[which_level]] == which_molecule, ]
    
    if (base::length(base::unique(polyA_table_unique$group)) == 2) {
      
      group_ctr_data <- polyA_table_unique$polyA_length[polyA_table_unique$group == control_group]
      group_trt_data <- polyA_table_unique$polyA_length[polyA_table_unique$group == treated_group]
      
      p_val <- base::suppressWarnings(stats::wilcox.test(stats::as.formula(formula), data = polyA_table_unique)$p.value)
      p_values[i] <- p_val
      
      mean_group_ctr[i] <- base::mean(group_ctr_data, na.rm = TRUE)
      mean_group_trt[i] <- base::mean(group_trt_data, na.rm = TRUE)
      
      fold_change[i] <- mean_group_trt[i] / mean_group_ctr[i]
      
      mean_ctr <- base::mean(group_ctr_data, na.rm = TRUE)
      mean_trt <- base::mean(group_trt_data, na.rm = TRUE)
      sd_ctr <- stats::sd(group_ctr_data, na.rm = TRUE)
      sd_trt <- stats::sd(group_trt_data, na.rm = TRUE)
      n_ctr <- base::length(group_ctr_data)
      n_trt <- base::length(group_trt_data)
      
      pooled_sd <- base::sqrt(((n_ctr - 1) * sd_ctr^2 + (n_trt - 1) * sd_trt^2) / (n_ctr + n_trt - 2))
      cohen_d[i] <- base::abs((mean_ctr - mean_trt) / pooled_sd)
      
      cohen_effect[i] <- base::ifelse(cohen_d[i] < 0.2, "small effect",
                                base::ifelse(cohen_d[i] < 0.5, "medium effect", "large effect"))
      
      diff_length[i] <- mean_group_trt[i] - mean_group_ctr[i]
      
    } else {
      p_values[i] <- NA
      mean_group_ctr[i] <- NA
      mean_group_trt[i] <- NA
      fold_change[i] <- NA
      cohen_d[i] <- NA
      cohen_effect[i] <- NA
      diff_length[i] <- NA
    }
  }
  
  how_molecules$p_value <- p_values
  how_molecules$mean_group_ctr <- mean_group_ctr
  how_molecules$mean_group_trt <- mean_group_trt
  how_molecules$diff_length <- diff_length
  how_molecules$fold_change <- fold_change
  how_molecules$cohen_d <- cohen_d
  how_molecules$cohen_effect <- cohen_effect
  how_molecules$padj <- stats::p.adjust(how_molecules$p_value, method = "fdr")
  how_molecules$Log2FC <- log2(how_molecules$fold_change)
  
  end_time <- base::Sys.time()
  
  base::message("Processing complete. Time taken: ", base::round(difftime(end_time, start_time, units = "mins"), 2), " minutes")
  base::message("Statistics have been calculated successfully.")
  
  return(how_molecules)
}
