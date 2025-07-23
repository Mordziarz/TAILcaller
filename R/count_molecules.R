#' Count Molecules and Calculate Statistics for Grouped Data
#'
#' This function performs statistical analysis on polyA tail length data by 
#' grouping the data according to specified factors and calculating summary 
#' statistics. It counts the number of molecules per group and computes the 
#' mean and median polyA tail lengths, providing comprehensive insights into 
#' the distribution of polyA tails across different experimental conditions.
#'
#' @param polyA_table A data frame containing polyA tail length information, 
#'   typically the output from the `get_gene_id()` function. The table must 
#'   contain columns for grouping factors, identification levels, and a 
#'   "polyA_length" column with numeric values.
#' @param grouping_factor A character string specifying the name of the column 
#'   in `polyA_table` that contains the grouping variable used to divide the 
#'   experiment into different groups or conditions. This column will be used 
#'   as the primary grouping factor for statistical calculations.
#' @param which_level A character string specifying the name of the column 
#'   that defines the level at which statistics should be calculated. This 
#'   parameter determines whether the analysis is performed at the transcript 
#'   level (e.g., "transcript_id") or gene level (e.g., "gene_id"), allowing 
#'   for flexible aggregation of molecular data.
#'
#' @author Mateusz Mazdziarz
#'
#' @importFrom dplyr group_by summarise n
#' @importFrom rlang sym
#' @importFrom magrittr %>%
#' @importFrom base message Sys.time round difftime mean
#' @importFrom stats median

count_molecules <- function(polyA_table=get_gene_id_out,grouping_factor="group",which_level="transcript_id"){

  base::message("Starting to process the data and calculate transcript statistics...")
  
  start_time <- base::Sys.time()
  
  count_table <- polyA_table %>%

    dplyr::group_by(
      !!sym(grouping_factor), 
      !!sym(which_level)
      
      ) %>%

    dplyr::summarise(
  
      count = dplyr::n(),  
  
      avg_polyA_length = base::mean(polyA_length),  
  
      median_polyA_length = stats::median(polyA_length)
  )
  
  end_time <- base::Sys.time()

  base::message("Processing complete. Time taken: ", base::round(difftime(end_time, start_time, units = "mins"), 2), " minutes")
  
  base::message("Transcript statistics have been calculated successfully.")
  
return(count_table)

}
