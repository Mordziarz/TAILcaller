#' Counting molecules, calculating the mean and median for the groups
#'
#' @param polyA_table the table created using the get_polyA function.
#' @param grouping_factor the name of the column in the table that divides the experiment into groups.
#' @param which_level the name of the column by which the statistics should be grouped, either transcripts or genes.
#' @return a table object.
#' @export
#'

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
