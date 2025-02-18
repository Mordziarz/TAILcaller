#' Creating a polyA table based on the bams from dorado
#'
#' @param get_polyA_out the table was created using the get_polyA_out function
#' @param transcript_column_gtf the name of the column that contains the transcript_id
#' @return A table object.
#' @export
#'
get_gene_id <- function(get_polyA_out=get_polyA_out,transcript_column_gtf="column_name") {
  
  if (!base::is.data.frame(get_polyA_out)) {
    stop("Argument 'get_polyA_out' must be a data frame.")
  }
  
  if (!base::is.character(transcript_column_gtf) || !(transcript_column_gtf %in% colnames(gtf_df))) {
    stop(paste("Argument 'transcript_column_gtf' must be a column name in the gtf data frame. Please check if the column", transcript_column_gtf, "exists."))
  }
  
  base::message("Starting to process the data...")
  
  start_time <- base::Sys.time()

  get_gene_id_out <- base::merge(get_polyA_out,gtf_df,by.x="transcript_id",by.y=transcript_column_gtf)
  
  end_time <- base::Sys.time()

  base::message("Execution time: ", base::round(difftime(end_time, start_time, units = "mins"),2), " minutes")
  
  return(get_gene_id_out)
}
