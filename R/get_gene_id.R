#' Merging a polyA table with gene_id from gtf file
#'
#' @param get_polyA_out the table was created using the get_polyA_out function
#' @param gtf_file the gtf file
#' @param transcript_column_gtf the name of the column that contains the transcript_id
#' @param gene_column_gtf the name of the column that contains the gene_id
#' @return A table object.
#' @export
#'
get_gene_id <- function(polyA_table=get_polyA_out,gtf_file= gtf, transcript_column_gtf="column_name",gene_column_gtf="column_name") {
  
  if (!base::is.data.frame(polyA_table)) {
    stop("Argument 'polyA_table' must be a data frame.")
  }
  
  if (!base::is.character(transcript_column_gtf) || !(transcript_column_gtf %in% colnames(gtf_file))) {
    stop(paste("Argument 'transcript_column_gtf' must be a column name in the gtf data frame. Please check if the column", transcript_column_gtf, "exists."))
  }
  
    if (!base::is.character(gene_column_gtf) || !(gene_column_gtf %in% colnames(gtf_file))) {
    stop(paste("Argument 'gene_column_gtf' must be a column name in the gtf data frame. Please check if the column", gene_column_gtf, "exists."))
  }

  base::message("Starting to process the data...")
  
  start_time <- base::Sys.time()

  gtf <- base::as.data.frame(gtf_file)

  gtf <- gtf[!base::duplicated(gtf[transcript_column_gtf]),]
  
  gtf <- gtf[,c(gene_column_gtf,transcript_column_gtf)]

  get_gene_id_out <- base::merge(polyA_table,gtf,by.x="transcript_id",by.y=transcript_column_gtf)
  
  end_time <- base::Sys.time()

  base::message("Execution time: ", base::round(difftime(end_time, start_time, units = "mins"),2), " minutes")
  
  return(get_gene_id_out)
}
