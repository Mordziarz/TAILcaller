#' Merge PolyA Table with Gene IDs from GTF File
#'
#' This function merges a polyA table (created by the get_polyA function) with 
#' gene identifiers extracted from a GTF (Gene Transfer Format) file. The merge 
#' is performed by matching transcript IDs between the two datasets, allowing 
#' for gene-level analysis of polyA tail length data.
#'
#' @param polyA_table A data frame containing polyA information, typically 
#'   the output from the `get_polyA()` function. Must contain a column named 
#'   "transcript_id" for merging purposes.
#' @param gtf_file A data frame containing GTF file information that has been 
#'   read into R. This should be a parsed GTF file with transcript and gene 
#'   identifier columns.
#' @param transcript_column_gtf A character string specifying the name of the 
#'   column in `gtf_file` that contains transcript identifiers. This column 
#'   will be used to match with the "transcript_id" column in `polyA_table`.
#' @param gene_column_gtf A character string specifying the name of the column 
#'   in `gtf_file` that contains gene identifiers. This column will be added 
#'   to the merged output.
#'
#' @author Mateusz Mazdziarz
#' 
#' @importFrom base is.data.frame is.character message Sys.time as.data.frame 

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
