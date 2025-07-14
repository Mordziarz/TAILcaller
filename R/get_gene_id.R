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
#' @return A data frame containing the merged polyA table with additional gene 
#'   information. The returned table includes all columns from the original 
#'   `polyA_table` plus the gene identifier column from the GTF file. Only 
#'   rows with matching transcript IDs between the two datasets will be retained.
#'
#' @export
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Validates input parameters and column existence
#'   \item Removes duplicate transcript entries from the GTF data
#'   \item Extracts only the necessary columns (gene ID and transcript ID) from GTF
#'   \item Merges the datasets using transcript IDs as the matching key
#'   \item Reports processing time upon completion
#' }
#' 
#' The merge operation is performed using `base::merge()` with the transcript_id 
#' column from the polyA table matched against the specified transcript column 
#' from the GTF file. This is an inner join, meaning only records with matching 
#' transcript IDs in both datasets will be included in the result.
#'
#' @section Input Requirements:
#' \itemize{
#'   \item `polyA_table` must be a data frame with a "transcript_id" column
#'   \item `gtf_file` must be a data frame with the specified transcript and gene columns
#'   \item Column names specified in `transcript_column_gtf` and `gene_column_gtf` 
#'     must exist in the GTF data frame
#' }
#'
#' @section Error Handling:
#' The function performs validation checks and will stop execution with 
#' informative error messages if:
#' \itemize{
#'   \item `polyA_table` is not a data frame
#'   \item `transcript_column_gtf` is not a character string or doesn't exist in GTF
#'   \item `gene_column_gtf` is not a character string or doesn't exist in GTF
#' }
#'
#' @seealso 
#' \code{\link{get_polyA}} for creating the polyA table input
#' \code{\link[base]{merge}} for the underlying merge functionality
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
