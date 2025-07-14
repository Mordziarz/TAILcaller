#' Create a PolyA Length Table from BAM Files
#'
#' This function processes BAM files containing nanopore sequencing data from 
#' Dorado to extract polyA tail length information. It iterates through 
#' multiple BAM files specified in a samples table and consolidates the 
#' polyA tail data into a single data frame.
#'
#' @param samples_table A data frame containing BAM file information. Must 
#'   include the following columns:
#'   \describe{
#'     \item{bam_path}{Character vector with full paths to BAM files}
#'     \item{sample_name}{Character vector with sample identifiers}
#'     \item{group}{Character vector with group classifications}
#'   }
#'
#' @return A data frame with the following columns:
#'   \describe{
#'     \item{read_id}{Character vector containing read identifiers}
#'     \item{transcript_id}{Character vector containing transcript identifiers}
#'     \item{polyA_length}{Numeric vector containing polyA tail lengths}
#'     \item{sample_name}{Character vector containing sample names}
#'     \item{group}{Character vector containing group classifications}
#'   }
#'   Only returns rows where both polyA_length and transcript_id are not NA.
#'
#' @export
#'
#' @details 
#' The function uses the Rsamtools package to scan BAM files and extract 
#' specific information including query names (qname), reference names (rname), 
#' and the "pt" tag which contains polyA tail length information from Dorado.
#' 
#' Progress messages are displayed during processing, showing the number of 
#' samples processed. The function also reports the total execution time upon 
#' completion.
#' 
#' @section Requirements:
#' - Input BAM files must be indexed
#' - BAM files must contain the "pt" tag (polyA tail length from Dorado)
#' - The samples_table must contain all required columns
#'
#' @section Error Handling:
#' The function stops execution if the required columns (bam_path, sample_name, 
#' group) are not present in the samples_table.
#'
#' @seealso 
#' \code{\link[Rsamtools]{scanBam}} for BAM file reading functionality
#' \code{\link[Rsamtools]{ScanBamParam}} for BAM scanning parameters
#'
#' @author Mateusz Mazdziarz
#'
#' @importFrom Rsamtools scanBam ScanBamParam
#' @importFrom base data.frame nrow rbind Sys.time message paste round difftime

get_polyA <- function(samples_table=samples_table) {

  message("Starting to process the data...")

  if(!all(c("bam_path", "sample_name") %in% colnames(samples_table))) {
    stop("The samples_table must contain the columns 'bam_path', 'sample_name' and 'group'.")
  }

  start_time <- base::Sys.time()

  param <- Rsamtools::ScanBamParam(what = c("qname", "rname"),tag="pt")

  table_ogony_bam <- base::data.frame()

  for (i in 1:base::nrow(samples_table)) {

    bam_data <- Rsamtools::scanBam(samples_table$bam_path[i],param=param)

    table_bam <- base::as.data.frame(bam_data)

    table_bam$samples_name <- samples_table$sample_name[i]

    table_bam$group <- samples_table$group[i]

    table_ogony_bam <- base::rbind(table_ogony_bam,table_bam)

    if (i %% 1 == 0) {
      message(paste("Processed samples:", i, "out of", base::nrow(samples_table)))
    }
  }

  colnames(table_ogony_bam) <- c("read_id", "transcript_id", "polyA_length", "sample_name", "group")

  table_ogony_bam <- table_ogony_bam[!is.na(table_ogony_bam$polyA_length),]

  table_ogony_bam <- table_ogony_bam[!is.na(table_ogony_bam$transcript_id),]

  end_time <- base::Sys.time()

  base::message("Execution time: ", base::round(difftime(end_time, start_time, units = "mins"),2), " minutes")

  return(table_ogony_bam)

}
