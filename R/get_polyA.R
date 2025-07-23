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
