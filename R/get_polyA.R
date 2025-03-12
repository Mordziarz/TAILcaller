#' Creating a polyA table based on the bams from dorado
#'
#' @param samples_table the table with bam description
#' @return a table object.
#' @export
#'

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
