#' Creating a polyA table based on the bams from dorado
#'
#' @param samples_table the table was created using the get_RSCU function
#' @return A table object.
#' @export
#'

get_polyA <- function(samples_table=samples_table) {

  message("Starting to process the data...")

  if(!all(c("bam_path", "sample_name") %in% colnames(samples_table))) {
    stop("The samples_table must contain the columns 'bam_path', 'sample_name' and 'group'.")
  }

  start_time <- base::Sys.time()

  param <- Rsamtools::ScanBamParam(what = c("qname", "rname"),tag="pt")

  tabela_ogony_bam <- data_frame()

  for (i in 1:base::nrow(samples_table)) {

    bam_data <- Rsamtools::scanBam(samples_table$bam_path[i],param=param)

    tabela_bam <- base::as.data.frame(bam_data)

    tabela_bam$samples_name <- samples_table$sample_name[i]

    tabela_bam$group <- samples_table$group[i]

    tabela_ogony_bam <- base::rbind(tabela_ogony_bam,tabela_bam)

    if (i %% 1 == 0) {
      message(paste("Processed samples:", i, "out of", base::nrow(samples_table)))
    }
  }

  colnames(tabela_ogony_bam) <- c("read_id", "transcript_id", "polyA_length", "sample_name", "group")

  tabela_ogony_bam <- tabela_ogony_bam[!is.na(tabela_ogony_bam$polyA_length),]

  tabela_ogony_bam <- tabela_ogony_bam[!is.na(tabela_ogony_bam$transcript_id),]

  end_time <- base::Sys.time()

  base::message("Execution time: ", base::round(difftime(end_time, start_time, units = "mins"),2), " minutes")

  return(tabela_ogony_bam)

}
