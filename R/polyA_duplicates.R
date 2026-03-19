#' Process polyA tail length data to identify and optionally remove duplicates.
#'
#' This function takes a data frame of polyA tail length measurements, identifies
#' duplicate `read_id` entries, and provides summaries and visualizations.
#' Optionally, it can remove duplicate rows based on a combined index of key columns.
#'
#' @param polyA_table A data frame containing polyA tail length data.
#'   It must include the columns: "read_id", "transcript_id", "polyA_length",
#'   "sample_name", and "group".
#' @param delete_duplicates A logical value. If `TRUE` (default), duplicate rows
#'   (defined by a combined index of `read_id`, `transcript_id`, `polyA_length`,
#'   `sample_name`, and `group`) will be removed from the `polyA_table`.
#'   If `FALSE`, the original `polyA_table` will be returned in the `processed_table`
#'   element, along with the duplicate summaries.
#' @param gene_column_gtf A character string specifying the name of the column 
#'   in `gtf_file` that contains gene identifiers. This column will be added 
#'   to the merged output.
#'
#' @author Mateusz Mazdziarz
#'
#' @import data.table
#' @import ggplot2

polyA_duplicates <- function(polyA_table, delete_duplicates = TRUE, gene_column_gtf = "gene_id") {
  
  if (missing(polyA_table)) stop("'polyA_table' must be defined.")
  
  dt <- as.data.table(polyA_table)
  
  if (!(gene_column_gtf %in% colnames(dt))) {
    stop(paste("Column", gene_column_gtf, "not found!"))
  }
  
  message("Starting processing for polyA duplicates...")
  

  summary_read_id <- dt[, .(count = .N), by = read_id][order(-count)]
  duplicated_reads_summary <- summary_read_id[count > 1]
  
  message(paste0("Found ", nrow(duplicated_reads_summary), " read_ids with duplicates."))
  
  read_categories <- summary_read_id[, .(
    num_reads = .N
  ), by = .(category = ifelse(count == 1, "Unique Single", "Duplicated or Multiple"))]
  
  read_categories[, percentage := num_reads / sum(num_reads) * 100]
  
  plot_categories <- ggplot(read_categories, aes(x = category, y = num_reads, fill = category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(num_reads, " (", round(percentage, 1), "%)")),
              vjust = -0.5) +
    labs(x = "Read category", y = "Number of reads") +
    theme_bw() + theme(legend.position = "none")
  
  if (delete_duplicates) {
    message("Deleting duplicates using internal indexing...")
    rows_before <- nrow(dt)
    
    cols_to_check <- c("read_id", gene_column_gtf, "polyA_length", "sample_name", "group")
    
    dt <- unique(dt, by = cols_to_check)
    
    rows_after <- nrow(dt)
    message(paste0("Removed ", rows_before - rows_after, " duplicate rows."))
  }
  
  return(list(
    processed_table = as.data.frame(dt),
    read_id_counts = as.data.frame(summary_read_id), 
    duplicated_read_ids = as.data.frame(duplicated_reads_summary),
    read_categories_summary = as.data.frame(read_categories),
    categories_plot = plot_categories
  ))
}