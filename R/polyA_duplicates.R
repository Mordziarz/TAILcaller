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
#' @importFrom dplyr group_by summarise n arrange filter mutate select
#' @importFrom ggplot2 ggplot aes geom_bar geom_text labs theme_bw theme element_text

polyA_duplicates <- function(polyA_table, delete_duplicates = TRUE,gene_column_gtf="gene_id") {
  
  if (missing(polyA_table)) {
    stop("'polyA_table' must be defined.")
  }

      if (!base::is.character(gene_column_gtf) || !(gene_column_gtf %in% colnames(polyA_table))) {
    stop(paste("Argument 'gene_column_gtf' must be a column name in the gtf data frame. Please check if the column", gene_column_gtf, "exists."))
  }

  
  message("Starting processing for polyA duplicates...")
  
  summary_read_id <- polyA_table %>% 
    dplyr::group_by(read_id) %>% 
    dplyr::summarise(count = n(), .groups = 'drop') %>%
    dplyr::arrange(desc(count))
  
  duplicated_reads_summary <- summary_read_id %>%
    dplyr::filter(count > 1)
  
  message(paste0("Found ", nrow(duplicated_reads_summary), " read_ids with duplicates."))
  
  read_categories <- summary_read_id %>%
    dplyr::mutate(category = ifelse(count == 1, "Unique Single", "Duplicated or Multiple")) %>%
    dplyr::group_by(category) %>%
    dplyr::summarise(num_reads = n(), .groups = 'drop') %>%
    dplyr::mutate(percentage = num_reads / sum(num_reads) * 100)
  
  plot_categories <- ggplot(read_categories, aes(x = category, y = num_reads, fill = category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(num_reads, " (", round(percentage, 1), "%)")),
              vjust = -0.5, color = "black") +
    labs(x = "Read category", y = "Number of reads") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),legend.position = "none")
  
  
  output_table <- polyA_table
  
  output_table$gene_id2 <- output_table[[gene_column_gtf]] 

  if (delete_duplicates == TRUE) {
    message("Deleting duplicates based on combined index...")
    output_table$index <- paste0(output_table$read_id, "_",
                                 output_table$gene_id2, "_",
                                 output_table$polyA_length, "_",
                                 output_table$sample_name, "_",
                                 output_table$group)
    
    rows_before_deduplication <- nrow(output_table)
    output_table <- output_table[!duplicated(output_table$index),]
    rows_after_deduplication <- nrow(output_table)
    
    message(paste0("Removed ", rows_before_deduplication - rows_after_deduplication, " duplicate rows."))
    
    output_table <- dplyr::select(output_table, -index, -gene_id2)
  }
  
  return(list(
    processed_table = output_table,
    read_id_counts = summary_read_id, 
    duplicated_read_ids = duplicated_reads_summary,
    read_categories_summary = read_categories,
    categories_plot = plot_categories
  ))
}