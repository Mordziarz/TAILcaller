#' Create a Wide-Format Matrix from Grouped Molecular Statistics
#'
#' This function transforms a long-format table containing grouped molecular 
#' statistics into a wide-format matrix suitable for downstream analysis, 
#' visualization, or statistical modeling. The function reshapes data by 
#' pivoting specified statistics across different experimental groups, with 
#' transcript or gene identifiers as row names and group conditions as columns.
#'
#' @param count_molecules_out A data frame containing grouped molecular 
#'   statistics, typically the output from the `count_molecules()` function. 
#'   The table must contain columns for grouping factors, identification 
#'   levels, and statistical measures (count, avg_polyA_length, 
#'   median_polyA_length).
#' @param grouping_factor A character string specifying the name of the column 
#'   in `count_molecules_out` that contains the grouping variable. Values from 
#'   this column will become the column names in the resulting matrix, 
#'   representing different experimental conditions or groups.
#' @param which_level A character string specifying the name of the column 
#'   that contains identifiers (transcript_id or gene_id). Values from this 
#'   column will become the row names in the resulting matrix, defining the 
#'   molecular entities being analyzed.
#' @param statistic A character string specifying which statistical measure 
#'   to use for matrix values. Must be one of the columns from the 
#'   `count_molecules_out` table, typically "count", "avg_polyA_length", 
#'   or "median_polyA_length".
#'
#' @author Mateusz Mazdziarz
#'
#' @importFrom dplyr select
#' @importFrom rlang sym
#' @importFrom magrittr %>%
#' @importFrom tidyr pivot_wider
#' @importFrom base message Sys.time missing colnames paste stop as.data.frame rownames as.matrix round difftime list

get_matrix <- function(count_molecules_out = count_molecules_out,
                       grouping_factor = "group",
                       which_level = "transcript_id",
                       statistic = "avg_polyA_length") {

  base::message("Starting to process the data and calculate transcript statistics...")

  start_time <- base::Sys.time()

  if (base::missing(count_molecules_out)) {
    stop("'count_molecules_out' must be defined.")
  }

  if (!(grouping_factor %in% base::colnames(count_molecules_out))) {
    stop(base::paste("Column", grouping_factor, "not found in the data frame."))
  }

  if (!(which_level %in% base::colnames(count_molecules_out))) {
    stop(base::paste("Column", which_level, "not found in the data frame."))
  }

  if (!(statistic %in% base::colnames(count_molecules_out))) {
    stop(base::paste("Column", statistic, "not found in the data frame."))
  }

  count_table_wide <- count_molecules_out %>%
    dplyr::select(!!rlang::sym(grouping_factor), !!rlang::sym(which_level), !!rlang::sym(statistic)) %>%
    tidyr::pivot_wider(
      names_from = !!rlang::sym(grouping_factor),
      values_from = !!rlang::sym(statistic),
      values_fill = base::list(statistic = NA)
    )

  count_table_wide <- base::as.data.frame(count_table_wide)

  base::rownames(count_table_wide) <- count_table_wide[[which_level]]

  count_table_wide[[which_level]] <- NULL

  count_table_wide <- base::as.matrix(count_table_wide)

  count_table_wide[is.na(count_table_wide)] <- 0

  end_time <- base::Sys.time()

  base::message("Processing complete. Time taken: ", base::round(difftime(end_time, start_time, units = "mins"), 2), " minutes")

  base::message("Transcript statistics have been calculated successfully.")

  return(count_table_wide)
}
