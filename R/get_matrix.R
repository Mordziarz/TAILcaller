#' Creating a ma plot
#'
#' @param count_molecules_out the table was created using the count_molecules_out function
#' @param grouping_factor the name of the column in the table that divides the experiment into groups
#' @param which_level the name of the column by which the statistics should be grouped, either transcripts or genes
#' @param statistic the column from table which was created using the count_molecules function
#' @return A table object.
#' @export
#'

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
    dplyr::select(!!sym(grouping_factor), !!rlang::sym(which_level), !!rlang::sym(statistic)) %>%
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
