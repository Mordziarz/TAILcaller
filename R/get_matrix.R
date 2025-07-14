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
#' @return A numeric matrix where:
#'   \describe{
#'     \item{Rows}{Represent individual transcripts or genes (from `which_level`)}
#'     \item{Columns}{Represent different experimental groups (from `grouping_factor`)}
#'     \item{Values}{Contain the specified statistic for each transcript/gene-group combination}
#'     \item{Missing values}{Are replaced with 0 for computational convenience}
#'   }
#'   Row names correspond to transcript or gene identifiers, and column names 
#'   correspond to group identifiers.
#'
#' @export
#'
#' @details
#' The function performs the following data transformation steps:
#' \itemize{
#'   \item Validates that all required columns exist in the input data
#'   \item Selects only the necessary columns (grouping factor, level, and statistic)
#'   \item Pivots the data from long to wide format using `tidyr::pivot_wider()`
#'   \item Converts the result to a data frame and then to a matrix
#'   \item Sets row names to the specified level identifiers
#'   \item Replaces NA values with 0 for missing group-level combinations
#' }
#' 
#' The transformation uses `tidyr::pivot_wider()` with non-standard evaluation 
#' through `rlang::sym()` to handle column names passed as character strings. 
#' This approach maintains flexibility while ensuring proper data reshaping.
#'
#' @section Input Validation:
#' The function performs comprehensive validation checks and will stop 
#' execution with informative error messages if:
#' \itemize{
#'   \item `count_molecules_out` parameter is missing or undefined
#'   \item `grouping_factor` column does not exist in the input data
#'   \item `which_level` column does not exist in the input data
#'   \item `statistic` column does not exist in the input data
#' }
#'
#' @section Data Requirements:
#' The input data frame must contain:
#' \itemize{
#'   \item Column specified by `grouping_factor` with group identifiers
#'   \item Column specified by `which_level` with transcript/gene identifiers
#'   \item Column specified by `statistic` with numeric statistical values
#'   \item No duplicate combinations of grouping factor and level values
#' }
#'
#' @section Matrix Properties:
#' The resulting matrix has the following characteristics:
#' \itemize{
#'   \item **Dimensions**: Rows = unique levels, Columns = unique groups
#'   \item **Data Type**: Numeric matrix suitable for mathematical operations
#'   \item **Missing Values**: Converted to 0 for computational stability
#'   \item **Row Names**: Set to transcript/gene identifiers for easy reference
#' }
#'
#' @section Performance:
#' The function reports processing status including:
#' \itemize{
#'   \item Start notification when processing begins
#'   \item Processing completion time in minutes
#'   \item Success confirmation message
#' }
#'
#' @seealso 
#' \code{\link{count_molecules}} for creating the input statistics table
#' \code{\link[tidyr]{pivot_wider}} for data reshaping functionality
#' \code{\link[base]{matrix}} for matrix creation and manipulation
#' \code{\link[base]{as.matrix}} for data frame to matrix conversion
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
