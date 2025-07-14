#' Count Molecules and Calculate Statistics for Grouped Data
#'
#' This function performs statistical analysis on polyA tail length data by 
#' grouping the data according to specified factors and calculating summary 
#' statistics. It counts the number of molecules per group and computes the 
#' mean and median polyA tail lengths, providing comprehensive insights into 
#' the distribution of polyA tails across different experimental conditions.
#'
#' @param polyA_table A data frame containing polyA tail length information, 
#'   typically the output from the `get_gene_id()` function. The table must 
#'   contain columns for grouping factors, identification levels, and a 
#'   "polyA_length" column with numeric values.
#' @param grouping_factor A character string specifying the name of the column 
#'   in `polyA_table` that contains the grouping variable used to divide the 
#'   experiment into different groups or conditions. This column will be used 
#'   as the primary grouping factor for statistical calculations.
#' @param which_level A character string specifying the name of the column 
#'   that defines the level at which statistics should be calculated. This 
#'   parameter determines whether the analysis is performed at the transcript 
#'   level (e.g., "transcript_id") or gene level (e.g., "gene_id"), allowing 
#'   for flexible aggregation of molecular data.
#'
#' @return A grouped data frame (tibble) with summary statistics containing 
#'   the following columns:
#'   \describe{
#'     \item{grouping_factor}{The values from the specified grouping column}
#'     \item{which_level}{The values from the specified level column}
#'     \item{count}{Integer count of molecules/observations in each group}
#'     \item{avg_polyA_length}{Numeric mean of polyA tail lengths for each group}
#'     \item{median_polyA_length}{Numeric median of polyA tail lengths for each group}
#'   }
#'   Each row represents a unique combination of grouping factor and level values.
#'
#' @export
#'
#' @details
#' The function performs the following statistical operations:
#' \itemize{
#'   \item Groups the data by the specified grouping factor and level columns
#'   \item Counts the number of observations (molecules) in each group using \code{n()}
#'   \item Calculates the arithmetic mean of polyA tail lengths using \code{mean()}
#'   \item Calculates the median of polyA tail lengths using \code{median()}
#'   \item Reports processing time and completion status
#' }
#' 
#' The function uses `dplyr` operations with non-standard evaluation (NSE) 
#' through the `sym()` function to handle column names passed as character 
#' strings. This allows for flexible column specification while maintaining 
#' the functionality of `dplyr` grouping and summarization operations.
#'
#' @section Performance:
#' The function displays progress messages during execution, including:
#' \itemize{
#'   \item Start notification when processing begins
#'   \item Processing completion time in minutes
#'   \item Success confirmation message
#' }
#' 
#' @section Statistical Measures:
#' \itemize{
#'   \item **Count**: Total number of molecules/reads in each group
#'   \item **Mean**: Average polyA tail length, sensitive to outliers
#'   \item **Median**: Middle value of polyA tail lengths, robust to outliers
#' }
#' 
#' @section Data Requirements:
#' The input `polyA_table` must contain:
#' \itemize{
#'   \item A column named "polyA_length" with numeric values
#'   \item The column specified in `grouping_factor` parameter
#'   \item The column specified in `which_level` parameter
#' }
#'
#' @seealso 
#' \code{\link{get_gene_id}} for creating the input polyA table
#' \code{\link[dplyr]{group_by}} for data grouping functionality
#' \code{\link[dplyr]{summarise}} for data summarization operations
#' \code{\link[base]{mean}} and \code{\link[stats]{median}} for statistical calculations
#'
#' @author Mateusz Mazdziarz
#'
#' @importFrom dplyr group_by summarise n
#' @importFrom rlang sym
#' @importFrom magrittr %>%
#' @importFrom base message Sys.time round difftime mean
#' @importFrom stats median

count_molecules <- function(polyA_table=get_gene_id_out,grouping_factor="group",which_level="transcript_id"){

  base::message("Starting to process the data and calculate transcript statistics...")
  
  start_time <- base::Sys.time()
  
  count_table <- polyA_table %>%

    dplyr::group_by(
      !!sym(grouping_factor), 
      !!sym(which_level)
      
      ) %>%

    dplyr::summarise(
  
      count = dplyr::n(),  
  
      avg_polyA_length = base::mean(polyA_length),  
  
      median_polyA_length = stats::median(polyA_length)
  )
  
  end_time <- base::Sys.time()

  base::message("Processing complete. Time taken: ", base::round(difftime(end_time, start_time, units = "mins"), 2), " minutes")
  
  base::message("Transcript statistics have been calculated successfully.")
  
return(count_table)

}
