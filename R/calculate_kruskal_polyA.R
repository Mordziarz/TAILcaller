#' Calculate Kruskal–Wallis p-values per entity for polyA lengths
#'
#' \code{calculate_kruskal_polyA} computes, for each unique entity in a
#' polyA length table, the Kruskal–Wallis test p-value comparing
#' \code{polyA_length} distributions across groups. Only entities with
#' more than two groups are tested; others receive \code{NA}.
#'
#' @param polyA_table A \code{data.frame} containing at least the columns
#'   \code{"polyA_length"}, the grouping factor, and the entity identifier.
#'   Defaults to \code{get_gene_id_out}.
#' @param grouping_factor A string naming the column in \code{polyA_table}
#'   that defines the groups for comparison. Default \code{"group"}.
#' @param which_level A string naming the column in \code{polyA_table}
#'   that defines the entity (e.g., gene, transcript). Default
#'   \code{"gene_id"}.
#' @param padj_method A string specifying the p-value adjustment method
#'   passed to \code{\link[stats]{p.adjust}}. Default \code{"fdr"}.
#'
#' @author Mateusz Mazdziarz
#'
#' @importFrom stats kruskal.test p.adjust

calculate_kruskal_polyA <- function(polyA_table=get_gene_id_out, grouping_factor="group", which_level="gene_id",padj_method="fdr") {
  
  if (missing(polyA_table)) {
    stop("'polyA_table' must be defined.")
  }
  
  base::message("Starting to process the data and calculate statistics...")
  
  start_time <- base::Sys.time()
  
  how_molecules <- polyA_table[!duplicated(polyA_table[[which_level]]), ]
  how_molecules <- how_molecules[, which_level, drop=FALSE]
  
  p_values <- base::numeric(base::nrow(how_molecules))
  
  formula <- base::paste0("polyA_length ~ ", grouping_factor)
  
  for (i in 1:base::nrow(how_molecules)) {
    
    which_molecule <- how_molecules[[which_level]][i]
    
    polyA_table_unique <- polyA_table[polyA_table[[which_level]] == which_molecule, ]
    
    if (base::nrow(base::unique(polyA_table_unique[grouping_factor])) > 2) {
      
      p_val <- base::suppressWarnings(stats::kruskal.test(stats::as.formula(formula), data = polyA_table_unique)$p.value)
      p_values[i] <- p_val
      
      
    } else {
      p_values[i] <- NA
    }
  }
  
  how_molecules$p_value <- p_values
  how_molecules$padj <- stats::p.adjust(how_molecules$p_value, method = padj_method)
  
  end_time <- base::Sys.time()
  
  base::message("Processing complete. Time taken: ", base::round(difftime(end_time, start_time, units = "mins"), 2), " minutes")
  base::message("Statistics have been calculated successfully.")
  
  return(how_molecules)
}