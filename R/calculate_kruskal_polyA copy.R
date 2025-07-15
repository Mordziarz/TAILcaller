#' Calculation of statistics and creation of the table
#'
#' @param polyA_table the table created using the get_polyA function.
#' @param grouping_factor the name of the column in the table that divides the experiment into groups.
#' @param which_level the name of the column by which the statistics should be grouped, either transcripts or genes.
#' @return a table object.
#' @export
#'

calculate_polyA_stat_n3 <- function(polyA_table    = get_gene_id_out, grouping_factor= "group", which_level    = "gene_id",padj_method    = "fdr") {
  if (missing(polyA_table)) {
    stop("'polyA_table' must be defined.")
  }
  required_cols <- c("polyA_length", grouping_factor, which_level)
  if (!all(required_cols %in% colnames(polyA_table))) {
    stop(sprintf(
      "polyA_table must contain columns: %s",
      paste(required_cols, collapse = ", ")
    ))
  }
  
  message("Starting calculations...")
  start_time <- Sys.time()
  
  units <- unique(polyA_table[[which_level]])
  n_units <- length(units)

  out_df <- data.frame(
    !!which_level := units,
    test_used     = character(n_units),
    p_value       = numeric(n_units),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(units)) {
    unit_id <- units[i]
    subdf <- polyA_table[polyA_table[[which_level]] == unit_id, , drop=FALSE]
    n_obs <- nrow(subdf)
    groups_present <- unique(subdf[[grouping_factor]])
    n_groups <- length(groups_present)
    
    p_val <- NA
    chosen_test <- NA_character_
    
    if (n_groups >= 2) {
      if (n_obs <= 5000) {
        p_norm <- stats::shapiro.test(subdf$polyA_length)$p.value
      } else {
        p_norm <- nortest::lillie.test(subdf$polyA_length)$p.value
      }
      lev_p <- tryCatch({
        car::leveneTest(
          as.formula(paste0("polyA_length ~ ", grouping_factor)),
          data = subdf
        )[["Pr(>F)"]][1]
      }, error = function(e) NA)
      
      if (n_groups > 2 && !is.na(p_norm) && p_norm > .05 && !is.na(lev_p) && lev_p > .05) {
        aov_res <- stats::aov(
          as.formula(paste0("polyA_length ~ ", grouping_factor)),
          data = subdf
        )
        p_val <- summary(aov_res)[[1]][["Pr(>F)"]][1]
        chosen_test <- "ANOVA"
      } else {
        p_val <- suppressWarnings(
          stats::kruskal.test(
            as.formula(paste0("polyA_length ~ ", grouping_factor)),
            data = subdf
          )$p.value
        )
        chosen_test <- "Kruskal–Wallis"
      }
    }
    
    out_df$test_used[i] <- chosen_test
    out_df$p_value[i]   <- p_val
  }
  
  out_df$padj <- stats::p.adjust(out_df$p_value, method = padj_method)
  
  end_time <- Sys.time()
  message(
    "Done in ",
    round(difftime(end_time, start_time, units = "mins"), 2),
    " mins."
  )
  
  return(out_df)
}
