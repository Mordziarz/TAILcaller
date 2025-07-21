#' Calculate per-entity group comparison p-values for polyA lengths
#'
#' \code{calculate_polyA_stat_n3} performs, for each unique entity in a 
#' polyA length table, a hypothesis test comparing distributions of 
#' \code{polyA_length} across groups. For each entity it:
#'   1. Checks if there are at least two groups and ≥3 total observations.
#'   2. Tests normality (Shapiro–Wilk or Lilliefors).
#'   3. Tests variance homogeneity (Levene’s test) if each group has ≥2 obs.
#'   4. Chooses one of:
#'      - One-way ANOVA (normal + homogeneous variances);
#'      - Welch ANOVA (normal + unequal variances);
#'      - Kruskal–Wallis (non-normal or insufficient obs).
#'   5. Records the p-value and test used.
#'   6. Adjusts all p-values for multiple testing.
#'
#' @param polyA_table A \code{data.frame} containing at least the columns 
#'   \code{"polyA_length"}, the grouping factor, and the entity identifier.
#'   Defaults to \code{get_gene_id_out}.
#' @param grouping_factor A string naming the column in \code{polyA_table} 
#'   that defines groups for comparison. Default \code{"group"}.
#' @param which_level A string naming the column in \code{polyA_table} 
#'   that defines the entity (e.g., gene, transcript). Default \code{"gene_id"}.
#' @param padj_method A string specifying the p-value adjustment method 
#'   passed to \code{\link[stats]{p.adjust}}. Default \code{"fdr"}.
#'
#' @return A \code{data.frame} with one row per entity and columns:
#'   \describe{
#'     \item{<which_level>}{Entity identifier (named as supplied).}
#'     \item{test_used}{The test applied: \code{"ANOVA"}, \code{"Welch ANOVA"}, 
#'       or \code{"Kruskal–Wallis"}.}
#'     \item{p_value}{Raw p-value from the chosen test.}
#'     \item{padj}{Adjusted p-value across all entities.}
#'   }
#'
#' @details
#' - Entities with fewer than two groups yield \code{NA} p-values.
#' - Entities with <3 total observations default to Kruskal–Wallis.
#' - Normality assessed by Shapiro–Wilk for n ≤ 5000, otherwise by Lilliefors.
#' - Variance homogeneity assessed by Levene’s test when each group has ≥2 obs.
#' - ANOVA used if normal + homogeneous; Welch ANOVA if normal + heterogeneous; 
#'   otherwise Kruskal–Wallis.
#'
#' @author Mateusz Mazdziarz
#'
#' @importFrom stats shapiro.test kruskal.test aov oneway.test p.adjust
#' @importFrom nortest lillie.test
#' @importFrom car leveneTest
#' @export

calculate_polyA_stat_n3 <- function(polyA_table = get_gene_id_out,grouping_factor = "group",which_level = "gene_id",padj_method = "fdr") {
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

  units   <- unique(polyA_table[[which_level]])
  n_units <- length(units)

  out_df <- data.frame(
    test_used = character(n_units),
    p_value   = numeric(n_units),
    stringsAsFactors = FALSE
  )
  out_df[[which_level]] <- units

  for (i in seq_along(units)) {
    unit_id <- units[i]
    subdf   <- polyA_table[polyA_table[[which_level]] == unit_id, , drop = FALSE]
    n_obs   <- nrow(subdf)
    grp_tbl <- table(subdf[[grouping_factor]])
    n_groups <- length(grp_tbl)

    p_val       <- NA_real_
    chosen_test <- NA_character_

    if (n_groups >= 2) {
      if (n_obs < 3 || any(grp_tbl < 2)) {
        p_val       <- suppressWarnings(
                         stats::kruskal.test(
                           as.formula(paste0("polyA_length ~ ", grouping_factor)),
                           data = subdf
                         )$p.value
                       )
        chosen_test <- "Kruskal–Wallis"
      } else {
        normality_p_values <- tapply(subdf$polyA_length, subdf[[grouping_factor]], function(x) {
          if (length(x) >= 3) {
            if (length(unique(x)) == 1) { # Added opening bracket
              return(0)
            } # Added closing bracket for 'if (length(unique(x)) == 1)'
            tryCatch({
              if (length(x) <= 5000) {
                stats::shapiro.test(x)$p.value
              } else {
                nortest::lillie.test(x)$p.value
              }
            }, error = function(e) {
              warning(paste("Error in normality test for group:", e$message))
              return(0)
            })
          } else {
            0
          }
        })
        all_groups_normal <- all(normality_p_values > 0.05, na.rm = TRUE)

        if (all(grp_tbl >= 2)) {
          lev_p <- tryCatch({
            car::leveneTest(
              as.formula(paste0("polyA_length ~ ", grouping_factor)),
              data = subdf
            )[["Pr(>F)"]][1]
          }, error = function(e) NA_real_)
        } else {
          lev_p <- NA_real_
        }

        if (all_groups_normal && !is.na(lev_p) && lev_p > 0.05) {
          aov_res <- stats::aov(
            as.formula(paste0("polyA_length ~ ", grouping_factor)),
            data = subdf
          )
          p_val       <- summary(aov_res)[[1]][["Pr(>F)"]][1]
          chosen_test <- "ANOVA"
        } else if (all_groups_normal && all(grp_tbl >= 2)) {
          welch_res <- stats::oneway.test(
            as.formula(paste0("polyA_length ~ ", grouping_factor)),
            data = subdf,
            var.equal = FALSE
          )
          p_val       <- welch_res$p.value
          chosen_test <- "Welch ANOVA"
        } else {
          p_val       <- suppressWarnings(
                           stats::kruskal.test(
                             as.formula(paste0("polyA_length ~ ", grouping_factor)),
                             data = subdf
                           )$p.value
                         )
          chosen_test <- "Kruskal–Wallis"
        }
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

  out_df <- out_df[c(which_level, "test_used", "p_value", "padj")]
  return(out_df)
}