#' Plot Density of Poly(A) Tail Lengths with Statistical Testing
#'
#' Generates a normalized density plot of poly(A) tail length distributions by group,
#' overlays group-wise median or mean lines, and performs:
#'  - No test if only one group  
#'  - Wilcoxon rank‐sum or Student’s t-test if two groups  
#'  - Kruskal–Wallis or one‐way ANOVA if more than two groups  
#'
#' @param polyA_table A data frame with per-read poly(A) tail lengths and grouping,  
#'   as produced by \code{\link{get_gene_id}()}. Must include columns  
#'   \code{polyA_length} (numeric), the grouping column, and \code{sample_name}.  
#' @param stats A string: \code{"median"} or \code{"mean"} for dashed lines.  
#' @param grouping_column Name of the column in \code{polyA_table} defining groups.  
#' @return A list with components:  
#'   \describe{  
#'     \item{test}{A message or an \code{htest} object (t-test/Wilcoxon/ANOVA/Kruskal–Wallis).}  
#'     \item{plot}{A \code{ggplot2} density plot with dashed lines at group means/medians.}  
#'   }  
#' @export
#'
plot_density <- function(polyA_table = get_gene_id_out,
                         stats = "median",
                         grouping_column = "group") {
                          
  if (missing(polyA_table) || !is.data.frame(polyA_table)) {
    stop("'polyA_table' must be defined as a data frame.")
  }
  if (missing(stats) || !(stats %in% c("median", "mean"))) {
    stop("'stats' must be 'median' or 'mean'.")
  }
  required_cols <- c("polyA_length", grouping_column, "sample_name")
  if (!all(required_cols %in% colnames(polyA_table))) {
    stop(sprintf(
      "The polyA_table must contain columns: %s.",
      paste(required_cols, collapse = ", ")
    ))
  }

  max_x <- stats::quantile(polyA_table$polyA_length, 0.99, na.rm = TRUE)

  density_plot <- ggplot2::ggplot(polyA_table,
    ggplot2::aes(
      x = polyA_length,
      y = ..ndensity..,
      color = !!rlang::sym(grouping_column)
    )
  ) +
    ggplot2::geom_density() +
    ggplot2::stat_density(geom = "line", position = "identity", size = 1) +
    ggplot2::labs(
      title = NULL,
      x = "PolyA length",
      y = "Density (normalized)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::coord_cartesian(xlim = c(0, max_x))

  if (stats == "median") {
    medians <- polyA_table %>%
      dplyr::group_by(!!rlang::sym(grouping_column)) %>%
      dplyr::summarise(median_polyA_length = stats::median(polyA_length, na.rm = TRUE))
    density_plot <- density_plot +
      ggplot2::geom_vline(
        data = medians,
        ggplot2::aes(xintercept = median_polyA_length, color = !!rlang::sym(grouping_column)),
        linetype = "dashed", size = 1
      )
  } else {
    means <- polyA_table %>%
      dplyr::group_by(!!rlang::sym(grouping_column)) %>%
      dplyr::summarise(mean_polyA_length = base::mean(polyA_length, na.rm = TRUE))
    density_plot <- density_plot +
      ggplot2::geom_vline(
        data = means,
        ggplot2::aes(xintercept = mean_polyA_length, color = !!rlang::sym(grouping_column)),
        linetype = "dashed", size = 1
      )
  }

  formula <- as.formula(paste0("polyA_length ~ ", grouping_column))
  group_count <- length(unique(polyA_table[[grouping_column]]))
  test_result <- NULL

  if (group_count == 1) {
    test_result <- "No test: only one group present."
  } else if (group_count == 2) {
    norm_p <- polyA_table %>%
      dplyr::group_by(!!rlang::sym(grouping_column)) %>%
      dplyr::summarise(p = stats::shapiro.test(polyA_length)$p.value) %>%
      dplyr::pull(p)
    if (all(norm_p > 0.05)) {
      test_result <- stats::t.test(formula, data = polyA_table)
    } else {
      test_result <- stats::wilcox.test(formula, data = polyA_table)
    }
  } else {
    aov_mod <- stats::aov(formula, data = polyA_table)
    res <- stats::residuals(aov_mod)
    bart_p <- stats::bartlett.test(formula, data = polyA_table)$p.value
    if (stats::shapiro.test(res)$p.value > 0.05 && bart_p > 0.05) {
      test_result <- summary(aov_mod)
    } else {
      test_result <- stats::kruskal.test(formula, data = polyA_table)
    }
  }

  return(list(test = test_result, plot = density_plot))
}