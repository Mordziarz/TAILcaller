#' Plot Density of Poly(A) Tail Lengths with Summary Statistics
#'
#' This function generates a density plot of poly(A) tail length distributions  
#' for each sample group, and optionally overlays group-wise mean or median  
#' lines. It also performs an appropriate nonparametric test (Wilcoxon or  
#' Kruskal–Wallis) based on the number of groups.
#'
#' @param polyA_table A data frame containing per-read poly(A) tail lengths,  
#'   as produced by \code{\link{get_gene_id}()}. Must include columns  
#'   \code{polyA_length} (numeric), the specified \code{grouping_column}, and  
#'   \code{sample_name}.
#' @param stats A character string selecting the summary statistic for vertical  
#'   lines: \code{"median"} or \code{"mean"}. Lines are drawn at each group’s  
#'   median or mean poly(A) length.
#' @param grouping_column A character string naming the column in  
#'   \code{polyA_table} that defines sample groups. Used for color mapping,  
#'   density grouping, and statistical testing.
#'
#' @return A list with components:  
#'   \describe{  
#'     \item{wilcox_test}{The result of a Wilcoxon rank–sum test (two groups)  
#'       or Kruskal–Wallis test (more than two groups), or a message if testing  
#'       is not possible.}  
#'     \item{plot}{A \code{ggplot} object showing the density curves and dashed  
#'       lines at group means or medians.}  
#'   }
#'
#' @export
#'
#' @details  
#' The function performs these steps:  
#' \itemize{  
#'   \item Validates that \code{polyA_table}, \code{stats}, and required columns  
#'     exist.  
#'   \item Computes the 99th percentile of \code{polyA_length} to set the x-axis  
#'     limit.  
#'   \item Plots density curves (normalized) for each group using  
#'     \code{geom_density()} and \code{stat_density()}.  
#'   \item Overlays dashed vertical lines at group-specific means or medians,  
#'     computed via \code{dplyr::summarise()} when \code{stats} is  
#'     \code{"mean"} or \code{"median"}.  
#'   \item Determines the number of unique groups:  
#'     \describe{  
#'       \item{n = 1}{Returns a message that testing cannot be performed.}  
#'       \item{n = 2}{Performs a Wilcoxon rank–sum test.}  
#'       \item{n > 2}{Performs a Kruskal–Wallis test.}  
#'     }  
#' }
#'
#' @section Statistical Testing:  
#' Uses \code{wilcox.test()} for two groups and \code{kruskal.test()} for more  
#' than two groups. If the grouping column is missing or contains only one level,  
#' an explanatory message is returned instead of a test result.
#'
#' @seealso  
#' \code{\link{get_gene_id}} for preparing \code{polyA_table};  
#' \code{\link[dplyr]{group_by}}, \code{\link[dplyr]{summarise}} for summary stats;  
#' \code{\link[ggplot2]{geom_density}}, \code{\link[ggplot2]{geom_vline}} for plotting;  
#' \code{\link[stats]{wilcox.test}}, \code{\link[stats]{kruskal.test}} for tests.
#'
#' @importFrom ggplot2 ggplot aes geom_density stat_density geom_vline labs theme_bw coord_cartesian
#' @importFrom rlang sym
#' @importFrom dplyr group_by summarise
#' @importFrom stats quantile median wilcox.test kruskal.test
#' @importFrom base stop missing unique nrow

plot_density <- function(polyA_table=get_gene_id_out,stats="median",grouping_column="group"){
  
  
  if (missing(polyA_table)) {
    stop("'polyA_table' must be defined.")
  }

  if (missing(stats)) {
    stop("'stats' must be defined. median or mean")
  }
  
  if(!all(c("polyA_length", "group", "sample_name") %in% colnames(polyA_table))) {
    stop("The polyA_table must contain the columns 'polyA_length', 'group' and 'sample_name'.")
  }
  
max_density <- stats::quantile(polyA_table$polyA_length, probs = c(0.99))[1]

density_plot <- ggplot2::ggplot(polyA_table, ggplot2::aes(x = polyA_length, y=..ndensity.., color = !!rlang::sym(grouping_column))) +
                ggplot2::geom_density() +
                ggplot2::stat_density(geom = "line", position = "identity", size = 1) +
                ggplot2::labs(title = "Density plot of polyA lengths", x = "PolyA length", y = "Density (normalized)") +
                ggplot2::theme_bw() + 
                ggplot2::coord_cartesian(xlim = c(0, max_density))

if (stats=="median") {
  
  medians <- polyA_table %>%
    dplyr::group_by(!!rlang::sym(grouping_column)) %>%
    dplyr::summarise(median_polyA_length = stats::median(polyA_length))
  
  density_plot <- density_plot + ggplot2::geom_vline(data = medians, ggplot2::aes(xintercept = median_polyA_length, color = !!rlang::sym(grouping_column)), 
                                            linetype = "dashed", size = 1)
}

if (stats=="mean") {

means <- polyA_table %>%
  dplyr::group_by(!!rlang::sym(grouping_column)) %>%
  dplyr::summarise(mean_polyA_length = base::mean(polyA_length))

density_plot <- density_plot + ggplot2::geom_vline(data = means, ggplot2::aes(xintercept = mean_polyA_length, color = !!rlang::sym(grouping_column)), 
           linetype = "dashed", size = 1)

}

  normality <- polyA_table %>%
    dplyr::group_by(!!rlang::sym(grouping_column)) %>%
    dplyr::summarise(
      n         = dplyr::n(),
      p_value   = if (n() <= 5000) {
                    stats::shapiro.test(polyA_length)$p.value
                  } else {
                    nortest::lillie.test(polyA_length)$p.value
                  }
    )
  
  ngroups <- dplyr::n_distinct(polyA_table[[grouping_column]])
  if (ngroups >= 2) {
    levene_res <- car::leveneTest(
      as.formula(paste0("polyA_length ~ ", grouping_column)),
      data = polyA_table
    )
  } else {
    levene_res <- NULL
  }
  
  if (ngroups < 2) {
    group_test <- "Only one group; no comparison performed."
  } else if (ngroups == 2) {
    use_t_test <- all(normality$p_value > 0.05) &&
                  (!is.null(levene_res) && levene_res[["Pr(>F)"]][1] > 0.05)
    if (use_t_test) {
      group_test <- stats::t.test(
        as.formula(paste0("polyA_length ~ ", grouping_column)),
        data = polyA_table
      )
    } else {
      group_test <- stats::wilcox.test(
        as.formula(paste0("polyA_length ~ ", grouping_column)),
        data = polyA_table
      )
    }
  } else {
    if (all(normality$p_value > 0.05) &&
        !is.null(levene_res) && levene_res[["Pr(>F)"]][1] > 0.05) {
      aov_res <- stats::aov(
        as.formula(paste0("polyA_length ~ ", grouping_column)),
        data = polyA_table
      )
      anova_res <- summary(aov_res)
      tukey_res <- stats::TukeyHSD(aov_res)
      group_test <- list(anova = anova_res, tukey = tukey_res)
    } else {
      kruskal_res <- stats::kruskal.test(
        as.formula(paste0("polyA_length ~ ", grouping_column)),
        data = polyA_table
      )
      dunn_res <- dunn.test::dunn.test(
        x = polyA_table$polyA_length,
        g = polyA_table[[grouping_column]],
        method = "bonferroni"
      )
      group_test <- list(kruskal = kruskal_res, dunn = dunn_res)
    }
  }
  
  return(list(
    plot        = density_plot,
    normality   = normality,
    variance    = levene_res,
    comparison  = group_test
  ))
}

