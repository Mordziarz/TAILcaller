#’ Plot Density of Poly(A) Tail Lengths with Summary Statistics and Group Comparisons
#’
#’ Generates a normalized density plot of poly(A) tail length distributions for each sample group, optionally overlays group‐specific mean or median lines, and conducts appropriate group‐comparison tests (t-test or Wilcoxon for two groups; ANOVA or Kruskal–Wallis plus post-hoc for more than two).
#’
#’ @param polyA_table Data frame of per‐read poly(A) tail lengths with at least these columns:
#’ \itemize{
#’ \item{\code{polyA_length} (numeric)}{—Tail length of each read.}
#’ \item{\code{sample_name} (character or factor)}{—Identifier for each sample.}
#’ \item{\code{<grouping_column>} (character or factor)}{—Grouping factor for comparisons.}
#’ }
#’ @param stats Character; summary statistic for vertical reference lines. Must be \code{"median"} or \code{"mean"}.
#’ @param grouping_column Character; name of the column in \code{polyA_table} to use for grouping, density coloring, and tests.
#’
#’ @return A named list with components:
#’ \describe{
#’ \item{\code{plot}}{A \code{ggplot} object showing density curves per group with dashed lines at group means or medians.}
#’ \item{\code{normality}}{Data frame of Shapiro–Wilk (if n ≤ 5000) or Lilliefors (if n > 5000) p-values per group.}
#’ \item{\code{variance}}{Result of Levene’s test for homogeneity of variance, or \code{NULL} if only one group.}
#’ \item{\code{test}}{For two groups, a t-test or Wilcoxon rank-sum test; for more than two, ANOVA + Tukey HSD or Kruskal–Wallis + Dunn’s test; or a message if only one group.}
#’ }
#’
#’ @details
#’ The function performs these steps:
#’ \itemize{
#’ \item Validates that \code{polyA_table}, \code{stats}, and required columns exist.
#’ \item Computes the 99th percentile of \code{polyA_length} to set the x‐axis limit.
#’ \item Plots normalized density curves per group with \code{geom_density()} and \code{stat_density()}.
#’ \item Adds dashed vertical lines at group means or medians using \code{geom_vline()}.
#’ \item Tests normality per group (Shapiro–Wilk or Lilliefors) and homogeneity of variance (Levene’s test).
#’ \item Determines number of groups (\code{ngroups}):
#’ \describe{
#’ \item{\code{ngroups < 2}}{Returns a message; no comparisons.}
#’ \item{\code{ngroups == 2}}{If normality and equal variances hold, performs two-sample t-test; otherwise Wilcoxon test.}
#’ \item{\code{ngroups > 2}}{If assumptions hold, performs one-way ANOVA with Tukey HSD post-hoc; otherwise Kruskal–Wallis with Dunn’s test.}
#’ }
#’ }
#’
#’ @section Statistical Tests:
#’ \itemize{
#’ \item Two groups: \code{t.test()} or \code{wilcox.test()}.
#’ \item More than two: \code{aov()} + \code{TukeyHSD()} or \code{kruskal.test()} + \code{dunn.test()}.
#’ }
#’
#’ @seealso
#’ \code{\link{get_gene_id}} for generating \code{polyA_table};
#’ \code{\link[dplyr]{group_by}}, \code{\link[dplyr]{summarise}} for summary stats;
#’ \code{\link[ggplot2]{geom_density}}, \code{\link[ggplot2]{geom_vline}};
#’ \code{\link[stats]{shapiro.test}}, \code{\link[nortest]{lillie.test}}, \code{\link[car]{leveneTest}}, \code{\link[stats]{t.test}}, \code{\link[stats]{wilcox.test}}, \code{\link[stats]{aov}}, \code{\link[stats]{TukeyHSD}}, \code{\link[stats]{kruskal.test}};
#’ \code{\link[dunn.test]{dunn.test}}.
#’
#’ @importFrom ggplot2 ggplot aes geom_density stat_density geom_vline labs theme_bw coord_cartesian
#’ @importFrom rlang sym
#’ @importFrom dplyr group_by summarise n_distinct n
#’ @importFrom stats quantile median t.test wilcox.test aov TukeyHSD kruskal.test
#’ @importFrom nortest lillie.test
#’ @importFrom car leveneTest
#’ @importFrom dunn.test dunn.test
#’ @export

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
    test        = group_test
  ))
}

