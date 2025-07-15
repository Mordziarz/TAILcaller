#' Plot Density of Poly(A) Tail Lengths with Summary Statistics and Group Comparisons
#'
#' Generates a normalized density plot of poly(A) tail length distributions for each sample group,
#' overlays group-specific mean or median reference lines, and conducts appropriate group-comparison tests:
#' \itemize{
#'   \item \strong{Two groups}: Student's t-test (if normality and equal variances), Welch's t-test (if normality, unequal variances), or Wilcoxon rank-sum test.
#'   \item \strong{More than two groups}: One-way ANOVA + Tukey HSD (if normality and equal variances), Welch's ANOVA (if normality, unequal variances), or Kruskal–Wallis + Dunn's test.
#' }
#'
#' @param polyA_table A \code{data.frame} with per-read poly(A) tail lengths. Must contain:
#'   \describe{
#'     \item{\code{polyA_length} (numeric)}{Tail length of each read.}
#'     \item{\code{sample_name} (character or factor)}{Sample identifier (checked but not directly used for plot/tests in this version).}
#'     \item{\code{<grouping_column>} (character or factor)}{Grouping variable for comparison.}
#'   }
#' @param stats Character; summary statistic for vertical reference lines. One of \code{"median"} or \code{"mean"}. Defaults to \code{"median"}.
#' @param grouping_column Character; name of the column in \code{polyA_table} to use for grouping, coloring, and statistical tests. Defaults to \code{"group"}.
#'
#' @return A named \code{list} with components:
#'   \describe{
#'     \item{\code{plot}}{A \code{ggplot} object showing normalized density curves per group with dashed mean/median lines.}
#'     \item{\code{normality}}{A \code{data.frame} of per-group normality test p-values (Shapiro–Wilk if n ≤ 5000, Lilliefors if n > 5000) and group sizes (\code{n}).}
#'     \item{\code{variance}}{Result of Levene’s test for homogeneity of variance (from \code{car::leveneTest}) or \code{NULL} if fewer than two groups or groups lack sufficient observations (min 2 per group).}
#'     \item{\code{test}}{Statistical comparison results, which can be:
#'       \itemize{
#'         \item A \code{character} message if fewer than two groups ("Only one group; no comparison performed.").
#'         \item For two groups: a \code{htest} object from \code{stats::t.test()} or \code{stats::wilcox.test()}.
#'         \item For more than two groups: a \code{list} containing \code{anova} and \code{tukey} results (from \code{stats::aov()} and \code{stats::TukeyHSD()}) for parametric ANOVA, or \code{kruskal} and \code{dunn} results (from \code{stats::kruskal.test()} and \code{dunn.test::dunn.test()}) for non-parametric tests, or \code{oneway.test} object for Welch's ANOVA.
#'       }
#'     }
#'   }
#'
#' @details
#' The function executes the following steps:
#'   \enumerate{
#'     \item Validates inputs and required columns: \code{polyA_length}, \code{sample_name}, and the specified \code{grouping_column}.
#'     \item Computes the 99th percentile of \code{polyA_length} to set a reasonable x-axis limit for the density plot, preventing outliers from stretching the plot excessively.
#'     \item Generates a \code{ggplot2} density plot:
#'       \itemize{
#'         \item Plots normalized density curves for each group using \code{geom_density()} and \code{stat_density()}.
#'         \item Adds dashed vertical lines at group means or medians (as specified by \code{stats}) using \code{geom_vline()}.
#'       }
#'     \item Performs per-group normality tests (Shapiro–Wilk for n ≤ 5000, Lilliefors for n > 5000) to assess the normal distribution assumption.
#'     \item Conducts Levene’s test for homogeneity of variance across groups using \code{car::leveneTest(center = median)} (more robust to non-normality than the mean).
#'     \item Determines the number of distinct groups (\code{ngroups}) in the \code{grouping_column} and selects the appropriate statistical test:
#'       \itemize{
#'         \item \code{ngroups < 2}: No comparison is performed, and a message is returned.
#'         \item \code{ngroups == 2}:
#'           \itemize{
#'             \item \strong{Student's t-test}: If both groups pass normality tests (\code{p > 0.05}) and Levene’s test (\code{p > 0.05}).
#'             \item \strong{Welch's t-test}: If both groups pass normality tests, but Levene’s test fails (\code{p <= 0.05}). This test does not assume equal variances.
#'             \item \strong{Wilcoxon rank-sum test}: If normality assumption is violated for any group.
#'           }
#'         \item \code{ngroups > 2}:
#'           \itemize{
#'             \item \strong{One-way ANOVA + Tukey HSD}: If all groups pass normality tests and Levene’s test passes.
#'             \item \strong{Welch's ANOVA}: If all groups pass normality tests, but Levene’s test fails. Post-hoc comparisons for Welch's ANOVA (e.g., Games-Howell) are not included in this function but should be considered if this test is chosen.
#'             \item \strong{Kruskal–Wallis + Dunn’s test}: If normality assumption is violated for any group, or if Levene's test could not be performed (e.g., due to insufficient observations per group for variance calculation).
#'           }
#'       }
#'   \end{enumerate}
#'
#' @section Statistical Tests Used:
#'   \itemize{
#'     \item \strong{Normality}: \code{stats::shapiro.test()} (for N <= 5000) or \code{nortest::lillie.test()} (for N > 5000).
#'     \item \strong{Homogeneity of Variance}: \code{car::leveneTest()}.
#'     \item \strong{Two groups}:
#'       \itemize{
#'         \item Parametric: \code{stats::t.test()} (Student's and Welch's variants).
#'         \item Non-parametric: \code{stats::wilcox.test()}.
#'       }
#'     \item \strong{More than two groups}:
#'       \itemize{
#'         \item Parametric: \code{stats::aov()} followed by \code{stats::TukeyHSD()} (for equal variances), or \code{stats::oneway.test()} with \code{var.equal = FALSE} (Welch's ANOVA for unequal variances).
#'         \item Non-parametric: \code{stats::kruskal.test()} followed by \code{dunn.test::dunn.test()}.
#'       }
#'   }
#'
#' @seealso
#'   \code{\link[dplyr]{group_by}}, \code{\link[dplyr]{summarise}}, \code{\link[dplyr]{n_distinct}}, \code{\link[dplyr]{n}} for data manipulation;
#'   \code{\link[ggplot2]{ggplot}}, \code{\link[ggplot2]{aes}}, \code{\link[ggplot2]{geom_density}}, \code{\link[ggplot2]{stat_density}}, \code{\link[ggplot2]{geom_vline}}, \code{\link[ggplot2]{labs}}, \code{\link[ggplot2]{theme_bw}}, \code{\link[ggplot2]{coord_cartesian}} for plotting;
#'   \code{\link[stats]{quantile}}, \code{\link[stats]{median}}, \code{\link[base]{mean}}, \code{\link[stats]{shapiro.test}}, \code{\link[stats]{t.test}}, \code{\link[stats]{wilcox.test}}, \code{\link[stats]{aov}}, \code{\link[stats]{TukeyHSD}}, \code{\link[stats]{kruskal.test}}, \code{\link[stats]{oneway.test}} for statistical functions;
#'   \code{\link[nortest]{lillie.test}} for normality testing;
#'   \code{\link[car]{leveneTest}} for homogeneity of variance;
#'   \code{\link[dunn.test]{dunn.test}} for post-hoc non-parametric tests;
#'   \code{\link[rlang]{sym}} for tidy evaluation.
#'
#' @author Mateusz Mazdziarz
#'
#' @importFrom ggplot2 ggplot aes geom_density stat_density geom_vline labs theme_bw coord_cartesian
#' @importFrom rlang sym
#' @importFrom dplyr group_by summarise n_distinct n
#' @importFrom stats quantile median mean t.test wilcox.test aov TukeyHSD kruskal.test oneway.test
#' @importFrom nortest lillie.test
#' @importFrom car leveneTest
#' @importFrom dunn.test dunn.test
#' @export

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

  ngroups <- dplyr::n_distinct(polyA_table[[grouping_column]])

  group_counts <- polyA_table %>%
    dplyr::group_by(!!rlang::sym(grouping_column)) %>%
    dplyr::summarise(count = dplyr::n(), .groups = 'drop')

  force_non_parametric <- any(group_counts$count < 3)

  normality <- NULL
  levene_res <- NULL
  group_test <- NULL

  if (ngroups < 2) {
    group_test <- "Only one group; no comparison performed."
  } else {
    if (!force_non_parametric) {
      normality <- polyA_table %>%
        dplyr::group_by(!!rlang::sym(grouping_column)) %>%
        dplyr::summarise(
          n = dplyr::n(),
          p_value = if (n() <= 5000) {
            stats::shapiro.test(polyA_length)$p.value
          } else {
            nortest::lillie.test(polyA_length)$p.value
          },
          .groups = 'drop'
        )

      if (all(group_counts$count >= 2)) {
        levene_res <- tryCatch({
          car::leveneTest(
            as.formula(paste0("polyA_length ~ ", grouping_column)),
            data = polyA_table,
            center = median 
          )
        }, error = function(e) {
          warning(paste0("Levene's test failed: ", e$message, ". Proceeding with non-parametric test."))
          NULL 
        })
      } else {
        levene_res <- NULL
      }
    }

    if (ngroups == 2) {
      if (force_non_parametric ||
          !all(normality$p_value > 0.05, na.rm = TRUE) || 
          is.null(levene_res) || levene_res[["Pr(>F)"]][1] <= 0.05) { 
        group_test <- stats::wilcox.test(
          as.formula(paste0("polyA_length ~ ", grouping_column)),
          data = polyA_table
        )
      } else {
        group_test <- stats::t.test(
          as.formula(paste0("polyA_length ~ ", grouping_column)),
          data = polyA_table,
          var.equal = TRUE
        )
      }
    } else { 
      if (force_non_parametric ||
          !all(normality$p_value > 0.05, na.rm = TRUE) || 
          is.null(levene_res) || levene_res[["Pr(>F)"]][1] <= 0.05) { 
        kruskal_res <- stats::kruskal.test(
          as.formula(paste0("polyA_length ~ ", grouping_column)),
          data = polyA_table
        )
        dunn_res <- dunn.test::dunn.test(
          x = polyA_table$polyA_length,
          g = polyA_table[[grouping_column]],
          method = "bonferroni",
          altp = TRUE 
        )
        group_test <- list(kruskal = kruskal_res, dunn = dunn_res)
      } else {
        aov_res <- stats::aov(
          as.formula(paste0("polyA_length ~ ", grouping_column)),
          data = polyA_table
        )
        anova_sum <- summary(aov_res)
        tukey_res <- stats::TukeyHSD(aov_res)
        group_test <- list(anova_summary = anova_sum, tukey_hsd = tukey_res)
      }
    }
  }

  return(list(
    plot = density_plot,
    normality = normality, 
    variance = levene_res, 
    test = group_test
  ))
}


