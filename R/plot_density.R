#' Plot Density of Poly(A) Tail Lengths with Summary Statistics and Group Comparisons
#'
#' Generates a normalized density plot of poly(A) tail length distributions for each sample group,
#' overlays group-specific mean or median reference lines, and conducts appropriate group-comparison tests:
#' \itemize{
#'   \item \strong{Two groups}: Student's t-test (if normality and equal variances), or Wilcoxon rank-sum test.
#'   \item \strong{More than two groups}: One-way ANOVA + Tukey HSD (if normality and equal variances), or Kruskal–Wallis + Dunn's test.
#' }
#' This function automatically determines the most suitable test based on data characteristics (number of groups, normality, homogeneity of variances) and gracefully handles cases with very few observations per group.
#'
#' @param polyA_table A \code{data.frame} with per-read poly(A) tail lengths. Must contain:
#'   \describe{
#'     \item{\code{polyA_length} (numeric)}{Tail length of each read.}
#'     \item{\code{sample_name} (character or factor)}{Sample identifier (checked but not directly used for plotting or statistical tests in this version).}
#'     \item{\code{<grouping_column>} (character or factor)}{Grouping variable for comparison, e.g., experimental condition.}
#'   }
#' @param stats Character; summary statistic for vertical reference lines on the plot. Must be one of \code{"median"} or \code{"mean"}. Defaults to \code{"median"}.
#' @param grouping_column Character; name of the column in \code{polyA_table} to use for defining groups, coloring the density curves, and performing statistical comparisons. Defaults to \code{"group"}.
#'
#' @return A named \code{list} with the following components:
#'   \describe{
#'     \item{\code{plot}}{A \code{ggplot} object displaying normalized density curves for each group, with dashed lines indicating group-specific means or medians.}
#'     \item{\code{normality}}{A \code{data.frame} containing per-group normality test p-values (Shapiro–Wilk for n ≤ 5000, Lilliefors for n > 5000) and group sizes (\code{n}). Returns \code{NULL} if any group has fewer than 3 observations (forcing a non-parametric test).}
#'     \item{\code{variance}}{Result of Levene’s test for homogeneity of variance (a \code{data.frame} from \code{car::leveneTest}). Returns \code{NULL} if fewer than two groups, if any group has fewer than 2 observations, or if a non-parametric test was forced due to small group sizes.}
#'     \item{\code{test}}{Statistical comparison results, which can be:
#'       \itemize{
#'         \item A \code{character} message ("Only one group; no comparison performed.") if \code{ngroups < 2}.
#'         \item For two groups: a \code{htest} object from \code{stats::t.test()} (Student's t-test) or \code{stats::wilcox.test()} (Wilcoxon rank-sum test).
#'         \item For more than two groups: a \code{list} containing \code{anova_summary} (from \code{summary(stats::aov())}) and \code{tukey_hsd} (from \code{stats::TukeyHSD()}) for parametric ANOVA, or a \code{list} containing \code{kruskal} (from \code{stats::kruskal.test()}) and \code{dunn} (from \code{dunn.test::dunn.test()}) for non-parametric tests.
#'       }
#'     }
#'   }
#'
#' @details
#' The function operates through the following steps:
#'   \enumerate{
#'     \item **Input Validation**: Ensures that `polyA_table` is defined and contains the necessary columns (`polyA_length`, `sample_name`, and the specified `grouping_column`). Also validates `stats` argument.
#'     \item **X-axis Limit Calculation**: Determines the 99th percentile of `polyA_length` to set a suitable upper limit for the plot's x-axis, preventing extreme outliers from distorting the visual representation.
#'     \item **Density Plot Generation**:
#'       \itemize{
#'         \item Creates normalized density curves for each group using `ggplot2::geom_density()`.
#'         \item Overlays dashed vertical lines representing either the mean or median (based on the `stats` parameter) for each group using `ggplot2::geom_vline()`.
#'       }
#'     \item **Group Size Check and Test Selection Logic**:
#'       \itemize{
#'         \item Calculates the number of unique groups (`ngroups`) and counts observations per group.
#'         \item If **any group has fewer than 3 observations**, or if Levene's test cannot be performed (e.g., fewer than 2 observations in any group), the function automatically defaults to a **non-parametric test path** (`force_non_parametric` is set to `TRUE`). In such cases, normality and homogeneity of variance tests are skipped, and their results (in `normality` and `variance` components of the return list) will be `NULL`.
#'       }
#'     \item **Statistical Test Execution**:
#'       \itemize{
#'         \item \code{ngroups < 2}: No statistical comparison is performed.
#'         \item \code{ngroups == 2}:
#'           \itemize{
#'             \item \strong{Student's t-test}: Chosen if `force_non_parametric` is `FALSE`, all groups pass normality tests (`p > 0.05`), and Levene’s test for homogeneity of variance passes (`p > 0.05`).
#'             \item \strong{Wilcoxon rank-sum test}: Performed if `force_non_parametric` is `TRUE`, or if normality/homogeneity of variance assumptions are violated or cannot be reliably assessed.
#'           }
#'         \item \code{ngroups > 2}:
#'           \itemize{
#'             \item \strong{One-way ANOVA + Tukey HSD}: Selected if `force_non_parametric` is `FALSE`, all groups pass normality tests (`p > 0.05`), and Levene’s test for homogeneity of variance passes (`p > 0.05`).
#'             \item \strong{Kruskal–Wallis + Dunn’s test}: Performed if `force_non_parametric` is `TRUE`, or if normality/homogeneity of variance assumptions are violated or cannot be reliably assessed.
#'           }
#'       }
#'   \end{enumerate}
#'
#' @section Statistical Tests Used:
#'   \itemize{
#'     \item \strong{Normality}: `stats::shapiro.test()` (for N <= 5000) or `nortest::lillie.test()` (for N > 5000). These are skipped if `force_non_parametric` is `TRUE`.
#'     \item \strong{Homogeneity of Variance}: `car::leveneTest(center = median)`. This is skipped if `force_non_parametric` is `TRUE`.
#'     \item \strong{Two groups}:
#'       \itemize{
#'         \item Parametric: `stats::t.test()` (Student's variant with `var.equal = TRUE`).
#'         \item Non-parametric: `stats::wilcox.test()`.
#'       }
#'     \item \strong{More than two groups}:
#'       \itemize{
#'         \item Parametric: `stats::aov()` followed by `stats::TukeyHSD()`.
#'         \item Non-parametric: `stats::kruskal.test()` followed by `dunn.test::dunn.test()`.
#'       }
#'   }
#'
#' @seealso
#'   `\link[dplyr]{group_by}`, `\link[dplyr]{summarise}`, `\link[dplyr]{n_distinct}`, `\link[dplyr]{n}` for data manipulation;
#'   `\link[ggplot2]{ggplot}`, `\link[ggplot2]{aes}`, `\link[ggplot2]{geom_density}`, `\link[ggplot2]{stat_density}`, `\link[ggplot2]{geom_vline}`, `\link[ggplot2]{labs}`, `\link[ggplot2]{theme_bw}`, `\link[ggplot2]{coord_cartesian}` for plotting;
#'   `\link[stats]{quantile}`, `\link[stats]{median}`, `\link[base]{mean}`, `\link[stats]{shapiro.test}`, `\link[stats]{t.test}`, `\link[stats]{wilcox.test}`, `\link[stats]{aov}`, `\link[stats]{TukeyHSD}`, `\link[stats]{kruskal.test}` for statistical functions;
#'   `\link[nortest]{lillie.test}` for normality testing;
#'   `\link[car]{leveneTest}` for homogeneity of variance;
#'   `\link[dunn.test]{dunn.test}` for post-hoc non-parametric tests;
#'   `\link[rlang]{sym}` for tidy evaluation.
#'
#' @author Mateusz Mazdziarz
#'
#' @importFrom ggplot2 ggplot aes geom_density stat_density geom_vline labs theme_bw coord_cartesian
#' @importFrom rlang sym
#' @importFrom dplyr group_by summarise n_distinct n
#' @importFrom stats quantile median mean t.test wilcox.test aov TukeyHSD kruskal.test
#' @importFrom nortest lillie.test
#' @importFrom car leveneTest
#' @importFrom dunn.test dunn.test
#' @export


plot_density <- function(polyA_table=get_gene_id_out,stats="median",grouping_column="group") {
  
  
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
      if (force_non_parametric || !all(normality$p_value > 0.05, na.rm = TRUE)) {
        group_test <- stats::wilcox.test(
          as.formula(paste0("polyA_length ~ ", grouping_column)),
          data = polyA_table
        )
      } else {
        if (!is.null(levene_res) && levene_res[["Pr(>F)"]][1] <= 0.05) { 
          group_test <- stats::t.test(
            as.formula(paste0("polyA_length ~ ", grouping_column)),
            data = polyA_table,
            var.equal = FALSE # Korekta Welcha
          )
        } else {
          group_test <- stats::t.test(
            as.formula(paste0("polyA_length ~ ", grouping_column)),
            data = polyA_table,
            var.equal = TRUE 
          )
        }
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