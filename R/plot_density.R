#' Plot Density of Poly(A) Tail Lengths with Summary Statistics and Group Comparisons
#'
#' Generates a normalized density plot of poly(A) tail length distributions for each sample group,
#' overlays group-specific mean or median reference lines, and conducts appropriate group-comparison tests:
#' \itemize{
#'   \item \strong{Two groups}: Student’s t-test (if normality and equal variances), Welch’s t-test (if normality and unequal variances), or Wilcoxon rank-sum test.
#'   \item \strong{More than two groups}: One-way ANOVA + Tukey HSD (if normality and equal variances), Welch’s ANOVA + Games-Howell (if normality and unequal variances), or Kruskal–Wallis + Dunn’s test.
#' }
#' The function automatically determines the most suitable test based on data characteristics
#' (number of groups, normality, homogeneity of variances) and gracefully handles cases with very few
#' observations per group.
#'
#' @param polyA_table A \code{data.frame} with per-read poly(A) tail lengths. Must contain:
#'   \describe{
#'     \item{\code{polyA_length} (numeric)}{Tail length of each read.}
#'     \item{\code{sample_name} (character or factor)}{Sample identifier (checked but not directly used for plotting or statistical tests in this version).}
#'     \item{\code{<grouping_factor>} (character or factor)}{Grouping variable for comparison, e.g., experimental condition.}
#'   }
#' @param stats Character; summary statistic for vertical reference lines on the plot. Must be one of
#'   \code{"median"} or \code{"mean"}. Defaults to \code{"median"}.
#' @param grouping_factor Character; name of the column in \code{polyA_table} to use for defining groups,
#'   coloring the density curves, and performing statistical comparisons. Defaults to \code{"group"}.
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
#' @importFrom rstatix games_howell_test

plot_density <- function(polyA_table=get_gene_id_out,stats="median",grouping_factor="group") {
  
  
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

density_plot <- ggplot2::ggplot(polyA_table, ggplot2::aes(x = polyA_length, y=..ndensity.., color = !!rlang::sym(grouping_factor))) +
                ggplot2::geom_density() +
                ggplot2::stat_density(geom = "line", position = "identity", size = 1) +
                ggplot2::labs(title = "Density plot of polyA lengths", x = "PolyA length", y = "Density (normalized)") +
                ggplot2::theme_bw() + 
                ggplot2::coord_cartesian(xlim = c(0, max_density))

if (stats=="median") {
  
  medians <- polyA_table %>%
    dplyr::group_by(!!rlang::sym(grouping_factor)) %>%
    dplyr::summarise(median_polyA_length = stats::median(polyA_length))
  
  density_plot <- density_plot + ggplot2::geom_vline(data = medians, ggplot2::aes(xintercept = median_polyA_length, color = !!rlang::sym(grouping_factor)), 
                                            linetype = "dashed", size = 1)
}

if (stats=="mean") {

means <- polyA_table %>%
  dplyr::group_by(!!rlang::sym(grouping_factor)) %>%
  dplyr::summarise(mean_polyA_length = base::mean(polyA_length))

density_plot <- density_plot + ggplot2::geom_vline(data = means, ggplot2::aes(xintercept = mean_polyA_length, color = !!rlang::sym(grouping_factor)), 
           linetype = "dashed", size = 1)

}

  ngroups <- dplyr::n_distinct(polyA_table[[grouping_factor]])

  group_counts <- polyA_table %>%
    dplyr::group_by(!!rlang::sym(grouping_factor)) %>%
    dplyr::summarise(count = dplyr::n(), .groups = 'drop')

  force_non_parametric <- any(group_counts$count < 3)

  normality <- NULL
  levene_res <- NULL
  group_test <- NULL
  all_normal <- FALSE
  homogeneous_variances <- FALSE


  if (ngroups < 2) {
  group_test <- "Only one group; no comparison performed."
} else {
  if (!force_non_parametric) {
    normality <- polyA_table %>%
      dplyr::group_by(!!rlang::sym(grouping_factor)) %>%
      dplyr::summarise(
        n = dplyr::n(),
        p_value = (function() {
          current_group_name <- as.character(dplyr::cur_group()[[1]])
          data_vector <- polyA_length

          if (dplyr::n() < 3) {
            warning(paste0("Not enough observations (n=", dplyr::n(), ") in group '", current_group_name, "' for normality test. Assuming non-normal."))
            return(0)
          }
          if (length(unique(data_vector)) == 1) {
            warning(paste0("All 'polyA_length' values are identical in group '", current_group_name, "'. Normality test cannot be performed. Assuming non-normal."))
            return(0) 
          }

          tryCatch({
            if (dplyr::n() <= 5000) {
              stats::shapiro.test(data_vector)$p.value
            } else {
              if (!requireNamespace("nortest", quietly = TRUE)) {
                stop("Package 'nortest' is required for Lilliefor's test when n > 5000. Please install it.")
              }
              nortest::lillie.test(data_vector)$p.value
            }
          }, error = function(e) {
            warning(paste0("Error in normality test for group '", current_group_name, "' (n=", dplyr::n(), "): ", e$message, ". Assuming non-normal."))
            return(0) 
          })
        })(), 
        .groups = 'drop'
      )
    all_normal <- all(normality$p_value > 0.05, na.rm = TRUE)

    if (all(group_counts$count >= 2)) {
      levene_res <- tryCatch(
        {
          car::leveneTest(
            as.formula(paste0("polyA_length ~ ", grouping_factor)),
            data = polyA_table,
            center = median
          )
        },
        error = function(e) {
          warning(paste0("Levene's test failed: ", e$message, ". Homogeneity of variance cannot be assessed."))
          NULL
        }
      )
      homogeneous_variances <- !is.null(levene_res) && levene_res[["Pr(>F)"]][1] > 0.05
    } else {
      warning("Not enough observations in some groups to perform Levene's test for homogeneity of variance. Assuming unequal variances for parametric tests if normality holds, or forcing non-parametric test if group sizes are too small (<3).")
      homogeneous_variances <- FALSE
    }
  }


    if (ngroups == 2) {
      if (force_non_parametric || !all_normal) {
        group_test <- stats::wilcox.test(
          as.formula(paste0("polyA_length ~ ", grouping_factor)),
          data = polyA_table
        )
      } else if (all_normal && homogeneous_variances) {
        group_test <- stats::t.test(
          as.formula(paste0("polyA_length ~ ", grouping_factor)),
          data = polyA_table,
          var.equal = TRUE
        )
      } else if (all_normal && !homogeneous_variances) {
        group_test <- stats::t.test(
          as.formula(paste0("polyA_length ~ ", grouping_factor)),
          data = polyA_table,
          var.equal = FALSE
        )
      }
    } else {
      if (force_non_parametric || !all_normal) {
        kruskal_res <- stats::kruskal.test(
          as.formula(paste0("polyA_length ~ ", grouping_factor)),
          data = polyA_table
        )
        dunn_res <- dunn.test::dunn.test(
          x = polyA_table$polyA_length,
          g = polyA_table[[grouping_factor]],
          method = "bonferroni",
          altp = TRUE
        )
        group_test <- list(kruskal = kruskal_res, dunn = dunn_res)
      } else if (all_normal && homogeneous_variances) {
        aov_res <- stats::aov(
          as.formula(paste0("polyA_length ~ ", grouping_factor)),
          data = polyA_table
        )
        anova_sum <- summary(aov_res)
        tukey_res <- stats::TukeyHSD(aov_res)
        group_test <- list(anova_summary = anova_sum, tukey_hsd = tukey_res)
      } else if (all_normal && !homogeneous_variances) {
        welch_res <- stats::oneway.test(
          as.formula(paste0("polyA_length ~ ", grouping_factor)),
          data = polyA_table,
          var.equal = FALSE
        )
        games_howell_res <- rstatix::games_howell_test(
          data = polyA_table,
          formula = as.formula(paste0("polyA_length ~ ", grouping_factor)),
          conf.level = 0.95,
          detailed = FALSE
        )
        group_test <- list(welch_anova = welch_res, games_howell = games_howell_res)
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