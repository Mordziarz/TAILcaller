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
#' @importFrom stats quantile median mean t.test wilcox.test aov TukeyHSD kruskal.test oneway.test
#' @importFrom nortest lillie.test
#' @importFrom car leveneTest
#' @importFrom dunn.test dunn.test
#' @importFrom rstatix games_howell_test
#' @import data.table

plot_density <- function(polyA_table = get_gene_id_out, stats = "median", grouping_factor = "group") {
  
  if (missing(polyA_table)) stop("'polyA_table' must be defined.")
  
  dt <- as.data.table(polyA_table)
  dt[, (grouping_factor) := as.factor(get(grouping_factor))]
  
  if(!all(c("polyA_length", grouping_factor, "sample_name") %in% colnames(dt))) {
    stop(paste("The polyA_table must contain the columns 'polyA_length',", grouping_factor, "and 'sample_name'."))
  }
  
  max_density <- quantile(dt$polyA_length, probs = 0.99, na.rm = TRUE)[[1]]
  ref_stats <- dt[, .(stat_val = if(stats == "median") as.numeric(median(polyA_length)) else mean(polyA_length)), 
                  by = grouping_factor]
  
  density_plot <- ggplot(dt, aes(x = polyA_length, y = after_stat(ndensity), color = .data[[grouping_factor]])) +
    geom_density(linewidth = 1) +
    geom_vline(data = ref_stats, aes(xintercept = stat_val, color = .data[[grouping_factor]]), 
               linetype = "dashed", linewidth = 1) +
    labs(title = "Density plot of polyA lengths", 
         x = "PolyA length", 
         y = "Density (normalized)",
         color = "Group") +
    theme_bw() + 
    coord_cartesian(xlim = c(0, max_density)) +
    guides(color = guide_legend(override.aes = list(linetype = "solid", linewidth = 1))) +
    theme(legend.key = element_blank())
  
  group_info <- dt[, .(n = .N), by = grouping_factor]
  ngroups <- nrow(group_info)
  force_non_parametric <- any(group_info$n < 3)
  
  normality <- NULL
  levene_res <- NULL
  group_test <- NULL
  
  if (ngroups >= 2) {
    if (!force_non_parametric) {
      normality <- dt[, {
        val <- polyA_length
        p_val <- tryCatch({
          if (.N <= 5000) stats::shapiro.test(val)$p.value
          else nortest::lillie.test(val)$p.value
        }, error = function(e) 0)
        .(n = .N, p_value = p_val)
      }, by = grouping_factor]
      
      all_normal <- all(normality$p_value > 0.05, na.rm = TRUE)
      form <- as.formula(paste0("polyA_length ~ ", grouping_factor))
      levene_res <- tryCatch(car::leveneTest(form, data = dt, center = median), error = function(e) NULL)
      homogeneous_variances <- !is.null(levene_res) && levene_res[["Pr(>F)"]][1] > 0.05
    }
    
    if (ngroups == 2) {
      if (force_non_parametric || !all_normal) {
        group_test <- stats::wilcox.test(form, data = dt)
      } else {
        group_test <- stats::t.test(form, data = dt, var.equal = homogeneous_variances)
      }
    } else {
      if (force_non_parametric || !all_normal) {
        kruskal_res <- stats::kruskal.test(form, data = dt)
        dunn_res <- dunn.test::dunn.test(x = dt$polyA_length, g = dt[[grouping_factor]], method = "bonferroni", altp = TRUE)
        group_test <- list(kruskal = kruskal_res, dunn = dunn_res)
      } else if (homogeneous_variances) {
        aov_res <- stats::aov(form, data = dt)
        group_test <- list(anova_summary = summary(aov_res), tukey_hsd = stats::TukeyHSD(aov_res))
      } else {
        group_test <- list(
          welch_anova = stats::oneway.test(form, data = dt, var.equal = FALSE),
          games_howell = rstatix::games_howell_test(data = dt, formula = form)
        )
      }
    }
  } else {
    group_test <- "Only one group; no comparison performed."
  }
  
  return(list(plot = density_plot, normality = as.data.frame(normality), variance = levene_res, test = group_test))
}