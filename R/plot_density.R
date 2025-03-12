#' Creating a density plot
#'
#' @param polyA_table the table was created using the get_polyA function
#' @param stats adopts a value of mean or median
#' @param grouping_column grouping factor that splits the plot by the given column
#' @return a plot object.
#' @export
#'

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

formula <- base::paste0("polyA_length ~ ", grouping_column)
res_stats <- stats::wilcox.test(stats::as.formula(formula), data=polyA_table)

return(list(wilcox_test = res_stats, plot = density_plot))

}
