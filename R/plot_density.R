#' Creating a density plot
#'
#' @param calculate_statistics_out the table was created using the calculate_statistics function
#' @return A table object.
#' @export
#'

plot_density <- function(polyA_table=get_gene_id_out,stats="median"){
  
  
  if (missing(polyA_table)) {
    stop("'polyA_table' must be defined.")
  }

  if (missing(stats)) {
    stop("'stats' must be defined. median or mean")
  }
  
  if(!all(c("polyA_length", "group") %in% colnames(polyA_table))) {
    stop("The polyA_table must contain the columns 'polyA_length' and 'group'.")
  }
  
max_density <- quantile(polyA_table$polyA_length, probs = c(0.99))[1]

density_plot <- ggplot(polyA_table, aes(x = polyA_length, y=..ndensity.., color = group)) +
                geom_density() +
                stat_density(geom = "line", position = "identity", size = 1) +
                labs(title = "Density plot of polyA lengths", x = "PolyA length", y = "Density (normalized)") +
                theme_bw() + 
                coord_cartesian(xlim = c(0, max_density)) 

if (stats=="median") {
  
  medians <- polyA_table %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(median_polyA_length = stats::median(polyA_length))
  
  density_plot <- density_plot + geom_vline(data = medians, aes(xintercept = median_polyA_length, color = group), 
                                            linetype = "dashed", size = 1)
}

if (stats=="mean") {

means <- polyA_table %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(mean_polyA_length = base::mean(polyA_length))

density_plot <- density_plot + geom_vline(data = means, aes(xintercept = mean_polyA_length, color = group), 
           linetype = "dashed", size = 1)

}

return(density_plot)

}
