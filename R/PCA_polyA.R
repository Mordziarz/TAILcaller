#' Creating a PCA based on the output from get_matrix function
#'
#' @param get_matrix_out the table was created using the get_RSCU function
#' @param samples_table the table with bam description
#' @param grouping_factor the name of the column in the table that divides the experiment into groups.
#' @return A plot object.
#' @export
#'

PCA_polyA <- function(get_matrix_out=get_matrix_out,samples_table=samples_table,grouping_factor="group"){
  
  if (base::missing(get_matrix_out)) {
    stop("The get_matrix_out predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }
  samples_table$names <- samples_table$sample_name
  get_matrix_out <- base::t(get_matrix_out)
  pca <- stats::prcomp(get_matrix_out, scale. = FALSE)
  pca_data <- base::as.data.frame(pca$x)
  pca_data$names <- base::rownames(get_matrix_out)
  pca_data <- base::merge(pca_data,samples_table,by.x="names",by.y="names",all.x = T)
  
  PCA_plot <- ggplot2::ggplot(pca_data, aes(x = PC1, y = PC2, label = names, color = !!rlang::sym(grouping_factor))) +
    ggplot2::geom_point(size = 3) +
    ggplot2::xlab(base::paste0("PC1 (", base::round(base::summary(pca)$importance[2,1]*100, 1), "%)")) +
    ggplot2::ylab(base::paste0("PC2 (", base::round(base::summary(pca)$importance[2,2]*100, 1), "%)")) +
    ggplot2::ggtitle("PCA") +
    ggplot2::theme_bw()+
    ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  
  return(PCA_plot)
  
  base::message(base::paste0("Success"))
  
}