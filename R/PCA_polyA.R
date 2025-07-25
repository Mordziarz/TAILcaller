#' Perform PCA on Poly(A) Count Matrix and Plot the First Two Components
#'
#' This function conducts principal component analysis (PCA) on a matrix of 
#' poly(A) length statistics (e.g., counts or summary measures) and produces a 
#' scatter plot of the first two principal components colored by a grouping factor.
#'
#' @param get_matrix_out A numeric matrix where rows represent samples and 
#'   columns represent features (e.g., transcript or gene bins). Typically 
#'   the output of \code{\link{get_matrix}()}. Rows are transposed internally 
#'   so that samples become observations in PCA.
#' @param samples_table A data frame containing sample metadata. Must include 
#'   a column named \code{sample_name} matching the row names of 
#'   \code{get_matrix_out} prior to transposition, and the grouping column.
#' @param grouping_factor A character string specifying the name of the column 
#'   in \code{samples_table} used to color points in the PCA plot (e.g., 
#'   \code{"group"}).
#'
#' @author Mateusz Mazdziarz
#'
#' @importFrom stats prcomp summary
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab ggtitle theme_bw theme element_blank
#' @importFrom base merge round paste0 t

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