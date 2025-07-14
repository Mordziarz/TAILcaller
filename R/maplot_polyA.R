#' Create a MA-Style Plot of Poly(A) Tail Changes
#'
#' This function generates an MA-style scatter plot where the x-axis is the  
#' average poly(A) tail length across control and treatment groups, and the  
#' y-axis is the log₂ fold change. Points are colored by whether tails are  
#' “Collapsed,” “Expansion,” or “No significant” based on significance criteria.
#'
#' @param calculate_statistics_out A data frame produced by \code{\link{calculate_statistics}()},  
#'   containing at minimum:  
#'   \describe{  
#'     \item{mean_group_ctr}{Numeric mean poly(A) length in control group.}  
#'     \item{mean_group_trt}{Numeric mean poly(A) length in treatment group.}  
#'     \item{Log2FC}{Numeric log₂ fold change.}  
#'     \item{padj}{Numeric adjusted p-value.}  
#'   }
#' @param collapsed_color A single color name or hex code for points classified as “Collapsed.”  
#' @param expansion_color A single color name or hex code for points classified as “Expansion.”  
#'
#' @return A \code{ggplot} object representing an MA plot of poly(A) tail changes,  
#'   with points colored and shaded by significance category.
#'
#' @export
#'
#' @details  
#' Steps performed by the function:  
#' \itemize{  
#'   \item Validates that \code{calculate_statistics_out} is provided and has the required columns.  
#'   \item Converts \code{Log2FC} and \code{padj} to numeric and filters out NA values.  
#'   \item Classifies each feature as:  
#'     \describe{  
#'       \item{Collapsed}{\code{Log2FC < 0} & \code{padj < 0.05}}  
#'       \item{Expansion}{\code{Log2FC > 0} & \code{padj < 0.05}}  
#'       \item{No significant}{all others}  
#'     }  
#'   \item Computes \code{mean = (mean_group_ctr + mean_group_trt)/2} for each feature.  
#'   \item Orders points so significant changes plot on top.  
#'   \item Builds a \code{ggplot2} scatter plot with a horizontal line at \code{y = 0},  
#'     custom axis labels, themes, and legends.  
#'   \item Applies manual scales for color, fill, and alpha using the provided  
#'     \code{collapsed_color} and \code{expansion_color}, with grey for non-significant points.  
#' }
#'
#' @seealso  
#' \code{\link{calculate_statistics}} for input data;  
#' \code{\link[ggplot2]{geom_point}}, \code{\link[ggplot2]{scale_color_manual}},  
#' \code{\link[ggplot2]{scale_fill_manual}} for styling.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline xlab ylab theme element_text theme_bw guides guide_legend scale_color_manual scale_fill_manual scale_alpha_manual
#' @importFrom stats as.numeric
#' @importFrom base ifelse paste0 stop missing order

maplot_polyA <- function(calculate_statistics_out=calculate_statistics_out,collapsed_color = "green",expansion_color = "red"){
  
  if (missing(calculate_statistics_out)) {
    stop("'calculate_statistics_out' must be defined.")
  }
  
  '%!in%' <- function(x,y)!('%in%'(x,y))

  calculate_statistics_out$padj <- base::as.numeric(calculate_statistics_out$padj)
  calculate_statistics_out <- calculate_statistics_out[calculate_statistics_out$padj %!in% NA,]
  
  calculate_statistics_out$Log2FC <- base::as.numeric(calculate_statistics_out$Log2FC)
  calculate_statistics_out <- calculate_statistics_out[calculate_statistics_out$Log2FC %!in% NA,]
  
  calculate_statistics_out$PolyA_tail_length <- "No significant"
  calculate_statistics_out$PolyA_tail_length <- base::ifelse(calculate_statistics_out$Log2FC < 0 & calculate_statistics_out$padj <0.05,  "Collapsed" , calculate_statistics_out$PolyA_tail_length)
  calculate_statistics_out$PolyA_tail_length <- base::ifelse(calculate_statistics_out$Log2FC > 0 & calculate_statistics_out$padj <0.05,  "Expansion" , calculate_statistics_out$PolyA_tail_length)
  calculate_statistics_out$Polya_difference <- "No significant"
  calculate_statistics_out$Polya_difference <- base::ifelse(calculate_statistics_out$Log2FC > 0 & calculate_statistics_out$padj <0.05,  "Significant" , calculate_statistics_out$Polya_difference)
  calculate_statistics_out$Polya_difference <- base::ifelse(calculate_statistics_out$Log2FC < 0 & calculate_statistics_out$padj <0.05,  "Significant" , calculate_statistics_out$Polya_difference)
  calculate_statistics_out$order <- base::ifelse(calculate_statistics_out$Polya_difference %in% "No significant", 2, 1)
  
  calculate_statistics_out$mean <- NA
  calculate_statistics_out$mean <- (calculate_statistics_out$mean_group_ctr + calculate_statistics_out$mean_group_trt)/2
  calculate_statistics_out <- calculate_statistics_out[order(calculate_statistics_out$Polya_difference),]
  
  polya_maplot <- ggplot2::ggplot(calculate_statistics_out, 
                                    ggplot2::aes(x = mean, y = Log2FC, 
                                                 color = PolyA_tail_length, 
                                                 fill = PolyA_tail_length,
                                                 alpha = PolyA_tail_length)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept=0,linetype="dotted", size=0.5,col="red") +
    ggplot2::xlab(expression("mean(polyA tail length)")) + ggplot2::ylab("log2FoldChange") +
    ggplot2::theme(
      plot.title = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 20),
      axis.title.y = ggplot2::element_text(size = 20, angle = 90),
      axis.text.y = ggplot2::element_text(size = 20),
      axis.text.x = ggplot2::element_text(size = 20)) +
    ggplot2::theme(legend.position = "right", plot.title = ggplot2::element_text(hjust = 0.5, size = 27),
                   legend.text = ggplot2::element_text(size = 10),
                   legend.title = ggplot2::element_text(size = 12),
                   legend.text.align = 0,
                   legend.title.align = 0,
                   axis.title = ggplot2::element_text(size = ggplot2::rel(4.0))) +
    ggplot2::scale_color_manual(values = c("Collapsed" = base::paste0(collapsed_color,2), "Expansion" = base::paste0(expansion_color,2), "No significant" = "grey35"),
                                name = "PolyA tail length") +  
    ggplot2::scale_fill_manual(values = c("Collapsed" = base::paste0(collapsed_color,4), "Expansion" = base::paste0(expansion_color,4), "No significant" = "grey50"),
                               name = "PolyA tail length") +
    ggplot2::scale_alpha_manual(values = c("Collapsed" = 0.5, "Expansion" = 0.5, "No significant" = 0.2),
                                name = "PolyA tail length") +  
    ggplot2::theme_bw() +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1, override.aes = list(shape = 15, size = 4)),
      fill = ggplot2::guide_legend(order = 1),
      alpha = ggplot2::guide_legend(order = 1)
    )
  
  base::message("Done !!!!")
  
  return(polya_maplot)
}
