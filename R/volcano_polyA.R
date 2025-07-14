#' Create a Volcano Plot of Poly(A) Tail Length Changes
#'
#' This function generates a volcano plot to visualize differential poly(A) tail
#' length changes between conditions. It classifies transcripts as “Collapsed,”
#' “Expansion,” or “No significant” based on log₂ fold change and adjusted
#' p-value thresholds, then plots log₂ fold change versus –log₁₀(padj) with
#' custom colors and transparency.
#'
#' @param calculate_statistics_out A data frame produced by \code{\link{calculate_statistics}()}, 
#'   containing at minimum the columns \code{Log2FC} (numeric) and \code{padj} (numeric).
#' @param collapsed_color A single color name or hex code used for points classified as “Collapsed.”
#' @param expansion_color A single color name or hex code used for points classified as “Expansion.”
#'
#' @return A \code{ggplot} object representing the volcano plot, with points colored
#'   and shaded by poly(A) tail length change category.
#'
#' @export
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Validates that \code{calculate_statistics_out} is provided and contains 
#'     \code{Log2FC} and \code{padj} columns.
#'   \item Converts \code{Log2FC} and \code{padj} to numeric and removes any NA values.
#'   \item Classifies each transcript into:
#'     \describe{
#'       \item{Collapsed}{\code{Log2FC < 0} and \code{padj < 0.05}}
#'       \item{Expansion}{\code{Log2FC > 0} and \code{padj < 0.05}}
#'       \item{No significant}{all others}
#'     }
#'   \item Sets plotting aesthetics: vertical line at \code{x = 0}, horizontal line at 
#'     \code{y = –log10(0.05)}, axis labels, theme, and legend formatting.
#'   \item Applies manual scales for color, fill, and alpha using the provided 
#'     \code{collapsed_color} and \code{expansion_color}, with greys for non-significant points.
#' }
#'
#' @section Color Mapping:
#' The \code{collapsed_color} and \code{expansion_color} arguments are used to build  
#' two opacity variants (alpha = 0.5 and 0.25) for fill and color aesthetics.
#'
#' @seealso
#' \code{\link{calculate_statistics}} for generating the input data;
#' \code{\link[ggplot2]{geom_point}}, \code{\link[ggplot2]{scale_color_manual}},
#' \code{\link[ggplot2]{scale_fill_manual}} for plot customization.
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline xlab ylab theme element_text
#' @importFrom ggplot2 scale_color_manual scale_fill_manual scale_alpha_manual theme_bw guides guide_legend
#' @importFrom stats as.numeric
#' @importFrom base ifelse paste0 stop missing

volcano_polyA <- function(calculate_statistics_out=calculate_statistics_out,collapsed_color = "green",expansion_color = "red"){
  
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


polya_volcano <- ggplot2::ggplot(calculate_statistics_out, 
                                  ggplot2::aes(x = Log2FC, y = -log10(padj), 
                                               color = PolyA_tail_length, 
                                               fill = PolyA_tail_length,
                                               alpha = PolyA_tail_length)) +
  ggplot2::geom_point() +
  ggplot2::geom_vline(xintercept = 0, linetype = "dotted", size = .5, col = "red") +
  ggplot2::geom_hline(yintercept = 1.30103, linetype = "dotted", size = 0.5, col = "red") +
  ggplot2::xlab(expression("log2FoldChange")) + ggplot2::ylab("-log10(padj)") +
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

return(polya_volcano)
}