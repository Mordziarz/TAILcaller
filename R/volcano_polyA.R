#' Creating a volcano plot
#'
#' @param calculate_statistics_out the table was created using the calculate_statistics function
#' @return A table object.
#' @export
#'

volcano_polyA <- function(calculate_statistics_out=calculate_statistics_out){
  
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
  ggplot2::theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 27),
        legend.text = ggplot2::element_text(size = 10),
        legend.title = ggplot2::element_text(size = 12),
        legend.text.align = 0,
        legend.title.align = 0,
        axis.title = ggplot2::element_text(size = rel(4.0))) +
  ggplot2::scale_color_manual(values = c("Collapsed" = "green2", "Expansion" = "red2", "No significant" = "grey35"),
                     name = "PolyA tail length") +  
  ggplot2::scale_fill_manual(values = c("Collapsed" = "green4", "Expansion" = "red4", "No significant" = "grey50"),
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