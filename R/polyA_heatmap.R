#' Create a Heatmap and Dendrogram of PolyA Tail Distributions
#'
#' This function computes a binned count matrix of poly(A) tail lengths for 
#' each sample group, then generates a heatmap and corresponding row 
#' dendrogram. Supports raw counts or normalized percentages.
#'
#' @param polyA_table A data frame containing per-read poly(A) tail lengths, 
#'   as produced by \code{\link{get_polyA}()}. Must include columns 
#'   \code{polyA_length} (numeric) and the grouping column.
#' @param grouping_factor A character string giving the name of the column in 
#'   \code{polyA_table} that defines sample groups. Values will become row 
#'   names and heatmap rows.
#' @param frame A single positive numeric value specifying the bin width for 
#'   counting poly(A) lengths. Defines the intervals 
#'   \code{[1, 1+frame)}, \code{[1+frame, 1+2*frame)}, etc.
#' @param select A character string, either \code{"base"} for raw count values 
#'   or \code{"normalized"} for row-wise percent normalization.
#' @param heatmap_color A character string selecting one of six color palettes 
#'   for the heatmap gradient: \code{"green_red"}, \code{"red_green"}, 
#'   \code{"blue_green"}, \code{"green_blue"}, \code{"blue_red"} or 
#'   \code{"red_blue"}.
#'
#' @return A list with components:
#'   \describe{
#'     \item{matrix}{A numeric matrix of counts or percentages (rows = groups, columns = length bins).}
#'     \item{heatmap}{A \code{Heatmap} object from ComplexHeatmap, drawn and ready for display.}
#'     \item{tree}{A \code{ggtree} object representing the row dendrogram.}
#'   }
#'
#' @export
#'
#' @details
#' The function performs these steps:
#' \itemize{
#'   \item Validates input data frame, grouping column, and parameters.
#'   \item Computes the 0.99 quantile of \code{polyA_length} to define the 
#'     maximum bin range.
#'   \item Creates equally spaced bins of width \code{frame} and counts reads 
#'     whose tail lengths fall into each bin for each group.
#'   \item Optionally normalizes each row to percentages if 
#'     \code{select = "normalized"}.
#'   \item Calculates decile quantiles of the matrix values to set color breaks.
#'   \item Constructs a heatmap with \code{ComplexHeatmap::Heatmap()}, using 
#'     a user-selected color ramp from \code{circlize::colorRamp2()}.
#'   \item Draws the heatmap and extracts the row dendrogram, converting it 
#'     to a \code{ggtree} plot.
#' }
#'
#' @section Color Palettes:
#' Available palettes map low values to the first color and high values to the 
#' last. For example, \code{"green_red"} transitions from dark green to red.
#'
#' @section Error Handling:
#' The function stops with an informative message if:
#' \itemize{
#'   \item \code{polyA_table} is not a data frame or lacks \code{polyA_length}.
#'   \item \code{grouping_factor} is missing or not a factor/character column.
#'   \item \code{frame} is not a single positive number.
#'   \item \code{select} is not one of \code{"base"} or \code{"normalized"}.
#'   \item \code{heatmap_color} is not one of the six valid palette names.
#' }
#'
#' @seealso
#' \code{\link{get_polyA}} for generating \code{polyA_table};  
#' \code{\link[circlize]{colorRamp2}} for palette creation;  
#' \code{\link[ComplexHeatmap]{Heatmap}} for heatmap generation;  
#' \code{\link[ggtree]{ggtree}} for dendrogram plotting.
#'
#' @importFrom stats quantile
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap draw row_dend
#' @importFrom ggtree ggtree
#' @importFrom grid unit

polyA_heatmap <- function(polyA_table=polyA_table,grouping_factor = "group", frame=10, select = "base", heatmap_color ="green_red"){

  if (!is.data.frame(polyA_table)) {
    stop("polyA_table must be a data frame.")
  }
  
  if (!"polyA_length" %in% names(polyA_table)) {
    stop("polyA_table must contain the column 'polyA_length'.")
  }
  
  if (!grouping_factor %in% names(polyA_table)) {
    stop(paste("Column", grouping_factor, "does not exist in polyA_table."))
  }

  if (!is.numeric(frame) || frame <= 0 || length(frame) != 1) {
    stop("frame must be a single, positive number.")
  }
  
  if (!is.numeric(polyA_table$polyA_length)) {
    stop("The column 'polyA_length' must contain numeric values.")
  }
  
  if (!is.factor(polyA_table[[grouping_factor]]) && !is.character(polyA_table[[grouping_factor]])) {
    stop(paste("Column", grouping_factor, "must be a factor or character string."))
  }

  if (base::missing(select)) {
    stop('The select predictions are required. Please provide a valid argument - "base" or "normalized".',
         call. = FALSE)
  }

if (select == "base") {

max_density <- stats::quantile(polyA_table$polyA_length, probs = c(0.99))[1]
breaks <- base::seq(1, base::max(max_density) + frame, by = frame)
matrix_data <- base::matrix(0, nrow = base::length(base::unique(polyA_table[[grouping_factor]])), ncol = base::length(breaks) - 1)
base::rownames(matrix_data) <- base::unique(polyA_table[[grouping_factor]])
base::colnames(matrix_data) <- base::paste0(breaks[-base::length(breaks)], "-", breaks[-1] - 1)

for (i in 1:nrow(polyA_table)) {
  sample_name <- polyA_table[[grouping_factor]][i]
  polyA_length <- polyA_table$polyA_length[i]
  
  for (j in 1:(base::length(breaks) - 1)) {
    if (polyA_length >= breaks[j] && polyA_length < breaks[j + 1]) {
      matrix_data[sample_name, j] <- matrix_data[sample_name, j] + 1
      break
    }
  }
}

quantiles <- stats::quantile(matrix_data, probs = base::seq(0, 1, length.out = 10))


valid_colors <- c("red_green",
                    "green_red",
                    "blue_green",
                    "green_blue",
                    "blue_red",
                    "red_blue")
  
  if (!(heatmap_color %in% valid_colors)) {
    stop("Invalid heatmap_color argument. Please choose from red_green, green_red, blue_green, green_blue, blue_red, red_blue valid options.")
  }
  
  if (heatmap_color=="green_red") {
    colors <- circlize::colorRamp2(quantiles,
                                   c("#006d2c","#2ca25f","#66c2a4","#99d8c9","#ccece6","#fdd49e","#fdbb84", "#fc8d59","#e34a33","#b30000"))
  }
  
  if (heatmap_color=="red_green") {
    colors <- circlize::colorRamp2(quantiles,
                                   c("#b30000","#e34a33","#fc8d59","#fdbb84","#fdd49e","#ccece6","#99d8c9","#66c2a4","#2ca25f","#006d2c"))
  }
  
  if (heatmap_color=="green_blue") {
    colors <- circlize::colorRamp2(quantiles,
                                   c("#006d2c","#2ca25f","#66c2a4","#99d8c9","#ccece6","#c6dbef","#9ecae1","#6baed6","#3182bd","#08519c"))
  }
  
  if (heatmap_color=="blue_green") {
    colors <- circlize::colorRamp2(quantiles,
                                   c("#08519c","#3182bd","#6baed6","#9ecae1","#c6dbef","#ccece6","#99d8c9","#66c2a4", "#2ca25f","#006d2c"))
  }
  
  if (heatmap_color=="blue_red") {
    colors <- circlize::colorRamp2(quantiles,
                                   c("#08519c","#3182bd","#6baed6","#9ecae1","#c6dbef","#fdd49e","#fdbb84","#fc8d59","#e34a33","#b30000"))
  }
  
  if (heatmap_color=="red_blue") {
    colors <- circlize::colorRamp2(quantiles,
                                   c("#b30000","#e34a33","#fc8d59","#fdbb84","#fdd49e","#c6dbef","#9ecae1","#6baed6","#3182bd","#08519c"))
  }


heatmap_polyA <- ComplexHeatmap::Heatmap(matrix_data,name = "poly(A)\ncounts",
                                        col = colors,
                                        column_dend_height = unit(3, "cm"),
                                        row_dend_width = unit(3, "cm"),
                                        cluster_rows = T,
                                        cluster_columns = F)
heatmap_polyA <- ComplexHeatmap::draw(heatmap_polyA)
tree <- ComplexHeatmap::row_dend(heatmap_polyA)
tree <- ggtree::ggtree(tree)

base::message("Success !!!")

return(list(matrix = matrix_data, heatmap = heatmap_polyA, tree = tree))
}

if (select == "normalized") {

max_density <- stats::quantile(polyA_table$polyA_length, probs = c(0.99))[1]
breaks <- base::seq(1, base::max(max_density) + frame, by = frame)
matrix_data <- base::matrix(0, nrow = base::length(base::unique(polyA_table[[grouping_factor]])), ncol = base::length(breaks) - 1)
base::rownames(matrix_data) <- base::unique(polyA_table[[grouping_factor]])
base::colnames(matrix_data) <- base::paste0(breaks[-base::length(breaks)], "-", breaks[-1] - 1)

for (i in 1:nrow(polyA_table)) {
  sample_name <- polyA_table[[grouping_factor]][i]
  polyA_length <- polyA_table$polyA_length[i]
  
  for (j in 1:(base::length(breaks) - 1)) {
    if (polyA_length >= breaks[j] && polyA_length < breaks[j + 1]) {
      matrix_data[sample_name, j] <- matrix_data[sample_name, j] + 1
      break
    }
  }
}

matrix_data_transform <- matrix_data

for (i in 1:base::nrow(matrix_data_transform)) {
  for (j in 1:base::ncol(matrix_data_transform)) {
    matrix_data_transform[i,j] <- (base::sum(matrix_data[i,j])/base::sum(matrix_data[i,])) * 100
  }
}

quantiles <- stats::quantile(matrix_data_transform, probs = base::seq(0, 1, length.out = 10))


valid_colors <- c("red_green",
                    "green_red",
                    "blue_green",
                    "green_blue",
                    "blue_red",
                    "red_blue")
  
  if (!(heatmap_color %in% valid_colors)) {
    stop("Invalid heatmap_color argument. Please choose from red_green, green_red, blue_green, green_blue, blue_red, red_blue valid options.")
  }
  
  if (heatmap_color=="green_red") {
    colors <- circlize::colorRamp2(quantiles,
                                   c("#006d2c","#2ca25f","#66c2a4","#99d8c9","#ccece6","#fdd49e","#fdbb84", "#fc8d59","#e34a33","#b30000"))
  }
  
  if (heatmap_color=="red_green") {
    colors <- circlize::colorRamp2(quantiles,
                                   c("#b30000","#e34a33","#fc8d59","#fdbb84","#fdd49e","#ccece6","#99d8c9","#66c2a4","#2ca25f","#006d2c"))
  }
  
  if (heatmap_color=="green_blue") {
    colors <- circlize::colorRamp2(quantiles,
                                   c("#006d2c","#2ca25f","#66c2a4","#99d8c9","#ccece6","#c6dbef","#9ecae1","#6baed6","#3182bd","#08519c"))
  }
  
  if (heatmap_color=="blue_green") {
    colors <- circlize::colorRamp2(quantiles,
                                   c("#08519c","#3182bd","#6baed6","#9ecae1","#c6dbef","#ccece6","#99d8c9","#66c2a4", "#2ca25f","#006d2c"))
  }
  
  if (heatmap_color=="blue_red") {
    colors <- circlize::colorRamp2(quantiles,
                                   c("#08519c","#3182bd","#6baed6","#9ecae1","#c6dbef","#fdd49e","#fdbb84","#fc8d59","#e34a33","#b30000"))
  }
  
  if (heatmap_color=="red_blue") {
    colors <- circlize::colorRamp2(quantiles,
                                   c("#b30000","#e34a33","#fc8d59","#fdbb84","#fdd49e","#c6dbef","#9ecae1","#6baed6","#3182bd","#08519c"))
  }

heatmap_polyA <- ComplexHeatmap::Heatmap(matrix_data_transform,name = "poly(A)\npercent",
                                        col = colors,
                                        column_dend_height = unit(3, "cm"),
                                        row_dend_width = unit(3, "cm"),
                                        cluster_rows = T,
                                        cluster_columns = F)
heatmap_polyA <- ComplexHeatmap::draw(heatmap_polyA)
tree <- ComplexHeatmap::row_dend(heatmap_polyA)
tree <- ggtree::ggtree(tree)

base::message("Success !!!")

return(list(matrix = matrix_data_transform, heatmap = heatmap_polyA, tree = tree))
}
}