polyA_heatmap <- function(polyA_table=polyA_table,grouping_factor = "group", jump = 5, frame=10){

  if (!is.data.frame(polyA_table)) {
    stop("polyA_table must be a data frame.")
  }
  
  if (!"polyA_length" %in% names(polyA_table)) {
    stop("polyA_table must contain the column 'polyA_length'.")
  }
  
  if (!grouping_factor %in% names(polyA_table)) {
    stop(paste("Column", grouping_factor, "does not exist in polyA_table."))
  }
  
  if (!is.numeric(jump) || jump <= 0 || length(jump) != 1) {
    stop("jump must be a single, positive number.")
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

max_density <- stats::quantile(polyA_table$polyA_length, probs = c(0.99))[1]
breaks <- base::seq(1, base::max(max_density) + frame, by = jump)
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

colors <- circlize::colorRamp2(quantiles,
                               c("#08519c","#3182bd","#6baed6","#9ecae1","#c6dbef","#fdd49e","#fdbb84","#fc8d59","#e34a33","#b30000"))


heatmap_rscu <- ComplexHeatmap::Heatmap(matrix_data,name = "poly(A) count",
                                        col = colors,
                                        column_dend_height = unit(3, "cm"),
                                        row_dend_width = unit(3, "cm"),
                                        cluster_rows = T,
                                        cluster_columns = F)
heatmap_rscu <- ComplexHeatmap::draw(heatmap_rscu)
tree <- ComplexHeatmap::row_dend(heatmap_rscu)
tree <- ggtree::ggtree(tree)

base::message("Success !!!")

return(list(matrix = matrix_data, heatmap = heatmap_rscu, tree = tree))
}