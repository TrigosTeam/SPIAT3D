plot_grid_metrics_discrete3D <- function(grid_metrics, metric_colname) {
  
  ## Check input parameters
  if (!(is.character(metric_colname) && metric_colname %in% c("proportion", "entropy"))) {
    stop("`metric_colname` is not 'proportion' or 'entropy'.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  
  ## Define low, medium and high categories
  # Low: between 0 and 1/3
  # Medium: between 1/3 and 2/3
  # High: between 2/3 and 1
  
  grid_metrics$rank <- ifelse(is.na(grid_metrics[[metric_colname]]), "na",
                              ifelse(grid_metrics[[metric_colname]] < 1/3, "low",
                                     ifelse(grid_metrics[[metric_colname]] < 2/3, "medium", "high")))
  grid_metrics$rank <- factor(grid_metrics$rank, c("low", "medium", "high", "na"))
  
  fig <- plot_ly(grid_metrics,
                 type = "scatter3d",
                 mode = 'markers',
                 x = ~x_coord,
                 y = ~y_coord,
                 z = ~z_coord,
                 color = ~rank,
                 colors = c("#AEB6E5", "#BC6EB9", "#A93154", "gray"),
                 symbol = 1,
                 symbols = "square",
                 marker = list(size = 4))
  
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                     yaxis = list(title = 'y'),
                                     zaxis = list(title = 'z')))
  return(fig)
}
