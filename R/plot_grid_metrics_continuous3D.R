plot_grid_metrics_continuous3D <- function(grid_metrics, 
                                           metric_colname) {
  
  ## Check input parameters
  if (!(is.character(metric_colname) && metric_colname %in% c("proportion", "entropy"))) {
    stop("`metric_colname` is not 'proportion' or 'entropy'.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  
  ## Color of each dot is related to its entropy
  pal <- colorRampPalette(hcl.colors(n = 5, palette = "Red-Blue", rev = TRUE))
  
  ## Add size column and for NA entropy values, make the size small
  grid_metrics$size <- ifelse(is.na(grid_metrics[[metric_colname]]), 3, 10)
  
  fig <- plot_ly(grid_metrics,
                 type = "scatter3d",
                 mode = 'markers',
                 x = ~x_coord,
                 y = ~y_coord,
                 z = ~z_coord,
                 color = as.formula(paste0('~', metric_colname)),
                 colors = pal(nrow(grid_metrics)),
                 marker = list(size = ~size),
                 symbol = 1,
                 symbols = "square")
  
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                     yaxis = list(title = 'y'),
                                     zaxis = list(title = 'z')))
  
  return(fig)
}
