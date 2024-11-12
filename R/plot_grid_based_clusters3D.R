plot_grid_based_clusters3D <- function(spe_with_grid, 
                                       plot_cell_types = NULL,
                                       plot_colours = NULL,
                                       feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spe_with_grid) != "SpatialExperiment") {
    stop("`spe_with_grid` is not a SpatialExperiment object.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe_with_grid[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  
  ## If no cell types chosen, use all cell types found in data frame
  if (is.null(plot_cell_types)) plot_cell_types <- unique(spe_with_grid[[feature_colname]])
  
  ## If cell types have been chosen, check they are found in the spe object
  unknown_cell_types <- setdiff(plot_cell_types, spe_with_grid[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following plot_cell_types are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  
  ## If no colours inputted, use rainbow palette
  if (is.null(plot_colours)) plot_colours <- rainbow(length(plot_cell_types))
  
  ## User inputs mismatching cell types and colours
  if (length(plot_cell_types) != length(plot_colours)) stop("Length of plot_cell_types is not equal to length of plot_colours")
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe_with_grid), colData(spe_with_grid))
  
  ## Factor for feature column
  df[[feature_colname]] <- factor(df[[feature_colname]], levels = plot_cell_types)
  
  ## Add points to fig
  fig <- plot_ly() %>%
    add_trace(
      data = df,
      type = "scatter3d",
      mode = 'markers',
      x = ~Cell.X.Position,
      y = ~Cell.Y.Position,
      z = ~Cell.Z.Position,
      marker = list(size = 2),
      color = ~.data[[feature_colname]],
      colors = plot_colours
    ) %>% 
    layout(scene = list(xaxis = list(title = 'x'),
                        yaxis = list(title = 'y'),
                        zaxis = list(title = 'z')))
  
  # Get number of grid-based clusters
  n_grid_based_clusters <- length(spe_with_grid@metadata[["grid_prisms"]])
  
  faces <- data.frame(edge1 = c(1, 1, 1, 1, 1, 1, 8, 8, 8, 8, 8, 8),
                      edge2 = c(2, 5, 2, 3, 3, 5, 6, 4 ,7, 6, 7, 4),
                      edge3 = c(6, 6, 4, 4, 7, 7, 2, 2, 5, 5, 3, 3))
  grid_based_colours <- rainbow(n_grid_based_clusters)
  
  ## Add grid-based clusters to fig, one by one  
  for (i in seq(n_grid_based_clusters)) {
    
    grid_based_cluster <- spe_with_grid@metadata[["grid_prisms"]][[i]]
    
    for (j in seq(nrow(grid_based_cluster))) {
      
      x <- grid_based_cluster$x[j]
      y <- grid_based_cluster$y[j]
      z <- grid_based_cluster$z[j]
      l <- grid_based_cluster$l[j]
      w <- grid_based_cluster$w[j]
      h <- grid_based_cluster$h[j]
      vertices <- data.frame(x = c(x, x + l, x, x + l, x, x + l, x, x + l),
                             y = c(y, y, y + w, y + w, y, y, y + w, y + w),
                             z = c(z, z, z, z, z + h, z + h, z + h, z + h))
      
      fig <- fig %>%
        add_trace(
          type = 'mesh3d',
          x = vertices[, 1], 
          y = vertices[, 2], 
          z = vertices[, 3],
          i = faces[, 1] - 1, 
          j = faces[, 2] - 1, 
          k = faces[, 3] - 1,
          opacity = 0.2,
          facecolor = rep(grid_based_colours[i], 12) # Always 12 faces per grid prism
        )      
    }
  }
  
  return(fig)
}
