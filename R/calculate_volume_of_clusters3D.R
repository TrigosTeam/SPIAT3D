calculate_volume_of_clusters3D <- function(spe, cluster_colname) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!is.character(cluster_colname)) {
    stop("`cluster_colname` is not a character. This should be 'alpha_hull_cluster', 'dbscan_cluster', or 'grid_based_cluster', depending on the chosen method.")
  }
  if (is.null(spe[[cluster_colname]])) {
    stop(paste("No column called", cluster_colname, "found in spe object."))
  }
  
  # Get number of clusters
  n_clusters <- max(spe[[cluster_colname]])
  
  ### 1. Estimate volume of each cluster by density of the window. ------------
  
  ## For each cluster, determine the number of cells in each cluster of each cluster
  result <- data.frame(matrix(nrow = n_clusters, ncol = 2))
  colnames(result) <- c("cluster_number", "n_cells")
  
  for (i in seq(n_clusters)) {
    result[i, "n_cells"] <- sum(spe[[cluster_colname]] == i)
  }
  result$cluster_number <- as.character(seq(n_clusters))
  
  ## Assume window is a rectangular prism
  spe_coords <- data.frame(spatialCoords(spe))
  
  length <- round(max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- round(max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  height <- round(max(spe_coords$Cell.Z.Position) - min(spe_coords$Cell.Z.Position))
  
  window_volume <- length * width * height
  
  result$volume_by_density <- (result$n_cells / ncol(spe)) * window_volume
  
  
  ### 2. If cluster_colname == "alpha_hull_cluster", use the volume method found in the alphashape3d package
  if (cluster_colname == "alpha_hull_cluster") {
    result$volume_by_alpha_hull <- volume_ashape3d(spe@metadata$alpha_hull$ashape3d_object, byComponents = T)
  }
  
  
  ### 3. If cluster_colname == "grid_based_cluster", sum the volume of each grid prism to get volume of each cluster
  if (cluster_colname == "grid_based_cluster") {
    result$volume_by_grid <- 0
    i <- 1
    for (grid_cluster in spe@metadata$grid_prisms) {
      result[i, "volume_by_grid"] <- sum(grid_cluster$l * grid_cluster$w * grid_cluster$h)
      i <- i + 1
    }
  }
  
  return(result)
}
