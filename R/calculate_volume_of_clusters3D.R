#' @title Algorithm to determine the volume of cell clusters in 3D spatial data.
#'
#' @description This function finds the volume of cell clusters in a 3D 
#'     SpatialExperiment Object, where the cell clusters have already been 
#'     identified using an existing SPIAT-3D cell clustering algorithm function.
#' 
#' @param spe A SpatialExperiment object containing 3D spatial information for 
#'     the cells. Naming of spatial coordinates MUST be "Cell.X.Position", 
#'     "Cell.Y.Position", "Cell.Z.Position" for the x-coordinate, y-coordinate 
#'     and z-coordinate of each cell. It must also contain the cell clustering 
#'     information, obtained by passing the SpatialExperiment object through one 
#'     of the cell clustering algorithm functions in SPIAT-3D 
#'     (alpha_hull_clustering3D, grid_based_clustering3D, dbscan_clustering3D).
#' @param cluster_colname A string specifying the name of the column in the 
#'     `colData` slot of the SpatialExperiment object that contains the cell 
#'     clustering information. Should be 'alpha_hull_cluster', 'dbscan_cluster', 
#'     or 'grid_based_cluster'.
#'
#' @return A data frame containing the volume of each cluster, using different
#'     volume calculation methods, depending on the clustering method used.
#'
#' @examples
#' volume_of_clusters <- calculate_volume_of_clusters3D(
#'     spe = SPIAT-3D::simulated_spe_with_alpha_hull_clustering,
#'     cluster_colname = "alpha_hull_cluster"
#' )
#' 
#' @export

calculate_volume_of_clusters3D <- function(spe, 
                                           cluster_colname) {
  
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
