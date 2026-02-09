#' @title Algorithm to determine cells which form the border of cell clusters in 3D spatial data.
#'
#' @description This function finds which cells form the border of cell clusters on in a 3D SpatialExperiment Object, 
#'     where the cell clusters have already been identified using an existing SPIAT-3D cell clustering algorithm function.
#' 
#' @param spe A SpatialExperiment object containing 3D spatial information for the cells. 
#'     Naming of spatial coordinates MUST be "Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position" 
#'     for the x-coordinate, y-coordinate and z-coordinate of each cell. It must also contain the
#'     cell clustering information, obtained by passing the SpatialExperiment object through one of the cell
#'     clustering algorithm functions in SPIAT-3D (alpha_hull_clustering3D, grid_based_clustering3D, dbscan_clustering3D)
#' @param cluster_colname A string specifying the name of the column in the `colData` slot of the SpatialExperiment object 
#'     that contains the cell clustering information. Should be 'alpha_hull_cluster', 'dbscan_cluster', or 'grid_based_cluster'
#' @param feature_colname A string specifying the name of the column in the `colData` slot of the SpatialExperiment
#'    object that contains the cell type information. Defaults to "Cell.Type"
#' @param plot_image A logical indicating whether to plot cell clusters with bordering cells. Defaults to TRUE.
#'
#' @return The same 3D SpatialExperiment object used as input for spe, with an added column in the `colData` slot
#'     called 'cluster_border' to specify which cells form the border of the clusters, 
#'     and which cells are outside or inside the cluster
#'
#' @examples
#' spe_with_border_of_clusters <- calculate_border_of_clusters3D(
#'     spe = SPIAT-3D::simulated_spe_with_alpha_hull_clustering,
#'     cluster_colname = "alpha_hull_cluster",
#'     feature_colname = "Cell.Type",
#'     plot_image = T
#' )
#' 
#' @export


calculate_border_of_clusters3D <- function(spe, 
                                           radius,
                                           cluster_colname, 
                                           feature_colname = "Cell.Type", 
                                           plot_image = T) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (!(is.numeric(radius) && length(radius) == 1 && radius > 0)) {
    stop("`radius` is not a positive numeric.")
  }
  if (!is.character(cluster_colname)) {
    stop("`cluster_colname` is not a character. This should be 'alpha_hull_cluster', 'dbscan_cluster', or 'grid_based_cluster', depending on the chosen method.")
  }
  if (is.null(spe[[cluster_colname]])) {
    stop(paste("No column called", cluster_colname, "found in spe object."))
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  ## Get spatial coords of spe
  spe_coords <- data.frame(spatialCoords(spe))
  
  ## Get coords of non-cluster cells
  non_cluster_coords <- spe_coords[spe[[cluster_colname]] == 0, ]
  
  # New column for spe object: 'cluster_border'. Default is 'outside'
  spe$cluster_border <- "outside"
  
  # Label cells part of a cluster (e.g. 'cluster1')
  spe$cluster_border[spe[[cluster_colname]] != 0] <- paste("inside_C", spe[[cluster_colname]][spe[[cluster_colname]] != 0], sep = "")
  
  ## Iterate for each cluster
  n_clusters <- max(spe[[cluster_colname]])
  
  for (i in seq_len(n_clusters)) {
    
    ## Subset for cells in the current cluster of interest
    cluster_coords <- spe_coords[spe[[cluster_colname]] == i, ]
    
    # For each cell in the current cluster, check how many other cells in the cluster are in its radius
    cluster_to_cluster_interactions <- dbscan::frNN(cluster_coords, radius)
    
    # Determine the median minimum number of cluster cells found in the radius of cluster cell. Use this as the threshold for non-cluster cells.
    non_cluster_threshold <- quantile(unlist(lapply(cluster_to_cluster_interactions$dist, length)), 0.5)
    
    # For each non-cluster cell, check how many cluster cells are in its radius.
    non_cluster_to_cluster_interactions <- dbscan::frNN(cluster_coords, radius, non_cluster_coords)
    
    # If number of cluster cells found in the radius of non-cluster cells is greater than threshold, non-cluster cell has probably infiltrated cluster too
    n_cluster_cells_in_non_cluster_cell_radius <- unlist(lapply(non_cluster_to_cluster_interactions$id, length))
    
    spe$cluster_border[as.numeric(names(non_cluster_to_cluster_interactions$id)[n_cluster_cells_in_non_cluster_cell_radius > non_cluster_threshold])] <- paste("infiltrated_C", i, sep = "")
    
    # If number of cluster cells found in the radius of non-cluster cells is less than threshold, but greater than 0, non-cluster cell is probably on the border
    spe$cluster_border[as.numeric(names(non_cluster_to_cluster_interactions$id)[n_cluster_cells_in_non_cluster_cell_radius > 0 & n_cluster_cells_in_non_cluster_cell_radius < non_cluster_threshold])] <- paste("border_C", i, sep = "")
  }
  
  ## Plot
  if (plot_image) {
    fig <- plot_cells3D(spe, feature_colname = "cluster_border")
    methods::show(fig)
  }
  
  return(spe)
}
