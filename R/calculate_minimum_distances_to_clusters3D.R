calculate_minimum_distances_to_clusters3D <- function(spe, 
                                                      cell_types_inside_cluster, 
                                                      cell_types_outside_cluster, 
                                                      cluster_colname, 
                                                      feature_colname = "Cell.Type", 
                                                      plot_image = T) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!is.character(cell_types_inside_cluster)) {
    stop("`cell_types_inside_cluster` is not a character vector.")
  }
  if (!is.character(cell_types_outside_cluster)) {
    stop("`cell_types_outside_cluster` is not a character vector.")
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
  
  ## Add Cell.ID column
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.Id column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }
  
  ## For each cell type outside clusters, get their set of coords. These exclude cell types in clusters
  spe_coords <- spatialCoords(spe)
  
  # Cells outside cluster have a cluster number of 0 (i.e. they are not in a cluster)
  spe_outside_cluster <- spe[ , spe[[cluster_colname]] == 0]
  
  cell_types_outside_cluster_coords <- list()
  for (cell_type in cell_types_outside_cluster) {
    cell_types_outside_cluster_coords[[cell_type]] <- spatialCoords(spe_outside_cluster)[spe_outside_cluster[[feature_colname]] == cell_type, ]
  }
  
  ## For each cluster, determine the minimum distance of each outside_cell_type  
  result <- vector()
  
  # Get number of clusters
  n_clusters <- max(spe[[cluster_colname]])
  
  for (i in seq(n_clusters)) {
    cluster_coords <- spe_coords[spe[[cluster_colname]] == i & spe[[feature_colname]] %in% cell_types_inside_cluster, ]
    cluster_cell_types <- spe[["Cell.Type"]][spe[[cluster_colname]] == i & spe[[feature_colname]] %in% cell_types_inside_cluster]
    cluster_cell_ids <- spe[["Cell.ID"]][spe[[cluster_colname]] == i & spe[[feature_colname]] %in% cell_types_inside_cluster]
    
    for (outside_cell_type in cell_types_outside_cluster) {
      curr_cell_type_coords <- cell_types_outside_cluster_coords[[outside_cell_type]]
      
      all_closest <- RANN::nn2(data = cluster_coords, 
                               query = curr_cell_type_coords, 
                               k = 1) 
      
      local_dist_mins <- data.frame(
        cluster_number = i,
        outside_cell_id = as.character(spe_outside_cluster$Cell.ID[spe_outside_cluster[["Cell.Type"]] == outside_cell_type]),
        outside_cell_type = outside_cell_type,
        inside_cell_id = cluster_cell_ids[c(all_closest$nn.idx)],
        inside_cell_type = cluster_cell_types[c(all_closest$nn.idx)],
        distance = all_closest$nn.dists
      )
      ## Remove any 0 distance rows
      local_dist_mins <- local_dist_mins[local_dist_mins$distance != 0, ]
      result <- rbind(result, local_dist_mins)
    }
    
    
    ## Plot
    if (plot_image) {
      
      cluster_number_labs <- paste("cluster_", seq(n_clusters), sep = "")
      names(cluster_number_labs) <- seq(n_clusters)
      
      fig <- ggplot(result, aes(x = outside_cell_type, y = distance, fill = outside_cell_type)) + 
        geom_violin() +
        facet_grid(cluster_number~., scales="free_x", labeller = labeller(cluster_number = cluster_number_labs)) +
        theme_bw() +
        theme(axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5), legend.position = "none") +
        labs(title="Minimum cell distances to clusters", x = "Cell type", y = "Distance") +
        stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
      
      methods::show(fig)
    }
    
  }
  return(result)
}
