library(alphashape3d)

alpha_hull_clustering3D <- function(spe, 
                                    cell_types_of_interest, 
                                    alpha, 
                                    minimum_cells_in_cluster,
                                    feature_colname = "Cell.Type", 
                                    plot_image = T) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  ## Check cell types of interst are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  if (!(is.numeric(alpha) && length(alpha) == 1 && alpha > 0)) {
    stop("`alpha` is not a positive numeric.")
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
  
  ## Subset for the chosen cell_types_of_interest
  spe_subset <- spe[ , spe[[feature_colname]] %in% cell_types_of_interest]
  spe_subset_coords <- spatialCoords(spe_subset)
  
  ## Get the alpha hull
  alpha_hull <- ashape3d(as.matrix(spe_subset_coords), alpha = alpha)
  
  if (sum(alpha_hull$triang[, 9]) == 0) stop("alpha value is too small? No alpha hulls identified")
  
  ## Determine which alpha hull cluster each cell_type_of_interest belongs to
  alpha_hull_clusters <- components_ashape3d(alpha_hull)
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe), colData(spe))
  
  df_cell_types_of_interest <- df[df[[feature_colname]] %in% cell_types_of_interest, ]
  df_other_cell_types <- df[!(df[[feature_colname]] %in% cell_types_of_interest), ]
  
  df_cell_types_of_interest$alpha_hull_cluster <- alpha_hull_clusters
  df_other_cell_types$alpha_hull_cluster <- 0
  
  ## Ignore cell_types_of_interest which belong to an alpha hull cluster with less than minimum_cells_in_cluster
  alpha_hull_clusters_table <- table(alpha_hull_clusters)
  maximium_alpha_hull_cluster <- Position(function(x) x < minimum_cells_in_cluster, alpha_hull_clusters_table)
  maximium_alpha_hull_cluster <- as.numeric(names(alpha_hull_clusters_table[maximium_alpha_hull_cluster]))

  if (!is.na(maximium_alpha_hull_cluster) && maximium_alpha_hull_cluster != -1) {
    spe_subset_coords <- spe_subset_coords[alpha_hull_clusters >= 1 & alpha_hull_clusters < maximium_alpha_hull_cluster, ]
    
    df_cell_types_of_interest$alpha_hull_cluster <- ifelse(alpha_hull_clusters >= 1 & alpha_hull_clusters < maximium_alpha_hull_cluster, 
                                                           alpha_hull_clusters, 0)
  
    ## Get the alpha hull again...
    alpha_hull <- ashape3d(as.matrix(spe_subset_coords), alpha = alpha)
  }
  
  ## Convert data frame to spe object
  df <- rbind(df_cell_types_of_interest, df_other_cell_types)
  
  spe <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = spe@metadata)
  
  ## Get the information of the vertices and faces of the alpha hull (what 3 vertices make up each face triangle?)
  vertices <- alpha_hull$x
  faces <- alpha_hull$triang[alpha_hull$triang[, 9] == 2, c("tr1", "tr2", "tr3")]
  spe@metadata$alpha_hull <- list(vertices = vertices, faces = faces, ashape3d_object = alpha_hull)
  
  ## Plot
  if (plot_image) {
    fig <- plot_alpha_hull_clusters3D(spe, feature_colname = feature_colname)
    methods::show(fig)
  }
  
  return(spe)
}
