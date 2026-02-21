#' @title Find cell clusters in 3D spatial data using dbscan clustering 
#'     algorithm.
#'
#' @description This function finds cell clusters in a 3D SpatialExperiment 
#'     object using the dbscan clustering algorithm. 
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for 
#'     the cells. Naming of spatial coordinates MUST be "Cell.X.Position", 
#'     "Cell.Y.Position", "Cell.Z.Position" for the x-coordinate, y-coordinate 
#'     and z-coordinate of each cell.
#' @param cell_types_of_interest A character vector specifying the cell types of 
#'     interest.
#' @param radius A positive numeric. Spheres of specified radius are drawn
#'     around each cell type of interest and the number of other cell types of
#'     interest are counted.
#' @param minimum_cells_in_radius A positive numeric. If the number of cells 
#'     types of interest within the sphere of another cell type of interest 
#'     surpasses this specified value, they form a cluster.
#' @param minimum_cells_in_cluster A positive numeric. Clusters identified with
#'     dbscan which have less than this specified value are relabelled as not a
#'     cluster.
#' @param feature_colname A string specifying the name of the column in the 
#'     `colData` slot of the SpatialExperiment object that contains the cell 
#'     type information. Defaults to "Cell.Type"
#' @param plot_image A logical indicating whether to plot 3D spatial data with 
#'     alpha hull clusters. Defaults to TRUE.
#'
#' @return The same 3D SpatialExperiment object used as input for spe, with an 
#'     added column in the `colData` slot to specify which dbscan cluster each 
#'     cell belongs to.
#'
#' @examples
#' dbscan_spe <- dbscan_clustering3D(
#'     spe = SPIAT-3D::simulated_spe,
#'     cell_types_of_interest = c("Tumour", "Immune"),
#'     radius = 30,
#'     minimum_cells_in_radius = 10,
#'     minimum_cells_in_cluster = 30,
#'     feature_colname = "Cell.Type",
#'     plot_image = TRUE
#' )
#' 
#' @export

dbscan_clustering3D <- function(spe,
                                cell_types_of_interest,
                                radius,
                                minimum_cells_in_radius,
                                minimum_cells_in_cluster,
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
  ## Check cell types of interst are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  if (!(is.numeric(radius) && length(radius) == 1 && radius > 0)) {
    stop("`radius` is not a positive numeric.")
  }
  if (!(is.integer(minimum_cells_in_radius) && length(minimum_cells_in_radius) == 1 || (is.numeric(minimum_cells_in_radius) && length(minimum_cells_in_radius) == 1 && minimum_cells_in_radius > 0 && minimum_cells_in_radius%%1 == 0))) {
    stop("`minimum_cells_in_radius` is not a positive integer.")
  }
  if (!(is.integer(minimum_cells_in_cluster) && length(minimum_cells_in_cluster) == 1 || (is.numeric(minimum_cells_in_cluster) && length(minimum_cells_in_cluster) == 1 && minimum_cells_in_cluster > 0 && minimum_cells_in_cluster%%1 == 0))) {
    stop("`minimum_cells_in_cluster` is not a positive integer.")
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
  
  # Focus on cell types of interest
  spe_subset <- spe[ , spe[[feature_colname]] %in% cell_types_of_interest]
  spe_subset_coords <- spatialCoords(spe_subset)
  
  db <- dbscan::dbscan(spe_subset_coords, eps = radius, minPts = minimum_cells_in_radius, borderPoints = F)
  n_clusters <- max(db$cluster)
  if (n_clusters == 0) {
    stop("No clusters identified. Consider increasing `radius` and/or decreasing `minimum_cells_in_radius`.")
  }
  
  ## Cell types of interest have a 'cluster' value of 0 if they are noise, 1 if they belong to cluster 1, ...
  ## Check if number of cells in cluster 1, cluster 2, ... is larger than minimum_cells_in_cluster, if they don't, these cells are also assigned a value of 0.
  for (i in seq_len(n_clusters)) {
    if (sum(db$cluster == i) < minimum_cells_in_cluster) {
      db$cluster[db$cluster == i] <- 0
      
      ## Re-number the clusters
      db$cluster[db$cluster > i] <- db$cluster[db$cluster > i] - 1
    }
  }
  n_clusters <- max(db$cluster)
  if (n_clusters == 0) {
    stop("All clusters identified do not meet the `minimum_cells_in_cluster` threshold. Consider lowering the `minimum_cells_in_cluster` parameter.")
  }
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe), colData(spe))
  
  df_cell_types_of_interest <- df[df[[feature_colname]] %in% cell_types_of_interest, ]
  df_other_cell_types <- df[!(df[[feature_colname]] %in% cell_types_of_interest), ]
  
  df_cell_types_of_interest$dbscan_cluster <- db$cluster
  df_other_cell_types$dbscan_cluster <- 0
  
  ## Convert data frame to spe object
  df <- rbind(df_cell_types_of_interest, df_other_cell_types)
  
  spe <- SpatialExperiment(
    assay = matrix(data = NA, nrow = nrow(df), ncol = nrow(df)),
    colData = df,
    spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position", "Cell.Z.Position"),
    metadata = spe@metadata)
  
  ## Plot
  if (plot_image) {
    df$dbscan_cluster <- ifelse(df$dbscan_cluster == 0, "non_cluster", paste("cluster_", df$dbscan_cluster, sep = ""))
    
    fig <- plot_ly(df,
                   type = "scatter3d",
                   mode = 'markers',
                   x = ~Cell.X.Position,
                   y = ~Cell.Y.Position,
                   z = ~Cell.Z.Position,
                   color = ~dbscan_cluster,
                   colors = rainbow(length(unique(df$dbscan_cluster))),
                   marker = list(size = 2)) %>% 
      layout(scene = list(xaxis = list(title = 'x'),
                          yaxis = list(title = 'y'),
                          zaxis = list(title = 'z')))
    
    methods::show(fig)
  }
  
  return(spe)
}
