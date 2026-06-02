#' @title Calculate the centre of clusters in 3D spatial data.
#'
#' @description This function finds the centre of cell clusters in a 3D
#'     SpatialExperiment Object, where the cell clusters have already been
#'     identified using an existing SPIAT-3D cell clustering algorithm function.
#'     We assume the cell clusters have uniform density and that the centre of
#'     each cell cluster is defined by its centre of mass, which can be
#'     estimated by taking the average of the x, y, and z coordinates of
#'     all cells within a cluster.
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
#' @return A data frame containing the x, y, z coordinates of the centre of each
#'     cell cluster.
#'
#' @examples
#' # Get simulated SpatialExperiment object to use as an example for analysis
#' simulated_spe <- readRDS(system.file("extdata", "simulated_spe.rds", package = "SPIAT3D")
#'
#' # Get dbscan clusters
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
#' # Get centre of dbscan clusters
#' cluster_centres <- calculate_center_of_clusters3D(
#'     spe = dbscan_spe
#'     cluster_colname = "dbscan_cluster"
#' )
#'
#' @export

calculate_center_of_clusters3D <- function(spe,
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

  # Get spe coords
  spe_coords <- spatialCoords(spe)

  ## For each cluster, determine the number of cells in each cluster of each cluster
  result <- data.frame(matrix(nrow = n_clusters, ncol = 4))
  colnames(result) <- c("cluster_number", "Centre.X.Position", "Centre.Y.Position", "Centre.Z.Position")

  result$cluster_number <- as.character(seq(n_clusters))
  for (i in seq(n_clusters)) {
    spe_cluster_coords <- spe_coords[spe[[cluster_colname]] == i, ]
    result[i, c("Centre.X.Position", "Centre.Y.Position", "Centre.Z.Position")] <-
      apply(spe_cluster_coords, 2, mean)
  }

  return(result)
}
