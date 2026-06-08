#' @title Calculate minimum distances to cell clusters in 3D spatial data.
#'
#' @description This function calculates the minimum distances between
#'     non-cluster cells to cluster cells in a  3D SpatialExperiment object,
#'     where the cell clusters have already been identified using an existing
#'     SPIAT3D cell clustering algorithm function.
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for
#'     the cells. Naming of spatial coordinates MUST be "Cell.X.Position",
#'     "Cell.Y.Position", "Cell.Z.Position" for the x-coordinate, y-coordinate
#'     and z-coordinate of each cell. It must also contain the cell clustering
#'     information, obtained by passing the SpatialExperiment object through one
#'     of the cell clustering algorithm functions in SPIAT3D
#'     (alpha_hull_clustering3D, grid_based_clustering3D, dbscan_clustering3D).
#' @param cell_types_inside_cluster A character vector specifying the cell types
#'     inside the cell clusters that are of interest. Minimum distances will be
#'     calculated TO these cell types.
#' @param cell_types_outside_cluster A character vector specifying the cell
#'     types outside the cell clusters that are of interest. Minimum distances
#'     will be calculated FROM these cell types.
#' @param cluster_colname A string specifying the name of the column in the
#'     `colData` slot of the SpatialExperiment object that contains the cell
#'     clustering information. Should be 'alpha_hull_cluster', 'dbscan_cluster',
#'     or 'grid_based_cluster'.
#' @param feature_colname A string specifying the name of the column in the
#'     `colData` slot of the SpatialExperiment object that contains the cell
#'     type information.
#' @param plot_image A logical indicating whether to plot violin plots of the
#'     minimum distances each non-cluster cell to each cluster. Defaults to
#'     TRUE.
#'
#' @return A data frame containing information about the minimum distances
#'     between FROM non-clusters TO cluster cells, for each cluster.
#'
#' @importFrom ggplot2 ggplot aes geom_violin facet_grid labeller labs stat_summary theme_bw theme element_blank element_text
#'
#' @examples
#' # Get simulated SpatialExperiment object to use as an example for analysis
#' simulated_spe <- readRDS(system.file("extdata", "simulated_spe.rds", package = "SPIAT3D"))
#'
#' # Get dbscan clusters
#' dbscan_spe <- dbscan_clustering3D(
#'     spe = simulated_spe,
#'     cell_types_of_interest = c("Tumour", "Immune"),
#'     radius = 8,
#'     minimum_cells_in_radius = 10,
#'     minimum_cells_in_cluster = 30,
#'     feature_colname = "Cell.Type",
#'     plot_image = TRUE
#' )
#'
#' # Get minimum distances to dbscan clusters
#' minimum_distances_to_clusters <- calculate_minimum_distances_to_clusters3D(
#'     spe = dbscan_spe,
#'     cell_types_inside_cluster = c("Tumour"),
#'     cell_types_outside_cluster = c("Tumour", "Immune"),
#'     cluster_colname = "dbscan_cluster",
#'     feature_colname = "Cell.Type",
#'     plot_image = TRUE
#' )
#'
#' @export

calculate_minimum_distances_to_clusters3D <- function(spe,
                                                      cell_types_inside_cluster,
                                                      cell_types_outside_cluster,
                                                      cluster_colname,
                                                      feature_colname,
                                                      plot_image = T) {

  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (!is.character(cell_types_inside_cluster)) {
    stop("`cell_types_inside_cluster` is not a character vector.")
  }
  if (!is.character(cell_types_outside_cluster)) {
    stop("`cell_types_outside_cluster` is not a character vector.")
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
  spe_coords <- SpatialExperiment::spatialCoords(spe)

  # Cells outside cluster have a cluster number of 0 (i.e. they are not in a cluster)
  spe_outside_cluster <- spe[ , spe[[cluster_colname]] == 0]

  cell_types_outside_cluster_coords <- list()
  for (cell_type in cell_types_outside_cluster) {
    cell_types_outside_cluster_coords[[cell_type]] <- SpatialExperiment::spatialCoords(spe_outside_cluster)[spe_outside_cluster[[feature_colname]] == cell_type, ]
  }

  ## For each cluster, determine the minimum distance of each outside_cell_type
  result <- vector()

  # Get number of clusters
  n_clusters <- max(spe[[cluster_colname]])

  for (i in seq(n_clusters)) {
    # Get coordinates of cells inside current cluster
    cluster_coords <- spe_coords[spe[[cluster_colname]] == i & spe[[feature_colname]] %in% cell_types_inside_cluster, ]
    # Get cell types of cells inside current cluster
    cluster_cell_types <- spe[["Cell.Type"]][spe[[cluster_colname]] == i & spe[[feature_colname]] %in% cell_types_inside_cluster]
    # Get ids of cells inside current cluster
    cluster_cell_ids <- spe[["Cell.ID"]][spe[[cluster_colname]] == i & spe[[feature_colname]] %in% cell_types_inside_cluster]

    for (outside_cell_type in cell_types_outside_cluster) {
      # Get coordinates of cells outside of any cluster
      curr_cell_type_coords <- cell_types_outside_cluster_coords[[outside_cell_type]]

      # Get minimum distances of each cell outside of any cluster to current cluster
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
        stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1), colour = "red")

      methods::show(fig)
    }

  }
  return(result)
}
