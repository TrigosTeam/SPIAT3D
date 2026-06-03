#' @title Calculate cell proportion of clusters in 3D spatial data.
#'
#' @description This function finds the proportion of cells which make up cell
#'     clusters in a 3D SpatialExperiment Object, where the cell clusters have
#'     already been identified using an existing SPIAT-3D cell clustering
#'     algorithm function.
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
#' @param feature_colname A string specifying the name of the column in the
#'     `colData` slot of the SpatialExperiment object that contains the cell
#'     type information. Defaults to "Cell.Type".
#' @param plot_image A logical indicating whether to plot cell proportion of
#'     the cell clusters in bar plot form. Defaults to TRUE.
#'
#' @return A data frame containing the cell proportion of all cell types in
#'     each cell cluster.
#'
#' @importFrom ggplot2 ggplot aes geom_bar labs scale_x_discrete guides theme_bw theme guide_legend element_text
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
#' # Get cell proportions that make up dbscan clusters
#' cluster_cell_props <- calculate_cell_proportions_of_clusters3D(
#'     spe = dbscan_spe,
#'     cluster_colname = "dbscan_cluster",
#'     feature_colname = "Cell.Type",
#'     plot_image = T
#' )
#'
#' @export

calculate_cell_proportions_of_clusters3D <- function(spe,
                                                     cluster_colname,
                                                     feature_colname = "Cell.Type",
                                                     plot_image = T) {

  # Get number of clusters
  n_clusters <- max(spe[[cluster_colname]])

  ## Get different cell types found in the clusters (alphabetical for consistency)
  cell_types <- unique(spe[[feature_colname]][spe[[cluster_colname]] != 0])
  cell_types <- cell_types[order(cell_types)]

  ## For each cluster, determine the size and cell proportion of each cluster
  result <- data.frame(matrix(nrow = n_clusters, ncol = 2 + length(cell_types)))
  colnames(result) <- c("cluster_number", "n_cells", cell_types)
  result$cluster_number <- as.character(seq(n_clusters))

  for (i in seq(n_clusters)) {
    cells_in_cluster <- spe[[feature_colname]][spe[[cluster_colname]] == i]
    result[i, "n_cells"] <- length(cells_in_cluster)

    for (cell_type in cell_types) {
      result[i, cell_type] <- sum(cells_in_cluster == cell_type) / result[i, "n_cells"]
    }
  }

  ## Plot
  if (plot_image) {
    plot_result <- reshape2::melt(result, id.vars = c("cluster_number", "n_cells"))
    fig <- ggplot(plot_result, aes(cluster_number, value, fill = variable)) +
      geom_bar(stat = "identity") +
      labs(title = "Cell proportions of each cluster", x = "", y = "Cell proportion") +
      scale_x_discrete(labels = paste("cluster_", result$cluster_number, ", n = ", result$n_cells, sep = "")) +
      guides(fill = guide_legend(title="Cell type")) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))

    methods::show(fig)
  }

  return(result)
}
