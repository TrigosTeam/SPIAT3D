#' @title Calculate cells in neighbourhood (gradient version) on 3D spatial
#'     data.
#'
#' @description This function calculates the cells in neighbourhood (gradient
#'     version) on a 3D SpatialExperiment Object. This metric finds the
#'     proportion of target cells around each reference cell, for each target
#'     cell type, and calculates the average proportion across all radii.
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for
#'     the cells. Naming of spatial coordinates MUST be "Cell.X.Position",
#'     "Cell.Y.Position", "Cell.Z.Position" for the x-coordinate, y-coordinate
#'     and z-coordinate of each cell.
#' @param reference_cell_type A string specifying the reference cell type.
#' @param target_cell_types A character vector specifying the target cell types.
#' @param radii A positive, ascending numeric vector specifying the set of
#'     radius values used to calculate cells in neighbourhood over a gradient.
#' @param feature_colname A string specifying the name of the column in the
#'     `colData` slot of the SpatialExperiment object that contains the cell
#'     type information. Defaults to "Cell.Type".
#' @param plot_image A logical indicating whether to plot cells in neighbourhood
#'     gradient as a line graph. Defaults to TRUE.
#'
#' @return A data frame containing the cells in neighbourhood values for each
#'     target cell type (columns) across each radii (rows).
#'
#' @examples
#' # Get simulated SpatialExperiment object to use as an example for analysis
#' simulated_spe <- readRDS(system.file("extdata", "simulated_spe.rds", package = "SPIAT3D")
#'
#' result <- calculate_cells_in_neighbourhood_gradient3D(
#'     spe = simulated_spe,
#'     reference_cell_type = "Tumour",
#'     target_cell_types = c("Tumour", "Immune"),
#'     radii = seq(20, 100, 10),
#'     feature_colname = "Cell.Type",
#'     plot_image = TRUE
#' )
#'
#' @export

calculate_cells_in_neighbourhood_gradient3D <- function(spe,
                                                        reference_cell_type,
                                                        target_cell_types,
                                                        radii,
                                                        feature_colname = "Cell.Type",
                                                        plot_image = TRUE) {

  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }

  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types

  for (i in seq(length(radii))) {
    cell_proportions_neighbourhood_proportions_df <- calculate_cells_in_neighbourhood3D(spe,
                                                                                        reference_cell_type,
                                                                                        target_cell_types,
                                                                                        radii[i],
                                                                                        feature_colname)

    if (is.null(cell_proportions_neighbourhood_proportions_df)) return(NULL)

    result[i, ] <- apply(cell_proportions_neighbourhood_proportions_df[ , paste(target_cell_types, "_prop", sep = "")], 2, mean, na.rm = T)
  }

  # Add a radius column to the result
  result$radius <- radii

  # Plot
  if (plot_image) {
    fig <- plot_cells_in_neighbourhood_gradient3D(result, reference_cell_type)
    methods::show(fig)
  }

  return(result)
}
