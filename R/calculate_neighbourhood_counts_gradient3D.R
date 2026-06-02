#' @title Calculate neighbourhood counts gradient on 3D spatial data.
#'
#' @description This function calculates the neighbourhood counts gradient on a
#'     3D SpatialExperiment Object. This metric finds the number of target cells
#'     around each reference cell, for each target cell type, and calculates the
#'     average number across all radii.
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for
#'     the cells. Naming of spatial coordinates MUST be "Cell.X.Position",
#'     "Cell.Y.Position", "Cell.Z.Position" for the x-coordinate, y-coordinate
#'     and z-coordinate of each cell.
#' @param reference_cell_type A string specifying the reference cell type.
#' @param target_cell_types A character vector specifying the target cell types.
#' @param radii A positive, ascending numeric vector specifying the set of
#'     radius values used.
#' @param feature_colname A string specifying the name of the column in the
#'     `colData` slot of the SpatialExperiment object that contains the cell
#'     type information. Defaults to "Cell.Type".
#' @param plot_image A logical indicating whether to plot neighbourhood counts
#'     gradient as a line graph. Defaults to TRUE.
#'
#' @return A data frame containing the neighbourhood counts values for each
#'     target cell type (columns) across each radii (rows).
#'
#' @examples
#' # Get simulated SpatialExperiment object to use as an example for analysis
#' simulated_spe <- readRDS(system.file("extdata", "simulated_spe.rds", package = "SPIAT3D"))
#'
#' result <- calculate_neighbourhood_counts_gradient3D(
#'     spe = simulated_spe,
#'     reference_cell_type = "Tumour",
#'     target_cell_types = c("Tumour", "Immune"),
#'     radii = seq(20, 100, 10),
#'     feature_colname = "Cell.Type",
#'     plot_image = TRUE
#' )
#'
#' @export

calculate_neighbourhood_counts_gradient3D <- function(spe,
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
    neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spe,
                                                                reference_cell_type,
                                                                target_cell_types,
                                                                radii[i],
                                                                feature_colname,
                                                                FALSE,
                                                                FALSE)

    if (is.null(neighbourhood_counts_df)) return(NULL)

    neighbourhood_counts_df$ref_cell_id <- NULL
    result[i, ] <- apply(neighbourhood_counts_df, 2, mean)
  }
  # Add a radius column to the result
  result$radius <- radii

  if (plot_image) {
    fig <- plot_neighbourhood_counts_gradient3D(result, reference_cell_type)
    methods::show(fig)
  }

  return(result)
}
