#' @title Calculate cross G gradient on 3D spatial data.
#'
#' @description This function calculates the cross G gradient on a 3D
#'     SpatialExperiment Object. This metric finds the proportion of reference
#'     cells which have at least one interaction with a target cell as well as
#'     the expected cross G values if the 3D spatial data followed a complete
#'     spatial randomness pattern. This is calculated across all radii.
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for
#'     the cells. Naming of spatial coordinates MUST be "Cell.X.Position",
#'     "Cell.Y.Position", "Cell.Z.Position" for the x-coordinate, y-coordinate
#'     and z-coordinate of each cell.
#' @param reference_cell_type A string specifying the reference cell type.
#' @param target_cell_type A string specifying the target cell type.
#' @param radii A positive, ascending numeric vector specifying the set of
#'     radius values used.
#' @param feature_colname A string specifying the name of the column in the
#'     `colData` slot of the SpatialExperiment object that contains the cell
#'     type information. Defaults to "Cell.Type".
#' @param plot_image A logical indicating whether to plot cross G gradient as a
#'     line graph showing both the observed and expected cross G values.
#'     Defaults to TRUE.
#'
#' @return A data frame containing observed and expected cross G values across
#'     each radii (rows).
#'
#' @examples
#' # Get simulated SpatialExperiment object to use as an example for analysis
#' simulated_spe <- readRDS(system.file("extdata", "simulated_spe.rds", package = "SPIAT3D"))
#'
#' result <- calculate_cross_G_gradient3D(
#'     spe = simulated_spe,
#'     reference_cell_type = "Tumour",
#'     target_cell_type = "Immune",
#'     radii = seq(20, 100, 10),
#'     feature_colname = "Cell.Type",
#'     plot_image = TRUE
#' )
#'
#' @export

calculate_cross_G_gradient3D <- function(spe,
                                         reference_cell_type,
                                         target_cell_type,
                                         radii,
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {

  # Check if radii input is valid
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }

  # Set up output
  result <- data.frame(matrix(nrow = length(radii), ncol = 2))
  colnames(result) <- c("observed_cross_G",
                        "expected_cross_G")

  # Iterate through each radii
  for (i in seq(length(radii))) {
    cross_G_df <- calculate_cross_G3D(spe,
                                      reference_cell_type,
                                      target_cell_type,
                                      radii[i],
                                      feature_colname)

    result[i, ] <- cross_G_df
  }

  # Add a radius column to the result
  result$radius <- radii

  # Plot
  if (plot_image) {
    fig <- plot_cross_G_gradient3D(result, reference_cell_type, target_cell_type)
    methods::show(fig)
  }

  return(result)
}
