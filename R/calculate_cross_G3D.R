#' @title Calculate cross G on 3D spatial data.
#'
#' @description This function calculates the cross G on a 3D SpatialExperiment
#'     Object. This metric finds the proportion of reference cells which have at
#'     least one interaction with a target cell as well as the expected cross G
#'     values if the 3D spatial data followed a complete spatial randomness
#'     pattern. This is calculated for a single radius value.
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for
#'     the cells. Naming of spatial coordinates MUST be "Cell.X.Position",
#'     "Cell.Y.Position", "Cell.Z.Position" for the x-coordinate, y-coordinate
#'     and z-coordinate of each cell.
#' @param reference_cell_type A string specifying the reference cell type.
#' @param target_cell_type A string specifying the target cell type.
#' @param radius A positive numeric specifying the radius value.
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
#' result <- calculate_cross_G3D(
#'     spe = simulated_spe,
#'     reference_cell_type = "Tumour",
#'     target_cell_type = "Immune",
#'     radius = 5,
#'     feature_colname = "Cell.Type"
#' )
#'
#' @export

calculate_cross_G3D <- function(spe,
                                reference_cell_type,
                                target_cell_type,
                                radius,
                                feature_colname = "Cell.Type") {

  ### Calculate the observed cross_G
  # Get the number of target cells in the radius around each reference cell
  neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spe,
                                                              reference_cell_type,
                                                              target_cell_type,
                                                              radius,
                                                              feature_colname)

  reference_target_interactions <- neighbourhood_counts_df[[target_cell_type]]

  # cross_G: essentially the proportion of reference cells with at least 1 target cell within the chosen radius.
  observed_cross_G <- sum(reference_target_interactions != 0) / length(reference_target_interactions)

  ### Calculate the expected cross_G
  # Get rough dimensions of the window the points are in
  spe_coords <- data.frame(SpatialExperiment::spatialCoords(spe))

  length <- round(max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- round(max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  height <- round(max(spe_coords$Cell.Z.Position) - min(spe_coords$Cell.Z.Position))

  # Get volume of the window the cells are in
  volume <- length * width * height

  # Get the number of target cells
  n_target_cells <- sum(spe[[feature_colname]] == target_cell_type)

  # Get target_cell_type intensity (density)
  target_cell_type_intensity <- n_target_cells / volume

  # Apply formula
  expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * (4 / 3) * pi * radius^3)

  result <- data.frame(observed_cross_G = observed_cross_G,
                       expected_cross_G = expected_cross_G)

  return(result)
}
