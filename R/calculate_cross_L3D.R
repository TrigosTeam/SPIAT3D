#' @title Calculate cross L-function on 3D spatial data.
#'
#' @description This function calculates the cross L-function on a 3D
#'     SpatialExperiment Object, which is the linearised transformation of the
#'     cross K-function.
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for
#'     the cells. Naming of spatial coordinates MUST be "Cell.X.Position",
#'     "Cell.Y.Position", "Cell.Z.Position" for the x-coordinate, y-coordinate
#'     and z-coordinate of each cell.
#' @param reference_cell_type A string specifying the reference cell type.
#' @param target_cell_types A character vector specifying the target cell types.
#' @param radius A positive numeric specifying the radius value.
#' @param feature_colname A string specifying the name of the column in the
#'     `colData` slot of the SpatialExperiment object that contains the cell
#'     type information.
#'
#' @return A data frame containing observed cross L-function values for each
#'     target cell type and expected cross L-function values.
#'
#' @examples
#' # Get simulated SpatialExperiment object to use as an example for analysis
#' simulated_spe <- readRDS(system.file("extdata", "simulated_spe.rds", package = "SPIAT3D"))
#'
#' result <- calculate_cross_L3D(
#'     spe = simulated_spe,
#'     reference_cell_type = "Tumour",
#'     target_cell_types = c("Tumour", "Immune"),
#'     radius = 30,
#'     feature_colname = "Cell.Type"
#' )
#'
#' @export

calculate_cross_L3D <- function(spe,
                                reference_cell_type,
                                target_cell_types,
                                radius,
                                feature_colname) {

  # Calculate cross K first
  result <- calculate_cross_K3D(spe = spe,
                                reference_cell_type = reference_cell_type,
                                target_cell_types = target_cell_types,
                                radius = radius,
                                feature_colname = feature_colname)

  # Linearise cross K to get cross L
  result[ , c("expected", target_cell_types)] <- (result[ , c("expected", target_cell_types)] / (4 * pi / 3)) ^ (1/3)

  return(result)
}
