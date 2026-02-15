#' @title Calculate the background entropy of 3D spatial data.
#'
#' @description This function calculates the background entropy of a 3D 
#'     SpatialExperiment object across a subset of cell types chosen by the
#'     user.
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for 
#'     the cells. Naming of spatial coordinates MUST be "Cell.X.Position", 
#'     "Cell.Y.Position", "Cell.Z.Position" for the x-coordinate, y-coordinate 
#'     and z-coordinate of each cell.
#' @param cell_types_of_interest A character vector specifying the cell types of 
#'     interest. If NULL, all cell types in the `feature_colname` column will be 
#'     considered.
#' @param feature_colname A string specifying the name of the column in the 
#'     `colData` slot of the SpatialExperiment object that contains the cell 
#'     type information. Defaults to "Cell.Type".
#'
#' @return A numeric representing the background entropy.
#'
#' @examples
#' background_entropy <- calculate_entropy_background3D(
#'     spe = SPIAT-3D::simulated_spe,
#'     cell_types_of_interest = c("Tumour", "Immune"),
#'     feature_colname = "Cell.Type"
#' )
#' 
#' @export

calculate_entropy_background3D <- function(spe,
                                           cell_types_of_interest, 
                                           feature_colname = "Cell.Type") {
  
  # NULL case: entropy is undefined
  if (is.null(cell_types_of_interest)) return(NA)
  
  # One cell type case: entropy is 0
  if (is.character(cell_types_of_interest) && length(cell_types_of_interest) == 1) return(0)
  
  cell_proportions_data <- calculate_cell_proportions3D(spe, cell_types_of_interest, feature_colname, FALSE)
  
  # Calculate entropy of the entire image
  entropy <- -1 * sum(cell_proportions_data$proportion * log(cell_proportions_data$proportion, length(cell_proportions_data$proportion)))
  
  return(entropy) 
}
