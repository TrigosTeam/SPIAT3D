#' @title Calculate co-occurrence gradient on 3D spatial data.
#'
#' @description This function calculates the co-occurrence gradient on a 3D 
#'     SpatialExperiment Object. This metric finds the proportion of target
#'     cells around reference cells relative to the proportion of target cells
#'     in the SpatialExperiment Object, for each target cell type, and across 
#'     all radii values.
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
#' @param plot_image A logical indicating whether to plot co-occurrence gradient
#'     as a line graph. Defaults to TRUE.
#'
#' @return A data frame containing the co-occurrence values for each target cell 
#'     type (columns) across each radii (rows).
#'
#' @examples
#' result <- calculate_co_occurrence_gradient3D(
#'     spe = SPIAT-3D::simulated_spe,
#'     reference = "Tumour",
#'     target = c("Tumour", "Immune"),
#'     radii = seq(20, 100, 10),
#'     feature_colname = "Cell.Type",
#'     plot_image = TRUE
#' )
#' 
#' @export

calculate_co_occurrence_gradient3D <- function(spe, 
                                               reference_cell_type, 
                                               target_cell_types, 
                                               radii, 
                                               feature_colname = "Cell.Type",
                                               plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types) + 1))
  colnames(result) <- c("reference", target_cell_types)
  
  for (i in seq(length(radii))) {
    co_occurrence_df <- calculate_co_occurrence3D(spe,
                                                  reference_cell_type,
                                                  target_cell_types,
                                                  radii[i],
                                                  feature_colname)
    
    result[i, ] <- co_occurrence_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_co_occurrence_gradient3D(result)
    methods::show(fig)
  }
  
  return(result)
}
