#' @title Calculate cross K gradient on 3D spatial data.
#'
#' @description This function calculates the cross K gradient on a 3D 
#'     SpatialExperiment Object. See paper on theory behind cross K gradient (I
#'     ain't explaining it here).
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
#' @param plot_image A logical indicating whether to plot cross K gradient as a 
#'     line graph showing both the observed cross K values for each cell type 
#'     and expected cross K values. Defaults to TRUE.
#'
#' @return A data frame containing observed cross K values for each target cell 
#'     type and expected cross K values (columns), all across each radii (rows).
#'
#' @examples
#' result <- calculate_cross_K_gradient3D(
#'     spe = SPIAT-3D::simulated_spe,
#'     reference_cell_type = "Tumour",
#'     target_cell_types = c("Tumour", "Immune"),
#'     radii = seq(20, 100, 10),
#'     feature_colname = "Cell.Type",
#'     plot_image = TRUE
#' )
#' 
#' @export

calculate_cross_K_gradient3D <- function(spe, 
                                         reference_cell_type, 
                                         target_cell_types, 
                                         radii, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  # Check if radii input is valid
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  # Set up output
  result <- data.frame(matrix(nrow = length(radii), ncol = 2 + length(target_cell_types)))
  colnames(result) <- c("reference", "expected", target_cell_types)
  
  # Iterate through each radii
  for (i in seq(length(radii))) {
    cross_K_df <- calculate_cross_K3D(spe,
                                      reference_cell_type,
                                      target_cell_types,
                                      radii[i],
                                      feature_colname)
    
    result[i, ] <- cross_K_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  # Plot
  if (plot_image) {
    fig <- plot_cross_K_gradient3D(result)
    methods::show(fig)
  }
  
  return(result)
}
