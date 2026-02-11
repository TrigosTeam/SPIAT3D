#' @title Calculate cells in neighbourhood on 3D spatial data.
#'
#' @description This function calculates the cells in neighbourhood on a 3D 
#'     SpatialExperiment Object. This metric finds the proportion of target 
#'     cells around each reference cell, for each target cell type.
#' 
#' @param spe A SpatialExperiment object containing 3D spatial information for 
#'     the cells. Naming of spatial coordinates MUST be "Cell.X.Position", 
#'     "Cell.Y.Position", "Cell.Z.Position" for the x-coordinate, y-coordinate 
#'     and z-coordinate of each cell.
#' @param reference_cell_type A string specifying the reference cell type.
#' @param target_cell_types A character vector specifying the target cell types.
#' @param radius A positive numeric specifying the radius values.
#' @param feature_colname A string specifying the name of the column in the 
#'     `colData` slot of the SpatialExperiment object that contains the cell 
#'     type information. Defaults to "Cell.Type".
#'
#' @return A data frame containing the cells in neighbourhood values for each
#'     reference cell (rows) and for each target cell type (columns).
#'
#' @examples
#' result <- calculate_cells_in_neighbourhood3D(
#'     spe = SPIAT-3D::simulated_spe,
#'     reference = "Tumour",
#'     target = c("Tumour", "Immune"),
#'     radius = 30,
#'     feature_colname = "Cell.Type"
#' )
#' 
#' @export

calculate_cells_in_neighbourhood3D <- function(spe, 
                                               reference_cell_type, 
                                               target_cell_types, 
                                               radius, 
                                               feature_colname = "Cell.Type") {
  
  ## Get cells in neighbourhood df
  neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spe,
                                                              reference_cell_type,
                                                              c(reference_cell_type, target_cell_types),
                                                              radius,
                                                              feature_colname,
                                                              FALSE,
                                                              FALSE)
  
  if (is.null(neighbourhood_counts_df)) return(NULL)
  
  neighbourhood_counts_df[ , paste(target_cell_types, "_prop", sep = "")] <- 
    neighbourhood_counts_df[ , target_cell_types] / (neighbourhood_counts_df[ , target_cell_types] + neighbourhood_counts_df[ , reference_cell_type])
  
  # If reference cell type is in target cell types, proportion should be 1
  if (reference_cell_type %in% target_cell_types) {
    neighbourhood_counts_df[neighbourhood_counts_df[[reference_cell_type]] != 0, paste(reference_cell_type, "_prop", sep = "")] <- 1
  }
  
  return(neighbourhood_counts_df)
}
