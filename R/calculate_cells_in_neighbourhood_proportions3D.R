calculate_cells_in_neighbourhood_proportions3D <- function(spe, 
                                                           reference_cell_type, 
                                                           target_cell_types, 
                                                           radius, 
                                                           feature_colname = "Cell.Type") {
  
  ## Get cells in neighbourhood df
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spe,
                                                                  reference_cell_type,
                                                                  c(reference_cell_type, target_cell_types),
                                                                  radius,
                                                                  feature_colname,
                                                                  FALSE,
                                                                  FALSE)
  
  if (is.null(cells_in_neighbourhood_df)) return(NULL)
  
  cells_in_neighbourhood_df[ , paste(target_cell_types, "_prop", sep = "")] <- 
    cells_in_neighbourhood_df[ , target_cell_types] / (cells_in_neighbourhood_df[ , target_cell_types] + cells_in_neighbourhood_df[ , reference_cell_type])
  
  # If reference cell type is in target cell types, proportion should be 1
  if (reference_cell_type %in% target_cell_types) {
    cells_in_neighbourhood_df[cells_in_neighbourhood_df[[reference_cell_type]] != 0, paste(reference_cell_type, "_prop", sep = "")] <- 1
  }
  
  return(cells_in_neighbourhood_df)
}
