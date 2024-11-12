calculate_cells_in_neighbourhood_proportions3D <- function(spe, 
                                                           reference_cell_type, 
                                                           target_cell_types, 
                                                           radius, 
                                                           feature_colname = "Cell.Type") {
  
  ## Get cells in neighbourhood df
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spe,
                                                                  reference_cell_type,
                                                                  target_cell_types,
                                                                  radius,
                                                                  feature_colname,
                                                                  FALSE,
                                                                  FALSE)
  
  if (is.null(cells_in_neighbourhood_df)) return(NULL)
  
  ## Get total number of target cells for each row (first column is the reference cell id column, so we exclude it)
  cells_in_neighbourhood_df$total <- apply(cells_in_neighbourhood_df[ , c(-1)], 1, sum)
  
  cells_in_neighbourhood_df[ , paste(target_cell_types, "_prop", sep = "")] <- cells_in_neighbourhood_df[ , target_cell_types] / cells_in_neighbourhood_df$total
  
  return(cells_in_neighbourhood_df)
}
