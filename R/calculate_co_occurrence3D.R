calculate_co_occurrence3D <- function(spe, 
                                      reference_cell_type, 
                                      target_cell_types, 
                                      radius, 
                                      feature_colname = "Cell.Type") {
  
  # Get all cell types in spe
  all_cell_types <- unique(spe[[feature_colname]])
  
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spe,
                                                                  reference_cell_type,
                                                                  all_cell_types,
                                                                  radius,
                                                                  feature_colname,
                                                                  F,
                                                                  F)
  
  cells_in_neighbourhood_df$total <- rowSums(cells_in_neighbourhood_df[, -1], na.rm = TRUE)
  
  result <- data.frame(reference = reference_cell_type)
  
  # Get total number of cells in spe
  n_cells_in_spe <- length(spe[[feature_colname]])
  
  # Get total number of cells in radius around reference cell type
  n_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_df$total)
  
  for (target_cell_type in target_cell_types) {
    
    # Get total number of target cells in radius around reference cell type
    n_target_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_df[[target_cell_type]])
    
    # Get proportion of target cells in radius around reference cell type
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    
    # Get proportion of target cell type in spe
    n_target_cells_in_spe <- sum(spe[[feature_colname]] == target_cell_type)
    target_cell_type_proportion_in_spe <- n_target_cells_in_spe / n_cells_in_spe
    
    # Get co-occurence value for taget cell type
    target_cell_type_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_spe
    
    # Add to result data frame
    result[[target_cell_type]] <- target_cell_type_co_occurrence
  }
  
  return(result)
}
