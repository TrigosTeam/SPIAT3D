calculate_entropy3D <- function(spe,
                                reference_cell_type,
                                target_cell_types,
                                radius,
                                feature_colname = "Cell.Type") {
  
  # Check target_cell_types
  if (!(is.character(target_cell_types) && length(target_cell_types) >= 2)) {
    stop("`target_cell_types` is not a character vector with at least 2 cell types.")
  }
  
  ## Users should ensure include the reference_cell_type as one of the target_cell_types
  cells_in_neighbourhood_proportion_df <- calculate_cells_in_neighbourhood_proportions3D(spe,
                                                                                         reference_cell_type,
                                                                                         target_cell_types,
                                                                                         radius,
                                                                                         feature_colname)
  
  if (is.null(cells_in_neighbourhood_proportion_df)) return(NULL)
  
  ## Get entropy for target_cell_type
  cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_entropy", sep = "")] <- 
    -1 * 
    (cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_prop", sep = "")] * log(cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_prop", sep = "")], 2) +
       (1 - cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_prop", sep = "")]) * log(1 - (cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_prop", sep = "")]), 2))
  
  return(cells_in_neighbourhood_proportion_df)
}
