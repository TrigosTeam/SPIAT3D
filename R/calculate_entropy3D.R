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

  ## Get entropy for each row
  cells_in_neighbourhood_proportion_df$entropy <- apply(cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_prop", sep = "")],
                                                        1,
                                                        function(x) -1 * sum(x * log(x, length(target_cell_types))))
  cells_in_neighbourhood_proportion_df$entropy <- ifelse(cells_in_neighbourhood_proportion_df$total > 0 & is.nan(cells_in_neighbourhood_proportion_df$entropy), 
                                                         0,
                                                         cells_in_neighbourhood_proportion_df$entropy)
  
  return(cells_in_neighbourhood_proportion_df)
}
