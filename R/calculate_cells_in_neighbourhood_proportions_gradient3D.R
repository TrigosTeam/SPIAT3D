calculate_cells_in_neighbourhood_proportions_gradient3D <- function(spe, 
                                                                    reference_cell_type, 
                                                                    target_cell_types, 
                                                                    radii, 
                                                                    feature_colname = "Cell.Type",
                                                                    plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (i in seq(length(radii))) {
    cell_proportions_neighbourhood_proportions_df <- calculate_cells_in_neighbourhood_proportions3D(spe,
                                                                                                    reference_cell_type,
                                                                                                    target_cell_types,
                                                                                                    radii[i],
                                                                                                    feature_colname)
    
    if (is.null(cell_proportions_neighbourhood_proportions_df)) return(NULL)
    
    result[i, ] <- apply(cell_proportions_neighbourhood_proportions_df[ , paste(target_cell_types, "_prop", sep = "")], 2, mean, na.rm = T)
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  # Plot
  if (plot_image) plot_cells_in_neighbourhood_proportions_gradient3D(result, reference_cell_type)
  
  return(result)
}
