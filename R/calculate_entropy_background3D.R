calculate_entropy_background3D <- function(spe,
                                           cell_types_of_interest, 
                                           feature_colname = "Cell.Type") {
  
  # NULL case: entropy is undefined
  if (is.null(cell_types_of_interest) == 0) return(NA)
  
  # One cell type case: entropy is 0
  if (is.character(cell_types_of_interest) && length(cell_types_of_interest) == 1) return(0)
  
  cell_proportions_data <- calculate_cell_proportions3D(spe, cell_types_of_interest, feature_colname, FALSE)
  
  # Calculate entropy of the entire image
  entropy <- -1 * sum(cell_proportions_data$proportion * log(cell_proportions_data$proportion, length(cell_proportions_data$proportion)))
  
  return(entropy) 
}
