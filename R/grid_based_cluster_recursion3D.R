grid_based_cluster_recursion3D <- function(df,  # Using a df is much faster than using a spe
                                           cell_types_of_interest,
                                           threshold_cell_proportion,
                                           x, y, z, l, w, h,
                                           feature_colname,
                                           answer) {
  
  # Look at cells only in the current grid prism
  df <- df[df$Cell.X.Position >= x &
             df$Cell.X.Position < (x + l) &
             df$Cell.Y.Position >= y &
             df$Cell.Y.Position < (y + w) &
             df$Cell.Z.Position >= z &
             df$Cell.Z.Position < (z + h), ]
  
  # Get cell types from spe grid prism
  cell_types <- df[[feature_colname]]
  
  # Number of cells in prism is getting too small
  if (length(cell_types) <= 2) return(data.frame())
  
  # Get total cell proportion for chosen cell_types_of_interest
  cell_proportion <- mean(cell_types %in% cell_types_of_interest)
  
  # Keep grid prism if cell proportion is above the threshold cell proportion
  if (cell_proportion >= threshold_cell_proportion) {
    return(data.frame(x, y, z, l, w, h))
  }
  
  # some cell_types_of_interest still in the grid prism, check sub-grid prisms (8 to check)
  else if (cell_proportion > 0) {
    # (0, 0, 0)
    answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x, y, z, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    
    # (0.5, 0, 0)
    answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x + l/2, y, z, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    
    # (0, 0.5, 0)
    answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x, y + w/2, z, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    # (0.5, 0.5, 0)
    answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x + l/2, y + w/2, z, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    
    # (0, 0, 0.5)
    answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x, y, z + h/2, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    
    # (0.5, 0, 0.5)
    answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x + l/2, y, z + h/2, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    
    # (0, 0.5, 0.5)
    answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x, y + w/2, z + h/2, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    # (0.5, 0.5, 0.5)
    answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                           cell_types_of_interest,
                                                           threshold_cell_proportion,
                                                           x + l/2, y + w/2, z + h/2, l/2, w/2, h/2,
                                                           feature_colname,
                                                           data.frame()))
    
    return(answer)
  }
  
  # cell proportion is zero
  else {
    return(data.frame())
  }
}
