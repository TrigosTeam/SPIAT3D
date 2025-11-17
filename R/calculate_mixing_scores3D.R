calculate_mixing_scores3D <- function(spe, 
                                      reference_cell_types, 
                                      target_cell_types, 
                                      radius, 
                                      feature_colname = "Cell.Type") {
  
  # Define result
  result <- data.frame()
  
  for (reference_cell_type in reference_cell_types) {
    
    for (target_cell_type in target_cell_types) {
      
      # No point getting mixing scores if comparing the same cell type
      if (reference_cell_type == target_cell_type) {
        next
      }
      
      # Get number of reference cells and target cells
      n_ref <- sum(spe[[feature_colname]] == reference_cell_type)
      n_tar <- sum(spe[[feature_colname]] == target_cell_type)
      
      
      # Can't get mixing scores if there are 0 or 1 reference cells
      if (n_ref == 0 || n_ref == 1) {
        result <-  rbind(result, 
                         c(reference_cell_type, 
                           target_cell_type, 
                           n_ref, 
                           n_tar, 
                           0, 
                           0, 
                           NA, 
                           NA))
      }
      
      
      ## Get cells in neighbourhood df
      cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spe,
                                                                      reference_cell_type,
                                                                      c(reference_cell_type, target_cell_type),
                                                                      radius,
                                                                      feature_colname,
                                                                      FALSE,
                                                                      FALSE)
      
      # Get number of ref-ref interactions
      # Halve it to avoid counting each ref-ref interaction twice
      n_ref_ref_interactions <- 0.5 * sum(cells_in_neighbourhood_df[[reference_cell_type]]) 
      
      # Get number of ref-tar interactions
      n_ref_tar_interactions <- sum(cells_in_neighbourhood_df[[target_cell_type]]) 
      
      
      # Can't get mixing scores if there are no target cells
      if (n_tar == 0) {
        
        result <-  rbind(result, 
                         c(reference_cell_type, 
                           target_cell_type, 
                           n_ref, 
                           0, 
                           0, 
                           n_ref_ref_interactions, 
                           NA, 
                           NA))
      }
      
      # Generic case: We have reference cells and target cells
      else {
        
        if (n_ref_ref_interactions != 0) {
          mixing_score <- n_ref_tar_interactions / n_ref_ref_interactions
          normalised_mixing_score <- 0.5 * mixing_score * n_ref / n_tar
        }
        else {
          mixing_score <- 0
          normalised_mixing_score <- 0
          methods::show(paste("There are no reference to reference interactions for", target_cell_type, "in the specified radius, cannot calculate mixing score"))
        }
        
        result <-  rbind(result, 
                         c(reference_cell_type, 
                           target_cell_type, 
                           n_ref, 
                           n_tar, 
                           n_ref_tar_interactions, 
                           n_ref_ref_interactions, 
                           mixing_score, 
                           normalised_mixing_score))
      }
    }
  }
  
  # Required column names of our output data frame
  colnames(result) <- c("ref_cell_type", 
                        "tar_cell_type", 
                        "n_ref_cells",
                        "n_tar_cells", 
                        "n_ref_tar_interactions",
                        "n_ref_ref_interactions", 
                        "mixing_score", 
                        "normalised_mixing_score")
  
  # Turn numeric data into numeric type
  result[ , 3:8] <- apply(result[ , 3:8], 2, as.numeric)
  
  return(result)
}
