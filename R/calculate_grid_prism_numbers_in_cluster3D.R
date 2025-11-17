### Start from the grid_prism with the maximum cell proportion.
## Look left, right, forward, back, up and down and see if that grid_prism has at least threshold cell proportion value
## If it does, add it to the answer
## Keep doing this until adjacent grid prisms don't have above threshold, or if you hit a boundary, or it has already been removed
## Return a vector containing all the grid prism numbers which COULD be part of the cluster
calculate_grid_prism_numbers_in_cluster3D <- function(curr_grid_prism_number, 
                                                      grid_prism_cell_proportions, 
                                                      threshold_cell_proportion,
                                                      n_splits,
                                                      answer) {
  
  ## If answer already has curr_grid_prism_number, go back
  if (as.character(curr_grid_prism_number) %in% answer) return(answer)
  
  grid_prism_numbers <- names(grid_prism_cell_proportions)
  
  ## If curr_grid_prism_number has already been removed from grid_prism_numbers, go back
  if (!(as.character(curr_grid_prism_number) %in% grid_prism_numbers)) return(answer)
  
  
  if (grid_prism_cell_proportions[as.character(curr_grid_prism_number)] > threshold_cell_proportion) {
    
    answer <- c(answer, as.character(curr_grid_prism_number))
    
    ### CHECK RIGHT, LEFT, FORWARD, BACKWARD, UP, DOWN
    ## Need to check if going right, left, forward, backward, up or down is possible
    
    # Right
    if (curr_grid_prism_number%%n_splits != 0) {
      answer <- calculate_grid_prism_numbers_in_cluster3D(curr_grid_prism_number + 1,
                                                          grid_prism_cell_proportions,
                                                          threshold_cell_proportion,
                                                          n_splits,
                                                          answer)
    }
    
    # Left
    if (curr_grid_prism_number%%n_splits != 1) {
      answer <- calculate_grid_prism_numbers_in_cluster3D(curr_grid_prism_number - 1,
                                                          grid_prism_cell_proportions,
                                                          threshold_cell_proportion,
                                                          n_splits,
                                                          answer)
    }
    
    # Forward
    if ((curr_grid_prism_number - 1)%%(n_splits^2) < n_splits^2 - n_splits) {
      answer <- calculate_grid_prism_numbers_in_cluster3D(curr_grid_prism_number + n_splits,
                                                          grid_prism_cell_proportions,
                                                          threshold_cell_proportion,
                                                          n_splits,
                                                          answer)
    }
    
    # Backward
    if (curr_grid_prism_number%%(n_splits^2) > n_splits) {
      answer <- calculate_grid_prism_numbers_in_cluster3D(curr_grid_prism_number - n_splits,
                                                          grid_prism_cell_proportions,
                                                          threshold_cell_proportion,
                                                          n_splits,
                                                          answer)
    }
    
    # Up
    if (curr_grid_prism_number <= n_splits^3 - n_splits^2) {
      answer <- calculate_grid_prism_numbers_in_cluster3D(curr_grid_prism_number + n_splits^2,
                                                          grid_prism_cell_proportions,
                                                          threshold_cell_proportion,
                                                          n_splits,
                                                          answer)
    }
    
    # Down
    if (curr_grid_prism_number > n_splits^2) {
      answer <- calculate_grid_prism_numbers_in_cluster3D(curr_grid_prism_number - n_splits^2,
                                                          grid_prism_cell_proportions,
                                                          threshold_cell_proportion,
                                                          n_splits,
                                                          answer)
    }
  }
  
  return(answer)
}
