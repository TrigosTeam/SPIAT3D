calculate_cross_L3D <- function(spe, 
                               reference_cell_type, 
                               target_cell_type, 
                               radius, 
                               feature_colname = "Cell.Type") {
  
  result <- calculate_cross_K3D(spe = spe,
                                reference_cell_type = reference_cell_type,
                                target_cell_type = target_cell_type,
                                radius = radius,
                                feature_colname = feature_colname)
  
  colnames(result) <- c("observed_cross_L", "expected_cross_L", "cross_L_ratio")
  result <- (result / (4 * pi / 3)) ^ (1/3)
  
  return(result)
}
