calculate_cross_K_gradient3D <- function(spe, 
                                         reference_cell_type, 
                                         target_cell_type, 
                                         radii, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 3))
  colnames(result) <- c("observed_cross_K", 
                        "expected_cross_K",
                        "cross_K_ratio")
  
  for (i in seq(length(radii))) {
    cross_K_df <- calculate_cross_K3D(spe,
                                      reference_cell_type,
                                      target_cell_type,
                                      radii[i],
                                      feature_colname)
    
    result[i, ] <- cross_K_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig1 <- plot_cross_K_gradient3D(result, reference_cell_type, target_cell_type)
    fig2 <- plot_cross_K_gradient_ratio3D(result, reference_cell_type, target_cell_type)
    
    combined_fig <- plot_grid(fig1, fig2, nrow = 2)
    methods::show(combined_fig)
  }
  
  return(result)
}