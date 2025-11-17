plot_cells_in_neighbourhood_gradient3D <- function(cells_in_neighbourhood_gradient_df, 
                                                   reference_cell_type = NULL) {
  
  plot_result <- reshape2::melt(cells_in_neighbourhood_gradient_df, "radius")
  
  fig <- ggplot(plot_result, aes(radius, value, color = variable)) + 
    geom_line() + 
    labs(title = "Average cells in neighbourhood gradient", x = "Radius", y = "Average cells in neighbourhood") + 
    scale_color_discrete(name = "Cell type") +
    theme_bw()
  
  if (!is.null(reference_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", paste(colnames(cells_in_neighbourhood_gradient_df)[seq(ncol(cells_in_neighbourhood_gradient_df) - 1)], collapse = ", "), sep = ""))
  }
  
  return(fig)
}
