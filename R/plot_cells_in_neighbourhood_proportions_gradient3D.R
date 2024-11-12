plot_cells_in_neighbourhood_proportions_gradient3D <- function(cells_in_neighbourhood_proportions_gradient_df, reference_cell_type = NULL) {
  
  plot_result <- reshape2::melt(cells_in_neighbourhood_proportions_gradient_df, id.vars = c("radius"))
  fig <- ggplot(plot_result, aes(radius, value, color = variable)) +
    geom_point() +
    geom_line() +
    labs(title = "Average cells in neighbourhood proportions gradient", x = "Radius", y = "Cell proportion", color = "Cell type") +
    theme_bw() +
    ylim(0, 1)
  
  if (!is.null(reference_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", paste(colnames(cells_in_neighbourhood_proportions_gradient_df)[seq(ncol(cells_in_neighbourhood_proportions_gradient_df) - 1)], collapse = ", "), sep = ""))
  }
  
  methods::show(fig)
  
  return(fig)
}
