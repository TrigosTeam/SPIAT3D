plot_cross_G_gradient3D <- function(cross_G_gradient_df, 
                                    reference_cell_type = NULL, 
                                    target_cell_type = NULL) {
  
  plot_result <- reshape2::melt(cross_G_gradient_df, "radius", c("observed_cross_G", "expected_cross_G"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross G-function gradient", x = "Radius", y = "Cross G-function value") +
    scale_colour_discrete(name = "", labels = c("Observed cross G", "Expected CSR cross G")) +
    theme_bw()
  
  if (!is.null(reference_cell_type) && !is.null(target_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", target_cell_type, sep = ""))
  }
  
  return(fig) 
}
