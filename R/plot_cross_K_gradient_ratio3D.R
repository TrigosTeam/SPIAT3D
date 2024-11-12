plot_cross_K_gradient_ratio3D <- function(cross_K_gradient_df, reference_cell_type = NULL, target_cell_type = NULL) {
  
  plot_result <- data.frame(radius = cross_K_gradient_df$radius,
                            observed_cross_K_gradient_ratio = cross_K_gradient_df$cross_K_ratio,
                            expected_cross_K_gradient_ratio = 1)
  
  plot_result <- reshape2::melt(plot_result, "radius", c("observed_cross_K_gradient_ratio", "expected_cross_K_gradient_ratio"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross K-function ratio gradient", x = "Radius", y = "Cross K-function ratio") +
    scale_colour_discrete(name = "", labels = c("Observed cross K ratio", "Expected CSR cross K ratio")) +
    theme_bw()
  
  if (!is.null(reference_cell_type) && !is.null(target_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", target_cell_type, sep = ""))
  }
  
  methods::show(fig)
 
  return(fig) 
}
