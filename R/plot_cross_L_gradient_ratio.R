plot_cross_L_gradient_ratio3D <- function(cross_L_gradient_df, reference_cell_type = NULL, target_cell_type = NULL) {
  
  plot_result <- data.frame(radius = cross_L_gradient_df$radius,
                            observed_cross_L_gradient_ratio = cross_L_gradient_df$cross_L_ratio,
                            expected_cross_L_gradient_ratio = 1)
  
  plot_result <- reshape2::melt(plot_result, "radius", c("observed_cross_L_gradient_ratio", "expected_cross_L_gradient_ratio"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross L-function ratio gradient", x = "Radius", y = "Cross L-function ratio") +
    scale_colour_discrete(name = "", labels = c("Observed cross L ratio", "Expected CSR cross L ratio")) +
    theme_bw()
  
  if (!is.null(reference_cell_type) && !is.null(target_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", target_cell_type, sep = ""))
  }
  
  methods::show(fig)
  
  return(fig) 
}
