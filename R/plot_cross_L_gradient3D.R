plot_cross_L_gradient3D <- function(cross_L_gradient_df, reference_cell_type = NULL, target_cell_type = NULL) {
  
  plot_result <- reshape2::melt(cross_L_gradient_df, "radius", c("observed_cross_L", "expected_cross_L", "cross_L_ratio"))
  plot_result <- plot_result[plot_result$variable != "cross_L_ratio", ]
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross L-function gradient", x = "Radius", y = "Cross L-function value") +
    scale_colour_discrete(name = "", labels = c("Observed cross L", "Expected CSR cross L")) +
    theme_bw()
  
  if (!is.null(reference_cell_type) && !is.null(target_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", target_cell_type, sep = ""))
  }
  
  methods::show(fig)
  
  return(fig) 
}
