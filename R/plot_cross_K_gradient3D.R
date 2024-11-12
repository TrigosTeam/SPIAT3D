plot_cross_K_gradient3D <- function(cross_K_gradient_df, reference_cell_type = NULL, target_cell_type = NULL) {
  
  plot_result <- reshape2::melt(cross_K_gradient_df, "radius", c("observed_cross_K", "expected_cross_K", "cross_K_ratio"))
  plot_result <- plot_result[plot_result$variable != "cross_K_ratio", ]
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross K-function gradient", x = "Radius", y = "Cross K-function value") +
    scale_colour_discrete(name = "", labels = c("Observed cross K", "Expected CSR cross K")) +
    theme_bw()
  
  if (!is.null(reference_cell_type) && !is.null(target_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", target_cell_type, sep = ""))
  }
  
  methods::show(fig)
 
  return(fig) 
}
