plot_cross_K_gradient_ratio3D <- function(cross_K_gradient_df) {
  
  target_cell_types <- colnames(cross_K_gradient_df)[!colnames(cross_K_gradient_df) %in% c("reference", "expected", "radius")]
  
  # Normalize columns by 'expected'
  for (target_cell_type in target_cell_types) {
    cross_K_gradient_df[[target_cell_type]] <- cross_K_gradient_df[[target_cell_type]] / cross_K_gradient_df$expected
  }
  cross_K_gradient_df$expected <- 1
  
  plot_result <- reshape2::melt(cross_K_gradient_df, "radius", c(target_cell_types, "expected"))
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross K-function gradient ratio", x = "Radius", y = "Cross K-function ratio") +
    scale_colour_discrete(name = "") +
    theme_bw()
  
  return(fig) 
}
