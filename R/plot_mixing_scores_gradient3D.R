plot_mixing_scores_gradient3D <- function(mixing_scores_gradient_df) {
  
  plot_result1 <- mixing_scores_gradient_df
  plot_result1$expected_normalised_mixing_score <- 1
  plot_result1 <- reshape2::melt(plot_result1, "radius", c("normalised_mixing_score", "expected_normalised_mixing_score"))
  
  fig1 <- ggplot(plot_result1, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Normalised mixing score (NMS) gradient", 
         subtitle = paste("Reference: ", mixing_scores_gradient_df$ref_cell_type[1], ", Target: ", mixing_scores_gradient_df$tar_cell_type[1], sep = ""), 
         x = "Radius", y = "NMS") +
    scale_colour_discrete(name = "", labels = c("Observed NMS", "Expected CSR NMS")) +
    theme_bw()
  
  
  plot_result2 <- mixing_scores_gradient_df
  n_tar_cells <- plot_result2$n_tar_cells[1]
  n_ref_cells <- plot_result2$n_ref_cells[1]
  plot_result2$expected_mixing_score <- n_tar_cells * n_ref_cells / ((n_ref_cells - 1) * n_ref_cells / 2)
  plot_result2 <- reshape2::melt(plot_result2, "radius", c("mixing_score", "expected_mixing_score"))
  
  fig2 <- ggplot(plot_result2, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Mixing score (MS) gradient", 
         subtitle = paste("Reference: ", mixing_scores_gradient_df$ref_cell_type[1], ", Target: ", mixing_scores_gradient_df$tar_cell_type[1], sep = ""), 
         x = "Radius", y = "MS") +
    scale_colour_discrete(name = "", labels = c("Observed MS", "Expected CSR MS  ")) +
    theme_bw()
  
  combined_fig <- plot_grid(fig1, fig2, nrow = 2)
  
  methods::show(combined_fig)
  
  return(combined_fig)
}
