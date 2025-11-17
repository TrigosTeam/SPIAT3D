plot_mixing_scores_gradient3D <- function(mixing_scores_gradient_df, 
                                          metric = "MS") {
  
  if (!metric %in% c("MS", "NMS")) {
    stop("'metric' should be 'MS' or 'NMS', for mixing score and normalised mixing score respectively.")
  }
  
  if (metric == "NMS") {
    plot_result <- mixing_scores_gradient_df
    plot_result$expected_normalised_mixing_score <- 1
    plot_result <- reshape2::melt(plot_result, "radius", c("normalised_mixing_score", "expected_normalised_mixing_score"))
    
    fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
      geom_line() +
      labs(title = "Normalised mixing score (NMS) gradient", 
           subtitle = paste("Reference: ", mixing_scores_gradient_df$ref_cell_type[1], ", Target: ", mixing_scores_gradient_df$tar_cell_type[1], sep = ""), 
           x = "Radius", y = "NMS") +
      scale_colour_discrete(name = "", labels = c("Observed NMS", "Expected CSR NMS")) +
      theme_bw() 
  }
  else if (metric == "MS") {
    plot_result <- mixing_scores_gradient_df
    n_tar_cells <- plot_result$n_tar_cells[1]
    n_ref_cells <- plot_result$n_ref_cells[1]
    plot_result$expected_mixing_score <- n_tar_cells * n_ref_cells / ((n_ref_cells - 1) * n_ref_cells / 2)
    plot_result <- reshape2::melt(plot_result, "radius", c("mixing_score", "expected_mixing_score"))
    
    fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
      geom_line() +
      labs(title = "Mixing score (MS) gradient", 
           subtitle = paste("Reference: ", mixing_scores_gradient_df$ref_cell_type[1], ", Target: ", mixing_scores_gradient_df$tar_cell_type[1], sep = ""), 
           x = "Radius", y = "MS") +
      scale_colour_discrete(name = "", labels = c("Observed MS", "Expected CSR MS  ")) +
      theme_bw()  
  }
  return(fig)
}
