plot_entropy_gradient3D <- function(entropy_gradient_df, expected_entropy = NULL, reference_cell_type = NULL, target_cell_types = NULL) {
  
  plot_result <- entropy_gradient_df

  if (!is.null(expected_entropy)) {
    if (!is.numeric(expected_entropy) || length(expected_entropy) != 1) stop("Please enter a single number for expected_entropy")
    plot_result$expected_entropy <- expected_entropy
    plot_result <- reshape2::melt(plot_result, "radius", c("entropy", "expected_entropy"))
    labels <- c("Observed entropy", "Expected CSR entropy")
  }
  else {
    plot_result <- reshape2::melt(plot_result, "radius", c("entropy"))
    labels <- c("Observed entropy")
  }
    
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Average entropy gradient", x = "Radius", y = "Entropy") +
    scale_colour_discrete(name = "", labels = labels) +
    theme_bw()
    
  if (!is.null(reference_cell_type) && !is.null(target_cell_types)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", target_cell_types, sep = ""))
  }
  
  methods::show(fig)
  
  return(fig)
}
