plot_co_occurrence_gradient3D <- function(co_occurrence_gradient_df) {
  
  target_cell_types <- colnames(co_occurrence_gradient_df)
  target_cell_types <- target_cell_types[!target_cell_types %in% c("reference", "radius")]
  
  plot_result <- reshape2::melt(co_occurrence_gradient_df, "radius", target_cell_types)
  
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Co-occurrence gradient", x = "Radius", y = "Co-occurrence value") +
    scale_colour_discrete(name = "", labels = target_cell_types) +
    theme_bw()
  
  methods::show(fig)
  
  return(fig) 
}
