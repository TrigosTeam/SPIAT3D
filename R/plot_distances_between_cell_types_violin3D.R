## For scales parameter, use "free_x" or "free". "free_y" looks silly
plot_distances_between_cell_types_violin3D <- function(distances_df, scales = "free_x") {
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  pair <- distance <- NULL
  
  fig <- ggplot(distances_df, aes(x = pair, y = distance)) + 
    geom_violin() +
    facet_wrap(~pair, scales=scales, strip.position="bottom") +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
    labs(title="Cell distances", x = "Reference/Target pair", y = "Distance") +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
  
  message("Plots show mean ± sd")
  
  return(fig)
}
