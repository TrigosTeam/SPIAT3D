## For scales parameter, use "free_x" or "free". "free_y" looks silly
plot_cells_in_neighbourhood_violin3D <- function(cells_in_neighbourhood_df, reference_cell_type, scales = "free_x") {
  
  ## Target cell types will be all the columns except the first column
  target_cell_types <- colnames(cells_in_neighbourhood_df)[c(-1)]
  
  df <- reshape2::melt(cells_in_neighbourhood_df, measure.vars = target_cell_types)
  colnames(df) <- c("ref_cell_id", "tar_cell_type", "count")
  
  # setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  tar_cell_type <- count <- NULL
  
  fig <- ggplot(df, aes(x = tar_cell_type, y = count)) + 
    geom_violin() +
    facet_wrap(~tar_cell_type, scales=scales, strip.position="bottom") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(title=paste("Cells in neighbourhood of", reference_cell_type, "cells"), x = "Target cell type", y = "Number of cells") +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
  
  message("Plots show mean ± sd")
  
  return(fig)
}
