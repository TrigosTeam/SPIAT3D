#' @title Function to plot neighbourhood counts data as a set of violin plots.
#'
#' @description This function plots the neighbourhood counts data as violin 
#'     plots, showing the neighbourhood counts for each reference and target
#'     cell type pair.
#' 
#' @param neighbourhood_counts_df A data frame obtained from the output of the 
#'     calculate_neighbourhood_counts3D function.
#'     
#' @param reference_cell_type A string specifying the reference cell type.
#'     
#' @param scales A string to control the facet wrap scales of the ggplot object.
#'     Either "fixed", "free_x", "free_y" or "free". Would not recommend 
#'     "free_y". Defaults to "free_x".
#'
#' @return A ggplot object containing the violin plots.
#'
#' @examples
#' result <- calculate_neighbourhood_counts3D(
#'     spe = SPIAT-3D::simulated_spe,
#'     reference_cell_type = "Tumour",
#'     target_cell_types = c("Tumour", "Immune"),
#'     radius = 30,
#'     feature_colname = "Cell.Type",
#'     show_summary = TRUE,
#'     plot_image = TRUE
#' )
#' 
#' fig <- plot_neighbourhood_counts_violin3D(
#'     neighbourhood_counts_df = result,
#'     reference_cell_type = "Tumour",
#'     scales = "free_x"
#' )
#' 
#' @export

plot_neighbourhood_counts_violin3D <- function(neighbourhood_counts_df, 
                                               reference_cell_type, 
                                               scales = "free_x") {
  
  # Target cell types will be all the columns except the first column
  target_cell_types <- colnames(neighbourhood_counts_df)[c(-1)]
  
  # Re-format input data frame
  df <- reshape2::melt(neighbourhood_counts_df, measure.vars = target_cell_types)
  colnames(df) <- c("ref_cell_id", "tar_cell_type", "count")
  
  # Setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  tar_cell_type <- count <- NULL
  
  # Plot
  fig <- ggplot(df, aes(x = tar_cell_type, y = count)) + 
    geom_violin() +
    facet_wrap(~tar_cell_type, scales=scales, strip.position="bottom") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(title = paste("Neighbourhood counts of", reference_cell_type, "cells"), x = "Target cell type", y = "Neighbourhood counts") +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
  
  message("Plots show mean ± sd")
  
  return(fig)
}
