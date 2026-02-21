#' @title Function to plot distances between cell types as a set of violin 
#'     plots.
#'
#' @description This function plots the distances between cell types data as
#'     violin plots, showing the distances for each cell type pair.
#' 
#' @param distances_df A data frame obtained from the output of the 
#'     calculate_pairwise_distances_between_cell_types3D or
#'     calculate_minimum_distances_between_cell_types3D functions.
#'     
#' @param scales A string to control the facet wrap scales of the ggplot object.
#'     Either "fixed", "free_x", "free_y" or "free". Would not recommend 
#'     "free_y". Defaults to "free_x".
#'
#' @return A ggplot object containing the violin plots.
#'
#' @examples
#' minimum_distances <- calculate_minimum_distances_between_cell_types3D(
#'     spe = SPIAT-3D::simulated_spe,
#'     cell_types_of_interest = NULL,
#'     feature_colname = "Cell.Type",
#'     show_summary = TRUE,
#'     plot_image = FALSE
#' )
#' 
#' fig <- plot_distances_between_cell_types_violin3D(
#'     distances_df = minimum_distances,
#'     scales = "free_x"
#' )
#' 
#' methods::show(fig)
#' 
#' @export

plot_distances_between_cell_types_violin3D <- function(distances_df, 
                                                       scales = "free_x") {
  
  # Setting these variables to NULL as otherwise get "no visible binding for global variable" in R check
  pair <- distance <- NULL
  
  # Plot
  fig <- ggplot(distances_df, aes(x = pair, y = distance)) + 
    geom_violin() +
    facet_wrap(~pair, scales = scales, strip.position = "bottom") +
    theme_bw() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.title = element_text(hjust = 0.5)) +
    labs(title="Cell distances", x = "Reference/Target pair", y = "Distance") +
    stat_summary(fun.data = "mean_sdl", fun.args = list(mult= 1), colour = "red")
  
  message("Plots show mean ± sd")
  
  return(fig)
}
