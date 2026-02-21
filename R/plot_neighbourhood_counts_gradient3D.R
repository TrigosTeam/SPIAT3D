#' @title Function to plot neighbourhood counts gradient data.
#'
#' @description This function plots the neighbourhood counts gradient data as a 
#'     line graph, showing average neighbourhood counts vs radius.
#' 
#' @param neighbourhood_counts_gradient_df A data frame obtained from the output 
#'     of calculate_neighbourhood_counts3D function.
#'     
#' @param reference_cell_type A string specifying the reference cell type. If 
#'     NULL, no reference cell type label will be showed on the plot. Defaults 
#'     to NULL.
#'
#' @return A ggplot object containing the line graph.
#'
#' @examples
#' result <- calculate_neighbourhood_counts_gradient3D(
#'     spe = SPIAT-3D::simulated_spe,
#'     reference_cell_type = "Tumour",
#'     target_cell_types = c("Tumour", "Immune"),
#'     radii = seq(20, 100, 10),
#'     feature_colname = "Cell.Type",
#'     plot_image = FALSE
#' )
#' 
#' fig <- plot_neighbourhood_counts_gradient3D(
#'     neighbourhood_counts_gradient_df = result,
#'     reference_cell_type = "Tumour"
#' )
#' 
#' methods::show(fig)
#' 
#' @export

plot_neighbourhood_counts_gradient3D <- function(neighbourhood_counts_gradient_df, 
                                                 reference_cell_type = NULL) {
  
  # Re-format input data frame
  plot_result <- reshape2::melt(neighbourhood_counts_gradient_df, "radius")
  
  # Plot
  fig <- ggplot(plot_result, aes(radius, value, color = variable)) + 
    geom_line() + 
    labs(title = "Average neighbourhood counts gradient", x = "Radius", y = "Average neighbourhood counts") + 
    scale_color_discrete(name = "Cell type") +
    theme_bw()
  
  # Add reference cell type label
  if (!is.null(reference_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, sep = ""))
  }
  
  return(fig)
}
