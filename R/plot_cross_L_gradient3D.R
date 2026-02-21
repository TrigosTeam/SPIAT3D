#' @title Function to plot cross L gradient data.
#'
#' @description This function plots the cross L gradient data as a line graph, 
#'     showing cross L vs radius.
#' 
#' @param cross_L_gradient_df A data frame obtained from the output of the 
#'     calculate_cross_L_gradient3D function.
#'
#' @return A ggplot object containing the line graph.
#'
#' @examples
#' result <- calculate_cross_L_gradient3D(
#'     spe = SPIAT-3D::simulated_spe,
#'     reference_cell_type = "Tumour",
#'     target_cell_types = c("Tumour", "Immune"),
#'     radii = seq(20, 100, 10),
#'     feature_colname = "Cell.Type",
#'     plot_image = FALSE
#' )
#' 
#' fig <- plot_cross_L_gradient3D(
#'     cross_L_gradient_df = result
#' )
#' 
#' methods::show(fig)
#' 
#' @export

plot_cross_L_gradient3D <- function(cross_L_gradient_df) {
  
  # Get target cell types
  target_cell_types <- colnames(cross_L_gradient_df)[!colnames(cross_L_gradient_df) %in% c("reference", "expected", "radius")]
  
  # Re-format input data frame
  plot_result <- reshape2::melt(cross_L_gradient_df, "radius", c(target_cell_types, "expected"))
  
  # Plot
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross L-function gradient", x = "Radius", y = "Cross L-function value") +
    scale_colour_discrete(name = "") +
    theme_bw()
  
  return(fig) 
}
