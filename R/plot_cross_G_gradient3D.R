#' @title Function to plot cross G gradient data.
#'
#' @description This function plots the cross G gradient data as a line graph, 
#'     showing cross G vs radius.
#' 
#' @param cross_G_gradient_df A data frame obtained from the output of the 
#'     calculate_cross_G_gradient3D function.
#'     
#' @param reference_cell_type A string specifying the reference cell type. If 
#'     NULL or if `target_cell_type` is NULL, no reference cell type label will 
#'     be showed on the plot. Defaults to NULL.
#'     
#' @param target_cell_type A string specifying the target cell type. If NULL or
#'     if `reference_cell_type` is NULL, no target cell type label will be 
#'     showed on the plot. Defaults to NULL.
#'
#' @return A ggplot object containing the line graph.
#'
#' @examples
#' result <- calculate_cross_G_gradient3D(
#'     spe = SPIAT-3D::simulated_spe,
#'     reference_cell_type = "Tumour",
#'     target_cell_type = "Immune",
#'     radii = seq(20, 100, 10),
#'     feature_colname = "Cell.Type",
#'     plot_image = FALSE
#' )
#' 
#' fig <- plot_cross_G_gradient3D(
#'     cross_G_gradient_df = result,
#'     reference_cell_type = "Tumour",
#'     target_cell_type = "Immune"
#' )
#' 
#' methods::show(fig)
#' 
#' @export

plot_cross_G_gradient3D <- function(cross_G_gradient_df, 
                                    reference_cell_type = NULL, 
                                    target_cell_type = NULL) {
  
  # Re-format input data frame
  plot_result <- reshape2::melt(cross_G_gradient_df, "radius", c("observed_cross_G", "expected_cross_G"))
  
  # Plot
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Cross G-function gradient", x = "Radius", y = "Cross G-function value") +
    scale_colour_discrete(name = "", labels = c("Observed cross G", "Expected CSR cross G")) +
    theme_bw()
  
  # Add reference cell type and target cell type labels
  if (!is.null(reference_cell_type) && !is.null(target_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", target_cell_type, sep = ""))
  }
  
  return(fig) 
}
