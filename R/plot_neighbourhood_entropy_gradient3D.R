#' @title Function to plot neighbourhood entropy gradient data.
#'
#' @description This function plots the neighbourhood entropy gradient data as a
#'     line graph, showing average neighbourhood entropy vs radius.
#' 
#' @param neighbourhood_entropy_gradient_df A data frame obtained from the
#'     output of calculate_neighbourhood_entropy_gradient3D function.
#'     
#' @param reference_cell_type A string specifying the reference cell type. If 
#'     NULL, no reference cell type label will be showed on the plot. Defaults 
#'     to NULL.
#'
#' @return A ggplot object containing the line graph.
#'
#' @examples
#' result <- calculate_neighbourhood_entropy_gradient3D(
#'     spe = SPIAT-3D::simulated_spe,
#'     reference_cell_type = "Tumour",
#'     target_cell_types = c("Tumour", "Immune"),
#'     radii = seq(20, 100, 10),
#'     feature_colname = "Cell.Type",
#'     plot_image = FALSE
#' )
#' 
#' fig <- plot_neighbourhood_entropy_gradient3D(
#'     neighbourhood_entropy_gradient_df = result,
#'     reference_cell_type = "Tumour"
#' )
#' 
#' methods::show(fig)
#' 
#' @export

plot_neighbourhood_entropy_gradient3D <- function(neighbourhood_entropy_gradient_df, 
                                                  reference_cell_type = NULL) {
  
  # Re-format input data frame
  plot_result <- reshape2::melt(neighbourhood_entropy_gradient_df, id.vars = c("radius"))
  
  # Plot
  fig <- ggplot(plot_result, aes(radius, value, color = variable)) +
    geom_point() +
    geom_line() +
    labs(title = "Neighbourhood entropy gradient", x = "Radius", y = "Neighbourhood entropy", color = "Cell type") +
    theme_bw() +
    ylim(0, 1)
  
  # Add reference cell type label
  if (!is.null(reference_cell_type)) {
    fig <- fig + labs(subtitle = paste("Reference: ", reference_cell_type, ", Target: ", paste(colnames(neighbourhood_entropy_gradient_df)[seq(ncol(neighbourhood_entropy_gradient_df) - 1)], collapse = ", "), sep = ""))
  }
  
  return(fig)
}
