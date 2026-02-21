#' @title Function to plot co-occurrence gradient data.
#'
#' @description This function plots the co-occurrence gradient data as a line 
#'     graph, showing co-occurrence vs radius.
#' 
#' @param co_occurrence_gradient_df A data frame obtained from the output of the
#'     calculate_co_occurrence_gradient3D function.
#'
#' @return A ggplot object containing the line graph.
#'
#' @examples
#' result <- calculate_co_occurrence_gradient3D(
#'     spe = SPIAT-3D::simulated_spe,
#'     reference_cell_type = "Tumour",
#'     target_cell_types = c("Tumour", "Immune"),
#'     radii = seq(20, 100, 10),
#'     feature_colname = "Cell.Type",
#'     plot_image = FALSE
#' )
#' 
#' fig <- plot_co_occurrence_gradient3D(
#'     co_occurrence_gradient_df = result
#' )
#' 
#' methods::show(fig)
#' 
#' @export

plot_co_occurrence_gradient3D <- function(co_occurrence_gradient_df) {
  
  # Get target cell types
  target_cell_types <- colnames(co_occurrence_gradient_df)
  target_cell_types <- target_cell_types[!target_cell_types %in% c("reference", "radius")]
  
  # Add expected co-occurrence column (expected value is 1)
  co_occurrence_gradient_df$expected <- 1
  
  # Re-format input data frame
  plot_result <- reshape2::melt(co_occurrence_gradient_df, "radius", c(target_cell_types, "expected"))
  
  # Plot
  fig <- ggplot(plot_result, aes(x = radius, y = value, color = variable)) +
    geom_line() +
    labs(title = "Co-occurrence gradient", x = "Radius", y = "Co-occurrence value") +
    scale_colour_discrete(name = "") +
    theme_bw()
  
  return(fig) 
}
