#' @title Plot 3D grid metrics (continuous method).
#'
#' @description This functions plots the 3D grid metrics, showing the
#'     proportion/entropy values of each rectangular prism using a continuous 
#'     color scale for the proportion/entropy values. 
#' 
#' @param grid_metrics A data frame containing the proportion/entropy and 
#'     spatial information for each rectangular prism. Obtained from the output
#'     of the calculate_cell_proportion_grid_metrics3D and 
#'     calculate_entropy_grid_metrics3D functions.
#' @param metric_colname A string specifying the name of the column in 
#'     `grid_metrics` containing the proportion/entropy information. Should be
#'     'proportion' or 'entropy'.
#'
#' @return A Plotly object containing the 3D grid metrics plot.
#'
#' @examples
#' cell_prop_grid_metrics <- calculate_cell_proportion_grid_metrics3D(
#'     spe = SPIAT-3D::simulated_spe,
#'     n_splits = 10,
#'     reference_cell_types = c("Tumour"),
#'     target_cell_types = c("Immune"),
#'     feature_colname = "Cell.Type",
#'     plot_image = FALSE
#' )
#' 
#' fig <- plot_grid_metrics_continuous3D(
#'     grid_metrics = cell_prop_grid_metrics,
#'     metric_colname = "proportion"
#' )
#' 
#' methods::show(fig)
#' 
#' @export

plot_grid_metrics_continuous3D <- function(grid_metrics, 
                                           metric_colname) {
  
  ## Check input parameters
  if (!(is.character(metric_colname) && metric_colname %in% c("proportion", "entropy"))) {
    stop("`metric_colname` is not 'proportion' or 'entropy'.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  
  ## Color of each dot is related to its proportion/entropy
  pal <- colorRampPalette(hcl.colors(n = 5, palette = "Red-Blue", rev = TRUE))
  
  ## Add size column and for NA proportion/entropy values, make the size small
  grid_metrics$size <- ifelse(is.na(grid_metrics[[metric_colname]]), 3, 10)
  
  # Plot
  fig <- plot_ly(grid_metrics,
                 type = "scatter3d",
                 mode = 'markers',
                 x = ~x_coord,
                 y = ~y_coord,
                 z = ~z_coord,
                 color = as.formula(paste0('~', metric_colname)),
                 colors = pal(nrow(grid_metrics)),
                 marker = list(size = ~size),
                 symbol = 1,
                 symbols = "square")
  
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'x'),
                                     yaxis = list(title = 'y'),
                                     zaxis = list(title = 'z')))
  
  return(fig)
}
