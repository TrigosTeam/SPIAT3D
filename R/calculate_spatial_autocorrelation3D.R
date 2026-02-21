#' @title Calculate spatial autocorrelation of 3D grid metrics.
#'
#' @description This functions calculates the spatial autocorrelation of 3D
#'     grid metrics, which quantifies the interaction of rectangular prisms with
#'     similar proportion/entropy values in the 3D grid metrics.
#' 
#' @param grid_metrics A data frame containing the proportion/entropy and 
#'     spatial information for each rectangular prism. Obtained from the output
#'     of the calculate_cell_proportion_grid_metrics3D and 
#'     calculate_entropy_grid_metrics3D functions.
#' @param metric_colname A string specifying the name of the column in 
#'     `grid_metrics` containing the proportion/entropy information. Should be
#'     'proportion' or 'entropy'.
#' @param weight_method A string specifying the method to use to determining 
#'     whether rectangular prisms are interacting or not, via a weight. A weight
#'     of 1 means the rectangular prisms are interacting, while a weight 0 means
#'     they are not. Methods include "IDW" (inverse distance weighting), "rook"
#'     (adjacent rectangular prisms have a weight of 1, else 0), and "queen"
#'     (adjacent rectangular prisms, including diagonally, have a weight of 1,
#'     else 0).
#'
#' @return A numeric representing the spatial autocorrelation value.
#'
#' @examples
#' cell_prop_grid_metrics <- calculate_cell_proportion_grid_metrics3D(
#'     spe = SPIAT-3D::simulated_spe,
#'     n_splits = 10,
#'     reference_cell_types = c("Tumour"),
#'     target_cell_types = c("Immune"),
#'     feature_colname = "Cell.Type",
#'     plot_image = T
#' )
#' 
#' spatial_autocorrelation <- calculate_spatial_autocorrelation3D(
#'     grid_metrics = cell_prop_grid_metrics,
#'     metric_colname = "proportion",
#'     weight_method = "queen"
#' )
#' 
#' @export

calculate_spatial_autocorrelation3D <- function(grid_metrics,
                                                metric_colname,
                                                weight_method) {
  
  ## Check input parameters
  if (!(is.character(metric_colname))) {
    stop("`metric_colname` is not a character. This should be 'proportion' or 'entropy', depending on the chosen method.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  if (!(is.character(weight_method) && weight_method %in% c("IDW", "rook", "queen"))) {
    stop("`weight_method` is not 'IDW', 'rook' or 'queen'.")
  }
  
  ## Get number of grid prisms
  n_grid_prisms <- nrow(grid_metrics)
  
  ## Get splitting number (should be the cube root of n_grid_prisms)
  n_splits <- (n_grid_prisms)^(1/3)
  
  ## Find the coordinates of each grid prism
  x <- ((seq(n_grid_prisms) - 1) %% n_splits)
  y <- (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits))
  z <- (floor((seq(n_grid_prisms) - 1) / (n_splits^2)))
  grid_prism_coords <- data.frame(x = x, y = y, z = z)
  
  ## Subset for non NA rows
  grid_prism_coords <- grid_prism_coords[!is.na(grid_metrics[[metric_colname]]), ]
  grid_metrics <- grid_metrics[!is.na(grid_metrics[[metric_colname]]), ]
  
  weight_matrix <- -1 * apcluster::negDistMat(grid_prism_coords)
  ## Use the inverse distance between two points as the weight (IDW is 'inverse distance weighting')
  if (weight_method == "IDW") {
    weight_matrix <- 1 / weight_matrix
  }
  ## Use rook method: adjacent points get a weight of 1, otherwise, weight of 0
  ## Adjacent points are within 1 unit apart. e.g. (0, 0, 0) vs (0, 0, 1)
  else if (weight_method == "rook") {
    weight_matrix <- ifelse(weight_matrix > 1, 0, 1)  
  }
  ## Use queen method: adjacent points get a weight of 1, otherwise, weight of 0
  ## Adjacent points are within sqrt(3) unit apart. e.g. (0, 0, 0) vs (0, 0, 1)
  else if (weight_method == "queen") {
    weight_matrix <- ifelse(weight_matrix > sqrt(3), 0, 1)  
  }
  
  ## Points along the diagonal are comparing the same point so its weight is zero
  diag(weight_matrix) <- 0
  
  n <- nrow(grid_metrics)
  
  # Center the data
  data_centered <- grid_metrics[[metric_colname]] - mean(grid_metrics[[metric_colname]])
  
  # Calculate numerator using matrix multiplication
  numerator <- sum(data_centered * (weight_matrix %*% data_centered))
  
  # Calculate denominator
  denominator <- sum(data_centered^2) * sum(weight_matrix)
  
  # Moran's I
  I <- (n * numerator) / denominator
  
  return(I)
}
