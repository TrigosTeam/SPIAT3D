#' @title Calculate grid metrics in 3D spatial data.
#'
#' @description This functions divides a 3D SpatialExperiment Object into a 3D 
#'     grid of rectangular prisms. The number and types of cells in each 
#'     rectangular prism and the coordinate of each rectangular prism is 
#'     calculated. The is a supporting function used in 
#'     calculate_cell_proportion_grid_metrics3D and 
#'     calculate_entropy_grid_metrics3D.
#' 
#' @param spe A SpatialExperiment object containing 3D spatial information for 
#'     the cells. Naming of spatial coordinates MUST be "Cell.X.Position", 
#'     "Cell.Y.Position", "Cell.Z.Position" for the x-coordinate, y-coordinate 
#'     and z-coordinate of each cell. 
#' @param n_splits A positive numeric integer specifying the number splits used
#'     to divide the x-axis, y-axis and z-axis.
#' @param feature_colname A string specifying the name of the column in the 
#'     `colData` slot of the SpatialExperiment object that contains the cell 
#'     type information. Defaults to "Cell.Type".
#'
#' @return A list containing the number and type of cells in each rectangular 
#'     prism, and the coordinate of each rectangular prism.
#' 
#' @export

get_spe_grid_metrics3D <- function(spe, 
                                   n_splits, 
                                   feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  
  spe_coords <- spatialCoords(spe)
  
  ## Get dimensions of the window
  min_x <- min(spe_coords[ , "Cell.X.Position"])
  min_y <- min(spe_coords[ , "Cell.Y.Position"])
  min_z <- min(spe_coords[ , "Cell.Z.Position"])
  
  max_x <- max(spe_coords[ , "Cell.X.Position"])
  max_y <- max(spe_coords[ , "Cell.Y.Position"])
  max_z <- max(spe_coords[ , "Cell.Z.Position"])
  
  length <- round(max_x - min_x)
  width  <- round(max_y - min_y)
  height <- round(max_z - min_z)
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  d_lay <- height / n_splits
  
  # Shift spe_coords so they begin at the origin
  spe_coords[, "Cell.X.Position"] <- spe_coords[, "Cell.X.Position"] - min_x
  spe_coords[, "Cell.Y.Position"] <- spe_coords[, "Cell.Y.Position"] - min_y
  spe_coords[, "Cell.Z.Position"] <- spe_coords[, "Cell.Z.Position"] - min_z
  
  ## Figure out which 'grid prism number' each cell is inside
  spe$grid_prism_num <- floor(spe_coords[ , "Cell.X.Position"] / d_row) +
    floor(spe_coords[ , "Cell.Y.Position"] / d_col) * n_splits + 
    floor(spe_coords[ , "Cell.Z.Position"] / d_lay) * n_splits^2 + 1
  
  ## Determine the cell types found in each grid prism
  n_grid_prisms <- n_splits^3
  grid_prism_cell_matrix <- as.data.frame.matrix(table(spe[[feature_colname]], factor(spe$grid_prism_num, levels = seq(n_grid_prisms))))
  grid_prism_cell_matrix <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       t(grid_prism_cell_matrix), check.names = FALSE)
  
  ## Determine centre coordinates of each grid prism
  grid_prism_coordinates <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       x_coord = ((seq(n_grid_prisms) - 1) %% n_splits + 0.5) * d_row + round(min_x),
                                       y_coord = (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits) + 0.5) * d_col + round(min_y),
                                       z_coord = (floor((seq(n_grid_prisms) - 1) / (n_splits^2)) + 0.5) * d_lay + round(min_z))
  
  spe@metadata[["grid_metrics"]] <- list("grid_prism_cell_matrix" = grid_prism_cell_matrix,
                                         "grid_prism_coordinates" = grid_prism_coordinates)
  
  return(spe)
}
