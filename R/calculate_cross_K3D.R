calculate_cross_K3D <- function(spe, 
                                reference_cell_type, 
                                target_cell_type, 
                                radius, 
                                feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!(is.character(reference_cell_type) && length(reference_cell_type) == 1)) {
    stop("`reference_cell_type` is not a character.")
  }
  if (!is.character(target_cell_types)) {
    stop("`target_cell_types` is not a character vector.")
  }
  if (!(is.numeric(radius) && length(radius) == 1 && radius > 0)) {
    stop("`radius` is not a positive numeric.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }  
  
  
  ## Get expected cross K-function
  expected_cross_K <- (4/3) * pi * radius^3
  
  # No reference or target cell in spe object: observed_cross_K is undefined.
  if (!(reference_cell_type %in% spe[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type, "is not found in the spe object"))
    
    result <- data.frame(observed_cross_K = NA,
                         expected_cross_K = expected_cross_K,
                         cross_K_ratio = NA)
    
    return(result)
  }
  if (!(target_cell_type %in% spe[[feature_colname]])) {
    warning(paste("The target_cell_type", target_cell_type, "is not found in the spe object"))
    
    result <- data.frame(observed_cross_K = NA,
                         expected_cross_K = expected_cross_K,
                         cross_K_ratio = NA)
    
    return(result)
  }
  
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood3D(spe,
                                                                  reference_cell_type,
                                                                  target_cell_type,
                                                                  radius,
                                                                  feature_colname,
                                                                  show_summary = FALSE,
                                                                  plot_image = FALSE)
  
  n_ref_tar_interactions <- sum(cells_in_neighbourhood_df[[target_cell_type]])
  n_ref_cells <- sum(spe[[feature_colname]] == reference_cell_type)
  n_tar_cells <- sum(spe[[feature_colname]] == target_cell_type)
  
  ## Get rough dimensions of the window the points are in
  spe_coords <- data.frame(spatialCoords(spe))
  
  length <- round(max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- round(max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  height <- round(max(spe_coords$Cell.Z.Position) - min(spe_coords$Cell.Z.Position))
  
  ## Get volume of the window the cells are in
  volume <- length * width * height
  
  ## Get observed cross K-function
  observed_cross_K <- (((volume * n_ref_tar_interactions) / n_ref_cells) / n_tar_cells)
  
  result <- data.frame(observed_cross_K = observed_cross_K,
                       expected_cross_K = expected_cross_K,
                       cross_K_ratio = observed_cross_K / expected_cross_K)
  
  return(result)
}
