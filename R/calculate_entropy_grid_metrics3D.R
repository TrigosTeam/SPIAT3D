calculate_entropy_grid_metrics3D <- function(spe, 
                                             n_splits,
                                             cell_types_of_interest,
                                             feature_colname = "Cell.Type",
                                             plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  ## Check cell_types_of_interest are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## If cell types have been chosen, check they are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, unique(spe[[feature_colname]]))
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Add grid metrics to spe
  spe <- get_spe_grid_metrics3D(spe, n_splits, feature_colname)
  
  # Get grid_prism_cell_matrix from spe
  grid_prism_cell_matrix <- spe@metadata$grid_metrics$grid_prism_cell_matrix
  
  ## Define data frame which contains all results
  n_grid_prisms <- n_splits^3
  result <- data.frame(row.names = seq(n_grid_prisms))
  
  for (cell_type in cell_types_of_interest) {
    result[[cell_type]] <- grid_prism_cell_matrix[[cell_type]]
  }
  result$total <- rowSums(result)
  
  ## Get data frame containing proportions for cell_types_of_interest
  df_props <- result[ , cell_types_of_interest] / result$total
  
  ## Use proportion data frame to get entropy
  calculate_entropy <- function(x) {
    entropy <- -1 * sum(x * ifelse(is.infinite(log(x, length(x))), 0, log(x, length(x))))
    return(entropy)
  }
  result$entropy <- apply(df_props, 1, calculate_entropy)
  
  # Add grid_prism_coordinates info to result
  result <- cbind(result, spe@metadata$grid_metrics$grid_prism_coordinates)
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous3D(result, "entropy")
    methods::show(fig)
  }
  
  return(result)
}
