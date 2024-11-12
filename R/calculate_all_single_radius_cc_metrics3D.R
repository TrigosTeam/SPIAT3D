### Calculate all single radius cell-colocalisation metrics
# If a function only requires one target cell type, iterate through each cell type in target_cell_types, else use all target_cell_types

calculate_all_single_radius_cc_metrics3D <- function(spe, 
                                                     reference_cell_type, 
                                                     target_cell_types, 
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
    stop(paste(radius, " is not a positive numeric."))
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  
  ## For reference_cell_type, check it is found in the spe object
  if (!(reference_cell_type %in% spe[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spe object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spe object
  unknown_cell_types <- setdiff(target_cell_types, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  
  # Define result
  result <- list("cells_in_neighbourhood" = list(),
                 "cells_in_neighbourhood_proportion" = list(),
                 "entropy" = list(),
                 "mixing_score" = list(),
                 "cross_K" = list())
  
  # Define other constants
  mixing_score_df_colnames <- c("ref_cell_type", 
                                "tar_cell_type", 
                                "n_ref_cells",
                                "n_tar_cells", 
                                "n_ref_tar_interactions",
                                "n_ref_ref_interactions", 
                                "mixing_score", 
                                "normalised_mixing_score")
  cross_K_df_colnames <- c("ref_cell_type",
                           "tar_cell_type",
                           "observed_cross_K",
                           "expected_cross_K",
                           "cross_K_ratio")

  # Get rough dimensions of window for cross_K
  spe_coords <- data.frame(spatialCoords(spe))
  length <- round(max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- round(max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  height <- round(max(spe_coords$Cell.Z.Position) - min(spe_coords$Cell.Z.Position))
  ## Get volume of the window the cells are in
  volume <- length * width * height
  


  # All single radius cc metrics stem from calculate_entropy3D function
  entropy_df <- calculate_entropy3D(spe, 
                                    reference_cell_type, 
                                    target_cell_types, 
                                    radius, 
                                    feature_colname)  
  
  ## Cells in neighbourhood ----------
  result[["cells_in_neighbourhood"]] <- entropy_df[ , c("ref_cell_id", target_cell_types)]
  
  ## Cells in neighbourhood proportion ----------
  result[["cells_in_neighbourhood_proportion"]] <- entropy_df[ , c("ref_cell_id", target_cell_types, paste(target_cell_types, "_prop", sep = ""))]
  
  ## Entropy --------------
  result[["entropy"]] <- entropy_df
  
  
  ## These metrics focus on a particular cell type 
  for (target_cell_type in target_cell_types) {
    mixing_score_df <- data.frame(matrix(nrow = 1, ncol = length(mixing_score_df_colnames)))
    colnames(mixing_score_df) <- mixing_score_df_colnames
    mixing_score_df$ref_cell_type <- reference_cell_type
    
    cross_K_df <- data.frame(matrix(nrow = 1, ncol = length(cross_K_df_colnames)))
    colnames(cross_K_df) <- cross_K_df_colnames
    cross_K_df$ref_cell_type <- reference_cell_type
    
    ## Mixing score -----------------
    # No need to fill in mixing_score_df if the reference and target cell is the same
    if (reference_cell_type != target_cell_type) {
      mixing_score_df$tar_cell_type <- target_cell_type
      mixing_score_df$n_ref_cells <- sum(spe[[feature_colname]] == reference_cell_type)
      mixing_score_df$n_tar_cells <- sum(spe[[feature_colname]] == target_cell_type)
      mixing_score_df$n_ref_tar_interactions <- sum(entropy_df[[target_cell_type]])
      mixing_score_df$n_ref_ref_interactions <- sum(entropy_df[[reference_cell_type]])
      mixing_score_df$mixing_score <- mixing_score_df$n_ref_tar_interactions / (0.5 * mixing_score_df$n_ref_ref_interactions)
      mixing_score_df$normalised_mixing_score <- 0.5 * mixing_score_df$mixing_score * mixing_score_df$n_ref_cells / mixing_score_df$n_tar_cell
      if (is.infinite(mixing_score_df$mixing_score)) mixing_score_df$mixing_score <- NA
      if (is.infinite(mixing_score_df$normalised_mixing_score)) mixing_score_df$normalised_mixing_score <- NA
      result[["mixing_score"]][[target_cell_type]] <- mixing_score_df
    }
    
    ## Cross_K ---------------------
    cross_K_df$tar_cell_type <- target_cell_type
    cross_K_df$observed_cross_K <- (((volume * sum(entropy_df[[target_cell_type]])) / sum(df[[feature_colname]] == reference_cell_type)) / sum(df[[feature_colname]] == target_cell_type))
    cross_K_df$expected_cross_K <- (4/3) * pi * radius^3
    cross_K_df$cross_K_ratio <- cross_K_df$observed_cross_K / cross_K_df$expected_cross_K
    result[["cross_K"]][[target_cell_type]] <- cross_K_df
  }
  
  return(result)
}
