calculate_all_gradient_cc_metrics3D <- function(spe, 
                                                reference_cell_type, 
                                                target_cell_types, 
                                                radii, 
                                                feature_colname = "Cell.Type", 
                                                plot_image = T) {
  

  ## Define result
  result <- list("mixing_score" = list(),
                 "cells_in_neighbourhood" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "cells_in_neighbourhood_proportion" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "entropy" = data.frame(matrix(nrow = length(radii), ncol = 1)),
                 "cross_K" = list())
  colnames(result[["cells_in_neighbourhood"]]) <- target_cell_types
  colnames(result[["cells_in_neighbourhood_proportion"]]) <- target_cell_types
  colnames(result[["entropy"]]) <- "entropy"

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
  
  # Define indiviudal data frames for mixing_score and cross_K
  for (target_cell_type in target_cell_types) {
    if (reference_cell_type != target_cell_type) {
      result[["mixing_score"]][[target_cell_type]] <- data.frame(matrix(nrow = length(radii), ncol = length(mixing_score_df_colnames)))
      colnames(result[["mixing_score"]][[target_cell_type]]) <- mixing_score_df_colnames
    }
    result[["cross_K"]][[target_cell_type]] <- data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames)))
    colnames(result[["cross_K"]][[target_cell_type]]) <- cross_K_df_colnames
  }
  
  # Get gradient results for each metric
  for (i in seq(length(radii))) {
    df <- calculate_all_single_radius_cc_metrics3D(spe,
                                                   reference_cell_type,
                                                   target_cell_types,
                                                   radii[i],
                                                   feature_colname)
    
    if (is.null(df)) return(NULL)
    
    df[["cells_in_neighbourhood"]]$ref_cell_id <- NULL
    
    result[["cells_in_neighbourhood"]][i, ] <- apply(df[["cells_in_neighbourhood"]], 2, mean)
    result[["cells_in_neighbourhood_proportion"]][i, ] <- apply(df[["cells_in_neighbourhood_proportion"]][ , paste(target_cell_types, "_prop", sep = "")], 2, mean, na.rm = T)
    result[["entropy"]][i, "entropy"] <- mean(df[["entropy"]]$entropy, na.rm = T)
    
    for (target_cell_type in names(df[["mixing_score"]])) {
      result[["mixing_score"]][[target_cell_type]][i, ] <- df[["mixing_score"]][[target_cell_type]]
    }
    
    for (target_cell_type in names(df[["cross_K"]])) {
      result[["cross_K"]][[target_cell_type]][i, ] <- df[["cross_K"]][[target_cell_type]]
    }
  }
  
  # Add radius column to each data frame
  result[["cells_in_neighbourhood"]]$radius <- radii
  result[["cells_in_neighbourhood_proportion"]]$radius <- radii
  result[["entropy"]]$radius <- radii
  for (target_cell_type in names(df[["mixing_score"]])) {
    result[["mixing_score"]][[target_cell_type]]$radius <- radii
  }
  
  for (target_cell_type in names(df[["cross_K"]])) {
    result[["cross_K"]][[target_cell_type]]$radius <- radii
  }
  
  
  ## Plot
  if (plot_image) {
    plot_cells_in_neighbourhood_gradient3D(result[["cells_in_neighbourhood"]], target_cell_types)
    plot_cells_in_neighbourhood_proportions_gradient3D(result[["cells_in_neighbourhood_proportion"]], target_cell_types)
    expected_entropy <- calculate_entropy_background3D(spe, target_cell_types, feature_colname)
    plot_entropy_gradient3D(result[["entropy"]], expected_entropy, reference_cell_type, target_cell_types)
    
    for (target_cell_type in names(df[["mixing_score"]])) {
      plot_mixing_scores_gradient3D(result[["mixing_score"]][[target_cell_type]])
    }
    
    for (target_cell_type in names(df[["cross_K"]])) {
      plot_cross_K_gradient3D(result[["cross_K"]][[target_cell_type]], reference_cell_type, target_cell_type)
    }
  }
  
  return(result)
}
