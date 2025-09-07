calculate_all_gradient_cc_metrics3D <- function(spe, 
                                                reference_cell_type, 
                                                target_cell_types, 
                                                radii, 
                                                feature_colname = "Cell.Type", 
                                                plot_image = T) {
  
  # Define constants
  cross_K_df_colnames <- c("reference",
                           "expected",
                           target_cell_types)
  mixing_score_df_colnames <- c("ref_cell_type", 
                                "tar_cell_type", 
                                "n_ref_cells",
                                "n_tar_cells", 
                                "n_ref_tar_interactions",
                                "n_ref_ref_interactions", 
                                "mixing_score", 
                                "normalised_mixing_score")
  cross_G_df_colnames <- c("observed_cross_G",
                           "expected_cross_G")
  co_occurrence_df_colnames <- c("reference",
                                 target_cell_types)
  
  ## Define result
  result <- list("mixing_score" = list(),
                 "cells_in_neighbourhood" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "cells_in_neighbourhood_proportion" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "entropy" = data.frame(matrix(nrow = length(radii), ncol = 1)),
                 "cross_K" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_L" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_G" = list(),
                 "co_occurrence" = data.frame(matrix(nrow = length(radii), ncol = length(co_occurrence_df_colnames))))
  colnames(result[["cells_in_neighbourhood"]]) <- target_cell_types
  colnames(result[["cells_in_neighbourhood_proportion"]]) <- target_cell_types
  colnames(result[["entropy"]]) <- "entropy"
  colnames(result[["cross_K"]]) <- cross_K_df_colnames
  colnames(result[["cross_L"]]) <- cross_K_df_colnames
  colnames(result[["co_occurrence"]]) <- co_occurrence_df_colnames
  
  # Define individual data frames for mixing_score and cross_G
  for (target_cell_type in target_cell_types) {
    if (reference_cell_type != target_cell_type) {
      result[["mixing_score"]][[target_cell_type]] <- data.frame(matrix(nrow = length(radii), ncol = length(mixing_score_df_colnames)))
      colnames(result[["mixing_score"]][[target_cell_type]]) <- mixing_score_df_colnames
    }
    result[["cross_G"]][[target_cell_type]] <- data.frame(matrix(nrow = length(radii), ncol = length(cross_G_df_colnames)))
    colnames(result[["cross_G"]][[target_cell_type]]) <- cross_G_df_colnames
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
    result[["cross_K"]][i, ] <- df[["cross_K"]]
    result[["cross_L"]][i, ] <- df[["cross_L"]]
    result[["co_occurrence"]][i, ] <- df[["co_occurrence"]]
    
    for (target_cell_type in names(df[["mixing_score"]])) {
      result[["mixing_score"]][[target_cell_type]][i, ] <- df[["mixing_score"]][[target_cell_type]]
    }
    for (target_cell_type in names(df[["cross_G"]])) {
      result[["cross_G"]][[target_cell_type]][i, ] <- df[["cross_G"]][[target_cell_type]]
    }
  }
  
  # Add radius column to each data frame
  result[["cells_in_neighbourhood"]]$radius <- radii
  result[["cells_in_neighbourhood_proportion"]]$radius <- radii
  result[["entropy"]]$radius <- radii
  result[["cross_K"]]$radius <- radii
  result[["cross_L"]]$radius <- radii
  result[["co_occurrence"]]$radius <- radii
  for (target_cell_type in names(df[["mixing_score"]])) {
    result[["mixing_score"]][[target_cell_type]]$radius <- radii
  }
  for (target_cell_type in names(df[["cross_G"]])) {
    result[["cross_G"]][[target_cell_type]]$radius <- radii
  }
  
  
  ## Plot
  if (plot_image) {
    fig_ACIN <- plot_cells_in_neighbourhood_gradient3D(result[["cells_in_neighbourhood"]], reference_cell_type)
    methods::show(fig_ACIN)
    
    fig_ACINP <- plot_cells_in_neighbourhood_proportions_gradient3D(result[["cells_in_neighbourhood_proportion"]], reference_cell_type)
    methods::show(fig_ACINP)
    
    expected_entropy <- calculate_entropy_background3D(spe, target_cell_types, feature_colname)
    fig_AE <- plot_entropy_gradient3D(result[["entropy"]], expected_entropy, reference_cell_type, target_cell_types)
    methods::show(fig_AE)
    
    for (target_cell_type in names(result[["mixing_score"]])) {
      fig_NMS <- plot_mixing_scores_gradient3D(result[["mixing_score"]][[target_cell_type]], "NMS")
      fig_MS <- plot_mixing_scores_gradient3D(result[["mixing_score"]][[target_cell_type]], "MS")
      fig_NMS_MS <- plot_grid(fig_NMS, fig_MS, nrow = 2)
      methods::show(fig_NMS_MS)
    }
    fig_CK <- plot_cross_K_gradient3D(result[["cross_K"]])
    fig_CKR <- plot_cross_K_gradient_ratio3D(result[["cross_K"]])
    fig_CK_CKR <- plot_grid(fig_CK, fig_CKR, nrow = 2)
    methods::show(fig_CK_CKR)
    
    fig_CL <- plot_cross_L_gradient3D(result[["cross_L"]])
    fig_CLR <- plot_cross_L_gradient_ratio3D(result[["cross_L"]])
    fig_CL_CLR <- plot_grid(fig_CL, fig_CLR, nrow = 2)
    methods::show(fig_CL_CLR)
    
    for (target_cell_type in names(result[["cross_G"]])) {
      fig_CG <- plot_cross_G_gradient3D(result[["cross_G"]][[target_cell_type]], reference_cell_type, target_cell_type)
      methods::show(fig_CG)
    }
    
    fig_co_occ <- plot_co_occurrence_gradient3D(result[["co_occurrence"]])
    methods::show(fig_co_occ)
  }
  
  return(result)
}
