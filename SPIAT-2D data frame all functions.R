library(SpatialExperiment)
library(dbscan)

library(apcluster)
library(plotly)
library(dplyr)
library(reshape2)
library(gtools)
library(cowplot)
library(Hmisc)


calculate_all_gradient_cc_metrics2D <- function(spatial_df, 
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
                 "entropy" = data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types))),
                 "cross_K" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_L" = data.frame(matrix(nrow = length(radii), ncol = length(cross_K_df_colnames))),
                 "cross_G" = list(),
                 "co_occurrence" = data.frame(matrix(nrow = length(radii), ncol = length(co_occurrence_df_colnames))))
  colnames(result[["cells_in_neighbourhood"]]) <- target_cell_types
  colnames(result[["cells_in_neighbourhood_proportion"]]) <- target_cell_types
  colnames(result[["entropy"]]) <- target_cell_types
  colnames(result[["cross_K"]]) <- cross_K_df_colnames
  colnames(result[["cross_L"]]) <- cross_K_df_colnames
  colnames(result[["co_occurrence"]]) <- co_occurrence_df_colnames
  
  # Define indiviudal data frames for mixing_score and cross_G
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
    df <- calculate_all_single_radius_cc_metrics2D(spatial_df,
                                                   reference_cell_type,
                                                   target_cell_types,
                                                   radii[i],
                                                   feature_colname)
    
    if (is.null(df)) return(NULL)
    
    df[["cells_in_neighbourhood"]]$ref_cell_id <- NULL
    
    result[["cells_in_neighbourhood"]][i, ] <- apply(df[["cells_in_neighbourhood"]], 2, mean)
    result[["cells_in_neighbourhood_proportion"]][i, ] <- apply(df[["cells_in_neighbourhood_proportion"]][ , paste(target_cell_types, "_prop", sep = "")], 2, mean, na.rm = T)
    result[["entropy"]][i, ] <- apply(df[["entropy"]][ , paste(target_cell_types, "_entropy", sep = "")], 2, mean, na.rm = T)
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
    fig_ACIN <- plot_cells_in_neighbourhood_gradient2D(result[["cells_in_neighbourhood"]], reference_cell_type)
    methods::show(fig_ACIN)
    
    fig_ACINP <- plot_cells_in_neighbourhood_proportions_gradient2D(result[["cells_in_neighbourhood_proportion"]], reference_cell_type)
    methods::show(fig_ACINP)
    
    expected_entropy <- calculate_entropy_background2D(spatial_df, target_cell_types, feature_colname)
    fig_AE <- plot_entropy_gradient2D(result[["entropy"]], expected_entropy, reference_cell_type, target_cell_types)
    methods::show(fig_AE)
    
    for (target_cell_type in names(result[["mixing_score"]])) {
      fig_NMS <- plot_mixing_scores_gradient2D(result[["mixing_score"]][[target_cell_type]], "NMS")
      fig_MS <- plot_mixing_scores_gradient2D(result[["mixing_score"]][[target_cell_type]], "MS")
      fig_NMS_MS <- plot_grid(fig_NMS, fig_MS, nrow = 2)
      methods::show(fig_NMS_MS)
    }
    fig_CK <- plot_cross_K_gradient2D(result[["cross_K"]])
    fig_CKR <- plot_cross_K_gradient_ratio2D(result[["cross_K"]])
    fig_CK_CKR <- plot_grid(fig_CK, fig_CKR, nrow = 2)
    methods::show(fig_CK_CKR)
    
    fig_CL <- plot_cross_L_gradient2D(result[["cross_L"]])
    fig_CLR <- plot_cross_L_gradient_ratio2D(result[["cross_L"]])
    fig_CL_CLR <- plot_grid(fig_CL, fig_CLR, nrow = 2)
    methods::show(fig_CL_CLR)
    
    for (target_cell_type in names(result[["cross_G"]])) {
      fig_CG <- plot_cross_G_gradient2D(result[["cross_G"]][[target_cell_type]], reference_cell_type, target_cell_type)
      methods::show(fig_CG)
    }
    
    fig_co_occ <- plot_co_occurrence_gradient2D(result[["co_occurrence"]])
    methods::show(fig_co_occ)
  }
  
  return(result)
}
### Calculate all single radius cell-colocalisation metrics
# If a function only requires one target cell type, iterate through each cell type in target_cell_types, else use all target_cell_types

calculate_all_single_radius_cc_metrics2D <- function(spatial_df, 
                                                     reference_cell_type, 
                                                     target_cell_types, 
                                                     radius, 
                                                     feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
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
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  
  ## For reference_cell_type, check it is found in the spatial_df object
  if (!(reference_cell_type %in% spatial_df[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spatial_df object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spatial_df object
  unknown_cell_types <- setdiff(target_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
  }
  # Define result
  result <- list("cells_in_neighbourhood" = list(),
                 "cells_in_neighbourhood_proportion" = list(),
                 "entropy" = list(),
                 "mixing_score" = list(),
                 "cross_K" = list(),
                 "cross_L" = list(),
                 "cross_G" = list(),
                 "co_occurrence" = list())
  
  # Define other constants
  mixing_score_df_colnames <- c("ref_cell_type", 
                                "tar_cell_type", 
                                "n_ref_cells",
                                "n_tar_cells", 
                                "n_ref_tar_interactions",
                                "n_ref_ref_interactions", 
                                "mixing_score", 
                                "normalised_mixing_score")
  cross_K_df_colnames <- c("reference",
                           "expected",
                           target_cell_types)
  cross_G_df_colnames <- c("observed_cross_G",
                           "expected_cross_G")
  co_occurrence_df_colnames <- c("reference",
                                 target_cell_types)
  
  # Get rough dimensions of window for cross_K
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position")]
  length <- round(max(spatial_df_coords$Cell.X.Position) - min(spatial_df_coords$Cell.X.Position))
  width  <- round(max(spatial_df_coords$Cell.Y.Position) - min(spatial_df_coords$Cell.Y.Position))

  ## Get volume of the window the cells are in
  volume <- length * width
  
  # All single radius cc metrics stem from calculate_entropy2D function
  entropy_df <- calculate_entropy2D(spatial_df, 
                                    reference_cell_type, 
                                    target_cell_types, 
                                    radius, 
                                    feature_colname)  

  ## Cells in neighbourhood ----------
  result[["cells_in_neighbourhood"]] <- entropy_df[ , c("ref_cell_id", target_cell_types)]
  
  ## Cells in neighbourhood proportion ----------
  result[["cells_in_neighbourhood_proportion"]] <- entropy_df[ , c("ref_cell_id", paste(target_cell_types, "_prop", sep = ""))]
  
  ## Entropy --------------
  result[["entropy"]] <- entropy_df[ , c("ref_cell_id", paste(target_cell_types, "_entropy", sep = ""))]
  
  ## Mixing score -----------------
  for (target_cell_type in target_cell_types) {
    mixing_score_df <- data.frame(matrix(nrow = 1, ncol = length(mixing_score_df_colnames)))
    colnames(mixing_score_df) <- mixing_score_df_colnames
    mixing_score_df$ref_cell_type <- reference_cell_type
    
    # No need to fill in mixing_score_df if the reference and target cell is the same
    if (reference_cell_type != target_cell_type) {
      mixing_score_df$tar_cell_type <- target_cell_type
      mixing_score_df$n_ref_cells <- sum(spatial_df[[feature_colname]] == reference_cell_type)
      mixing_score_df$n_tar_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
      mixing_score_df$n_ref_tar_interactions <- sum(entropy_df[[target_cell_type]])
      mixing_score_df$n_ref_ref_interactions <- sum(entropy_df[[reference_cell_type]])
      mixing_score_df$mixing_score <- mixing_score_df$n_ref_tar_interactions / (0.5 * mixing_score_df$n_ref_ref_interactions)
      mixing_score_df$normalised_mixing_score <- 0.5 * mixing_score_df$mixing_score * mixing_score_df$n_ref_cells / mixing_score_df$n_tar_cell
      if (is.infinite(mixing_score_df$mixing_score)) mixing_score_df$mixing_score <- NA
      if (is.infinite(mixing_score_df$normalised_mixing_score)) mixing_score_df$normalised_mixing_score <- NA
      result[["mixing_score"]][[target_cell_type]] <- mixing_score_df
    }
  }
  
  ## Cross_K ---------------------
  cross_K_df <- data.frame(matrix(nrow = 1, ncol = length(cross_K_df_colnames)))
  colnames(cross_K_df) <- cross_K_df_colnames
  cross_K_df$reference <- reference_cell_type
  cross_K_df$expected <- pi * radius^2
  
  for (target_cell_type in target_cell_types) {
    cross_K_df[[target_cell_type]] <- (((volume * sum(entropy_df[[target_cell_type]])) / sum(spatial_df[[feature_colname]] == reference_cell_type)) / sum(spatial_df[[feature_colname]] == target_cell_type)) 
  }
  result[["cross_K"]] <- cross_K_df
  
  ## Cross_L ---------------------
  cross_L_df <- cross_K_df
  cross_L_df[ , c("expected", target_cell_types)] <- (cross_L_df[ , c("expected", target_cell_types)] / (pi)) ^ (1/2)
  result[["cross_L"]] <- cross_L_df
  
  ## Cross_G ---------------------
  for (target_cell_type in target_cell_types) {
    cross_G_df <- data.frame(matrix(nrow = 1, ncol = length(cross_G_df_colnames)))
    colnames(cross_G_df) <- cross_G_df_colnames
    
    reference_target_interactions <- entropy_df[[target_cell_type]]
    n_target_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
    target_cell_type_intensity <- n_target_cells / volume
    observed_cross_G <- sum(reference_target_interactions != 0) / length(reference_target_interactions)
    expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * pi * radius^2)
    
    cross_G_df$observed_cross_G <- observed_cross_G
    cross_G_df$expected_cross_G <- expected_cross_G
    result[["cross_G"]][[target_cell_type]] <- cross_G_df
  } 
  
  ## Co_occurrence ---------------
  all_cell_types <- unique(spatial_df[[feature_colname]])
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spatial_df,
                                                                  reference_cell_type,
                                                                  all_cell_types,
                                                                  radius,
                                                                  feature_colname,
                                                                  F,
                                                                  F)
  
  cells_in_neighbourhood_df$total <- rowSums(cells_in_neighbourhood_df[, -1], na.rm = TRUE)
  
  
  co_occurrence_df <- data.frame(matrix(nrow = 1, ncol = length(co_occurrence_df_colnames)))
  colnames(co_occurrence_df) <- co_occurrence_df_colnames
  co_occurrence_df$reference <- reference_cell_type
  
  n_cells_in_spe <- length(spatial_df[[feature_colname]])
  n_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_df$total)
  
  for (target_cell_type in target_cell_types) {
    n_target_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_df[[target_cell_type]])
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    n_target_cells_in_spe <- sum(spatial_df[[feature_colname]] == target_cell_type)
    target_cell_type_proportion_in_spe <- n_target_cells_in_spe / n_cells_in_spe
    target_cell_type_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_spe
    
    co_occurrence_df[[target_cell_type]] <- target_cell_type_co_occurrence
  }
  result[["co_occurrence"]] <- co_occurrence_df
  
  return(result)
}


calculate_cell_proportion_grid_metrics2D <- function(spatial_df, 
                                                     n_splits,
                                                     reference_cell_types,
                                                     target_cell_types,
                                                     feature_colname = "Cell.Type",
                                                     plot_image = TRUE) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  ## Check reference_cell_types are found in the spatial_df object
  unknown_cell_types <- setdiff(reference_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in reference_cell_types are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## Check target_cell_types are found in the spatial_df object
  unknown_cell_types <- setdiff(target_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  # Check if there is intersection between reference_cell_types and target_cell_types
  if (length(intersect(reference_cell_types, target_cell_types)) > 0) {
    stop("Cannot have same cells in both reference_cell_types and target_cell_types")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Add grid metrics to spatial_df
  grid_metrics <- get_spatial_df_grid_metrics2D(spatial_df, n_splits, feature_colname)
  
  # Get grid_prism_cell_matrix from spatial_df
  grid_prism_cell_matrix <- grid_metrics$grid_prism_cell_matrix
  
  ## Define data frame which contains all results
  n_grid_prisms <- n_splits^2
  result <- data.frame(row.names = seq(n_grid_prisms))
  
  # Fill in the result data frame
  if (length(reference_cell_types) == 1) {
    result$reference <- grid_prism_cell_matrix[[reference_cell_types]]
  }
  else {
    result$reference <- rowSums(grid_prism_cell_matrix[ , reference_cell_types])
  }
  if (length(target_cell_types) == 1) {
    result$target <- grid_prism_cell_matrix[[target_cell_types]]
  }
  else {
    result$target <- rowSums(grid_prism_cell_matrix[ , target_cell_types])
  }
  result$total <- result$reference + result$target
  result$proportion <- result$target / result$total
  
  # Add grid_prism_coordinates info to result
  result <- cbind(result, grid_metrics$grid_prism_coordinates)
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous2D(result, "proportion")
    methods::show(fig)
  }
  
  return(result)
}


calculate_cell_proportions2D <- function(spatial_df,
                                         cell_types_of_interest = NULL, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (ncol(spatial_df) == 0) {
    stop("No cells found for calculating cell proportions.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Creates frequency/bar plot of all cell types in the entire image
  cell_proportions <- data.frame(table(spatial_df[[feature_colname]]))
  names(cell_proportions) <- c("cell_type", 'frequency')
  
  # Only include cell types the user has chosen
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spatial_df object
    unknown_cell_types <- setdiff(cell_types_of_interest, cell_proportions$cell_type)
    if (length(unknown_cell_types) != 0) {
      stop(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                 paste(unknown_cell_types, collapse = ", ")))
    }
    
    # Subset for cell types chosen by user
    cell_proportions <- cell_proportions[(cell_proportions$cell_type %in% cell_types_of_interest), ]
    
  }
  
  # Get frequency total for all cells
  cell_type_frequency_total <- sum(cell_proportions$frequency)
  
  # Get proportions and percentages
  cell_proportions$proportion <- cell_proportions$frequency / cell_type_frequency_total
  cell_proportions$percentage <- cell_proportions$proportion * 100
  
  # Order the cell types by proportion (highest cell proportion is first)
  cell_proportions <- cell_proportions[rev(order(cell_proportions$proportion)), ]
  rownames(cell_proportions) <- seq(nrow(cell_proportions))
  
  # Plot
  if (plot_image) {
    
    labels <- paste(round(cell_proportions$percentage, 1), "%", sep = "")
    
    fig <- ggplot(cell_proportions, aes(x = factor(cell_type, cell_type), y = percentage, fill = cell_type)) +
      geom_bar(stat='identity') + 
      theme_bw() +
      labs(title="Cell proportions", x = "Cell type", y = "Percentage") +
      theme(plot.title = element_text(hjust = 0.5), 
            legend.position = "none") +
      geom_text(aes(label = labels), vjust = 0)
    
    methods::show(fig)
  }
  
  return(cell_proportions)
}
calculate_cells_in_neighbourhood_gradient2D <- function(spatial_df, 
                                                        reference_cell_type, 
                                                        target_cell_types, 
                                                        radii, 
                                                        feature_colname = "Cell.Type",
                                                        plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (i in seq(length(radii))) {
    cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spatial_df,
                                                                    reference_cell_type,
                                                                    target_cell_types,
                                                                    radii[i],
                                                                    feature_colname,
                                                                    FALSE,
                                                                    FALSE)
    
    if (is.null(cells_in_neighbourhood_df)) return(NULL)
    
    cells_in_neighbourhood_df$ref_cell_id <- NULL
    result[i, ] <- apply(cells_in_neighbourhood_df, 2, mean)
  }
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_cells_in_neighbourhood_gradient2D(result, reference_cell_type)
    methods::show(fig)
  }
  
  return(result)
}
calculate_cells_in_neighbourhood_proportions_gradient2D <- function(spatial_df, 
                                                                    reference_cell_type, 
                                                                    target_cell_types, 
                                                                    radii, 
                                                                    feature_colname = "Cell.Type",
                                                                    plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types)))
  colnames(result) <- target_cell_types
  
  for (i in seq(length(radii))) {
    cell_proportions_neighbourhood_proportions_df <- calculate_cells_in_neighbourhood_proportions2D(spatial_df,
                                                                                                    reference_cell_type,
                                                                                                    target_cell_types,
                                                                                                    radii[i],
                                                                                                    feature_colname)
    
    if (is.null(cell_proportions_neighbourhood_proportions_df)) return(NULL)
    
    result[i, ] <- apply(cell_proportions_neighbourhood_proportions_df[ , paste(target_cell_types, "_prop", sep = "")], 2, mean, na.rm = T)
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  # Plot
  if (plot_image) {
    fig <- plot_cells_in_neighbourhood_proportions_gradient2D(result, reference_cell_type)
    methods::show(fig)
  }
  
  return(result)
}
calculate_cells_in_neighbourhood_proportions2D <- function(spatial_df, 
                                                           reference_cell_type, 
                                                           target_cell_types, 
                                                           radius, 
                                                           feature_colname = "Cell.Type") {
  
  ## Get cells in neighbourhood df
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spatial_df,
                                                                  reference_cell_type,
                                                                  c(reference_cell_type, target_cell_types),
                                                                  radius,
                                                                  feature_colname,
                                                                  FALSE,
                                                                  FALSE)
  
  if (is.null(cells_in_neighbourhood_df)) return(NULL)
  
  cells_in_neighbourhood_df[ , paste(target_cell_types, "_prop", sep = "")] <- 
    cells_in_neighbourhood_df[ , target_cell_types] / (cells_in_neighbourhood_df[ , target_cell_types] + cells_in_neighbourhood_df[ , reference_cell_type])
  
  # If reference cell type is in target cell types, proportion should be 1
  if (reference_cell_type %in% target_cell_types) {
    cells_in_neighbourhood_df[cells_in_neighbourhood_df[[reference_cell_type]] != 0, paste(reference_cell_type, "_prop", sep = "")] <- 1
  }
  
  return(cells_in_neighbourhood_df)
}
calculate_cells_in_neighbourhood2D <- function(spatial_df, 
                                               reference_cell_type, 
                                               target_cell_types, 
                                               radius, 
                                               feature_colname = "Cell.Type",
                                               show_summary = TRUE,
                                               plot_image = TRUE) {
  
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
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
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  ## For reference_cell_type, check it is found in the spatial_df object
  if (!(reference_cell_type %in% spatial_df[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spatial_df object"))
    return(NULL)
  }
  ## For target_cell_types, check they are found in the spatial_df object
  unknown_cell_types <- setdiff(target_cell_types, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in target_cell_types are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
  }
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }  
  
  # Get spatial_df coords
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position")]
  
  # Get reference_cell_type coords
  reference_cell_type_coords <- spatial_df_coords[spatial_df[[feature_colname]] == reference_cell_type, ]
  
  result <- data.frame(ref_cell_id = spatial_df$Cell.ID[spatial_df[[feature_colname]] == reference_cell_type])
  
  for (target_cell_type in target_cell_types) {
    
    if (sum(spatial_df[[feature_colname]] == target_cell_type) == 0) {
      result[[target_cell_type]] <- NA
      next
    }
    
    ## Get target_cell_type coords
    target_cell_type_coords <- spatial_df_coords[spatial_df[[feature_colname]] == target_cell_type, ]
    
    ## Determine number of target cells spatial_dfcified distance for each reference cell
    ref_tar_result <- dbscan::frNN(target_cell_type_coords, 
                                   eps = radius,
                                   query = reference_cell_type_coords, 
                                   sort = FALSE)
    
    n_targets <- rapply(ref_tar_result$id, length)
    
    
    # Don't want to include the reference cell as one of the target cells
    if (reference_cell_type == target_cell_type) n_targets <- n_targets - 1
    
    ## Add to data frame
    result[[target_cell_type]] <- n_targets
  }
  
  ## Print summary
  if (show_summary) {
    print(summarise_cells_in_neighbourhood2D(result))    
  }
  
  ## Plot
  if (plot_image) {
    fig <- plot_cells_in_neighbourhood_violin2D(result, reference_cell_type)
    methods::show(fig)
  }
  
  return(result)
}

calculate_co_occurrence_gradient2D <- function(spatial_df, 
                                               reference_cell_type, 
                                               target_cell_types, 
                                               radii, 
                                               feature_colname = "Cell.Type",
                                               plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = length(target_cell_types) + 1))
  colnames(result) <- c("reference", target_cell_types)
  
  for (i in seq(length(radii))) {
    co_occurrence_df <- calculate_co_occurrence2D(spatial_df,
                                                  reference_cell_type,
                                                  target_cell_types,
                                                  radii[i],
                                                  feature_colname)
    
    result[i, ] <- co_occurrence_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_co_occurrence_gradient2D(result)
    methods::show(fig)
  }
  
  return(result)
}
calculate_co_occurrence2D <- function(spatial_df, 
                                      reference_cell_type, 
                                      target_cell_types, 
                                      radius, 
                                      feature_colname = "Cell.Type") {
  
  # Get all cell types in spatial_df
  all_cell_types <- unique(spatial_df[[feature_colname]])
  
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spatial_df,
                                                                  reference_cell_type,
                                                                  all_cell_types,
                                                                  radius,
                                                                  feature_colname,
                                                                  F,
                                                                  F)
  
  cells_in_neighbourhood_df$total <- rowSums(cells_in_neighbourhood_df[, -1], na.rm = TRUE)

  result <- data.frame(reference = reference_cell_type)
  
  # Get total number of cells in spatial_df
  n_cells_in_spatial_df <- length(spatial_df[[feature_colname]])
  
  # Get total number of cells in radius around reference cell type
  n_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_df$total)
  
  for (target_cell_type in target_cell_types) {
    
    # Get total number of target cells in radius around reference cell type
    n_target_cells_in_reference_cell_type_radius <- sum(cells_in_neighbourhood_df[[target_cell_type]])
    
    # Get proportion of target cells in radius around reference cell type
    target_cell_type_proportion_in_reference_cell_type_radius <- n_target_cells_in_reference_cell_type_radius / n_cells_in_reference_cell_type_radius
    
    # Get proportion of target cell type in spatial_df
    n_target_cells_in_spatial_df <- sum(spatial_df[[feature_colname]] == target_cell_type)
    target_cell_type_proportion_in_spatial_df <- n_target_cells_in_spatial_df / n_cells_in_spatial_df
    
    # Get co-occurence value for taget cell type
    target_cell_type_co_occurrence <- target_cell_type_proportion_in_reference_cell_type_radius / target_cell_type_proportion_in_spatial_df
    
    # Add to result data frame
    result[[target_cell_type]] <- target_cell_type_co_occurrence
  }
  
  return(result)
}
calculate_cross_G_gradient2D <- function(spatial_df, 
                                         reference_cell_type, 
                                         target_cell_type, 
                                         radii, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 2))
  colnames(result) <- c("observed_cross_G", 
                        "expected_cross_G")
  
  for (i in seq(length(radii))) {
    cross_G_df <- calculate_cross_G2D(spatial_df,
                                      reference_cell_type,
                                      target_cell_type,
                                      radii[i],
                                      feature_colname)
    
    result[i, ] <- cross_G_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig <- plot_cross_G_gradient2D(result, reference_cell_type, target_cell_type)
    methods::show(fig)
  }
  
  return(result)
}
calculate_cross_G2D <- function(spatial_df,
                                reference_cell_type,
                                target_cell_type,
                                radius,
                                feature_colname = "Cell.Type") {
  
  ### Calculate the observed cross_G
  # Get the number of target cells in the radius around each reference cell
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spatial_df,
                                                                  reference_cell_type,
                                                                  target_cell_type,
                                                                  radius,
                                                                  feature_colname,
                                                                  show_summary = FALSE,
                                                                  plot_image = FALSE)
  
  reference_target_interactions <- cells_in_neighbourhood_df[[target_cell_type]]
  
  # cross_G: essentially the proportion of reference cells with at least 1 target cell within the chosen radius.
  observed_cross_G <- sum(reference_target_interactions != 0) / length(reference_target_interactions)
  
  ### Calculate the expected cross_G
  # Get rough dimensions of the window the points are in
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position")]
  
  length <- round(max(spatial_df_coords$Cell.X.Position) - min(spatial_df_coords$Cell.X.Position))
  width  <- round(max(spatial_df_coords$Cell.Y.Position) - min(spatial_df_coords$Cell.Y.Position))
  
  # Get volume of the window the cells are in
  volume <- length * width
  
  # Get the number of target cells
  n_target_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
  
  # Get target_cell_type intensity (density)
  target_cell_type_intensity <- n_target_cells / volume
  
  # Apply formula
  expected_cross_G <- 1 - exp(-1 * target_cell_type_intensity * pi * radius^2)
  
  result <- data.frame(observed_cross_G = observed_cross_G,
                       expected_cross_G = expected_cross_G)
  
  return(result)
}
calculate_cross_K_gradient2D <- function(spatial_df, 
                                         reference_cell_type, 
                                         target_cell_types, 
                                         radii, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 2 + length(target_cell_types)))
  colnames(result) <- c("reference", "expected", target_cell_types)
  
  for (i in seq(length(radii))) {
    cross_K_df <- calculate_cross_K2D(spatial_df,
                                      reference_cell_type,
                                      target_cell_types,
                                      radii[i],
                                      feature_colname)
    
    result[i, ] <- cross_K_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig1 <- plot_cross_K_gradient2D(result)
    fig2 <- plot_cross_K_gradient_ratio2D(result)
    
    combined_fig <- plot_grid(fig1, fig2, nrow = 2)
    methods::show(combined_fig)
  }
  
  return(result)
}
calculate_cross_K2D <- function(spatial_df, 
                                reference_cell_type, 
                                target_cell_types, 
                                radius, 
                                feature_colname = "Cell.Type") {
  
  if (is.null(spatial_df[[feature_colname]])) stop(paste("No column called", feature_colname, "found in spatial_df object"))
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }  
  
  
  ## Get expected cross K-function
  expected_cross_K <- pi * radius^2
  
  ## For reference_cell_type, check it is found in the spatial_df object
  if (!(reference_cell_type %in% spatial_df[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spatial_df object"))
    result <- data.frame(observed_cross_K = NA,
                         expected_cross_K = expected_cross_K,
                         cross_K_ratio = NA)
    return(result)
  }
  
  ## Get rough dimensions of the window the points are in
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position")]
  
  length <- round(max(spatial_df_coords$Cell.X.Position) - min(spatial_df_coords$Cell.X.Position))
  width  <- round(max(spatial_df_coords$Cell.Y.Position) - min(spatial_df_coords$Cell.Y.Position))

  ## Get volume of the window the cells are in
  volume <- length * width
  
  # Number of reference cell types is constant
  n_ref_cells <- sum(spatial_df[[feature_colname]] == reference_cell_type)
  
  # Define result data frame
  result <- data.frame(reference = reference_cell_type, expected = expected_cross_K)
  
  cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spatial_df,
                                                                  reference_cell_type,
                                                                  target_cell_types,
                                                                  radius,
                                                                  feature_colname,
                                                                  show_summary = FALSE,
                                                                  plot_image = FALSE)
  
  for (target_cell_type in target_cell_types) {
    
    n_ref_tar_interactions <- sum(cells_in_neighbourhood_df[[target_cell_type]])
    
    n_tar_cells <- sum(spatial_df[[feature_colname]] == target_cell_type)
    
    ## Get observed cross K-function
    if (n_tar_cells == 0) {
      observed_cross_K <- NA
    }
    else {
      observed_cross_K <- (volume * n_ref_tar_interactions) / (n_ref_cells * n_tar_cells)  
    }
    result[[target_cell_type]] <- observed_cross_K
  }
  
  return(result)
}
calculate_cross_L_gradient2D <- function(spatial_df, 
                                         reference_cell_type, 
                                         target_cell_types, 
                                         radii, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 2 + length(target_cell_types)))
  colnames(result) <- c("reference", "expected", target_cell_types)
  
  for (i in seq(length(radii))) {
    cross_L_df <- calculate_cross_L2D(spatial_df,
                                      reference_cell_type,
                                      target_cell_types,
                                      radii[i],
                                      feature_colname)
    
    result[i, ] <- cross_L_df
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig1 <- plot_cross_L_gradient2D(result)
    fig2 <- plot_cross_L_gradient_ratio2D(result)
    
    combined_fig <- plot_grid(fig1, fig2, nrow = 2)
    methods::show(combined_fig)
  }
  
  return(result)
}
calculate_cross_L2D <- function(spatial_df, 
                                reference_cell_type, 
                                target_cell_types, 
                                radius, 
                                feature_colname = "Cell.Type") {
  
  result <- calculate_cross_K2D(spatial_df = spatial_df,
                                reference_cell_type = reference_cell_type,
                                target_cell_types = target_cell_types,
                                radius = radius,
                                feature_colname = feature_colname)
  
  result[ , c("expected", target_cell_types)] <- (result[ , c("expected", target_cell_types)] / (pi)) ^ (1/2)
  
  return(result)
}
calculate_entropy_background2D <- function(spatial_df,
                                           cell_types_of_interest, 
                                           feature_colname = "Cell.Type") {
  
  # NULL case: entropy is undefined
  if (is.null(cell_types_of_interest)) return(NA)
  
  # One cell type case: entropy is 0
  if (is.character(cell_types_of_interest) && length(cell_types_of_interest) == 1) return(0)
  
  cell_proportions_data <- calculate_cell_proportions2D(spatial_df, cell_types_of_interest, feature_colname, FALSE)
  
  # Calculate entropy of the entire image
  entropy <- -1 * sum(cell_proportions_data$proportion * log(cell_proportions_data$proportion, length(cell_proportions_data$proportion)))
  
  return(entropy) 
}
calculate_entropy_gradient2D <- function(spatial_df,
                                         reference_cell_type,
                                         target_cell_types,
                                         radii,
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 1))
  colnames(result) <- "entropy"
  
  for (i in seq(length(radii))) {
    entropy_df <- calculate_entropy2D(spatial_df,
                                      reference_cell_type,
                                      target_cell_types,
                                      radii[i],
                                      feature_colname)
    
    if (is.null(entropy_df)) return(NULL)
    
    result[i, "entropy"] <- mean(entropy_df$entropy, na.rm = T)
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    expected_entropy <- calculate_entropy_background2D(spatial_df, target_cell_types, feature_colname)
    fig <- plot_entropy_gradient2D(result, expected_entropy, reference_cell_type, target_cell_types)
    methods::show(fig)
  }
  
  return(result)
}
calculate_entropy_grid_metrics2D <- function(spatial_df, 
                                             n_splits,
                                             cell_types_of_interest,
                                             feature_colname = "Cell.Type",
                                             plot_image = TRUE) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  ## Check cell_types_of_interest are found in the spatial_df object
  unknown_cell_types <- setdiff(cell_types_of_interest, spatial_df[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  ## If cell types have been chosen, check they are found in the spatial_df object
  unknown_cell_types <- setdiff(cell_types_of_interest, unique(spatial_df[[feature_colname]]))
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    return(NULL)
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Add grid metrics to spatial_df
  grid_metrics <- get_spatial_df_grid_metrics2D(spatial_df, n_splits, feature_colname)
  
  # Get grid_prism_cell_matrix from spatial_df
  grid_prism_cell_matrix <- grid_metrics$grid_prism_cell_matrix
  
  ## Define data frame which contains all results
  n_grid_prisms <- n_splits^2
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
  result <- cbind(result, grid_metrics$grid_prism_coordinates)
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_metrics_continuous2D(result, "entropy")
    methods::show(fig)
  }
  
  return(result)
}
calculate_entropy2D <- function(spatial_df,
                                reference_cell_type,
                                target_cell_types,
                                radius,
                                feature_colname = "Cell.Type") {
  
  # Check target_cell_types
  if (!(is.character(target_cell_types) && length(target_cell_types) >= 2)) {
    stop("`target_cell_types` is not a character vector with at least 2 cell types.")
  }
  
  ## Users should ensure include the reference_cell_type as one of the target_cell_types
  cells_in_neighbourhood_proportion_df <- calculate_cells_in_neighbourhood_proportions2D(spatial_df,
                                                                                         reference_cell_type,
                                                                                         target_cell_types,
                                                                                         radius,
                                                                                         feature_colname)
  
  if (is.null(cells_in_neighbourhood_proportion_df)) return(NULL)

  ## Get entropy for target_cell_type
  cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_entropy", sep = "")] <- 
    -1 * 
    (cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_prop", sep = "")] * log(cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_prop", sep = "")], 2) +
    (1 - cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_prop", sep = "")]) * log(1 - (cells_in_neighbourhood_proportion_df[ , paste(target_cell_types, "_prop", sep = "")]), 2))
  
  return(cells_in_neighbourhood_proportion_df)
}


calculate_minimum_distances_between_cell_types2D <- function(spatial_df,
                                                             cell_types_of_interest = NULL,
                                                             feature_colname = "Cell.Type",
                                                             show_summary = TRUE,
                                                             plot_image = TRUE) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (ncol(spatial_df) < 2) {
    stop("There must be at least two cells in spatial_df.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }  
  
  # De-factor feature column in spatial_df object
  spatial_df[[feature_colname]] <- as.character(spatial_df[[feature_colname]])
  
  # Subset spatial_df to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spatial_df object
    unknown_cell_types <- setdiff(cell_types_of_interest, spatial_df[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    spatial_df <- spatial_df[spatial_df[[feature_colname]] %in% cell_types_of_interest, ]
  }
  # If cell_types_of_interest is NULL, use all cells in spatial_df
  else {
    cell_types_of_interest <- unique(spatial_df[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(spatial_df$Cell.ID[spatial_df[[feature_colname]] == cell_type])
  }
  
  # Get spatial_df coords
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position")]
  
  # Get different possible cell type combinations
  # Each row represents a combination
  # If a row is [1 , 2], then we are comparing cell type 1 and cell type 2
  permu <- gtools::permutations(length(cell_types_of_interest), 2, repeats.allowed = TRUE)
  
  result <- data.frame()
  
  for (i in seq(nrow(permu))) {
    cell_type1 <- cell_types_of_interest[permu[i, 1]]
    cell_type2 <- cell_types_of_interest[permu[i, 2]]
    
    # Don't have one of the cells
    if (sum(spatial_df[[feature_colname]] == cell_type1) == 0 || sum(spatial_df[[feature_colname]] == cell_type2) == 0) {
      result <- rbind(result, data.frame(ref_cell_id = NA, ref_cell_type = cell_type1, nearest_cell_id = NA, nearest_cell_type = cell_type2, distance = NA))
      next
    }
    
    # Get x, y, z coords for all cells of cell_type1 and cell_type2
    cell_type1_coords <- spatial_df_coords[spatial_df[[feature_colname]] == cell_type1, ]
    cell_type2_coords <- spatial_df_coords[spatial_df[[feature_colname]] == cell_type2, ]
    
    # Find all of closest points
    # For each cell of cell_type1, find the closest cell of cell_type2
    if (cell_type1 != cell_type2) {
      nearest_neighbours <- RANN::nn2(data = cell_type2_coords, 
                                      query = cell_type1_coords, 
                                      k = 1)  
    }
    # If we are comparing the same cell_type, and there is only one of this cell type, move on
    else if (nrow(cell_type1_coords) == 1) {
      warning("There is only 1 '", cell_type1, "' cell in your data. It has no nearest neighbour of the same cell type.", sep = "")
      result <- rbind(result, data.frame(ref_cell_id = NA, ref_cell_type = cell_type1, nearest_cell_id = NA, nearest_cell_type = cell_type2, distance = NA))
      next
    }
    # If we are comparing the same cell_type, use the second closest neighbour
    else {
      nearest_neighbours <- RANN::nn2(data = cell_type2_coords, 
                                      query = cell_type1_coords, 
                                      k = 2)
      nearest_neighbours[['nn.idx']] <- nearest_neighbours[['nn.idx']][ , 2]
      nearest_neighbours[['nn.dists']] <- nearest_neighbours[['nn.dists']][ , 2]
    }
    
    # Create the data frame containing the chosen cells and their ids, as well as the nearest cell to them and their ids, and the distance between
    
    df <- data.frame(
      ref_cell_id = cell_type_ids[[cell_type1]],
      ref_cell_type = cell_type1,
      nearest_cell_id = cell_type_ids[[cell_type2]][c(nearest_neighbours$nn.idx)],
      nearest_cell_type = cell_type2,
      distance = nearest_neighbours$nn.dists
    )
    result <- rbind(result, df)
  }
  
  result$pair <- paste(result$ref_cell_type, result$nearest_cell_type,sep = "/")
  
  # Print summary
  if (show_summary) {
    print(summarise_distances_between_cell_types2D(result))  
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_distances_between_cell_types_violin2D(result)
    methods::show(fig)
  }
  
  return(result)
}

calculate_mixing_scores_gradient2D <- function(spatial_df, 
                                               reference_cell_type, 
                                               target_cell_type, 
                                               radii, 
                                               feature_colname = "Cell.Type",
                                               plot_image = TRUE) {
  
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }
  
  result <- data.frame(matrix(nrow = length(radii), ncol = 8))
  colnames(result) <- c("ref_cell_type", 
                        "tar_cell_type", 
                        "n_ref_cells",
                        "n_tar_cells", 
                        "n_ref_tar_interactions",
                        "n_ref_ref_interactions", 
                        "mixing_score", 
                        "normalised_mixing_score")
  
  for (i in seq(length(radii))) {
    mixing_scores <- calculate_mixing_scores2D(spatial_df,
                                               reference_cell_type,
                                               target_cell_type,
                                               radii[i],
                                               feature_colname)
    
    result[i, ] <- mixing_scores
  }
  
  # Add a radius column to the result
  result$radius <- radii
  
  if (plot_image) {
    fig1 <- plot_mixing_scores_gradient2D(result, "NMS")
    fig2 <- plot_mixing_scores_gradient2D(result, "MS")
    combined_fig <- plot_grid(fig1, fig2, nrow = 2)
    methods::show(combined_fig)
  }
  
  return(result)
}
calculate_mixing_scores2D <- function(spatial_df, 
                                      reference_cell_types, 
                                      target_cell_types, 
                                      radius, 
                                      feature_colname = "Cell.Type") {
  
  # Define result
  result <- data.frame()
  
  for (reference_cell_type in reference_cell_types) {
    
    for (target_cell_type in target_cell_types) {
      
      # No point getting mixing scores if comparing the same cell type
      if (reference_cell_type == target_cell_type) {
        next
      }
      
      # Get number of reference cells and target cells
      n_ref <- sum(spatial_df[[feature_colname]] == reference_cell_type)
      n_tar <- sum(spatial_df[[feature_colname]] == target_cell_type)
      
      
      # Can't get mixing scores if there are 0 or 1 reference cells
      if (n_ref == 0 || n_ref == 1) {
        result <-  rbind(result, 
                         c(reference_cell_type, 
                           target_cell_type, 
                           n_ref, 
                           n_tar, 
                           0, 
                           0, 
                           NA, 
                           NA))
      }
      
      
      ## Get cells in neighbourhood df
      cells_in_neighbourhood_df <- calculate_cells_in_neighbourhood2D(spatial_df,
                                                                      reference_cell_type,
                                                                      c(reference_cell_type, target_cell_type),
                                                                      radius,
                                                                      feature_colname,
                                                                      FALSE,
                                                                      FALSE)
      
      # Get number of ref-ref interactions
      # Halve it to avoid counting each ref-ref interaction twice
      n_ref_ref_interactions <- 0.5 * sum(cells_in_neighbourhood_df[[reference_cell_type]]) 
      
      # Get number of ref-tar interactions
      n_ref_tar_interactions <- sum(cells_in_neighbourhood_df[[target_cell_type]]) 
      
      
      # Can't get mixing scores if there are no target cells
      if (n_tar == 0) {
        
        result <-  rbind(result, 
                         c(reference_cell_type, 
                           target_cell_type, 
                           n_ref, 
                           0, 
                           0, 
                           n_ref_ref_interactions, 
                           NA, 
                           NA))
      }
      
      # Generic case: We have reference cells and target cells
      else {
        
        if (n_ref_ref_interactions != 0) {
          mixing_score <- n_ref_tar_interactions / n_ref_ref_interactions
          normalised_mixing_score <- 0.5 * mixing_score * n_ref / n_tar
        }
        else {
          mixing_score <- 0
          normalised_mixing_score <- 0
          methods::show(paste("There are no reference to reference interactions for", target_cell_type, "in the spatial_dfcified radius, cannot calculate mixing score"))
        }
        
        result <-  rbind(result, 
                         c(reference_cell_type, 
                           target_cell_type, 
                           n_ref, 
                           n_tar, 
                           n_ref_tar_interactions, 
                           n_ref_ref_interactions, 
                           mixing_score, 
                           normalised_mixing_score))
      }
    }
  }
  
  # Required column names of our output data frame
  colnames(result) <- c("ref_cell_type", 
                        "tar_cell_type", 
                        "n_ref_cells",
                        "n_tar_cells", 
                        "n_ref_tar_interactions",
                        "n_ref_ref_interactions", 
                        "mixing_score", 
                        "normalised_mixing_score")
  
  # Turn numeric data into numeric type
  result[ , 3:8] <- apply(result[ , 3:8], 2, as.numeric)
  
  return(result)
}
calculate_pairwise_distances_between_cell_types2D <- function(spatial_df,
                                                              cell_types_of_interest = NULL,
                                                              feature_colname = "Cell.Type",
                                                              show_summary = TRUE,
                                                              plot_image = TRUE) {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (ncol(spatial_df) < 2) {
    stop("There must be at least two cells in spatial_df.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  if (is.null(spatial_df[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spatial_df")
    spatial_df$Cell.ID <- paste("Cell", seq(nrow(spatial_df)), sep = "_")
  }
  
  # De-factor feature column in spatial_df object
  spatial_df[[feature_colname]] <- as.character(spatial_df[[feature_colname]])
  
  # Subset spatial_df to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spatial_df object
    unknown_cell_types <- setdiff(cell_types_of_interest, spatial_df[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spatial_df object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    spatial_df <- spatial_df[spatial_df[[feature_colname]] %in% cell_types_of_interest, ]
  }
  # If cell_types_of_interest is NULL, use all cells in spatial_df
  else {
    cell_types_of_interest <- unique(spatial_df[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(spatial_df$Cell.ID[spatial_df[[feature_colname]] == cell_type])
  }
  
  # Calculate cell to cell distances
  distance_matrix <- -1 * apcluster::negDistMat(spatial_df[ , c("Cell.X.Position", "Cell.Y.Position")])
  rownames(distance_matrix) <- spatial_df$Cell.ID
  colnames(distance_matrix) <- spatial_df$Cell.ID
  
  result <- data.frame()
  
  for (i in seq(length(cell_types_of_interest))) {
    
    for (j in i:length(cell_types_of_interest)) {
      
      # Get current cell types and cell ids
      cell_type1 <- names(cell_type_ids)[i]
      cell_type2 <- names(cell_type_ids)[j]
      
      cell_type1_ids <- cell_type_ids[[cell_type1]]
      cell_type2_ids <- cell_type_ids[[cell_type2]]
      
      ## Don't have a cell type, or the same cell type with only one cell
      if (length(cell_type1_ids) == 0 || length(cell_type2_ids) == 0) {
        result <- rbind(result, data.frame(Var1 = NA, Var2 = NA, value = NA, cell_type1 = cell_type1, cell_type2 = cell_type2, pair = paste(cell_type1, cell_type2, sep="/")))
        next
      }
      
      ## Same cell type only one cell
      if (cell_type1 == cell_type2 && length(cell_type1_ids) == 1) {
        warning("There is only 1 '", cell_type1, "' cell in your data. It has no pair of the same cell type.", sep = "")
        result <- rbind(result, data.frame(Var1 = NA, Var2 = NA, value = NA, cell_type1 = cell_type1, cell_type2 = cell_type2, pair = paste(cell_type1, cell_type2, sep="/")))
        next
      }
      
      # Subset distance_matrix for current cell types
      distance_matrix_subset <- distance_matrix[rownames(distance_matrix) %in% cell_type1_ids, 
                                                colnames(distance_matrix) %in% cell_type2_ids]
      
      ## Different cell types, each only has one cell
      if (length(cell_type1_ids) == 1 && length(cell_type2_ids) == 1) {
        distance_matrix_subset <- as.matrix(distance_matrix_subset)
        rownames(distance_matrix_subset) <- cell_type1_ids
        colnames(distance_matrix_subset) <- cell_type2_ids
      }    
      ## Different cell types, only one cell of cell_type1
      else if (length(cell_type1_ids) == 1) {
        distance_matrix_subset <- as.matrix(distance_matrix_subset)
        colnames(distance_matrix_subset) <- cell_type1_ids
      }
      ## Different cell types, only one cell of cell_type2
      else if (length(cell_type2_ids) == 1) {
        distance_matrix_subset <- as.matrix(distance_matrix_subset)
        colnames(distance_matrix_subset) <- cell_type2_ids
      }
      ## Same cell type, only need part of the matrix (make irrelevant part of matrix equal to NA)
      if (cell_type1 == cell_type2) distance_matrix_subset[upper.tri(distance_matrix_subset, diag = TRUE)] <- NA
      
      # Convert distance_matrix_subset to a data frame
      df <- reshape2::melt(distance_matrix_subset, na.rm = TRUE)
      df$cell_type1 <- cell_type1
      df$cell_type2 <- cell_type2
      df$pair <- paste(cell_type1, cell_type2, sep="/")
      
      result <- rbind(result, df)
    }
  }
  
  # Rearrange columns 
  colnames(result)[c(1, 2, 3)] <- c("cell_type1_id", "cell_type2_id", "distance")
  result <- result[ , c("cell_type1_id", "cell_type1", "cell_type2_id", "cell_type2", "distance", "pair")]
  
  # Print summary
  if (show_summary) {
    print(summarise_distances_between_cell_types2D(result))  
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_distances_between_cell_types_violin2D(result)
    methods::show(fig)
  }
  
  return(result)
}
calculate_prevalence_gradient_AUC2D <- function(prevalence_gradient_df) {
  
  return(sum(prevalence_gradient_df$prevalence) * 0.01)
}
calculate_prevalence_gradient2D <- function(grid_metrics,
                                            metric_colname,
                                            show_AUC = T,
                                            plot_image = T) {
  
  ## Check input parameters
  if (!(is.character(metric_colname))) {
    stop("`metric_colname` is not a character. This should be 'proportion' or 'entropy', depending on the chosen method.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  if (!is.logical(show_AUC)) {
    stop("`show_AUC` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  # Thresholds range from 0 to 1
  thresholds <- seq(0.01, 1, 0.01)
  
  # Define result
  result <- data.frame(threshold = thresholds)
  
  # Get prevalences for each threshold
  result$prevalence <- sapply(thresholds, function(threshold) { 
    calculate_prevalence2D(grid_metrics, metric_colname, threshold) 
  })
  
  # Show AUC of prevalence gradient graph
  if (show_AUC) {
    print(paste("AUC:", round(calculate_prevalence_gradient_AUC2D(result), 2)))
  }
  
  # Plot
  if (plot_image) {
    fig <- ggplot(result, aes(threshold, prevalence)) +
      geom_line() +
      theme_bw() +
      labs(x = "Threshold",
           y = "Prevalence",
           title = paste("Prevalence vs Threshold (", metric_colname, ")", sep = "")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      ylim(0, 100)
    methods::show(fig)
  }
  
  return(result)
}
calculate_prevalence2D <- function(grid_metrics,
                                   metric_colname,
                                   threshold,
                                   above = TRUE) {
  
  ## Check input parameters
  if (!(is.character(metric_colname))) {
    stop("`metric_colname` is not a character. This should be 'proportion' or 'entropy', depending on the chosen method.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  if (!(is.numeric(threshold) && length(threshold) == 1 && threshold >= 0 && threshold <= 1)) {
    stop("`threshold` is not a numeric between 0 and 1.")
  }
  if (!is.logical(above)) {
    stop("`above` is not a logical (TRUE or FALSE).")
  }
  
  
  ## Exclude rows with NA values
  grid_metrics <- grid_metrics[!is.na(grid_metrics[[metric_colname]]), ]
  
  if (above) {
    prevalence <- sum(grid_metrics[[metric_colname]] >= threshold) / nrow(grid_metrics) * 100
  }
  else {
    prevalence <- sum(grid_metrics[[metric_colname]] < threshold) / nrow(grid_metrics) * 100    
  }
  
  return(prevalence)
}
calculate_spatial_autocorrelation2D <- function(grid_metrics,
                                                metric_colname,
                                                weight_method = 0.1) {
  
  ## Check input parameters
  if (!(is.character(metric_colname))) {
    stop("`metric_colname` is not a character. This should be 'proportion' or 'entropy', depending on the chosen method.")
  }
  if (is.null(grid_metrics[[metric_colname]])) {
    stop("`metric_colname` is not a column in `grid_metrics`.")
  }
  if (!((is.numeric(weight_method) && length(weight_method) == 1 && weight_method > 0 && weight_method < 1) ||
        (is.character(weight_method) && weight_method %in% c("IDW", "rook", "queen")))) {
    stop("`weight_method` is not a numeric between 0 and 1 or either 'IDW', 'rook' or 'queen'.")
  }
  
  ## Get number of grid prisms
  n_grid_prisms <- nrow(grid_metrics)
  
  ## Get splitting number (should be the cube root of n_grid_prisms)
  n_splits <- (n_grid_prisms)^(1/3)
  
  ## Find the coordinates of each grid prism
  x <- ((seq(n_grid_prisms) - 1) %% n_splits)
  y <- (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits))
  grid_prism_coords <- data.frame(x = x, y = y)
  
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
  ## Adjacent points are within sqrt(2) unit apart. e.g. (0, 0) vs (0, 1)
  else if (weight_method == "queen") {
    weight_matrix <- ifelse(weight_matrix > sqrt(2), 0, 1)  
  }
  ## If a number (x) between 0 and 1 is supplied, set a threshold to be x quantile value of c(weight_matrix)
  ## Grid prisms within this spatial_dfcified threshold have a weight of 1, otherwise, weight of 0
  else if (as.numeric(weight_method) && 0 < weight_method && weight_method < 1) {
    threshold <- quantile(c(weight_matrix), weight_method)
    weight_matrix <- ifelse(weight_matrix > threshold, 0, 1)
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

get_spatial_df_grid_metrics2D <- function(spatial_df, 
                                          n_splits, 
                                          feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spatial_df) != "data.frame") {
    stop("`spatial_df` is not a data.frame object.")
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spatial_df[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spatial_df object."))
  }
  
  spatial_df_coords <- spatial_df[ , c("Cell.X.Position", "Cell.Y.Position")]
  
  ## Get dimensions of the window
  min_x <- min(spatial_df_coords[ , "Cell.X.Position"])
  min_y <- min(spatial_df_coords[ , "Cell.Y.Position"])
  
  max_x <- max(spatial_df_coords[ , "Cell.X.Position"])
  max_y <- max(spatial_df_coords[ , "Cell.Y.Position"])
  
  length <- round(max_x - min_x)
  width  <- round(max_y - min_y)
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  
  # Shift spatial_df_coords so they begin at the origin
  spatial_df_coords[, "Cell.X.Position"] <- spatial_df_coords[, "Cell.X.Position"] - min_x
  spatial_df_coords[, "Cell.Y.Position"] <- spatial_df_coords[, "Cell.Y.Position"] - min_y
  
  ## Figure out which 'grid prism number' each cell is inside
  spatial_df$grid_prism_num <- floor(spatial_df_coords[ , "Cell.X.Position"] / d_row) +
    floor(spatial_df_coords[ , "Cell.Y.Position"] / d_col) * n_splits
  
  ## Determine the cell types found in each grid prism
  n_grid_prisms <- n_splits^2
  grid_prism_cell_matrix <- as.data.frame.matrix(table(spatial_df[[feature_colname]], factor(spatial_df$grid_prism_num, levels = seq(n_grid_prisms))))
  grid_prism_cell_matrix <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       t(grid_prism_cell_matrix), check.names = FALSE)
  
  ## Determine centre coordinates of each grid prism
  grid_prism_coordinates <- data.frame(grid_prism_num = seq(n_grid_prisms),
                                       x_coord = ((seq(n_grid_prisms) - 1) %% n_splits + 0.5) * d_row + round(min_x),
                                       y_coord = (floor(((seq(n_grid_prisms) - 1) %% (n_splits)^2) / n_splits) + 0.5) * d_col + round(min_y))
  
  grid_metrics <- list("grid_prism_cell_matrix" = grid_prism_cell_matrix,
                       "grid_prism_coordinates" = grid_prism_coordinates)
  
  return(grid_metrics)
}



summarise_cells_in_neighbourhood2D <- function(cells_in_neighbourhood_df) {
  
  ## Target cell types will be all the columns except the first column
  target_cell_types <- colnames(cells_in_neighbourhood_df)[c(-1)]
  
  ## Set up data frame for summarised_results list
  df <- data.frame(row.names = c("mean", "min", "max", "median", "st_dev"))
  
  for (target_cell_type in target_cell_types) {
    
    ## Get statistical measures for each target cell type
    target_cell_type_values <- cells_in_neighbourhood_df[[target_cell_type]]
    df[[target_cell_type]] <- c(mean(target_cell_type_values),
                                min(target_cell_type_values),
                                max(target_cell_type_values),
                                median(target_cell_type_values),
                                sd(target_cell_type_values))
    
  }
  
  return(data.frame(t(df)))
}
summarise_distances_between_cell_types2D <- function(distances_df) {
  
  pair <- distance <- NULL
  
  # summarise the results
  distances_df_summarised <- distances_df %>% 
    dplyr::group_by(pair) %>%
    dplyr::summarise(mean(distance), 
                     min(distance), 
                     max(distance),
                     stats::median(distance), 
                     stats::sd(distance))
  
  distances_df_summarised <- data.frame(distances_df_summarised)
  
  colnames(distances_df_summarised) <- c("pair", 
                                         "mean", 
                                         "min", 
                                         "max", 
                                         "median", 
                                         "std_dev")
  
  for (i in seq(nrow(distances_df_summarised))) {
    # Get cell_types for each pair
    cell_types <- strsplit(distances_df_summarised[i,"pair"], "/")[[1]]
    
    distances_df_summarised[i, "reference"] <- cell_types[1]
    distances_df_summarised[i, "target"] <- cell_types[2]
  }
  
  return(distances_df_summarised)
}

