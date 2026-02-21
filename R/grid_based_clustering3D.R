#' @title Find cell clusters in 3D spatial data using grid based clustering 
#'     algorithm.
#'
#' @description This function finds cell clusters in a 3D SpatialExperiment 
#'     object using the grid based clustering algorithm. 
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for 
#'     the cells. Naming of spatial coordinates MUST be "Cell.X.Position", 
#'     "Cell.Y.Position", "Cell.Z.Position" for the x-coordinate, y-coordinate 
#'     and z-coordinate of each cell.
#' @param cell_types_of_interest A character vector specifying the cell types of 
#'     interest.
#' @param n_splits A positive numeric integer specifying the number splits used
#'     to divide the x-axis, y-axis and z-axis.
#' @param minimum_cells_in_cluster A positive numeric. Clusters identified which 
#'     have less than this specified value are relabelled as not a cluster.
#' @param feature_colname A string specifying the name of the column in the 
#'     `colData` slot of the SpatialExperiment object that contains the cell 
#'     type information. Defaults to "Cell.Type"
#' @param plot_image A logical indicating whether to plot 3D spatial data with 
#'     grid based clusters. Defaults to TRUE.
#'
#' @return The same 3D SpatialExperiment object used as input for spe, with an 
#'     added column in the `colData` slot to specify which grid based cluster 
#'     each cell belongs to, and added metadata containing information needed to 
#'     plot the grid based clusters.
#'
#' @examples
#' grid_based_spe <- grid_based_clustering3D(
#'     spe = SPIAT-3D::simulated_spe,
#'     cell_types_of_interest = c("Tumour", "Immune"),
#'     n_splits = 10,
#'     minimum_cells_in_cluster = 30,
#'     feature_colname = "Cell.Type",
#'     plot_image = TRUE
#' )
#' 
#' @export

grid_based_clustering3D <- function(spe,
                                    cell_types_of_interest,
                                    n_splits,
                                    minimum_cells_in_cluster,
                                    feature_colname = "Cell.Type",
                                    plot_image = TRUE) {
  
  # Other functions needed
  ### Start from the grid_prism with the maximum cell proportion.
  ## Look left, right, forward, back, up and down and see if that grid_prism has at least threshold cell proportion value
  ## If it does, add it to the answer
  ## Keep doing this until adjacent grid prisms don't have above threshold, or if you hit a boundary, or it has already been removed
  ## Return a vector containing all the grid prism numbers which COULD be part of the cluster
  calculate_grid_prism_numbers_in_cluster3D <- function(curr_grid_prism_number, 
                                                        grid_prism_cell_proportions, 
                                                        threshold_cell_proportion,
                                                        n_splits,
                                                        answer) {
    
    ## If answer already has curr_grid_prism_number, go back
    if (as.character(curr_grid_prism_number) %in% answer) return(answer)
    
    grid_prism_numbers <- names(grid_prism_cell_proportions)
    
    ## If curr_grid_prism_number has already been removed from grid_prism_numbers, go back
    if (!(as.character(curr_grid_prism_number) %in% grid_prism_numbers)) return(answer)
    
    
    if (grid_prism_cell_proportions[as.character(curr_grid_prism_number)] > threshold_cell_proportion) {
      
      answer <- c(answer, as.character(curr_grid_prism_number))
      
      ### CHECK RIGHT, LEFT, FORWARD, BACKWARD, UP, DOWN
      ## Need to check if going right, left, forward, backward, up or down is possible
      
      # Right
      if (curr_grid_prism_number%%n_splits != 0) {
        answer <- calculate_grid_prism_numbers_in_cluster3D(curr_grid_prism_number + 1,
                                                            grid_prism_cell_proportions,
                                                            threshold_cell_proportion,
                                                            n_splits,
                                                            answer)
      }
      
      # Left
      if (curr_grid_prism_number%%n_splits != 1) {
        answer <- calculate_grid_prism_numbers_in_cluster3D(curr_grid_prism_number - 1,
                                                            grid_prism_cell_proportions,
                                                            threshold_cell_proportion,
                                                            n_splits,
                                                            answer)
      }
      
      # Forward
      if ((curr_grid_prism_number - 1)%%(n_splits^2) < n_splits^2 - n_splits) {
        answer <- calculate_grid_prism_numbers_in_cluster3D(curr_grid_prism_number + n_splits,
                                                            grid_prism_cell_proportions,
                                                            threshold_cell_proportion,
                                                            n_splits,
                                                            answer)
      }
      
      # Backward
      if (curr_grid_prism_number%%(n_splits^2) > n_splits) {
        answer <- calculate_grid_prism_numbers_in_cluster3D(curr_grid_prism_number - n_splits,
                                                            grid_prism_cell_proportions,
                                                            threshold_cell_proportion,
                                                            n_splits,
                                                            answer)
      }
      
      # Up
      if (curr_grid_prism_number <= n_splits^3 - n_splits^2) {
        answer <- calculate_grid_prism_numbers_in_cluster3D(curr_grid_prism_number + n_splits^2,
                                                            grid_prism_cell_proportions,
                                                            threshold_cell_proportion,
                                                            n_splits,
                                                            answer)
      }
      
      # Down
      if (curr_grid_prism_number > n_splits^2) {
        answer <- calculate_grid_prism_numbers_in_cluster3D(curr_grid_prism_number - n_splits^2,
                                                            grid_prism_cell_proportions,
                                                            threshold_cell_proportion,
                                                            n_splits,
                                                            answer)
      }
    }
    
    return(answer)
  }
  
  # Divide each grid prism into smaller grid prisms, keeping only sub-divisions of each grid prism which satisfy the threshold
  grid_based_cluster_recursion3D <- function(df,  # Using a df is much faster than using a spe
                                             cell_types_of_interest,
                                             threshold_cell_proportion,
                                             x, y, z, l, w, h,
                                             feature_colname,
                                             answer) {
    
    # Look at cells only in the current grid prism
    df <- df[df$Cell.X.Position >= x &
               df$Cell.X.Position < (x + l) &
               df$Cell.Y.Position >= y &
               df$Cell.Y.Position < (y + w) &
               df$Cell.Z.Position >= z &
               df$Cell.Z.Position < (z + h), ]
    
    # Get cell types from spe grid prism
    cell_types <- df[[feature_colname]]
    
    # Number of cells in prism is getting too small
    if (length(cell_types) <= 2) return(data.frame())
    
    # Get total cell proportion for chosen cell_types_of_interest
    cell_proportion <- mean(cell_types %in% cell_types_of_interest)
    
    # Keep grid prism if cell proportion is above the threshold cell proportion
    if (cell_proportion >= threshold_cell_proportion) {
      return(data.frame(x, y, z, l, w, h))
    }
    
    # some cell_types_of_interest still in the grid prism, check sub-grid prisms (8 to check)
    else if (cell_proportion > 0) {
      # (0, 0, 0)
      answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                             cell_types_of_interest,
                                                             threshold_cell_proportion,
                                                             x, y, z, l/2, w/2, h/2,
                                                             feature_colname,
                                                             data.frame()))
      
      # (0.5, 0, 0)
      answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                             cell_types_of_interest,
                                                             threshold_cell_proportion,
                                                             x + l/2, y, z, l/2, w/2, h/2,
                                                             feature_colname,
                                                             data.frame()))
      
      # (0, 0.5, 0)
      answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                             cell_types_of_interest,
                                                             threshold_cell_proportion,
                                                             x, y + w/2, z, l/2, w/2, h/2,
                                                             feature_colname,
                                                             data.frame()))
      # (0.5, 0.5, 0)
      answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                             cell_types_of_interest,
                                                             threshold_cell_proportion,
                                                             x + l/2, y + w/2, z, l/2, w/2, h/2,
                                                             feature_colname,
                                                             data.frame()))
      
      # (0, 0, 0.5)
      answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                             cell_types_of_interest,
                                                             threshold_cell_proportion,
                                                             x, y, z + h/2, l/2, w/2, h/2,
                                                             feature_colname,
                                                             data.frame()))
      
      # (0.5, 0, 0.5)
      answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                             cell_types_of_interest,
                                                             threshold_cell_proportion,
                                                             x + l/2, y, z + h/2, l/2, w/2, h/2,
                                                             feature_colname,
                                                             data.frame()))
      
      # (0, 0.5, 0.5)
      answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                             cell_types_of_interest,
                                                             threshold_cell_proportion,
                                                             x, y + w/2, z + h/2, l/2, w/2, h/2,
                                                             feature_colname,
                                                             data.frame()))
      # (0.5, 0.5, 0.5)
      answer <- rbind(answer, grid_based_cluster_recursion3D(df,
                                                             cell_types_of_interest,
                                                             threshold_cell_proportion,
                                                             x + l/2, y + w/2, z + h/2, l/2, w/2, h/2,
                                                             feature_colname,
                                                             data.frame()))
      
      return(answer)
    }
    
    # cell proportion is zero
    else {
      return(data.frame())
    }
  }
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  ## Check cell types of interst are found in the spe object
  unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  if (!(is.integer(n_splits) && length(n_splits) == 1 || (is.numeric(n_splits) && length(n_splits) == 1 && n_splits > 0 && n_splits%%1 == 0))) {
    stop("`n_splits` is not a positive integer.")
  }
  if (!(is.integer(minimum_cells_in_cluster) && length(minimum_cells_in_cluster) == 1 || (is.numeric(minimum_cells_in_cluster) && length(minimum_cells_in_cluster) == 1 && minimum_cells_in_cluster > 0 && minimum_cells_in_cluster%%1 == 0))) {
    stop("`minimum_cells_in_cluster` is not a positive integer.")
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
  
  ## Calculate proportions for each grid prism
  if (length(cell_types_of_interest) == 1) {
    grid_prism_cell_proportions <- grid_prism_cell_matrix[ , cell_types_of_interest]
  }
  else {
    grid_prism_cell_proportions <- rowSums(grid_prism_cell_matrix[ , cell_types_of_interest])
  }
  grid_prism_cell_proportions <- grid_prism_cell_proportions / rowSums(grid_prism_cell_matrix[ , unique(spe[[feature_colname]])])
  n_grid_prisms <- n_splits^3
  names(grid_prism_cell_proportions) <- seq(n_grid_prisms)
  
  
  ## Create template for final result
  result <- list()
  n_clusters <- 1
  
  ## Get dimensions of the window
  spe_coords <- data.frame(spatialCoords(spe))
  
  min_x <- min(spe_coords$Cell.X.Position)
  min_y <- min(spe_coords$Cell.Y.Position)
  min_z <- min(spe_coords$Cell.Z.Position)
  
  max_x <- max(spe_coords$Cell.X.Position)
  max_y <- max(spe_coords$Cell.Y.Position)
  max_z <- max(spe_coords$Cell.Z.Position)
  
  length <- round(max_x - min_x)
  width  <- round(max_y - min_y)
  height <- round(max_z - min_z)
  
  ## Get distance of row, col and lay
  d_row <- length / n_splits
  d_col <- width / n_splits
  d_lay <- height / n_splits
  
  
  ### CLUSTER DETECTION RECURSIVE ALGORITHM LOOP ###
  
  # First, remove all 0s and NANs from grid_prism_cell_proportions
  grid_prism_cell_proportions <- grid_prism_cell_proportions[grid_prism_cell_proportions != 0 & !is.nan(grid_prism_cell_proportions)]
  
  while (length(grid_prism_cell_proportions) != 0) {
    # Get the maximum cell proportion and its corresponding grid prism number
    maximum_cell_proportion <- max(grid_prism_cell_proportions)
    maximum_cell_proportion_prism_number <- as.numeric(names(which.max(grid_prism_cell_proportions)))
    
    # Break out the loop if maximum cell proportion is less than 0.5
    if (maximum_cell_proportion < 0.5) break 
    
    # Else, find all the grid prisms adjacent to the maximum cell proportion grid prism. 
    # These are potentially apart of the cluster
    # Adjacent grid prisms must have cell proportion > 0.25 * max cell proportion
    grid_prisms_in_cluster <- calculate_grid_prism_numbers_in_cluster3D(maximum_cell_proportion_prism_number,
                                                                        grid_prism_cell_proportions,
                                                                        0.25 * maximum_cell_proportion,
                                                                        n_splits,
                                                                        c())
    
    # Perform the recursive algorithm on each grid prism potentially apart of the cluster to get a more precise shape of each cluster
    # Create data frame with spatial coords and cell types as columns. Use this as input
    result[[n_clusters]] <- data.frame()
    df <- spe_coords
    df[[feature_colname]] <- spe[[feature_colname]] 
    for (grid_prism in as.numeric(grid_prisms_in_cluster)) {
      result[[n_clusters]] <- rbind(result[[n_clusters]],
                                    grid_based_cluster_recursion3D(df,
                                                                   cell_types_of_interest,
                                                                   0.75 * maximum_cell_proportion,
                                                                   ((grid_prism - 1) %% n_splits) * d_row + round(min_x),
                                                                   (floor(((grid_prism - 1) %% n_splits^2) / n_splits)) * d_col + round(min_y),
                                                                   (floor((grid_prism - 1) / n_splits^2)) * d_lay + round(min_z),
                                                                   d_row, d_col, d_lay,
                                                                   feature_colname,
                                                                   data.frame()))
      
      
    }
    colnames(result[[n_clusters]]) <- c("x", "y", "z", "l", "w", "h")
    n_clusters <- n_clusters + 1
    
    # Remove grid prisms which have just been examined
    grid_prism_cell_proportions <- grid_prism_cell_proportions[setdiff((names(grid_prism_cell_proportions)), 
                                                                       grid_prisms_in_cluster)]
    
  }
  # Name each grid_based cluster
  names(result) <- paste("cluster", seq_len(length(result)), sep = "_")
  
  ## Add grid_based_cluster column to spe, indicating which cluster each cell belongs to
  spe$grid_based_cluster <- 0
  cluster_number <- 1
  
  for (i in seq_len(length(result))) {
    cluster_info <- result[[paste("cluster", i, sep = "_")]]
    for (j in seq(nrow(cluster_info))) {
      x <- cluster_info$x[j]
      y <- cluster_info$y[j]
      z <- cluster_info$z[j]
      l <- cluster_info$l[j]
      w <- cluster_info$w[j]
      h <- cluster_info$h[j]
      
      spe$grid_based_cluster <- ifelse(spe_coords$Cell.X.Position >= x &
                                         spe_coords$Cell.X.Position < (x + l) &
                                         spe_coords$Cell.Y.Position >= y &
                                         spe_coords$Cell.Y.Position < (y + w) &
                                         spe_coords$Cell.Z.Position >= z &
                                         spe_coords$Cell.Z.Position < (z + h) &
                                         spe[[feature_colname]] %in% cell_types_of_interest, 
                                       cluster_number, 
                                       spe$grid_based_cluster)
      
    }
    # Check if current cluster surpasses the minimum_cells_in_cluster threshold
    if (sum(spe$grid_based_cluster == cluster_number) < minimum_cells_in_cluster) {
      spe$grid_based_cluster[spe$grid_based_cluster == cluster_number] <- 0
      result[[paste("cluster", i, sep = "_")]] <- NULL
      
    }
    else {
      cluster_number <- cluster_number + 1 
    }
  }
  
  n_clusters <- max(spe$grid_based_cluster)
  if (n_clusters == 0) {
    stop("All clusters identified do not meet the `minimum_cells_in_cluster` threshold. Consider lowering the `minimum_cells_in_cluster` parameter.")
  }
  
  # re-name each grid_based cluster
  names(result) <- paste("cluster", seq_len(length(result)), sep = "_")
  
  # Add grid_clustering result to spe metadata
  spe@metadata[["grid_prisms"]] <- result
  
  ## Plot
  if (plot_image) {
    fig <- plot_grid_based_clusters3D(spe, feature_colname = feature_colname)
    methods::show(fig)
  }
  
  return(spe)
}
