calculate_pairwise_distances_between_cell_types3D <- function(spe,
                                                              cell_types_of_interest = NULL,
                                                              feature_colname = "Cell.Type",
                                                              show_summary = TRUE,
                                                              plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (ncol(spe) < 2) {
    stop("There must be at least two cells in spe.")
  }
  if (!(is.null(cell_types_of_interest) && is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  if (!is.logical(show_summary)) {
    stop("`show_summary` is not a logical (TRUE or FALSE).")
  }
  if (!is.logical(plot_image)) {
    stop("`plot_image` is not a logical (TRUE or FALSE).")
  }
  
  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }
  
  # De-factor feature column in spe object
  spe[[feature_colname]] <- as.character(spe[[feature_colname]])
  
  # Subset spe to only contain the cells of interest
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spe object
    unknown_cell_types <- setdiff(cell_types_of_interest, spe[[feature_colname]])
    if (length(unknown_cell_types) != 0) {
      warning(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
                    paste(unknown_cell_types, collapse = ", ")))
    }
    
    spe <- spe[ , spe[[feature_colname]] %in% cell_types_of_interest]
  }
  # If cell_types_of_interest is NULL, use all cells in spe
  else {
    cell_types_of_interest <- unique(spe[[feature_colname]])
  }
  
  # Create a list containing the cell IDs of each cell type
  cell_type_ids <- list()
  for (cell_type in cell_types_of_interest) {
    cell_type_ids[[cell_type]] <- as.character(spe$Cell.ID[spe[[feature_colname]] == cell_type])
  }
  
  # Calculate cell to cell distances
  distance_matrix <- -1 * apcluster::negDistMat(spatialCoords(spe))
  rownames(distance_matrix) <- spe$Cell.ID
  colnames(distance_matrix) <- spe$Cell.ID
  
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
    print(summarise_distances_between_cell_types3D(result))  
  }
  
  # Plot
  if (plot_image) {
    fig <- plot_distances_between_cell_types_violin3D(result)
    methods::show(fig)
  }
  
  return(result)
}
