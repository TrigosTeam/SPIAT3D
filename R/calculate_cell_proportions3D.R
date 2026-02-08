#' @title Calculate cell proportions in 3D spatial data.
#'
#' @description This function calculates the proportions of different cell types in a 3D SpatialExperiment object. 
#'    It can optionally plot a bar chart of the cell proportions.
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for the cells.
#' @param cell_types_of_interest A character vector specifying the cell types of interest.
#'    If NULL, all cell types in the `feature_colname` column will be considered.
#' @param feature_colname A string specifying the name of the column in the `colData` slot of the SpatialExperiment
#'    object that contains the cell type information.
#' @param plot_image A logical indicating whether to plot violin plots of the minimum distances 
#'    between cell type pairs. Defaults to TRUE.
#'
#' @return A data frame containing the cell types, their frequencies, proportions, and percentages.
#'
#' @examples
#' cell_proportions <- calculate_cell_proportions3D(
#'     spe = SPIAT-3D::simulated_spe, 
#'     cell_types_of_interest = NULL, 
#'     feature_colname = "Cell.Type", 
#'     plot_image = TRUE
#' )
#' 
#' @export


calculate_cell_proportions3D <- function(spe,
                                         cell_types_of_interest = NULL, 
                                         feature_colname = "Cell.Type",
                                         plot_image = TRUE) {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  # Check if there are empty strings or string of only spaces in 'cell_types_of_interest'
  if (length(spe[[feature_colname]][trimws(spe[[feature_colname]]) == ""]) > 0) {
    stop("spe cannot contain cell types that are an empty string or a string of only spaces.")
  }
  if (ncol(spe) == 0) {
    stop("No cells found for calculating cell proportions.")
  }
  if (!(is.null(cell_types_of_interest) || is.character(cell_types_of_interest))) {
    stop("`cell_types_of_interest` is not a character vector or NULL.")
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
  
  # Creates frequency/bar plot of all cell types in the entire image
  cell_proportions <- data.frame(table(spe[[feature_colname]]))
  names(cell_proportions) <- c("cell_type", 'frequency')
  
  # Only include cell types the user has chosen
  if (!is.null(cell_types_of_interest)) {
    
    ## If cell types have been chosen, check they are found in the spe object
    unknown_cell_types <- setdiff(cell_types_of_interest, cell_proportions$cell_type)
    if (length(unknown_cell_types) != 0) {
      stop(paste("The following cell types in cell_types_of_interest are not found in the spe object:\n   ",
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
