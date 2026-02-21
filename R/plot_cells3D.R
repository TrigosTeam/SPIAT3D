#' @title Function to plot cells in 3D spatial data.
#'
#' @description This functions plots the cells of 3D spatial data in a 3D 
#'     SpatialExperiment Object.
#' 
#' @param spe A SpatialExperiment object containing 3D spatial information for 
#'     the cells. Naming of spatial coordinates MUST be "Cell.X.Position", 
#'     "Cell.Y.Position", "Cell.Z.Position" for the x-coordinate, y-coordinate 
#'     and z-coordinate of each cell.
#' @param plot_cell_types A string vector specifying the cell types to plot. If 
#'     NULL, all cell types in the `feature_colname` column will be considered. 
#'     Defaults to NULL.
#' @param plot_colours A string vector specifying the colours of the cell types
#'     when plotting. Must match the number of cell types specified in 
#'     `plot_cell_types`. If NULL, the rainbow color pallete will be used. 
#'     Defaults to NULL.
#' @param feature_colname A string specifying the name of the column in the 
#'     `colData` slot of the SpatialExperiment object that contains the cell 
#'     type information. Defaults to "Cell.Type"
#'
#' @return A Plotly object plotting the cells of the 3D SpatialExperiment 
#'     Object.
#'
#' @examples
#' fig <- plot_cells3D(
#'     spe = SPIAT-3D::simulated_spe,
#'     plot_cell_types = NULL,
#'     plot_colours = NULL,
#'     feature_colname = "Cell.Type"
#' )
#' 
#' @export

plot_cells3D <- function(spe,
                         plot_cell_types = NULL,
                         plot_colours = NULL,
                         feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spe) != "SpatialExperiment") {
    stop("`spe` is not a SpatialExperiment object.")
  }
  if (!is.null(plot_cell_types)) {
    if(!is.character(plot_cell_types)) {
      stop("`plot_cell_types` is not a character vector.")
    }
  } 
  if (!is.null(plot_colours)) {
    non_colours <- plot_colours[which(!(sapply(plot_colours, function(X) {
      tryCatch(is.matrix(col2rgb(X)), 
               error = function(e) FALSE)
    })))]
    if (length(non_colours) > 0) {
      stop(paste("The following plot_colours are not colours:\n   ",
                 paste(non_colours, collapse = ", ")))
    } 
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe[[feature_colname]])) {
    stop(paste(feature_colname, "is not a valid column in your spe object."))
  }
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe), "Cell.Type" = spe[[feature_colname]])
  
  ## If no cell types chosen, use all cell types found in data frame
  if (is.null(plot_cell_types)) {
    warning("plot_cell_types not specified, all cell types found in the spe object will be used.")
    plot_cell_types <- unique(df[["Cell.Type"]])
  }
  ## If no colours inputted, use rainbow palette
  if (is.null(plot_colours)) {
    warning("plot_colours not specified, rainbow palette will be used.")
    plot_colours <- rainbow(length(plot_cell_types))
  }
  ## User inputs mismatching cell types and colours
  if (length(plot_cell_types) != length(plot_colours)) {
    stop("Length of plot_cell_types is not equal to length of plot_colours")
  }
  
  ## If cell types have been chosen, check they are found in the spe object
  spe_cell_types <- unique(spe[[feature_colname]])
  unknown_cell_types <- setdiff(plot_cell_types, spe_cell_types)
  
  if (length(unknown_cell_types) == length(plot_cell_types)) {
    stop("None of the plot_cell_types are found in the spe object")
  }
  
  if (length(unknown_cell_types) != 0) {
    warning(paste("The following plot_cell_types are not found in the spe object:\n   ",
                  paste(unknown_cell_types, collapse = ", ")))
    plot_colours <- plot_colours[which(plot_cell_types %in% spe_cell_types)]
    plot_cell_types <- intersect(plot_cell_types, spe_cell_types)
  }
  
  ## Factor for feature column
  df[, "Cell.Type"] <- factor(df[, "Cell.Type"],
                              levels = plot_cell_types)
  
  ## Plot
  fig <- plot_ly(df,
                 type = "scatter3d",
                 mode = 'markers',
                 x = ~Cell.X.Position,
                 y = ~Cell.Y.Position,
                 z = ~Cell.Z.Position,
                 color = ~Cell.Type,
                 colors = plot_colours,
                 marker = list(size = 2))
  
  fig <- fig %>% layout(scene = list(xaxis = list(title = 'x', showgrid = T, showaxeslabels = F, showticklabels = T, gridwidth = 5, 
                                                  titlefont = list(size = 20), tickfont = list(size = 15)),
                                     yaxis = list(title = 'y', showgrid = T, showaxeslabels = F, showticklabels = T, gridwidth = 5,
                                                  titlefont = list(size = 20), tickfont = list(size = 15)),
                                     zaxis = list(title = 'z', showgrid = T, showaxeslabels = F, showticklabels = T, gridwidth = 5,
                                                  titlefont = list(size = 20), tickfont = list(size = 15))))
  
  return(fig)
}
