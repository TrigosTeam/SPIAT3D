#' @title Function to plot alpha hull clusters in 3D spatial data.
#'
#' @description This functions plots the cells of 3D spatial data and also
#'     displays the alpha cluster clusters in a 3D SpatialExperiment Object
#'     where alpha hull clusters have already been identified.
#' 
#' @param spe A SpatialExperiment object containing 3D spatial information for 
#'     the cells. Naming of spatial coordinates MUST be "Cell.X.Position", 
#'     "Cell.Y.Position", "Cell.Z.Position" for the x-coordinate, y-coordinate 
#'     and z-coordinate of each cell. It must also contain the cell clustering 
#'     information, obtained by passing the SpatialExperiment object through the 
#'     alpha_hull_clustering3D function.
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
#' @return A Plotly object plotting the cells and alpha hull clusters of the 3D 
#'     SpatialExperiment Object.
#'
#' @examples
#' alpha_hull_spe <- alpha_hull_clustering3D(
#'     spe = SPIAT-3D::simulated_spe,
#'     cell_types_of_interest = c("Tumour", "Immune"),
#'     alpha = 8,
#'     feature_colname = "Cell.Type",
#'     plot_image = TRUE
#' )
#' 
#' fig <- plot_alpha_hull_clusters3D(
#'     spe_with_alpha_hull = alpha_hull_spe,
#'     plot_cell_types = NULL,
#'     plot_colours = NULL,
#'     feature_colname = "Cell.Type"
#' )
#' 
#' @export

plot_alpha_hull_clusters3D <- function(spe_with_alpha_hull, 
                                       plot_cell_types = NULL,
                                       plot_colours = NULL,
                                       feature_colname = "Cell.Type") {
  
  # Check input parameters
  if (class(spe_with_alpha_hull) != "SpatialExperiment") {
    stop("`spe_with_alpha_hull` is not a SpatialExperiment object.")
  }
  if (!is.character(feature_colname)) {
    stop("`feature_colname` is not a character.")
  }
  if (is.null(spe_with_alpha_hull[[feature_colname]])) {
    stop(paste("No column called", feature_colname, "found in spe object."))
  }
  
  ## If no cell types chosen, use all cell types found in data frame
  if (is.null(plot_cell_types)) plot_cell_types <- unique(spe_with_alpha_hull[[feature_colname]])
  
  ## If cell types have been chosen, check they are found in the spe object
  unknown_cell_types <- setdiff(plot_cell_types, spe_with_alpha_hull[[feature_colname]])
  if (length(unknown_cell_types) != 0) {
    stop(paste("The following plot_cell_types are not found in the spe object:\n   ",
               paste(unknown_cell_types, collapse = ", ")))
  }
  
  ## If no colours inputted, use rainbow palette
  if (is.null(plot_colours)) {
    plot_colours <- rainbow(length(plot_cell_types))
  }
  
  ## User inputs mismatching cell types and colours
  if (length(plot_cell_types) != length(plot_colours)) {
    stop("Length of plot_cell_types is not equal to length of plot_colours")
  }
  
  ## Convert spe object to data frame
  df <- data.frame(spatialCoords(spe_with_alpha_hull), "Cell.Type" = spe_with_alpha_hull[[feature_colname]])
  
  ## Factor for feature column
  df[["Cell.Type"]] <- factor(df[, "Cell.Type"],
                              levels = plot_cell_types)
  
  ## Add points to fig
  fig <- plot_ly() %>%
    add_trace(
      data = df,
      type = "scatter3d",
      mode = 'markers',
      x = ~Cell.X.Position,
      y = ~Cell.Y.Position,
      z = ~Cell.Z.Position,
      marker = list(size = 2),
      color = ~Cell.Type,
      colors = plot_colours
    ) %>% 
    layout(scene = list(xaxis = list(title = 'x'),
                        yaxis = list(title = 'y'),
                        zaxis = list(title = 'z')))
  
  
  ## Get alpha hull numbers (ignoring 0)
  alpha_hull_clusters <- spe_with_alpha_hull$alpha_hull_cluster[spe_with_alpha_hull$alpha_hull_cluster != 0]
  
  # Get number of alpha hulls
  n_alpha_hulls <- length(unique(alpha_hull_clusters))
  
  vertices <- spe_with_alpha_hull@metadata$alpha_hull$vertices
  faces <- data.frame(spe_with_alpha_hull@metadata$alpha_hull$faces)
  alpha_hull_colours <- rainbow(n_alpha_hulls)
  
  ## Add alpha hulls to fig, one by one  
  for (i in seq(n_alpha_hulls)) {
    faces_temp <- faces[faces[ , 1] %in% which(alpha_hull_clusters == i) , ]
    
    ## Ignore the weird cases where some cells represent clusters, but no faces are associated with them??
    if (nrow(faces_temp) == 0) next
    
    opacity_level <- 0.20
    
    fig <- fig %>%
      add_trace(
        type = 'mesh3d',
        x = vertices[, 1], 
        y = vertices[, 2], 
        z = vertices[, 3],
        i = faces_temp[, 1] - 1, 
        j = faces_temp[, 2] - 1, 
        k = faces_temp[, 3] - 1,
        opacity = opacity_level,
        facecolor = rep(alpha_hull_colours[i], nrow(faces_temp))
      )
  }
  
  return(fig)
}
