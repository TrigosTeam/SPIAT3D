#' @title Calculate cross K-function on 3D spatial data.
#'
#' @description This function calculates the cross K-function on a 3D
#'     SpatialExperiment Object. See paper on theory behind cross K-function (I
#'     ain't explaining it here).
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for
#'     the cells. Naming of spatial coordinates MUST be "Cell.X.Position",
#'     "Cell.Y.Position", "Cell.Z.Position" for the x-coordinate, y-coordinate
#'     and z-coordinate of each cell.
#' @param reference_cell_type A string specifying the reference cell type.
#' @param target_cell_types A character vector specifying the target cell types.
#' @param radius A positive numeric specifying the radius value.
#' @param feature_colname A string specifying the name of the column in the
#'     `colData` slot of the SpatialExperiment object that contains the cell
#'     type information. Defaults to "Cell.Type".
#'
#' @return A data frame containing observed cross K-function values for each
#'     target cell type and expected cross K-function values.
#'
#' @examples
#' # Get simulated SpatialExperiment object to use as an example for analysis
#' simulated_spe <- readRDS(system.file("extdata", "simulated_spe.rds", package = "SPIAT3D"))
#'
#' result <- calculate_cross_K3D(
#'     spe = simulated_spe,
#'     reference_cell_type = "Tumour",
#'     target_cell_types = c("Tumour", "Immune"),
#'     radius = 30,
#'     feature_colname = "Cell.Type"
#' )
#'
#' @export

calculate_cross_K3D <- function(spe,
                                reference_cell_type,
                                target_cell_types,
                                radius,
                                feature_colname = "Cell.Type") {

  # Check if inputs are valid
  if (is.null(spe[[feature_colname]])) stop(paste("No column called", feature_colname, "found in spe object"))

  if (is.null(spe[["Cell.ID"]])) {
    warning("Temporarily adding Cell.ID column to your spe")
    spe$Cell.ID <- paste("Cell", seq(ncol(spe)), sep = "_")
  }


  ## Get expected cross K-function
  expected_cross_K <- (4/3) * pi * radius^3

  ## For reference_cell_type, check it is found in the spe object
  if (!(reference_cell_type %in% spe[[feature_colname]])) {
    warning(paste("The reference_cell_type", reference_cell_type,"is not found in the spe object"))
    result <- data.frame(observed_cross_K = NA,
                         expected_cross_K = expected_cross_K,
                         cross_K_ratio = NA)
    return(result)
  }

  ## Get rough dimensions of the window the points are in
  spe_coords <- data.frame(SpatialExperiment::spatialCoords(spe))

  length <- round(max(spe_coords$Cell.X.Position) - min(spe_coords$Cell.X.Position))
  width  <- round(max(spe_coords$Cell.Y.Position) - min(spe_coords$Cell.Y.Position))
  height <- round(max(spe_coords$Cell.Z.Position) - min(spe_coords$Cell.Z.Position))

  ## Get volume of the window the cells are in
  volume <- length * width * height

  # Number of reference cell types is constant
  n_ref_cells <- sum(spe[[feature_colname]] == reference_cell_type)

  # Define result data frame
  result <- data.frame(reference = reference_cell_type, expected = expected_cross_K)

  neighbourhood_counts_df <- calculate_neighbourhood_counts3D(spe,
                                                              reference_cell_type,
                                                              target_cell_types,
                                                              radius,
                                                              feature_colname,
                                                              show_summary = FALSE,
                                                              plot_image = FALSE)

  # Calculate cross K-fucnction for each target cell type
  for (target_cell_type in target_cell_types) {

    n_ref_tar_interactions <- sum(neighbourhood_counts_df[[target_cell_type]])

    n_tar_cells <- sum(spe[[feature_colname]] == target_cell_type)

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
