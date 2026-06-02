#' @title Calculate mixing scores gradient on 3D spatial data.
#'
#' @description This function calculates the mixing scores gradient on a 3D
#'     SpatialExperiment Object. See paper on theory behind mixing scores
#'     gradient (I ain't explaining it here).
#'
#' @param spe A SpatialExperiment object containing 3D spatial information for
#'     the cells. Naming of spatial coordinates MUST be "Cell.X.Position",
#'     "Cell.Y.Position", "Cell.Z.Position" for the x-coordinate, y-coordinate
#'     and z-coordinate of each cell.
#' @param reference_cell_type A string specifying the reference cell type.
#' @param target_cell_type A string specifying the target cell type.
#' @param radii A positive, ascending numeric vector specifying the set of
#'     radius values used.
#' @param feature_colname A string specifying the name of the column in the
#'     `colData` slot of the SpatialExperiment object that contains the cell
#'     type information. Defaults to "Cell.Type".
#' @param plot_image A logical indicating whether to plot mixing scores gradient
#'     as a  line graph. Defaults to TRUE.
#'
#' @return A data frame containing mixing score values and associated
#'     information across each radii (rows).
#'
#' @examples
#' # Get simulated SpatialExperiment object to use as an example for analysis
#' simulated_spe <- readRDS(system.file("extdata", "simulated_spe.rds", package = "SPIAT3D"))
#'
#' result <- calculate_mixing_scores_gradient3D(
#'     spe = simulated_spe,
#'     reference_cell_type = "Tumour",
#'     target_cell_type = "Immune",
#'     radii = seq(20, 100, 10),
#'     feature_colname = "Cell.Type",
#'     plot_image = TRUE
#' )
#'
#' @export

calculate_mixing_scores_gradient3D <- function(spe,
                                               reference_cell_type,
                                               target_cell_type,
                                               radii,
                                               feature_colname = "Cell.Type",
                                               plot_image = TRUE) {

  # Check radii input
  if (!(is.numeric(radii) && length(radii) > 1)) {
    stop("`radii` is not a numeric vector with at least 2 values")
  }

  # Set up output
  result <- data.frame(matrix(nrow = length(radii), ncol = 8))
  colnames(result) <- c("ref_cell_type",
                        "tar_cell_type",
                        "n_ref_cells",
                        "n_tar_cells",
                        "n_ref_tar_interactions",
                        "n_ref_ref_interactions",
                        "mixing_score",
                        "normalised_mixing_score")

  # Iterate through each radius
  for (i in seq(length(radii))) {
    mixing_scores <- calculate_mixing_scores3D(spe,
                                               reference_cell_type,
                                               target_cell_type,
                                               radii[i],
                                               feature_colname)

    result[i, ] <- mixing_scores
  }

  # Add a radius column to the result
  result$radius <- radii

  if (plot_image) {
    fig1 <- plot_mixing_scores_gradient3D(result, "NMS")
    fig2 <- plot_mixing_scores_gradient3D(result, "MS")
    combined_fig <- plot_grid(fig1, fig2, nrow = 2)
    methods::show(combined_fig)
  }

  return(result)
}
