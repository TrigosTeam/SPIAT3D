#' @title Calculate prevalence of 3D grid metrics.
#'
#' @description This functions calculates the 'prevalence', or equivalently,
#'     the proportion of rectangular prisms of the 3D grid metrics has a value
#'     above/below the required threshold.
#'
#' @param grid_metrics A data frame containing the proportion/entropy and
#'     spatial information for each rectangular prism. Obtained from the output
#'     of the calculate_cell_proportion_grid_metrics3D and
#'     calculate_entropy_grid_metrics3D functions.
#' @param metric_colname A string specifying the name of the column in
#'     `grid_metrics` containing the proportion/entropy information. Should be
#'     'proportion' or 'entropy'.
#' @param threshold A numeric between 0 and 1 representing the threshold value.
#' @param above A logical indicating whether to calculate the proportion of
#'     rectangular prisms which surpass the threshold (TRUE) or do not (FALSE).
#'     Defaults to TRUE.
#'
#' @return A numeric representing the prevalence value.
#'
#' @examples
#' # Get simulated SpatialExperiment object to use as an example for analysis
#' simulated_spe <- readRDS(system.file("extdata", "simulated_spe.rds", package = "SPIAT3D")
#'
#' # Get cell proportion grid metrics
#' cell_prop_grid_metrics <- calculate_cell_proportion_grid_metrics3D(
#'     spe = simulated_spe,
#'     n_splits = 10,
#'     reference_cell_types = c("Tumour"),
#'     target_cell_types = c("Immune"),
#'     feature_colname = "Cell.Type",
#'     plot_image = T
#' )
#'
#' # Get prevalence from cell proportion grid metrics
#' prevalence <- calculate_prevalence3D(
#'     grid_metrics = cell_prop_grid_metrics,
#'     metric_colname = "proportion",
#'     threshold = 0.3,
#'     above = TRUE
#' )
#'
#' @export

calculate_prevalence3D <- function(grid_metrics,
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
