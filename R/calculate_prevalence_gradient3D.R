#' @title Calculate prevalence gradient of 3D grid metrics.
#'
#' @description This functions calculates the 'prevalence gradient', or
#'     equivalently, the proportion of rectangular prisms of the 3D grid metrics
#'     that have a value above the required threshold, for a gradient of
#'     thresholds.
#'
#' @param grid_metrics A data frame containing the proportion/entropy and
#'     spatial information for each rectangular prism. Obtained from the output
#'     of the calculate_cell_proportion_grid_metrics3D and
#'     calculate_entropy_grid_metrics3D functions.
#' @param metric_colname A string specifying the name of the column in
#'     `grid_metrics` containing the proportion/entropy information. Should be
#'     'proportion' or 'entropy'.
#' @param show_AUC A logical indicating whether to print out the prevalence
#'     gradient AUC value. Defaults to TRUE.
#' @param plot_image A logical indicating whether to plot a line graph showing
#'     prevalence vs threshold values. Defaults to TRUE.
#'
#' @return A data frame containing the prevalence values for each threshold.
#'
#' @examples
#' # Get simulated SpatialExperiment object to use as an example for analysis
#' simulated_spe <- readRDS(system.file("extdata", "simulated_spe.rds", package = "SPIAT3D"))
#'
#' # Get cell proportion grid metrics
#' cell_prop_grid_metrics <- calculate_cell_proportion_grid_metrics3D(
#'     spe = simulated_spe,
#'     n_splits = 10,
#'     reference_cell_types = c("Tumour"),
#'     target_cell_types = c("Immune"),
#'     feature_colname = "Cell.Type",
#'     plot_image = TRUE
#' )
#'
#' # Get prevalence gradient from cell proportion grid metrics
#' prevalence_gradient <- calculate_prevalence_gradient3D(
#'     grid_metrics = cell_prop_grid_metrics,
#'     metric_colname = "proportion",
#'     show_AUC = TRUE,
#'     plot_image = TRUE
#' )
#'
#' @export

calculate_prevalence_gradient3D <- function(grid_metrics,
                                            metric_colname,
                                            show_AUC = TRUE,
                                            plot_image = TRUE) {

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
    calculate_prevalence3D(grid_metrics, metric_colname, threshold)
  })

  # Show AUC of prevalence gradient graph
  if (show_AUC) {
    AUC_value <- sum(diff(result$threshold) * (head(result$prevalence, -1) + tail(result$prevalence, -1)) / 2)
    print(paste("AUC:", round(AUC_value, 2)))
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
