#' @title Summarise neighbourhood metric data.
#'
#' @description This function summarises the data from either neighbourhood
#'     counts, cells in neighbourhood or neighbourhood entropy, showing, min,
#'     max, median and standard deviation for each reference and target cell
#'     type pair.
#'
#' @param neighbourhood_metric_df A data frame obtained from the output of
#'     either the calculate_neighbourhood_counts3D,
#'     calculate_cells_in_neighbourhood3D or calculate_neighbourhood_entropy3D
#'     function.
#'
#' @return A data frame containing the summarised data.
#'
#' @examples
#' # Get simulated SpatialExperiment object to use as an example for analysis
#' simulated_spe <- readRDS(system.file("extdata", "simulated_spe.rds", package = "SPIAT3D"))
#'
#' # Calculate neighbourhood counts from simulated spe
#' result <- summarise_neighbourhood_metric3D(
#'     spe = simulated_spe,
#'     reference_cell_type = "Tumour",
#'     target_cell_types = c("Tumour", "Immune"),
#'     radius = 30,
#'     feature_colname = "Cell.Type",
#'     show_summary = FALSE,
#'     plot_image = FALSE
#' )
#'
#' # Get summary
#' neighbourhood_counts_summary <- summarise_neighbourhood_metric3D(
#'     neighbourhood_metric_df = result
#' )
#'
#' @export

summarise_neighbourhood_metric3D <- function(neighbourhood_metric_df) {

  ## Target cell types will be all the columns except the first column
  target_cell_types <- colnames(neighbourhood_metric_df)[c(-1)]

  ## Set up data frame for summarised_results list
  df <- data.frame(row.names = c("mean", "min", "max", "median", "st_dev"))

  for (target_cell_type in target_cell_types) {

    ## Get statistical measures for each target cell type
    target_cell_type_values <- neighbourhood_metric_df[[target_cell_type]]
    df[[target_cell_type]] <- c(mean(target_cell_type_values),
                                min(target_cell_type_values),
                                max(target_cell_type_values),
                                median(target_cell_type_values),
                                sd(target_cell_type_values))

  }

  return(data.frame(t(df)))
}
