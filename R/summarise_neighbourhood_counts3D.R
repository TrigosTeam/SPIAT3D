#' @title Summarise neighbourhood counts data.
#'
#' @description This function summarises the neighbourhood counts data, showing 
#'     mean, min, max, median and standard deviation for each reference and
#'     target cell type pair.
#' 
#' @param neighbourhood_counts_dfs A data frame obtained from the output of the 
#'     calculate_neighbourhood_counts3D function.
#'
#' @return A data frame containing the summarised data.
#'
#' @examples
#' result <- calculate_neighbourhood_counts3D(
#'     spe = SPIAT-3D::simulated_spe,
#'     reference_cell_type = "Tumour",
#'     target_cell_types = c("Tumour", "Immune"),
#'     radius = 30,
#'     feature_colname = "Cell.Type",
#'     show_summary = TRUE,
#'     plot_image = TRUE
#' )
#' 
#' neighbourhood_counts_summary <- summarise_neighbourhood_counts3D(
#'     neighbourhood_counts_df = result
#' )
#' 
#' @export

summarise_neighbourhood_counts3D <- function(neighbourhood_counts_df) {
  
  ## Target cell types will be all the columns except the first column
  target_cell_types <- colnames(neighbourhood_counts_df)[c(-1)]
  
  ## Set up data frame for summarised_results list
  df <- data.frame(row.names = c("mean", "min", "max", "median", "st_dev"))
  
  for (target_cell_type in target_cell_types) {
    
    ## Get statistical measures for each target cell type
    target_cell_type_values <- neighbourhood_counts_df[[target_cell_type]]
    df[[target_cell_type]] <- c(mean(target_cell_type_values),
                                min(target_cell_type_values),
                                max(target_cell_type_values),
                                median(target_cell_type_values),
                                sd(target_cell_type_values))
    
  }
  
  return(data.frame(t(df)))
}
