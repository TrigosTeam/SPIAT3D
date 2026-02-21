#' @title Summarise distances between cell types.
#'
#' @description This function summarises the distances between cell types data,
#'     showing mean, min, max, median and standard deviation for each cell type
#'     pair.
#' 
#' @param distances_df A data frame obtained from the output of the 
#'     calculate_pairwise_distances_between_cell_types3D or
#'     calculate_minimum_distances_between_cell_types3D functions.
#'
#' @return A data frame containing the summarised data.
#'
#' @examples
#' minimum_distances <- calculate_minimum_distances_between_cell_types3D(
#'     spe = SPIAT-3D::simulated_spe,
#'     cell_types_of_interest = NULL,
#'     feature_colname = "Cell.Type",
#'     show_summary = FALSE,
#'     plot_image = FALSE
#' )
#' 
#' distances_summary <- summarise_distances_between_cell_types3D(
#'     distances_df = minimum_distances
#' )
#' 
#' @export

summarise_distances_between_cell_types3D <- function(distances_df) {
  
  # Just in case?
  pair <- distance <- NULL
  
  # Summarise the results
  distances_df_summarised <- distances_df %>% 
    dplyr::group_by(pair) %>%
    dplyr::summarise(mean(distance), 
                     min(distance), 
                     max(distance),
                     stats::median(distance), 
                     stats::sd(distance))
  
  distances_df_summarised <- data.frame(distances_df_summarised)
  
  colnames(distances_df_summarised) <- c("pair", 
                                         "mean", 
                                         "min", 
                                         "max", 
                                         "median", 
                                         "std_dev")
  
  # Add columns for reference cell type and target cell type.
  for (i in seq(nrow(distances_df_summarised))) {
    # Get cell_types for each pair
    cell_types <- strsplit(distances_df_summarised[i,"pair"], "/")[[1]]
    
    distances_df_summarised[i, "reference"] <- cell_types[1]
    distances_df_summarised[i, "target"] <- cell_types[2]
  }
  
  return(distances_df_summarised)
}
