summarise_distances_between_cell_types3D <- function(distances_df) {
  
  pair <- distance <- NULL
  
  # summarise the results
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
  
  for (i in seq(nrow(distances_df_summarised))) {
    # Get cell_types for each pair
    cell_types <- strsplit(distances_df_summarised[i,"pair"], "/")[[1]]
    
    distances_df_summarised[i, "reference"] <- cell_types[1]
    distances_df_summarised[i, "target"] <- cell_types[2]
  }
  
  return(distances_df_summarised)
}
