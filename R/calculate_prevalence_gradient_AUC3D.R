calculate_prevalence_gradient_AUC3D <- function(prevalence_gradient_df) {
  
  return(
    sum(diff(prevalence_gradient_df$threshold) * (head(prevalence_gradient_df$prevalence, -1) + tail(prevalence_gradient_df$prevalence, -1)) / 2)
  )
  
}
