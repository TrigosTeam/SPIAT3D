calculate_prevalence_gradient_AUC3D <- function(prevalence_gradient_df) {
  
  return(sum(prevalence_gradient_df$prevalence) * 0.01)
}
