# 6 categories: ME, RE, SE, MN, RN, SN
# 10 simulations for each category: 60 simulations total.
# First 5 simulations, increasing arrangement parameter (e.g. increasing mixing)
# Second 5 simulations, increasing size parameter (e.g. increasing ellipsoid size)
library(tibble)
library(ggplot2)
library(gridExtra)

generate_simulation_parameters_for_metric_testing <- function() {
  
  simulation_parameters <- data.frame(arrangement = rep(rep(c("mixed", "ringed", "separated"), each = 10), 2),
                                      shape = rep(c("ellipsoid", "network"), each = 30))
  
  parameter_values <- list("mixed_cell_type_A_proportion" = list("fixed" = 0.5, "increasing" = seq(0.3, 0.7, 0.1), "type" = c("arrangement", "mixed")),
                           "ringed_ring_width_factor" = list("fixed" = 0.15, "increasing" = seq(0.1, 0.2, 0.025), "type" = c("arrangement", "ringed")),
                           "separated_cluster1_x_coordinate" = list("fixed" = 150, "increasing" = seq(125, 175, 12.5), "type" = c("arrangement", "separated")),
                           "ellipsoid_x_radius" = list("fixed" = 75, "increasing" = seq(50, 100, 12.5), "type" = c("shape", "ellipsoid")),
                           "ellipsoid_y_radius" = list("fixed" = 100, "increasing" = seq(75, 125, 12.5), "type" = c("shape", "ellipsoid")),
                           "ellipsoid_z_radius" = list("fixed" = 125, "increasing" = seq(100, 150, 12.5), "type" = c("shape", "ellipsoid")),
                           "network_width" = list("fixed" = 30, "increasing" = seq(25, 35, 2.5), "type" = c("shape", "network")))

  for (parameter in names(parameter_values)) {
    simulation_parameters[[parameter]] <- ifelse(simulation_parameters[[parameter_values[[parameter]]$type[1]]] == parameter_values[[parameter]]$type[2],
                                                                 parameter_values[[parameter]][["fixed"]], NA)
  }
  
  variable_order <- c("mixed_cell_type_A_proportion", "ellipsoid_x_radius",
                      "ringed_ring_width_factor", "ellipsoid_x_radius",
                      "separated_cluster1_x_coordinate", "ellipsoid_x_radius",
                      "mixed_cell_type_A_proportion", "network_width",
                      "ringed_ring_width_factor", "network_width",
                      "separated_cluster1_x_coordinate", "network_width")
  
  simulation_parameters$variable <- rep(variable_order, each = 5)
  
  prev_variable <- ""
  for (i in seq(nrow(simulation_parameters))) {
    curr_variable <- simulation_parameters$variable[i]
    
    if (curr_variable == prev_variable) {
      next
    }
    
    simulation_parameters[[curr_variable]][i:(i + 4)] <- parameter_values[[curr_variable]][["increasing"]]
    if (curr_variable == "ellipsoid_x_radius") {
      simulation_parameters[["ellipsoid_y_radius"]][i:(i + 4)] <- parameter_values[["ellipsoid_y_radius"]][["increasing"]]
      simulation_parameters[["ellipsoid_z_radius"]][i:(i + 4)] <- parameter_values[["ellipsoid_z_radius"]][["increasing"]]
    }
    
    prev_variable <- curr_variable
  }
  
  # Get ring width as well, depending on if it is an ellipsoid or network shape
  ringed_ring_width_factors <- simulation_parameters$ringed_ring_width_factor
  ellipsoid_x_radii <- simulation_parameters[["ellipsoid_x_radius"]]
  ellipsoid_y_radii <- simulation_parameters[["ellipsoid_y_radius"]]
  ellipsoid_z_radii <- simulation_parameters[["ellipsoid_z_radius"]]
  network_widths <- simulation_parameters[["network_width"]]

  ringed_ring_widths_ellipsoids <- ringed_ring_width_factors * 
    (ellipsoid_x_radii + ellipsoid_y_radii + ellipsoid_z_radii) / 3
  
  ringed_ring_widths_networks <- ringed_ring_width_factors * network_widths
  
  ringed_ring_widths <- ifelse(!is.na(ringed_ring_widths_ellipsoids), ringed_ring_widths_ellipsoids, ringed_ring_widths_networks)
  
  simulation_parameters <- add_column(simulation_parameters, ringed_ring_width = ringed_ring_widths, .after = 4)
  
  return(simulation_parameters)
}

generate_simulation_metadata_for_metric_testing <- function() {
  
  # Define fixed background metadata parameters
  background_parameters <- list(
    number_of_cells = 30000,
    length = 600, # Units: micrometers (um)
    width = 600,
    height = 300,
    minimum_distance_between_cells = 10,
    cell_types = c("A", "B", "O"),
    cell_proportions = c(0, 0, 1)
  )
  
  spe_metadata_background <- spe_metadata_background_template("random", NULL)
  spe_metadata_background$background$length <- background_parameters$length
  spe_metadata_background$background$width <-  background_parameters$width
  spe_metadata_background$background$height <- background_parameters$height
  spe_metadata_background$background$n_cells <- background_parameters$number_of_cells
  spe_metadata_background$background$minimum_distance_between_cells <- background_parameters$minimum_distance_between_cells
  spe_metadata_background$background$cell_types <- background_parameters$cell_types
  spe_metadata_background$background$cell_proportions <- background_parameters$cell_proportions
  
  simulation_parameters <- generate_simulation_parameters_for_metric_testing()
  
  simulation_metadata <- list()
  
  for (index in seq_len(nrow(simulation_parameters))) {
    
    shape <- simulation_parameters$shape[index]
    arrangement <- simulation_parameters$arrangement[index]
    
    # Determine cluster_type from arrangement
    if (arrangement == "mixed") {
      cluster_type <- "regular"
    }
    else if (arrangement == "ringed") {
      cluster_type <- "ring"
    }
    else if (arrangement == "separated") {
      cluster_type <- "regular"
    }
    
    # Add cluster metadata
    spe_metadata_cluster <- spe_metadata_cluster_template(cluster_type, shape, spe_metadata_background)
    spe_metadata_cluster$cluster_1$cluster_cell_types <- c('A', 'B')
    spe_metadata_cluster$cluster_1$centre_loc <- c(300, 300, 150)
    
    # Add sphere if arrangement is separated
    if (arrangement == "separated") {
      spe_metadata_cluster <- spe_metadata_cluster_template(cluster_type, "sphere", spe_metadata_cluster)
      spe_metadata_cluster$cluster_2$cluster_cell_types <- c('A', 'B')
      spe_metadata_cluster$cluster_2$cluster_cell_proportions <- c(0, 1)
      spe_metadata_cluster$cluster_2$radius <- 100
      spe_metadata_cluster$cluster_2$centre_loc <- c(450, 300, 150)
    }
    
    simulation_metadata[[index]] <- spe_metadata_cluster
    
    if (shape == "ellipsoid") {
      simulation_metadata[[index]][["cluster_1"]][["x_radius"]] <- simulation_parameters$ellipsoid_x_radius[index]
      simulation_metadata[[index]][["cluster_1"]][["y_radius"]] <- simulation_parameters$ellipsoid_y_radius[index]
      simulation_metadata[[index]][["cluster_1"]][["z_radius"]] <- simulation_parameters$ellipsoid_z_radius[index]
      simulation_metadata[[index]][["cluster_1"]][["y_z_rotation"]] <- runif(1, min = 0, max = 180) # Choose random angle
      simulation_metadata[[index]][["cluster_1"]][["x_z_rotation"]] <- runif(1, min = 0, max = 180) # Choose random angle
      simulation_metadata[[index]][["cluster_1"]][["x_y_rotation"]] <- runif(1, min = 0, max = 180) # Choose random angle
    }
    else if (shape == "network") {
      simulation_metadata[[index]][["cluster_1"]][["n_edges"]] <- 20
      simulation_metadata[[index]][["cluster_1"]][["width"]] <- simulation_parameters$network_width[index]
      simulation_metadata[[index]][["cluster_1"]][["radius"]] <- 125
    }
    
    
    ### Arrangements
    if (arrangement == "mixed") {
      simulation_metadata[[index]][["cluster_1"]][["cluster_cell_proportions"]] <-
        c(simulation_parameters$mixed_cell_type_A_proportion[index], 1 - simulation_parameters$mixed_cell_type_A_proportion[index])
      
    }
    else if (arrangement == "ringed") {
      simulation_metadata[[index]][["cluster_1"]][["cluster_cell_proportions"]] <- c(1, 0)
      simulation_metadata[[index]][["cluster_1"]][["ring_cell_types"]] <- c('A', 'B')
      simulation_metadata[[index]][["cluster_1"]][["ring_cell_proportions"]] <- c(0, 1)
      simulation_metadata[[index]][["cluster_1"]][["ring_width"]] <- simulation_parameters$ringed_ring_width[index]
    }
    else if (arrangement == "separated") {
      simulation_metadata[[index]][["cluster_1"]][["cluster_cell_proportions"]] <- c(1, 0)
      simulation_metadata[[index]][["cluster_1"]][["centre_loc"]] <-
        c(simulation_parameters$separated_cluster1_x_coordinate[index], 300, 150)
    }
  }

  return(simulation_metadata)
}

# In this function, change the metric function used for the metric you want to test
generate_and_analyse_simulations_for_metric_testing <- function() {
  
  simulation_metadata <- generate_simulation_metadata_for_metric_testing()
  
  results <- list()
  
  for (i in seq_along(simulation_metadata)) {
    print(i)
    curr_simulation <- simulate_spe_metadata3D(simulation_metadata[[i]], plot_image = FALSE)
    
    # Choose your metric here
    metric_result <- calculate_Gcross_gradient3D(curr_simulation, 
                                                 reference_cell_type = "A",
                                                 target_cell_type = "B",
                                                 radii = seq(10, 100, 10),
                                                 plot_image = FALSE)
    
    results[[i]] <- metric_result
  }
  return(results)
}


## Code it yourself
parameters <- generate_simulation_parameters_for_metric_testing()
simulation_metadata <- generate_simulation_metadata_for_metric_testing()
results <- generate_and_analyse_simulations_for_metric_testing()

# AUC
calculate_auc <- function(x, y) {
  # Check if x and y are of equal length
  if (length(x) != length(y)) {
    stop("x and y must have the same length!")
  }
  
  # Calculate the AUC using the trapezoidal rule
  auc <- sum(diff(x) * (y[-1] + y[-length(y)]) / 2)
  
  return(auc)
}

index <- 1
plots <- list()
while (index <= nrow(parameters)) {
  curr_variable <- parameters$variable[index]
  curr_simulation_type <- paste(parameters$arrangement[index], parameters$shape[index], sep = "-")
  
  parameter_values <- parameters[[curr_variable]][index:(index + 4)]
  metric_values <- c()
  for (curr_index in index:(index + 4)) {
    auc <- calculate_auc(results[[curr_index]][["radius"]], results[[curr_index]][["observed_Gcross"]])   
    metric_values <- c(metric_values, auc)
  }

  plot_df <- data.frame(parameter = parameter_values, metric = metric_values)
  
  plot <- ggplot(plot_df, aes(parameter, metric)) + 
    geom_point() + labs(
      title = curr_simulation_type,
      x = curr_variable,
      y = "Gcross"
    )
  
  plots[[length(plots) + 1]] <- plot
  
  index <- index + 5
}

grid.arrange(grobs = plots, nrow = 3, ncol = 4)


# As is (not AUC)
plots <- list()
for (index in 51:55) {

  curr_result <- results[[index]]
  
  plot <- plot_Gcross_gradient3D(curr_result)
  
  plots[[length(plots) + 1]] <- plot
}
grid.arrange(grobs = plots, nrow = 5, ncol = 1)
