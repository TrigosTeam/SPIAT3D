spe1 <- spe_cluster


### 1. Basic Metrics ----------------------------------------------------------

# Calculate Cell Proportions
cell_props1 <- calculate_cell_proportions3D(spe1,
                                            cell_types_of_interest = NULL,
                                            plot_image = TRUE)
print(cell_props1)

cell_props2 <- calculate_cell_proportions3D(spe1,
                                            cell_types_of_interest = c("Tumour", "Immune"),
                                            plot_image = TRUE)
print(cell_props2)


### 2. Colocalization metrics -------------------------------------------------

### Calculate Pairwise Distances between Cells
pairwise_distances <- calculate_pairwise_distances_between_cell_types3D(spe1,
                                                                        cell_types_of_interest = c("Tumour", "Immune"),
                                                                        plot_image = TRUE)


### Calculate Minimum Distances between cells
minimum_distances <- calculate_minimum_distances_between_cell_types3D(spe1,
                                                                      cell_types_of_interest = c("Tumour", "Immune"),
                                                                      plot_image = TRUE)





### Calculate Mixing Scores
mixing_scores <- calculate_mixing_scores3D(spe1,
                                           reference_cell_types = c("Tumour", "Immune"),
                                           target_cell_types = c("Tumour", "Immune"),
                                           radius = 20)
print(mixing_scores)

mixing_scores_gradient <- calculate_mixing_scores_gradient3D(spe1,
                                                             reference_cell_type = "Tumour",
                                                             target_cell_type = "Immune",
                                                             radii = seq(50))


### Calculate cells in the neighbourhood
neighbourhood_cells <- calculate_cells_in_neighbourhood3D(spe1,
                                                          reference_cell_type = "Tumour",
                                                          target_cell_types = c("Tumour", "Immune"),
                                                          radius = 30,
                                                          plot_image = F)

neighbourhood_cells_gradient <- calculate_cells_in_neighbourhood_gradient3D(spe1,
                                                                            reference_cell_type = "Tumour",
                                                                            target_cell_types = c("Tumour", "Immune"),
                                                                            radii = seq(1, 30, 2),
                                                                            plot_image = T)


## Calculate cell proportions in the neighbourhood
neighbourhood_cell_proportions <- calculate_cells_in_neighbourhood_proportions3D(spe1,
                                                                                 reference_cell_type = "Tumour",
                                                                                 target_cell_types = c("Tumour", "Immune"),
                                                                                 radius = 20)
print(neighbourhood_cell_proportions)

neighbourhood_cell_proportions_gradient <- calculate_cells_in_neighbourhood_proportions_gradient3D(spe1,
                                                                                                   reference_cell_type = "Tumour",
                                                                                                   target_cell_types = c("Tumour", "Immune"),
                                                                                                   radii = seq(1, 50, 3))


### Calculate cross-K function
cross_K <- calculate_cross_K3D(spe1,
                               reference_cell_type = "Tumour",
                               target_cell_type = "Immune",
                               radius = 20)
print(cross_K)


cross_K_gradient <- calculate_cross_K_gradient3D(spe1,
                                                 reference_cell_type = "Tumour",
                                                 target_cell_type = "Immune",
                                                 radii = seq(1, 50))


### Calculate entropy
entropy_background <- calculate_entropy_background3D(spe1,
                                                     cell_types_of_interest = c("Tumour", "Immune"))

print(entropy_background)

entropy_result <- calculate_entropy3D(spe1,
                                      radius = 20,
                                      reference_cell_type = "Tumour",
                                      target_cell_types = c("Tumour", "Immune"))





entropy_gradient <- calculate_entropy_gradient3D(spe1,
                                                 reference_cell_type = "Tumour",
                                                 target_cell_types = c("Tumour", "Immune"),
                                                 radii = seq(1, 50, 2),
                                                 plot_image = TRUE)


### Using all_single_radius and all_gradient functions

all_single_radius_result <- calculate_all_single_radius_cc_metrics3D(spe1, "Tumour", c("Tumour", "Immune"), 20)

all_gradient_result <- calculate_all_gradient_cc_metrics3D(spe1, "Tumour", c("Tumour", "Immune"), seq(1, 50, 2))


### 3. Spatial Heterogeneity metrics ------------------------------------------

### calculate entropy grid metrics
entropy_grid_metrics <- calculate_entropy_grid_metrics3D(spe1,
                                                         n_splits = 8,
                                                         cell_types_of_interest = c("Tumour", "Immune", "Immune1"),
                                                         plot_image = TRUE)
plot_grid_metrics_discrete3D(entropy_grid_metrics, "entropy")


### calculate entropy prevalence
entropy_prevalence <- calculate_prevalence3D(entropy_grid_metrics,
                                             metric_colname = "entropy",
                                             threshold = 0.5)
print(entropy_prevalence)

entropy_prevalence_gradient <- calculate_prevalence_gradient3D(entropy_grid_metrics,
                                                               "entropy")

### calculate spatial autocorrelation
entropy_spatial_autocorrelation <- calculate_spatial_autocorrelation3D(entropy_grid_metrics,
                                                                       metric_colname = "entropy",
                                                                       weight_method = "IDW")
print(entropy_spatial_autocorrelation)


### calculate cell proportion grid metrics
cell_proportion_grid_metrics <- calculate_cell_proportion_grid_metrics3D(spe1,
                                                                         n_splits = 10,
                                                                         reference_cell_types = c("Tumour"),
                                                                         target_cell_types = c("Immune"),
                                                                         plot_image = TRUE)
plot_grid_metrics_discrete3D(cell_proportion_grid_metrics, "proportion")


### calculate cell proportion prevalence
cell_proportion_prevalence <- calculate_prevalence3D(cell_proportion_grid_metrics,
                                                     metric_colname = "proportion",
                                                     threshold = 0.5)
print(cell_proportion_prevalence)

cell_proportion_prevalence_gradient <- calculate_prevalence_gradient3D(cell_proportion_grid_metrics,
                                                                       metric_colname = "proportion")

## calculate spatial autocorrelation for cell proportions
cell_proportion_spatial_autocorrelation <- calculate_spatial_autocorrelation3D(cell_proportion_grid_metrics,
                                                                               metric_colname = "proportion",
                                                                               weight_method = 0.10)
print(cell_proportion_spatial_autocorrelation)



### 4. Clustering algorithms --------------------------------------------------
spe_alpha_hull <- alpha_hull_clustering3D(spe1, c("Tumour", "Immune"), alpha = 4.2, minimum_cells_in_alpha_hull = 15)

plot_alpha_hull_clusters3D(spe_alpha_hull, c("Tumour", "Immune", "Others"), c("orange", "skyblue", "lightgray"))

alpha_hull_props <- calculate_cell_proportions_of_clusters3D(spe_alpha_hull, cluster_colname = "alpha_hull_cluster")

alpha_hull_min_distances <- calculate_minimum_distances_to_clusters3D(spe_alpha_hull, cluster_colname = "alpha_hull_cluster", 
                                                                      cell_types_inside_cluster = c("Tumour", "Immune"),
                                                                      cell_types_outside_cluster = c("Tumour"))

calculate_volume_of_clusters3D(spe_alpha_hull, cluster_colname = "alpha_hull_cluster")

calculate_center_of_clusters3D(spe_alpha_hull, "alpha_hull_cluster")

spe_alpha_hull <- calculate_border_of_clusters3D(spe_alpha_hull, 7, "alpha_hull_cluster")


spe_dbscan <- dbscan_clustering3D(spe1, c("Tumour", "Immune"), radius = 10, minimum_cells_in_radius = 20)

dbscan_props <- calculate_cell_proportions_of_clusters3D(spe_dbscan, cluster_colname = "dbscan_cluster")

dbscan_min_distances <- calculate_minimum_distances_to_clusters3D(spe_dbscan, cluster_colname = "dbscan_cluster", 
                                                                  cell_types_inside_cluster = c("Tumour", "Immune"),
                                                                  cell_types_outside_cluster = c("Tumour"))

calculate_volume_of_clusters3D(spe_dbscan, cluster_colname = "dbscan_cluster")

calculate_center_of_clusters3D(spe_dbscan, "dbscan_cluster")

spe_dbscan <- calculate_border_of_clusters3D(spe_dbscan, 6, "dbscan_cluster")


spe_grid <- grid_based_clustering3D(spe1, cell_types_of_interest = c("Tumour", "Immune"), n_splits = 10)

plot_grid_based_clusters3D(spe_grid, c("Tumour", "Immune", "Others"), c("orange", "skyblue", "lightgray"))

grid_props <- calculate_cell_proportions_of_clusters3D(spe_grid, cluster_colname = "grid_based_cluster")

grid_min_distances <- calculate_minimum_distances_to_clusters3D(spe_grid, cluster_colname = "grid_based_cluster", 
                                                                cell_types_inside_cluster = c("Tumour", "Immune"),
                                                                cell_types_outside_cluster = c("Tumour"))

calculate_volume_of_clusters3D(spe_grid, cluster_colname = "grid_based_cluster")

calculate_center_of_clusters3D(spe_grid, "grid_based_cluster")

spe_grid <- calculate_border_of_clusters3D(spe_grid, 8, "grid_based_cluster")


### 5. Presence of cluster metrics --------------------------------------------
ANNI_result <- average_nearest_neighbor_index3D(data,
                                                cell_types_of_interest = c("Tumour", "Immune"),
                                                n_simulations = 10)

library(ggplot2)
df <- data.frame(x = 0,
                 y = ANNI_result$ANNI_mean,
                 low = ANNI_result$ANNI_95confidence_interval[1],
                 up = ANNI_result$ANNI_95confidence_interval[2])

ggplot(df, aes(x, y)) + geom_point() + geom_errorbar(aes(ymin = low, ymax = up))



### 6. Plot data -------------------------------------------------------------
plot_cells3D(spe1,
             plot_cell_types = c("Tumour", "Immune", "Others"),
             plot_colours = c("orange", "skyblue", "lightgray"))
