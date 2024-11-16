# Cases to test: ------
# 1. Correct spe input
md1 <- spe_metadata_background_template("random")
md1 <- spe_metadata_cluster_template("regular", "sphere", md1)
md1 <- spe_metadata_cluster_template("regular", "ellipsoid", md1)
md1 <- spe_metadata_cluster_template("regular", "cylinder", md1)
spe1 <- simulate_spe_metadata3D(md1)
spe1$Cell.Type[spe1$Cell.Type == "Endothelial"] <- "C"
spe1$Cell.Type[spe1$Cell.Type == "Tumour"] <- "A"
spe1$Cell.Type[spe1$Cell.Type == "Immune"] <- "B"
spe1$Cell.Type[spe1$Cell.Type == "Others"] <- "O"

# 2. Not spe object input
spe2 <- list()

# 3. Spe object input, but with invalid feature colname
md3 <- spe_metadata_background_template("random")
spe3 <- simulate_spe_metadata3D(md3)
names(colData(spe3))[1] <- "dummy"

# 4. Zero cells
md4 <- spe_metadata_background_template("random")
md4$background$cell_types <- c("A", "B")
spe4 <- simulate_spe_metadata3D(md4)
spe4 <- spe4[ , 1]

# 5. One cell
md5 <- spe_metadata_background_template("random")
md5$background$n_cells <- 1
md5$background$cell_types <- c("A", "B")
spe5 <- simulate_spe_metadata3D(md5)

# 6. Two cells
md6 <- spe_metadata_background_template("random")
md6$background$n_cells <- 2
md6$background$cell_types <- c("A", "B")
spe6 <- simulate_spe_metadata3D(md6)

# 7. Many cells, but 0 cell type of interest/reference/target
md7 <- spe_metadata_background_template("random")
md7$background$cell_types <- c("C")
md7$background$cell_proportions <- c(1)
spe7 <- simulate_spe_metadata3D(md7)

# 8. Many cells, but 1 cell type of interest/reference/target
md8 <- spe_metadata_background_template("random")
md8$background$cell_types <- c("C")
md8$background$cell_proportions <- c(1)
spe8 <- simulate_spe_metadata3D(md8)
spe8[["Cell.Type"]][1] <- "A"

# 9. Many cells, but 2 cell type of interest/reference/target
md9 <- spe_metadata_background_template("random")
md9$background$cell_types <- c("C")
md9$background$cell_proportions <- c(1)
spe9 <- simulate_spe_metadata3D(md9)
spe9[["Cell.Type"]][1] <- "A"
spe9[["Cell.Type"]][1] <- "B"

# 10. One cell type only
md10 <- spe_metadata_background_template("random")
md10$background$cell_types <- c("C")
md10$background$cell_proportions <- c(1)
spe10 <- simulate_spe_metadata3D(md10)

# 11. Reference cell type and target cell type are separated.
md11 <- spe_metadata_background_template("random")
md11$background$cell_types <- c("C")
md11$background$cell_proportions <- c(1)
spe11 <- simulate_spe_metadata3D(md11)
spe11$Cell.Type[spatialCoords(spe11[ , 3]) < 20] <- "A"
spe11$Cell.Type[spatialCoords(spe11[ , 3]) > 280] <- "B"


# Testing functions -----
chosen_spe <- spe1


cell_props1 <- calculate_cell_proportions3D(chosen_spe,
                                            cell_types_of_interest = NULL,
                                            plot_image = TRUE)
print(cell_props1)

cell_props2 <- calculate_cell_proportions3D(chosen_spe,
                                            cell_types_of_interest = c("A", "B"),
                                            plot_image = TRUE)
print(cell_props2)

cell_props3 <- calculate_cell_proportions3D(chosen_spe,
                                            cell_types_of_interest = c("A"),
                                            plot_image = TRUE)

print(cell_props3)

cell_props4 <- calculate_cell_proportions3D(chosen_spe,
                                            cell_types_of_interest = c("A", "D"),
                                            plot_image = TRUE)

print(cell_props4)



### Calculate Pairwise Distances between Cells
pairwise_distances <- calculate_pairwise_distances_between_cell_types3D(chosen_spe,
                                                                        cell_types_of_interest = c("A", "B"),
                                                                        plot_image = TRUE)


### Calculate Minimum Distances between cells
minimum_distances <- calculate_minimum_distances_between_cell_types3D(chosen_spe,
                                                                      cell_types_of_interest = c("A", "B"),
                                                                      plot_image = TRUE)





### Calculate Mixing Scores
mixing_scores <- calculate_mixing_scores3D(chosen_spe,
                                           reference_cell_types = c("A", "B"),
                                           target_cell_types = c("A", "B"),
                                           radius = 20)
print(mixing_scores)

mixing_scores <- calculate_mixing_scores3D(chosen_spe,
                                           reference_cell_types = c("A", "B", "D"),
                                           target_cell_types = c("A", "B"),
                                           radius = 20)
print(mixing_scores)

mixing_scores_gradient <- calculate_mixing_scores_gradient3D(chosen_spe,
                                                             reference_cell_type = "A",
                                                             target_cell_type = "B",
                                                             radii = seq(10, 200, 10))

mixing_scores_gradient <- calculate_mixing_scores_gradient3D(chosen_spe,
                                                             reference_cell_type = "B",
                                                             target_cell_type = "A",
                                                             radii = seq(10, 200, 10))


### Calculate cells in the neighbourhood
neighbourhood_cells <- calculate_cells_in_neighbourhood3D(chosen_spe,
                                                          reference_cell_type = "A",
                                                          target_cell_types = c("A", "B"),
                                                          radius = 30,
                                                          plot_image = F)

neighbourhood_cells_gradient <- calculate_cells_in_neighbourhood_gradient3D(chosen_spe,
                                                                            reference_cell_type = "A",
                                                                            target_cell_types = c("A", "B"),
                                                                            radii = seq(1, 30, 2),
                                                                            plot_image = T)


## Calculate cell proportions in the neighbourhood
neighbourhood_cell_proportions <- calculate_cells_in_neighbourhood_proportions3D(chosen_spe,
                                                                                 reference_cell_type = "A",
                                                                                 target_cell_types = c("A", "B"),
                                                                                 radius = 20)
print(neighbourhood_cell_proportions)

neighbourhood_cell_proportions_gradient <- calculate_cells_in_neighbourhood_proportions_gradient3D(chosen_spe,
                                                                                                   reference_cell_type = "A",
                                                                                                   target_cell_types = c("A", "B"),
                                                                                                   radii = seq(1, 50, 3))


### Calculate cross-K function
cross_K <- calculate_cross_K3D(chosen_spe,
                               reference_cell_type = "A",
                               target_cell_type = "B",
                               radius = 20)
print(cross_K)


cross_K_gradient <- calculate_cross_K_gradient3D(chosen_spe,
                                                 reference_cell_type = "A",
                                                 target_cell_type = "B",
                                                 radii = seq(1, 50))


### Calculate entropy
entropy_background <- calculate_entropy_background3D(chosen_spe,
                                                     cell_types_of_interest = c("A", "B"))

print(entropy_background)

entropy_result <- calculate_entropy3D(chosen_spe,
                                      radius = 20,
                                      reference_cell_type = "A",
                                      target_cell_types = c("A", "B"))





entropy_gradient <- calculate_entropy_gradient3D(chosen_spe,
                                                 reference_cell_type = "A",
                                                 target_cell_types = c("A", "B"),
                                                 radii = seq(1, 50, 2),
                                                 plot_image = TRUE)


### Using all_single_radius and all_gradient functions

all_single_radius_result <- calculate_all_single_radius_cc_metrics3D(chosen_spe, "A", c("A", "B"), 20)

all_gradient_result <- calculate_all_gradient_cc_metrics3D(chosen_spe, "A", c("A", "B"), seq(1, 50, 2))



### Calculate entropy grid metrics
entropy_grid_metrics <- calculate_entropy_grid_metrics3D(chosen_spe,
                                                         n_splits = 8,
                                                         cell_types_of_interest = c("A", "B", "B1"),
                                                         plot_image = TRUE)
plot_grid_metrics_discrete3D(entropy_grid_metrics, "entropy")


### Calculate entropy prevalence
entropy_prevalence <- calculate_prevalence3D(entropy_grid_metrics,
                                             metric_colname = "entropy",
                                             threshold = 0.5)
print(entropy_prevalence)

entropy_prevalence_gradient <- calculate_prevalence_gradient3D(entropy_grid_metrics,
                                                               "entropy")

### Calculate spatial autocorrelation
entropy_spatial_autocorrelation <- calculate_spatial_autocorrelation3D(entropy_grid_metrics,
                                                                       metric_colname = "entropy",
                                                                       weight_method = "IDW")
print(entropy_spatial_autocorrelation)


### Calculate cell proportion grid metrics
cell_proportion_grid_metrics <- calculate_cell_proportion_grid_metrics3D(chosen_spe,
                                                                         n_splits = 10,
                                                                         reference_cell_types = c("A"),
                                                                         target_cell_types = c("B"),
                                                                         plot_image = TRUE)
plot_grid_metrics_discrete3D(cell_proportion_grid_metrics, "proportion")


### Calculate cell proportion prevalence
cell_proportion_prevalence <- calculate_prevalence3D(cell_proportion_grid_metrics,
                                                     metric_colname = "proportion",
                                                     threshold = 0.5)
print(cell_proportion_prevalence)

cell_proportion_prevalence_gradient <- calculate_prevalence_gradient3D(cell_proportion_grid_metrics,
                                                                       metric_colname = "proportion")

## Calculate spatial autocorrelation for cell proportions
cell_proportion_spatial_autocorrelation <- calculate_spatial_autocorrelation3D(cell_proportion_grid_metrics,
                                                                               metric_colname = "proportion",
                                                                               weight_method = 0.10)
print(cell_proportion_spatial_autocorrelation)




spe_alpha_hull <- alpha_hull_clustering3D(chosen_spe, c("A", "B"), alpha = 15, minimum_cells_in_alpha_hull = 30)

plot_alpha_hull_clusters3D(spe_alpha_hull, c("A", "B", "C"), c("orange", "skyblue", "lightgray"))

alpha_hull_props <- calculate_cell_proportions_of_clusters3D(spe_alpha_hull, cluster_colname = "alpha_hull_cluster")

alpha_hull_min_distances <- calculate_minimum_distances_to_clusters3D(spe_alpha_hull, cluster_colname = "alpha_hull_cluster", 
                                                                      cell_types_inside_cluster = c("A", "B"),
                                                                      cell_types_outside_cluster = c("A"))

calculate_volume_of_clusters3D(spe_alpha_hull, cluster_colname = "alpha_hull_cluster")

calculate_center_of_clusters3D(spe_alpha_hull, "alpha_hull_cluster")

spe_alpha_hull <- calculate_border_of_clusters3D(spe_alpha_hull, 15, "alpha_hull_cluster")


spe_dbscan <- dbscan_clustering3D(chosen_spe, c("A", "B"), radius = 30, minimum_cells_in_radius = 20, minimum_cells_in_cluster = 30)

dbscan_props <- calculate_cell_proportions_of_clusters3D(spe_dbscan, cluster_colname = "dbscan_cluster")

dbscan_min_distances <- calculate_minimum_distances_to_clusters3D(spe_dbscan, cluster_colname = "dbscan_cluster", 
                                                                  cell_types_inside_cluster = c("A", "B"),
                                                                  cell_types_outside_cluster = c("A"))

calculate_volume_of_clusters3D(spe_dbscan, cluster_colname = "dbscan_cluster")

calculate_center_of_clusters3D(spe_dbscan, "dbscan_cluster")

spe_dbscan <- calculate_border_of_clusters3D(spe_dbscan, 6, "dbscan_cluster")


spe_grid <- grid_based_clustering3D(chosen_spe, cell_types_of_interest = c("A", "B"), n_splits = 10, minimum_cells_in_cluster = 30)

plot_grid_based_clusters3D(spe_grid, c("A", "B", "C"), c("orange", "skyblue", "lightgray"))

grid_props <- calculate_cell_proportions_of_clusters3D(spe_grid, cluster_colname = "grid_based_cluster")

grid_min_distances <- calculate_minimum_distances_to_clusters3D(spe_grid, cluster_colname = "grid_based_cluster", 
                                                                cell_types_inside_cluster = c("A", "B"),
                                                                cell_types_outside_cluster = c("A"))

calculate_volume_of_clusters3D(spe_grid, cluster_colname = "grid_based_cluster")

calculate_center_of_clusters3D(spe_grid, "grid_based_cluster")

spe_grid <- calculate_border_of_clusters3D(spe_grid, 8, "grid_based_cluster")



plot_cells3D(chosen_spe,
             plot_cell_types = c("A", "B", "C"),
             plot_colours = c("orange", "skyblue", "lightgray"))
