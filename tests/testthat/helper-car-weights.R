# Neighbourhood graph for the CAR models (besag, bym, bym2, leroux): the 4 nearest
# centroids of the 42 areas of mys / sim_area. mys_proxmat links every pair of
# areas, which a CAR model reads as a complete graph.
mys_knn <- create_weights(cbind(sim_area$x_coord, sim_area$y_coord), method = "knn", k = 4,
                          style = "B")
