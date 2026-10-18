# Builds data/mys_proxmat.rda: the 4-nearest-neighbour graph of the 42 areas of
# mys / sim_area, made symmetric and row-standardised.
#
# The centroids are those of sim_area (x_coord, y_coord). Each area is linked to
# its 4 nearest areas, and a link in either direction counts, so every area has
# 4 or more neighbours and the graph is connected. Row-standardised weights suit
# the SAR models (eblup_sfh(), eblup_stfh(), hb_area(spatial = "slm")); the CAR
# models of hb_area() use the neighbour pattern W > 0.
#
# Until fastsae 1.0.0, mys_proxmat was a dense distance-decay matrix that linked
# every pair of areas (a complete graph for the CAR models).

load("data/sim_area.rda")
coords <- as.matrix(sim_area[, c("x_coord", "y_coord")])
D <- nrow(coords)
k <- 4

dist_mat <- as.matrix(stats::dist(coords))
adj <- matrix(0, D, D)
for (i in seq_len(D)) {
  adj[i, order(dist_mat[i, ])[seq_len(k) + 1L]] <- 1
}
adj <- pmax(adj, t(adj))

reach <- 1L
repeat {
  nxt <- union(reach, which(colSums(adj[reach, , drop = FALSE]) > 0))
  if (length(nxt) == length(reach)) break
  reach <- nxt
}
stopifnot(isSymmetric(adj), all(diag(adj) == 0), length(reach) == D)

mys_proxmat <- adj / rowSums(adj)
dimnames(mys_proxmat) <- list(as.character(seq_len(D)), as.character(seq_len(D)))

save(mys_proxmat, file = "data/mys_proxmat.rda", compress = "bzip2")
