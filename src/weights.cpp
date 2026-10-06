// src/weights.cpp
// C++ helpers for create_weights bottlenecks (ponytail: only for n>500, else R vectorized is fine)
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// knn binary adjacency (directed, before symmetrization)
// D: symmetric distance matrix (n x n), diag = 0
// [[Rcpp::export(.knn_adj_cpp)]]
arma::mat knn_adj_cpp(const arma::mat& D, int k) {
  int n = D.n_rows;
  if (k < 1) k = 1;
  if (k >= n) k = n - 1;
  arma::mat A(n, n, fill::zeros);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  // ponytail: pre-check interrupt for large n (n>5k non-interruptible otherwise)
  if (n > 5000) Rcpp::checkUserInterrupt();
  for (int i = 0; i < n; ++i) {
    // copy row — alloc per row is needed for sort_index; keep O(n) not O(n^2) per thread
    arma::vec d(n);
    for (int j = 0; j < n; ++j) d(j) = D(i, j);
    d(i) = arma::datum::inf;
    arma::uvec idx = arma::sort_index(d); // ascending
    for (int t = 0; t < k; ++t) {
      int j = idx(t);
      A(i, j) = 1.0;
    }
    if ((i % 2048 == 0) && n > 5000) Rcpp::checkUserInterrupt();
  }
  return A;
}

// distance band binary adjacency
// [[Rcpp::export(.dist_adj_cpp)]]
arma::mat dist_adj_cpp(const arma::mat& D, double d_max) {
  int n = D.n_rows;
  arma::mat A(n, n, fill::zeros);
  if (n > 5000) Rcpp::checkUserInterrupt();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) continue;
      double d = D(i, j);
      if (d > 0 && d <= d_max) A(i, j) = 1.0;
    }
  }
  return A;
}

// inverse distance weight matrix (vectorized equivalent, but in C++ for n>3000)
// [[Rcpp::export(.idw_mat_cpp)]]
arma::mat idw_mat_cpp(const arma::mat& D, double alpha) {
  int n = D.n_rows;
  arma::mat W(n, n, fill::zeros);
  if (n > 5000) Rcpp::checkUserInterrupt();
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      if (i == j) continue;
      double d = D(i, j);
      if (d > 0 && std::isfinite(d)) W(i, j) = 1.0 / std::pow(d, alpha);
    }
  }
  return W;
}
