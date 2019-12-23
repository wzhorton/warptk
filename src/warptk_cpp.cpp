#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export(".dmat_C")]]
arma::mat dmat (arma::vec x) {
  int n = x.n_elem;
  arma::mat dmat(n, n);
  double dist;

  for(int i = 0; i < n; i++){
    for(int j = 0; j <= i; j++){
      dist = abs(x(i) - x(j));
      dmat(i,j) = dist;
      dmat(j,i) = dist;
    }
  }
  return dmat;
}
