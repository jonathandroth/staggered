#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/Dense>
// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;

  return Rcpp::wrap(C);
}


// [[Rcpp::export]]
Eigen::MatrixXd solve_least_squares_svd(
    Eigen::MatrixXd A,
    Eigen::MatrixXd B
) {
  Eigen::MatrixXd x = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(B);
  return(x);
}


// [[Rcpp::export]]
Eigen::MatrixXd solve_least_squares_normal(
    Eigen::MatrixXd A,
    Eigen::MatrixXd B
) {
  Eigen::MatrixXd x = (A.transpose() * A).ldlt().solve(A.transpose() * B);
  return(x);
}

