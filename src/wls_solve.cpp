#include <RcppEigen.h>
# include <omp.h>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

using Eigen::VectorXd;
using Eigen::LeastSquaresConjugateGradient;
using Eigen::SparseMatrix;

// this code was produced by OpenAI, https://beta.openai.com/playground?model=code-davinci-002
// see also https://stackoverflow.com/questions/26160815/mappedsparsematrix-in-rcppeigen

// [[Rcpp::export]]
Rcpp::List wls_solve(Eigen::SparseMatrix<double> X, Eigen::VectorXd y, Eigen::VectorXd w) { // , int n_cores = 2
  // Eigen::setNbThreads(n_cores);
  Eigen::VectorXd weights = w.array().sqrt();
  Eigen::SparseMatrix<double> X_weighted = X;
  for (int k = 0; k < X.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(X, k); it; ++it) {
      it.valueRef() *= weights(it.row());
    }
  }
  Eigen::VectorXd y_weighted = y.cwiseProduct(weights);
  Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > solver;
  solver.setMaxIterations(100);
  solver.setTolerance(1E-20);
  solver.compute(X_weighted);
  VectorXd                        res = solver.solve(y_weighted);
  return Rcpp::List::create(Rcpp::Named("res") = res);
}

/*** R
library(Matrix)
set.seed(1234)
n <- 1000
p <- 500
X <- Matrix(rnorm(n * p), nrow = n)
# X[sample(1:n * p, size = n * p * 0.8)] <- 0
X <- as(X, "sparseMatrix")
y <- rnorm(n)
w <- 1/(y+0.1)
wls_solve(X, y, w)$res
*/