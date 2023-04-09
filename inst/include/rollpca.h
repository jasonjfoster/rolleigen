#ifndef ROLLPCA_H
#define ROLLPCA_H

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <roll.h>
using namespace Rcpp;
using namespace RcppParallel;

namespace rollpca {

// 'Worker' function for rolling eigenvalues and eigenvectors
struct RollEigenSlices : public Worker {
  
  const arma::cube arma_cov;      // source
  const int n_rows;
  const int n_cols;
  arma::mat& arma_eigen_values;   // destination (pass by reference)
  arma::cube& arma_eigen_vectors; // destination (pass by reference)
  
  // initialize with source and destination
  RollEigenSlices(const arma::cube arma_cov, const int n_rows,
                  const int n_cols, arma::mat& arma_eigen_values,
                  arma::cube& arma_eigen_vectors)
    : arma_cov(arma_cov), n_rows(n_rows),
      n_cols(n_cols), arma_eigen_values(arma_eigen_values),
      arma_eigen_vectors(arma_eigen_vectors) { }
  
  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      arma::mat x = arma_cov.slice(i);
      arma::mat A = x.submat(0, 0, n_cols - 1, n_cols - 1);
      arma::vec eigen_values(n_cols);
      arma::mat eigen_vectors(n_cols, n_cols);
      
      // check if missing value is present
      bool any_na = x.has_nan();
      
      // don't compute if missing value 
      if (!any_na) {
        
        // check if solution is found
        bool status = arma::eig_sym(eigen_values, eigen_vectors, A);
        
        // don't find approximate solution for rank deficient system
        if (status) {
          
          // reverse order for consistency with R's eigen
          std::reverse(eigen_values.begin(), eigen_values.end());
          eigen_vectors = arma::fliplr(eigen_vectors);
          
          arma_eigen_values.row(i) = trans(eigen_values);
          arma_eigen_vectors.slice(i) = eigen_vectors;
          
        } else if (!status) {
          
          arma::vec no_solution_row(n_cols);
          no_solution_row.fill(NA_REAL);
          
          arma::mat no_solution_slice(n_cols, n_cols);
          no_solution_slice.fill(NA_REAL);
          
          arma_eigen_values.row(i) = trans(no_solution_row);
          arma_eigen_vectors.slice(i) = no_solution_slice;
          
        }
        
      } else {
        
        arma::vec no_solution_row(n_cols);
        no_solution_row.fill(NA_REAL);
        
        arma::mat no_solution_slice(n_cols, n_cols);
        no_solution_slice.fill(NA_REAL);
        
        arma_eigen_values.row(i) = trans(no_solution_row);
        arma_eigen_vectors.slice(i) = no_solution_slice;
        
      }
      
    }
  }
  
};

}

#endif