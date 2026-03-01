#ifndef ROLLEIGEN_H
#define ROLLEIGEN_H

#define ARMA_WARN_LEVEL 0

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <roll.h>
#include "rolleigen_check.h"
#include "rolleigen_attr.h"
using namespace Rcpp;
using namespace RcppParallel;

namespace rolleigen {

// 'Worker' function for computing the rolling statistic using a standard algorithm
struct RollOrderSlices : public Worker {
  
  const int n_rows_x;
  const int n_cols_x;
  arma::mat& arma_eigen_values;
  arma::cube& arma_eigen_vectors;
  
  RollOrderSlices(const int n_rows_x, const int n_cols_x,
                  arma::mat& arma_eigen_values, arma::cube& arma_eigen_vectors)
    : n_rows_x(n_rows_x), n_cols_x(n_cols_x),
      arma_eigen_values(arma_eigen_values), arma_eigen_vectors(arma_eigen_vectors) { }
  
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      if (i > 0) {

        arma::vec eigen_values = trans(arma_eigen_values.row(i));
        arma::mat eigen_vectors = arma_eigen_vectors.slice(i);
        arma::mat eigen_vectors0 = arma_eigen_vectors.slice(i - 1);
        
        // check if missing value is present
        bool any_na = eigen_vectors.has_nan();
        bool any_na0 = eigen_vectors0.has_nan();
        
        // don't compute if missing value 
        if (!any_na && !any_na0) {
          
          arma::mat similarity = trans(eigen_vectors) * eigen_vectors0;
          arma::uvec order = arma::index_max(arma::abs(similarity), 1);
          
          eigen_vectors = eigen_vectors.cols(order);
          similarity = similarity.cols(order);
          
          arma::vec signs = arma::sign(similarity.diag());
          
          arma_eigen_values.row(i) = trans(eigen_values(order));
          arma_eigen_vectors.slice(i) = eigen_vectors.each_col() % signs;
        
        }

      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using a standard algorithm
struct RollEigenSlices : public Worker {
  
  const arma::cube arma_cov;      // source
  const int n_rows_x;
  const int n_cols_x;
  arma::mat& arma_eigen_values;   // destination (pass by reference)
  arma::cube& arma_eigen_vectors;
  
  // initialize with source and destination
  RollEigenSlices(const arma::cube arma_cov, const int n_rows_x,
                  const int n_cols_x, arma::mat& arma_eigen_values,
                  arma::cube& arma_eigen_vectors)
    : arma_cov(arma_cov), n_rows_x(n_rows_x),
      n_cols_x(n_cols_x), arma_eigen_values(arma_eigen_values),
      arma_eigen_vectors(arma_eigen_vectors) { }
  
  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {
      
      arma::mat sigma = arma_cov.slice(i);
      
      // check if missing value is present
      bool any_na = sigma.has_nan();
      
      // don't compute if missing value 
      if (!any_na) {

          arma::mat A = sigma.submat(0, 0, n_cols_x - 1, n_cols_x - 1);
          arma::vec eigen_values(n_cols_x);
          arma::mat eigen_vectors(n_cols_x, n_cols_x);
        
        // check if solution is found
        bool status = arma::eig_sym(eigen_values, eigen_vectors, A);
        
        // don't find approximate solution for rank deficient system
        if (status) {
          
          // reverse order for consistency with R's eigen
          std::reverse(eigen_values.begin(), eigen_values.end());
          eigen_vectors = arma::fliplr(eigen_vectors);
          
          arma_eigen_values.row(i) = trans(eigen_values);
          arma_eigen_vectors.slice(i) = eigen_vectors;
          
        } else {
          
          arma_eigen_values.row(i).fill(NA_REAL);
          arma_eigen_vectors.slice(i).fill(NA_REAL);
          
        }
        
      } else {
        
        arma_eigen_values.row(i).fill(NA_REAL);
        arma_eigen_vectors.slice(i).fill(NA_REAL);
        
      }
      
    }
  }
  
};

// 'Worker' function for computing the rolling statistic using a standard algorithm
struct RollPcrInterceptTRUE : public Worker {
  
  const arma::cube arma_cov;            // source
  const int n_rows_xy;
  const int n_cols_x;
  const arma::uvec arma_cols;
  const arma::uvec arma_comps;
  const arma::vec arma_n_obs;
  const arma::mat arma_mean;
  const arma::cube arma_eigen_vectors;
  arma::mat& arma_coef;                 // destination (pass by reference)
  arma::vec& arma_rsq;
  
  // initialize with source and destination
  RollPcrInterceptTRUE(const arma::cube arma_cov, const int n_rows_xy,
                       const int n_cols_x, const arma::uvec arma_cols,
                       const arma::uvec arma_comps, const arma::vec arma_n_obs,
                       const arma::mat arma_mean, const arma::cube arma_eigen_vectors,
                       arma::mat& arma_coef, arma::vec& arma_rsq)
    : arma_cov(arma_cov), n_rows_xy(n_rows_xy),
      n_cols_x(n_cols_x), arma_cols(arma_cols),
      arma_comps(arma_comps), arma_n_obs(arma_n_obs),
      arma_mean(arma_mean), arma_eigen_vectors(arma_eigen_vectors),
      arma_coef(arma_coef), arma_rsq(arma_rsq) { }
  
  // function call operator that iterates by slice
  void operator()(std::size_t begin_slice, std::size_t end_slice) {
    for (std::size_t i = begin_slice; i < end_slice; i++) {

      arma::mat sigma = arma_cov.slice(i);
      
      // check if missing value is present
      bool any_na = sigma.has_nan();
      
      // don't compute if missing value 
      if (!any_na) {

        arma::mat A = sigma.submat(0, 0, n_cols_x - 2, n_cols_x - 2);
        arma::mat b = sigma.submat(0, n_cols_x - 1, n_cols_x - 2, n_cols_x - 1);
        arma::vec gamma(n_cols_x - 1);
        arma::mat eigen_vectors = arma_eigen_vectors.slice(i);
        
        // check if solution is found
        bool status_solve = arma::solve(gamma, A * eigen_vectors, b, arma::solve_opts::no_approx);
        int df_fit = n_cols_x;
        
        // don't find approximate solution for rank deficient system,
        // and the width and current row must be greater than the
        // number of variables
        if (status_solve && (arma_n_obs[i] >= df_fit)) {
          
          // coefficients
          arma::vec gamma_subset = gamma(arma_comps - 1);
          arma::vec coef = eigen_vectors.submat(arma_cols, arma_comps - 1) * gamma_subset;
          arma_coef.submat(i, 1, i, n_cols_x - 1) = trans(coef);
          
          // intercept
          arma::vec mean_x = trans(arma_mean.submat(i, 0, i, n_cols_x - 2));
          arma_coef(i, 0) = arma_mean(i, n_cols_x - 1) - arma::dot(mean_x, coef);
          
          // r-squared
          long double var_y = sigma(n_cols_x - 1, n_cols_x - 1);
          if (var_y > arma::datum::eps) {    
            arma_rsq[i] = arma::dot(coef, A * coef) / var_y;
          } else {
            arma_rsq[i] = NA_REAL;
          }
          
        } else {
          
          arma_coef.row(i).fill(NA_REAL);
          arma_rsq[i] = NA_REAL;
          
        }
        
      } else {

        arma_coef.row(i).fill(NA_REAL);
        arma_rsq[i] = NA_REAL;
        
      }
      
    }
  }
  
};

}

#endif