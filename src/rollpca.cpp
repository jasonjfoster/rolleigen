#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <roll.h>
using namespace Rcpp;
using namespace RcppParallel;

void check_width(const int& width) {
  
  if (width < 1) {
    stop("value of 'width' must be greater than zero");
  }
  
}

void check_weights_x(const int& n_rows_x, const int& width,
                     const arma::vec& weights) {
  
  if ((int)weights.size() < std::min(width, n_rows_x)) {
    stop("length of 'weights' must equal either the number of rows in 'x' or 'width'");
  }
  
}

bool check_lambda(const arma::vec& weights, const int& n_rows_x,
                  const int& width, const bool& online) {
  
  // check if equal-weights
  bool status_eq = all(weights == weights[0]);
  bool status_exp = true;
  
  // check if exponential-weights
  if (!status_eq) {
    
    int i = 0;
    int n = weights.size();
    long double lambda = 0;
    long double lambda_prev = 0;
    
    // check if constant ratio
    while (status_exp && (i <= (n - 2))) {
      
      // ratio of weights
      lambda_prev = lambda;
      lambda = weights[n - i - 2] / weights[n - i - 1];
      
      // tolerance for consistency with R's all.equal
      if (((i > 0) && (std::abs(lambda - lambda_prev) > sqrt(arma::datum::eps))) ||
          ((weights[n - i - 2] > weights[n - i - 1]) && (width < n_rows_x)) ||
          (std::isnan(lambda) || (std::isinf(lambda)))) {
        
        status_exp = false;
        
      }
      
      i += 1;
      
    }
    
  }
  
  if (!status_exp && online) {
    warning("'online' is only supported for equal or exponential decay 'weights'");
  }
  
  return status_exp;
  
}

void check_min_obs(const int& min_obs) {
  
  if (min_obs < 1) {
    stop("value of 'min_obs' must be greater than zero");
  }
  
}

arma::uvec any_na_x(const NumericMatrix& x) {
  
  int n_rows_x = x.nrow();
  int n_cols_x = x.ncol();
  arma::uvec result(n_rows_x);
  
  for (int i = 0; i < n_rows_x; i++) {
    
    int any_na = 0;
    int j = 0;
    
    while ((any_na == 0) && (j < n_cols_x)) {
      if (std::isnan(x(i, j))) {
        any_na = 1;
      }
      j += 1;
    }
    
    result[i] = any_na;
    
  }
  
  return result;
  
}

arma::cube roll_cov_z(const NumericMatrix& x, const int& width,
                      const arma::vec& weights, const bool& center,
                      const bool& scale, const int& min_obs,
                      const bool& complete_obs, const bool& na_restore,
                      const bool& online) {
  
  int n = weights.size();
  int n_rows = x.nrow();
  int n_cols = x.ncol();
  arma::uvec arma_any_na(n_rows);
  arma::cube arma_cov(n_cols, n_cols, n_rows);
  
  // check 'width' argument for errors
  check_width(width);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  check_weights_x(n_rows, width, weights);
  bool status = check_lambda(weights, n_rows, width, online);
  
  // default 'min_obs' argument is 'width',
  // otherwise check argument for errors
  check_min_obs(min_obs);
  
  // default 'complete_obs' argument is 'true',
  // otherwise check argument for errors
  if (complete_obs) {
    arma_any_na = any_na_x(x);
  } else {
    arma_any_na.fill(0);
  }
  
  // compute rolling covariances
  if (status && online) {
    
    roll::RollCovOnlineMatXX roll_cov_online(x, n, n_rows, n_cols, width,
                                             weights, center, scale, min_obs,
                                             arma_any_na, na_restore,
                                             arma_cov);
    parallelFor(0, n_cols, roll_cov_online);
    
  } else {
    
    roll::RollCovOfflineMatXX roll_cov_offline(x, n, n_rows, n_cols, width,
                                               weights, center, scale, min_obs,
                                               arma_any_na, na_restore,
                                               arma_cov);
    parallelFor(0, n_rows * n_cols * (n_cols + 1) / 2, roll_cov_offline);
    
  }
  
  return arma_cov;
  
}

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

// [[Rcpp::export(.roll_eigen)]]
List roll_eigen(const NumericMatrix& x, const int& width,
                const arma::vec& weights, const bool& center,
                const bool& scale, const int& min_obs,
                const bool& complete_obs, const bool& na_restore,
                const bool& online) {
  
  int n_rows = x.nrow();  
  int n_cols = x.ncol();
  arma::cube arma_cov = roll_cov_z(x, width, weights, center, scale, min_obs, complete_obs, 
                                   na_restore, online);
  arma::mat arma_eigen_values(n_rows, n_cols);
  arma::cube arma_eigen_vectors(n_cols, n_cols, n_rows);
  
  // compute rolling eigenvalues and eigenvectors
  RollEigenSlices roll_eigen_slices(arma_cov, n_rows, n_cols,
                                    arma_eigen_values, arma_eigen_vectors);
  parallelFor(0, n_rows, roll_eigen_slices);
  
  // create and return a matrix or xts object for eigenvalues
  NumericMatrix eigen_values(wrap(arma_eigen_values));
  List dimnames = x.attr("dimnames");
  if (dimnames.size() > 1) {
    eigen_values.attr("dimnames") = List::create(dimnames[0], R_NilValue);
  }
  eigen_values.attr("index") = x.attr("index");
  eigen_values.attr(".indexCLASS") = x.attr(".indexCLASS");
  eigen_values.attr(".indexTZ") = x.attr(".indexTZ");
  eigen_values.attr("tclass") = x.attr("tclass");
  eigen_values.attr("tzone") = x.attr("tzone");
  eigen_values.attr("class") = x.attr("class");
  
  // create and return a cube for eigenvectors
  NumericVector eigen_vectors(wrap(arma_eigen_vectors));
  eigen_vectors.attr("dim") = IntegerVector::create(n_cols, n_cols, n_rows);
  if (dimnames.size() > 1) {
    eigen_vectors.attr("dimnames") = List::create(dimnames[1], R_NilValue);
  }
  
  // create and return a list
  List result = List::create(Named("values") = eigen_values,
                             Named("vectors") = eigen_vectors);
  
  return result;
  
}