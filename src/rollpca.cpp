#include "rollpca.h"

void check_width(const int& width) {
  
  if (width < 1) {
    stop("value of 'width' must be greater than zero");
  }
  
}

void check_comps(const arma::uvec& comps, const unsigned int& n_cols) {
  
  if (comps.max() > n_cols) {
    stop("maximum value of 'n_comps' must be less than or equal to number of columns in 'x'");
  }
  
  if (comps.min() < 1) {
    stop("minimum value of 'n_comps' must be greater than or equal to one");
  }
  
  if (comps.size() > n_cols) {
    stop("length of 'n_comps' must be less than or equal to number of columns in 'x'");
  }
  
}

void check_weights_x(const int& n_rows_x, const int& width,
                     const arma::vec& weights) {
  
  if ((int)weights.size() < std::min(width, n_rows_x)) {
    stop("length of 'weights' must equal either the number of rows in 'x' or 'width'");
  }
  
}

void check_weights_lm(const int& n_rows_xy, const int& width,
                      const arma::vec& weights) {
  
  if ((int)weights.size() < std::min(width, n_rows_xy)) {
    stop("length of 'weights' must equal either the number of rows in 'x' (and 'y') or 'width'");
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

void check_lm(const int& n_rows_x, const int& n_rows_y) {
  
  if (n_rows_x != n_rows_y) {
    stop("number of rows in 'x' must equal the number of rows in 'y'");
  }
  
}

List dimnames_lm_x(const List& input, const int& n_cols_x,
                   const bool& intercept) {
  
  if (intercept && (input.size() > 1)) {
    
    CharacterVector dimnames_cols = input[1];
    CharacterVector result(n_cols_x);
    result(0) = "(Intercept)";
    
    std::copy(dimnames_cols.begin(), dimnames_cols.end(), result.begin() + 1);
    
    return List::create(input[0], result);
    
  } else if (!intercept && (input.size() > 1)) {
    
    return List::create(input[0], input[1]);
    
  } else if (intercept) {
    
    CharacterVector result(n_cols_x);
    result(0) = "(Intercept)";
    
    for (int i = 1; i < n_cols_x; i++) {
      
      result[i] = "x";
      result[i] += i;
      
    }
    
    return List::create(R_NilValue, result);
    
  } else {
    
    CharacterVector result(n_cols_x);
    
    for (int i = 0; i < n_cols_x; i++) {
      
      result[i] = "x";
      result[i] += i + 1;
      
    }
    
    return List::create(R_NilValue, result);
    
  }
  
}

CharacterVector dimnames_lm_y(const List& input, const int& n_cols_y) {
  
  if (input.size() > 1) {
    
    return input[1];
    
  } else {
    
    CharacterVector result(n_cols_y);
    
    for (int i = 0; i < n_cols_y; i++) {
      
      result[i] = "y";
      result[i] += i + 1;
      
    }
    
    return result;
    
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

arma::uvec seq(const int& size) {
  
  arma::uvec result(size);
  
  for (int i = 0; i < size; i++) {
    result[i] = i;
  }
  
  return result;
  
}

arma::cube roll_cov_z(const NumericMatrix& x, const int& width,
                      const arma::vec& weights, const bool& center,
                      const bool& scale, const int& min_obs,
                      const bool& complete_obs, const bool& na_restore,
                      const bool& online) {
  
  int n = weights.size();
  int n_rows_x = x.nrow();
  int n_cols_x = x.ncol();
  arma::uvec arma_any_na(n_rows_x);
  arma::cube arma_cov(n_cols_x, n_cols_x, n_rows_x);
  
  // check 'width' argument for errors
  check_width(width);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  check_weights_x(n_rows_x, width, weights);
  bool status = check_lambda(weights, n_rows_x, width, online);
  
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
    
    roll::RollCovOnlineMatXX roll_cov_online(x, n, n_rows_x, n_cols_x, width,
                                             weights, center, scale, min_obs,
                                             arma_any_na, na_restore,
                                             arma_cov);
    parallelFor(0, n_cols_x, roll_cov_online);
    
  } else {
    
    roll::RollCovOfflineMatXX roll_cov_offline(x, n, n_rows_x, n_cols_x, width,
                                               weights, center, scale, min_obs,
                                               arma_any_na, na_restore,
                                               arma_cov);
    parallelFor(0, n_rows_x * n_cols_x * (n_cols_x + 1) / 2, roll_cov_offline);
    
  }
  
  return arma_cov;
  
}

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
  rollpca::RollEigenSlices roll_eigen_slices(arma_cov, n_rows, n_cols,
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

List roll_pcr_z(const NumericMatrix& x, const NumericVector& y,
                const int& width, const int& n_comps,
                const arma::vec& weights, const bool& intercept,
                const bool& center, const bool& scale,
                const int& min_obs, const bool& complete_obs,
                const bool& na_restore, const bool& online) {
  
  int n = weights.size();
  int n_rows_xy = x.nrow();
  int n_cols_x = x.ncol() + 1;
  arma::uvec arma_comps = seq(n_comps) + 1;
  arma::uvec arma_any_na(n_rows_xy);
  arma::vec arma_n_obs(n_rows_xy);
  arma::vec arma_sum_w(n_rows_xy);
  arma::mat arma_mean(n_rows_xy, n_cols_x);
  arma::cube arma_cov(n_cols_x, n_cols_x, n_rows_xy);
  arma::mat arma_eigen_values(n_rows_xy, n_cols_x - 1);
  arma::cube arma_eigen_vectors(n_cols_x - 1, n_cols_x - 1, n_rows_xy);
  arma::mat arma_coef(n_rows_xy, n_cols_x);
  arma::mat arma_rsq(n_rows_xy, 1);
  
  // check 'x' and 'y' arguments for errors
  check_lm(n_rows_xy, y.size());
  
  // check 'width' argument for errors
  check_width(width);
  
  // default 'n_comps' argument is all components
  check_comps(arma_comps, n_cols_x - 1);
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  check_weights_lm(n_rows_xy, width, weights);
  bool status = check_lambda(weights, n_rows_xy, width, online);
  
  // check 'intercept' argument for errors
  if (!intercept) {
    warning("'intercept = FALSE' is not supported");
  }
  
  // default 'min_obs' argument is 'width',
  // otherwise check argument for errors
  check_min_obs(min_obs);
  
  // cbind x and y variables
  NumericMatrix data(n_rows_xy, n_cols_x);
  std::copy(x.begin(), x.end(), data.begin());
  std::copy(y.begin(), y.end(), data.begin() + n_rows_xy * (n_cols_x - 1));
  
  // default 'complete_obs' argument is 'true',
  // otherwise check argument for errors
  if (complete_obs) {
    arma_any_na = any_na_x(data);
  } else {
    
    warning("'complete_obs = FALSE' is not supported");
    arma_any_na = any_na_x(data);
    
  }
  
  // compute rolling covariances
  if (status && online) {
    
    roll::RollCovOnlineMatLm roll_cov_online(data, n, n_rows_xy, n_cols_x, width,
                                             weights, true, min_obs,
                                             arma_any_na, na_restore,
                                             arma_n_obs, arma_sum_w, arma_mean,
                                             arma_cov);
    parallelFor(0, n_cols_x, roll_cov_online);
    
  } else {
    
    roll::RollCovOfflineMatLm roll_cov_offline(data, n, n_rows_xy, n_cols_x, width,
                                               weights, true, min_obs,
                                               arma_any_na, na_restore,
                                               arma_n_obs, arma_sum_w, arma_mean,
                                               arma_cov);
    parallelFor(0, n_rows_xy * n_cols_x * (n_cols_x + 1) / 2, roll_cov_offline);
    
  }
  
  // compute rolling eigenvalues and eigenvectors
  rollpca::RollEigenSlices roll_eigen_slices(arma_cov, n_rows_xy, n_cols_x - 1,
                                             arma_eigen_values, arma_eigen_vectors);
  parallelFor(0, n_rows_xy, roll_eigen_slices);
  
  // compute rolling principal component regressions
  arma::uvec arma_cols = seq(n_cols_x - 1);
  rollpca::RollPcrInterceptTRUE roll_pcr_slices(arma_cov, n_rows_xy, n_cols_x, arma_cols, arma_comps,
                                                arma_n_obs, arma_mean, arma_eigen_vectors,
                                                arma_coef, arma_rsq);
  parallelFor(0, n_rows_xy, roll_pcr_slices);
  
  // create and return a list
  List result = List::create(Named("coefficients") = arma_coef,
                             Named("r.squared") = arma_rsq);
  
  return result;
  
}

// [[Rcpp::export(.roll_pcr)]]
List roll_pcr(const SEXP& x, const SEXP& y,
              const int& width, const int& n_comps,
              const arma::vec& weights, const bool& intercept,
              const bool& center, const bool& scale,
              const int& min_obs, const bool& complete_obs,
              const bool& na_restore, const bool& online) {
  
  if (Rf_isMatrix(x) && Rf_isMatrix(y)) {
    
    NumericMatrix xx(x);
    NumericMatrix yy(y);
    int n_rows_xy = xx.nrow();  
    int n_cols_x = xx.ncol();
    int n_cols_y = yy.ncol();
    List result_coef(n_cols_y);
    List result_rsq(n_cols_y);
    // List result_se(n_cols_y);
    List result_z(2);
    List result(2);
    
    if (intercept) {
      n_cols_x += 1;
    }
    
    // create a list of matrices,
    // otherwise a list of lists
    if (n_cols_y == 1) {
      
      result_z = roll_pcr_z(xx, yy(_, 0), width, n_comps, 
                            weights, intercept, center, scale,
                            min_obs, complete_obs, na_restore,
                            online);
      
      arma::mat arma_coef_z = result_z[0];
      arma::mat arma_rsq_z = result_z[1];
      // arma::mat arma_se_z = result_z[2];
      
      // create and return a matrix or xts object for coefficients
      NumericVector coef(wrap(arma_coef_z));
      coef.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
      List x_dimnames = xx.attr("dimnames");
      coef.attr("dimnames") = dimnames_lm_x(x_dimnames, n_cols_x, intercept);
      coef.attr("index") = xx.attr("index");
      coef.attr(".indexCLASS") = xx.attr(".indexCLASS");
      coef.attr(".indexTZ") = xx.attr(".indexTZ");
      coef.attr("tclass") = xx.attr("tclass");
      coef.attr("tzone") = xx.attr("tzone");
      coef.attr("class") = xx.attr("class");
      
      // create and return a matrix or xts object for r-squareds
      NumericVector rsq(wrap(arma_rsq_z));
      rsq.attr("dim") = IntegerVector::create(n_rows_xy, 1);
      if (x_dimnames.size() > 1) {
        rsq.attr("dimnames") = List::create(x_dimnames[0], "R-squared");
      } else {
        rsq.attr("dimnames") = List::create(R_NilValue, "R-squared");
      }
      rsq.attr("index") = xx.attr("index");
      rsq.attr(".indexCLASS") = xx.attr(".indexCLASS");
      rsq.attr(".indexTZ") = xx.attr(".indexTZ");
      rsq.attr("tclass") = xx.attr("tclass");
      rsq.attr("tzone") = xx.attr("tzone");
      rsq.attr("class") = xx.attr("class");
      
      // // create and return a matrix or xts object for standard errors
      // NumericVector se(wrap(arma_se_z));
      // se.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
      // se.attr("dimnames") = coef.attr("dimnames");
      // se.attr("index") = xx.attr("index");
      // se.attr(".indexCLASS") = xx.attr(".indexCLASS");
      // se.attr(".indexTZ") = xx.attr(".indexTZ");
      // se.attr("tclass") = xx.attr("tclass");
      // se.attr("tzone") = xx.attr("tzone");
      // se.attr("class") = xx.attr("class");
      
      // create and return a list
      result = List::create(Named("coefficients") = coef,
                            Named("r.squared") = rsq);
      // Named("std.error") = se);
      
    } else {
      
      for (int z = 0; z < n_cols_y; z++) {
        
        result_z = roll_pcr_z(xx, yy(_, z), width, n_comps,
                              weights, intercept, center, scale,
                              min_obs, complete_obs, na_restore,
                              online);
        
        arma::mat arma_coef_z = result_z[0];
        arma::mat arma_rsq_z = result_z[1];
        // arma::mat arma_se_z = result_z[2];
        
        // create and return a matrix or xts object for coefficients
        NumericVector coef(wrap(arma_coef_z));
        coef.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
        List dimnames_x = xx.attr("dimnames");
        coef.attr("dimnames") = dimnames_lm_x(dimnames_x, n_cols_x, intercept);
        coef.attr("index") = xx.attr("index");
        coef.attr(".indexCLASS") = xx.attr(".indexCLASS");
        coef.attr(".indexTZ") = xx.attr(".indexTZ");
        coef.attr("tclass") = xx.attr("tclass");
        coef.attr("tzone") = xx.attr("tzone");
        coef.attr("class") = xx.attr("class");
        
        // create and return a matrix or xts object for r-squareds
        NumericVector rsq(wrap(arma_rsq_z));
        rsq.attr("dim") = IntegerVector::create(n_rows_xy, 1);
        if (dimnames_x.size() > 1) {
          rsq.attr("dimnames") = List::create(dimnames_x[0], "R-squared");
        } else {
          rsq.attr("dimnames") = List::create(R_NilValue, "R-squared");
        }
        rsq.attr("index") = xx.attr("index");
        rsq.attr(".indexCLASS") = xx.attr(".indexCLASS");
        rsq.attr(".indexTZ") = xx.attr(".indexTZ");
        rsq.attr("tclass") = xx.attr("tclass");
        rsq.attr("tzone") = xx.attr("tzone");
        rsq.attr("class") = xx.attr("class");
        
        // // create and return a matrix or xts object for standard errors
        // NumericVector se(wrap(arma_se_z));
        // se.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
        // se.attr("dimnames") = coef.attr("dimnames");
        // se.attr("index") = xx.attr("index");
        // se.attr(".indexCLASS") = xx.attr(".indexCLASS");
        // se.attr(".indexTZ") = xx.attr(".indexTZ");
        // se.attr("tclass") = xx.attr("tclass");
        // se.attr("tzone") = xx.attr("tzone");
        // se.attr("class") = xx.attr("class");
        
        result_coef(z) = coef;
        result_rsq(z) = rsq;
        // result_se(z) = se;
        
      }
      
      // add names to each list
      List dimnames_y = yy.attr("dimnames");
      result_coef.attr("names") = dimnames_lm_y(dimnames_y, n_cols_y);
      result_rsq.attr("names") = result_coef.attr("names");
      // result_se.attr("names") = result_coef.attr("names");
      
      // create and return a list
      result = List::create(Named("coefficients") = result_coef,
                            Named("r.squared") = result_rsq);
      // Named("std.error") = result_se);
      
    }
    
    return result;
    
  } else if (Rf_isMatrix(x)) {
    
    NumericMatrix xx(x);
    NumericVector yy(y);
    
    int n_rows_xy = xx.nrow();
    int n_cols_x = xx.ncol();
    List result_coef(1);
    List result_rsq(1);
    // List result_se(1);
    List result_z(2);
    List result(2);
    
    if (intercept) {
      n_cols_x += 1;
    }
    
    // create a list of matrices
    result_z = roll_pcr_z(xx, yy, width, n_comps,
                          weights, intercept, center, scale,
                          min_obs, complete_obs, na_restore,
                          online);
    
    arma::mat arma_coef_z = result_z[0];
    arma::mat arma_rsq_z = result_z[1];
    // arma::mat arma_se_z = result_z[2];
    
    // create and return a matrix or xts object for coefficients
    NumericVector coef(wrap(arma_coef_z));
    coef.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
    List dimnames_x = xx.attr("dimnames");
    coef.attr("dimnames") = dimnames_lm_x(dimnames_x, n_cols_x, intercept);
    coef.attr("index") = xx.attr("index");
    coef.attr(".indexCLASS") = xx.attr(".indexCLASS");
    coef.attr(".indexTZ") = xx.attr(".indexTZ");
    coef.attr("tclass") = xx.attr("tclass");
    coef.attr("tzone") = xx.attr("tzone");
    coef.attr("class") = xx.attr("class");
    
    // create and return a matrix or xts object for r-squareds
    NumericVector rsq(wrap(arma_rsq_z));
    rsq.attr("dim") = IntegerVector::create(n_rows_xy, 1);
    if (dimnames_x.size() > 1) {
      rsq.attr("dimnames") = List::create(dimnames_x[0], "R-squared");
    } else {
      rsq.attr("dimnames") = List::create(R_NilValue, "R-squared");
    }
    rsq.attr("index") = xx.attr("index");
    rsq.attr(".indexCLASS") = xx.attr(".indexCLASS");
    rsq.attr(".indexTZ") = xx.attr(".indexTZ");
    rsq.attr("tclass") = xx.attr("tclass");
    rsq.attr("tzone") = xx.attr("tzone");
    rsq.attr("class") = xx.attr("class");
    
    // // create and return a matrix or xts object for standard errors
    // NumericVector se(wrap(arma_se_z));
    // se.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
    // se.attr("dimnames") = coef.attr("dimnames");
    // se.attr("index") = xx.attr("index");
    // se.attr(".indexCLASS") = xx.attr(".indexCLASS");
    // se.attr(".indexTZ") = xx.attr(".indexTZ");
    // se.attr("tclass") = xx.attr("tclass");
    // se.attr("tzone") = xx.attr("tzone");
    // se.attr("class") = xx.attr("class");
    
    // create and return a list
    result = List::create(Named("coefficients") = coef,
                          Named("r.squared") = rsq);
    // Named("std.error") = se);
    
    return result;
    
  } else {
    
    stop("'x' is not a matrix");
    
    // return 0;
    
  }
  
}