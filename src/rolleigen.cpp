#include "rolleigen.h"

List dimnames_eigen(const int& n_cols) {
  
  CharacterVector result(n_cols);
  
  for (int i = 0; i < n_cols; i++) {
    
    result[i] = "pc";
    result[i] += i + 1;
    
  }
  
  return List::create(R_NilValue, result);
  
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
  rolleigen::check_pos_int(width, "width");
  
  // default 'weights' argument is equal-weighted,
  // otherwise check argument for errors
  rolleigen::check_weights(n_rows_x, width, weights, "'x'");
  bool status = rolleigen::check_lambda(weights, n_rows_x, width, online);
  
  // default 'min_obs' argument is 'width',
  // otherwise check argument for errors
  rolleigen::check_pos_int(min_obs, "min_obs");
  
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
List roll_eigen(const SEXP& x, const int& width,
                const arma::vec& weights, const bool& center,
                const bool& scale, const bool& order,
                const int& min_obs, const bool& complete_obs,
                const bool& na_restore, const bool& online) {
  
  if (Rf_isMatrix(x)) {
    
    NumericMatrix xx(x);
    int n_rows = xx.nrow();
    int n_cols = xx.ncol();
    
    arma::cube arma_cov = roll_cov_z(xx, width, weights, center, scale, min_obs, complete_obs,
                                     na_restore, online);
    arma::mat arma_eigen_values(n_rows, n_cols);
    arma::cube arma_eigen_vectors(n_cols, n_cols, n_rows);
    
    // compute rolling eigenvalues and eigenvectors
    rolleigen::RollEigenSlices roll_eigen_slices(arma_cov, n_rows, n_cols,
                                                 arma_eigen_values, arma_eigen_vectors);
    parallelFor(0, n_rows, roll_eigen_slices);
    
    if (order) {
      
      // order rolling eigenvalues and eigenvectors
      rolleigen::RollOrderSlices roll_order_slices(n_rows, n_cols,
                                                   arma_eigen_values, arma_eigen_vectors);
      roll_order_slices(0, n_rows);
      
    }
    
    // create and return a matrix or xts object for eigenvalues
    NumericMatrix eigen_values(wrap(arma_eigen_values));
    List dimnames = xx.attr("dimnames");
    rolleigen::xts_attr(eigen_values, xx, dimnames_eigen(n_cols));
    
    // create and return a cube for eigenvectors
    NumericVector eigen_vectors(wrap(arma_eigen_vectors));
    rolleigen::cube_attr(eigen_vectors, n_cols, n_cols, n_rows, dimnames, List());
    
    // create and return a list
    List result = List::create(Named("values") = eigen_values,
                               Named("vectors") = eigen_vectors);
    
    return result;
    
  } else {
    
    NumericVector xx(x);
    NumericMatrix xxx(xx.size(), 1, xx.begin());
    int n_rows = xxx.nrow();
    int n_cols = xxx.ncol();
    
    arma::cube arma_cov = roll_cov_z(xxx, width, weights, center, scale, min_obs, complete_obs,
                                     na_restore, online);
    arma::mat arma_eigen_values(n_rows, n_cols);
    arma::cube arma_eigen_vectors(n_cols, n_cols, n_rows);
    
    // compute rolling eigenvalues and eigenvectors
    rolleigen::RollEigenSlices roll_eigen_slices(arma_cov, n_rows, n_cols,
                                                 arma_eigen_values, arma_eigen_vectors);
    parallelFor(0, n_rows, roll_eigen_slices);
    
    if (order) {
      
      // order rolling eigenvalues and eigenvectors
      rolleigen::RollOrderSlices roll_order_slices(n_rows, n_cols,
                                                   arma_eigen_values, arma_eigen_vectors);
      roll_order_slices(0, n_rows);
      
    }
    
    // create and return a matrix or xts object for eigenvalues
    NumericMatrix eigen_values(wrap(arma_eigen_values));
    List dimnames = xx.attr("dimnames");
    rolleigen::xts_attr(eigen_values, xx, dimnames);
    
    // create and return a cube for eigenvectors
    NumericVector eigen_vectors(wrap(arma_eigen_vectors));
    rolleigen::cube_attr(eigen_vectors, n_cols, n_cols, n_rows, dimnames, List());
    
    // create and return a list
    List result = List::create(Named("values") = eigen_values,
                               Named("vectors") = eigen_vectors);
    
    return result;
    
  }
  
}

List roll_pcr_z(const SEXP& x, const NumericVector& y,
                const int& width, const int& n_comps,
                const arma::vec& weights, const bool& intercept,
                const bool& center, const bool& scale,
                const int& min_obs, const bool& complete_obs,
                const bool& na_restore, const bool& online) {
  
  if (Rf_isMatrix(x)) {
    
    NumericMatrix xx(x);
    int n = weights.size();
    int n_rows_xy = xx.nrow();
    int n_cols_x = xx.ncol() + 1;
    arma::uvec arma_comps = seq(n_comps) + 1;
    arma::uvec arma_any_na(n_rows_xy);
    arma::vec arma_n_obs(n_rows_xy);
    arma::vec arma_sum_w(n_rows_xy);
    arma::mat arma_mean(n_rows_xy, n_cols_x);
    arma::cube arma_cov(n_cols_x, n_cols_x, n_rows_xy);
    arma::mat arma_eigen_values(n_rows_xy, n_cols_x - 1);
    arma::cube arma_eigen_vectors(n_cols_x - 1, n_cols_x - 1, n_rows_xy);
    arma::mat arma_coef(n_rows_xy, n_cols_x);
    arma::vec arma_rsq(n_rows_xy);
    
    // check 'x' and 'y' arguments for errors
    rolleigen::check_rows_equal(n_rows_xy, y.size(), "x", "y");
    
    // check 'width' argument for errors
    rolleigen::check_pos_int(width, "width");
    
    // default 'n_comps' argument is all components
    rolleigen::check_comps(arma_comps, n_cols_x - 1, "n_comps", "x");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    rolleigen::check_weights(n_rows_xy, width, weights, "'x' (and 'y')");
    bool status = rolleigen::check_lambda(weights, n_rows_xy, width, online);
    
    // check 'intercept' argument for errors
    if (!intercept) {
      warning("'intercept = FALSE' is not supported");
    }
    
    // check 'center' argument for errors
    if (!center) {
      warning("'center = FALSE' is not supported");
    }
    
    // check 'scale' argument for errors
    if (scale) {
      warning("'scale = TRUE' is not supported");
    }
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    rolleigen::check_pos_int(min_obs, "min_obs");
    
    // cbind x and y variables
    NumericMatrix data(n_rows_xy, n_cols_x);
    std::copy(xx.begin(), xx.end(), data.begin());
    std::copy(y.begin(), y.end(), data.begin() + n_rows_xy * (n_cols_x - 1));
    
    // default 'complete_obs' argument is 'true',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_x(data);
    } else {
      
      warning("'complete_obs = FALSE' is not supported");
      arma_any_na = any_na_x(data);
      
    }
    
    // compute rolling crossproducts
    if (status && online) {
      
      roll::RollCrossProdOnlineMatXX roll_cov_online(data, n, n_rows_xy, n_cols_x, width,
                                                     weights, true, false, min_obs,
                                                     arma_any_na, na_restore,
                                                     arma_n_obs, arma_sum_w, arma_mean,
                                                     arma_cov);
      parallelFor(0, n_cols_x, roll_cov_online);
      
    } else {
      
      roll::RollCrossProdOfflineMatXX roll_cov_offline(data, n, n_rows_xy, n_cols_x, width,
                                                       weights, true, false, min_obs,
                                                       arma_any_na, na_restore,
                                                       arma_n_obs, arma_sum_w, arma_mean,
                                                       arma_cov);
      parallelFor(0, n_rows_xy * n_cols_x * (n_cols_x + 1) / 2, roll_cov_offline);
      
    }
    
    // compute rolling eigenvalues and eigenvectors
    rolleigen::RollEigenSlices roll_eigen_slices(arma_cov, n_rows_xy, n_cols_x - 1,
                                                 arma_eigen_values, arma_eigen_vectors);
    parallelFor(0, n_rows_xy, roll_eigen_slices);
    
    // compute rolling principal component regressions
    arma::uvec arma_cols = seq(n_cols_x - 1);
    rolleigen::RollPcrInterceptTRUE roll_pcr_slices(arma_cov, n_rows_xy, n_cols_x, arma_cols, arma_comps,
                                                    arma_n_obs, arma_mean, arma_eigen_vectors,
                                                    arma_coef, arma_rsq);
    parallelFor(0, n_rows_xy, roll_pcr_slices);
    
    // create and return a list
    List result = List::create(Named("coefficients") = arma_coef,
                               Named("r.squared") = arma_rsq);
    
    return result;
    
  } else {
    
    NumericVector xx(x);
    int n = weights.size();
    int n_rows_xy = xx.size();
    int n_cols_x = 1 + 1;
    arma::uvec arma_comps = seq(n_comps) + 1;
    arma::uvec arma_any_na(n_rows_xy);
    arma::vec arma_n_obs(n_rows_xy);
    arma::vec arma_sum_w(n_rows_xy);
    arma::mat arma_mean(n_rows_xy, n_cols_x);
    arma::cube arma_cov(n_cols_x, n_cols_x, n_rows_xy);
    arma::mat arma_eigen_values(n_rows_xy, n_cols_x - 1);
    arma::cube arma_eigen_vectors(n_cols_x - 1, n_cols_x - 1, n_rows_xy);
    arma::mat arma_coef(n_rows_xy, n_cols_x);
    arma::vec arma_rsq(n_rows_xy);
    
    // check 'x' and 'y' arguments for errors
    rolleigen::check_rows_equal(n_rows_xy, y.size(), "x", "y");
    
    // check 'width' argument for errors
    rolleigen::check_pos_int(width, "width");
    
    // default 'n_comps' argument is all components
    rolleigen::check_comps(arma_comps, n_cols_x - 1, "n_comps", "x");
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    rolleigen::check_weights(n_rows_xy, width, weights, "'x' (and 'y')");
    bool status = rolleigen::check_lambda(weights, n_rows_xy, width, online);
    
    // check 'intercept' argument for errors
    if (!intercept) {
      warning("'intercept = FALSE' is not supported");
    }
    
    // check 'center' argument for errors
    if (!center) {
      warning("'center = FALSE' is not supported");
    }
    
    // check 'scale' argument for errors
    if (scale) {
      warning("'scale = TRUE' is not supported");
    }
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    rolleigen::check_pos_int(min_obs, "min_obs");
    
    // cbind x and y variables
    NumericMatrix data(n_rows_xy, n_cols_x);
    std::copy(xx.begin(), xx.end(), data.begin());
    std::copy(y.begin(), y.end(), data.begin() + n_rows_xy * (n_cols_x - 1));
    
    // default 'complete_obs' argument is 'true',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_x(data);
    } else {
      
      warning("'complete_obs = FALSE' is not supported");
      arma_any_na = any_na_x(data);
      
    }
    
    // compute rolling crossproducts
    if (status && online) {
      
      roll::RollCrossProdOnlineMatXX roll_cov_online(data, n, n_rows_xy, n_cols_x, width,
                                                     weights, true, false, min_obs,
                                                     arma_any_na, na_restore,
                                                     arma_n_obs, arma_sum_w, arma_mean,
                                                     arma_cov);
      parallelFor(0, n_cols_x, roll_cov_online);
      
    } else {
      
      roll::RollCrossProdOfflineMatXX roll_cov_offline(data, n, n_rows_xy, n_cols_x, width,
                                                       weights, true, false, min_obs,
                                                       arma_any_na, na_restore,
                                                       arma_n_obs, arma_sum_w, arma_mean,
                                                       arma_cov);
      parallelFor(0, n_rows_xy * n_cols_x * (n_cols_x + 1) / 2, roll_cov_offline);
      
    }
    
    // compute rolling eigenvalues and eigenvectors
    rolleigen::RollEigenSlices roll_eigen_slices(arma_cov, n_rows_xy, n_cols_x - 1,
                                                 arma_eigen_values, arma_eigen_vectors);
    parallelFor(0, n_rows_xy, roll_eigen_slices);
    
    // compute rolling principal component regressions
    arma::uvec arma_cols = seq(n_cols_x - 1);
    rolleigen::RollPcrInterceptTRUE roll_pcr_slices(arma_cov, n_rows_xy, n_cols_x, arma_cols, arma_comps,
                                                    arma_n_obs, arma_mean, arma_eigen_vectors,
                                                    arma_coef, arma_rsq);
    parallelFor(0, n_rows_xy, roll_pcr_slices);
    
    // create and return a list
    List result = List::create(Named("coefficients") = arma_coef,
                               Named("r.squared") = arma_rsq);
    
    return result;
    
  }
  
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
      rolleigen::xts_attr(coef, xx, dimnames_lm_x(x_dimnames, n_cols_x, intercept));

      // create and return a matrix or xts object for r-squareds
      NumericVector rsq(wrap(arma_rsq_z));
      rolleigen::rsq_attr(rsq, n_rows_xy, xx, x_dimnames);

      // // create and return a matrix or xts object for standard errors
      // NumericVector se(wrap(arma_se_z));
      // se.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
      // rolleigen::xts_attr(se, xx, coef.attr("dimnames"));

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
        rolleigen::xts_attr(coef, xx, dimnames_lm_x(dimnames_x, n_cols_x, intercept));

        // create and return a matrix or xts object for r-squareds
        NumericVector rsq(wrap(arma_rsq_z));
        rolleigen::rsq_attr(rsq, n_rows_xy, xx, dimnames_x);

        // // create and return a matrix or xts object for standard errors
        // NumericVector se(wrap(arma_se_z));
        // se.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
        // rolleigen::xts_attr(se, xx, coef.attr("dimnames"));

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
    rolleigen::xts_attr(coef, xx, dimnames_lm_x(dimnames_x, n_cols_x, intercept));
    
    // create and return a matrix or xts object for r-squareds
    NumericVector rsq(wrap(arma_rsq_z));
    rolleigen::rsq_attr(rsq, n_rows_xy, xx, dimnames_x);
    
    // // create and return a matrix or xts object for standard errors
    // NumericVector se(wrap(arma_se_z));
    // se.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
    // rolleigen::xts_attr(se, xx, coef.attr("dimnames"));
    
    // create and return a list
    result = List::create(Named("coefficients") = coef,
                          Named("r.squared") = rsq);
    // Named("std.error") = se);
    
    return result;
    
  } else if (Rf_isMatrix(y)) {
    
    NumericVector xx(x);
    NumericMatrix yy(y);
    // xx.attr("dim") = IntegerVector::create(xx.size(), 1);
    // NumericMatrix xxx(wrap(xx));
    NumericMatrix xxx(xx.size(), 1, xx.begin());
    
    int n_rows_xy = xxx.nrow();
    int n_cols_x = xxx.ncol();
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
      
      result_z = roll_pcr_z(xxx, yy(_, 0), width, n_comps,
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
      rolleigen::xts_attr(coef, yy, dimnames_lm_x(dimnames_x, n_cols_x, intercept));
      
      // create and return a matrix or xts object for r-squareds
      NumericVector rsq(wrap(arma_rsq_z));
      rolleigen::rsq_attr(rsq, n_rows_xy, yy, dimnames_x);
      
      // // create and return a matrix or xts object for standard errors
      // NumericVector se(wrap(arma_se_z));
      // se.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
      // rolleigen::xts_attr(se, yy, coef.attr("dimnames"));
      
      // create and return a list
      result = List::create(Named("coefficients") = coef,
                            Named("r.squared") = rsq);
      // Named("std.error") = se);
      
    } else {
      
      for (int z = 0; z < n_cols_y; z++) {
        
        result_z = roll_pcr_z(xxx, yy(_, z), width, n_comps,
                              weights, intercept, center, scale,
                              min_obs, complete_obs, na_restore,
                              online);
        
        arma::mat arma_coef_z = result_z[0];
        arma::mat arma_rsq_z = result_z[1];
        // arma::mat arma_se_z = result_z[2];
        
        // create and return a matrix or xts object for coefficients
        NumericVector coef(wrap(arma_coef_z));
        coef.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
        List dimnames_x = xxx.attr("dimnames");
        rolleigen::xts_attr(coef, yy, dimnames_lm_x(dimnames_x, n_cols_x, intercept));
        
        // create and return a matrix or xts object for r-squareds
        NumericVector rsq(wrap(arma_rsq_z));
        rolleigen::rsq_attr(rsq, n_rows_xy, yy, dimnames_x);
        
        // // create and return a matrix or xts object for standard errors
        // NumericVector se(wrap(arma_se_z));
        // se.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
        // rolleigen::xts_attr(se, yy, coef.attr("dimnames"));
        
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
    
  } else {
    
    NumericVector xx(x);
    NumericVector yy(y);
    
    int n_rows_xy = xx.size();
    int n_cols_x = 1;
    int n_cols_y = 1;
    List result_coef(n_cols_y);
    List result_rsq(n_cols_y);
    // List result_se(n_cols_y);
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
    rolleigen::xts_attr(coef, xx, dimnames_lm_x(dimnames_x, n_cols_x, intercept));
    
    // create and return a matrix or xts object for r-squareds
    NumericVector rsq(wrap(arma_rsq_z));
    rolleigen::rsq_attr(rsq, n_rows_xy, xx, dimnames_x);
    
    // // create and return a matrix or xts object for standard errors
    // NumericVector se(wrap(arma_se_z));
    // se.attr("dim") = IntegerVector::create(n_rows_xy, n_cols_x);
    // rolleigen::xts_attr(se, xx, coef.attr("dimnames"));
    
    // create and return a list
    result = List::create(Named("coefficients") = coef,
                          Named("r.squared") = rsq);
                          // Named("std.error") = se);
    
    return result;
    
  }
  
}