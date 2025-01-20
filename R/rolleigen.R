##' Rolling Eigenvalues and Eigenvectors
##'
##' A function for computing the rolling and expanding eigenvalues and eigenvectors of time-series data.
##' 
##' @param x vector or matrix. Rows are observations and columns are variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param center logical. If \code{TRUE} then the weighted mean of each variable is used,
##' if \code{FALSE} then zero is used.
##' @param scale logical. If \code{TRUE} then the weighted standard deviation of each variable is used,
##' if \code{FALSE} then no scaling is done.
##' @param order logical. Change sign and order of the components.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then pairwise is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @return A list containing the following components:
##' \item{values}{An object of the same class and dimension as \code{x} with the rolling and expanding
##' eigenvalues.}
##' \item{vectors}{A cube with each slice the rolling and expanding eigenvectors.}
##' @examples
##' n <- 15
##' m <- 3
##' x <- matrix(rnorm(n * m), nrow = n, ncol = m)
##' weights <- 0.9 ^ (n:1)
##' 
##' # rolling eigenvalues and eigenvectors with complete windows
##' roll_eigen(x, width = 5)
##' 
##' # rolling eigenvalues and eigenvectors with partial windows
##' roll_eigen(x, width = 5, min_obs = 1)
##' 
##' # expanding eigenvalues and eigenvectors with partial windows
##' roll_eigen(x, width = n, min_obs = 1)
##' 
##' # expanding eigenvalues and eigenvectors with partial windows and weights
##' roll_eigen(x, width = n, min_obs = 1, weights = weights)
##' @export
roll_eigen <- function(x, width, weights = rep(1, width),
                       center = TRUE, scale = FALSE, order = TRUE,
                       min_obs = width, complete_obs = TRUE, na_restore = FALSE,
                       online = TRUE) {
  return(.Call(`_rolleigen_roll_eigen`,
               x,
               as.integer(width),
               as.numeric(weights),
               as.logical(center),
               as.logical(scale),
               as.logical(order),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}

##' Rolling Principal Component Regressions 
##'
##' A function for computing the rolling and expanding principal component regressions of time-series data.
##' 
##' @param x vector or matrix. Rows are observations and columns are the independent variables.
##' @param y vector or matrix. Rows are observations and columns are the dependent variables.
##' @param width integer. Window size.
##' @param n_comps integer. Number of principal components.
##' @param weights vector. Weights for each observation within a window.
##' @param intercept logical. Either \code{TRUE} to include or \code{FALSE} to remove the intercept.
##' @param center logical. If \code{TRUE} then the weighted mean of each variable is used,
##' if \code{FALSE} then zero is used.
##' @param scale logical. If \code{TRUE} then the weighted standard deviation of each variable is used,
##' if \code{FALSE} then no scaling is done.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then pairwise is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @return A list containing the following components:
##' \item{coefficients}{A list of objects with the rolling and expanding coefficients for each \code{y}.
##' An object is the same class and dimension (with an added column for the intercept) as \code{x}.}
##' \item{r.squared}{A list of objects with the rolling and expanding r-squareds for each \code{y}.
##' An object is the same class as \code{x}.}
##' @examples
##' n <- 15
##' x <- rnorm(n)
##' y <- rnorm(n)
##' weights <- 0.9 ^ (n:1)
##' 
##' # rolling regressions with complete windows
##' roll_pcr(x, y, width = 5, n_comps = 1)
##' 
##' # rolling regressions with partial windows
##' roll_pcr(x, y, width = 5, n_comps = 1, min_obs = 1)
##' 
##' # expanding regressions with partial windows
##' roll_pcr(x, y, width = n, n_comps = 1, min_obs = 1)
##' 
##' # expanding regressions with partial windows and weights
##' roll_pcr(x, y, width = n, n_comps = 1, min_obs = 1, weights = weights)
##' @export
roll_pcr <- function(x, y, width, n_comps = ncol(x), weights = rep(1, width),
                     intercept = TRUE, center = TRUE, scale = FALSE,
                     min_obs = width, complete_obs = TRUE, na_restore = FALSE,
                     online = TRUE) {
  return(.Call(`_rolleigen_roll_pcr`,
               x, y,
               as.integer(width),
               as.integer(n_comps),
               as.numeric(weights),
               as.logical(intercept),
               as.logical(center),
               as.logical(scale),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}