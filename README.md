# rolleigen

[![](https://github.com/jasonjfoster/rolleigen/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/jasonjfoster/rolleigen/actions/workflows/check-standard.yaml)
[![](https://codecov.io/gh/jasonjfoster/rolleigen/graph/badge.svg)](https://app.codecov.io/github/jasonjfoster/rolleigen)

## Overview

'rolleigen' provides fast and efficient computation of rolling and expanding eigenanalysis for time-series data.

The 'rolleigen' package decomposes the covariance matrix of the explanatory variables into eigenvalues and eigenvectors to perform principal component analysis (Pearson, 1901, <doi:10.1080/14786440109462720>; Hotelling, 1933, <doi:10.1037/h0071325>) and principal component regression (Massy, 1965, <doi:10.1080/01621459.1965.10480787>) over rolling and expanding windows. For each window, the eigenvalues and eigenvectors are computed from the covariance matrix and, optionally, ordered from largest to smallest to summarize the directions of greatest variation in the data. A subset of leading components is then used to fit a regression that mitigates collinearity in the explanatory variables. Use cases include:

* **Dimensionality reduction**: extracting a small number of components that capture the majority of variation across many correlated variables
* **Factor extraction**: identifying latent factors that drive co-movement in a set of explanatory variables through time
* **Principal component regression**: fitting a regression on leading components to mitigate collinearity in the explanatory variables

The package supports rolling and expanding windows, weights, and handling of missing values via the min_obs, complete_obs, and na_restore arguments. The implementation uses the online and offline algorithms from the 'roll' package to compute rolling and expanding cross-products efficiently, with parallelism across columns and windows provided by 'RcppParallel'.

## Installation

Install the development version from GitHub:

```r
# install.packages("pak")
pak::pak("jasonjfoster/rolleigen")
```

## Usage

Load the package and supply a dataset:

```r
library(rolleigen) # roll (>= 1.1.7)

n <- 15
m <- 3
x <- matrix(rnorm(n * m), nrow = n, ncol = m)
y <- rnorm(n)
weights <- 0.9 ^ (n:1)
```

Then, to compute rolling and expanding eigenvalues and eigenvectors, use the `roll_eigen()` function:

```r
# rolling eigenvalues and eigenvectors with complete windows
roll_eigen(x, width = 5)

# rolling eigenvalues and eigenvectors with partial windows
roll_eigen(x, width = 5, min_obs = 1)

# expanding eigenvalues and eigenvectors with partial windows
roll_eigen(x, width = n, min_obs = 1)

# expanding eigenvalues and eigenvectors with partial windows and weights
roll_eigen(x, width = n, min_obs = 1, weights = weights)
```

Or use the `roll_pcr()` function to compute rolling and expanding principal component regressions:

```r
# rolling regressions with complete windows
roll_pcr(x, y, width = 5, n_comps = 1)

# rolling regressions with partial windows
roll_pcr(x, y, width = 5, n_comps = 1, min_obs = 1)

# expanding regressions with partial windows
roll_pcr(x, y, width = n, n_comps = 1, min_obs = 1)

# expanding regressions with partial windows and weights
roll_pcr(x, y, width = n, n_comps = 1, min_obs = 1, weights = weights)
```

Handling of missing values is also supported (see the `min_obs`, `complete_obs`, and `na_restore` arguments).

## References

Hotelling, H. (1933). "Analysis of a Complex of Statistical Variables Into Principal Components." *Journal of Educational Psychology* 24 (6): 417-441. <doi:10.1037/h0071325>

Massy, W.F. (1965). "Principal Components Regression in Exploratory Statistical Research." *Journal of the American Statistical Association* 60 (309): 234-256. <doi:10.1080/01621459.1965.10480787>

Pearson, K. (1901). "On Lines and Planes of Closest Fit to Systems of Points in Space." *Philosophical Magazine* 2 (11): 559-572. <doi:10.1080/14786440109462720>
