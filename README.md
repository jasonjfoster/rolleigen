# rolleigen

## Overview

`rolleigen` is a package that provides analytical computation of rolling and expanding eigenvalues and eigenvectors for time-series data.

## Installation

Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("jjf234/rollpca") # roll (>= 1.1.7)
```

## Usage

Load the package and supply a dataset:

``` r
library(rolleigen)

n <- 15
m <- 3
x <- matrix(rnorm(n * m), nrow = n, ncol = m)
y <- rnorm(n)
weights <- 0.9 ^ (n:1)
```
Then, to compute rolling eigenvalues and eigenvectors, use the `roll_eigen` function:

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

Or use the `roll_pcr` function to compute rolling and expanding principal component regressions:

``` r
# rolling regressions with complete windows
roll_pcr(x, y, width = 5, n_comps = 1)

# rolling regressions with partial windows
roll_pcr(x, y, width = 5, n_comps = 1, min_obs = 1)

# expanding regressions with partial windows
roll_pcr(x, y, width = n, n_comps = 1, min_obs = 1)

# expanding regressions with partial windows and weights
roll_pcr(x, y, width = n, n_comps = 1, min_obs = 1, weights = weights)
```

Note that handling of missing values is supported as well (see the `min_obs`, `complete_obs`, and `na_restore` arguments).