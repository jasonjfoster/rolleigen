# rollpca

[![](https://github.com/jjf234/rollpca/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/jjf234/rollpca/actions/workflows/check-standard.yaml)

## Overview

`rollpca` is a package that provides analytical computation of rolling and expanding principal component analysis for time-series data.

## Installation

Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("jjf234/rollpca") # roll (>= 1.1.7)
```

## Usage

Load the package and supply a dataset:

``` r
library(roll)
library(rollpca)

n <- 15
m <- 3
x <- matrix(rnorm(n * m), nrow = n, ncol = m)
```
Then, to compute rolling eigenvalues and eigenvectors, use the `roll_eigen` function:

```r
# rolling eigenvalues and eigenvectors with complete windows
roll_eigen(x, width = 5)
```

Note that handling of missing values is supported as well (see the `min_obs`, `complete_obs`, and `na_restore` arguments).