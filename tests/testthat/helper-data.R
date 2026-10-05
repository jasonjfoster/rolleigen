set.seed(5640)
n_vars <- 5
n_obs <- 15 # 15 # 1500
n_size <- n_obs * n_vars
lambda <- 0.9 # 0.9 # 1 / 0.9
dates <- rev(seq(Sys.Date(), length.out = n_obs, by = "-1 day"))

# test arguments
test_width <- c(1, 5, 10, 15)
test_comps <- c(1, 2, 3)
test_intercept <- c(TRUE, FALSE)
test_center <- c(TRUE, FALSE)
test_scale <- c(TRUE, FALSE)
test_order <- c(TRUE, FALSE)
test_min_obs <- c(1, n_vars + 1)

# test data
test_ls <- list("random matrix with 0's" =
                  matrix(sample(c(0, rnorm(n_size - 1)), n_size, replace = TRUE,
                                prob = rep(1 / n_size, n_size)),
                         nrow = n_obs, ncol = n_vars),
                
                "random matrix with distinct eigenvalues" =
                  matrix(rnorm(n_size), nrow = n_obs, ncol = n_vars) %*%
                  diag(2 ^ (n_vars:1)))

if (requireNamespace("zoo", quietly = TRUE)) {
  test_ls <- lapply(test_ls, function(x) zoo::zoo(x, dates))
}

test_ls <- lapply(test_ls, setNames, paste0("x", rep(1:n_vars)))

# colnames(test_ls[[1]]) <- paste0("x", rep(1:n_vars))
