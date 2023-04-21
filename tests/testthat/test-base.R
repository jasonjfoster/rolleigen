test_that("equivalent to base::eigen", {
  
  # skip("long-running test")
  
  for (a in 1:(length(test_ls))) {
    for (b in 1:length(test_width)) {
      
      width <- test_width[b]     
      test_weights <- list(rep(1, width))
      
      test_roll <- roll_eigen(test_ls[[a]], width,
                              min_obs = test_min_obs[1])
      test_rollapplyr <- rollapplyr_eigen(test_ls[[a]], width)
      
      expect_equal(test_roll$values, test_rollapplyr$values,
                   check.attributes = FALSE)
      
      expect_equal(abs(test_roll$vectors)[ , 1, ], abs(test_rollapplyr$vectors)[ , 1, ],
                   check.attributes = FALSE)
      
    }
  }
  
})