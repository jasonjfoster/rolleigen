test_that("equivalent to base::eigen", {
  
  # skip("long-running test")
  
  if (!requireNamespace("zoo", quietly = TRUE)) {
    skip("zoo package required for this test")
  }
  
  # test data
  test_zoo <- c(test_ls, list("random vector with 0's" = test_ls[[1]][ , 1]))
  
  for (a in 1:(length(test_zoo))) {
    for (b in 1:length(test_width)) {
      
      width <- test_width[b]     
      test_weights <- list(rep(1, width))
      
      for (c in 1:length(test_center)) {
        for (d in 1:length(test_scale)) {
          
          test_roll <- roll_eigen(test_zoo[[a]], width,
                                  test_weights[[1]], test_center[c],
                                  test_scale[d], test_order[2],
                                  test_min_obs[1])
          test_rollapplyr <- rollapplyr_eigen(test_zoo[[a]], width,
                                              test_center[c], test_scale[d])
          
          expect_equal(test_roll$values, test_rollapplyr$values,
                       check.attributes = FALSE)
          
          expect_equal(abs(test_roll$vectors[ , 1, ]), abs(test_rollapplyr$vectors[ , 1, ]),
                       check.attributes = FALSE)
          
        }
      }
      
    }
  }
  
})