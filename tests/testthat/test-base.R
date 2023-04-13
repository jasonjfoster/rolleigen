test_that("equivalent to base::eigen", {
  
  # skip("long-running test")
  
  if (!requireNamespace("zoo", quietly = TRUE)) {
    skip("zoo package required for this test")
  }
  
  for (a in 1:(length(test_ls))) {
    for (b in 1:length(test_width)) {
      
      width <- test_width[b]     
      test_weights <- list(rep(1, width))
      
      test_roll <- roll_eigen(test_ls[[a]], width)
      test_rollapplyr <- rollapplyr_eigen(test_ls[[a]], width)
      
      expect_equal(test_roll$values, test_rollapplyr$values,
                   check.attributes = FALSE)
      
      expect_equal(abs(test_roll$vectors), abs(test_rollapplyr$vectors),
                   check.attributes = FALSE)
      
    }
  }
  
})