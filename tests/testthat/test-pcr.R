test_that("equivalent to pls::pcr", {
  
  # skip("long-running test")
  
  packages <- c("zoo", "pls")
  
  status <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  
  if (!all(status)) {
    
    skip(paste0(paste0("'", paste0(names(status[!status]), collapse = ", '", sep = "'")),
                "' package(s) required for this test"))
    
  }
  
  # test data
  test_data_x <- test_ls[[1]][ , 1:3]
  test_data_y <- test_ls[[1]][ , 4, drop = FALSE] # univariate 'y' for pls::pcr
  
  # for (ax in 1:(length(test_data_x))) {
  for (b in 1:length(test_width)) {
    
    width <- test_width[b]     
    test_weights <- list(rep(1, width))
    
    for (c in 1:length(test_comps)) {
      # for (ay in 1:(length(test_data_y))) {
      
      expect_equal(roll_pcr(test_data_x, test_data_y,
                            width, test_comps[c],
                            test_weights[[1]], test_intercept[1],
                            test_center[1], test_scale[2],
                            test_min_obs[1]),
                   rollapplyr_pcr(test_data_x, test_data_y,
                                  width, test_comps[c],
                                  test_intercept[1], test_center[1],
                                  test_scale[2]),
                   check.attributes = FALSE)
      
      # }
    }
    
  }
  # }
  
})