test_that("equivalent to pls::pcr", {
  
  packages <- c("zoo", "pls")
  
  status <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  
  if (!all(status)) {
    
    skip(paste0("'", paste0(names(status[!status]), collapse = "', '"),
                "' package(s) required for this test"))
    
  }
  
  # test data
  test_data_x <- c(lapply(test_ls[1:2], function(x){x[ , 1:3]}),
                   list("random vector with 0's" = test_ls[[1]][ , 1]))
  test_data_y <- c(lapply(test_ls[1:2], function(x){x[ , 4, drop = FALSE]}), # univariate 'y' for pls::pcr
                   list("random vector with 0's" = test_ls[[1]][ , 4]))
  
  for (ax in 1:length(test_data_x)) {
    for (b in 1:length(test_width)) {
      
      width <- test_width[b]
      test_weights <- list(rep(1, width))
      
      for (f in 1:length(test_weights)) {
        for (g in 1:length(test_comps)) {
          
          # "maximum value of 'n_comps' must be less than or equal to number of columns in 'x'"
          n_comps <- min(test_comps[g], NCOL(test_data_x[[ax]]))
          
          for (ay in 1:length(test_data_y)) {
            
            expect_equal(roll_pcr(test_data_x[[ax]], test_data_y[[ay]],
                                  width, n_comps,
                                  test_weights[[f]], test_intercept[1],
                                  test_center[1], test_scale[2],
                                  test_min_obs[1]),
                         rollapplyr_pcr(test_data_x[[ax]], test_data_y[[ay]],
                                        width, n_comps,
                                        test_intercept[1], test_center[1],
                                        test_scale[2]))
            
          }
          
        }
      }
      
    }
  }
  
})
