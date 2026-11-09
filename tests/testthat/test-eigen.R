test_that("equivalent to base::eigen", {
  
  packages <- c("zoo")
  
  status <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  
  if (!all(status)) {
    
    skip(paste0("'", paste0(names(status[!status]), collapse = "', '"),
                "' package(s) required for this test"))
    
  }
  
  # test data
  test_zoo <- c(test_ls[1:2], list("random vector with 0's" = test_ls[[1]][ , 1]))
  
  for (ax in 1:length(test_zoo)) {
    for (b in 1:length(test_width)) {
      
      width <- test_width[b]     
      # test_weights <- list(rep(1, width))
      test_weights <- list(lambda ^ (width:1))
      
      for (c in 1:length(test_center)) {
        for (d in 1:length(test_scale)) {
          for (e in 1:length(test_order)) {
            
            # ordering assumes full-rank windows
            if (test_order[e]) {
              min_obs <- test_min_obs[2]
            } else {
              min_obs <- test_min_obs[1]
            }
            
            if (width >= min_obs) {
              
              test_roll <- roll_eigen(test_zoo[[ax]], width,
                                      test_weights[[1]], test_center[c],
                                      test_scale[d], test_order[e],
                                      min_obs)
              test_rollapplyr <- rollapplyr_eigen(test_zoo[[ax]], width,
                                                  test_weights[[1]], test_center[c],
                                                  test_scale[d], test_order[e],
                                                  min_obs)
              
              expect_equal(test_roll$values, test_rollapplyr$values,
                           check.attributes = FALSE)
              
              if (test_order[e]) {
                
                test_signs <- sign(test_roll$vectors[1, , min_obs] *
                                     test_rollapplyr$vectors[1, , min_obs])
                
                expect_equal(sweep(test_roll$vectors[ , , min_obs:n_obs, drop = FALSE], 2, test_signs, "*"),
                             test_rollapplyr$vectors[ , , min_obs:n_obs, drop = FALSE],
                             check.attributes = FALSE)
                
              } else {
                
                expect_equal(abs(test_roll$vectors[ , 1, ]), abs(test_rollapplyr$vectors[ , 1, ]),
                             check.attributes = FALSE)
                
              }
              
            }
            
          }
        }
      }
      
    }
  }
  
})
