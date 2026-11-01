test_that("equal to online algorithm", {
  
  # test data
  test_roll_x <- c(lapply(test_ls, function(x){x[ , 1:3]}),
                   list("random vector with 0's and NA's" = test_ls[[3]][ , 1]))
  test_roll_y <- c(lapply(test_ls, function(x){x[ , 4:5]}),
                   list("random vector with 0's and NA's" = test_ls[[3]][ , 4]))
  
  if (requireNamespace("zoo", quietly = TRUE)) {
    names(test_roll_x[[4]]) <- zoo::index(test_roll_x[[1]])
    names(test_roll_y[[4]]) <- zoo::index(test_roll_y[[1]])
  }
  
  for (ax in 1:length(test_roll_x)) {
    for (b in 1:length(test_width)) {
      
      width <- test_width[b]
      # test_weights <- list(rep(1, width))
      test_weights <- list(lambda ^ ((2 * width):1))
      
      for (c in 1:length(test_min_obs)) {
        for (d in 1:length(test_complete_obs)) {
          for (e in 1:length(test_na_restore)) {
            for (f in 1:length(test_weights)) {
              for (g in 1:length(test_center)) {
                for (h in 1:length(test_scale)) {
                  for (j in 1:length(test_order)) {
                    
                    # ordering assumes full-rank windows
                    if (test_order[j]) {
                      min_obs <- test_min_obs[2]
                    } else {
                      min_obs <- test_min_obs[c]
                    }
                    
                    if (width >= min_obs) {
                      
                      result1 <- roll_eigen(test_roll_x[[ax]], width,
                                            test_weights[[f]], test_center[g],
                                            test_scale[h], test_order[j],
                                            min_obs, test_complete_obs[d],
                                            test_na_restore[e], test_online[1])
                      result2 <- roll_eigen(test_roll_x[[ax]], width,
                                            test_weights[[f]], test_center[g],
                                            test_scale[h], test_order[j],
                                            min_obs, test_complete_obs[d],
                                            test_na_restore[e], test_online[2])
                      
                      expect_equal(result1$values, result2$values)
                      
                      if (test_order[j]) {
                        
                        expect_equal(result1$vectors, result2$vectors)
                        
                      } else {
                        
                        # eigenvectors are unique (up to sign) only if eigenvalues are distinct
                        expect_equal(abs(result1$vectors[ , 1, ]), abs(result2$vectors[ , 1, ]))
                        
                      }
                      
                    }
                    
                  }
                }
              }
              
              for (g in 1:length(test_comps)) {
                
                # "maximum value of 'n_comps' must be less than or equal to number of columns in 'x'"
                n_comps <- min(test_comps[g], NCOL(test_roll_x[[ax]]))
                
                for (ay in 1:length(test_roll_y)) {
                  
                  # "'complete_obs = FALSE' is not supported"
                  expect_equal(suppressWarnings(roll_pcr(test_roll_x[[ax]], test_roll_y[[ay]],
                                                         width, n_comps,
                                                         test_weights[[f]], test_intercept[1],
                                                         test_center[1], test_scale[2],
                                                         test_min_obs[c], test_complete_obs[d],
                                                         test_na_restore[e], test_online[1])),
                               suppressWarnings(roll_pcr(test_roll_x[[ax]], test_roll_y[[ay]],
                                                         width, n_comps,
                                                         test_weights[[f]], test_intercept[1],
                                                         test_center[1], test_scale[2],
                                                         test_min_obs[c], test_complete_obs[d],
                                                         test_na_restore[e], test_online[2])))
                  
                }
              }
              
            }
          }
        }
      }
      
    }
  }
  
})
