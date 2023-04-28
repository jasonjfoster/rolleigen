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
  # test_data_y <- test_ls[[1]][ , 4:5]
  test_data_y <- test_ls[[1]][ , 4, drop = FALSE]
  
  # for (a in 1:(length(test_data_x))) {
  for (b in 1:length(test_width)) {
    
    width <- test_width[b]     
    test_weights <- list(rep(1, width))
    
    for (c in 1:length(test_comps)) {
      
      expect_equal(roll_pcr(test_data_x, test_data_y,
                            width = width, n_comps = test_comps[c],
                            min_obs = test_min_obs[1]),
                   rollapplyr_pcr(test_data_x, test_data_y,
                                  width = width, n_comps = test_comps[c]),
                   check.attributes = FALSE)
      
    }
    
  }
  # }
  
})