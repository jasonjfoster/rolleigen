test_that("input checks", {
  
  # test data
  width <- test_width[2]
  
  for (ax in 1:length(test_check_x)) {
    
    # check_rows_equal
    expect_error(roll_pcr(test_check_x[[ax]], test_check_y[-1], width),
                 "number of rows in 'x' must equal the number of rows in 'y'")
    
    # check_pos_int: 'width'
    expect_error(roll_eigen(test_check_x[[ax]], width = 0),
                 "value of 'width' must be greater than or equal to one")
    expect_error(roll_pcr(test_check_x[[ax]], test_check_y, width = 0),
                 "value of 'width' must be greater than or equal to one")
    
    # check_pos_int: 'min_obs'
    expect_error(roll_eigen(test_check_x[[ax]], width, min_obs = 0),
                 "value of 'min_obs' must be greater than or equal to one")
    expect_error(roll_pcr(test_check_x[[ax]], test_check_y, width, min_obs = 0),
                 "value of 'min_obs' must be greater than or equal to one")
    
    # check_pos_int: 'n_comps'
    expect_error(roll_pcr(test_check_x[[ax]], test_check_y, width, n_comps = 0),
                 "value of 'n_comps' must be greater than or equal to one")
    
    # check_weights
    expect_error(roll_eigen(test_check_x[[ax]], width, rep(1, width - 1)),
                 "length of 'weights' must be greater than or equal to the number of rows in 'x' or 'width'")
    expect_error(roll_pcr(test_check_x[[ax]], test_check_y, width, weights = rep(1, width - 1)),
                 "length of 'weights' must be greater than or equal to the number of rows in 'x' (and 'y') or 'width'",
                 fixed = TRUE)
    
    # check_comps: 'n_comps' (maximum only, since 'n_comps' is a scalar)
    expect_error(roll_pcr(test_check_x[[ax]], test_check_y, width, NCOL(test_check_x[[ax]]) + 1),
                 "maximum value of 'n_comps' must be less than or equal to number of columns in 'x'")
    
  }
  
})

test_that("input warnings", {
  
  # test data
  width <- test_width[2]
  
  for (ax in 1:length(test_check_x)) {
    
    # check_lambda
    expect_warning(roll_eigen(test_check_x[[ax]], width, 1:width),
                   "'online' is only supported for equal or exponential decay 'weights'")
    expect_silent(roll_eigen(test_check_x[[ax]], width, 1:width, online = FALSE))
    expect_warning(roll_pcr(test_check_x[[ax]], test_check_y, width, weights = 1:width),
                   "'online' is only supported for equal or exponential decay 'weights'")
    expect_silent(roll_pcr(test_check_x[[ax]], test_check_y, width, weights = 1:width, online = FALSE))
    
    # 'intercept = FALSE'
    expect_warning(roll_pcr(test_check_x[[ax]], test_check_y, width, intercept = FALSE),
                   "'intercept = FALSE' is not supported")
    
    # 'center = FALSE'
    expect_warning(roll_pcr(test_check_x[[ax]], test_check_y, width, center = FALSE),
                   "'center = FALSE' is not supported")
    
    # 'scale = TRUE'
    expect_warning(roll_pcr(test_check_x[[ax]], test_check_y, width, scale = TRUE),
                   "'scale = TRUE' is not supported")
    
    # 'complete_obs = FALSE' (roll_pcr only)
    expect_warning(roll_pcr(test_check_x[[ax]], test_check_y, width, complete_obs = FALSE),
                   "'complete_obs = FALSE' is not supported")
    expect_silent(roll_eigen(test_check_x[[ax]], width, complete_obs = FALSE))
    
  }
  
})
