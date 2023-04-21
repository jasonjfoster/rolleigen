rollapplyr_cube <- function(f, x, width) {
  
  n_rows_x <- nrow(x)
  n_cols_x <- ncol(x)
  result <- array(as.numeric(NA), c(n_cols_x, n_cols_x, n_rows_x))
  
  for (i in 1:n_rows_x) {
    result[ , , i] <- f(x[max(1, i - width + 1):i, , drop = FALSE])
  }
  
  return(result)
  
}

rollapplyr_eigen <- function(x, width) {
  
  n_rows_x <- nrow(x)
  n_cols_x <- ncol(x)
  result <- list("values" = matrix(as.numeric(NA), n_rows_x, n_cols_x),
                 "vectors" = array(as.numeric(NA), c(n_cols_x, n_cols_x, n_rows_x)))
  
    x_attr <- attributes(x)
  
  cov <- rollapplyr_cube(cov, x, width)
  
  for (i in 1:n_rows_x) {
    
    if (anyNA(cov[ , , i])) {
      
      result[["values"]][i , ] <- NA
      result[["vectors"]][ , , i] <- NA
      
    } else {
      
      eigen <- eigen(cov[ , , i])
      result[["values"]][i, ] <- eigen$values
      result[["vectors"]][ , , i] <- eigen$vectors
      
    }
    
  }
  
  if (exists("x_attr")) {
    
    attributes(result[["values"]]) <- x_attr
    dimnames(result[["vectors"]]) <- list(x_attr$dimnames[[2]], NULL)
    
  }
  
  return(result)
  
}