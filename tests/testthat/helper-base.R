rollapplyr_cube <- function(f, x, width) {
  
  n_rows <- nrow(x)
  n_cols <- ncol(x)
  result <- array(as.numeric(NA), c(n_cols, n_cols, n_rows))
  
  for (i in width:n_rows) {
    
    temp <- f(x[max(1, i - width + 1):i, , drop = FALSE])
    
    result[ , , i] <- temp
    
  }
  
  return(result)
  
}

rollapplyr_eigen <- function(x, width) {
  
  n_rows <- nrow(x)
  n_cols <- ncol(x)
  evals <- matrix(as.numeric(NA), n_rows, n_cols)
  evecs <- array(as.numeric(NA), c(n_cols, n_cols, n_rows))
  result <- list()
  
  if (zoo::is.zoo(x)) {
    x_attr <- attributes(x)
  }
  
  cov <- rollapplyr_cube(cov, x, width)
  
  for (i in 1:n_rows) {
    
    if (anyNA(cov[ , , i])) {
      
      evals[i , ] <- NA
      evecs[ , , i] <- NA
      
    } else {
      
      eigen <- eigen(cov[ , , i])
      evals[i, ] <- eigen$values
      evecs[ , , i] <- eigen$vectors
      
    }
    
  }
  
  if (exists("x_attr")) {
    
    attributes(evals) <- x_attr
    dimnames(evecs) <- list(x_attr$dimnames[[2]], NULL)
    
  }
  
  result <- list("values" = evals,
                 "vectors" = evecs)
  
  return(result)
  
}