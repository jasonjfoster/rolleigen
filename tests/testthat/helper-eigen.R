cov.wt_cov <- function(x, center) {
  
  result <- cov.wt(x, center = center)$cov
  
  return(result)
  
}

cov.wt_cor <- function(x, center) {
  
  result <- cov.wt(x, cor = TRUE, center = center)$cor
  
  return(result)
  
}

rollapplyr_cube <- function(f, x, width, ...) {
  
  n_rows_x <- nrow(x)
  n_cols_x <- ncol(x)
  result <- array(as.numeric(NA), c(n_cols_x, n_cols_x, n_rows_x))
  
  for (i in 1:n_rows_x) {
    result[ , , i] <- f(x[max(1, i - width + 1):i, , drop = FALSE], ...)
  }
  
  # attr(result, "dim") <- c(n_cols_x, n_cols_x, n_rows_x)
  # 
  # x_dimnames <- dimnames(x)
  # attr(result, "dimnames") <- list(x_dimnames[[2]], x_dimnames[[2]], NULL)
  
  return(result)
  
}

rollapplyr_eigen <- function(x, width, center, scale, order, min_obs) {
  
  if (!is.matrix(x)) {
    
    temp_attr <- attributes(x)
    x <- as.matrix(zoo::coredata(x))
    attr(x, "dimnames") <- NULL
    attr(x, "index") <- temp_attr[["index"]]
    attr(x, "class") <- temp_attr[["class"]]
    
  }
  
  if (zoo::is.zoo(x)) {
    
    x_attr <- attributes(x)
    # x_attr[["dim"]] <- NULL
    x_attr[["dimnames"]] <- NULL
    
  }
  
  n_rows_x <- nrow(x)
  n_cols_x <- ncol(x)
  result <- list("values" = matrix(as.numeric(NA), n_rows_x, n_cols_x),
                 "vectors" = array(as.numeric(NA), c(n_cols_x, n_cols_x, n_rows_x)))
    
  if (scale) {
    cov <- rollapplyr_cube(cov.wt_cor, x, width, center = center)
  } else {
    cov <- rollapplyr_cube(cov.wt_cov, x, width, center = center)
  }
  
  for (i in 1:n_rows_x) {
    
    if ((min(i, width) < min_obs) || anyNA(cov[ , , i]) || any(is.infinite(cov[ , , i]))) {
      
      result[["values"]][i , ] <- NA
      result[["vectors"]][ , , i] <- NA
      
    } else {
      
      eigen <- eigen(cov[ , , i])
      result[["values"]][i, ] <- eigen$values
      result[["vectors"]][ , , i] <- eigen$vectors
      
    }
    
  }
  
  if (order) {
    result <- rollapplyr_order(result)
  }
  
  if (exists("x_attr")) {
    attributes(result[["values"]]) <- x_attr
  }
  
  x_dimnames <- dimnames(x)
  
  # if (length(x_dimnames) > 1) {
  #   attr(result[["values"]], "dimnames") <- list(x_dimnames[[1]], NULL)
  # }
  attr(result[["values"]], "dimnames") <- list(x_dimnames[[1]], NULL)
  attr(result[["vectors"]], "dimnames") <- list(x_dimnames[[2]], NULL)
  
  return(result)
  
}

rollapplyr_order <- function(result) {
  
  n_rows_x <- dim(result[["vectors"]])[3]
  n_cols_x <- dim(result[["vectors"]])[1]
  
  for (i in 1:n_rows_x) {
    
    if (i > 1) {
      
      eigen_values <- result[["values"]][i, ]
      eigen_vectors <- matrix(result[["vectors"]][ , , i], n_cols_x, n_cols_x)
      eigen_vectors0 <- matrix(result[["vectors"]][ , , i - 1], n_cols_x, n_cols_x)
      
      if (!anyNA(eigen_vectors) && !anyNA(eigen_vectors0)) {
        
        similarity <- crossprod(eigen_vectors, eigen_vectors0)
        order <- apply(abs(similarity), 2, which.max)
        
        eigen_vectors <- eigen_vectors[ , order, drop = FALSE]
        similarity <- similarity[order, , drop = FALSE]
        
        signs <- sign(diag(similarity))
        
        result[["values"]][i, ] <- eigen_values[order]
        result[["vectors"]][ , , i] <- sweep(eigen_vectors, 2, signs, "*")
        
      }
      
    }
    
  }
  
  return(result)
  
}
