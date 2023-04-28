dimnames_lm_x <- function(dimnames, n_cols_x, intercept) {
  
  if (intercept && (length(dimnames) > 1)) {
    
    result <- list(dimnames[[1]], c("(Intercept)", dimnames[[2]]))
    
  } else if (!intercept && (length(dimnames) > 1)) {
    
    result <- list(dimnames[[1]], dimnames[[2]])
    
  } else if (intercept) {
    
    result <- list(NULL, c("(Intercept)", paste0("x", rep(1:(n_cols_x - 1)))))
    
  } else {
    
    result <- list(NULL, paste0("x", rep(1:n_cols_x)))
    
  }
  
  return(result)
  
}

rollapplyr_pcr <- function(x, y, width, n_comps) {
  
  # "could not find function 'mvrValstats'"
  require(pls)
  
  if (is.matrix(y)) {
    
    n_rows_xy <- nrow(x)
    n_cols_x <- ncol(x) + 1
    
    result <- list("coefficients" = matrix(as.numeric(NA), n_rows_xy, n_cols_x),
                   "r.squared" = matrix(as.numeric(NA), n_rows_xy, 1))
    
    if (zoo::is.zoo(x)) {
      
      x_attr <- attributes(x)
      x_attr[["dim"]] <- NULL
      x_attr[["dimnames"]] <- NULL
      
    } else if (zoo::is.zoo(y)) {
      
      x_attr <- attributes(y)
      x_attr[["dim"]] <- NULL
      x_attr[["dimnames"]] <- NULL
      
    }
    
    for (i in 1:n_rows_xy) {
      
      x_subset <- x[max(1, i - width + 1):i, , drop = FALSE]
      y_subset <- y[max(1, i - width + 1):i, , drop = FALSE]
      data <- as.data.frame(cbind(y_subset, x_subset))
      
      fit <- tryCatch(pls::pcr(reformulate(termlabels = ".", response = names(data)[1]), data = data),
                      error = function(x) NA)
      df_fit <- n_cols_x
      
      if (!all(is.na(fit)) & (i >= df_fit)) {
        
        fit_coef <- coef(fit, ncomp = n_comps, intercept = TRUE)[ , , 1]
        
        result[["coefficients"]][i, ] <- fit_coef
        result[["r.squared"]][i, ] <- pls:::R2.mvr(fit, ncomp = n_comps, intercept = TRUE)$val[ , , 2]
        
      }
      
    }
    
    if (exists("x_attr")) {
      
      attributes(result[["coefficients"]]) <- x_attr
      attributes(result[["r.squared"]]) <- x_attr
      
    }
    
    attr(result[["coefficients"]], "dim") <- c(n_rows_xy, n_cols_x)
    attr(result[["r.squared"]], "dim") <- c(n_rows_xy, 1)
    
    x_dimnames <- dimnames(x)
    y_dimnames <- dimnames(y)
    
    attr(result[["coefficients"]], "dimnames") <- dimnames_lm_x(x_dimnames, n_cols_x, TRUE)
    if (length(x_dimnames) > 1) {
      attr(result[["r.squared"]], "dimnames") <- list(x_dimnames[[1]], "R-squared")
    } else {
      attr(result[["r.squared"]], "dimnames") <- list(NULL, "R-squared")
    }
    
  } else {
    
    n_rows_xy <- length(x)
    n_cols_x <- 1
    
    result <- list("coefficients" = rep(as.numeric(NA), n_rows_xy),
                   "r.squared" = rep(as.numeric(NA), n_rows_xy))
    # "std.error" = rep(as.numeric(NA), n_rows_xy))
    
    if (zoo::is.zoo(x)) {
      
      x_attr <- attributes(x)
      x_attr[["dim"]] <- NULL
      x_attr[["dimnames"]] <- NULL
      
    } else if (zoo::is.zoo(y)) {
      
      x_attr <- attributes(y)
      x_attr[["dim"]] <- NULL
      x_attr[["dimnames"]] <- NULL
      
    }
    
    for (i in 1:n_rows_xy) {
      
      x_subset <- x[max(1, i - width + 1):i]
      y_subset <- y[max(1, i - width + 1):i]
      data <- as.data.frame(cbind(y_subset, x_subset))
      
      fit <- tryCatch(pls::pcr(reformulate(termlabels = ".", response = names(data)[1]), data = data),
                      error = function(x) NA)
      df_fit <- n_cols_x
      
      summary_fit <- summary(fit)
      summary_fit_coef <- coef(summary_fit)[ , "Estimate"]
      
      if (!all(is.na(fit)) & (i >= df_fit)) {
        
        fit_coef <- coef(fit, ncomp = n_comps, intercept = TRUE)[ , , 1]
        
        result[["coefficients"]][i, ] <- fit_coef
        result[["r.squared"]][i, ] <- pls:::R2.mvr(fit, ncomp = n_comps, intercept = TRUE)$val[ , , 2]
        
      }
      
    }
    
    if (exists("x_attr")) {
      
      attributes(result[["coefficients"]]) <- x_attr
      attributes(result[["r.squared"]]) <- x_attr
      
    }
    
  }
  
  return(result)
  
}