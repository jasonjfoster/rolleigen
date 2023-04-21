rollapplyr_pcr <- function(x, y, width, n_comps) {
  
  # could not find function "mvrValstats"
  require(pls)
  
  # if (!requireNamespace("pls", quietly = TRUE)) {
  #   stop("pls package required for this function")
  # }
  
  if (is.matrix(y)) {
    
    n_rows_xy <- nrow(x)
    n_cols_x <- ncol(x) + 1
    
    result <- list("coefficients" = matrix(as.numeric(NA), n_rows_xy, n_cols_x),
                   "r.squared" = matrix(as.numeric(NA), n_rows_xy, 1))
    
    x_attr <- attributes(x)
    x_attr[["dim"]] <- NULL
    x_attr[["dimnames"]] <- NULL
    
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
      
      attr(result[["coefficients"]], "dim") <- c(n_rows_xy, n_cols_x)
      attr(result[["r.squared"]], "dim") <- c(n_rows_xy, 1)
      
    }
    
  } # else {
  #   
  #   n_rows_xy <- length(x)
  #   n_cols_x <- 1
  #   
  #   result <- list("coefficients" = rep(as.numeric(NA), n_rows_xy),
  #                  "r.squared" = rep(as.numeric(NA), n_rows_xy),
  #                  "std.error" = rep(as.numeric(NA), n_rows_xy))
  #   
  #   
  #   if (zoo::is.zoo(x)) {
  #     
  #     x_attr <- attributes(x)
  #     x_attr[["dim"]] <- NULL
  #     x_attr[["dimnames"]] <- NULL
  #     
  #   } else if (zoo::is.zoo(y)) {
  #     
  #     x_attr <- attributes(y)
  #     x_attr[["dim"]] <- NULL
  #     x_attr[["dimnames"]] <- NULL
  #     
  #   }
  #   
  #   for (i in 1:n_rows_xy) {
  #     
  #     x_subset <- x[max(1, i - width + 1):i]
  #     y_subset <- y[max(1, i - width + 1):i]
  #     data <- as.data.frame(cbind(y_subset, x_subset))
  #     
  #     if (intercept) {
  #       fit <- lm(reformulate(termlabels = ".", response = names(data)[1]), data = data)
  #     } else {
  #       fit <- lm(reformulate(termlabels = ".-1", response = names(data)[1]), data = data)
  #     }
  #     
  #     summary_fit <- summary(fit)
  #     summary_fit_coef <- coef(summary_fit)[ , "Estimate"]
  #     
  #     if (nrow(coef(summary_fit)) == n_cols_x) {
  #       
  #       result[["coefficients"]][i] <- summary_fit_coef
  #       
  #       # "In summary.lm(fit) : essentially perfect fit: summary may be unreliable"
  #       if (!(isTRUE(all.equal(as.numeric(rep(summary_fit_coef[1], length(y_subset))), as.numeric(y_subset))) &&
  #             isTRUE(all.equal(as.numeric(summary_fit_coef[-1]), rep(0, length(summary_fit_coef[-1])))))) {
  #         
  #         result[["r.squared"]][i] <- summary_fit$r.squared
  #         result[["std.error"]][i] <- coef(summary_fit)[ , "Std. Error"]
  #         
  #       } else {
  #         
  #         result[["r.squared"]][i] <- as.numeric(NA)
  #         result[["std.error"]][i] <- as.numeric(NA)
  #         
  #       }
  #       
  #     }
  #     
  #   }
  #   
  #   if (exists("x_attr")) {
  #     
  #     attributes(result[["coefficients"]]) <- x_attr
  #     attributes(result[["r.squared"]]) <- x_attr
  #     attributes(result[["std.error"]]) <- x_attr
  #     
  #   }
  #   
  # }
  
  return(result)
  
}