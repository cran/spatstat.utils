#'
#'   utiloptim.R
#'
#'   Utilities for optimization
#'
#'  $Revision: 1.4 $  $Date: 2021/10/03 08:17:51 $
#'

optimizeWithTrace <- local({

  tracer <- function(x, ..., .TheFunction, .Enviro) {
    y <- .TheFunction(x, ...)
    xx <- get("xx", envir=.Enviro)
    yy <- get("yy", envir=.Enviro)
    assign("xx", c(xx, as.numeric(x)), envir=.Enviro)
    assign("yy", c(yy, y), envir=.Enviro)
    return(y)
  }
  
  optimizeWithTrace <- function(f, interval, ..., 
                                lower = min(interval), upper = max(interval)) {
    e <- new.env()
    assign("xx", numeric(0), envir=e)
    assign("yy", numeric(0), envir=e)
    result <- optimize(tracer, lower=lower, upper=upper,
                       ..., .TheFunction=f, .Enviro=e)
    result$x <- get("xx", envir=e)
    result$y <- get("yy", envir=e)
    rm(e)
    return(result)
  }

  optimizeWithTrace
})

which.min.fair <- function(x) {
  a <- min(x, na.rm=TRUE)
  i <- which(x == a)
  if(length(i) > 1)
    i <- sample(i, 1)
  return(i)
}

which.max.fair <- function(x) {
  a <- max(x, na.rm=TRUE)
  i <- which(x == a)
  if(length(i) > 1)
    i <- sample(i, 1)
  return(i)
}

