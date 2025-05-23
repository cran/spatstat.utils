#'
#'   utilseq.R
#'
#'  Utilities for sequences, vectors, ranges of values
#'
#'       $Revision: 1.25 $ $Date: 2025/04/04 07:00:04 $
#'
#'  ==>>  ORIGINAL FILE is in spatstat/develop/Spatstat/R  <<==

dropifsingle <- function(x) if(length(x) == 1) x[[1L]] else x

#' ...............  ordering ......................

# defines the current favorite algorithm for 'order' 
fave.order <- function(x) { sort.list(x, method="quick", na.last=NA) }

# order statistic (for use in lapply calls) 
orderstats <- function(x, k, decreasing=FALSE) {
  if(decreasing) k <- length(x) + 1L - k
  if(anyNA(x)) sort(x, na.last=!decreasing)[k] else sort(x, partial=k)[k]
}

# ranks (for use in lapply calls) 
orderwhich <- function(x, k, decreasing=FALSE) {
  sort.int(x, decreasing=decreasing, method="quick", na.last=NA,
            index.return=TRUE)$ix[k]
}

sortunique <- function(x) {
  ## WAS FASTER: rle(sort(x))$values
  sort(unique(x))
}

## ................ reverse cumulative sum .....................

revcumsum <- function(x) {
  #' equivalent to rev(cumsum(rev(x)))
  if(is.complex(x)) {
    a <- revcumsum(Re(x))
    b <- revcumsum(Im(x))
    return(complex(real=a, imaginary=b))
  }
  n <- length(x)
  if(identical(storage.mode(x), "integer")) {
    z <- .C(C_irevcumsum,
            x=as.integer(x),
            as.integer(n))
    return(z$x)
  } else {
    z <- .C(C_drevcumsum,
            x=as.double(x),
            as.integer(n))
    return(z$x)
  }
}

## ............. vectors of length 2 .................

as2vector <- function(x) {
  ## convert various wacky formats to numeric vector of length 2
  ## for use as coordinates of a single point.
  if(is.numeric(x)) {
    if(length(x) == 2)
      return(x)
    xname <- short.deparse(substitute(x))
    stop(paste(xname, "should have length 2"))
  }
  if(inherits(x, "ppp")) {
    #' coded so that it works even if spatstat is not loaded
    if(x$n == 1)
      return(c(x$x, x$y))
    xname <- short.deparse(substitute(x))
    stop(paste(xname, "should consist of exactly one point"))
  }
  if(is.list(x) && all(c("x", "y") %in% names(x))) {
    if(length(x$x) == 1 && length(x$y == 1))
      return(c(x$x, x$y))
    xname <- short.deparse(substitute(x))
    if(length(x$x) != 1) stop(paste0(xname, "$x should have length 1"))
    if(length(x$y) != 1) stop(paste0(xname, "$y should have length 1"))

  }
  stop(paste("Format of", sQuote(xname), "not understood"))
}

ensure2vector <- function(x) {
  if(!is.numeric(x)) {
    xname <- short.deparse(substitute(x))
    stop(paste(xname, "is not numeric"))
  }
  n <- length(x)
  if(n == 2) return(x)
  if(n == 1) return(rep(x,2))
  xname <- short.deparse(substitute(x))
  stop(paste(xname, "should be of length 1 or 2"))
}


## ............. sequences ...................

prolongseq <- function(x, newrange, step=NULL) {
  ## Extend a sequence x so that it covers the new range.
  stopifnot(length(newrange) == 2 && newrange[1L] < newrange[2L])
  ## Check 'x' is an evenly-spaced sequence
  if(length(x) > 1) {
    dx <- diff(x)
    if(any(dx <= 0))
      stop("x must be an increasing sequence")
    if(diff(range(dx)) > 0.01 * abs(mean(dx)))
      stop("x must be evenly spaced")
  }
  ## Infer step length
  if(!is.null(step)) {
    check.1.real(step)
    stopifnot(step > 0)
  } else if(length(x) > 1) {
    step <- mean(dx)
  } else stop("step is needed when x is a single value")

  ## 
  if(max(x) < newrange[1L] || min(x) > newrange[2L])
    stop("x lies entirely outside the desired range")
    
  ## add or trim data to left
  if(x[1L] > newrange[1L]) {
    leftbit <- seq(from=x[1L], to=newrange[1L], by= -step)[-1L]
    x <- c(rev(leftbit), x)
    nleft <- length(leftbit)
  } else {
    nx <- length(x)
    x <- x[x >= newrange[1L]]
    nleft <- length(x) - nx
  }

  # add or trim data to right
  nx <- length(x)
  if(newrange[2L] > x[nx]) {
    rightbit <- seq(from=x[nx], to=newrange[2L], by= step)[-1L]
    x <- c(x, rightbit)
    nright <- length(rightbit)
  } else {
    x <- x[x <= newrange[2L]]
    nright <- length(x) - nx
  }
  attr(x, "nleft") <- nleft
  attr(x, "nright") <- nright
  return(x)
}

## fill gaps in a sequence
fillseq <- function(x, step=NULL) {
  n <- length(x)
  if(n <= 1) return(x)
  rx <- range(x)
  dx <- diff(x)
  if(any(dx < 0)) {
    xname <- short.deparse(substitute(x))
    stop(paste(xname, "should be an increasing sequence"),
         call.=FALSE)
  }
  ## guess step length
  if(is.null(step)) {
    eps <- diff(rx)/1e7
    step <- min(dx[dx > eps])
  }
  ## make new sequence
  y <- seq(rx[1L], rx[2L], by=step)
  ny <- length(y)
  ## mapping from x to y
  i <- round((x - rx[1L])/step) + 1L
  i <- pmin(ny, pmax(1, i))
  return(list(xnew=y, i=i))
}

geomseq <- function(from, to, length.out) {
  check.1.real(from)
  check.1.real(to)
  if(from <= 0 || to <= 0) stop("range limits must be positive")
  y <- exp(seq(from=log(from), to=log(to), length.out=length.out))
  y[1L] <- from  #' avoid numerical error
  y[length.out] <- to
  return(y)
}

## ............. ranges ...................

intersect.ranges <- function(r, s, fatal=TRUE) {
  if(!is.null(r) && !is.null(s)) {
    lo <- max(r[1L],s[1L])
    hi <- min(r[2L],s[2L])
    if(lo <= hi)
      return(c(lo, hi))
  }
  if(fatal) stop("Intersection is empty")
  return(NULL)
}

inside.range <- function(x, r) {
  stopifnot(length(r) == 2 && r[1L] <= r[2L])
  return(x >= r[1L] & x <= r[2L])
}

check.in.range <- function(x, r, fatal=TRUE) {
  if(all(inside.range(x, r)))
    return(TRUE)
  if(!fatal)
    return(FALSE)
  xname <- deparse(substitute(x))
  stop(paste(xname, "should be a number between",
             r[1L], "and", r[2L]),
       call.=FALSE)
}

startinrange <- function(x0, dx, r) {
  ## find y = x0 + n * dx such that y \in r
  if(all(inside.range(x0, r))) return(x0)
  stopifnot(is.numeric(dx) && length(dx) == 1)
  y <- x0 + dx * round((mean(r) - x0)/dx)
  y[!inside.range(y, r)] <- NA
  return(y)
}

prettyinside <- function(x, ...) {
  r <- range(x, na.rm=TRUE)
  if(diff(r) == 0) return(r[1L])
  ## call 'pretty' after removing any NULL arguments
  p <- do.call(pretty, resolve.defaults(list(x=quote(x), ...),
                                        .MatchNull=FALSE,
                                        .StripNull=TRUE))
  ok <- inside.range(p, r)
  return(p[ok])
}

prettydiscrete <- function(x, n=10) {
  nx <- length(x)
  dx <- nx %/% n
  if(dx < 1) return(x)
  i <- 1 + (0:(n-1)) * dx
  return(x[i])
}


check.range <- function(r, fatal=TRUE) {
  if(is.numeric(r) && identical(r, range(r, na.rm=TRUE)))
    return(TRUE)
  if(!fatal)
    return(FALSE)
  rname <- deparse(substitute(r))
  stop(paste(rname, "should be a vector of length 2 giving (min, max)"))
}

evenly.spaced <- function(x, tol=1e-07) {
  ## test whether x is evenly spaced and increasing
  dx <- diff(x)
  if(any(dx <= .Machine$double.eps))
    return(FALSE)
  ## The following test for equal spacing is used in hist.default
  if(length(dx) && (diff(range(dx)) > tol * mean(dx)))
    return(FALSE)
  return(TRUE)
}

equispaced <- function(z, reltol=0.001) {
  .Deprecated("evenly.spaced")
  evenly.spaced(z, reltol)
}


adjustthinrange <- function(ur,vstep,vr) {
  if(diff(ur) >= vstep) return(ur)
  ur <- mean(ur) + c(-1,1) * vstep/2
  if(ur[1L] < vr[1L]) ur <- vr[1L] + c(0,1)*vstep
  if(ur[2L] > vr[2L]) ur <- vr[2L] - c(1,0)*vstep
  return(ur)
}

fastFindInterval <- function(x, b, labels=FALSE, reltol=0.001, dig.lab=3L,
                             left.open=TRUE) {
  nintervals <- length(b) - 1
  nx <- length(x)
  if(nx == 0) {
    y <- integer(0)
  } else if(evenly.spaced(b, reltol)) {
    ## breaks are equally spaced
    if(left.open) {
      ## intervals are left-open, right-closed ( ] except the first interval
      zz <- .C(C_fastCinterv,
               x          = as.double(x),
               n          = as.integer(nx),
               brange     = as.double(range(b)),
               nintervals = as.integer(nintervals),
               y          = as.integer(integer(nx))
               )
    } else {
      ## intervals are left-closed, right-open [ ) except the last interval
      zz <- .C(C_fastFinterv,
               x          = as.double(x),
               n          = as.integer(nx),
               brange     = as.double(range(b)),
               nintervals = as.integer(nintervals),
               y          = as.integer(integer(nx))
               )
    }
    y <- zz$y
  } else {
    ## use R's interval search algorithm
    y <- findInterval(x, b, rightmost.closed=TRUE, left.open=left.open)
  }
  if(labels) {
    #' digits in labels code copied from base::cut.default (with adaptations)
    for(dig in dig.lab:max(12L, dig.lab)) {
      ch.br <- formatC(0 + b, digits = dig, width = 1L)
      if(all(ch.br[-1L] != ch.br[1L:nintervals]))
        break
    }
    if(left.open) {
      ## left-open except the first one
      blab <- paste0(c("[", rep("(", nintervals-1)), 
                     ch.br[1:nintervals],
                     ",",
                     ch.br[-1L],
                     "]")
    } else {
      ## right-open except the last one
      blab <- paste0("[",
                     ch.br[1:nintervals],
                     ",",
                     ch.br[-1L],
                     c(rep(")", nintervals-1), "]"))
    }
    y <- as.integer(y)
    levels(y) <- as.character(blab)
    class(y) <- "factor"
  }
  return(y)
}

# ...................................................
# efficient replacements for ifelse()
# 'a' and 'b' are single values
# 'x' and 'y' are vectors of the same length as 'test'

# ifelse(test, a, b)
ifelseAB <- function(test,  a, b) {
  y <- rep.int(b, length(test))
  y[test] <- a
  return(y)
}

# ifelse(test, a, x)
ifelseAX <- function(test, a, x) {
  y <- x
  y[test] <- a
  return(y)
}

# ifelse(test, x, b)
ifelseXB <- function(test, x, b) {
  y <- rep.int(b, length(test))
  y[test] <- x[test]
  return(y)
}
  
# ifelse(test, x, y)
ifelseXY <- function(test, x, y) {
  z <- y
  z[test] <- x[test]
  return(z)
}

#.... very special cases ......

# ifelse(test, 1, NA)
ifelse1NA <- function(test) {
  y <- as.integer(test)
  y[!test] <- NA
  return(y)
}

# ifelse(test, 0, NA)
ifelse0NA <- function(test) {
  nyet <- !test
  y <- as.integer(nyet)
  y[nyet] <- NA
  return(y)
}

# ifelse(test, -x, x)
ifelseNegPos <- function(test, x) {
  y <- x
  y[test] <- -x[test]
  return(y)
}


ratiotweak <- function(a, b, overzero=NA, zerozero=overzero) {
  # map x/0 to 'overzero' and 0/0 to 'zerozero'
  result <- a/b
  bzero <- (b == 0)
  result[ bzero ] <- overzero
  if(!missing(zerozero)) {
    abzero <- bzero & (a == 0)
    result[ abzero ] <- zerozero
  }
  return(result)
}

natozero <- function(x) {
  #' map NA to zero (e.g. in tapply)
  x[is.na(x)] <- 0
  return(x)
}

insertinlist <- function(x, i, y) {
  ## insert a possibly longer or shorter list 'y'
  ## into serial position 'i' in list 'x'
  n <- length(x)
  if(n == 0) return(y)
  m <- seq_len(n)
  names(m) <- names(x)
  i <- m[[i]] # convert 'i' to integer index
  stopifnot(length(i) == 1)
  if(n == 1) return(y)
  xleft <- x[seq_len(i-1L)]
  xright <- x[i + seq_len(n-i)]
  z <- c(xleft, y, xright)
  return(z)
}

#' ..... rounding ..............................

dround <- function(x) {
  round(x, getOption('digits'))
}

niceround <- function(x, m=c(1,2,5,10)) {
  expo <- 10^as.integer(floor(log10(x)))
  y <- m * expo
  z <- y[which.min(abs(y - x))]
  return(z)
}

## ..............................................

exactCutBreaks <- function(x, breaks) {
  ## determine the exact breakpoints used in cut.default
  ## This code is extracted from base::cut.default
  stopifnot(is.numeric(x))
  if(length(breaks) > 1L) {
    ## numeric vector of breaks
    breaks <- sort.int(as.double(breaks))
    if(anyDuplicated(breaks)) 
      stop("'breaks' are not unique")
  } else if(length(breaks) == 1L) {
    ## number of breaks
    if (is.na(breaks) || breaks < 2L) 
      stop("invalid number of intervals")
    nb <- as.integer(breaks + 1)
    dx <- diff(rx <- range(x, na.rm = TRUE))
    if(dx == 0) {
      dx <- if(rx[1L] != 0) abs(rx[1L]) else 1
      breaks <- seq.int(rx[1L] - dx/1000,
                        rx[2L] + dx/1000, 
                        length.out = nb)
    } else {
      breaks <- seq.int(rx[1L], rx[2L], length.out = nb)
      breaks[c(1L, nb)] <- c(rx[1L] - dx/1000,
                             rx[2L] + dx/1000)
    }
  } else stop("breaks must be specified")
  return(breaks)
}


## ....  Harmonic mean ..........

harmonicmean <- function(x, na.rm=TRUE) {
  ## harmonic mean, robust against zeroes and small values
  if(anyNA(x)) {
    if(na.rm) x <- x[!is.na(x)] else return(NA)
  }
  if(length(x) == 0) return(NaN) # consistent with mean()
  mx <- min(abs(x))
  if(mx == 0) return(mx)
  return(mx/mean(mx/x))
}

harmonicsum <- function(x, na.rm=TRUE) {
  # computes 1/sum(1/x) robustly
  if(anyNA(x)) {
    if(na.rm) x <- x[!is.na(x)] else return(NA)
  }
  if(length(x) == 0) return(0) # consistent with sum()
  mx <- min(abs(x))
  if(mx == 0) return(mx)
  return(mx/sum(mx/x))
}
