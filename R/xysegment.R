#
#      xysegment.S
#
#     $Revision: 1.22 $    $Date: 2022/04/28 03:25:59 $
#
# Low level utilities for analytic geometry for line segments
#
# author: Adrian Baddeley 2001
#         from an original by Rob Foxall 1997
#
# distpl(p, l) 
#       distance from a single point p  = (xp, yp)
#       to a single line segment l = (x1, y1, x2, y2)
#
# distppl(p, l) 
#       distances from each of a list of points p[i,]
#       to a single line segment l = (x1, y1, x2, y2)
#       [uses only vector parallel ops]
#
# distppll(p, l) 
#       distances from each of a list of points p[i,]
#       to each of a list of line segments l[i,] 
#       [interpreted code uses large matrices and 'outer()']
#       [Fortran implementation included!]
#
# NNdist2segs
#       distance to nearest line segment, from each point in a list.



distpl <- function(p, l) {
  xp <- p[1]
  yp <- p[2]
  dx <- l[3]-l[1]
  dy <- l[4]-l[2]
  leng <- sqrt(dx^2 + dy^2)
  # vector from 1st endpoint to p
  xpl <- xp - l[1]
  ypl <- yp - l[2]
  # distance from p to 1st & 2nd endpoints
  d1 <- sqrt(xpl^2 + ypl^2)
  d2 <- sqrt((xp-l[3])^2 + (yp-l[4])^2)
  dmin <- min(d1,d2)
  # test for zero length
  if(leng < .Machine$double.eps)
    return(dmin)
  # rotation sine & cosine
  co <- dx/leng
  si <- dy/leng
  # back-rotated coords of p
  xpr <- co * xpl + si * ypl
  ypr <-  - si * xpl + co * ypl
  # test
  if(xpr >= 0 && xpr <= leng)
    dmin <- min(dmin, abs(ypr))
  return(dmin)
}

distppl <- function(p, l) {
  xp <- p[,1]
  yp <- p[,2]
  dx <- l[3]-l[1]
  dy <- l[4]-l[2]
  leng <- sqrt(dx^2 + dy^2)
  # vector from 1st endpoint to p
  xpl <- xp - l[1]
  ypl <- yp - l[2]
  # distance from p to 1st & 2nd endpoints
  d1 <- sqrt(xpl^2 + ypl^2)
  d2 <- sqrt((xp-l[3])^2 + (yp-l[4])^2)
  dmin <- pmin.int(d1,d2)
  # test for zero length
  if(leng < .Machine$double.eps)
    return(dmin)
  # rotation sine & cosine
  co <- dx/leng
  si <- dy/leng
  # back-rotated coords of p
  xpr <- co * xpl + si * ypl
  ypr <-  - si * xpl + co * ypl
  # ypr is perpendicular distance to infinite line
  # Applies only when xp, yp in the middle
  middle <- (xpr >= 0 & xpr <= leng)
  if(any(middle))
    dmin[middle] <- pmin.int(dmin[middle], abs(ypr[middle]))
  
  return(dmin)
}

distppll <- function(p, l, mintype=0,
                     method=c("C", "Fortran", "interpreted"), listit=FALSE) {
  np <- nrow(p)
  nl <- nrow(l)
  xp <- p[,1]
  yp <- p[,2]
  if(is.na(match(mintype,0:2)))
    stop(paste("Argument", sQuote("mintype"), "must be 0, 1 or 2.\n"))
  method <- match.arg(method)
  if(method == "Fortran") {
    warning("method='Fortran' is no longer supported; method='C' was used")
    method <- "C"
  }
  switch(method,
         interpreted={
           dx <- l[,3]-l[,1]
           dy <- l[,4]-l[,2]
           # segment lengths
           leng <- sqrt(dx^2 + dy^2)
           # rotation sines & cosines
           co <- dx/leng
           si <- dy/leng
           co <- matrix(co, nrow=np, ncol=nl, byrow=TRUE)
           si <- matrix(si, nrow=np, ncol=nl, byrow=TRUE)
           # matrix of squared distances from p[i] to 1st endpoint of segment j
           xp.x1 <- outer(xp, l[,1], "-")
           yp.y1 <- outer(yp, l[,2], "-")
           d1 <- xp.x1^2 + yp.y1^2
           # ditto for 2nd endpoint
           xp.x2 <- outer(xp, l[,3], "-")
           yp.y2 <- outer(yp, l[,4], "-")
           d2 <- xp.x2^2 + yp.y2^2
           # for each (i,j) rotate p[i] around 1st endpoint of segment j
           # so that line segment coincides with x axis
           xpr <- xp.x1 * co + yp.y1 * si
           ypr <-  - xp.x1 * si + yp.y1 * co
           d3 <- ypr^2
           # test
           lenf <- matrix(leng, nrow=np, ncol=nl, byrow=TRUE)
           zero <- (lenf < .Machine$double.eps) 
           outside <- (zero | xpr < 0 | xpr > lenf) 
           if(any(outside))
             d3[outside] <- Inf

           dsq <- matrix(pmin.int(d1, d2, d3),nrow=np, ncol=nl)
           d <- sqrt(dsq)
           if(mintype >= 1)
             min.d <- apply(d, 1, min)
           if(mintype == 2)
             min.which <- apply(d, 1, which.min)
         },
         C = {
           eps <- .Machine$double.eps
           temp <- .C(C_prdist2segs,
                      x=as.double(xp),
                      y=as.double(yp),
                      npoints =as.integer(np),
                      x0=as.double(l[,1]),
                      y0=as.double(l[,2]),
                      x1=as.double(l[,3]),
                      y1=as.double(l[,4]),
                      nsegments=as.integer(nl),
                      epsilon=as.double(eps),
                      dist2=as.double(numeric(np * nl)))
           d <- matrix(sqrt(temp$dist2), nrow=np, ncol=nl)
           if(mintype == 1) {
             min.d <- apply(d, 1, min)
           } else if(mintype == 2) {
             min.which <- apply(d, 1, which.min)
             min.d <- d[cbind(1:np, min.which)]
           }
         })
  ###### end switch #####
  if(mintype==0)
    return(if(listit) list(d=d) else d)
  else if(mintype==1)
    return(list(d=d, min.d=min.d))
  else if(mintype==2) 
    return(list(d=d, min.d=min.d, min.which=min.which))
}

# (distance to) nearest segment 

distppllmin <- function(p, l, big=NULL) {
  if(is.null(big)) {
    xdif <- diff(range(c(p[,1],l[, c(1,3)])))
    ydif <- diff(range(c(p[,2],l[, c(2,4)])))
    big <- 2 * (xdif^2 + ydif^2)
  }
  z <- NNdist2segments(p[,1], p[,2],
                       l[,1], l[,2], l[,3], l[,4],
                       big)
  return(list(min.d=sqrt(z$dist2), min.which=z$index))
}

NNdist2segments <- function(xp, yp, x0, y0, x1, y1, bigvalue, wantindex=TRUE) {
  np <- length(xp)
  ns <- length(x0)
  dist2 <- rep(bigvalue, np)
  if(wantindex) {
    z <- .C(C_nndist2segs,
            xp=as.double(xp),
            yp=as.double(yp),
            npoints=as.integer(np),
            x0=as.double(x0),
            y0=as.double(y0),
            x1=as.double(x1),
            y1=as.double(y1),
            nsegments=as.integer(ns),
            epsilon=as.double(.Machine$double.eps),
            dist2=as.double(dist2),
            index=as.integer(integer(np)))
    return(list(dist2=z$dist2,
                index=z$index + 1L))
  } else {
    z <- .C(C_nnd2segs,
            xp=as.double(xp),
            yp=as.double(yp),
            npoints=as.integer(np),
            x0=as.double(x0),
            y0=as.double(y0),
            x1=as.double(x1),
            y1=as.double(y1),
            nsegments=as.integer(ns),
            epsilon=as.double(.Machine$double.eps),
            dist2=as.double(dist2))
    return(list(dist2=z$dist2))
  }
}
