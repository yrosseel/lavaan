# this file is from the 'numDeriv' package
# it is included to avoid an extra package dependency that will
# most likely disappear again in the next release

jacobian <- function (func, x, method="Richardson",
                              method.args=list(), ...) UseMethod("jacobian")

jacobian.default <- function(func, x, method="Richardson",
      method.args=list(eps=1e-4, d=0.0001,
      zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE), ...){
  f <- func(x, ...)
  n <- length(x)         #number of variables.
  if(method=="simple"){
    #  very simple numerical approximation
    args <- list(eps=1e-4) # default
    args[names(method.args)] <- method.args
    eps <- args$eps
    df <-matrix(NA, length(f), n)
    for (i in 1:n) {
      dx <- x
      dx[i] <- dx[i] +eps
      df[,i] <- (func(dx, ...)-f)/eps
     }
    return(df)
    } else
  if(method=="Richardson"){
    args <- list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE) # default
    args[names(method.args)] <- method.args
    eps <- args$eps
    d <- args$d
    r <- args$r
    v <- args$v
    a <- array(NA, c(length(f),r, n) )

    h <- abs(d*x)+eps*(abs(x) < args$zero.tol)
    for(k in 1:r)  { # successively reduce h                
       for(i in 1:n)  {
         a[,k,i] <- (func(x + h*(i==seq(n)), ...) -
                     func(x - h*(i==seq(n)), ...))/(2*h[i])
         #if((k != 1)) a[,(abs(a[,(k-1),i]) < 1e-20)] <- 0 #some func are unstable near zero
         }
       h <- h/v     # Reduced h by 1/v.
       }
   for(m in 1:(r - 1)) {
       a <- (a[,2:(r+1-m),,drop=FALSE]*(4^m)-a[,1:(r-m),,drop=FALSE])/(4^m-1)
     }
  # drop second dim of a, which is now 1 (but not other dim's even if they are 1
  return(array(a, dim(a)[c(1,3)]))
  } else stop("indicated method ", method, "not supported.")
}

