# numerical derivatives using complex numbers
# see Squire & Trapp 1998, siam rev 40(1) 110-112

# it would seem that you can choose h to be fairly small, without
# sacrifycing accuracy due to rounding errors??

# YR 17 July 2012

lavGradientC <- function(func, x, h=1e-12, ..., check=FALSE) {

    # check current point, see if it is a scalar function
    if(check) {
        f0 <- func(x, ...)
        stopifnot(length(f0) == 1L)
    }

    nvar <- length(x)

    # determine 'h' per element of x
    h <- pmax(h, abs(h*x))

    # simple 'forward' method
    if(nvar == 1L) {
        dx <- Im(func(x + h*1i))/h
    } else {
        dx <- rep(as.numeric(NA), nvar)
        for(p in seq_len(nvar)) 
            dx[p] <- Im(func(x + h*1i*(seq.int(nvar) == p),...))/h[p]
    }

    dx
}

lavJacobianC <- function(func, x, h=(.Machine$double.eps)^(2/3), ...) {

    f0 <- func(x, ...)
    nres <- length(f0)
    nvar <- length(x)

    # determine 'h' per element of x
    h <- pmax(h, abs(h*x))

    # get exact h, per x
    tmp <- x + h
    h <- (tmp - x)

    # simple 'forward' method
    dx <- matrix(as.numeric(NA), nres, nvar)
    for(p in seq_len(nvar))
        dx[,p] <- Im(func(x + h*1i*(seq.int(nvar) == p), ...))/h[p]

    dx
}

lavJacobianD <- function(func, x, h=sqrt(.Machine$double.eps), ...) {

    f0 <- func(x, ...)
    nres <- length(f0)
    nvar <- length(x)

    # determine 'h' per element of x
    h <- pmax(h, abs(h*x))

    # get exact h, per x
    tmp <- x + h
    h <- (tmp - x)

    # simple 'forward' method
    dx <- matrix(as.numeric(NA), nres, nvar)
    for(p in seq_len(nvar))
        dx[,p] <- (func(x + h*(seq.int(nvar) == p), ...) - func(x,...))/h[p]

    dx
}

# quick and dirty (FIXME!!!) way to get
# surely there must be a more elegant way?
# dCor/dCov
# vech.idx <- lavaan:::vech.idx; diag.idx <- lavaan:::diag.idx
lav_deriv_cov2cor <- function(COV = NULL, num.idx = NULL) {

    # dCor/dvar1 = - cov / (2*var1 * sqrt(var1) * sqrt(var2))
    # dCor/dvar2 = - cov / (2*var2 * sqrt(var1) * sqrt(var2))
    # dCor/dcov  =  1/(sqrt(var1) * sqrt(var2))

    # diagonal: diag(vech(tcrossprod(1/delta)))

    nvar <- ncol(COV);  pstar <- nvar*(nvar+1)/2
    delta <- sqrt(diag(COV))
    if(length(num.idx) > 0L) {
        delta[num.idx] <- 1.0
    }

    A <- COV * -1/( 2*delta^2*tcrossprod(delta) )
    if(length(num.idx) > 0L) {
        A[num.idx,] <- 0; A[cbind(num.idx, num.idx)] <- 1
    }
    A2 <- diag(nvar) %x% t(A)

    OUT <- diag( pstar )
    diag(OUT) <- vech(tcrossprod(1/delta))
    var.idx <- which(!vech.idx(nvar) %in% vech.idx(nvar, diagonal=FALSE))
    DUP <- duplicationMatrix(nvar)
    OUT[,var.idx] <- t(DUP) %*% A2[,diag.idx(nvar)]

    if(length(num.idx) > 0L) {
        var.idx <- var.idx[-num.idx]
    }
    OUT[var.idx, var.idx] <- 0

    OUT
}


lav_deriv_cov2cor_numerical <- function(COV, num.idx=integer(0)) {

    compute.R <- function(x) {
        S <- vech.reverse(x)
        diagS <- diag(S); delta <- 1/sqrt(diagS)
        if(length(num.idx) > 0L) {
            delta[num.idx] <- 1.0
        }
        R <- diag(delta) %*% S %*% diag(delta)
        #R <- cov2cor(S)
        R.vec <- lavaan:::vech(R, diagonal = TRUE)
        R.vec
    }

    x <- vech(COV, diagonal = TRUE)
    dx <- lavJacobianC(func=compute.R, x=x)

    dx
}

