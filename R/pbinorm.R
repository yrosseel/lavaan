# bivariate standard normal cdf  VECTORIZED!!
# 
# upper.xx and upper.y are the UPPER limits of the bivariate integral
# optionally, lower.x and lowerly can be provided
#
# this is using the MVTDST routine from package mvtnorm, but avoiding
# the overhead of pmvnorm
#
# YR 26 June 2012

pbinorm <- function(upper.x=NULL, upper.y=NULL, rho=0.0, 
                    lower.x=-Inf, lower.y=-Inf, check=FALSE) {

    p2_i <- function(lower.x, lower.y, upper.x, upper.y, rho, infin) {
        .Fortran("MVTDST", N=2L, NU=0L,
                 lower=c(lower.x, lower.y), upper=c(upper.x, upper.y),
                 infin=infin, correl=c(rho),
                 delta=c(0,0), maxpts=25000L, abseps=0.001, releps=0,
                 error=0, value=0, inform=0L, PACKAGE="mvtnorm")$value
    }
    
    if(check) {
        lower.x <- as.double(lower.x); lower.y <- as.double(lower.y)
        upper.x <- as.double(upper.x); upper.y <- as.double(upper.y)
        
        rho <- as.double(rho)
        stopifnot(any(abs(rho) < 1), length(lower.y) == N,
              upper.x == N, upper.y == N)
    }

    N <- length(lower.x)
    if(N > 1L && length(rho) == 1L)
        rho <- rep(rho, N)

    # set up infin
    #INFIN  INTEGER, array of integration limits flags:
    #         if INFIN(I) < 0, Ith limits are (-infinity, infinity);
    #         if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
    #         if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
    #         if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
    infin.x <- rep(2L, N); infin.y <- rep(2L, N)
    infin.x[ lower.x < -10.0 & upper.x > +10.0 ] <- -1L
    infin.x[ lower.x < -10.0 & upper.x < +10.0 ] <-  0L
    infin.x[ lower.x > -10.0 & upper.x > +10.0 ] <-  1L
    infin.y[ lower.y < -10.0 & upper.y > +10.0 ] <- -1L
    infin.y[ lower.y < -10.0 & upper.y < +10.0 ] <-  0L
    infin.y[ lower.y > -10.0 & upper.y > +10.0 ] <-  1L
    lower.x[ !is.finite(lower.x) ] <- 0
    lower.y[ !is.finite(lower.y) ] <- 0
    upper.x[ !is.finite(upper.x) ] <- 0
    upper.y[ !is.finite(upper.y) ] <- 0

    # vectorize (this would be faster if the loop is in the fortran code!) 
    res <- sapply(seq_len(N), function(i) 
                      p2_i(lower.x[i], lower.y[i], 
                           upper.x[i], upper.y[i], rho[i],
                           infin=c(infin.x[i], infin.y[i])))

    res
}

pbinorm2 <- function(x, y, rho) {

    N <- length(x)

    # check input
    stopifnot(any(abs(rho) < 1), length(y) == N)
    if(N > 1L && length(rho) == 1L)
        rho <- rep(rho, N)

    # vectorize 
    res <- sapply(seq_len(N), function(i) bivnor(ah=-x[i], ak=-y[i], r=rho[i]))

    res
}

# the bivnor function below is an R translation of MATLAB code written
# by John Burkardt: 
# http://people.sc.fsu.edu/~jburkardt/c_src/toms462/toms462.html
#
# R port by Yves Rosseel 26 June 2012

#*****************************************************************************80
#
## BIVNOR computes the bivariate normal CDF.
#
#  Discussion:
#
#    BIVNOR computes the probability for two normal variates X and Y
#    whose correlation is R, that AH <= X and AK <= Y.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    13 April 2012
#
#  Author:
#
#    Original FORTRAN77 version by Thomas Donnelly.
#    MATLAB version by John Burkardt.
#
#  Reference:
#
#    Thomas Donnelly,
#    Algorithm 462: Bivariate Normal Distribution,
#    Communications of the ACM,
#    October 1973, Volume 16, Number 10, page 638.
#
#  Parameters:
#
#    Input, real AH, AK, the lower limits of integration.
#
#    Input, real R, the correlation between X and Y.
#
#    Output, real VALUE, the bivariate normal CDF.
#
#  Local Parameters:
#
#    Local, integer IDIG, the number of significant digits
#    to the right of the decimal point desired in the answer.
#

bivnor <- function(ah, ak, r, idig=15) {
    b <- 0.0

    gh <- pnorm(-ah)/2.0
    gk <- pnorm(-ak)/2.0

    if(r == 0.0) {
        b <- 4.00 * gh * gk
        b <- max(b, 0.0)
        b <- min(b, 1.0)
        return(b)
    }

    rr <- (1.0 + r) * (1.0 - r)
    if(rr < 0.0) 
        stop("1 < |R|")

    if(rr == 0.0) {
        if(r < 0.0) {
            if(ah + ak < 0.0) {
                b <- 2.0 * (gh + gk) - 1.0
            } else if ( ah - ak < 0.0 ) {
                b <- 2.0 * gk
            } else {
                b <- 2.0 * gh
            }
        }
        b <- max(b, 0.0)
        b <- min(b, 1.0)
        return(b)
    }

    sqr <- sqrt(rr)

    if(idig == 15) {
        con <- 2.0 * pi * 1.0E-15 / 2.0
    } else {
        con <- pi
        for(i in 1:idig) {
          con <- con / 10.0
        }
    }
#
#  (0,0)
#
    if(ah == 0.0 && ak == 0.0) {
        b <- 0.25 + 0.5 * asin(r) / pi
        b <- max(b, 0.0)
        b <- min(b, 1.0)
        return(b)
    }
#
#  (0,nonzero)
#
    if(ah == 0.0 && ak != 0.0) {
        b <- gk
        wh <- -ak
        wk <- (ah / ak - r) / sqr
        gw <- 2.0 * gk
        is <- 1
#
#  (nonzero,0)
#
    } else if(ah != 0.0 && ak == 0.0) {
        b <- gh
        wh <- -ah
        wk <- (ak / ah - r) / sqr
        gw <- 2.0 * gh
        is <- -1
#
#  (nonzero,nonzero)
#
    } else if(ah != 0.0 && ak != 0.0) {
        b <- gh + gk
        if(ah*ak < 0.0) {
            b <- b - 0.5
        }
        wh <- - ah
        wk <- (ak / ah - r) / sqr
        gw <- 2.0 * gh
        is <- -1
    }

    while(TRUE) {
        sgn <- -1.0
        t <- 0.0

        if(wk != 0.0) {
            if(abs(wk) == 1.0) {
                t <- wk * gw * ( 1.0 - gw ) / 2.0
                b <- b + sgn * t
            } else {
                if (1.0 < abs(wk)) {
                    sgn <- -sgn
                    wh <- wh * wk
                    g2 <- pnorm(wh)
                    wk <- 1.0 / wk

                    if(wk < 0.0) {
                        b <- b + 0.5
                    }
                    
                    b <- b - ( gw + g2 ) / 2.0 + gw * g2
                }

                h2 <- wh * wh
                a2 <- wk * wk
                h4 <- h2 / 2.0
                ex <- exp(-h4)
                w2 <- h4 * ex
                ap <- 1.0
                s2 <- ap - ex
                sp <- ap
                s1 <- 0.0
                sn <- s1
                conex <- abs(con / wk)

                while(TRUE) {

                    cn <- ap * s2 / ( sn + sp )
                    s1 <- s1 + cn

                    if(abs(cn) <= conex)
                        break

                    sn <- sp
                    sp <- sp + 1.0
                    s2 <- s2 - w2
                    w2 <- w2 * h4 / sp
                    ap <- - ap * a2
                }

               t <- 0.5 * ( atan(wk) - wk * s1 ) / pi
               b <- b + sgn * t
            }
        }      

        if(0 <= is)
            break

        if(ak == 0.0)
            break

        wh <- -ak
        wk <- ( ah / ak - r ) / sqr
        gw <- 2.0 * gk
        is <- 1
    }

    b <- max(b, 0.0)
    b <- min(b, 1.0)

    b
}
