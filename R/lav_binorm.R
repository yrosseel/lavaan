# functions to deal with bivariate normal distributions
# YR
# TODO: better handling of rho=1.0

# density of a bivariate standard normal
dbinorm <- function(u, v, rho) {
    # dirty hack to handle extreme large values for rho
    # note that u, v, and rho are vectorized!
    RHO.limit <- 0.9999
    abs.rho <- abs(rho); idx <- which(abs.rho > RHO.limit)
    if(length(idx) > 0L) rho[idx] <- sign(rho[idx]) * RHO.limit

    R <- 1 - rho*rho
    1/(2*pi*sqrt(R)) * exp( - 0.5*(u*u - 2*rho*u*v + v*v)/R )
}

# partial derivative - rho
dbinorm_drho <- function(u, v, rho) {
    R <- 1 - rho*rho
    dbinorm(u,v,rho) * (u*v*R -rho*(u*u - 2*rho*u*v + v*v) + rho*R )/R*R
}

# partial derivative - u
dbinorm_du <- function(u, v, rho) {
    R <- 1 - rho*rho
    -dbinorm(u,v,rho) * (u - rho*v)/R
}

# partial derivative - v
dbinorm_dv <- function(u, v, rho) {
    R <- 1 - rho*rho
    -dbinorm(u,v,rho) * (v - rho*u)/R
}

# CDF of bivariate standard normal
# function pbinorm(upper.x, upper.y, rho)

# partial derivative pbinorm - upper.x
pbinorm_dupper.x <- function(upper.x, upper.y, rho=0.0) {
    R <- 1 - rho*rho
    dnorm(upper.x) * pnorm( (upper.y - rho*upper.x)/R )
}

pbinorm_dupper.y <- function(upper.x, upper.y, rho=0.0) {
    R <- 1 - rho*rho
    dnorm(upper.y) * pnorm( (upper.x - rho*upper.y)/R )
}

pbinorm_drho <- function(upper.x, upper.y, rho=0.0) {
    dbinorm(upper.x, upper.y, rho)
}




# switch between pbivnorm, mnormt, ...
pbinorm <- function(upper.x=NULL, upper.y=NULL, rho=0.0,
                    lower.x=-Inf, lower.y=-Inf, check=FALSE) {

    pbinorm2(upper.x=upper.x, upper.y=upper.y, rho=rho,
             lower.x=lower.x, lower.y=lower.y, check=check)

}

# using vectorized version (a la pbivnorm)
pbinorm2 <- function(upper.x=NULL, upper.y=NULL, rho=0.0,
                     lower.x=-Inf, lower.y=-Inf, check=FALSE) {

    N <- length(upper.x)
    stopifnot(length(upper.y) == N)
    if(N > 1L) {
        if(length(rho) == 1L)
            rho <- rep(rho, N)
        if(length(lower.x) == 1L)
            lower.x <- rep(lower.x, N)
        if(length(lower.y) == 1L)
            lower.y <- rep(lower.y, N)
    }

   upper.only <- all(lower.x == -Inf & lower.y == -Inf)
   if(upper.only) {
        upper.x[upper.x == +Inf] <-  exp(10) # better pnorm?
        upper.y[upper.y == +Inf] <-  exp(10)
        upper.x[upper.x == -Inf] <- -exp(10)
        upper.y[upper.y == -Inf] <- -exp(10)
        res <- pbivnorm(upper.x, upper.y, rho=rho)
    } else {
        # pbivnorm does not handle -Inf well...
        lower.x[lower.x == -Inf] <- -exp(10)
        lower.y[lower.y == -Inf] <- -exp(10)
        res <- pbivnorm(upper.x, upper.y, rho=rho) -
               pbivnorm(lower.x, upper.y, rho=rho) -
               pbivnorm(upper.x, lower.y, rho=rho) +
               pbivnorm(lower.x, lower.y, rho=rho)
    }

    res
}


# using non-vectorized version
#pbinorm1 <- function(upper.x=NULL, upper.y=NULL, rho=0.0,
#                     lower.x=-Inf, lower.y=-Inf, check=FALSE) {
#
#    p2_i <- function(lower.x, lower.y, upper.x, upper.y, rho) {
#        # MVTNORM
#        #pmvnorm(lower=c(lower.x, lower.y),
#        #        upper=c(upper.x, upper.y),
#        #        corr=matrix(c(1,rho,rho,1),2L,2L))
#
#        # MNORMT
#        biv.nt.prob(df=0,
#                    lower=c(lower.x, lower.y),
#                    upper=c(upper.x, upper.y),
#                    mean=c(0,0),
#                    S=matrix(c(1,rho,rho,1),2L,2L))
#
#        # PBIVNORM
#
#    }
#
#    N <- length(upper.x)
#    stopifnot(length(upper.y) == N)
#    if(N > 1L) {
#        if(length(rho) == 1L)
#            rho <- rep(rho, N)
#        if(length(lower.x) == 1L)
#            lower.x <- rep(lower.x, N)
#        if(length(lower.y) == 1L)
#            lower.y <- rep(lower.y, N)
#    }
#    # biv.nt.prob does not handle +Inf well for upper
#    upper.x[upper.x == +Inf] <- exp(10) # better pnorm?
#    upper.y[upper.y == +Inf] <- exp(10) # better pnorm?
#    # biv.nt.prob does allow abs(rho) > 1
#    stopifnot(all(abs(rho) <= 1))
#
#    # vectorize (this would be faster if the loop is in the fortran code!)
#    res <- sapply(seq_len(N), function(i)
#                      p2_i(lower.x[i], lower.y[i],
#                           upper.x[i], upper.y[i],
#                           rho[i]))
#    res
#}
