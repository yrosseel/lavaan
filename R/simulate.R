# simulate data starting from a user-specified model
#
# initial version: YR 24 jan 2011
# revision for 0.4-11: YR 21 okt 2011
simulateData <- function(
                         # user-specified model
                         model           = NULL,
                         model.type      = "sem",

                         # model modifiers
                         meanstructure   = FALSE,
                         int.ov.free     = TRUE,
                         int.lv.free     = FALSE,
                         fixed.x         = FALSE,
                         orthogonal      = FALSE,
                         std.lv          = TRUE,

                         auto.fix.first  = FALSE,
                         auto.fix.single = FALSE,
                         auto.var        = TRUE,
                         auto.cov.lv.x   = TRUE,
                         auto.cov.y      = TRUE,
                         ...,

                         # data properties
                         sample.nobs     = 500L,
                         ov.var          = NULL,
                         group.label     = paste("G", 1:ngroups, sep=""),
                         skewness        = NULL,
                         kurtosis        = NULL,

                         # control
                         seed = NULL,
                         empirical = FALSE,

                         return.type = "data.frame"
                        )
{
    if(!is.null(seed)) set.seed(seed)
    if(!exists(".Random.seed", envir = .GlobalEnv))
        runif(1)               # initialize the RNG if necessary
    RNGstate <- .Random.seed

    # lavaanify
    lav <- lavaanify(model = model, 
                     meanstructure=meanstructure,
                     int.ov.free=int.ov.free, 
                     int.lv.free=int.lv.free,
                     fixed.x=fixed.x,
                     orthogonal=orthogonal,
                     std.lv=std.lv,
                     auto.fix.first=auto.fix.first,
                     auto.fix.single=auto.fix.single,
                     auto.var=auto.var,
                     auto.cov.lv.x=auto.cov.lv.x,
                     auto.cov.y=auto.cov.y)

    # unstandardize 
    if(!is.null(ov.var)) {
        # FIXME: if ov.var is named, check the order of the elements

        # 1. unstandardize observed variables
        lav$ustart <- unstandardize.est.ov(partable=lav, ov.var=ov.var)

        # 2. unstandardized latent variables
    }

    # run lavaan to set up the model matrices
    fit <- lavaan(model=lav, sample.nobs=sample.nobs, ...)

    # the model-implied moments for the population
    Sigma.hat <- computeSigmaHat(fit@Model)
       Mu.hat <- computeMuHat(fit@Model)

    # ngroups
    ngroups <- length(sample.nobs)

    # prepare 
    X <- vector("list", length=ngroups)
    out <- vector("list", length=ngroups)

    for(g in 1:ngroups) {
        # FIXME: change to rmvnorm once we include the library?
        if(is.null(skewness) && is.null(kurtosis)) {
            X[[g]] <- MASS.mvrnorm(n = sample.nobs[g],
                                   mu = Mu.hat[[g]],
                                   Sigma = Sigma.hat[[g]],
                                   empirical = empirical)
        } else {
            # first generate Z
            Z <- ValeMaurelli1983(n        = sample.nobs[g], 
                                  COR      = cov2cor(Sigma.hat[[g]]), 
                                  skewness = skewness,  # FIXME: per group?
                                  kurtosis = kurtosis)
            # rescale
            X[[g]] <- scale(Z, center = Mu.hat[[g]],
                               scale  = 1/sqrt(diag(Sigma.hat[[g]])))
        }

        if(return.type == "data.frame") X[[g]] <- as.data.frame(X[[g]])
    }

    if(return.type == "matrix") {
        if(ngroups == 1L) {
            return(X[[1L]])
        } else {
            return(X)
        }

    } else if (return.type == "data.frame") {
        Data <- X[[1L]]

        # if multiple groups, add group column
        if(ngroups > 1L) {
            for(g in 2:ngroups) {
                Data <- rbind(Data, X[[g]])
            }
            Data$group <- rep(1:ngroups, times=sample.nobs)
        }
        var.names <- vnames(fit@ParTable, type="ov", group=1L)
        if(ngroups > 1L) var.names <- c(var.names, "group")
        names(Data) <- var.names
        return(Data)

    } else if (return.type == "cov") {
        if(ngroups == 1L) {
            return(cov(X[[1L]]))
        } else {
            cov.list <- lapply(X, cov)
            return(cov.list)
        }
    }
}

Skewness <- function(x., N1=TRUE) {
    x <- x.; x <- x[!is.na(x)]; N <- length(x)
    mean.x <- mean(x); xc <- x - mean.x; var.x <- var(x)
    if(!N1) var.x <- var.x * (N-1)/N
    sd.x <- sqrt(var.x)
    sk <- sum(xc^3)/sd.x^3
    skewness <- N*sk/((N-1)*(N-2))
    skewness
}

Kurtosis <- function(x., N1=TRUE) {
    x <- x.; x <- x[!is.na(x)]; N <- length(x)
    mean.x <- mean(x); xc <- x - mean.x; var.x <- var(x)
    if(!N1) var.x <- var.x * (N-1)/N
    k <- sum(xc^4)/var.x^2
    kurtosis <- N*(N+1)*k/((N-1)*(N-2)*(N-3))-3*(N-1)^2/((N-2)*(N-3))
    kurtosis
}

fleishman1978 <- function(n=100, skewness=0, kurtosis=0, verbose=FALSE) {

    system.function <- function(x, skewness, kurtosis) {
        b=x[1L]; c=x[2L]; d=x[3L]
        eq1 <- b^2 + 6*b*d + 2*c^2 + 15*d^2 - 1
        eq2 <- 2*c*(b^2 + 24*b*d + 105*d^2 + 2) - skewness
        eq3 <- 24*(b*d + c^2*(1 + b^2 + 28*b*d) +
                   d^2*(12 + 48*b*d + 141*c^2 + 225*d^2)) - kurtosis
        eq <- c(eq1,eq2,eq3)
        sum(eq^2) ## SS
    }

    out <- nlminb(start=c(1,0,0), objective=system.function,
                  scale=10,
                  control=list(trace=ifelse(verbose,1,0)),
                  skewness=skewness, kurtosis=kurtosis)
    if(out$convergence != 0) warning("no convergence")
    b <- out$par[1L]; c <- out$par[2L]; d <- out$par[3L]; a <- -c

    Z <- rnorm(n=n)
    Y <- a + b*Z + c*Z^2 + d*Z^3
    Y
}

ValeMaurelli1983 <- function(n=100L, COR, skewness, kurtosis) {

    fleishman1978_abcd <- function(skewness, kurtosis) {
        system.function <- function(x, skewness, kurtosis) {
            b.=x[1L]; c.=x[2L]; d.=x[3L]
            eq1 <- b.^2 + 6*b.*d. + 2*c.^2 + 15*d.^2 - 1
            eq2 <- 2*c.*(b.^2 + 24*b.*d. + 105*d.^2 + 2) - skewness
            eq3 <- 24*(b.*d. + c.^2*(1 + b.^2 + 28*b.*d.) +
                       d.^2*(12 + 48*b.*d. + 141*c.^2 + 225*d.^2)) - kurtosis
            eq <- c(eq1,eq2,eq3)
            sum(eq^2) ## SS
        }

        out <- nlminb(start=c(1,0,0), objective=system.function,
                      scale=10,
                      control=list(trace=0),
                      skewness=skewness, kurtosis=kurtosis)
        if(out$convergence != 0) warning("no convergence")
        b. <- out$par[1L]; c. <- out$par[2L]; d. <- out$par[3L]; a. <- -c.
        c(a.,b.,c.,d.)
    }

    getICOV <- function(b1, c1, d1, b2, c2, d2, R) {
        objectiveFunction <- function(x, b1, c1, d1, b2, c2, d2, R) {
            rho=x[1L]
            eq <- rho*(b1*b2 + 3*b1*d2 + 3*d1*b2 + 9*d1*d2) +
                  rho^2*(2*c1*c2) + rho^3*(6*d1*d2) - R
            eq^2
        }

        #gradientFunction <- function(x, bcd1, bcd2, R) {
        #
        #}

        out <- nlminb(start=R, objective=objectiveFunction,
                      scale=10, control=list(trace=0),
                      b1=b1, c1=c1, d1=d1, b2=b2, c2=c2, d2=d2, R=R)
        if(out$convergence != 0) warning("no convergence")
        rho <- out$par[1L]
        rho
    }

    # number of variables
    nvar <- ncol(COR)
    # check skewness
    if(length(skewness) == nvar) {
        SK <- skewness
    } else if(length(skewness == 1L)) {
        SK <- rep(skewness, nvar)
    } else {
        stop("skewness has wrong length")
    }

    if(length(kurtosis) == nvar) {
        KU <- kurtosis
    } else if(length(skewness == 1L)) {
        KU <- rep(kurtosis, nvar)
    } else {
        stop("kurtosis has wrong length")
    }

    # create Fleishman table
    FTable <- matrix(0, nvar, 4L)
    for(i in 1:nvar) {
        FTable[i,] <- fleishman1978_abcd(skewness=SK[i], kurtosis=KU[i])
    }

    # compute intermediate correlations between all pairs
    ICOR <- diag(nvar)
    for(j in 1:(nvar-1L)) {
        for(i in (j+1):nvar) {
            if(COR[i,j] == 0) next
            ICOR[i,j] <- ICOR[j,i] <-
                getICOV(FTable[i,2], FTable[i,3], FTable[i,4],
                        FTable[j,2], FTable[j,3], FTable[j,4], R=COR[i,j])
        }
    }

    ICOR

    # generate Z ## FIXME: replace by rmvnorm once we use that package
    X <- Z <- MASS.mvrnorm(n=n, mu=rep(0,nvar), Sigma=ICOR)

    # transform Z using Fleishman constants
    for(i in 1:nvar) {
        X[,i] <- FTable[i,1L] + FTable[i,2L]*Z[,i] + FTable[i,3L]*Z[,i]^2 +
                 FTable[i,4L]*Z[,i]^3
    }

    X
}

