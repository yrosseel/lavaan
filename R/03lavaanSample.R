# constructor for the 'Sample' class
#
# initial version: YR 25/03/2009
# major revision: YR 5/11/2011: separate data.obs and sample statistics

# extract the data we need for this particular model
getData <- function(data        = NULL, 
                    ov.names    = character(0),

                    # standardize?
                    std.ov      = FALSE,

                    # multiple groups?
                    group       = NULL,
                    group.label = character(0)
                   ) 
{
    # number of observed variables
    nvar  <- length(ov.names)

    # number of groups
    ngroups <- 1L
    if(length(group.label) > 1L) ngroups <- length(group.label) 

    # prepare empty list for data.matrix per group
    X <- vector("list", length=ngroups)

    # does the data contain all the observed variables
    # needed in the user-specified model for this group
    idx.missing <- which(!(ov.names %in% names(data)))
    if(length(idx.missing)) {
        stop("lavaan ERROR: missing observed variables in dataset: ",
             paste(ov.names[idx.missing], collapse=" "))
    }

    # for each group
    for(g in 1:ngroups) {

        # extract variables in correct order
        var.idx <- match(ov.names, names(data))
        if(ngroups > 1L) {
            case.idx <- data[, group] == group.label[g]
            data.tmp <- data[case.idx, var.idx]
        } else {
            data.tmp <- data[, var.idx]
        }

        # check if we have enough observations
        if(nrow(data.tmp) < nvar)
            stop("lavaan ERROR: too few observations (nobs < nvar)")

        # data should contain numeric values only
        data.tmp <- data.matrix(data.tmp)

        # standardize observed variables?
        if(std.ov) {
            data.tmp <- scale(data.tmp)[,]
        }

        X[[g]] <- data.tmp

    } # ngroups

    X
}

getMissingPatterns <- function(X       = NULL, 
                               missing = "listwise",
                               warn    = TRUE,
                               verbose = FALSE) {

    # number of groups
    ngroups <- length(X)

    # set missing flag (can be overriden later)
    if(missing == "ml") {
        missing.flag <- rep(TRUE, ngroups)
    } else {
        missing.flag <- rep(FALSE, ngroups)
    }

    # prepare empty list for missing data
    missing <- vector("list", length=ngroups)

    for(g in 1:ngroups) {
        missing[[g]] <- list()
        
        if(!missing.flag[g]) {
            # listwise deletion
            keep.idx <- complete.cases(X[[g]])
            nobs <- sum(keep.idx)
            # check again if we have enough observations
            if(nobs == 0L)
                stop("lavaan ERROR: no cases left after listwise deletion")
            if(nobs < ncol(X[[g]]))
                stop("lavaan ERROR: too few observations (nobs < nvar)")
            # fill in some basic information in missing
            missing[[g]] <- list(npatterns=0L, flag=FALSE)
        } else {
            # do more (but only if missing!)
            #   - get missing patterns
            #   - store sufficient statistics (per missing pattern group)
            #   - compute pairwise coverage
            missing[[g]] <- missing.patterns(X[[g]], warn=warn)
            if(missing[[g]]$npatterns > 1L) {
                # estimate moments
                #out <- estimate.moments.fiml(X=data.obs, M=missing[[g]],
                #                             verbose=verbose)
                out <- estimate.moments.EM(X=X[[g]], M=missing[[g]],
                                           verbose=verbose)
                missing[[g]]$sigma <- out$sigma
                missing[[g]]$mu    <- out$mu
                missing[[g]]$h1    <- out$fx
                missing[[g]]$flag  <- TRUE
            } else {
                # data is complete after all (for this group)
                missing[[g]]$flag  <- FALSE
            }
        }

    } # ngroups

    missing
}

getSampleStats <- function(X           = NULL,
                             M           = NULL,
                             boot.idx    = NULL,
                             rescale     = FALSE) {

    # number of groups
    ngroups <- length(X)

    # sample statistics per group
    cov         <- vector("list", length=ngroups)
    var         <- vector("list", length=ngroups)
    mean        <- vector("list", length=ngroups)
    nobs        <- vector("list", length=ngroups)
    norig       <- vector("list", length=ngroups)
    ov.names    <- vector("list", length=ngroups)

    for(g in 1:ngroups) {

        # get variable names for this group
        ov.names[[g]] <- colnames(X[[g]])
        X[[g]] <- unname(X[[g]])

        # bootstrap sample?
        if(!is.null(boot.idx)) {
            X[[g]] <- X[[g]][boot.idx[[g]]]
        }

        # listwise deletion?
        norig[[g]] <- nrow(X[[g]])
        if(is.null(M)) {
            keep.idx <- complete.cases(X[[g]])
            X[[g]] <- X[[g]][keep.idx,,drop=TRUE]
        }
        nobs[[g]] <- nrow(X[[g]])

        # fill in the other slots
        cov[[g]]  <-   cov(X[[g]], use="pairwise") # must be pairwise
        #var[[g]]  <- apply(X[[g]], 2,  var, na.rm=TRUE)
        mean[[g]] <- apply(X[[g]], 2, mean, na.rm=TRUE)

        # rescale cov by (N-1)/N?
        if(rescale) {
            # we 'transform' the sample cov (divided by n-1) 
            # to a sample cov divided by 'n'
            cov[[g]] <- (nobs[[g]]-1)/nobs[[g]] * cov[[g]]
        }

    } # ngroups

    # construct SampleStats object
    SampleStats <- new("SampleStats",

                         # sample moments
                         mean           = mean,
                         cov            = cov,
                         #var            = var,

                         # convenience
                         nobs           = nobs,
                         norig          = norig,
                         ngroups        = ngroups,
                         ntotal         = sum(unlist(nobs)),
                         ov.names       = ov.names,

                         # missingness
                         missing        = M
                        )

    SampleStats
}


getSampleStatsExtra <- function(X             = NULL, 
                                sample        = NULL,

                                estimator     = "ML",
                                mimic         = "lavaan",
                                meanstructure = FALSE) {

    # number of groups
    ngroups <- sample@ngroups

    # extra sample statistics per group
    icov        <- vector("list", length=ngroups)
    cov.log.det <- vector("list", length=ngroups)
    cov.vecs    <- vector("list", length=ngroups)
    WLS.V       <- vector("list", length=ngroups)

    # icov and cov.log.det
    for(g in 1:ngroups) {
        tmp <- try(inv.chol(sample@cov[[g]], logdet=TRUE))
        if(inherits(tmp, "try-error")) {
            if(ngroups > 1) {
                stop("sample covariance can not be inverted in group", g)
            } else {
                stop("sample covariance can not be inverted")
            }
        } else {
            cov.log.det[[g]] <- attr(tmp, "logdet")
            attr(tmp, "logdet") <- NULL
            icov[[g]]        <- tmp
        }
    }

    # cov.vecs
    if(estimator == "GLS" || estimator == "WLS") {
        for(g in 1:ngroups) {
            cov.vecs[[g]] <- vech(sample@cov[[g]])
        }
    }

    # WLS.V (for GLS and WLS only)
    if(estimator == "GLS") {
        for(g in 1:ngroups) {
            if(meanstructure) {
                V11 <- icov[[g]]
                if(mimic == "Mplus") { # is this a bug in Mplus?
                    V11 <- V11 * sample@nobs[[g]]/(sample@nobs[[g]]-1)
                }
                V22 <- 0.5 * D.pre.post(icov[[g]] %x% icov[[g]])
                WLS.V[[g]] <- bdiag(V11,V22)
            } else {
                WLS.V[[g]] <-
                    0.5 * D.pre.post(icov[[g]] %x% icov[[g]])
            }
        }
    } else if(estimator == "WLS") {
        for(g in 1:ngroups) {
            # sample size large enough?
            pstar <- nvar*(nvar+1)/2
            if(meanstructure) pstar <- pstar + nvar
            if(d.nobs[g] < pstar) {
                if(g > 1L) cat("in group: ", g, ":\n", sep="")
                stop("lavaan ERROR: cannot compute Gamma: number of observations too small")
            }

            Gamma <- compute.Gamma(X[[g]], meanstructure=meanstructure,
                                   Mplus.WLS=(mimic=="Mplus"))
            # Gamma should be po before we invert
            ev <- eigen(Gamma, symmetric=FALSE, only.values=TRUE)$values
            if(is.complex(ev) || any(Re(ev) < 0)) {
               stop("lavaan ERROR: Gamma (weight) matrix is not positive-definite")
            }
            #d.WLS.V[[g]] <- MASS.ginv(Gamma) # can we avoid ginv?
            WLS.V[[g]] <- inv.chol(Gamma)
        }
    }



    # construct Sample object
    SampleStatsExtra <- new("SampleStatsExtra",

                            icov           = icov,
                            cov.log.det    = cov.log.det,
                            cov.vecs       = cov.vecs,
                            WLS.V          = WLS.V
                           )

    SampleStatsExtra
}

