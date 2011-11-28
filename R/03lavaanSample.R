# constructor for the 'Sample' class
#
# initial version: YR 25/03/2009
# major revision: YR 5/11/2011: separate data.obs and sample statistics

# extract the data we need for this particular model
getData <- function(data        = NULL, 
                    ov.names    = character(0),

                    # transform? (for bollen.stine bootstrap)
                    model.cov   = NULL,
                    model.mean  = NULL,
                    
                    # standardize?
                    std.ov      = FALSE,

                    # multiple groups?
                    group       = NULL,
                    ngroups     = 1L,
                    group.label = character(0),
    
                    # how to deal with missing data?
                    missing     = "listwise",

                    debug       = FALSE,
                    warn        = TRUE,
                    verbose     = FALSE
                   ) 
{
    # set missing flag (can be overriden later)
    if(missing == "ml") {
        missing.flag <- rep(TRUE, ngroups)
    } else {
        missing.flag <- rep(FALSE, ngroups)
    }

    # number of observed variables
    nvar  <- length(ov.names)

    # prepare empty list for complete data
    X <- vector("list", length=ngroups)
    norig <- integer(ngroups)
    nobs  <- integer(ngroups)

    # prepare empty list for missing data
    d.missing <- vector("list", length=ngroups)

    # does the data contain all the observed variables
    # needed in the user-specified model for this group
    idx.missing <- which(!(ov.names %in% names(data)))
    if(length(idx.missing)) {
        stop("lavaan ERROR: missing observed variables in dataset: ",
             paste(ov.names[idx.missing], collapse=" "))
    }

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
        # strip dimnames and coerce to matrix
        data.tmp <- data.matrix(data.tmp); dimnames(data.tmp) <- NULL

        # transform observed variables (eg. for bollen-stine bootstrap)
        # simple transform is no missing data!
        if(!is.null(model.cov) && !missing.flag[g]) {
            cov.g <- cov(data.tmp, use="pairwise") # full dataset
            sigma.sqrt <- sqrtSymmetricMatrix( model.cov[[g]] )
            S.inv.sqrt <- sqrtSymmetricMatrix( solve(cov.g) )

            # center
            data.tmp <- scale(data.tmp, center=TRUE, scale=FALSE)[,]

            # transform
            data.tmp <- data.tmp %*% S.inv.sqrt %*% sigma.sqrt

            # add model.mean[[g]]
            if(!is.null(model.mean)) {
                data.tmp <- scale(data.tmp, center=(-1*model.mean[[g]]), 
                                  scale=FALSE)[,]
            }
        }

        # standardize observed variables?
        if(std.ov) {
            data.tmp <- scale(data.tmp)[,]
        }

        # number of observations
        norig[g] <- nobs[g] <- nrow(data.tmp)

        # missing data?
        d.missing[[g]] <- list()
        if(!missing.flag[g]) {
            # listwise deletion
            keep.idx <- complete.cases(data.tmp)
            if(length(keep.idx) > 0L) {
                data.tmp <- data.tmp[keep.idx,,drop=FALSE]
                nobs[g] <- nrow(data.tmp)
            }
            # check again if we have enough observations
            if(nrow(data.tmp) == 0L)
                stop("lavaan ERROR: no cases left after listwise deletion")
            if(nrow(data.tmp) < nvar)
                stop("lavaan ERROR: too few observations (nobs < nvar)")
            # fill in some basic information in d.missing
            d.missing[[g]] <- list(npatterns=0L)
        } else {
            # do more (but only if missing!)
            #   - get missing patterns
            #   - store sufficient statistics (per missing pattern group)
            #   - compute pairwise coverage
            d.missing[[g]] <- missing.patterns(data.tmp, warn=warn)
            if(d.missing[[g]]$npatterns > 1L) {
                # estimate moments
                #out <- estimate.moments.fiml(X=data.obs, M=d.missing[[g]],
                #                             verbose=verbose)
                out <- estimate.moments.EM(X=data.tmp, M=d.missing[[g]],
                                           verbose=verbose)
                d.missing[[g]]$sigma <- out$sigma
                d.missing[[g]]$mu    <- out$mu
                d.missing[[g]]$h1    <- out$fx
            } else {
                # data is complete after all (for this group)
                missing.flag[g] <- FALSE
            }
        }

        X[[g]] <- data.tmp

    } # ngroups

    DataObject <- new("FullData",
                      X           = X,
                      ov.names    = ov.names,
                      nvar        = nvar,
                      ngroups     = ngroups,
                      group.label = group.label,
                      norig       = norig,
                      nobs        = nobs,
                      Missing     = d.missing,
                      missing.flag = missing.flag
                     )

    DataObject
}

getSampleStatsFromDataObject <- function(data          = NULL, 
                                         boot.idx      = NULL,

                                         estimator     = "ML",
                                         likelihood    = "normal",
                                         mimic         = "lavaan",
                                         meanstructure = FALSE,

                                         debug         = FALSE,
                                         warn          = TRUE,
                                         verbose       = FALSE)

{

    # number of variables
    nvar  <- ncol(data@X[[1L]])
    ngroups <- length(data@X)

    # sample statistics per group
    cov         <- vector("list", length=ngroups)
    icov        <- vector("list", length=ngroups)
    cov.log.det <- vector("list", length=ngroups)
    cov.vecs    <- vector("list", length=ngroups)
    var         <- vector("list", length=ngroups)
    mean        <- vector("list", length=ngroups)
    nobs        <- vector("list", length=ngroups)
    missing     <- vector("list", length=ngroups)
    WLS.V       <- vector("list", length=ngroups)

    for(g in 1:ngroups) {

        # bootstrap sample?
        #if(!is.null(boot.idx)) {
        #    if(ngroups == 1L) {
        #        data.obs <- data.obs[boot.idx,]
        #    } else {
        #        CASE.idx <- which(case.idx)
        #        gboot.idx <- boot.idx[which(boot.idx %in% CASE.idx)] 
        #        in.idx <- match(gboot.idx, CASE.idx)
        #        data.obs <- data.obs[in.idx,]
        #    }
        #}


        # fill in the other slots
        cov[[g]]  <-   cov(data@X[[g]], use="pairwise") # must be pairwise
        var[[g]]  <- apply(data@X[[g]], 2,  var, na.rm=TRUE)
        mean[[g]] <- apply(data@X[[g]], 2, mean, na.rm=TRUE)
        nobs[[g]] <- data@nobs[g]

    } # ngroups

    # rescale d.cov? only if ML and likelihood == "normal"
    if((estimator == "ML") && likelihood == "normal") {
        for(g in 1:ngroups) {
            # we 'transform' the sample cov (divided by n-1) 
            # to a sample cov divided by 'n'
            cov[[g]] <- (nobs[[g]]-1)/nobs[[g]] * cov[[g]]
        }
    }

    # icov and cov.log.det
    for(g in 1:ngroups) {
        tmp <- try(inv.chol(cov[[g]], logdet=TRUE))
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
    for(g in 1:ngroups) {
        cov.vecs[[g]] <- vech(cov[[g]])
    }

    # WLS.V (for GLS and WLS only)
    if(estimator == "GLS") {
        for(g in 1:ngroups) {
            if(meanstructure) {
                V11 <- icov[[g]]
                if(mimic == "Mplus") { # is this a bug in Mplus?
                    V11 <- V11 * nobs[[g]]/(nobs[[g]]-1)
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

            Gamma <- compute.Gamma(data@X[[g]], meanstructure=meanstructure,
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
    SampleStats <- new("SampleStats",

                      # sample moments
                      mean           = mean,
                      cov            = cov,
                      var            = var,
                      nobs           = nobs,
                      nvar           = nvar,
                      ntotal         = sum(unlist(nobs)),

                      # convenience
                      ov.names       = data@ov.names,
                      ngroups        = data@ngroups,
                      group.label    = data@group.label,

                      # missing data information
                      missing        = data@Missing,
                      missing.flag   = data@missing.flag,

                      # extra
                      icov           = icov,
                      cov.log.det    = cov.log.det,
                      cov.vecs       = cov.vecs,
                      WLS.V          = WLS.V
                   )
}

