# constructor for the 'Sample' class
#
# initial version: YR 25/03/2009

Sample <- function(data=NULL,
                   group=NULL,
                   sample.cov=NULL,
                   sample.mean=NULL,
                   sample.nobs=NULL,
                   std.ov=FALSE,

                   ov.names=character(0),
                   data.type="unknown",
                   ngroups=1L,
                   group.label=character(0),
                   estimator="ML",
                   likelihood="normal",
                   mimic="Mplus",
                   meanstructure=FALSE,
                   missing="listwise",

                   warn=TRUE,
                   verbose=FALSE)
{

    # number of observed variables
    nvar  <- length(ov.names)

    # set missing flag (can be overriden later)
    if(missing == "ml") {
        missing.flag <- rep(TRUE, ngroups)
    } else {
        missing.flag <- rep(FALSE, ngroups)
    }

    # sample statistics per group
    d.cov         <- vector("list", length=ngroups)
    d.icov        <- vector("list", length=ngroups)
    d.cov.log.det <- vector("list", length=ngroups)
    d.cov.vecs    <- vector("list", length=ngroups)
    d.var     <- vector("list", length=ngroups)
    d.mean    <- vector("list", length=ngroups)
    d.nobs    <- vector("list", length=ngroups)
    d.missing <- vector("list", length=ngroups)
    d.data    <- vector("list", length=ngroups)
    d.WLS.V   <- vector("list", length=ngroups)

    #### FULL DATA FRAME ####
    if(data.type == "full") {

        # does the data frame contain all the observed variables
        # needed in the user-specified model?
        idx.missing <- which(!(ov.names %in% names(data)))
        if(length(idx.missing)) {
            stop("missing observed variables in dataset: ",
                 paste(ov.names[idx.missing], collapse=" "))
        }
     
        for(g in 1:ngroups) {

            # extract variables in correct order
            var.idx <- match(ov.names, names(data))
            if(ngroups > 1L) {
                case.idx <- data[, group] == group.label[g]
                data.obs <- data[case.idx, var.idx]
            } else {
                data.obs <- data[, var.idx]
            }

            # check if we have enough observations
            if(nrow(data.obs) < nvar)
                stop("lavaan ERROR: too few observations (nobs < nvar)")

            # data.obs should contain numeric values only
            # strip dimnames and coerce to matrix
            data.obs <- data.matrix(data.obs); dimnames(data.obs) <- NULL

            # standardize observed variables?
            if(std.ov) {
                data.obs <- scale(data.obs)[,]
            }

            # missing data?
            d.missing[[g]] <- list()
            if(!missing.flag[g]) {
                # store original nobs per group
                d.missing[[g]]$norig <- nrow(data.obs)
                data.obs <- na.omit(data.obs)
                d.missing[[g]]$nobs <- nrow(data.obs)
                # check again if we have enough observations
                if(nrow(data.obs) == 0) 
                    stop("lavaan ERROR: no cases left after listwise deletion")
                if(nrow(data.obs) < nvar)
                    stop("lavaan ERROR: too few observations (nobs < nvar)")
            } else {
                # do more (but only if missing!)
                #   - get missing patterns
                #   - store sufficient statistics (per missing pattern group)
                #   - compute pairwise coverage
                d.missing[[g]] <- missing.patterns(data.obs, warn=warn)
                if(d.missing[[g]]$npatterns > 1L) {
                    # ok, we have missing values - check the estimator
                    if(estimator == "GLS" || estimator == "WLS") {
                        stop("estimator ", estimator, 
                        " not allowed when data contains missing values\n")
                    }
                    # estimate moments
                    out <- estimate.moments.fiml(X=data.obs, M=d.missing[[g]],
                                                 verbose=verbose)
                    d.missing[[g]]$sigma <- out$sigma
                    d.missing[[g]]$mu    <- out$mu
                    d.missing[[g]]$h1    <- out$fx
                } else {
                    # data is complete after all (for this group)
                    missing.flag[g] <- FALSE
                    d.missing[[g]]$norig <- nrow(data.obs)
                    d.missing[[g]]$nobs  <- nrow(data.obs)
                }
            }
  
            # fill in the other slots
            d.cov[[g]] <- cov(data.obs, use="pairwise") # must be pairwise
            d.var[[g]]  <- apply(data.obs, 2, var, na.rm=TRUE)
            d.mean[[g]] <- apply(data.obs, 2, mean, na.rm=TRUE)
            d.nobs[[g]] <- d.missing[[g]]$nobs

            # for convenience (and for predict and residual), we store
            # a local copy of the data inside
            # FIXME!!
            # for large datasets, this may be overkill!
            if(!missing.flag[g]) {
                d.data[[g]] <- data.obs
            }
        } # ngroups

    } 

    if(data.type == "moment") {

        sample.nobs <- as.list(as.integer(sample.nobs))
        if(!is.list(sample.cov)) sample.cov  <- list(sample.cov)
        if(!is.null(sample.mean) && !is.list(sample.mean)) 
            sample.mean <- list(sample.mean)

        for(g in 1:ngroups) {

            tmp.cov <- sample.cov[[g]]

            # is it a not a correlation matrix?
            if(g == 1 && ( sum(diag(tmp.cov)) == ncol(tmp.cov) ) ) {
                # ok, we do not stop here, but spit out a big fat warning
                text <- 
paste("  \nsample covariance matrix looks like a correlation matrix!\n\n",
      "  lavaan currently does not support the analysis of correlation\n",
      "  matrices; the standard errors in the summary output will be most \n",
      "  likely wrong; see the following reference:\n\n",
      "  Cudeck, R. (1989). Analysis of correlation matrices using covariance\n",
      "  structure models. Psychological Bulletin, 105, 317-327.\n\n", sep="")
                warning(text)
            }
 
            # lower triangle must be filled
            if(all(tmp.cov[lower.tri(tmp.cov)] == 0)) {
                stop("please provide a lower-triangular covariance matrix!\n")
            }
            # make sure that the matrix is fully symmetric
            T <- t(tmp.cov)
            tmp.cov[upper.tri(tmp.cov)] <- T[upper.tri(T)]

            # check dimnames
            if(is.null(rownames(tmp.cov))) {
                stop("please provide row names for the covariance matrix!\n")
            }
   
            # extract only the part we need (using ov.names)
            idx <- match(ov.names, rownames(tmp.cov))
            if(any(is.na(idx))) {
                cat("found: ", rownames(tmp.cov)[idx], "\n")
                cat("expected: ", ov.names, "\n")
                stop("rownames of covariance matrix do not match the model!\n")
            } else {
                tmp.cov <- tmp.cov[idx,idx]
            }

            # strip dimnames
            dimnames(tmp.cov) <- NULL
        
            if(is.null(sample.mean)) {
                # assume zero mean vector
                tmp.mean <- numeric(nvar)
            } else {
                # extract only the part we need
                tmp.mean <- sample.mean[[g]][idx]
                names(tmp.mean) <- NULL
            }

            d.cov[[g]]  <- tmp.cov
            d.var[[g]]  <- diag(tmp.cov)
            d.mean[[g]] <- tmp.mean
            d.nobs[[g]] <- sample.nobs[[g]]
            d.missing[[g]] <- list(norig=sample.nobs[[g]])
        } # g
    } # moment

    # rescale d.cov? only if ML and likelihood == "normal"
    if((estimator == "ML") && likelihood == "normal") {
        for(g in 1:ngroups) {
            # we 'transform' the sample cov (divided by n-1) 
            # to a sample cov divided by 'n'
            d.cov[[g]] <- (d.nobs[[g]]-1)/d.nobs[[g]] * d.cov[[g]]
        }
    }

    # icov and cov.log.det
    for(g in 1:ngroups) {
        tmp <- try(inv.chol(d.cov[[g]], logdet=TRUE))
        if(inherits(tmp, "try-error")) {
            if(ngroups > 1) {
                stop("sample covariance can not be inverted in group", g)
            } else {
                stop("sample covariance can not be inverted")
            }
        } else {
            d.cov.log.det[[g]] <- attr(tmp, "logdet")
            attr(tmp, "logdet") <- NULL
            d.icov[[g]]        <- tmp
        }
    }

    # cov.vecs
    for(g in 1:ngroups) {
        d.cov.vecs[[g]] <- vech(d.cov[[g]])
    }

    # WLS.V (for GLS and WLS only)
    if(estimator == "GLS") {
        for(g in 1:ngroups) {
            if(meanstructure) {
                V11 <- d.icov[[g]]
                if(mimic == "Mplus") { # is this a bug in Mplus?
                    V11 <- V11 * d.nobs[[g]]/(d.nobs[[g]]-1)
                }
                V22 <- 0.5 * D.pre.post(d.icov[[g]] %x% d.icov[[g]])
                d.WLS.V[[g]] <- bdiag(V11,V22)
            } else {
                d.WLS.V[[g]] <-
                    0.5 * D.pre.post(d.icov[[g]] %x% d.icov[[g]])
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

            Gamma <- compute.Gamma(d.data[[g]], meanstructure=meanstructure,
                                   Mplus.WLS=(mimic=="Mplus"))
            # Gamma should be po before we invert
            ev <- eigen(Gamma, symmetric=FALSE, only.values=TRUE)$values
            if(is.complex(ev) || any(Re(ev) < 0)) {
                stop("lavaan ERROR: Gamma (weight) matrix is not positive-definite")
            }
            #d.WLS.V[[g]] <- MASS.ginv(Gamma) # can we avoid ginv?
            d.WLS.V[[g]] <- inv.chol(Gamma)
        }
    }
   
    # construct Sample object
    Sample <- new("Sample",
                  ov.names=ov.names,
                  nvar=nvar,
                  ngroups=ngroups,
                  group.label=group.label,
                  missing.flag=missing.flag,

                  ntotal=sum(unlist(d.nobs)),
                  nobs=d.nobs,

                  mean=d.mean,
                  cov=d.cov,
                  icov=d.icov,
                  cov.log.det=d.cov.log.det,
                  cov.vecs=d.cov.vecs,
                  WLS.V=d.WLS.V,
                  var=d.var,
                  data.obs=d.data,
                  
                  missing=d.missing

                 )

    Sample
}

