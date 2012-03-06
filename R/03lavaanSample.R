# constructor for the 'Sample' class
#
# initial version: YR 25/03/2009
# major revision: YR 5/11/2011: separate data.obs and sample statistics

# extract the data we need for this particular model
getData <- function(data          = NULL,          # data.frame
                    group         = NULL,          # multiple groups?
                    ov.names      = NULL,          # variables needed in model
                    std.ov        = FALSE,         # standardize ov's?
                    missing       = "listwise",    # remove missings?
                    warn          = TRUE           # produce warnings?
                   ) 
{
    # number of groups
    ngroups <- 1L; group.label <- character(0)
    if(!is.null(group)) {
        if(!(group %in% names(data))) {
            stop("lavaan ERROR: grouping variable ", sQuote(group),
                 " not found;\n  ",
                 "variable names found in data frame are:\n  ", 
                 paste(names(data), collapse=" "))
        }
        # note: we use the order as in the data; not as in levels(data)
        group.label <- unique(as.character(data[,group]))
        if(warn && any(is.na(group.label))) {
            cat("lavaan WARNING: group variable ", sQuote(group), 
                " contains missing values\n", sep="")
        }
        group.label <- group.label[!is.na(group.label)]
        ngroups     <- length(group.label)
    }

    # ov.names
    if(ngroups > 1L) {
        if(is.list(ov.names)) {
            if(length(ov.names) != ngroups)
                stop("lavaan ERROR: ov.names assumes ", length(ov.names),
                     " groups; data contains ", ngroups, " groups")
        } else {
            tmp <- ov.names
            ov.names <- vector("list", length=ngroups)
            ov.names[1:ngroups] <- list(tmp)
        }
    } else {
        if(is.list(ov.names)) {
            if(length(ov.names) > 1L)
                stop("lavaan ERROR: model syntax defines multiple groups; data suggests a single group")
        } else {
            ov.names <- list(ov.names)
        }
    }

    # prepare empty list for data.matrix per group
    ov.idx   <- vector("list", length=ngroups)
    case.idx <- vector("list", length=ngroups)
    nobs     <- vector("list", length=ngroups)
    norig    <- vector("list", length=ngroups)
    X        <- vector("list", length=ngroups)

    # for each group
    for(g in 1:ngroups) {

        # does the data contain all the observed variables
        # needed in the user-specified model for this group
        idx.missing <- which(!(ov.names[[g]] %in% names(data)))
        if(length(idx.missing)) {
            stop("lavaan ERROR: missing observed variables in dataset: ",
                 paste(ov.names[[g]][idx.missing], collapse=" "))
        }

        # extract variables in correct order
        ov.idx[[g]] <- match(ov.names[[g]], names(data))

        # extract cases per group
        if(ngroups > 1L) {
            if(missing == "listwise") {
                case.idx[[g]] <- which(data[, group] == group.label[g] &
                                       complete.cases(data[,ov.idx[[g]]]))
                nobs[[g]] <- length(case.idx[[g]])
                norig[[g]] <- length(which(data[, group] == group.label[g]))
            } else {
                case.idx[[g]] <- which(data[, group] == group.label[g])
                nobs[[g]] <- norig[[g]] <- length(case.idx[[g]])
            }
        } else {
            if(missing == "listwise") {
                case.idx[[g]] <- which(complete.cases(data[,ov.idx[[g]]]))
                nobs[[g]] <- length(case.idx[[g]])
                norig[[g]] <- nrow(data)
            } else {
                case.idx[[g]] <- 1:nrow(data)
                nobs[[g]] <- norig[[g]] <- length(case.idx[[g]])
            }
        }

        # check if we have enough observations
        if( nobs[[g]] < (nvar <- length(ov.idx[[g]])) ) {
            txt <- ""
            if(ngroups > 1L) txt <- paste(" in group ", g, sep="")
            stop("lavaan ERROR: too few observations (nobs < nvar)", txt,
                 "\n  nobs = ", nobs[[g]], " nvar = ", nvar)
        }


        ### ONLY if we wish to store X ####
     
        # extract data
        X[[g]] <- data.matrix( data[case.idx[[g]], ov.idx[[g]]] )
        #print( tracemem(X[[g]]) )
 
        # get rid of row names, but keep column names
        #rownames(X[[g]]) <- NULL ### WHY are two copies made here? 
                                  ### answer: rownames is NOT a primitive
        dimnames(X[[g]]) <- list(NULL, ov.names[[g]]) # only 1 copy

        # standardize observed variables?
        if(std.ov) {
            X[[g]] <- scale(X[[g]])[,] # three copies are made!
        }

    } # ngroups


    lavaanData <- new("lavaanData",
                      ngroups         = ngroups,
                      group.label     = group.label,
                      nobs            = nobs,
                      norig           = norig,
                      ov.names        = ov.names,
                      ov.idx          = ov.idx,
                      case.idx        = case.idx,
                      X               = X,
                      isComplete      = TRUE,
                      missingPatterns = list()
                     )
    lavaanData                     
}

getMissingPatterns <- function(Data    = NULL, 
                               missing = "listwise",
                               warn    = TRUE,
                               verbose = FALSE) {

    # get X
    X <- Data@X

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
            missing[[g]] <- list(npatterns=0L, flag=FALSE, nobs=nobs)
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

getSampleStatsFromData <- function(Data        = NULL,
                                   M           = NULL,
                                   boot.idx    = NULL,
                                   rescale     = FALSE,
                                   WLS.V       = list()) {

    # get X
    X <- Data@X

    # number of groups
    ngroups <- Data@ngroups
    nobs <- Data@nobs
   
    # sample statistics per group
    cov         <- vector("list", length=ngroups)
    var         <- vector("list", length=ngroups)
    mean        <- vector("list", length=ngroups)
    # extra sample statistics per group
    icov        <- vector("list", length=ngroups)
    cov.log.det <- vector("list", length=ngroups)
    cov.vecs    <- vector("list", length=ngroups)

    for(g in 1:ngroups) {

        # bootstrap sample?
        if(!is.null(boot.idx)) {
            X[[g]] <- X[[g]][boot.idx[[g]],,drop=FALSE]
        }

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

        # icov and cov.log.det
        tmp <- try(inv.chol(cov[[g]], logdet=TRUE))
        if(inherits(tmp, "try-error")) {
            if(ngroups > 1) {
                stop("lavaan ERROR: sample covariance can not be inverted in group: ", g)
            } else {
                stop("lavaan ERROR: sample covariance can not be inverted")
            }
        } else {
            cov.log.det[[g]] <- attr(tmp, "logdet")
            attr(tmp, "logdet") <- NULL
            icov[[g]]        <- tmp
        }

        # cov.vecs
        cov.vecs[[g]] <- vech(cov[[g]])

    } # ngroups

    # construct SampleStats object
    SampleStats <- new("SampleStats",

                       # sample moments
                       mean        = mean,
                       cov         = cov,
                       #var        = var,

                       # convenience
                       nobs        = nobs,
                       ntotal      = sum(unlist(nobs)),
                       ngroups     = ngroups,

                       # extra sample statistics
                       icov        = icov,
                       cov.log.det = cov.log.det,
                       cov.vecs    = cov.vecs,
                       WLS.V       = WLS.V,                     

                       # missingness
                       missing     = M
                      )

    SampleStats
}


getSampleStatsFromMoments <- function(sample.cov  = NULL,
                                      sample.mean = NULL,
                                      sample.nobs = NULL,
                                      rescale     = FALSE,
                                      ov.names    = NULL,
                                      WLS.V       = list()) {

    # matrix -> list
    if(!is.list(sample.cov)) sample.cov  <- list(sample.cov)
        if(!is.null(sample.mean) && !is.list(sample.mean))
            sample.mean <- list(sample.mean)

    # number of groups
    ngroups <- length(sample.cov)
   
    # sample statistics per group
    cov         <- vector("list", length=ngroups)
    #var         <- vector("list", length=ngroups)
    mean        <- vector("list", length=ngroups)
    nobs        <- as.list(as.integer(sample.nobs))
    # extra sample statistics per group
    icov        <- vector("list", length=ngroups)
    cov.log.det <- vector("list", length=ngroups)
    cov.vecs    <- vector("list", length=ngroups)
    
    # prepare empty list for missing data
    missing <- vector("list", length=ngroups)

    for(g in 1:ngroups) {

        tmp.cov <- sample.cov[[g]]

        # make sure that the matrix is fully symmetric (NEEDED?)
        T <- t(tmp.cov)
        tmp.cov[upper.tri(tmp.cov)] <- T[upper.tri(T)]

        # check dimnames
        if(!is.null(rownames(tmp.cov))) {
            cov.names <- rownames(tmp.cov)
        } else if(!is.null(colnames(tmp.cov))) {
            cov.names <- colnames(tmp.cov)
        } else {
            stop("lavaan ERROR: please provide row/col names ",
                "for the covariance matrix!\n")
        }

        # extract only the part we need (using ov.names)
        idx <- match(ov.names[[g]], cov.names)
        if(any(is.na(idx))) {
            cat("found: ", cov.names, "\n")
            cat("expected: ", ov.names[[g]], "\n")
            stop("lavaan ERROR: rownames of covariance matrix do not match ",
                 "the model!\n", 
                 "  found: ", paste(cov.names, collapse=" "), "\n",
                 "  expected: ", paste(ov.names[[g]], collapse=" "), "\n")
        } else {
            tmp.cov <- tmp.cov[idx,idx]
        }

        # strip dimnames
        dimnames(tmp.cov) <- NULL

        if(is.null(sample.mean)) {
            # assume zero mean vector
            tmp.mean <- numeric(ncol(tmp.cov))
        } else {
            # extract only the part we need
            tmp.mean <- sample.mean[[g]][idx]
            names(tmp.mean) <- NULL
        }

        cov[[g]]  <- tmp.cov
        #var[[g]]  <- diag(tmp.cov)
        mean[[g]] <- tmp.mean

        # rescale cov by (N-1)/N?
        if(rescale) {
            # we 'transform' the sample cov (divided by n-1) 
            # to a sample cov divided by 'n'
            cov[[g]] <- (nobs[[g]]-1)/nobs[[g]] * cov[[g]]
        }

        # icov and cov.log.det
        tmp <- try(inv.chol(cov[[g]], logdet=TRUE))
        if(inherits(tmp, "try-error")) {
            if(ngroups > 1) {
                stop("lavaan ERROR: sample covariance can not be inverted in group: ", g)
            } else {
                stop("lavaan ERROR: sample covariance can not be inverted")
            }
        } else {
            cov.log.det[[g]] <- attr(tmp, "logdet")
            attr(tmp, "logdet") <- NULL
            icov[[g]]        <- tmp
        }

        # cov.vecs
        cov.vecs[[g]] <- vech(cov[[g]])

        # missing
        missing[[g]] <- list(npatterns=0L, flag=FALSE)

    } # ngroups

    # construct SampleStats object
    SampleStats <- new("SampleStats",

                       # sample moments
                       mean        = mean,
                       cov         = cov,
                       #var        = var,

                       # convenience
                       nobs        = nobs,
                       ntotal      = sum(unlist(nobs)),
                       ngroups     = ngroups,

                       # extra sample statistics
                       icov        = icov,
                       cov.log.det = cov.log.det,
                       cov.vecs    = cov.vecs,
                       WLS.V       = WLS.V,                     

                       # missingness
                       missing     = missing
                      )

    SampleStats
}

getWLS.V <- function(Data          = NULL, 
                     sample        = NULL,
                     boot.idx      = NULL,
                     estimator     = "ML",
                     mimic         = "lavaan",
                     meanstructure = FALSE) {

    # number of groups
    if(is.null(sample)) {
        X <- Data@X
        ngroups <- length(X)
    } else {
        ngroups <- sample@ngroups
    }

    # bootstrap sample?
    if(!is.null(boot.idx)) {
        for(g in 1:ngroups) {
            X[[g]] <- X[[g]][boot.idx[[g]],,drop=FALSE]
        }
    }

    WLS.V       <- vector("list", length=ngroups)

    # WLS.V (for GLS and WLS only)
    if(estimator == "GLS") {
        if(is.null(sample)) {
            # FIXME: maybe we should avoid sample = NULL alltogether...
            # tmp cov and icov, assuming data is complete
            cov         <- lapply(X, cov, use="pairwise")
            icov        <- lapply(cov, inv.chol, logdet=FALSE)
            nobs        <- lapply(X, nrow)
        } else {
            cov <- sample@cov
            icov <- sample@icov
            nobs <- sample@nobs
        }
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
            nvar <- ncol(X[[g]])
            pstar <- nvar*(nvar+1)/2
            if(meanstructure) pstar <- pstar + nvar
            if(nrow(X[[g]]) < pstar) {
                if(ngroups > 1L) { 
                    txt <- cat(" in group: ", g, "\n", sep="")
                } else {
                    txt <- "\n"
                }
                stop("lavaan ERROR: cannot compute Gamma: ",
                     "number of observations (", nrow(X[[g]]), ") too small",
                     txt)
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

    WLS.V
}

