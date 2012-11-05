# constructor for the 'lavData' class
#
# the lavData class describes how the data looks like
#  - do we have a full data frame, or only sample statistics
#  - variable type ("numeric", "ordered", ...)
#  - how many groups, how many observations, ...
#  - what about missing patterns?
#
# initial version: YR 14 April 2012

# extract the data we need for this particular model
lavData <- function(data          = NULL,          # data.frame
                    group         = NULL,          # multiple groups?
                    group.label   = NULL,          # custom group labels?
                    ov.names      = names(data),   # variables needed in model
                    ordered       = NULL,          # ordered variables
                    ov.names.x    = character(0),  # exo variables
                    std.ov        = FALSE,         # standardize ov's?
                    missing       = "listwise",    # remove missings?
                    sample.cov    = NULL,          # sample covariance(s)
                    sample.mean   = NULL,          # sample mean vector(s)
                    sample.nobs   = NULL,          # sample nobs
                    warn          = TRUE           # produce warnings?
                   ) 
{

    # three scenarios:
    #    1) data is full data.frame
    #    2) data are sample statistics only
    #    3) no data at all

    # 1) full data
    if(!is.null(data)) {
        stopifnot(is.data.frame(data)) ## FIXME!! we should also allow matrices
        lavData <- getDataFull(data        = data,
                               group       = group,
                               group.label = group.label,
                               ov.names    = ov.names,
                               ordered     = ordered,
                               ov.names.x  = ov.names.x,
                               std.ov      = std.ov,
                               missing     = missing,
                               warn        = warn)
    } 
    
    
    # 2) sample moments
    else if(!is.null(sample.cov)) {
        
        # we also need the number of observations (per group)
        if(is.null(sample.nobs))
            stop("lavaan ERROR: please specify number of observations")

        # if meanstructure=TRUE, we need sample.mean
        #if(meanstructure == TRUE && is.null(sample.mean))
        #    stop("lavaan ERROR: please provide sample.mean if meanstructure=TRUE")
        # if group.equal contains "intercepts", we need sample.mean
        #if("intercepts" %in% group.equal && is.null(sample.mean))
        #    stop("lavaan ERROR: please provide sample.mean if group.equal contains \"intercepts\"")

        # list?
        if(is.list(sample.cov)) {
            # multiple groups, multiple cov matrices
            if(!is.null(sample.mean)) {
                stopifnot(length(sample.mean) == length(sample.cov))
            }
            # multiple groups, multiple cov matrices
            ngroups     <- length(sample.cov)
            LABEL <- names(sample.cov)
            if(is.null(group.label)) {
                if(is.null(LABEL))
                    group.label <- paste("Group ", 1:ngroups, sep="")
                else
                    group.label <- LABEL
            } else {
                if(is.null(LABEL)) {
                    stopifnot(length(group.label) == ngroups)
                } else {
                    # FIXME!!!!
                    # check if they match
                }   
            }
        } else {
            ngroups <- 1L; group.label <- character(0)
            if(!is.matrix(sample.cov))
                stop("lavaan ERROR: sample.cov must be a matrix or a list of matrices")
        }

        # handle ov.names
        if(!is.list(ov.names)) {
            tmp <- ov.names; ov.names <- vector("list", length=ngroups)
            ov.names[1:ngroups] <- list(tmp)
        } else {
            if(length(ov.names) != ngroups)
                stop("lavaan ERROR: ov.names assumes ", length(ov.names),
                     " groups; data contains ", ngroups, " groups")
        }
        # handle ov.names.x
        if(!is.list(ov.names.x)) {
            tmp <- ov.names.x; ov.names.x <- vector("list", length=ngroups)
            ov.names.x[1:ngroups] <- list(tmp)
        } else {
            if(length(ov.names.x) != ngroups)
                stop("lavaan ERROR: ov.names.x assumes ", length(ov.names.x),
                     " groups; data contains ", ngroups, " groups")
        }

        ov <- list()
        ov$name <- unique(unlist(c(ov.names,ov.names.x)))
        nvar    <- length(ov$name)
        ov$idx  <- rep(NA, nvar)
        ov$nobs <- rep(sample.nobs, nvar)
        ov$type <- rep("numeric", nvar)

        # construct lavData object
        lavData <- new("lavData",
                       data.type   = "moment",
                       ngroups     = ngroups, 
                       group       = character(0L),
                       group.label = group.label,
                       nobs        = as.list(sample.nobs),
                       norig       = as.list(sample.nobs),
                       ov.names    = ov.names, 
                       ov.names.x  = ov.names.x,
                       ov          = ov,
                       missing     = "listwise")



    # 3) data.type = "none":  both data and sample.cov are NULL
    } else {
        if(is.null(sample.nobs)) sample.nobs <- 0L
        sample.nobs <- as.list(sample.nobs)
        ngroups <- length(unlist(sample.nobs))
        if(ngroups > 1L)
            group.label <- paste("Group ", 1:ngroups, sep="")
        else
            group.label <- character(0)

        # handle ov.names
        if(!is.list(ov.names)) {
            tmp <- ov.names; ov.names <- vector("list", length=ngroups)
            ov.names[1:ngroups] <- list(tmp)
        }
        # handle ov.names.x
        if(!is.list(ov.names.x)) {
            tmp <- ov.names.x; ov.names.x <- vector("list", length=ngroups)
            ov.names.x[1:ngroups] <- list(tmp)
        }

        ov <- list()
        ov$name <- unique(unlist(c(ov.names,ov.names.x)))
        nvar    <- length(ov$name)
        ov$idx  <- rep(NA, nvar)
        ov$nobs <- rep(0L, nvar)
        ov$type <- rep("numeric", nvar)

        # construct lavData object
        lavData <- new("lavData",
                       data.type   = "none",
                       ngroups     = ngroups,
                       group       = character(0L),
                       group.label = group.label,
                       nobs        = sample.nobs,
                       norig       = sample.nobs,
                       ov.names    = ov.names, 
                       ov.names.x  = ov.names.x,
                       ov          = ov,
                       missing     = "listwise")
    }

    lavData
}


# handle full data
getDataFull <- function(data          = NULL,          # data.frame
                        group         = NULL,          # multiple groups?
                        group.label   = NULL,          # custom group labels?
                        ov.names      = names(data),   # variables needed in model
                        ordered       = NULL,          # ordered variables
                        ov.names.x    = character(0),  # exo variables
                        std.ov        = FALSE,         # standardize ov's?
                        missing       = "listwise",    # remove missings?
                        warn          = TRUE           # produce warnings?
                       )
{
    # number of groups and group labels
    if(!is.null(group) && length(group) > 0L) {
        if(!(group %in% names(data))) {
            stop("lavaan ERROR: grouping variable ", sQuote(group),
                 " not found;\n  ",
                 "variable names found in data frame are:\n  ", 
                 paste(names(data), collapse=" "))
        }
        # note: by default, we use the order as in the data; 
        # not as in levels(data[,group])
        if(length(group.label) == 0L) {
            group.label <- unique(as.character(data[,group]))
            if(warn && any(is.na(group.label))) {
                warning("lavaan WARNING: group variable ", sQuote(group), 
                        " contains missing values\n", sep="")
            }
            group.label <- group.label[!is.na(group.label)]
        } else {
            group.label <- unique(as.character(group.label))
            # check if user-provided group labels exist
            LABEL <- unique(as.character(data[,group]))
            idx <- match(group.label, LABEL)
            if(warn && any(is.na(idx))) {
                warning("lavaan WARNING: some group.labels do not appear ",
                        "in the grouping variable: ",  
                        paste(group.label[which(is.na(idx))], collapse=" "))
            }
            group.label <- group.label[!is.na(idx)]
            # any groups left?
            if(length(group.label) == 0L)
                stop("lavaan ERROR: no group levels left; check the group.label argument")
        }
        ngroups     <- length(group.label)
    } else {
        if(warn && length(group.label) > 0L)
            warning("lavaan WARNING: `group.label' argument",
                    " will be ignored if `group' argument is missing")
        ngroups <- 1L
        group.label <- character(0L)
        group <- character(0L)
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
        if(is.list(ov.names.x)) {
            if(length(ov.names.x) != ngroups)
                stop("lavaan ERROR: ov.names assumes ", length(ov.names.x),
                     " groups; data contains ", ngroups, " groups")
        } else {
            tmp <- ov.names.x
            ov.names.x <- vector("list", length=ngroups)
            ov.names.x[1:ngroups] <- list(tmp)
        }
    } else {
        if(is.list(ov.names)) {
            if(length(ov.names) > 1L)
                stop("lavaan ERROR: model syntax defines multiple groups; data suggests a single group")
        } else {
            ov.names <- list(ov.names)
        }
        if(is.list(ov.names.x)) {
            if(length(ov.names.x) > 1L)
                stop("lavaan ERROR: model syntax defines multiple groups; data suggests a single group")
        } else {
            ov.names.x <- list(ov.names.x)
        }
    }

    # construct OV list -- FIXME: surely, this can be done more elegantly??
    for(g in 1:ngroups) {
        # does the data contain all the observed variables
        # needed in the user-specified model for this group
        idx.missing <- which(!(ov.names[[g]] %in% names(data)))
        if(length(idx.missing)) {
            stop("lavaan ERROR: missing observed variables in dataset: ",
                 paste(ov.names[[g]][idx.missing], collapse=" "))
        }
    }

    # here, we now for sure all ov.names exist in the data.frame
    # create varTable
    ov <- varTable(data, ov.names = ov.names, ov.names.x = ov.names.x, 
                   as.data.frame. = FALSE)

    # do some checking
    # check for unordered factors
    if("factor" %in%  ov$type) {
        f.names <- ov$name[ov$type == "factor"]
        if(any(f.names %in% unlist(ov.names)))
            warning(paste("lavaan WARNING: unordered factor(s) detected in data:", paste(f.names, collapse=" ")))
    }
    # check for zero-cases
    idx <- which(ov$nobs == 0L || ov$var == 0)
    if(length(idx) > 0L) {
        OV <- as.data.frame(ov)
        rn <- rownames(OV)
        rn[idx] <- paste(rn[idx], "***", sep="")
        rownames(OV) <- rn
        print(OV)
        stop("lavaan ERROR: some variables have no values (only missings) or no variance")
    }

    # prepare empty list for data.matrix per group
    case.idx <- vector("list", length=ngroups)
    nobs     <- vector("list", length=ngroups)
    norig    <- vector("list", length=ngroups)
    Mp       <- vector("list", length=ngroups)
    X        <- vector("list", length=ngroups)
    eXo      <- vector("list", length=ngroups)

    # for each group
    for(g in 1:ngroups) {

        # extract variables in correct order
        ov.idx  <- ov$idx[match(ov.names[[g]],   ov$name)]
        exo.idx <- ov$idx[match(ov.names.x[[g]], ov$name)] 

        # extract cases per group
        if(ngroups > 1L || length(group.label) > 0L) {
            if(missing == "listwise") {
                case.idx[[g]] <- which(data[, group] == group.label[g] &
                                       complete.cases(data[,ov.idx]))
                nobs[[g]] <- length(case.idx[[g]])
                norig[[g]] <- length(which(data[, group] == group.label[g]))
            } else {
                case.idx[[g]] <- which(data[, group] == group.label[g])
                nobs[[g]] <- norig[[g]] <- length(case.idx[[g]])
            }
        } else {
            if(missing == "listwise") {
                case.idx[[g]] <- which(complete.cases(data[,ov.idx]))
                nobs[[g]] <- length(case.idx[[g]])
                norig[[g]] <- nrow(data)
            } else {
                case.idx[[g]] <- 1:nrow(data)
                nobs[[g]] <- norig[[g]] <- length(case.idx[[g]])
            }
        }

        # extract data
        X[[g]] <- data.matrix( data[case.idx[[g]], ov.idx, drop=FALSE] )
        dimnames(X[[g]]) <- NULL
        if(length(exo.idx) > 0L) {
            eXo[[g]] <- data.matrix( data[case.idx[[g]], exo.idx, drop=FALSE] )
            dimnames(eXo[[g]]) <- NULL
        } else {
            eXo[g] <- list(NULL)
        }
        #print( tracemem(X[[g]]) )

        # standardize observed variables?
        if(std.ov) {
            X[[g]]  <- scale(X[[g]])[,] # three copies are made!
            if(length(exo.idx) > 0L)
                eXo[[g]] <- scale(eXo[[g]])[,]
        }

        # missing data
        if(missing != "listwise") {
            # get missing patterns
            Mp[[g]] <- getMissingPatterns(X[[g]])
            # checking!
            if(length(Mp[[g]]$empty.idx) > 0L) {
                X[[g]] <- X[[g]][-Mp[[g]]$empty.idx,,drop=FALSE]
                warning("lavaan WARNING: some cases are empty and will be removed:\n  ", paste(Mp[[g]]$empty.idx, collapse=" "))
            }
            if(any(Mp[[g]]$coverage < 0.1)) {
                warning("lavaan WARNING: due to missing values, some pairwise combinations have less than 10% coverage")
            }
            # in case we had observations with only missings
            nobs[[g]] <- Mp[[g]]$nobs
        }

        # warn if we have a small number of observations (but NO error!)
        if( nobs[[g]] < (nvar <- length(ov.idx)) ) {
            txt <- ""
            if(ngroups > 1L) txt <- paste(" in group ", g, sep="")
            warning("lavaan WARNING: small number of observations (nobs < nvar)", txt,
                    "\n  nobs = ", nobs[[g]], " nvar = ", nvar)
        }

    } # ngroups

    lavData <- new("lavData",
                      data.type       = "full",
                      ngroups         = ngroups,
                      group           = group,
                      group.label     = group.label,
                      std.ov          = std.ov,
                      nobs            = nobs,
                      norig           = norig,
                      ov.names        = ov.names,
                      ov.names.x      = ov.names.x,
                      #ov.types        = ov.types,
                      #ov.idx          = ov.idx,
                      ov              = ov,
                      case.idx        = case.idx,
                      missing         = missing,
                      X               = X,
                      eXo             = eXo,
                      Mp              = Mp
                     )
    lavData                     
}
