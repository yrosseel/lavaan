# constructor for the 'lavData' class
#
# the lavData class describes how the data looks like
#  - do we have a full data frame, or only sample statistics?
#    (TODO: allow for patterns + freq, if data is categorical)
#  - variable type ("numeric", "ordered", ...)
#  - how many groups, how many observations, ...
#  - what about missing patterns?
#
# initial version: YR 14 April 2012

# YR 23 feb 2017: blocks/levels/groups, but everything is group-based!

# FIXME: if nlevels > 1L, and ngroups > 1L, we should check that
# group is at the upper-level

# extract the data we need for this particular model
lavData <- function(data              = NULL,          # data.frame
                    group             = NULL,          # multiple groups?
                    cluster           = NULL,          # clusters?
                    ov.names          = NULL,          # variables in model
                    ov.names.x        = character(0),  # exo variables
                    ov.names.l        = list(),        # names per level
                    ordered           = NULL,          # ordered variables
                    sampling.weights  = NULL,          # sampling weights
                    sample.cov        = NULL,          # sample covariance(s)
                    sample.mean       = NULL,          # sample mean vector(s)
                    sample.th         = NULL,          # sample thresholds
                    sample.nobs       = NULL,          # sample nobs

                    lavoptions        = lavOptions(),  # lavoptions
                    allow.single.case = FALSE          # for newdata in predict
                   )
{

    # get info from lavoptions

    # group.labels
    group.label <- lavoptions$group.label
    if(is.null(group.label)) {
        group.label <- character(0L)
    }

    # level.labels
    level.label <- lavoptions$level.label
    if(is.null(level.label)) {
        level.label <- character(0L)
    }

    # std.ov?
    std.ov <- lavoptions$std.ov
    if(is.null(std.ov)) {
        std.ov <- FALSE
    }

    # missing?
    missing <- lavoptions$missing
    if(is.null(missing) || missing == "default") {
        missing <- "listwise"
    }

    # warn?
    warn <- lavoptions$warn
    if(is.null(warn)) {
        warn <- TRUE
    }

    # four scenarios:
    #    0) data is already a lavData object: do nothing
    #    1) data is full data.frame (or a matrix)
    #    2) data are sample statistics only
    #    3) no data at all

    # 1) full data
    if(!is.null(data)) {

       # catch lavaan/lavData objects
        if(inherits(data, "lavData")) {
            return(data)
        } else if(inherits(data, "lavaan")) {
            return(data@Data)
        }

        # catch matrix
        if(!is.data.frame(data)) {
            # is it a matrix?
            if(is.matrix(data)) {
                if(nrow(data) == ncol(data)) {
                    # perhaps it is a covariance matrix?
                    if(data[2,1] == data[1,2] && warn) { # not perfect...
                        warning("lavaan WARNING: data argument looks like a covariance matrix; please use the sample.cov argument instead")
                    }
                }
                # or perhaps it is a data matrix?
                ### FIXME, we should avoid as.data.frame() and handle
                ### data matrices directly
                data <- as.data.frame(data, stringsAsFactors = FALSE)
            } else {
                stop("lavaan ERROR: data object of class ", class(data))
            }
        }

        # no ov.names?
        if(is.null(ov.names)) {
            ov.names <- names(data)
            # remove group variable, if provided
            if(length(group) > 0L) {
                group.idx <- which(ov.names == group)
                ov.names <- ov.names[-group.idx]
            }
            # remove cluster variable, if provided
            if(length(cluster) > 0L) {
                cluster.idx <- which(ov.names == cluster)
                ov.names <- ov.names[-cluster.idx]
            }
        }

        lavData <- lav_data_full(data              = data,
                                 group             = group,
                                 cluster           = cluster,
                                 group.label       = group.label,
                                 level.label       = level.label,
                                 ov.names          = ov.names,
                                 ordered           = ordered,
                                 sampling.weights  = sampling.weights,
                                 ov.names.x        = ov.names.x,
                                 ov.names.l        = ov.names.l,
                                 std.ov            = std.ov,
                                 missing           = missing,
                                 warn              = warn,
                                 allow.single.case = allow.single.case)
        sample.cov <- NULL # not needed, but just in case
    }


    # 2) sample moments
    if(is.null(data) && !is.null(sample.cov)) {

        # for now: no levels!!
        nlevels <- 1L

        # we also need the number of observations (per group)
        if(is.null(sample.nobs)) {
            stop("lavaan ERROR: please specify number of observations")
        }

        # list?
        if(is.list(sample.cov)) {
            # multiple groups, multiple cov matrices
            if(!is.null(sample.mean)) {
                stopifnot(length(sample.mean) == length(sample.cov))
            }
            if(!is.null(sample.th)) {
                stopifnot(length(sample.th) == length(sample.cov))
            }
            # multiple groups, multiple cov matrices
            ngroups     <- length(sample.cov)
            LABEL <- names(sample.cov)
            if(is.null(group.label) || length(group.label) == 0L) {
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
            sample.cov <- list(sample.cov)
        }

        # get ov.names
        if (is.null(ov.names)) {
            ov.names <- lapply(sample.cov, row.names)
        } else if (!is.list(ov.names)) {
            # duplicate ov.names for each group
            tmp <- ov.names; ov.names <- vector("list", length = ngroups)
            ov.names[1:ngroups] <- list(tmp)
        } else {
            if (length(ov.names) != ngroups)
                stop("lavaan ERROR: ov.names assumes ", length(ov.names),
                     " groups; data contains ", ngroups, " groups")
            # nothing to do
        }

        # handle ov.names.x
        if(!is.list(ov.names.x)) {
            tmp <- ov.names.x; ov.names.x <- vector("list", length = ngroups)
            ov.names.x[1:ngroups] <- list(tmp)
        } else {
            if(length(ov.names.x) != ngroups)
                stop("lavaan ERROR: ov.names.x assumes ", length(ov.names.x),
                     " groups; data contains ", ngroups, " groups")
        }

        ov <- list()
        ov$name <- unique( unlist(c(ov.names, ov.names.x)) )
        nvar    <- length(ov$name)
        ov$idx  <- rep(NA, nvar)
        ov$nobs <- rep(sum(unlist(sample.nobs)), nvar)
        ov$type <- rep("numeric", nvar)
        ov$nlev <- rep(0, nvar)
        # check for categorical
        if(!is.null(sample.th)) {
            th.idx <- attr(sample.th, "th.idx")
            if(is.list(th.idx)) {
                th.idx <- th.idx[[1]] ## FIRST group only (assuming same ths!)
            }
            if(any(th.idx > 0)) {
                TAB <- table(th.idx[th.idx > 0])
                ord.idx <- as.numeric(names(TAB))
                nlev <- as.integer(unname(TAB) + 1)
                ov$type[ord.idx ] <- "ordered"
                ov$nlev[ord.idx ] <- nlev
            }
        }

        # if std.ov = TRUE, give a warning (suggested by Peter Westfall)
        if(std.ov && warn) {
            warning("lavaan WARNING: std.ov argument is ignored if only sample statistics are provided.")
        }

        # construct lavData object
        lavData <- new("lavData",
                       data.type   = "moment",
                       ngroups     = ngroups,
                       group       = character(0L),
                       nlevels     = 1L, # for now
                       cluster     = character(0L),
                       group.label = group.label,
                       level.label = character(0L),
                       nobs        = as.list(sample.nobs),
                       norig       = as.list(sample.nobs),
                       ov.names    = ov.names,
                       ov.names.x  = ov.names.x,
                       ov.names.l  = ov.names.l,
                       ordered     = as.character(ordered),
                       weights     = vector("list", length = ngroups),
                       sampling.weights = character(0L),
                       ov          = ov,
                       std.ov      = FALSE,
                       missing     = "listwise",
                       case.idx    = vector("list", length = ngroups),
                       Mp          = vector("list", length = ngroups),
                       Rp          = vector("list", length = ngroups),
                       Lp          = vector("list", length = ngroups),
                       X           = vector("list", length = ngroups),
                       eXo         = vector("list", length = ngroups)
                      )

    }

    # 3) data.type = "none":  both data and sample.cov are NULL
    if(is.null(data) && is.null(sample.cov)) {

        # clustered/multilevel? --> ov.names.l should be filled in
        if(length(ov.names.l) > 0L) {
            nlevels <- length(ov.names.l[[1]]) # we assume the same number
                                               # of levels in each group!

            # do we have a cluster argument? if not, create one
            if(is.null(cluster)) {
                if(nlevels == 2L) {
                    cluster <- "cluster"
                } else {
                    cluster <- paste0("cluster", seq_len(nlevels - 1L))
                }
            }

            # default level.labels
            if(length(level.label) == 0L) {
                level.label <- c("within", cluster)
            } else {
                # check if length(level.label) = 1 + length(cluster)
                if(length(level.label) != length(cluster) + 1L) {
                    stop("lavaan ERROR: length(level.label) != length(cluster) + 1L")
                }
                # nothing to do
            }
        } else {
            nlevels <- 1L
            cluster <- character(0L)
            level.label <- character(0L)
        }

        # ngroups: ov.names (when group: is used), or sample.nobs
        if(is.null(ov.names)) {
            warning("lavaan WARNING: ov.names is NULL")
            ov.names <- character(0L)
            if(is.null(sample.nobs)) {
                ngroups <- 1L
                sample.nobs <- rep(list(0L), ngroups)
            } else {
                sample.nobs <- as.list(sample.nobs)
                ngroups <- length(sample.nobs)
            }
        } else if(!is.list(ov.names)) {
            if(is.null(sample.nobs)) {
                ngroups <- 1L
                sample.nobs <- rep(list(0L), ngroups)
            } else {
                sample.nobs <- as.list(sample.nobs)
                ngroups <- length(sample.nobs)
            }
            ov.names <- rep(list(ov.names), ngroups)
        } else if(is.list(ov.names)) {
            ngroups <- length(ov.names)
            if(is.null(sample.nobs)) {
                sample.nobs <- rep(list(0L), ngroups)
            } else {
                sample.nobs <- as.list(sample.nobs)
                if(length(sample.nobs) != ngroups) {
                    stop("lavaan ERROR: length(sample.nobs) = ",
                          length(sample.nobs),
                         " but syntax implies ngroups = ", ngroups)
                }
            }
        }


        # group.label
        if(ngroups > 1L) {
            if(is.null(group)) {
                group <- "group"
            }
            group.label <- paste("Group", 1:ngroups, sep="")
        } else {
            group <- character(0L)
            group.label <- character(0L)
        }

        # handle ov.names.x
        if(!is.list(ov.names.x)) {
            ov.names.x <- rep(list(ov.names.x), ngroups)
        }

        ov <- list()
        ov$name <- unique( unlist(c(ov.names, ov.names.x)) )
        nvar    <- length(ov$name)
        ov$idx  <- rep(NA, nvar)
        ov$nobs <- rep(0L, nvar)
        ov$type <- rep("numeric", nvar)
        ov$nlev <- rep(0L, nvar)

        # collect information per upper-level group
        Lp <- vector("list", length = ngroups)
        for(g in 1:ngroups) {
            if(nlevels > 1L) {
                Lp[[g]] <- lav_data_cluster_patterns(Y = NULL, clus = NULL,
                                                 cluster = cluster,
                                                 multilevel = TRUE,
                                                 ov.names = ov.names[[g]],
                                                 ov.names.l = ov.names.l[[g]])
            }
        } # g

        # construct lavData object
        lavData <- new("lavData",
                       data.type   = "none",
                       ngroups     = ngroups,
                       group       = group,
                       nlevels     = nlevels,
                       cluster     = cluster,
                       group.label = group.label,
                       level.label = level.label,
                       nobs        = sample.nobs,
                       norig       = sample.nobs,
                       ov.names    = ov.names,
                       ov.names.x  = ov.names.x,
                       ov.names.l  = ov.names.l,
                       ordered     = as.character(ordered),
                       weights     = vector("list", length = ngroups),
                       sampling.weights = character(0L),
                       ov          = ov,
                       missing     = "listwise",
                       case.idx    = vector("list", length = ngroups),
                       Mp          = vector("list", length = ngroups),
                       Rp          = vector("list", length = ngroups),
                       Lp          = Lp,
                       X           = vector("list", length = ngroups),
                       eXo         = vector("list", length = ngroups)
                      )
    }

    lavData
}


# handle full data
lav_data_full <- function(data          = NULL,          # data.frame
                          group         = NULL,          # multiple groups?
                          cluster       = NULL,          # clustered?
                          group.label   = NULL,          # custom group labels?
                          level.label   = NULL,
                          ov.names      = NULL,          # variables needed
                                                         # in model
                          ordered       = NULL,          # ordered variables
                          sampling.weights = NULL,       # sampling weights
                          ov.names.x    = character(0L), # exo variables
                          ov.names.l    = list(),        # var per level
                          std.ov        = FALSE,         # standardize ov's?
                          missing       = "listwise",    # remove missings?
                          warn          = TRUE,          # produce warnings?
                          allow.single.case = FALSE      # allow single case?
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
            group.label <- unique(as.character(data[[group]]))
            if(warn && any(is.na(group.label))) {
                warning("lavaan WARNING: group variable ", sQuote(group),
                        " contains missing values\n", sep="")
            }
            group.label <- group.label[!is.na(group.label)]
        } else {
            group.label <- unique(as.character(group.label))
            # check if user-provided group labels exist
            LABEL <- unique(as.character(data[[group]]))
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

    # sampling weights
    if(!is.null(sampling.weights)) {
        if(is.character(sampling.weights)) {
            if(!(sampling.weights %in% names(data))) {
                stop("lavaan ERROR: sampling weights variable ",
                     sQuote(sampling.weights),
                     " not found;\n  ",
                     "variable names found in data frame are:\n  ",
                     paste(names(data), collapse=" "))
            }
            # check for missing values in sampling weight variable
            if(any(is.na(data[[sampling.weights]]))) {
                stop("lavaan ERROR: sampling.weights variable ",
                        sQuote(sampling.weights),
                        " contains missing values\n", sep = "")
            }
        } else {
            stop("lavaan ERROR: sampling weights argument should be a variable name in the data.frame")
        }
    }

    # clustered?
    if(!is.null(cluster) && length(cluster) > 0L) {

        # cluster variable in data?
        if(!all(cluster %in% names(data))) {
            # which one did we not find?
            not.ok <- which(!cluster %in% names(data))

            stop("lavaan ERROR: cluster variable(s) ", sQuote(cluster[not.ok]),
                 " not found;\n  ",
                 "variable names found in data frame are:\n  ",
                 paste(names(data), collapse = " "))
        }

        # check for missing values in cluster variable(s)
        for(cl in 1:length(cluster)) {
            if(warn && anyNA(data[[cluster[cl]]])) {
                warning("lavaan WARNING: cluster variable ",
                        sQuote(cluster[cl]),
                        " contains missing values\n", sep = "")
            }
        }

        # multilevel?
        if(length(ov.names.l) > 0L) {
            # default level.labels
            if(length(level.label) == 0L) {
                level.label <- c("within", cluster)
            } else {
                # check if length(level.label) = 1 + length(cluster)
                if(length(level.label) != length(cluster) + 1L) {
                    stop("lavaan ERROR: length(level.label) != length(cluster) + 1L")
                }
                # nothing to do
            }
            nlevels <- length(level.label)
        } else {
            # just clustered data, but no random effects
            nlevels <- 1L
            level.label <- character(0L)
        }
    } else {
        if(warn && length(level.label) > 0L)
            warning("lavaan WARNING: `level.label' argument",
                    " will be ignored if `cluster' argument is missing")
        nlevels <- 1L
        level.label <- character(0L)
        cluster <- character(0L)
    }

    # check ov.names vs ngroups
    if(ngroups > 1L) {
        if(is.list(ov.names)) {
            if(length(ov.names) != ngroups)
                stop("lavaan ERROR: ov.names assumes ", length(ov.names),
                     " groups; data contains ", ngroups, " groups")
        } else {
            tmp <- ov.names
            ov.names <- vector("list", length = ngroups)
            ov.names[1:ngroups] <- list(tmp)
        }
        if(is.list(ov.names.x)) {
            if(length(ov.names.x) != ngroups)
                stop("lavaan ERROR: ov.names assumes ", length(ov.names.x),
                     " groups; data contains ", ngroups, " groups")
        } else {
            tmp <- ov.names.x
            ov.names.x <- vector("list", length = ngroups)
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

    # check if all ov.names can be found in the data.frame
    for(g in 1:ngroups) {
        # does the data contain all the observed variables
        # needed in the user-specified model for this group
        ov.all <- unique(ov.names[[g]], ov.names.x[[g]]) # no overlap if categ

        # handle interactions
        ov.int.names <- ov.all[ grepl(":", ov.all) ]
        n.int <- length(ov.int.names)
        if(n.int > 0L) {
            ov.names.noint <- ov.all[!ov.all %in% ov.int.names]
            for(iv in seq_len(n.int)) {
                NAMES <- strsplit(ov.int.names[iv], ":", fixed = TRUE)[[1L]]
                if(all(NAMES %in% ov.names.noint)) {
                    # add this interaction term to the data.frame, unless
                    # it already exists
                    if(is.null(data[[ ov.int.names[iv] ]])) {
                        data[[ ov.int.names[iv] ]] <-
                            data[[NAMES[1L]]] * data[[NAMES[2L]]]
                    }
                }
            }
        }

        # check for missing observed variables
        idx.missing <- which(!(ov.all %in% names(data)))

        if(length(idx.missing)) {
            stop("lavaan ERROR: missing observed variables in dataset: ",
                 paste(ov.all[idx.missing], collapse=" "))
        }
    }


    # here, we know for sure all ov.names exist in the data.frame
    # create varTable
    # FIXME: should we add the 'group'/'cluster' variable (no for now)
    ov <- lav_dataframe_vartable(frame = data, ov.names = ov.names,
                                 ov.names.x = ov.names.x, ordered = ordered,
                                 as.data.frame. = FALSE)

    # do some checking
    # check for unordered factors (but only if nlev > 2)
    if("factor" %in%  ov$type) {
        f.names     <- ov$name[ov$type == "factor" & ov$nlev > 2L]
        f.names.all <- ov$name[ov$type == "factor"]
        OV.names <- unlist(ov.names)
        OV.names.x <- unlist(ov.names.x)
        OV.names.nox <- OV.names[! OV.names %in% OV.names.x]
        if(any(f.names %in% OV.names.x)) {
            stop(paste("lavaan ERROR: unordered factor(s) with more than 2 levels detected as exogenous covariate(s):", paste(f.names, collapse=" ")))
        } else if(any(f.names.all %in% OV.names.nox)) {
            stop(paste("lavaan ERROR: unordered factor(s) detected; make them numeric or ordered:", paste(f.names.all, collapse=" ")))
        }
    }
    # check for ordered exogenous variables
    if("ordered" %in% ov$type[ov$name %in% unlist(ov.names.x)]) {
        f.names <- ov$name[ov$type == "ordered" &
                           ov$name %in% unlist(ov.names.x)]
        if(warn && any(f.names %in% unlist(ov.names.x)))
            warning(paste("lavaan WARNING: exogenous variable(s) declared as ordered in data:", paste(f.names, collapse=" ")))
    }
    # check for zero-cases
    idx <- which(ov$nobs == 0L | ov$var == 0)
    if(length(idx) > 0L) {
        OV <- as.data.frame(ov)
        rn <- rownames(OV)
        rn[idx] <- paste(rn[idx], "***", sep="")
        rownames(OV) <- rn
        print(OV)
        stop("lavaan ERROR: some variables have no values (only missings) or no variance")
    }
    # check for single cases (no variance!)
    idx <- which(ov$nobs == 1L | (ov$type == "numeric" & !is.finite(ov$var)))
    if(!allow.single.case && length(idx) > 0L) {
        OV <- as.data.frame(ov)
        rn <- rownames(OV)
        rn[idx] <- paste(rn[idx], "***", sep="")
        rownames(OV) <- rn
        print(OV)
        stop("lavaan ERROR: some variables have only 1 observation or no finite variance")
    }
    # check for ordered variables with only 1 level
    idx <- which(ov$type == "ordered" & ov$nlev == 1L)
    if(length(idx) > 0L) {
        OV <- as.data.frame(ov)
        rn <- rownames(OV)
        rn[idx] <- paste(rn[idx], "***", sep="")
        rownames(OV) <- rn
        print(OV)
        stop("lavaan ERROR: ordered variable(s) has/have only 1 level")
    }
    # check for mix small/large variances (NOT including exo variables)
    if(!std.ov && !allow.single.case && warn && any(ov$type == "numeric")) {
        num.idx <- which(ov$type == "numeric" & ov$exo == 0L)
        if(length(num.idx) > 0L) {
            min.var <- min(ov$var[num.idx])
            max.var <- max(ov$var[num.idx])
            rel.var <- max.var/min.var
            if(warn && rel.var > 1000) {
                warning("lavaan WARNING: some observed variances are (at least) a factor 1000 times larger than others; use varTable(fit) to investigate")
            }
        }
    }
    # check for really large variances (perhaps -999999 for missing?)
    if(!std.ov && warn && any(ov$type == "numeric")) {
        num.idx <- which(ov$type == "numeric" & ov$exo == 0L)
        if(length(num.idx) > 0L) {
            max.var <- max(ov$var[num.idx])
            if(warn && max.var > 1000000) {
                warning("lavaan WARNING: some observed variances are larger than 1000000\n", "  lavaan NOTE: use varTable(fit) to investigate")
            }
        }
    }
    # check for all-exogenous variables (eg in f <~ x1 + x2 + x3)
    if(warn && all(ov$exo == 1L)) {
        warning("lavaan WARNING: all observed variables are exogenous; model may not be identified")
    }

    # prepare empty lists

    # group-based
    case.idx <- vector("list", length = ngroups)
    Mp       <- vector("list", length = ngroups)
    Rp       <- vector("list", length = ngroups)
    norig    <- vector("list", length = ngroups)
    nobs     <- vector("list", length = ngroups)
    X        <- vector("list", length = ngroups)
    eXo      <- vector("list", length = ngroups)
    Lp       <- vector("list", length = ngroups)
    weights  <- vector("list", length = ngroups)

    # collect information per upper-level group
    for(g in 1:ngroups) {

        # extract variables in correct order
        ov.idx  <- ov$idx[match(ov.names[[g]],   ov$name)]
        exo.idx <- ov$idx[match(ov.names.x[[g]], ov$name)]
        all.idx <- unique(c(ov.idx, exo.idx))

        # extract cases per group
        if(ngroups > 1L || length(group.label) > 0L) {
            if(missing == "listwise") {
                case.idx[[g]] <- which(data[[group]] == group.label[g] &
                                           complete.cases(data[all.idx]))
                nobs[[g]] <- length(case.idx[[g]])
                norig[[g]] <- length(which(data[[group]] == group.label[g]))
            #} else if(missing == "pairwise" && length(exo.idx) > 0L) {
            #    case.idx[[g]] <- which(data[[group]] == group.label[g] &
            #                           complete.cases(data[exo.idx]))
            #    nobs[[g]] <- length(case.idx[[g]])
            #    norig[[g]] <- length(which(data[[group]] == group.label[g]))
            } else if(length(exo.idx) > 0L && missing != "ml.x") {
                case.idx[[g]] <- which(data[[group]] == group.label[g] &
                                       complete.cases(data[exo.idx]))
                nobs[[g]] <- length(case.idx[[g]])
                norig[[g]] <- length(which(data[[group]] == group.label[g]))
                if(warn && (nobs[[g]] < norig[[g]])) {
                    warning("lavaan WARNING: ", (norig[[g]] - nobs[[g]]),
                        " cases were deleted in group ", group.label[g],
                        " due to missing values in ",
                        "\n\t\t  exogenous variable(s), while fixed.x = TRUE.")
                }
            } else {
                case.idx[[g]] <- which(data[[group]] == group.label[g])
                nobs[[g]] <- norig[[g]] <- length(case.idx[[g]])
            }
        } else {
            if(missing == "listwise") {
                case.idx[[g]] <- which(complete.cases(data[all.idx]))
                nobs[[g]] <- length(case.idx[[g]])
                norig[[g]] <- nrow(data)
            #} else if(missing == "pairwise" && length(exo.idx) > 0L) {
            #    case.idx[[g]] <- which(complete.cases(data[exo.idx]))
            #    nobs[[g]] <- length(case.idx[[g]])
            #    norig[[g]] <- nrow(data)
            } else if(length(exo.idx) > 0L && missing != "ml.x") {
                case.idx[[g]] <- which(complete.cases(data[exo.idx]))
                nobs[[g]] <- length(case.idx[[g]])
                norig[[g]] <- nrow(data)
                if(warn && (nobs[[g]] < norig[[g]])) {
                    warning("lavaan WARNING: ", (norig[[g]] - nobs[[g]]),
                        " cases were deleted due to missing values in ",
                        "\n\t\t  exogenous variable(s), while fixed.x = TRUE.")
                }
            } else {
                case.idx[[g]] <- 1:nrow(data)
                nobs[[g]] <- norig[[g]] <- length(case.idx[[g]])
            }
        }

        # extract data
        X[[g]] <- data.matrix( data[case.idx[[g]], ov.idx, drop = FALSE] )
        dimnames(X[[g]]) <- NULL ### copy?

        if(!is.null(sampling.weights)) {
            WT <- data[[sampling.weights]][case.idx[[g]]]
            if(any(WT < 0)) {
                stop("lavaan ERROR: some sampling weights are negative")
            }

            # check for missing values in sampling weight variable
            if(any(is.na(WT))) {
                stop("lavaan ERROR: sampling.weights variable ",
                        sQuote(sampling.weights),
                        " contains missing values\n", sep = "")
            }

            # rescale, so sum equals sample size in this group
            WT2 <- WT / sum(WT) * nobs[[g]]
            weights[[g]] <- WT2
        }

        # construct integers for user-declared 'ordered' factors
        # FIXME: is this really (always) needed???
        #  (but still better than doing lapply(data[,idx], ordered) which
        #   generated even more copies)
        user.ordered.names <- ov$name[ov$type == "ordered" & ov$user == 1L]
        user.ordered.idx <- which(ov.names[[g]] %in% user.ordered.names)
        if(length(user.ordered.idx) > 0L) {
            for(i in user.ordered.idx) {
                X[[g]][,i] <- as.numeric(as.factor(X[[g]][,i]))
            }
        }

        ## FIXME:
        ## - why also in X? (for samplestats, for now)
        if(length(exo.idx) > 0L) {
            eXo[[g]] <- data.matrix(data[case.idx[[g]], exo.idx, drop = FALSE])
            dimnames(eXo[[g]]) <- NULL
        } else {
            eXo[g] <- list(NULL)
        }

        # standardize observed variables? numeric only!
        if(std.ov) {
            num.idx <- which(ov$name %in% ov.names[[g]] &
                             ov$type == "numeric" & ov$exo == 0L)
            if(length(num.idx) > 0L) {
                X[[g]][,num.idx] <-
                   scale(X[[g]][,num.idx,drop = FALSE])[,,drop = FALSE]
                # three copies are made!!!!!
            }
            if(length(exo.idx) > 0L) {
                eXo[[g]] <- scale(eXo[[g]])[,,drop = FALSE]
            }
        }

        # missing data
        if(missing != "listwise") {
            # get missing patterns
            Mp[[g]] <- lav_data_missing_patterns(X[[g]],
                           sort.freq = TRUE, coverage = TRUE)
            # checking!
            if(length(Mp[[g]]$empty.idx) > 0L) {
                empty.case.idx <- Mp[[g]]$empty.idx
                if(warn) {
                    warning("lavaan WARNING: some cases are empty and will be ignored:\n  ", paste(empty.case.idx, collapse=" "))
                }
            }
            if(warn && any(Mp[[g]]$coverage < 0.1)) {
                warning("lavaan WARNING: due to missing values, some pairwise combinations have less than 10% coverage")
            }
            # in case we had observations with only missings
            nobs[[g]] <- NROW(X[[g]]) - length(Mp[[g]]$empty.idx)
        }

        # response patterns (ordered variables only)
        ord.idx <- which(ov.names[[g]] %in% ov$name[ov$type == "ordered"])
        if(length(ord.idx) > 0L) {
            Rp[[g]] <- lav_data_resp_patterns(X[[g]][,ord.idx, drop = FALSE])
        }

        # warn if we have a small number of observations (but NO error!)
        if( !allow.single.case && warn &&
            nobs[[g]] < (nvar <- length(ov.idx)) ) {
            txt <- ""
            if(ngroups > 1L) txt <- paste(" in group ", g, sep="")
            warning("lavaan WARNING: small number of observations (nobs < nvar)", txt,
                    "\n  nobs = ", nobs[[g]], " nvar = ", nvar)
        }

        # cluster information
        if(length(cluster) > 0L) {
            # extract cluster variable(s), for this group
            clus <- data.matrix(data[case.idx[[g]], cluster])
            if(nlevels > 1L) {
                multilevel <- TRUE
            } else {
                multilevel <- FALSE
            }
            Lp[[g]] <- lav_data_cluster_patterns(Y = X[[g]], clus = clus,
                                                 cluster = cluster,
                                                 multilevel = multilevel,
                                                 ov.names = ov.names[[g]],
                                                 ov.names.l = ov.names.l[[g]])

            # new in 0.6-4
            # check for 'level-1' variables with zero within variance
            l1.idx <- c(Lp[[g]]$within.idx[[2]], # within only
                        Lp[[g]]$both.idx[[2]])
            l1.names <- c(Lp[[g]]$within.names[[2]],
                          Lp[[g]]$both.names[[2]])
            for(v in l1.idx) {
                within.var <- tapply(X[[g]][,v], Lp[[g]]$cluster.idx[[2]],
                                     FUN = var, na.rm = TRUE)
                # ignore singeltons
                singleton.idx <- which( Lp[[g]]$cluster.size[[2]] == 1L )
                if(length(singleton.idx) > 0L) {
                    within.var[singleton.idx] <- 10 # non-zero variance
                }
                zero.var <- which(within.var < .Machine$double.eps)
                if(length(zero.var) == 0L) {
                    # all is good
                } else if(length(zero.var) == length(within.var)) {
                    # all zero! possibly a between-level variable
                    gtxt <- if(ngroups > 1L) {
                                paste(" in group ", g, ".", sep = "")
                            } else { "." }
                    txt <- c("Level-1 variable ", dQuote(l1.names[v]),
                             " has no variance at the within level", gtxt,
                             " The variable appears to be a between-level
                             variable. Please remove this variable from
                             the level 1 section in the model syntax.")
                    warning(lav_txt2message(txt))
                } else {
                    # some zero variances!
                    gtxt <- if(ngroups > 1L) {
                                paste(" in group ", g, ".", sep = "")
                            } else { "." }
                    txt <- c("Level-1 variable ", dQuote(l1.names[v]),
                       " has no variance within some clusters", gtxt,
                       " The cluster ids with zero within variance are:\n",
                       paste( Lp[[g]]$cluster.id[[2]][zero.var],
                              collapse = " "))
                    warning(lav_txt2message(txt))
                }
            }

        } # clustered data

    } # groups, at first level

    if(is.null(sampling.weights)) {
        sampling.weights <- character(0L)
    }

    lavData <- new("lavData",
                   data.type       = "full",
                   ngroups         = ngroups,
                   group           = group,
                   nlevels         = nlevels,
                   cluster         = cluster,
                   group.label     = group.label,
                   level.label     = level.label,
                   std.ov          = std.ov,
                   nobs            = nobs,
                   norig           = norig,
                   ov.names        = ov.names,
                   ov.names.x      = ov.names.x,
                   ov.names.l      = ov.names.l,
                   #ov.types        = ov.types,
                   #ov.idx          = ov.idx,
                   ordered         = as.character(ordered),
                   weights         = weights,
                   sampling.weights = sampling.weights,
                   ov              = ov,
                   case.idx        = case.idx,
                   missing         = missing,
                   X               = X,
                   eXo             = eXo,
                   Mp              = Mp,
                   Rp              = Rp,
                   Lp              = Lp
                  )
    lavData
}

# get missing patterns
lav_data_missing_patterns <- function(Y, sort.freq = FALSE, coverage = FALSE) {

    # construct TRUE/FALSE matrix: TRUE if value is observed
    OBS <- !is.na(Y)

    # empty cases
    empty.idx <- which(rowSums(OBS) == 0L)

    # this is what we did in < 0.6
    #if(length(empty.idx) > 0L) {
    #    OBS <- OBS[-empty.idx,,drop = FALSE]
    #}

    # pattern of observed values per observation
    case.id <- apply(1L * OBS, 1L, paste, collapse = "")

    # remove empty patterns
    if(length(empty.idx)) {
        case.id.nonempty <- case.id[-empty.idx]
    } else {
        case.id.nonempty <- case.id
    }

    # sort non-empty patterns (from high occurence to low occurence)
    if(sort.freq) {
        TABLE <- sort(table(case.id.nonempty), decreasing = TRUE)
    } else {
        TABLE <- table(case.id.nonempty)
    }

    # unique pattern ids
    pat.id <- names(TABLE)

    # number of patterns
    pat.npatterns  <- length(pat.id)

    # case idx per pattern
    pat.case.idx <- lapply(seq_len(pat.npatterns),
                           function(p) which(case.id == pat.id[p]))

    # unique pattern frequencies
    pat.freq <- as.integer(TABLE)

    # first occurrence of each pattern
    pat.first <- match(pat.id, case.id)

    # TRUE/FALSE for each pattern
    pat.obs <- OBS[pat.first,,drop = FALSE] # observed per pattern

    Mp <- list(npatterns = pat.npatterns, id = pat.id, freq = pat.freq,
               case.idx = pat.case.idx, pat = pat.obs, empty.idx = empty.idx)

    if(coverage) {
        # FIXME: if we have empty cases, include them in N?
        # no for now
        Mp$coverage <- crossprod(OBS) / sum(pat.freq)
        #Mp$coverage <- crossprod(OBS) / NROW(Y)
    }

    Mp
}

# get response patterns (ignore empty cases!)
lav_data_resp_patterns <- function(Y) {

    # construct TRUE/FALSE matrix: TRUE if value is observed
    OBS <- !is.na(Y)

    # empty cases
    empty.idx <- which(rowSums(OBS) == 0L)

    # removeYempty cases
    if(length(empty.idx) > 0L) {
        Y <- Y[-empty.idx,,drop = FALSE]
    }

    ntotal <- nrow(Y); nvar <- ncol(Y)

    # identify, label and sort response patterns
    id <- apply(Y, MARGIN = 1, paste, collapse = "")

    # sort patterns (from high occurence to low occurence)
    TABLE <- sort(table(id), decreasing = TRUE)
    order <- names(TABLE)
    npatterns <- length(TABLE)
    pat <- Y[match(order, id), , drop = FALSE]
    row.names(pat) <- as.character(TABLE)

    # handle NA?
    Y[is.na(Y)] <- -9
    total.patterns <- prod(apply(Y, 2, function(x) length(unique(x))))
    empty.patterns <- total.patterns - npatterns
    # return a list
    #out <- list(nobs=ntotal, nvar=nvar,
    #            id=id, npatterns=npatterns,
    #            order=order, pat=pat)

    # only return pat
    out <- list(npatterns=npatterns, pat=pat, total.patterns=total.patterns,
                empty.patterns=empty.patterns)

    out
}

# get cluster information
# - cluster can be a vector!
# - clus can contain multiple columns!
lav_data_cluster_patterns <- function(Y = NULL,
                                      clus = NULL,    # the cluster ids
                                      cluster = NULL, # the cluster 'names'
                                      multilevel = FALSE,
                                      ov.names, ov.names.l) {

    # how many levels?
    nlevels <- length(cluster) + 1L

    # did we get any data (or is this just for simulateData)
    if(!is.null(Y) && !is.null(clus)) {
        haveData <- TRUE
    } else {
        haveData <- FALSE
    }

    # check clus
    if(haveData) {
        stopifnot(ncol(clus) == (nlevels - 1L), nrow(Y) == nrow(clus))
    }

    cluster.size    <- vector("list", length = nlevels)
    cluster.id      <- vector("list", length = nlevels)
    cluster.idx     <- vector("list", length = nlevels)
    nclusters       <- vector("list", length = nlevels)
    cluster.sizes   <- vector("list", length = nlevels)
    ncluster.sizes  <- vector("list", length = nlevels)
    cluster.size.ns <- vector("list", length = nlevels)
    ov.idx          <- vector("list", length = nlevels)
    both.idx        <- vector("list", length = nlevels)
    within.idx      <- vector("list", length = nlevels)
    between.idx     <- vector("list", length = nlevels)
    both.names      <- vector("list", length = nlevels)
    within.names    <- vector("list", length = nlevels)
    between.names   <- vector("list", length = nlevels)

    # level-1 is special
    if(haveData) {
        nclusters[[1]] <- NROW(Y)
    }

    if(multilevel) {
        ov.idx[[1]] <- match(ov.names.l[[1]], ov.names)
    }

    # for the remaining levels...
    for(l in 2:nlevels) {
        if(haveData) {
            CLUS <- clus[,(l-1L)]
            cluster.id[[l]]      <- unique(CLUS)
            cluster.idx[[l]]     <- match(CLUS, cluster.id[[l]])
            cluster.size[[l]]    <- tabulate(cluster.idx[[l]])
            nclusters[[l]]       <- length(cluster.size[[l]])
            cluster.sizes[[l]]   <- unique(cluster.size[[l]])
            ncluster.sizes[[l]]  <- length(cluster.sizes[[l]])
            cluster.size.ns[[l]] <- as.integer(table(factor(cluster.size[[l]],
                                     levels = as.character(cluster.sizes[[l]]))))
        } else {
            cluster.id[[l]]      <- integer(0L)
            cluster.idx[[l]]     <- integer(0L)
            cluster.size[[l]]    <- integer(0L)
            nclusters[[l]]       <- integer(0L)
            cluster.sizes[[l]]   <- integer(0L)
            ncluster.sizes[[l]]  <- integer(0L)
            cluster.size.ns[[l]] <- integer(0L)
        }

        if(multilevel) {
            # index of ov.names for this level
            ov.idx[[l]]         <- match(ov.names.l[[l]], ov.names)

            both.idx[[l]]       <- which( ov.names %in% ov.names.l[[1]] &
                                          ov.names %in% ov.names.l[[2]])
            within.idx[[l]]     <- which( ov.names %in% ov.names.l[[1]] &
                                         !ov.names %in% ov.names.l[[2]])
            between.idx[[l]]    <- which(!ov.names %in% ov.names.l[[1]] &
                                          ov.names %in% ov.names.l[[2]])

            # names
            both.names[[l]]     <- ov.names[ ov.names %in% ov.names.l[[1]] &
                                             ov.names %in% ov.names.l[[2]] ]
            within.names[[l]]   <- ov.names[ ov.names %in% ov.names.l[[1]] &
                                            !ov.names %in% ov.names.l[[2]] ]
            between.names[[l]]  <- ov.names[!ov.names %in% ov.names.l[[1]] &
                                             ov.names %in% ov.names.l[[2]] ]
        }
    }

    out <- list(cluster = cluster, # clus = clus,
                # per level
                nclusters = nclusters,
                cluster.size = cluster.size, cluster.id = cluster.id,
                cluster.idx = cluster.idx, cluster.sizes = cluster.sizes,
                ncluster.sizes = ncluster.sizes,
                cluster.size.ns = cluster.size.ns,
                ov.idx = ov.idx, both.idx = both.idx, within.idx = within.idx,
                between.idx = between.idx,
                both.names = both.names, within.names = within.names,
                between.names = between.names)

    out
}

setMethod("show", "lavData",
function(object) {
    # print 'lavData' object
    lav_data_print_short(object)
})

lav_data_print_short <- function(object) {

    lavdata <- object

    # listwise deletion?
    listwise <- FALSE
    for(g in 1:lavdata@ngroups) {
       if(lavdata@nobs[[1L]] != lavdata@norig[[1L]]) {
           listwise <- TRUE
           break
       }
    }

    #cat("Data information:\n\n")
    if(lavdata@ngroups == 1L) {
        if(listwise) {
            cat(sprintf("  %-40s", ""), sprintf("  %10s", "Used"),
                                        sprintf("  %10s", "Total"),
                "\n", sep="")
        }
        t0.txt <- sprintf("  %-40s", "Number of observations")
        t1.txt <- sprintf("  %10i", lavdata@nobs[[1L]])
        t2.txt <- ifelse(listwise,
                  sprintf("  %10i", lavdata@norig[[1L]]), "")
        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

        if( (.hasSlot(lavdata, "nlevels")) && # in case we have an old obj
            (lavdata@nlevels > 1L) ) {
            #cat("\n")
            for(l in 2:lavdata@nlevels) {
                t0.txt <- sprintf("  %-40s",
                    paste("Number of clusters [", lavdata@cluster[l-1], "]",
                          sep = ""))
                t1.txt <- sprintf("  %10i", lavdata@Lp[[1]]$nclusters[[l]])
                #t2.txt <- ifelse(listwise,
                #          sprintf("  %10i", lavdata@norig[[1L]]), "")
                t2.txt <- ""
                cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
            }
        } else if( (.hasSlot(lavdata, "cluster")) &&
                   (length(lavdata@cluster) > 0L) ) {
            t0.txt <- sprintf("  %-40s",
            paste("Number of clusters [", lavdata@cluster, "]", sep = ""))
            t1.txt <- sprintf("  %10i", lavdata@Lp[[1]]$nclusters[[2]])
            t2.txt <- ""
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }
    } else {
        if(listwise) {
            cat(sprintf("  %-40s", ""), sprintf("  %10s", "Used"),
                                       sprintf("  %10s", "Total"),
               "\n", sep="")
        }
        t0.txt <- sprintf("  %-40s", "Number of observations per group")
        cat(t0.txt, "\n")
        for(g in 1:lavdata@ngroups) {
            t.txt <- sprintf("  %-40s  %10i", lavdata@group.label[[g]],
                                              lavdata@nobs[[g]])
            t2.txt <- ifelse(listwise,
                      sprintf("  %10i", lavdata@norig[[g]]), "")
            cat(t.txt, t2.txt, "\n", sep="")

            if( (.hasSlot(lavdata, "nlevels")) &&
                (lavdata@nlevels > 1L) ) {
                #cat("\n")
                for(l in 2:lavdata@nlevels) {
                    t0.txt <- sprintf("  %-40s",
                        paste("Number of clusters [", lavdata@cluster[l-1], "]",
                              sep = ""))
                    t1.txt <- sprintf("  %10i",
                                      lavdata@Lp[[g]]$nclusters[[l]])
                    #t2.txt <- ifelse(listwise,
                    #          sprintf("  %10i", lavdata@norig[[1L]]), "")
                    t2.txt <- ""
                    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
                }
            } else if( (.hasSlot(lavdata, "cluster")) &&
                   (length(lavdata@cluster) > 0L) ) {
                t0.txt <- sprintf("  %-40s",
                paste("Number of clusters [", lavdata@cluster, "]", sep = ""))
                t1.txt <- sprintf("  %10i", lavdata@Lp[[g]]$nclusters[[2]])
                t2.txt <- ""
                cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
            }
        } # g
    }

    # missing patterns?
    if(!is.null(lavdata@Mp[[1L]])) {
        if(lavdata@ngroups == 1L) {
            t0.txt <- sprintf("  %-40s", "Number of missing patterns")
            t1.txt <- sprintf("  %10i",
                              lavdata@Mp[[1L]]$npatterns)
            cat(t0.txt, t1.txt, "\n", sep="")
        } else {
            t0.txt <- sprintf("  %-40s", "Number of missing patterns per group")
            cat(t0.txt, "\n")
            for(g in 1:lavdata@ngroups) {
                t.txt <- sprintf("  %-40s  %10i", lavdata@group.label[[g]],
                                 lavdata@Mp[[g]]$npatterns)
                cat(t.txt, "\n", sep="")
            }
        }
    }

    # sampling weights?
    if( (.hasSlot(lavdata, "weights")) && # in case we have an old object
        (!is.null(lavdata@weights[[1L]])) ) {
        t0.txt <- sprintf("  %-30s", "Sampling weights variable")
        t1.txt <- sprintf("  %20s", lavdata@sampling.weights)
        cat(t0.txt, t1.txt, "\n", sep="")
    }

    cat("\n")
}

