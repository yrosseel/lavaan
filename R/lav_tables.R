# construct pairwise frequency tables
# YR. 10 April 2013
# Notes:
# - we do NOT make a distinction here between unordered and ordered categorical
#   variables
# - object can be a matrix (most likely with integers), a full data frame, 
#   a fitted lavaan object, a varTable, or a lavData object
# - 11 May 2013: added collapse=TRUE, min.std.resid options (suggested
#   by Myrsini Katsikatsou
# - 11 June 2013: added dimension, to get one-way and two-way (three-way?)
#   tables

lavTables <- function(object, dimension=2L, categorical=NULL,
                      std.resid=TRUE, min.std.resid = 0.0,
                      average = FALSE, collapse = FALSE) {

    # catch matrix
    if(is.matrix(object)) {
        object <- as.data.frame(object, stringsAsFactors = FALSE)
        if(is.null(categorical)) {
            object[,] <- lapply(object, base::factor)
        }
    }

    # check object class
    if(!class(object) %in% c("data.frame", "lavData", "lavaan")) {
        stop("lavaan ERROR: object must either be a matrix or a data.frame, or an object of class lavaan or class lavData")
    }

    if(inherits(object, "data.frame")) {
        ov.names <- names(object)
        vartable <- varTable(object, ov.names=ov.names, factor=categorical,
                             as.data.frame.=FALSE)
        X <- data.matrix(object)
        # manually construct integers for user-declared categorical variables
        user.categorical.names <-
            vartable$name[vartable$type %in% c("ordered","factor") &
                          vartable$user == 1L]
        user.categorical.idx <- which(ov.names %in% user.categorical.names)
        for(i in seq_len(length(user.categorical.idx))) {
            X[,i] <- as.numeric(as.factor(X[,i]))
        }
        X <- list(X); ov.names <- list(ov.names)
    } else if(inherits(object, "lavData")) {
        vartable <- object@ov
        X <- object@X
        ov.names <- object@ov.names
    } else if(inherits(object, "lavaan")) {
        vartable <- object@Data@ov
        X <- object@Data@X
        ov.names <- object@Data@ov.names
    }

    if(dimension == 1L) {
        out <- lav_oneway_tables(object = object, 
                   vartable = vartable, X = X, ov.names = ov.names)
    } else if(dimension == 2L) {
        out <- lav_pairwise_tables(object = object,
                   vartable = vartable, X = X, ov.names = ov.names,
                   std.resid = std.resid, min.std.resid = min.std.resid,
                   average = average, collapse = collapse)
    } else {
        stop("lavaan ERROR: dimension must be 1 or 2 of one-way or two-way tables")
    }

    out
}


lav_oneway_tables <- function(object, vartable=NULL, X=NULL, ov.names=NULL) {

    # identify 'categorical' variables
    cat.idx <- which(vartable$type %in% c("ordered","factor"))
    ncat <- length(cat.idx)

    # do we have any categorical variables?
    if(length(cat.idx) == 0L) {
        stop("lavaan ERROR: no categorical variables are found")
    } else {
        labels <- strsplit(vartable$lnam[cat.idx], "\\|")
    }

    # ok, we have an overview of all categorical variables in the data
    ngroups <- length(X)

    # for each group, for each categorical variable, collect information
    TABLES <- vector("list", length=ngroups)
    for(g in 1:ngroups) {
        TABLES[[g]] <- lapply(seq_len(ncat),
            FUN=function(x) {
                idx <- cat.idx[x]
                nrow <- vartable$nlev[idx]
                ncell<- nrow

                # compute observed frequencies
                FREQ <- tabulate(X[[g]][,idx])

                list(   id = rep.int(x, ncell),
                       lhs = rep.int(vartable$name[idx], ncell),
                       # op = rep.int("freq", ncell), 
                       rhs = labels[[x]],
                     group = rep.int(g, ncell),
                      nobs = rep.int(sum(FREQ), ncell),
                      obs.freq = FREQ,
                      obs.prop = FREQ/sum(FREQ) 
                    )
            })
    }

    for(g in 1:ngroups) {
        TABLE <- TABLES[[g]]
        TABLE <- lapply(TABLE, as.data.frame, stringAsFactors=FALSE)
        if(g == 1L) {
            out <- do.call(rbind, TABLE)
            # remove group column 
            out$group <- NULL
        } else {
            out <- rbind(out, do.call(rbind, TABLE))
        }
    }

    class(out) <- c("lavaan.data.frame", "data.frame")
    out
}

lav_pairwise_tables <- function(object, 
                                vartable=NULL, X=NULL, ov.names=NULL,
                                std.resid=TRUE, min.std.resid = 0.0,
                                average = FALSE, collapse = FALSE) {

    showAsMatrix <- FALSE
    if(! (is.logical(collapse) && !collapse) ) {
        average <- TRUE
        if(is.character(collapse) && tolower(collapse) == "matrix") {
            collapse <- TRUE
            showAsMatrix <- TRUE
        }
    }

    out <- lav_pairwise_tables_freq(vartable = vartable,
                                    X = X,
                                    ov.names = ov.names,
                                    as.data.frame. = TRUE)

    if(inherits(object, "data.frame") || inherits(object, "lavData")) {
        # compute chisq test for independence for each table + p.value
        # note: 'independent' tables (non-significant) are problematic here
        if(collapse) {
            ntables <- length(unique(out$id))
            ncells <- tabulate(out$id)
            out$ncells <- unlist(lapply(seq_len(ntables),
                    function(x) rep(ncells[x], each=ncells[x])))
            tableList <- lapply(seq_len(ntables),
                function(id) { matrix(out$obs.freq[     out$id == id ],
                                      nrow=max(out$row[ out$id == id ]),
                                      ncol=max(out$col[ out$id == id ]))  })
            tableChisqTest <- lapply(tableList, stats::chisq.test)
            table.chisq <- unlist(lapply(tableChisqTest, 
                function(x) as.numeric(x$statistic)))
            table.df <- unlist(lapply(tableChisqTest, 
                function(x) as.numeric(x$parameter)))
            table.pvalue <- unlist(lapply(tableChisqTest,   
                function(x) as.numeric(x$p.value)))
            out$chisq <- unlist(lapply(seq_len(ntables),
                    function(x) rep(table.chisq[x], each=ncells[x])))
            out$df <- unlist(lapply(seq_len(ntables),
                    function(x) rep(table.df[x], each=ncells[x])))
            out$pvalue <- unlist(lapply(seq_len(ntables),
                    function(x) rep(table.pvalue[x], each=ncells[x])))
        }
    } else if(inherits(object, "lavaan")) {
        # add predicted frequencies, fit indices

        #out$pi <- unlist(lav_pairwise_tables_pi(object))
        #out$freq.est <- out$pi * out$nobs
        obs.prop <- out$obs.freq/out$nobs
        est.prop <- unlist(lav_pairwise_tables_pi(object))
        out$est.freq <- est.prop * out$nobs
 
        # Joreskog & Moustaki equation 35
        if(std.resid) {
            out$std.resid <- out$nobs*(obs.prop-est.prop)^2/est.prop
            #if(check) {
            #    out$check <- ifelse(out$std.resid > 4.0, "***", "")
            #}

            if(average) {
                # fit per table
                table.average <- tapply(out$std.resid, INDEX=out$id, FUN=mean)
                table.percAboveMin <- tapply(out$std.resid, INDEX=out$id,
                    FUN=function(x) {sum(x > min.std.resid)/length(x)})
                table.numAboveMin <- tapply(out$std.resid, INDEX=out$id,
                    FUN=function(x) {sum(x > min.std.resid)})
                ntables <- length(table.average)
                ncells <- tabulate(out$id)
                out$ncells <- unlist(lapply(seq_len(ntables),
                    function(x) rep(ncells[x], each=ncells[x])))
                out$str.average <- unlist(lapply(seq_len(ntables),
                    function(x) rep(table.average[x], each=ncells[x])))
                out$str.min <- rep(min.std.resid, length(out$id))
                out$str.plarge <- unlist(lapply(seq_len(ntables),
                    function(x) rep(table.percAboveMin[x], each=ncells[x])))
                out$str.nlarge <- unlist(lapply(seq_len(ntables),
                    function(x) rep(table.numAboveMin[x], each=ncells[x])))

            } else {
                if(min.std.resid > 0.0) {
                    # select only rows where std.resid >= min.std.resid
                    idx <- which(out$std.resid >= min.std.resid)
                    out <- out[idx,]
                }
            }
        }
    }

    if(collapse) {
        # only 1 row per table
        row.idx <- which(!duplicated(out$id))
        out <- out[row.idx,,drop=FALSE]

        # remove some cell-specific columns
        out$row <- NULL; out$col <- NULL
        out$obs.freq <- NULL; out$est.freq <- NULL
        out$std.resid <- NULL; # out$ncells <- NULL
    }

    if(collapse && showAsMatrix) {
        tmp <- out
        RN <- unique(out$rhs)
        if(inherits(object, "lavaan")) {
            out <- lavaan::getCov(out$str.average, lower=FALSE, 
                                  names=unique(out$lhs))
            rownames(out) <- RN
        } else {
            out <- lavaan::getCov(out$pvalue, lower=FALSE, 
                                  names=unique(out$lhs))
            rownames(out) <- RN
        }
        class(out) <- c("lavaan.matrix.symmetric", "matrix")
    } else {
        class(out) <- c("lavaan.data.frame", "data.frame")
    }
    out
}

# internal function - X is a list of matrices!
lav_pairwise_tables_freq <- function(vartable = NULL, X = NULL, ov.names = NULL,
                                     as.data.frame. = TRUE) {

    stopifnot(is.list(X))

    # construct data.frame
    vartable <- as.data.frame(vartable, stringsAsFactors = FALSE)

    # identify 'categorical' variables
    cat.idx <- which(vartable$type %in% c("ordered","factor"))

    # do we have any categorical variables?
    if(length(cat.idx) == 0L) {
        stop("lavaan ERROR: no categorical variables are found")
    } else if(length(cat.idx) == 1L) {
        stop("lavaan ERROR: at least two categorical variables are needed")
    }

    # ok, we have an overview of all categorical variables in the data
    ngroups <- length(X)

    # pairwise tables
    pairwise.tables <- utils::combn(vartable$name[cat.idx], m=2L)
    pairwise.tables <- rbind(pairwise.tables, seq_len(ncol(pairwise.tables)))
    ntables <- ncol(pairwise.tables)

    # for each group, for each pairwise table, collect information
    TABLES <- vector("list", length=ngroups)
    for(g in 1:ngroups) {
        TABLES[[g]] <- apply(pairwise.tables, MARGIN=2,
            FUN=function(x) {
                idx1 <- which(vartable$name == x[1])
                idx2 <- which(vartable$name == x[2])
                id <- (g-1)*ntables + as.numeric(x[3]) 
                nrow <- vartable$nlev[idx1]
                ncol <- vartable$nlev[idx2]
                ncell <- nrow*ncol

                # compute two-way observed frequencies
                Y1 <- X[[g]][,idx1]
                Y2 <- X[[g]][,idx2]
                # FREQ <- table(Y1, Y2) # we loose missings; useNA is ugly
                FREQ <- pc_freq(Y1, Y2)
 
                list(   id = rep.int(id, ncell),
                       lhs = rep.int(x[1], ncell), 
                       # op = rep.int("table", ncell), 
                       rhs = rep.int(x[2], ncell),
                     group = rep.int(g, ncell),
                      nobs = rep.int(sum(FREQ), ncell),
                       row = rep.int(seq_len(ncol), times=nrow),
                       col = rep(seq_len(nrow), each=ncol),
                      obs.freq = vec(FREQ) # col by col!
                    )
            })
    }

    if(as.data.frame.) {
        for(g in 1:ngroups) {
            TABLE <- TABLES[[g]]
            TABLE <- lapply(TABLE, as.data.frame, stringAsFactors=FALSE)
            if(g == 1) {
                out <- do.call(rbind, TABLE)
                # remove group column 
                out$group <- NULL
            } else {
                out <- rbind(out, do.call(rbind, TABLE))
            }
        }
    } else {
        if(ngroups == 1L) {
            out <- TABLES[[1]]
        } else {
            out <- TABLES
        }
    }

    out
}


# low-level function to compute expected proportions per
lav_pairwise_tables_pi <- function(object) {

    stopifnot(class(object) == "lavaan", object@Model@categorical)
    ngroups <- object@Data@ngroups
    ov.types <- object@Data@ov$type

    th.idx <- object@Model@th.idx
    num.idx <- object@Model@num.idx
    Sigma.hat <- computeSigmaHat(object@Model)
    TH <- computeTH(object@Model)

    PI <- vector("list", length=ngroups)
    for(g in 1:ngroups) {
        Sigmahat <- Sigma.hat[[g]]
        cors <- Sigmahat[lower.tri(Sigmahat)]
        if(any(abs(cors) > 1)) {
            warning("lavaan WARNING: some model-implied correlations are larger than 1.0")
        }
        nvar <- nrow(Sigmahat)

        # shortcut for all ordered - tablewise
        if(all(ov.types == "ordered") && !is.null(object@Cache[[g]]$LONG)) {
            #FREQ.OBS <- c(FREQ.OBS, object@Cache[[g]]$bifreq)
            LONG2 <- LongVecTH.Rho(no.x               = nvar,
                                   all.thres          = TH[[g]],
                                   index.var.of.thres = th.idx[[g]],
                                   rho.xixj           = cors)
            # get expected probability per table, per pair
            PI[[g]] <- pairwiseExpProbVec(ind.vec = object@Cache[[g]]$LONG, 
                                          th.rho.vec=LONG2)
        } else {
            PI.group <- integer(0)
            # order! first i, then j, vec(table)!
            for(i in seq_len(nvar-1L)) {
                for(j in (i+1L):nvar) {
                    if(ov.types[i] == "ordered" && ov.types[j] == "ordered") {
                        PI.table <- pc_PI(rho   = Sigmahat[i,j],
                                          th.y1 = TH[[g]][ th.idx[[g]] == i ],
                                          th.y2 = TH[[g]][ th.idx[[g]] == j ])
                        PI.group <- c(PI.group, vec(PI.table))
                    }
                }
            }
            PI[[g]] <- PI.group
        }

    } # g
    
    PI
}
