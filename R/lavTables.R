# construct pairwise frequency tables
# Notes:
# - we do NOT make a distinction here between unordered and ordered categorical
#   variables
# - object can be a matrix (most likely with integers), a full data frame, 
#   a fitted lavaan object, a varTable, or a lavData object
# YR. 10 April 2013
lavTables <- function(object, categorical=NULL, as.data.frame.=TRUE,
                      fit=TRUE, check=TRUE, fit.average=FALSE) {

    if(is.matrix(object)) {
        object <- as.data.frame(object, stringsAsFactors = FALSE)
        if(is.null(categorical)) {
            object[,] <- lapply(object, base::factor)
        }
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
    } else {
        stop("lavaan ERROR: object must either be a matrix or a data.frame, or an object of class lavaan or class lavData")
    }

    out <- lav_pairwise_tables_freq(vartable = vartable,
                                    X = X,
                                    ov.names = ov.names,
                                    as.data.frame. = as.data.frame.)

    # add predicted frequencies, fit indices
    if(inherits(object, "lavaan")) {
        #out$pi <- unlist(lav_pairwise_tables_pi(object))
        #out$freq.est <- out$pi * out$nobs
        prop <- out$freq/out$nobs
        pi <- unlist(lav_pairwise_tables_pi(object))
        out$freq.est <- unlist(lav_pairwise_tables_pi(object)) * out$nobs
 
        # Joreskog & Moustaki equation 35
        if(fit) {
            out$fit <- out$nobs*(prop-pi)^2/pi
            if(check) {
                out$check <- ifelse(out$fit > 4.0, "***", "")
            }

            if(fit.average) {
                # fit per table
                table.average <- tapply(out$fit, INDEX=out$id, FUN=mean)
                ntables <- length(table.average)
                ncells <- tabulate(out$id)
                out$fit.average <- unlist(lapply(seq_len(ntables), function(x)
                                                 rep(table.average[x], 
                                                     each=ncells[x])))
            }
        }
    }

    class(out) <- c("lavaan.data.frame", "data.frame")
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
                        op = rep.int("table", ncell), 
                       rhs = rep.int(x[2], ncell),
                     group = rep.int(g, ncell),
                      nobs = rep.int(sum(FREQ), ncell),
                       row = rep.int(seq_len(ncol), times=nrow),
                       col = rep(seq_len(nrow), each=ncol),
                      freq = vec(FREQ) # col by col!
                    )
            })
    }

    if(as.data.frame.) {
        for(g in 1:ngroups) {
            TABLE <- TABLES[[g]]
            TABLE <- lapply(TABLE, as.data.frame, stringAsFactors=FALSE)
            if(g == 1) {
                out <- do.call(rbind, TABLE)
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
