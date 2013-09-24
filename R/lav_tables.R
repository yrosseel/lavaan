# construct 1D or 2D frequency tables
# YR. 10 April 2013
# Notes:
# - we do NOT make a distinction here between unordered and ordered categorical
#   variables
# - object can be a matrix (most likely with integers), a full data frame, 
#   a fitted lavaan object, or a lavData object
# - 11 May 2013: added collapse=TRUE, min.std.resid options (suggested
#   by Myrsini Katsikatsou
# - 11 June 2013: added dimension, to get one-way and two-way (three-way?)
#   tables
# - 20 Sept 2013: allow for sample-based or model-based cell probabilities
#                 re-organize/re-name to provide a more consistent interface
#                 rows in the output can be either: cells, tables or patterns

lavTables <- function(object,
                      # what type of table?
                      dimension    = 2L, 
                      type         = "cells",
                      # if raw data, additional attributes
                      categorical  = NULL,
                      group        = NULL,
                      # which statistics / fit indices?
                      statistic    = "default",
                      # pvalues for statistics?
                      p.value      = FALSE,
                      # output format
                      showAsMatrix = FALSE) {

    # check input
    if(! (dimension == 1L || dimension == 2L) ) {
        stop("lavaan ERROR: dimension must be 1 or 2 for one-way or two-way tables")
    }
    stopifnot(type %in% c("cells", "table", "pattern"))
    statistic  <- toupper(statistic)

    # showAsMatrix
    if(showAsMatrix) {
       stopifnot(dimension == 2L, length(statistic) == 1L, type == "table")
    }

    # extract or create lavdata
    lavdata <- lavData(object, ordered = categorical, group = group)

    # identify 'categorical' variables
    cat.idx <- which(lavdata@ov$type %in% c("ordered","factor"))

    # is 'object' a lavaan object?
    lavobject <- NULL
    if(inherits(object, "lavaan")) {
        lavobject <- object
    }

    if(dimension == 1L) {
        # only cells
        if(length(cat.idx) == 0L) {
            stop("lavaan ERROR: no categorical variables are found")
        }

        # first, basic freq tables
        out <- lav_tables_oneway(lavdata = lavdata)

        # default statistics
        if(inherits(object, "lavaan")) {
            if(length(statistic) == 1L && statistic == "DEFAULT") {
                statistic <- c("TH", "TH.H0")
            } else {
                stopifnot(statistic %in% c("TH", "TH.H0", "LR", "GF"))
            }
            
            if(any(c("TH", "TH.H0", "LR", "GF") %in% statistic)) {

                # sample based
                th <- unlist(lapply(1:lavdata@ngroups, function(x) {
                      unname(unlist(tapply(lavobject@SampleStats@th[[x]],
                                    INDEX=lavobject@SampleStats@th.idx[[x]],
                                    function(x) c(x,Inf)))) }))
                # overwrite obs.prop
                # NOTE: if we have exogenous variables, obs.prop will NOT
                #       correspond with qnorm(th)
                out$obs.prop <- unname(unlist(tapply(th, INDEX=out$id,   
                                   FUN=function(x) (pnorm(c(x,Inf)) -    
                                      pnorm(c(-Inf,x)))[-(length(x)+1)]  )))

                if("TH" %in% statistic) {
                    out$th <- th
                }

                # model based
                th.h0 <- unlist(lapply(1:lavdata@ngroups, function(x) {
                        unname(unlist(tapply(lavobject@Fit@TH[[x]],
                                      INDEX=lavobject@SampleStats@th.idx[[x]],
                                      function(x) c(x,Inf)))) }))

                est.prop <- unname(unlist(tapply(th.h0, INDEX=out$id, 
                                FUN=function(x) (pnorm(c(x,Inf)) - 
                                    pnorm(c(-Inf,x)))[-(length(x)+1)]  )))
                out$est.prop <- est.prop
         
                if("TH.H0" %in% statistic) {
                    out$th.h0 <- th.h0
                }
                if("LR" %in% statistic) {
                    out$LR <- lav_tables_stat_LR(out$obs.prop, out$est.prop,
                                                 out$nobs)
                }
                if("GF" %in% statistic) {
                    out$GF <- lav_tables_stat_GF(out$obs.prop, out$est.prop,
                                                 out$nobs)
                }
            }
        } else {
            if(length(statistic) == 1L && statistic == "DEFAULT") {
                # if data, none by default
                statistic <- character(0L)
            } else {
                stopifnot(statistic %in% c("TH"))
            }

            if("TH" %in% statistic) {
                out$th <- unlist(tapply(out$obs.prop, INDEX=out$id,
                                 FUN=function(x) qnorm(cumsum(x))))
            }
        }

    } else if(dimension == 2L) {
        if(type == "cells") {
            if(length(cat.idx) == 0L) {
                stop("lavaan ERROR: no categorical variables are found")
            } else if(length(cat.idx) == 1L) {
                stop("lavaan ERROR: at least two categorical variables are needed")
            }

            # default statistics
            if(inherits(object, "lavaan")) {
                if(length(statistic) == 1L && statistic == "DEFAULT") {
                    statistic <- c("GF", "LR")
                } else {
                    stopifnot(statistic %in% c("GF","LR","GF.BVN","LR.BVN"))
                }
            } else {
                if(length(statistic) == 1L && statistic == "DEFAULT") {
                    # if data, none by default
                    statistic <- character(0L)
                } else {
                    stopifnot(statistic %in% c("GF.BVN","LR.BVN"))
                }
            }

            out <- lav_tables_pairwise_cells(lavobject = lavobject, 
                                             lavdata   = lavdata,
                                             statistic = statistic)
        } else if(type == "table" || type == "tables") {
            # default statistics
            if(inherits(object, "lavaan")) {
                if(length(statistic) == 1L && statistic == "DEFAULT") {
                    statistic <- c("LR")
                } else {
                    stopifnot(statistic %in% c("GF","LR","GF.BVN","LR.BVN",
                                               "RMSEA.BVN", "RMSEA"))
                }
            } else {
                if(length(statistic) == 1L && statistic == "DEFAULT") {
                    # if data, none by default
                    statistic <- character(0L)
                } else {
                    stopifnot(statistic %in% c("GF.BVN","LR.BVN","RMSEA.BVN"))
                }
            }

            out <- lav_tables_pairwise_table(lavobject = lavobject, 
                                             lavdata   = lavdata,
                                             statistic = statistic,
                                             p.value   = p.value)
        } else if(type == "pattern") {
            stop("not implemented yet")
        }
    }

    # showAsMatrix?
    if(showAsMatrix) {
        # determine column we need
        NAMES <- names(out)
        stat.idx <- which(NAMES %in% c("LR",       "LR.BVN",
                                       "GF",       "GF.BVN",
                                       "RMSEA", "RMSEA.BVN"))
        if(length(stat.idx) == 0) {
            stop("lavaan ERROR: no statistic found in table for showAsMatrix")
        } else if(length(stat.idx) > 1) {
            stop("lavaan ERROR: more than one statistic found for showAsMatrix")
        } 
        OUT <- vector("list", length=lavdata@ngroups)
        for(g in 1:lavdata@ngroups) {
            if(lavdata@ngroups == 1L) { # no group column
                STAT <- out[[stat.idx]]
                RN <- unique(out$lhs)
            } else {
                STAT <- out[[stat.idx]][ out$group == g ]
                RN <- unique(out$lhs[ out$group == g ])
            }
            OUT[[g]] <- getCov(STAT, lower = FALSE, names = RN)
            class(OUT[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
        if(lavdata@ngroups > 1L) {
            names(OUT) <- lavdata@group.label
            out <- OUT
        } else {
            out <- OUT[[1]]      
        }
    } else {
        class(out) <- c("lavaan.data.frame", "data.frame")
    }

    out
}

# shortcut
lavTablesFit <- function(object,
                         # what type of table?
                         dimension    = 2L,
                         # if raw data, additional attributes
                         categorical  = NULL,
                         group        = NULL,
                         # which statistics / fit indices?
                         statistic    = "default",
                         # pvalues for statistics?
                         p.value      = FALSE,
                         # output format
                         showAsMatrix = FALSE) {

    lavTables(object = object, dimension = dimension, 
              categorical = categorical, group = group,
              statistic = statistic, p.value = p.value, 
              showAsMatrix = showAsMatrix)
}

lavTables1D <- function(object,
                        # if raw data, additional attributes
                        categorical  = NULL,
                        group        = NULL,
                        # which statistics / fit indices?
                        statistic    = "default") {

    lavTables(object = object, dimension = 1L,
              categorical = categorical, group = group,
              statistic = statistic, p.value = FALSE,
              showAsMatrix = FALSE)
}


# pairwise tables, rows = table cells
lav_tables_pairwise_cells <- function(lavobject = NULL, lavdata = NULL, 
                                      statistic = character(0L)) {

    # initial table, observed cell frequencies
    out <- lav_tables_pairwise_freq_cell(lavdata = lavdata, 
                                         as.data.frame. = TRUE)
    out$obs.prop <- out$obs.freq/out$nobs
    if(any(c("GF.BVN", "LR.BVN") %in% statistic)) {
        PI <- lav_tables_pairwise_sample_pi(lavobject = lavobject,
                                            lavdata   = lavdata)
        out$est.prop.bvn <- unlist(PI)
        if("LR.BVN" %in% statistic) {
            out$LR.BVN <- lav_tables_stat_LR(out$obs.prop, out$est.prop.bvn, 
                                             out$nobs)
        }
        if("GF.BVN" %in% statistic) {
            out$GF.BVN <- lav_tables_stat_GF(out$obs.prop, out$est.prop.bvn, 
                                             out$nobs)
        }
    }
    if(any(c("GF", "LR") %in% statistic)) {
        PI <- lav_tables_pairwise_model_pi(lavobject = lavobject)
        out$est.prop <- unlist(PI)
        if("LR" %in% statistic) {
            out$LR <- lav_tables_stat_LR(out$obs.prop, out$est.prop, 
                                         out$nobs)
        }
        if("GF" %in% statistic) {
            out$GF <- lav_tables_stat_GF(out$obs.prop, out$est.prop,
                                         out$nobs)
        }
    }

    out
}

# LR statistic
lav_tables_stat_LR <- function(obs.prop = NULL, est.prop = NULL, nobs = NULL) {
    # not defined if out$obs.prop is (close to) zero
    # we use freq=0.5 for these empty cells
    zero.idx <- which(obs.prop < .Machine$double.eps)
    obs.prop[zero.idx] <- 0.5/nobs[zero.idx]
    LR <- 2*nobs*(obs.prop*log(obs.prop/est.prop))
    LR
}

# GF (aka X2) statistic
lav_tables_stat_GF <- function(obs.prop = NULL, est.prop = NULL, nobs = NULL) {
    GF <- nobs*(obs.prop-est.prop)^2/est.prop
    GF
}

# pairwise tables, rows = tables
lav_tables_pairwise_table <- function(lavobject = NULL, lavdata = NULL,
                                      statistic = character(0L),
                                      p.value = FALSE) {

    # identify 'categorical' variables
    #cat.idx <- which(lavdata@ov$type %in% c("ordered","factor"))

    # pairwise tables
    #pairwise.tables <- utils::combn(vartable$name[cat.idx], m=2L)
    #pairwise.tables <- rbind(seq_len(ncol(pairwise.tables)), 
    #                         pairwise.tables)
    #ntables <- ncol(pairwise.tables)

    # initial table, observed cell frequencies
    #out <- as.data.frame(t(pairwise.tables))
    #names(out) <- c("id", "lhs", "rhs")

    # collapse approach
    stat.cell <- character(0)
    if("GF" %in% statistic) {
        stat.cell <- c(stat.cell, "GF")
    }
    if("LR" %in% statistic || "RMSEA" %in% statistic) {
        stat.cell <- c(stat.cell, "LR")
    }
    if("GF.BVN" %in% statistic) {
        stat.cell <- c(stat.cell, "GF.BVN")
    }
    if("LR.BVN" %in% statistic || "RMSEA.BVN" %in% statistic) {
        stat.cell <- c(stat.cell, "LR.BVN")
    }

    # get table with table cells
    out.cell <- lav_tables_pairwise_cells(lavobject = lavobject,
                                          lavdata   = lavdata,
                                          statistic = stat.cell)
    # only 1 row per table
    row.idx <- which(!duplicated(out.cell$id))
    out <- out.cell[row.idx,c("lhs","rhs","nobs"),drop=FALSE]

    # df
    nrow <- tapply(out.cell$row, INDEX=out.cell$id, FUN=max)
    ncol <- tapply(out.cell$col, INDEX=out.cell$id, FUN=max)
    out$df <- nrow*ncol - nrow - ncol

    # GF
    if("GF" %in% statistic) {
        out$GF <- tapply(out.cell$GF, INDEX=out.cell$id, FUN=sum)
        if(p.value) {
            out$GF.pval <- pchisq(out$GF, df=out$df, lower.tail=FALSE)
        }
    }
    if("GF.BVN" %in% statistic) {
        out$GF.BVN <- tapply(out.cell$GF.BVN, INDEX=out.cell$id, FUN=sum)
        if(p.value) {
            out$GF.BVN.pval <- pchisq(out$GF.BVN, df=out$df, lower.tail=FALSE)
        }
    }
 
    # LR
    if("LR" %in% statistic) {
        out$LR <- tapply(out.cell$LR, INDEX=out.cell$id, FUN=sum)
        if(p.value) {
            out$LR.pval <- pchisq(out$LR, df=out$df, lower.tail=FALSE)
        }
    }
    if("LR.BVN" %in% statistic) {
        out$LR.BVN <- tapply(out.cell$LR.BVN, INDEX=out.cell$id, FUN=sum)
        if(p.value) {
            out$LR.BVN.pval <- pchisq(out$LR.BVN, df=out$df, lower.tail=FALSE)
        }
    }

    if("RMSEA" %in% statistic) {
        LR <- tapply(out.cell$LR, INDEX=out.cell$id, FUN=sum)
        out$RMSEA <- sqrt( pmax(0, (LR - out$df)/ (2*out$nobs*out$df) ) )
        if(p.value) {
            # note: MUST use 1 - pchisq (instead of lower.tail = FALSE)
            # because for ncp > 80, routine only computes lower tail
            out$RMSEA.pval <- 1.0 - pchisq(LR,
                                           ncp = 0.1^2*out$nobs*out$df,
                                           df=out$df, lower.tail = TRUE)
        }
    }
    if("RMSEA.BVN" %in% statistic) {
        LR <- tapply(out.cell$LR.BVN, INDEX=out.cell$id, FUN=sum)
        out$RMSEA.BVN <- sqrt( pmax(0, (LR - out$df)/ (2*out$nobs*out$df) ) )
        if(p.value) {
            # note: MUST use 1 - pchisq (instead of lower.tail = FALSE)
            # because for ncp > 80, routine only computes lower tail
            out$RMSEA.BVN.pval <- 1.0 - pchisq(LR,
                                               ncp = 0.1^2*out$nobs*out$df,
                                               df=out$df, lower.tail = TRUE)
        }
    }
    
    out
}


lav_tables_oneway <- function(lavdata = NULL) {

    # shortcuts
    vartable <- lavdata@ov
    X        <- lavdata@X

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
                nvar <- length(lavdata@ov.names[[g]])
                id <- (g-1)*nvar + x

                # compute observed frequencies
                FREQ <- tabulate(X[[g]][,idx])

                list(   id = rep.int(id, ncell),
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
        } else {
            out <- rbind(out, do.call(rbind, TABLE))
        }
    }
    if(g == 1) {
        # remove group column 
        out$group <- NULL
    }


    class(out) <- c("lavaan.data.frame", "data.frame")
    out
}

# YR: - 20 Sept 2013, argument PI must be given
lav_tables_pairwise_OLD <- function(lavdata = NULL,
                                PI      = NULL, # fitted cell probabilities
                                min.std.resid = 0.0,
                                average = FALSE, collapse = FALSE,
                                method="GF") {

    # check input
    showAsMatrix <- FALSE
    if(! (is.logical(collapse) && !collapse) ) {
        average <- TRUE
        if(is.character(collapse) && tolower(collapse) == "matrix") {
            collapse <- TRUE
            showAsMatrix <- TRUE
        }
    }

    out <- lav_tables_pairwise_freq_cell(lavdata = lavdata,
                                    as.data.frame. = TRUE)

    if(is.null(PI)) {
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
            table.p.value <- unlist(lapply(tableChisqTest,   
                function(x) as.numeric(x$p.value)))
            out$chisq <- unlist(lapply(seq_len(ntables),
                    function(x) rep(table.chisq[x], each=ncells[x])))
            out$df <- unlist(lapply(seq_len(ntables),
                    function(x) rep(table.df[x], each=ncells[x])))
            out$p.value <- unlist(lapply(seq_len(ntables),
                    function(x) rep(table.p.value[x], each=ncells[x])))
        }
    } else {
        # add predicted frequencies, fit indices
        obs.prop <- out$obs.freq/out$nobs
        est.prop <- PI
        out$est.freq <- est.prop * out$nobs
 
        # Joreskog & Moustaki equation 34/35
        if(method == "GF") {
            out$GF <- out$nobs*(obs.prop-est.prop)^2/est.prop
        } else if(method == "LR") {
            # not defined if out$obs.prop is (close to) zero
            # we use freq=0.5 for these empty cells
            zero.idx <- which(obs.prop < .Machine$double.eps)
            obs.prop[zero.idx] <- 0.5/out$nobs[zero.idx]
            out$LR <- 2*out$nobs*(obs.prop*log(obs.prop/est.prop))
        }

        if(average) {
            # fit per table
            table.average <- tapply(out$std.resid, INDEX=out$id, FUN=mean)
            table.sum <- tapply(out$std.resid, INDEX=out$id, FUN=sum)
            table.percAboveMin <- tapply(out$std.resid, INDEX=out$id,
                FUN=function(x) {sum(x > min.std.resid)/length(x)})
            table.numAboveMin <- tapply(out$std.resid, INDEX=out$id,
                FUN=function(x) {sum(x > min.std.resid)})
            ntables <- length(table.average)
            ncells <- tabulate(out$id)
            nrow <- tapply(out$row, INDEX=out$id, FUN=max)
            ncol <- tapply(out$col, INDEX=out$id, FUN=max)

            out$nrow <- unlist(lapply(seq_len(ntables),
                function(x) rep(nrow[x], each=ncells[x])))
            out$ncol <- unlist(lapply(seq_len(ntables),
                function(x) rep(ncol[x], each=ncells[x])))
            #out$ncells <- unlist(lapply(seq_len(ntables),
            #    function(x) rep(ncells[x], each=ncells[x])))
            out$str.average <- unlist(lapply(seq_len(ntables),
                function(x) rep(table.average[x], each=ncells[x])))
            out$str.sum <- unlist(lapply(seq_len(ntables),
                function(x) rep(table.sum[x], each=ncells[x])))
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

    if(collapse) {
        # only 1 row per table
        row.idx <- which(!duplicated(out$id))
        out <- out[row.idx,,drop=FALSE]

        # remove some cell-specific columns
        out$row <- NULL; out$col <- NULL
        out$obs.freq <- NULL; out$est.freq <- NULL
    }

    if(collapse && showAsMatrix) {
        tmp <- out
        RN <- unique(out$rhs)
        if(inherits(object, "lavaan")) {
            out <- lavaan::getCov(out$str.average, lower=FALSE, 
                                  names=unique(out$lhs))
            rownames(out) <- RN
        } else {
            out <- lavaan::getCov(out$p.value, lower=FALSE, 
                                  names=unique(out$lhs))
            rownames(out) <- RN
        }
        class(out) <- c("lavaan.matrix.symmetric", "matrix")
    } else {
        class(out) <- c("lavaan.data.frame", "data.frame")
    }
    out
}

# compute pairwise (two-way) frequency tables
lav_tables_pairwise_freq_cell <- function(lavdata = NULL, 
                                          as.data.frame. = TRUE) {

    # shortcuts
    vartable <- as.data.frame(lavdata@ov, stringsAsFactors = FALSE)
    X        <- lavdata@X
    ov.names <- lavdata@ov.names
    ngroups  <- lavdata@ngroups

    # identify 'categorical' variables
    cat.idx <- which(vartable$type %in% c("ordered","factor"))

    # do we have any categorical variables?
    if(length(cat.idx) == 0L) {
        stop("lavaan ERROR: no categorical variables are found")
    } else if(length(cat.idx) == 1L) {
        stop("lavaan ERROR: at least two categorical variables are needed")
    }

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
            } else {
                out <- rbind(out, do.call(rbind, TABLE))
            }
        }
        if(g == 1) {
            # remove group column 
            out$group <- NULL
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


# low-level function to compute expected proportions per cell
# object
lav_tables_pairwise_model_pi <- function(lavobject = NULL) {

    stopifnot(lavobject@Model@categorical)

    # shortcuts
    ngroups   <- lavobject@Data@ngroups
    ov.types  <- lavobject@Data@ov$type
    th.idx    <- lavobject@Model@th.idx
    num.idx   <- lavobject@Model@num.idx
    Sigma.hat <- lavobject@Fit@Sigma.hat
    TH        <- lavobject@Fit@TH

    PI <- vector("list", length=ngroups)
    for(g in 1:ngroups) {
        Sigmahat <- Sigma.hat[[g]]
        cors <- Sigmahat[lower.tri(Sigmahat)]
        if(any(abs(cors) > 1)) {
            warning("lavaan WARNING: some model-implied correlations are larger than 1.0")
        }
        nvar <- nrow(Sigmahat)

        # shortcut for all ordered - tablewise
        if(all(ov.types == "ordered") && !is.null(lavobject@Cache[[g]]$LONG)) {
            #FREQ.OBS <- c(FREQ.OBS, lavobject@Cache[[g]]$bifreq)
            LONG2 <- LongVecTH.Rho(no.x               = nvar,
                                   all.thres          = TH[[g]],
                                   index.var.of.thres = th.idx[[g]],
                                   rho.xixj           = cors)
            # get expected probability per table, per pair
            PI[[g]] <- pairwiseExpProbVec(ind.vec = lavobject@Cache[[g]]$LONG, 
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

# low-level function to compute expected proportions per cell
# using sample-based correlations + thresholds
#
# object can be either lavData or lavaan class
lav_tables_pairwise_sample_pi <- function(lavobject = NULL, lavdata = NULL) {

    # get COR, TH and th.idx
    if(!is.null(lavobject)) {
        COR    <- lavobject@SampleStats@cov
        TH     <- lavobject@SampleStats@th
        th.idx <- lavobject@SampleStats@th.idx
    } else if(!is.null(lavdata)) {
        COR <- lav_cor(lav.data = lavdata, WLS.W = FALSE, details = TRUE,
                       labels = FALSE, verbose = FALSE)
        # relist
        if(!is.list(COR)) {
            COR <- list(COR)
        }
        TH <- lapply(lapply(COR, attr, "TH"), unlist)
        th.idx <- lapply(lapply(COR, attr, "TH.IDX"), unlist)
    } else {
        stop("lavaan ERROR: both lavobject and lavdata are NULL")
    }

    lav_tables_pairwise_sample_pi_cor(COR = COR, TH = TH,
                                      th.idx = th.idx)
}

# low-level function to compute expected proportions per cell
lav_tables_pairwise_sample_pi_cor <- function(COR = NULL, TH = NULL, 
                                              th.idx = NULL) {

    ngroups <- length(COR)

    PI <- vector("list", length=ngroups)
    for(g in 1:ngroups) {
        Sigmahat <- COR[[g]]
        cors <- Sigmahat[lower.tri(Sigmahat)]
        if(any(abs(cors) > 1)) {
            warning("lavaan WARNING: some model-implied correlations are larger than 1.0")
        }
        nvar <- nrow(Sigmahat)

        # reconstruct ov.types
        ov.types <- rep("numeric", nvar)
        ord.idx <- unique(which(th.idx[[g]] > 0))
        ov.types[ord.idx] <- "ordered"

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
    } # g
    
    PI
}
