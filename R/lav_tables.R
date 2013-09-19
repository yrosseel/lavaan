# construct 1D or 2D frequency tables
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
                      model.based = TRUE,
                      std.resid=TRUE, min.std.resid = 0.0, method = "GF",
                      average = FALSE, collapse = FALSE) {
    
    # extract data
    lav_data <- lav_data_extract(object = object, categorical = categorial)
    vartable <- lav_data$vartable
    X        <- lav_data$X
    ov.names <- lav_data$ov.names

    if(dimension == 1L) {
        out <- lav_oneway_tables(object = object, 
                   vartable = vartable, X = X, ov.names = ov.names)
    } else if(dimension == 2L) {
        out <- lav_pairwise_tables(object = object,
                   vartable = vartable, X = X, ov.names = ov.names,
                   average = average, model.based = model.based,
                   std.resid = std.resid, collapse = collapse, method = method)
    } else {
        stop("lavaan ERROR: dimension must be 1 or 2 of one-way or two-way tables")
    }

    out
}

# always collapse, 1 statistic per table
lavTablesFit <- function(object, dimension=2L, statistic="LR", 
                         model.based = TRUE,
                         p.value = FALSE, showAsMatrix = FALSE) {

    statistic <- toupper(statistic)

    # check object class
    if(!class(object) %in% c("lavaan")) {
        stop("lavaan ERROR: object must either an object of class lavaan")
    }

    vartable <- object@Data@ov
    X <- object@Data@X
    ov.names <- object@Data@ov.names

    if(dimension == 1L) {
        out <- lav_oneway_tables(object = object,
                   vartable = vartable, X = X, ov.names = ov.names)
    } else if(dimension == 2L) {
        if(statistic == "LR") {
            std.resid <- TRUE; method <- "LR"
        } else if(statistic == "GF") {
            std.resid <- TRUE; method <- "GF"
        } else if(statistic == "RMSEA") {
            std.resid <- TRUE; method <- "LR"
        } else {
            stop("lavaan ERROR: can not handle statistic ", statistic)
        }
        out <- lav_pairwise_tables(object = object,
                   vartable = vartable, X = X, ov.names = ov.names,
                   average = TRUE, model.based = model.based,
                   std.resid = std.resid, method = method, collapse = TRUE)

        # df
        out$df <- out$nrow*out$ncol - out$nrow - out$ncol

        # stat
        if(statistic == "LR") {
            out$LR <- out$str.sum
            STAT <- out$LR
            if(p.value) {
                out$p.value <- pchisq(out$LR, df=out$df, lower.tail = FALSE)
            }
        } else if(statistic == "GF") {
            out$GF <- out$str.sum
            STAT <- out$GF
            if(p.value) {
                out$p.value <- pchisq(out$GF, df=out$df, lower.tail = FALSE)
            }
        } else if(statistic == "RMSEA") {
            out$RMSEA <- sqrt( pmax(0, (out$str.sum - out$df)/
                                       (2*out$nobs*out$df)     ) )
            STAT <- out$RMSEA

            if(p.value) {
                # note: MUST use 1 - pchisq (instead of lower.tail = FALSE)
                # because for ncp > 80, routine only computes lower tail
                out$p.value <- 1.0 - pchisq(out$str.sum, 
                                            ncp = 0.1^2*out$nobs*out$df, 
                                            df=out$df, lower.tail = TRUE)
            }
        }

        # remove columns
        out$id <- out$nobs <- out$nrow <- out$ncol <- NULL
        out$str.average <- out$str.sum <- NULL
        out$str.min <- out$str.plarge <- out$str.nlarge <- NULL

        if(showAsMatrix) {
            RN <- unique(out$rhs)
            out <- lavaan::getCov(STAT, lower=FALSE,
                                  names=unique(out$lhs))
            rownames(out) <- RN
            class(out) <- c("lavaan.matrix.symmetric", "matrix")
        }

    } else {
        stop("lavaan ERROR: dimension must be 1 or 2 of one-way or two-way tables")
    }

    out
}

# Mariska Barendse Cp statistic
lavTablesFitCp <- function(object, alpha = 0.05) {

    # check object class
    if(!class(object) %in% c("lavaan")) {
        stop("lavaan ERROR: object must either an object of class lavaan")
    }

    vartable <- object@Data@ov
    X <- object@Data@X
    ov.names <- object@Data@ov.names

    std.resid <- TRUE; method <- "LR"
    out <- lav_pairwise_tables(object = object,
               vartable = vartable, X = X, ov.names = ov.names,
               average = TRUE,
               std.resid = std.resid, method = method, collapse = TRUE)
    out$LR <- out$str.sum
    out$id <- out$nobs <- NULL
    out$str.average <- out$str.sum <- NULL
    out$str.min <- out$str.plarge <- out$str.nlarge <- NULL

    # df 
    out$df <- out$nrow*out$ncol - out$nrow - out$ncol

    # p-value
    ntests <- length(out$p.value)
    out$p.value.adj <- pchisq(out$LR, 
                              df=out$df,
                              lower.tail = FALSE) * ntests

    # Bonferroni alpha
    #ntests <- length(out$p.value)
    #out$alpha.star <- rep(alpha / ntests, length(out$p.value))

    out
}

lavTablesFitCpMax <- function(object, alpha = 0.05) {
    out <- lavTablesFitCp(object = object, alpha = alpha)

    # find largest LR
    max.idx <- which(out$LR == max(out$LR))

    list(LR=out$LR[max.idx], df=out$df[max.idx], 
         p.value=out$p.value[max.idx], alpha.star=out$alpha.star[max.idx],
         p.value.Bonferroni=out$p.value[max.idx]*length(out$LR))
}

# Mariska Barendse CF statistic
lavTablesFitCF <- function(object) {

    # check object class
    if(!class(object) %in% c("lavaan")) {
        stop("lavaan ERROR: object must either an object of class lavaan")
    }

    ngroups <- length( object@Data@X )

    CF.group <- numeric(ngroups)
    DF.group <- numeric(ngroups)

    for(g in 1:ngroups) {
        logLik.group <- estimator.FML(Sigma.hat = object@Fit@Sigma.hat[[g]],
                                      TH        = object@Fit@TH[[g]],
                                      th.idx    = object@Model@th.idx[[g]],
                                      num.idx   = object@Model@num.idx[[g]],
                                      X         = object@Model@X[[g]],
                                      cache     = object@Cache[[g]])

        freq <- as.numeric( rownames(object@Data@Rp[[g]]$pat) )
        CF.group[g] <- 2*logLik.group + 2*sum(freq*log(freq/sum(freq)))


        # ord var in this group
        ov.ord <- object@pta$vnames$ov.ord[[g]]
        ov.idx <- which(ov.ord %in% object@Data@ov$name)
        ov.nlev <- object@Data@ov$nlev[ ov.idx ]

        DF.group[g] <- prod(ov.nlev) - object@Fit@npar - 1L
    }

    # check for negative values
    CF.group[CF.group < 0] <- 0.0

    # global test statistic
    CF <- sum(CF.group)

    attr(CF, "CF.group") <- CF.group
    attr(CF, "DF.group") <- DF.group

    CF
}

lavTablesFitCF.h1 <- function(object) {

    # check object class
    if(!class(object) %in% c("lavaan")) {
        stop("lavaan ERROR: object must either an object of class lavaan")
    }

    ngroups <- length( object@Data@X )

    CF.group <- numeric(ngroups)
    DF.group <- numeric(ngroups)

    for(g in 1:ngroups) {
        logLik.group <- estimator.FML(Sigma.hat = object@SampleStats@cov[[g]],
                                      TH        = object@SampleStats@th[[g]],
                                      th.idx    = object@Model@th.idx[[g]],
                                      num.idx   = object@Model@num.idx[[g]],
                                      X         = object@Model@X[[g]],
                                      cache     = object@Cache[[g]])

        freq <- as.numeric( rownames(object@Data@Rp[[g]]$pat) )
        CF.group[g] <- 2*logLik.group + 2*sum(freq*log(freq/sum(freq)))


        # ord var in this group
        ov.ord <- object@pta$vnames$ov.ord[[g]]
        ov.idx <- which(ov.ord %in% object@Data@ov$name)
        ov.nlev <- object@Data@ov$nlev[ ov.idx ]

        DF.group[g] <- prod(ov.nlev) - getNDAT(object@ParTable) - 1L
    }

    # check for negative values
    CF.group[CF.group < 0] <- 0.0

    # global test statistic
    CF <- sum(CF.group)

    attr(CF, "CF.group") <- CF.group
    attr(CF, "DF.group") <- DF.group

    CF
}

lavTablesFitCM <- function(object) {

   CF.h0 <- lavTablesFitCF(object)
   CF.h1 <- lavTablesFitCF.h1(object)

   CF.h0.group <- attr(CF.h0, "CF.group")
   CF.h1.group <- attr(CF.h1, "CF.group")
   DF.h0.group <- attr(CF.h0, "DF.group")
   DF.h1.group <- attr(CF.h1, "DF.group")

   attributes(CF.h0) <- NULL
   attributes(CF.h1) <- NULL

   CM <- CF.h0 - CF.h1
   attr(CM, "CM.group") <- CF.h0.group - CF.h1.group
   attr(CM, "DF.group") <- DF.h0.group - DF.h1.group

   CM
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

# YR: - 23 Aug 2013: added Mariska Barendse Cp fit index
lav_pairwise_tables <- function(object, 
                                vartable=NULL, X=NULL, ov.names=NULL,
                                std.resid=TRUE, min.std.resid = 0.0,
                                average = FALSE, collapse = FALSE,
                                model.based = TRUE,
                                method="GF") {

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
            table.p.value <- unlist(lapply(tableChisqTest,   
                function(x) as.numeric(x$p.value)))
            out$chisq <- unlist(lapply(seq_len(ntables),
                    function(x) rep(table.chisq[x], each=ncells[x])))
            out$df <- unlist(lapply(seq_len(ntables),
                    function(x) rep(table.df[x], each=ncells[x])))
            out$p.value <- unlist(lapply(seq_len(ntables),
                    function(x) rep(table.p.value[x], each=ncells[x])))
        }
    } else if(inherits(object, "lavaan")) {
        # add predicted frequencies, fit indices

        #out$pi <- unlist(lav_pairwise_tables_pi(object))
        #out$freq.est <- out$pi * out$nobs
        obs.prop <- out$obs.freq/out$nobs
        if(model.based) {
            est.prop <- unlist(lav_pairwise_tables_model_pi(object))
        } else {
            COR <- object@SampleStats@cov
            TH  <- object@SampleStats@th
            th.idx <- object@SampleStats@th.idx
            est.prop <- unlist(lav_pairwise_tables_sample_pi(COR = COR, TH = TH,
                                                             th.idx = th.idx))
        }
        out$est.freq <- est.prop * out$nobs
 
        # Joreskog & Moustaki equation 34/35
        if(std.resid) {
            if(method == "GF") {
                out$std.resid <- out$nobs*(obs.prop-est.prop)^2/est.prop
            } else if(method == "LR") {
                # not defined if out$obs.prop is (close to) zero
                # we use freq=0.5 for these empty cells
                zero.idx <- which(obs.prop < .Machine$double.eps)
                obs.prop[zero.idx] <- 0.5/out$nobs[zero.idx]
                out$std.resid <- 2*out$nobs*(obs.prop*log(obs.prop/est.prop))
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
lav_pairwise_tables_model_pi <- function(object) {

    stopifnot(class(object) == "lavaan", object@Model@categorical)
    ngroups <- object@Data@ngroups
    ov.types <- object@Data@ov$type

    th.idx <- object@Model@th.idx
    num.idx <- object@Model@num.idx

    #if(model.based) {
        Sigma.hat <- computeSigmaHat(object@Model)
        TH <- computeTH(object@Model)
    #} else {
    #    Sigma.hat <- object@SampleStats@cov
    #    TH <- object@SampleStats@th
    #}

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

# low-level function to compute expected proportions per cell
# no object!
lav_pairwise_tables_sample_pi <- function(COR = NULL, TH = NULL, 
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
