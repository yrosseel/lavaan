# construct 1D, 2D or pattern-based frequency tables
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
# - 20 Sept 2013: - allow for sample-based or model-based cell probabilities
#                   re-organize/re-name to provide a more consistent interface
#                   rows in the output can be either: cells, tables or patterns
#                 - dimension=0 equals type="pattern
#                 - collapse=TRUE is replaced by type="table"
#                 - changed names of statistics: std.resid is now GR.average
#                 - added many more statistics; some based on the model, some
#                   on the unrestricted model
# - 8 Nov 2013:   - skip empty cells for LR, instead of adding 0.5 to obs

lavTables <- function(object,
                      # what type of table?
                      dimension    = 2L, 
                      type         = "cells",
                      # if raw data, additional attributes
                      categorical  = NULL,
                      group        = NULL,
                      # which statistics / fit indices?
                      statistic    = "default",
                      LR.min       = 3.0,        # needed for LR.{p/n}large
                      GF.min       = 3.0,        # needed for GF.{p/n}large
                      # pvalues for statistics?
                      p.value      = FALSE,
                      # Bonferonni
                      # alpha.adj    = FALSE,
                      # output format
                      output       = "data.frame",
                      patternAsString = TRUE) {

    # check input
    if(! (dimension == 0L || dimension == 1L || dimension == 2L) ) {
        stop("lavaan ERROR: dimension must be 0, 1 or 2 for pattern, one-way or two-way tables")
    }
    stopifnot(type %in% c("cells", "table", "pattern"))
    if(type == "pattern") {
        dimension <- 0L
    }

    # extract or create lavdata
    lavdata <- lavData(object, ordered = categorical, group = group)

    # is 'object' a lavaan object?
    lavobject <- NULL
    if(inherits(object, "lavaan")) {
        lavobject <- object
    }

    # case 1: response patterns
    if(dimension == 0L) {
        out <- lav_tables_pattern(lavobject = lavobject, lavdata = lavdata,
                                  statistic = statistic, 
                                  patternAsString = patternAsString)
        # output format
        if(output == "data.frame") {
            class(out) <- c("lavaan.data.frame", "data.frame")
        } else {
            warning("lavaan WARNING: output option `", output, "' is not available; ignored.")
        }

    # case 2: one-way/univariate
    } else if(dimension == 1L) {
        out <- lav_tables_oneway(lavobject = lavobject, lavdata = lavdata,
                                 statistic = statistic)

        # output format
        if(output == "data.frame") {
            class(out) <- c("lavaan.data.frame", "data.frame")
        } else {
            warning("lavaan WARNING: output option `", output, "' is not available; ignored.")
        }

    # case 3a: two-way/pairwise/bivariate + cells
    } else if(dimension == 2L && type == "cells") {
        out <- lav_tables_pairwise_cells(lavobject = lavobject, 
                                         lavdata   = lavdata,
                                         statistic = statistic)
        # output format
        if(output == "data.frame") {
            class(out) <- c("lavaan.data.frame", "data.frame")
        } else if(output == "table") {
            out <- lav_tables_cells_format(out, lavdata = lavdata)
        } else {
            warning("lavaan WARNING: output option `", output, "' is not available; ignored.")           
        }

    #  case 3b: two-way/pairwise/bivariate + collapsed table
    } else if(dimension == 2L && (type == "table" || type == "tables")) {
        out <- lav_tables_pairwise_table(lavobject = lavobject, 
                                         lavdata   = lavdata,
                                         statistic = statistic,
                                         LR.min    = LR.min,
                                         GF.min    = GF.min,
                                         p.value   = p.value)
        # output format
        if(output == "data.frame") {
            class(out) <- c("lavaan.data.frame", "data.frame")
        } else if(output == "table") {
            out <- lav_tables_table_format(out, lavdata = lavdata, 
                                           lavobject = lavobject)
        } else {
            warning("lavaan WARNING: output option `", output, "' is not available; ignored.")
        }
    }

    out
}

# shortcut, always dim=2, type="cells"
lavTablesFit <- function(object,
                         # if raw data, additional attributes
                         categorical  = NULL,
                         group        = NULL,
                         # which statistics / fit indices?
                         statistic    = "default",
                         LR.min       = 3.0,
                         GF.min       = 3.0,
                         # pvalues for statistics?
                         p.value      = FALSE,
                         # output format
                         output       = "data.frame") {

    lavTables(object = object, dimension = 2L, type = "table",
              categorical = categorical, group = group,
              statistic = statistic, 
              LR.min = LR.min, GF.min = GF.min, p.value = p.value, 
              output = output, patternAsString = FALSE)
}

lavTables1D <- function(object,
                        # if raw data, additional attributes
                        categorical  = NULL,
                        group        = NULL,
                        # which statistics / fit indices?
                        statistic    = "default",
                        # output format
                        output = "data.frame") {

    lavTables(object = object, dimension = 1L,
              categorical = categorical, group = group,
              statistic = statistic, p.value = FALSE,
              output = output, patternAsString = FALSE)
}


lav_tables_pattern <- function(lavobject = NULL, lavdata = NULL,
                               statistic = NULL, patternAsString = TRUE) {

    # this only works if we have 'categorical' variables
    cat.idx <- which(lavdata@ov$type %in% c("ordered","factor"))
    if(length(cat.idx) == 0L) {
        stop("lavaan ERROR: no categorical variables are found")
    }

    # default statistics
    statistic <- toupper(statistic)
    if(!is.null(lavobject)) {
        if(length(statistic) == 1L && statistic == "DEFAULT") {
            statistic <- c("LR", "GF")
        } else {
            stopifnot(statistic %in% c("LR.BVN", "GF.BVN", "LR", "GF"))
        }  
    } else {
        if(length(statistic) == 1L && statistic == "DEFAULT") {
            # if data, none by default
            statistic <- character(0L)
        } else {
            stopifnot(statistic %in% c("LR.BVN", "GF.BVN"))
        }
    }

    # first, create basic table with response patterns
    for(g in 1:lavdata@ngroups) {
        pat <- lav_data_resppatterns(lavdata@X[[g]])$pat
        obs.freq <- as.integer( rownames(pat) )
        if(patternAsString) {
            pat <- data.frame(pattern = apply(pat, 1, paste, collapse=""),
                              stringsAsFactors = FALSE)
        } else {
            pat <- as.data.frame(pat, stringsAsFactors = FALSE)
            names(pat) <- lavdata@ov.names[[g]]
        }
        #pat$id <- 1:nrow(pat)
        if(lavdata@ngroups > 1L) {
            pat$group <- rep(g, nrow(pat))
        }
        pat$nobs <- rep(nrow(lavdata@X[[g]]), nrow(pat))
        pat$obs.freq <- obs.freq
        rownames(pat) <- NULL
        if(g == 1L) {
            out <- pat
        } else {
            out <- rbind(out, pat)
        }
    }

    out$obs.prop <- out$obs.freq/out$nobs

    if(any(c("GF.BVN", "LR.BVN") %in% statistic)) {
        # not a good statistic... we only have uni+bivariate information
        warning("lavaan WARNING: limited information used for thresholds and correlations; but GF/LR assumes full information")
        PI <- lav_tables_resp_pi(lavobject = lavobject, lavdata = lavdata,
                                 est = "h1")

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
        if(lavobject@Options$estimator == "FML") {
            # ok, nothing to say
        } else if(lavobject@Options$estimator %in% c("WLS","DWLS","PML")) {
            warning("lavaan WARNING: estimator ", lavobject@Options$estimator,
                    " is not using full information while est.prop is using full information")
        } else {
            stop("lavaan ERROR: estimator ", lavobject@Options$estimator,
                 " is not supported.")
        }

        PI <- lav_tables_resp_pi(lavobject = lavobject, lavdata = lavdata,
                                 est = "h0")

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
 
    # remove nobs?
    # out$nobs <- NULL

    out
}

# pairwise tables, rows = table cells
lav_tables_pairwise_cells <- function(lavobject = NULL, lavdata = NULL, 
                                      statistic = character(0L)) {

    # this only works if we have at least two 'categorical' variables
    cat.idx <- which(lavdata@ov$type %in% c("ordered","factor"))
    if(length(cat.idx) == 0L) {
        stop("lavaan ERROR: no categorical variables are found")
    } else if(length(cat.idx) == 1L) {
        stop("lavaan ERROR: at least two categorical variables are needed")
    }

    # default statistics
    statistic <- toupper(statistic)
    if(!is.null(lavobject)) {
        if(length(statistic) == 1L && statistic == "DEFAULT") {
            statistic <- c("LR")
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
    zero.idx <- which(obs.prop < .Machine$double.eps)
    if(length(zero.idx)) {
        obs.prop[zero.idx] <- as.numeric(NA)
    }
    # the usual LR formula
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
                                      LR.min = 3.0,
                                      GF.min = 3.0,
                                      p.value = FALSE) {

    # default statistics
    statistic <- toupper(statistic)
    if(!is.null(lavobject)) {
        if(length(statistic) == 1L && statistic == "DEFAULT") {
            statistic <- c("LR", "LR.AVERAGE")
        } else {
            stopifnot(statistic %in% c("GF","LR","GF.BVN","LR.BVN",
                                       "RMSEA.BVN", "RMSEA",
                                       "LR.AVERAGE",
                                       "LR.NLARGE",
                                       "LR.PLARGE",
                                       "GF.AVERAGE",
                                       "GF.NLARGE",
                                       "GF.PLARGE"))
        }
    } else {
        if(length(statistic) == 1L && statistic == "DEFAULT") {
            # if data, none by default
            statistic <- character(0L)
        } else {
            stopifnot(statistic %in% c("GF.BVN","LR.BVN","RMSEA.BVN"))
        }
    }

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
    if(any(c("LR","LR.AVERAGE","LR.PLARGE","LR.NLARGE") %in% statistic)) {
        stat.cell <- c(stat.cell, "LR")
    }
    if(any(c("GF","GF.AVERAGE","GF.PLARGE","GF.NLARGE") %in% statistic)) {
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
    if(is.null(out.cell$group)) {
        out <- out.cell[row.idx,c("lhs","rhs","nobs"),drop=FALSE]
    } else {
        out <- out.cell[row.idx,c("lhs","rhs","group", "nobs"),drop=FALSE]
    }

    # df
    if(length(statistic) > 0L) {
        nrow <- tapply(out.cell$row, INDEX=out.cell$id, FUN=max)
        ncol <- tapply(out.cell$col, INDEX=out.cell$id, FUN=max)
        out$df <- nrow*ncol - nrow - ncol
    }

    # GF
    if("GF" %in% statistic) {
        out$GF <- tapply(out.cell$GF, INDEX=out.cell$id, FUN=sum, 
                         na.rm=TRUE)
        if(p.value) {
            out$GF.pval <- pchisq(out$GF, df=out$df, lower.tail=FALSE)
        }
    }
    if("GF.BVN" %in% statistic) {
        out$GF.BVN <- tapply(out.cell$GF.BVN, INDEX=out.cell$id, FUN=sum,
                             na.rm=TRUE)
        if(p.value) {
            out$GF.BVN.pval <- pchisq(out$GF.BVN, df=out$df, lower.tail=FALSE)
        }
    }
 
    # LR
    if("LR" %in% statistic) {
        out$LR <- tapply(out.cell$LR, INDEX=out.cell$id, FUN=sum,
                         na.rm=TRUE)
        if(p.value) {
            out$LR.pval <- pchisq(out$LR, df=out$df, lower.tail=FALSE)
        }
    }
    if("LR.BVN" %in% statistic) {
        out$LR.BVN <- tapply(out.cell$LR.BVN, INDEX=out.cell$id, FUN=sum,
                             na.rm=TRUE)
        if(p.value) {
            out$LR.BVN.pval <- pchisq(out$LR.BVN, df=out$df, lower.tail=FALSE)
        }
    }

    if("RMSEA" %in% statistic) {
        LR <- tapply(out.cell$LR, INDEX=out.cell$id, FUN=sum, na.rm=TRUE)
        # note: there seems to be a mistake in Appendix 1 eqs 43/44 of Joreskog
        # SSI paper (2005) 'SEM with ordinal variables using LISREL'
        # 2*N*d should N*d
        out$RMSEA <- sqrt( pmax(0, (LR - out$df)/ (out$nobs*out$df) ) )
        if(p.value) {
            # note: MUST use 1 - pchisq (instead of lower.tail = FALSE)
            # because for ncp > 80, routine only computes lower tail
            out$RMSEA.pval <- 1.0 - pchisq(LR,
                                           ncp = 0.1^2*out$nobs*out$df,
                                           df=out$df, lower.tail = TRUE)
        }
    }
    if("RMSEA.BVN" %in% statistic) {
        LR <- tapply(out.cell$LR.BVN, INDEX=out.cell$id, FUN=sum,
                     na.rm=TRUE)
        # note: there seems to be a mistake in Appendix 1 eqs 43/44 of Joreskog
        # SSI paper (2005) 'SEM with ordinal variables using LISREL'
        # 2*N*d should N*d
        out$RMSEA.BVN <- sqrt( pmax(0, (LR - out$df)/ (out$nobs*out$df) ) )
        if(p.value) {
            # note: MUST use 1 - pchisq (instead of lower.tail = FALSE)
            # because for ncp > 80, routine only computes lower tail
            out$RMSEA.BVN.pval <- 1.0 - pchisq(LR,
                                               ncp = 0.1^2*out$nobs*out$df,
                                               df=out$df, lower.tail = TRUE)
        }
    }

    if("LR.AVERAGE" %in% statistic) {
        out$LR.average <- tapply(out.cell$LR, INDEX=out.cell$id, FUN=mean,
                                 na.rm=TRUE)
    }

    if("LR.NLARGE" %in% statistic) {
         out$LR.min <- rep(LR.min, length(out$lhs))
         out$LR.nlarge <- tapply(out.cell$LR, INDEX=out.cell$id,
                                 FUN=function(x) sum(x > LR.min, na.rm=TRUE) )
    }

    if("LR.PLARGE" %in% statistic) {
         out$LR.min <- rep(LR.min, length(out$lhs))
         out$LR.plarge <- tapply(out.cell$LR, INDEX=out.cell$id,
                                 FUN=function(x) sum(x > LR.min, na.rm=TRUE)/length(x) )
    }

    if("GF.AVERAGE" %in% statistic) {
        out$GF.average <- tapply(out.cell$GF, INDEX=out.cell$id, FUN=mean,
                                 na.rm=TRUE)

    }

    if("GF.NLARGE" %in% statistic) {
         out$GF.min <- rep(GF.min, length(out$lhs))
         out$GF.nlarge <- tapply(out.cell$GF, INDEX=out.cell$id, 
                                 FUN=function(x) sum(x > GF.min, na.rm=TRUE) )   
    }

    if("GF.PLARGE" %in% statistic) {
         out$GF.min <- rep(GF.min, length(out$lhs))
         out$GF.plarge <- tapply(out.cell$GF, INDEX=out.cell$id,
                                 FUN=function(x) sum(x > GF.min, na.rm=TRUE)/length(x) )
    }
    
    out
}


lav_tables_oneway <- function(lavobject = NULL, lavdata = NULL, 
                              statistic = NULL) {

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
        TABLE <- lapply(TABLE, as.data.frame, stringsAsFactors=FALSE)
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

    # default statistics
    statistic <- toupper(statistic)
    if(!is.null(lavobject)) {
        if(length(statistic) == 1L && statistic == "DEFAULT") {
            statistic <- c("LR")
        } else {
            stopifnot(statistic %in% c("TH.UVN", 
                                       "TH", "LR", "GF"))
        }

        # sample based
        # note, there is no LR.UVN or GF.UVN: always saturated!
        if("TH.UVN" %in% statistic) {
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

            out$th.uvn <- th
        }
        
        # model based
        if(any(c("TH","LR","GF") %in% statistic)) {

            # model based
            th.h0 <- unlist(lapply(1:lavdata@ngroups, function(x) {
                    unname(unlist(tapply(lavobject@Fit@TH[[x]],
                                  INDEX=lavobject@SampleStats@th.idx[[x]],
                                  function(x) c(x,Inf)))) }))

            est.prop <- unname(unlist(tapply(th.h0, INDEX=out$id, 
                            FUN=function(x) (pnorm(c(x,Inf)) - 
                                pnorm(c(-Inf,x)))[-(length(x)+1)]  )))
            out$est.prop <- est.prop
     
            if("TH" %in% statistic) {
                out$th <- th.h0
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
            stopifnot(statistic %in% c("TH.UVN"))
        }

        if("TH.UVN" %in% statistic) {
            out$th.uvn <- unlist(tapply(out$obs.prop, INDEX=out$id,
                                 FUN=function(x) qnorm(cumsum(x))))
        }
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
            TABLE <- lapply(TABLE, as.data.frame, stringsAsFactors=FALSE)
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
        TH.IDX <- lavobject@SampleStats@th.idx
    } else if(!is.null(lavdata)) {
        COR <- lav_cor(lav.data = lavdata, WLS.W = FALSE, details = TRUE,
                       labels = FALSE, verbose = FALSE)
        # relist
        if(!is.list(COR)) {
            COR <- list(COR)
        }
        TH <- lapply(lapply(COR, attr, "TH"), unlist)
        TH.IDX <- lapply(lapply(COR, attr, "TH.IDX"), unlist)
    } else {
        stop("lavaan ERROR: both lavobject and lavdata are NULL")
    }

    lav_tables_pairwise_sample_pi_cor(COR = COR, TH = TH,
                                      TH.IDX = TH.IDX)
}

# low-level function to compute expected proportions per cell
lav_tables_pairwise_sample_pi_cor <- function(COR = NULL, TH = NULL, 
                                              TH.IDX = NULL) {

    ngroups <- length(COR)

    PI <- vector("list", length=ngroups)
    for(g in 1:ngroups) {
        Sigmahat <- COR[[g]]
        cors <- Sigmahat[lower.tri(Sigmahat)]
        if(any(abs(cors) > 1)) {
            warning("lavaan WARNING: some model-implied correlations are larger than 1.0")
        }
        nvar <- nrow(Sigmahat)
        th.idx <- TH.IDX[[g]]

        # reconstruct ov.types
        ov.types <- rep("numeric", nvar)
        ord.idx <- unique(th.idx[th.idx > 0])
        ov.types[ord.idx] <- "ordered"

        PI.group <- integer(0)
        # order! first i, then j, vec(table)!
        for(i in seq_len(nvar-1L)) {
            for(j in (i+1L):nvar) {
                if(ov.types[i] == "ordered" && ov.types[j] == "ordered") {
                    PI.table <- pc_PI(rho   = Sigmahat[i,j],
                                      th.y1 = TH[[g]][ th.idx == i ],
                                      th.y2 = TH[[g]][ th.idx == j ])
                    PI.group <- c(PI.group, vec(PI.table))
                }
            }
        }
        PI[[g]] <- PI.group
    } # g
    
    PI
}

# low-level function to compute expected proportions per PATTERN
# using sample-based correlations + thresholds
#
# object can be either lavData or lavaan class
#
# only valid if estimator = FML, POM or NOR
#
lav_tables_resp_pi <- function(lavobject = NULL, lavdata = NULL,
                               est = "h0") {

    # shortcuts
    ngroups <- lavdata@ngroups

    # h0 or unrestricted?
    if(est == "h0") {
        Sigma.hat <- lavobject@Fit@Sigma.hat
        TH        <- lavobject@Fit@TH
        TH.IDX    <- lavobject@SampleStats@th.idx
    } else {
        if(is.null(lavobject)) {
            Sigma.hat <- lav_cor(lav.data = lavdata, WLS.W = FALSE, 
                                 details = TRUE, labels = FALSE, 
                                 verbose = FALSE)
            # relist
            if(!is.list(Sigma.hat)) {
                Sigma.hat <- list(Sigma.hat)
            }
            TH <- lapply(lapply(Sigma.hat, attr, "TH"), unlist)
            TH.IDX <- lapply(lapply(Sigma.hat, attr, "TH.IDX"), unlist)
        } else {
            Sigma.hat <- lavobject@SampleStats@cov
            TH        <- lavobject@SampleStats@th
            TH.IDX    <- lavobject@SampleStats@th.idx
        }
    }

    PI <- vector("list", length=ngroups)
    for(g in 1:ngroups) {
        Sigmahat <- Sigma.hat[[g]]
        cors <- Sigmahat[lower.tri(Sigmahat)]
        if(any(abs(cors) > 1)) {
            warning("lavaan WARNING: some model-implied correlations are larger than 1.0")
        }
        nvar <- nrow(Sigmahat)
        th.idx <- TH.IDX[[g]]
        MEAN <- rep(0, nvar)

        # reconstruct ov.types
        ov.types <- rep("numeric", nvar)
        ord.idx <- unique(th.idx[th.idx > 0])
        ov.types[ord.idx] <- "ordered"

        if(all(ov.types == "ordered")) {
            # get patterns ## FIXME GET it
            if(!is.null(lavdata@Rp[[g]]$pat)) {
                PAT <- lavdata@Rp[[g]]$pat
            } else {
                PAT <- lav_data_resppatterns( lavdata@X[[g]] )
            }
            npatterns <- nrow(PAT)
            freq <- as.numeric( rownames(PAT) )
            PI.group <- numeric(npatterns)
            TH.VAR <- lapply(1:nvar, 
                function(x) c(-Inf, TH[[g]][th.idx==x], +Inf))
            # FIXME!!! ok to set diagonal to 1.0?
            diag(Sigmahat) <- 1.0
            for(r in 1:npatterns) {
                # compute probability for each pattern
                lower <- sapply(1:nvar, function(x) TH.VAR[[x]][ PAT[r,x]      ])
                upper <- sapply(1:nvar, function(x) TH.VAR[[x]][ PAT[r,x] + 1L ])
                PI.group[r] <- sadmvn(lower, upper, mean=MEAN, varcov=Sigmahat)
            }
        } else { # case-wise
            stop("not implemented")
        }  

        PI[[g]] <- PI.group
    } # g

    PI
}

lav_tables_table_format <- function(out, lavdata = lavdata,
                                    lavobject = lavobject) {

    # determine column we need
    NAMES <- toupper(names(out))
    stat.idx <- which(NAMES %in% c("LR",       "LR.BVN",
                                   "GF",       "GF.BVN",
                                   "RMSEA", "RMSEA.BVN",
                                   "LR.AVERAGE", "LR.PLARGE", "LR.NLARGE",
                                   "GF.AVERAGE", "GF.PLARGE", "GF.NLARGE"))
    if(length(stat.idx) == 0) {
        if(!is.null(out$obs.freq)) {
            stat.idx <- which(NAMES == "OBS.FREQ")
        } else if(!is.null(out$nobs)) {
            stat.idx <- which(NAMES == "NOBS")
        }
        UNI <- NULL
    } else if(length(stat.idx) > 1) {
        stop("lavaan ERROR: more than one statistic for table output: ", 
              paste(NAMES[stat.idx], collapse=" "))
    } else {
        # univariate version of same statistic
        if(NAMES[stat.idx] == "LR.AVERAGE") {
            UNI <- lavTables1D(lavobject, statistic="LR")
        } else if(NAMES[stat.idx] == "GF.AVERAGE") {
            UNI <- lavTables1D(lavobject, statistic="GF")
        } else {
            UNI <- NULL
        }
    }

    OUT <- vector("list", length=lavdata@ngroups)
    for(g in 1:lavdata@ngroups) {
        if(lavdata@ngroups == 1L) { # no group column
            STAT <- out[[stat.idx]]
        } else {
            STAT <- out[[stat.idx]][ out$group == g ]
        }
        RN <- lavdata@ov.names[[g]]
        OUT[[g]] <- getCov(STAT, diagonal = FALSE, lower = FALSE, names = RN)
        # change diagonal elements: replace by univariate stat
        # if possible
        diag(OUT[[g]]) <- as.numeric(NA)
        if(!is.null(UNI)) {
            if(!is.null(UNI$group)) {
                idx <- which( UNI$group == g )
            } else {
                idx <- 1:length(UNI$lhs)
            }
            if(NAMES[stat.idx] == "LR.AVERAGE") { 
                diag(OUT[[g]]) <- tapply(UNI$LR[idx], INDEX=UNI$id[idx], 
                                         FUN=mean)
            } else if(NAMES[stat.idx] == "GF.AVERAGE") {
                diag(OUT[[g]]) <- tapply(UNI$GF[idx], INDEX=UNI$id[idx], 
                                         FUN=mean)
            }
        }
        class(OUT[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
    }
    if(lavdata@ngroups > 1L) {
        names(OUT) <- lavdata@group.label
        out <- OUT
    } else {
        out <- OUT[[1]]      
    }

    out
}

lav_tables_cells_format <- function(out, lavdata = lavdata) {

    OUT <- vector("list", length=lavdata@ngroups)
    if(is.null(out$group)) {
        out$group <- rep(1L, length(out$lhs))
    }
    # do we have a statistic?
    # determine column we need
    NAMES <- toupper(names(out))
    stat.idx <- which(NAMES %in% c("LR",       "LR.BVN",
                           "GF",       "GF.BVN",
                           "RMSEA", "RMSEA.BVN",
                           "LR.AVERAGE", "LR.PLARGE", "LR.NLARGE",
                           "GF.AVERAGE", "GF.PLARGE", "GF.NLARGE"))
    if(length(stat.idx) == 0) {
        statistic <- "obs.freq"
    } else if(length(stat.idx) > 1) {
         stop("lavaan ERROR: more than one statistic for table output: ",
              paste(NAMES[stat.idx], collapse=" "))
    } else {
        statistic <- NAMES[stat.idx]
    }

    for(g in 1:lavdata@ngroups) {
        case.idx <- which( out$group == g )
        ID.group <- unique( out$id[ out$group == g] )
        TMP <-lapply(ID.group, function(x) {
                  Tx <- out[out$id == x,]
                  M <- matrix(Tx[,statistic],
                              max(Tx$row), max(Tx$col))
                  rownames(M) <- unique(Tx$row)
                  colnames(M) <- unique(Tx$col)
                  class(M) <- c("lavaan.matrix", "matrix")
                  M })
        names(TMP) <- unique(paste(out$lhs[case.idx], out$rhs[case.idx], 
                                   sep="_"))
        OUT[[g]] <- TMP
    }
    if(lavdata@ngroups > 1L) {
        out <- OUT
        names(out) <- lavdata@group.label
    } else {
        out <- OUT[[1L]]
    }

    out
}
