# constructor for the lavParTable model description
#
# initial version: YR 22/05/2009
#  major revision: YR 02/11/2010: - FLATTEN the model syntax and turn it into a
#                                   data.frame, with a "modifiers" attribute
#                                 - add default elements here
#                                 - check for duplicate elements
#                                 - allow for every possible model...
#                                 - since 0.4-5
#                                 - the end result is a full description of
#                                   a model (but no matrix representation)
# - 14 Jan 2014: merge 02lavaanUser.R with lav_partable.R
#                move syntax-based code to lav_syntax.R

lavaanify <- lavParTable <- function(

                      model            = NULL, 
                      meanstructure    = FALSE,
                      int.ov.free      = FALSE,
                      int.lv.free      = FALSE,
                      orthogonal       = FALSE, 
                      std.lv           = FALSE,
                      fixed.x          = TRUE,
                      parameterization = "delta",
                      constraints      = NULL,

                      auto             = FALSE,
                      model.type       = "sem",
                      auto.fix.first   = FALSE,
                      auto.fix.single  = FALSE,
                      auto.var         = FALSE,
                      auto.cov.lv.x    = FALSE,
                      auto.cov.y       = FALSE,
                      auto.th          = FALSE,
                      auto.delta       = FALSE,

                      varTable         = NULL,
                      ngroups          = 1L,
                      group.equal      = NULL,
                      group.partial    = NULL,
                      group.w.free     = FALSE,
                      debug            = FALSE,
                      warn             = TRUE,
                      
                      as.data.frame.   = TRUE) {


    # check if model is already FLAT or a full parameter table
    if(is.list(model) && !is.null(model$lhs)) {
        if(is.null(model$mod.idx)) {
            warning("lavaan WARNING: input already looks like a parameter table; returning as is")
            return(model)
        } else {
            FLAT <- model
        }
    } else {
        # parse the model syntax and flatten the user-specified model
        # return a data.frame, where each line is a model element (rhs, op, lhs)
        FLAT <- lavParseModelString(model.syntax=model, warn=warn, debug=FALSE)
    }
    # user-specified *modifiers* are returned as an attribute
    MOD  <- attr(FLAT, "modifiers"); attr(FLAT, "modifiers") <- NULL
    # user-specified *constraints* are returned as an attribute
    CON  <- attr(FLAT, "constraints"); attr(FLAT, "constraints") <- NULL

    # extra constraints?
    if(!is.null(constraints) && nchar(constraints) > 0L) {
        FLAT2 <- lavParseModelString(model.syntax=constraints, warn=warn)
        CON2 <- attr(FLAT2, "constraints"); rm(FLAT2)
        CON <- c(CON, CON2)
    }
   
    if(debug) {
        cat("[lavaan DEBUG]: FLAT (flattened user model):\n")
        print(FLAT)
        cat("[lavaan DEBUG]: MOD (modifiers):\n")
        print( str(MOD) )
        cat("[lavaan DEBUG]: CON (constraints):\n")
        print( str(CON) )
    }

    # check for wrongly specified variances/covariances/intercepts
    # of exogenous variables in model syntax (if fixed.x=TRUE)
    if(fixed.x) { # we ignore the groups here!
        # we only call this function for the warning message
        tmp <- vnames(FLAT, "ov.x", warn=TRUE); rm(tmp)
    }

    # auto=TRUE?
    if(auto && model.type == "sem") { # mimic sem/cfa auto behavior
        if(model.type == "sem") {
            int.ov.free     = TRUE
            int.lv.free     = FALSE
            auto.fix.first  = !std.lv
            auto.fix.single = TRUE
            auto.var        = TRUE
            auto.cov.lv.x   = TRUE
            auto.cov.y      = TRUE
            auto.th         = TRUE
            auto.delta      = TRUE
        } else 
        
        if(model.type == "growth") {
            model.type      = "growth"
            int.ov.free     = FALSE
            int.lv.free     = TRUE
            auto.fix.first  = !std.lv
            auto.fix.single = TRUE
            auto.var        = TRUE
            auto.cov.lv.x   = TRUE
            auto.cov.y      = TRUE
            auto.th         = TRUE    
            auto.delta      = TRUE
        }
    }

    # check for meanstructure
    if(any(FLAT$op == "~1")) meanstructure <- TRUE

    # check for group identifiers in the syntax (op = ":")
    n.group.flat <- length(which(FLAT$op == ":"))
    if(n.group.flat > 0L && n.group.flat != ngroups) {
        stop("lavaan ERROR: syntax defines ", n.group.flat, " groups; ",
             "data suggests ", ngroups, " groups")
    }

    # for each `group' in n.group.flat, produce a USER
    if(n.group.flat > 0L) {
        # split the FLAT data.frame per `group', create LIST
        # for each `group', and bind them together
        FLAT <- as.data.frame(FLAT, stringsAsFactors=FALSE)
        group.op.idx <- c(which(FLAT$op == ":"), nrow(FLAT)+1L)
        for(g in 1:n.group.flat) {
            FLAT.group <- FLAT[(group.op.idx[g]+1L):(group.op.idx[g+1]-1L),]    
            LIST.group <- lav_partable_flat(FLAT.group, meanstructure = meanstructure, 
                int.ov.free = int.ov.free, int.lv.free = int.lv.free,
                orthogonal = orthogonal, std.lv = std.lv, fixed.x = fixed.x,
                parameterization = parameterization,
                auto.fix.first = auto.fix.first, 
                auto.fix.single = auto.fix.single,
                auto.var = auto.var, auto.cov.lv.x = auto.cov.lv.x,
                auto.cov.y = auto.cov.y, auto.th = auto.th, 
                auto.delta = auto.delta, 
                varTable = varTable, group.equal = NULL, 
                group.w.free = group.w.free, ngroups = 1L)
            LIST.group <- as.data.frame(LIST.group, stringsAsFactors=FALSE)
            if(g == 1L) {
                LIST <- LIST.group
            } else {
                LIST.group$group <- rep(g, length(LIST.group$lhs))
                LIST.group$id <- LIST.group$id + max(LIST$id)
                LIST <- rbind(LIST, LIST.group)
            }
        }
        LIST <- as.list(LIST)
    } else {
        LIST <- lav_partable_flat(FLAT, meanstructure = meanstructure, 
            int.ov.free = int.ov.free, int.lv.free = int.lv.free,
            orthogonal = orthogonal, std.lv = std.lv, fixed.x = fixed.x,
            parameterization = parameterization,
            auto.fix.first = auto.fix.first, auto.fix.single = auto.fix.single,
            auto.var = auto.var, auto.cov.lv.x = auto.cov.lv.x,
            auto.cov.y = auto.cov.y, auto.th = auto.th, 
            auto.delta = auto.delta,
            varTable = varTable, group.equal = group.equal, 
            group.w.free = group.w.free,
            ngroups = ngroups)
    }        
    if(debug) {
        cat("[lavaan DEBUG]: parameter LIST without MODIFIERS:\n")
        print( as.data.frame(LIST, stringsAsFactors=FALSE) )
    }

    # apply user-specified modifiers
    if(length(MOD)) {
        for(el in 1:length(MOD)) {
            idx <- which(LIST$mod.idx == el) # for each group
            MOD.fixed <- MOD[[el]]$fixed  
            MOD.start <- MOD[[el]]$start
            MOD.label <- MOD[[el]]$label 
            MOD.prior <- MOD[[el]]$prior

            # check for single argument if multiple groups
            if(ngroups > 1L && length(idx) > 1L) {
                # Ok, this is not very consistent:
                # A) here we force same behavior across groups
                if(length(MOD.fixed) == 1L) MOD.fixed <- rep(MOD.fixed, ngroups)
                if(length(MOD.start) == 1L) MOD.start <- rep(MOD.start, ngroups)
                if(length(MOD.prior) == 1L) MOD.prior <- rep(MOD.prior, ngroups)
                # B) here we do NOT! otherwise, it would imply an equality 
                #                    constraint...
                #    except if group.equal="loadings"!
                if(length(MOD.label) == 1L) {
                    if("loadings" %in% group.equal) {
                        MOD.label <- rep(MOD.label, ngroups)
                    } else {
                        MOD.label <- c(MOD.label, rep("", (ngroups-1L)) )
                    }
                }
            }

            # check for wrong number of arguments if multiple groups
            nidx <- length(idx)
            if( (!is.null(MOD.fixed) && nidx != length(MOD.fixed)) ||
                (!is.null(MOD.start) && nidx != length(MOD.start)) ||
                (!is.null(MOD.prior) && nidx != length(MOD.prior)) ||
                (!is.null(MOD.label) && nidx != length(MOD.label)) ) {
                el.idx <- which(LIST$mod.idx == el)[1L]
                stop("lavaan ERROR: wrong number of arguments in modifier (",
                    paste(MOD.label, collapse=","), ") of element ", 
                    LIST$lhs[el.idx], LIST$op[el.idx], LIST$rhs[el.idx])
            }

            # apply modifiers
            if(!is.null(MOD.fixed)) {
                # two options: constant or NA
                na.idx <- which(is.na(MOD.fixed))
                not.na.idx <- which(!is.na(MOD.fixed))

                # constant
                LIST$ustart[idx][not.na.idx] <- MOD.fixed[not.na.idx]
                LIST$free[  idx][not.na.idx] <- 0L

                # NA* modifier
                LIST$free[  idx][na.idx] <- 1L # eg factor loading
                LIST$ustart[idx][na.idx] <- as.numeric(NA)
            }
            if(!is.null(MOD.start)) {
                LIST$ustart[idx] <- MOD.start
            }
            if(!is.null(MOD.prior)) {
                # do we already have a `prior' column? if not, create one
                if(is.null(LIST$prior)) {
                    LIST$prior <- character( length(LIST$lhs) )
                }
                LIST$prior[idx] <- MOD.prior
            }
            if(!is.null(MOD.label)) {
                LIST$label[idx] <- MOD.label
            }
        }
    }
    # remove mod.idx column
    LIST$mod.idx <- NULL

    if(debug) {
        cat("[lavaan DEBUG]: parameter LIST with MODIFIERS:\n")
        print( as.data.frame(LIST, stringsAsFactors=FALSE) )
    }

    # get 'virtual' parameter labels
    LABEL <- lav_partable_labels(partable=LIST, group.equal=group.equal,
                                group.partial=group.partial)

    if(debug) {
        cat("[lavaan DEBUG]: parameter LIST with LABELS:\n")
        tmp <- LIST; tmp$LABEL <- LABEL
        print( as.data.frame(tmp, stringsAsFactors=FALSE) )
    }

    # handle user-specified equality constraints
    # lavaan 0.5-18
    # - rewrite 'LABEL-based' equality constraints as == constraints
    # - create plabel: internal labels, based on id
    # - create CON entries, using these internal labels
    LIST$plabel <- paste(".p", LIST$id, ".", sep="")
    idx.eq.label <- which(duplicated(LABEL))
    if(length(idx.eq.label) > 0L) {
        CON.idx <- length(CON)
        # add 'user' column
        CON <- lapply(CON, function(x) {x$user <- 1L; x} )
        for(idx in idx.eq.label) {
            eq.label <- LABEL[idx]
            all.idx <- which(LABEL == eq.label) # all same-label parameters   
            ref.idx <- all.idx[1L]              # the first one only 

            # two possibilities: 
            # 1. all.idx contains a fixed parameter: in this case,
            #    we fix them all (hopefully to the same value)
            # 2. all.idx contains only free parameters

            # 1. fixed?
            if(any(LIST$free[all.idx] == 0L)) {

                # which one is fixed? (only pick the first!)
                fixed.all <- all.idx[ LIST$free[all.idx] == 0L ]
                fixed.idx <- fixed.all[1] # only pick the first!

                # sanity check: are all ustart values equal?
                ustart1 <- LIST$ustart[ fixed.idx ]
                if(! all(ustart1 == LIST$ustart[fixed.all]) ) {
                    warning("lavaan WARNING: equality constraints involve fixed parameters with different values; only the first one will be used")
                }


                fixed.idx <- fixed.all[1] # only pick the first!

                # fix current 'idx'
                LIST$ustart[idx] <- LIST$ustart[fixed.idx]
                LIST$free[idx] <- 0L  # not free anymore, since it must
                                      # be equal to the 'fixed' parameter
                                      # (Note: Mplus ignores this)

                # just in case: if ref.idx is not equal to fixed.idx, 
                # fix this one too
                LIST$ustart[ref.idx] <- LIST$ustart[fixed.idx]
                LIST$free[ref.idx] <- 0L
            } else {
            # 2. ref.idx is a free parameter
                # user-label?
                #if(nchar(LIST$label[ref.idx])  > 0) {
                #    lhs.lab <- LIST$label[ref.idx]
                #} else {
                #    lhs.lab <- PLABEL[ref.idx]
                #}
                CON.idx <- CON.idx + 1L
                CON[[CON.idx]] <- list(op   = "==", 
                                       lhs  = LIST$plabel[ref.idx], 
                                       rhs  = LIST$plabel[idx],
                                       user = 2L)

                # just to trick semTools, also add something in the label
                # colum, *if* it is empty
                for(i in all.idx) {
                    if(nchar(LIST$label[i]) == 0L) {
                        LIST$label[i] <- LIST$plabel[ ref.idx ]
                    }
                }
            }
        }
    }
    if(debug) {
        print(CON)
    }


    # count free parameters
    idx.free <- which(LIST$free > 0)
    LIST$free[idx.free] <- seq_along(idx.free)
    if(!is.null(LIST$unco)) {
         LIST$unco[idx.free] <- seq_along(idx.free)
    }

    # 2. add free counter to this element
    #idx.equal <- which(LIST$eq.id > 0)
    #LIST$free[idx.equal] <- LIST$free[ LIST$eq.id[idx.equal] ]

    # 3. which parameters would be free without equality constraints?
    #idx.unco <- which(LIST$free > 0)
    #LIST$unco[idx.unco] <- seq_along(idx.unco)


    # handle constraints (if any) (NOT per group, but overall - 0.4-11)
    if(length(CON) > 0L) {
        #cat("DEBUG:\n"); print(CON)
        lhs  = unlist(lapply(CON, "[[", "lhs"))
         op  = unlist(lapply(CON, "[[",  "op"))
        rhs  = unlist(lapply(CON, "[[", "rhs"))
        user = unlist(lapply(CON, "[[", "user"))
        LIST$id         <- c(LIST$id,         length(LIST$id) + seq_along(lhs) )
        LIST$lhs        <- c(LIST$lhs,        lhs)
        LIST$op         <- c(LIST$op,         op)
        LIST$rhs        <- c(LIST$rhs,        rhs)
        LIST$user       <- c(LIST$user,       user)
        LIST$group      <- c(LIST$group,      rep(0L, each=length(CON)))
        LIST$free       <- c(LIST$free,       rep(0L, length(lhs)) )
        LIST$ustart     <- c(LIST$ustart,     rep(as.numeric(NA), length(lhs)))
        LIST$exo        <- c(LIST$exo,        rep(0L, length(lhs)) )
        LIST$label      <- c(LIST$label,      rep("",  length(lhs)) )
        if(!is.null(LIST$prior)) {
            LIST$prior      <- c(LIST$prior,      rep("",  length(lhs)) )
        }
        LIST$plabel     <- c(LIST$plabel,     rep("",  length(lhs)) )
        if(!is.null(LIST$eq.id)) {
            LIST$eq.id      <- c(LIST$eq.id,      rep(0L,  length(lhs)) )
        }
        if(!is.null(LIST$unco)) {
            LIST$unco       <- c(LIST$unco,       rep(0L,  length(lhs)) )
        }
    }

    # put lhs of := elements in label column
    def.idx <- which(LIST$op == ":=")
    LIST$label[def.idx] <- LIST$lhs[def.idx]


    if(debug) { 
        cat("[lavaan DEBUG] lavParTable\n")
        print( as.data.frame(LIST) )
    }

    # data.frame?
    if(as.data.frame.) {
        LIST <- as.data.frame(LIST, stringsAsFactors = FALSE)

        # order? first by 'op', then by 'user'
        #LIST[with(LIST, order(op, -user)),]
    }

    LIST
}


# lav_partable (old name: utils-user.R)
#
# functions to generate/compute/extract information from the lavaan
# `parameter table'
#
# YR. 29 june 2013 (as lav_partable)

# user visible function to add 'matrix' entries in the parameter table
lavMatrixRepresentation <- function(partable, representation = "LISREL",
                                    as.data.frame. = TRUE) {

    # get model matrices
    if(representation == "LISREL") {
        REP <- representation.LISREL(partable, target=NULL, extra=FALSE)
    } else {
        stop("lavaan ERROR: only representation \"LISREL\" has been implemented.")
    }

    partable$mat <- REP$mat
    partable$row <- REP$row
    partable$col <- REP$col

    if(as.data.frame.) {
        partable <- as.data.frame(partable, stringsAsFactors=FALSE)
        class(partable) <- c("lavaan.data.frame", "data.frame")
    }

    partable
}


# return 'attributes' of a lavaan partable -- generate a new set if necessary
lav_partable_attributes <- function(partable, pta=NULL) {

    if(is.null(pta)) {
        # attached to partable?
        pta <- attributes(partable)
        if(!is.null(pta$vnames) && !is.null(pta$ngroups)) {
            # looks like a pta
            return(pta)
        } else {
            pta <- list()
        }
    }

    # vnames
    pta$vnames <- lav_partable_vnames(partable, type="all", group="list")

    # ngroups
    pta$ngroups <- length(pta$vnames$ov)

    # nvar
    pta$nvar <- lapply(pta$vnames$ov, length)

    # nfac
    pta$nfac <- lapply(pta$vnames$lv, length)

    # nfac.nonnormal - for numerical integration
    pta$nfac.nonnormal <- lapply(pta$vnames$lv.nonnormal, length)

    pta
}


lav_partable_ov_idx <- function(partable, type="th", group=NULL) {

    stopifnot(is.list(partable), !missing(type),
              type %in% c("th"))

    if(type == "th") {
        ovn <- lav_partable_vnames(partable, type="ov.nox", group=group)
        ov.num <- lav_partable_vnames(partable, type="ov.num", group=group)
        th  <- lav_partable_vnames(partable, type="th.mean", group=group)
        th[th %in% ov.num] <- "__NUM__"
        th1 <- gsub("\\|t[0-9]*","",th)
        out <- match(th1, ovn)
        out[is.na(out)] <- 0
    }

    out
}

lav_partable_ndat <- function(partable, group=NULL) {

    ngroups <- max(partable$group)
    meanstructure <- any(partable$op == "~1")
    fixed.x <- any(partable$exo > 0L & partable$free == 0L)
    categorical <- any(partable$op == "|")
    if(categorical) meanstructure <- TRUE

    if(categorical) {
        ov.names <- lav_partable_vnames(partable, "ov.nox", group=group)
    } else {
        ov.names <- lav_partable_vnames(partable, "ov", group=group)
    }
    nvar <- length(ov.names)

    pstar <- nvar*(nvar+1)/2; if(meanstructure) pstar <- pstar + nvar
    ndat  <- ngroups*pstar

    # correction for fixed.x?
    if(!categorical && fixed.x) {
        ov.names.x <- lav_partable_vnames(partable, "ov.x", group=group)
        nvar.x <- length(ov.names.x)
        pstar.x <- nvar.x * (nvar.x + 1) / 2
        if(meanstructure) pstar.x <- pstar.x + nvar.x
        ndat <- ndat - (ngroups * pstar.x)
    }

    # correction for ordinal data?
    if(categorical) {
        ov.names.x <- lav_partable_vnames(partable, "ov.x", group=group)
        nexo     <- length(ov.names.x)
        ov.ord   <- lav_partable_vnames(partable, "ov.ord", group=group)
        nvar.ord <- length(ov.ord)
        th       <- lav_partable_vnames(partable, "th", group=group)
        nth      <- length(th)
        # no variances
        ndat <- ndat - (ngroups * nvar.ord)
        # no means
        ndat <- ndat - (ngroups * nvar.ord)
        # but additional thresholds
        ndat <- ndat + (ngroups * nth)
        # add slopes
        ndat <- ndat + (ngroups * nvar * nexo)
    }

    # correction for group proportions?
    group.idx <- which(partable$lhs == "group" & 
                       partable$op == "%")
    if(length(group.idx) > 0L) {
        # ndat <- ndat + (length(group.idx) - 1L) # G - 1 (sum to one)
        ndat <- ndat + length(group.idx) # poisson: each cell a parameter
    }

    ndat
}

lav_partable_npar <- function(partable) {
    npar <- max(partable$free)
    npar
}

lav_partable_df <- function(partable, group=NULL) {

    npar <- lav_partable_npar(partable)
    ndat <- lav_partable_ndat(partable, group=group)

    # degrees of freedom
    df <- ndat - npar

    as.integer(df)
}

lav_partable_labels <- function(partable, group.equal="", group.partial="", 
                               type="user") {

    # catch empty partable
    if(length(partable$lhs) == 0L) return(character(0L))

    # default labels
    label <- paste(partable$lhs, partable$op, partable$rhs, sep="")
    
    # handle multiple groups
    ngroups <- max(partable$group)
    if(ngroups > 1L) {
        for(g in 2:ngroups) {
            label[partable$group == g] <- 
                paste(label[partable$group == g], ".g", g, sep="")
        }
    }
 
    #cat("DEBUG: label start:\n"); print(label); cat("\n")
    #cat("group.equal = ", group.equal, "\n")
    #cat("group.partial = ", group.partial, "\n")

    # use group.equal so that equal sets of parameters get the same label
    if(ngroups > 1L && length(group.equal) > 0L) {

        if("intercepts" %in% group.equal ||
           "residuals"  %in%  group.equal ||
           "residual.covariances" %in%  group.equal) {
            ov.names.nox <- vector("list", length=ngroups)
            for(g in 1:ngroups)
                ov.names.nox[[g]] <- lav_partable_vnames(partable, "ov.nox", group=g)
        }
        if("thresholds" %in% group.equal) {
            ov.names.ord <- vector("list", length=ngroups)
            for(g in 1:ngroups)
                ov.names.ord[[g]] <- lav_partable_vnames(partable, "ov.ord", group=g)
        }
        if("means" %in% group.equal ||
           "lv.variances" %in% group.equal ||
           "lv.covariances" %in% group.equal) {
            lv.names <- vector("list", length=ngroups)
            for(g in 1:ngroups)
                lv.names[[g]] <- lav_partable_vnames(partable, "lv", group=g)
        }

        # g1.flag: TRUE if included, FALSE if not
        g1.flag <- logical(length(which(partable$group == 1L)))

        # LOADINGS
        if("loadings" %in% group.equal)
            g1.flag[ partable$op == "=~" & partable$group == 1L  ] <- TRUE
        # INTERCEPTS (OV)
        if("intercepts" %in% group.equal)
            g1.flag[ partable$op == "~1"  & partable$group == 1L  &
                     partable$lhs %in% ov.names.nox[[1L]] ] <- TRUE
        # THRESHOLDS (OV-ORD)
        if("thresholds" %in% group.equal)
            g1.flag[ partable$op == "|"  & partable$group == 1L  &
                     partable$lhs %in% ov.names.ord[[1L]] ] <- TRUE
        # MEANS (LV)
        if("means" %in% group.equal)
            g1.flag[ partable$op == "~1" & partable$group == 1L &
                     partable$lhs %in% lv.names[[1L]] ] <- TRUE
        # REGRESSIONS
        if("regressions" %in% group.equal)
            g1.flag[ partable$op == "~" & partable$group == 1L ] <- TRUE
        # RESIDUAL variances (FIXME: OV ONLY!)
        if("residuals" %in% group.equal)
            g1.flag[ partable$op == "~~" & partable$group == 1L &
                     partable$lhs %in% ov.names.nox[[1L]] &
                     partable$lhs == partable$rhs ] <- TRUE
        # RESIDUAL covariances (FIXME: OV ONLY!)
        if("residual.covariances" %in% group.equal)
            g1.flag[ partable$op == "~~" & partable$group == 1L &
                     partable$lhs %in% ov.names.nox[[1L]] &
                     partable$lhs != partable$rhs ] <- TRUE
        # LV VARIANCES
        if("lv.variances" %in% group.equal)
            g1.flag[ partable$op == "~~" & partable$group == 1L &
                     partable$lhs %in% lv.names[[1L]] &
                     partable$lhs == partable$rhs ] <- TRUE
        # LV COVARIANCES
        if("lv.covariances" %in% group.equal)
            g1.flag[ partable$op == "~~" & partable$group == 1L &
                     partable$lhs %in% lv.names[[1L]] &
                     partable$lhs != partable$rhs ] <- TRUE

        # if group.partial, set corresponding flag to FALSE
        if(length(group.partial) > 0L) {
            g1.flag[ label %in% group.partial &
                     partable$group == 1L ] <- FALSE
        }

        # for each (constrained) parameter in 'group 1', find a similar one
        # in the other groups (we assume here that the models need
        # NOT be the same across groups!
        g1.idx <- which(g1.flag)
        for(i in 1:length(g1.idx)) {
            ref.idx <- g1.idx[i]
            idx <- which(partable$lhs == partable$lhs[ref.idx] &
                         partable$op  == partable$op[ ref.idx] &
                         partable$rhs == partable$rhs[ref.idx] &
                         partable$group > 1L)
            label[idx] <- label[ref.idx]
        }
    }

    #cat("DEBUG: g1.idx = ", g1.idx, "\n")
    #cat("DEBUG: label after group.equal:\n"); print(label); cat("\n")

    # user-specified labels -- override everything!!
    user.idx <- which(nchar(partable$label) > 0L)
    label[user.idx] <- partable$label[user.idx]

    #cat("DEBUG: user.idx = ", user.idx, "\n")
    #cat("DEBUG: label after user.idx:\n"); print(label); cat("\n")

    # which labels do we need?
    if(type == "user") {
        idx <- 1:length(label)
    } else if(type == "free") {
        idx <- which(partable$free > 0L & !duplicated(partable$free))
    #} else if(type == "unco") {
    #    idx <- which(partable$unco > 0L & !duplicated(partable$unco))
    } else {
        stop("argument `type' must be one of free or user")
    }

    label[idx]
}

# only for simsem ....
getParameterLabels <- lav_partable_labels


lav_partable_full <- function(partable = NULL, group = NULL,
                              strict.exo = FALSE,
                              free = FALSE, start = FALSE) {

    # check minimum requirements: lhs, op, rhs
    stopifnot( !is.null(partable$lhs), 
               !is.null(partable$op), 
               !is.null(partable$rhs) )

    # meanstructure
    meanstructure <- any(partable$op == "~1")

    # number of groups
    if(!is.null(partable$group)) {
        ngroups <- max(partable$group)
    } else {
        ngroups <- 1L
    }

    # extract `names' of various types of variables:
    lv.names     <- lav_partable_vnames(partable, type="lv",  group=group)   # latent variables
    ov.names     <- lav_partable_vnames(partable, type="ov",  group=group)   # observed variables
    ov.names.x   <- lav_partable_vnames(partable, type="ov.x",group=group)   # exogenous x covariates
    ov.names.nox <- lav_partable_vnames(partable, type="ov.nox",group=group) # ov's without exo's
    lv.names.x   <- lav_partable_vnames(partable, type="lv.x",group=group)   # exogenous lv
    ov.names.y   <- lav_partable_vnames(partable, type="ov.y",group=group)   # dependent ov
    lv.names.y   <- lav_partable_vnames(partable, type="lv.y",group=group)   # dependent lv
    lvov.names.y <- c(ov.names.y, lv.names.y)
    ov.names.ord <- lav_partable_vnames(partable, type="ov.ord", group=group)


    # 1 "=~"
    l.lhs <- r.rhs <- op <- character(0)
    l.lhs <- rep(lv.names, each=length(ov.names.nox))
    l.rhs <- rep(ov.names.nox, times=length(lv.names))
    l.op  <- rep("=~", length(l.lhs))

    # 2a. "~~" ov ## FIXME: ov.names.nox or ov.names??
    ov.lhs <- ov.rhs <- ov.op <- character(0)
    if(strict.exo) {
        OV <- ov.names.nox
    } else {
        OV <- ov.names
    }
    nx <- length(OV)
    idx <- lower.tri(matrix(0, nx, nx), diag=TRUE)
    ov.lhs <- rep(OV,  each=nx)[idx] # fill upper.tri
    ov.rhs <- rep(OV, times=nx)[idx]
    ov.op  <- rep("~~", length(ov.lhs))

    # 2b. "~~" lv
    lv.lhs <- lv.rhs <- lv.op <- character(0)
    nx <- length(lv.names)
    idx <- lower.tri(matrix(0, nx, nx), diag=TRUE)
    lv.lhs <- rep(lv.names,  each=nx)[idx] # fill upper.tri
    lv.rhs <- rep(lv.names, times=nx)[idx]
    lv.op  <- rep("~~", length(lv.lhs))

    # 3 regressions?
    r.lhs <- r.rhs <- r.op <- character(0)
    if(any(partable$op == "~")) {

        if(strict.exo) {
            eqs.y <- unique( partable$lhs[ partable$op == "~" &
                                           !partable$lhs %in% ov.names.x ] )
        } else {
            eqs.y <- unique( partable$lhs[ partable$op == "~" ] )
        }
        eqs.x <- unique( partable$rhs[ partable$op == "~" ] )

        r.lhs <- rep(eqs.y, each=length(eqs.x))
        r.rhs <- rep(eqs.x, times=length(eqs.y))

        # remove self-arrows
        idx <- which(r.lhs == r.rhs)
        r.lhs <- r.lhs[-idx]
        r.rhs <- r.rhs[-idx]
        r.op <- rep("~", length(r.rhs))
    }

    # 4. intercepts
    int.lhs <- int.rhs <- int.op <- character(0)
    if(meanstructure) {
        if(strict.exo) {
            int.lhs <- c(ov.names.nox, lv.names)
        } else {
            int.lhs <- c(ov.names, lv.names)
        }
        int.rhs <- rep("",   length(int.lhs))
        int.op  <- rep("~1", length(int.lhs))
    }

    # 5. thresholds
    th.lhs <- th.rhs <- th.op <- character(0)
    if(length(ov.names.ord) > 0L) {
        tmp <- strsplit(lav_partable_vnames(partable, "th", group=group), "\\|")
        th.lhs <- sapply(tmp, function(x) x[1])
        th.rhs <- sapply(tmp, function(x) x[2])
        th.op  <- rep("|", length(th.lhs))
    }

    # 6. scaling parameters
    delta.lhs <- delta.rhs <- delta.op <- character(0)
    if(ngroups > 1L && length(ov.names.ord) > 0L) {
        delta.lhs <- ov.names.ord
        delta.rhs <- ov.names.ord
        delta.op  <- rep("~*~", length(delta.lhs))
    }

    # combine
    lhs <- c(l.lhs, ov.lhs, lv.lhs, r.lhs, int.lhs, th.lhs, delta.lhs)
    rhs <- c(l.rhs, ov.rhs, lv.rhs, r.rhs, int.rhs, th.rhs, delta.rhs)
     op <- c(l.op,  ov.op,  lv.op,  r.op,  int.op,  th.op,  delta.op)


    # multiple groups!
    group <- 1L
    if(ngroups > 1) {
        group   <- rep(1:ngroups, each=length(lhs))
        lhs     <- rep(lhs,     times=ngroups)
        op      <- rep(op,      times=ngroups)
        rhs     <- rep(rhs,     times=ngroups)
    }

    LIST <- data.frame(lhs=lhs, op=op, rhs=rhs, group=group,
                       stringsAsFactors=FALSE)

    if(free) {
        LIST$free <- rep(0L, nrow(LIST))
    }

    if(start) {
        LIST$start <- rep(0, nrow(LIST))
    }

    LIST
}

lav_partable_flat <- function(FLAT = NULL,
                              meanstructure    = FALSE,
                              int.ov.free      = FALSE,
                              int.lv.free      = FALSE,
                              orthogonal       = FALSE,
                              std.lv           = FALSE,
                              fixed.x          = TRUE,
                              parameterization = "delta",
                              auto.fix.first   = FALSE,
                              auto.fix.single  = FALSE,
                              auto.var         = FALSE,
                              auto.cov.lv.x    = FALSE,
                              auto.cov.y       = FALSE,
                              auto.th          = FALSE,
                              auto.delta       = FALSE,
                              varTable         = NULL,
                              group.equal      = NULL,
                              group.w.free     = FALSE,
                              ngroups          = 1L) {

    categorical <- FALSE

    ### DEFAULT elements: parameters that are typically not specified by
    ###                   users, but should typically be considered, 
    ###                   either free or fixed

    # extract `names' of various types of variables:
    lv.names     <- lav_partable_vnames(FLAT, type="lv")     # latent variables
    #lv.names.r   <- lav_partable_vnames(FLAT, type="lv.regular") # regular latent variables
    lv.names.f   <- lav_partable_vnames(FLAT, type="lv.formative") # formative latent variables
    ov.names     <- lav_partable_vnames(FLAT, type="ov")     # observed variables
    ov.names.x   <- lav_partable_vnames(FLAT, type="ov.x")   # exogenous x covariates 
    ov.names.nox <- lav_partable_vnames(FLAT, type="ov.nox")
    lv.names.x   <- lav_partable_vnames(FLAT, type="lv.x")   # exogenous lv
    ov.names.y   <- lav_partable_vnames(FLAT, type="ov.y")   # dependent ov
    lv.names.y   <- lav_partable_vnames(FLAT, type="lv.y")   # dependent lv
    #lvov.names.y <- c(ov.names.y, lv.names.y)
    lvov.names.y <- c(lv.names.y, ov.names.y)


    # get 'ordered' variables, either from FLAT or varTable
    ov.names.ord1 <- lav_partable_vnames(FLAT, type="ov.ord")
    # check if we have "|" for exogenous variables
    if(length(ov.names.ord1) > 0L) {
        idx <- which(ov.names.ord1 %in% ov.names.x)
        if(length(idx) > 0L) {
            warning("lavaan WARNING: thresholds are defined for exogenous variables: ", paste(ov.names.ord1[idx], collapse=" "))
        }
    }
 
    if(!is.null(varTable)) {
        ov.names.ord2 <- as.character(varTable$name[ varTable$type == "ordered" ])
        # remove fixed.x variables
        idx <- which(ov.names.ord2 %in% ov.names.x)
        if(length(idx) > 0L)
            ov.names.ord2 <- ov.names.ord2[-idx]
    } else {
        ov.names.ord2 <- character(0)
    }
    #### FIXME!!!!! ORDER!
    ov.names.ord <- unique(c(ov.names.ord1, ov.names.ord2))

    # if we have the "|" in the model syntax, check the number of thresholds
    if(!is.null(varTable) && length(ov.names.ord1) > 0L) {
        for(o in ov.names.ord1) {
            nth <- varTable$nlev[ varTable$name == o ] - 1L
            nth.in.partable <- sum(FLAT$op == "|" & FLAT$lhs == o)
            if(nth != nth.in.partable) {
                stop("lavaan ERROR: expected ", nth, 
                     " threshold(s) for variable ",
                     sQuote(o), "; syntax contains ", nth.in.partable, "\n")
            }
        }
    }

    if(length(ov.names.ord) > 0L)
        categorical <- TRUE

    lhs <- rhs <- character(0)

    # 1. THRESHOLDS (based on varTable)
    #    NOTE: new in 0.5-18: ALWAYS include threshold parameters in partable,
    #    but only free them if auto.th = TRUE
    nth <- 0L
    #if(auto.th && length(ov.names.ord2) > 0L) {
    if(length(ov.names.ord2) > 0L) {
        for(o in ov.names.ord2) {
            nth  <- varTable$nlev[ varTable$name == o ] - 1L
            if(nth < 1L) next
            lhs <- c(lhs, rep(o, nth))
            rhs <- c(rhs, paste("t", seq_len(nth), sep=""))
        }
        nth <- length(lhs)
    }

    # 2. default (residual) variances and covariances

    # a) (residual) VARIANCES (all ov's except exo, and all lv's)
    # NOTE: change since 0.5-17: we ALWAYS include the vars in the 
    #       parameter table; but only if auto.var = TRUE, we set them free
    #if(auto.var) {
        ov.var <- ov.names.nox
        # auto-remove ordinal variables
        #idx <- match(ov.names.ord, ov.var)
        #if(length(idx)) ov.var <- ov.var[-idx]
        lhs <- c(lhs, ov.var, lv.names)
        rhs <- c(rhs, ov.var, lv.names)
    #}

    # b) `independent` latent variable COVARIANCES (lv.names.x)
    if(auto.cov.lv.x && length(lv.names.x) > 1L) {
        tmp <- utils::combn(lv.names.x, 2)
        lhs <- c(lhs, tmp[1,]) # to fill upper.tri
        rhs <- c(rhs, tmp[2,])
    }

    # c) `dependent` latent variables COVARIANCES (lv.y.idx + ov.y.lv.idx)
    if(auto.cov.y && length(lvov.names.y) > 1L) {
        tmp <- utils::combn(lvov.names.y, 2L)
        lhs <- c(lhs, tmp[1,]) # to fill upper.tri
        rhs <- c(rhs, tmp[2,])
    }

    # d) exogenous x covariates: VARIANCES + COVARIANCES
    if(!categorical && (nx <- length(ov.names.x)) > 0L) {
        idx <- lower.tri(matrix(0, nx, nx), diag=TRUE)
        lhs <- c(lhs, rep(ov.names.x,  each=nx)[idx]) # fill upper.tri
        rhs <- c(rhs, rep(ov.names.x, times=nx)[idx])
    }
 
    # create 'op' (thresholds come first, then variances)
    op <- rep("~~", length(lhs)); op[seq_len(nth)] <- "|"

    # LATENT RESPONSE SCALES (DELTA)
    if(auto.delta && auto.th && length(ov.names.ord) > 0L && 
       # length(lv.names) > 0L &&
       (ngroups > 1L || any(FLAT$op == "~*~") || parameterization == "theta")) {
        lhs <- c(lhs, ov.names.ord)
        rhs <- c(rhs, ov.names.ord)
         op <- c(op,  rep("~*~", length(ov.names.ord)))
    }

    # 3. INTERCEPTS
    if(meanstructure) {
        if(categorical) {
            ov.int <- ov.names.nox
        } else {
            ov.int <- ov.names
        }
        # auto-remove ordinal variables
        #idx <- which(ov.int %in% ov.names.ord)
        #if(length(idx)) ov.int <- ov.int[-idx]

        int.lhs <- c(ov.int, lv.names)
        lhs <- c(lhs, int.lhs)
        rhs <- c(rhs, rep("",   length(int.lhs)))
        op  <- c(op,  rep("~1", length(int.lhs)))
    }

    # free group weights
    if(group.w.free) {
        lhs <- c(lhs, "group")
        rhs <- c(rhs, "w")
         op <- c(op,  "%") 
    }

    DEFAULT <- data.frame(lhs=lhs, op=op, rhs=rhs,
                          mod.idx=rep(0L, length(lhs)),
                          stringsAsFactors=FALSE)


    # 4. USER: user-specified elements
    lhs     <- FLAT$lhs
     op     <- FLAT$op
    rhs     <- FLAT$rhs
    mod.idx <- FLAT$mod.idx

    lv.names     <- lav_partable_vnames(FLAT, type="lv")     # latent variables
    ov.names     <- lav_partable_vnames(FLAT, type="ov")     # observed variables

    # check order of covariances: we only fill the upper.tri!
    cov.idx <- which(op == "~~" & lhs != rhs)
    for(i in cov.idx) {
        lv.ov.names <- c(lv.names, ov.names) ### FIXME!!! OK??
        lv.idx <- match(c(lhs[i], rhs[i]), lv.ov.names)
        if(lv.idx[1] > lv.idx[2]) { # swap!
            tmp <- lhs[i]; lhs[i] <- rhs[i]; rhs[i] <- tmp
        }
        if(lhs[i] %in% lv.names && rhs[i] %in% lv.names) {
            lv.idx <- match(c(lhs[i], rhs[i]), lv.names)
            if(lv.idx[1] > lv.idx[2]) { # swap!
                tmp <- lhs[i]; lhs[i] <- rhs[i]; rhs[i] <- tmp
            }
        } else if(lhs[i] %in% ov.names && rhs[i] %in% ov.names) {
            ov.idx <- match(c(lhs[i], rhs[i]), ov.names)
            if(ov.idx[1] > ov.idx[2]) { # swap!
                tmp <- lhs[i]; lhs[i] <- rhs[i]; rhs[i] <- tmp
            }
        } else { # mixed!! # we allow this since 0.4-10
            lv.ov.names <- c(lv.names, ov.names) ### FIXME!!! OK??
            lv.idx <- match(c(lhs[i], rhs[i]), lv.ov.names)
            if(lv.idx[1] > lv.idx[2]) { # swap!
                tmp <- lhs[i]; lhs[i] <- rhs[i]; rhs[i] <- tmp
            }
        }
    }

    USER <- data.frame(lhs=lhs, op=op, rhs=rhs, mod.idx=mod.idx,
                       stringsAsFactors=FALSE)

    # check for duplicated elements in USER
    TMP <- USER[,1:3]
    idx <- which(duplicated(TMP))
    if(length(idx) > 0L) {
        warning("duplicated elements in model syntax have been ignored: ", TMP[idx,])
        USER <- USER[-idx,]
    }

    # check for duplicated elements in DEFAULT
    # - FIXME: can we not avoid this somehow??
    # - for example, if the user model includes 'x1 ~~ x1'
    #   or 'x1 ~ 1' 
    # - remove them from DEFAULT
    TMP <- rbind(DEFAULT[,1:3], USER[,1:3])
    idx <- which(duplicated(TMP, fromLast=TRUE)) # idx should be in DEFAULT
    if(length(idx)) {
        for(i in idx) {
            flat.idx <- which(USER$lhs   == DEFAULT$lhs[i] &
                              USER$op    == DEFAULT$op[i]  &
                              USER$rhs   == DEFAULT$rhs[i])
            if(length(flat.idx) != 1L) {
                cat("[lavaan DEBUG] idx in TMP: i = ", i, "\n"); print(TMP[i,])
                cat("[lavaan DEBUG] idx in DEFAULT: i = ", i, "\n"); print(DEFAULT[i,])
               cat("[lavaan DEBUG] flat.idx:"); print(flat.idx)
            }
        }
        DEFAULT <- DEFAULT[-idx,]
    }

    # now that we have removed all duplicated elements, we can construct
    # the LIST for a single group
    lhs     <- c(USER$lhs, DEFAULT$lhs)
    op      <- c(USER$op,  DEFAULT$op)
    rhs     <- c(USER$rhs, DEFAULT$rhs)
    user    <- c(rep(1L, length(USER$lhs)),
                 rep(0L, length(DEFAULT$lhs)))
    mod.idx <- c(USER$mod.idx, DEFAULT$mod.idx)
    free    <- rep(1L,  length(lhs))
    ustart  <- rep(as.numeric(NA), length(lhs))
    #label   <- paste(lhs, op, rhs, sep="")
    label   <- rep(character(1), length(lhs))
    exo     <- rep(0L, length(lhs))

    # 0a. if auto.th = FALSE, set fix the thresholds
    if(!auto.th) {
        th.idx <- which(op == "|" & user == 0L)
        free[th.idx] <- 0L
    }

    # 0b. if auto.var = FALSE, set the unspecified variances to zero
    if(!auto.var) {
        var.idx <- which(op == "~~" &
                         lhs == rhs &
                         user == 0L)
        ustart[var.idx] <- 0.0
          free[var.idx] <- 0L
    } else {
        # 'formative' (residual) variances are set to zero by default
        var.idx <- which(op == "~~" &
                         lhs == rhs &
                         lhs %in% lv.names.f &
                         user == 0L)
        ustart[var.idx] <- 0.0
          free[var.idx] <- 0L
    }

    # 1. fix metric of regular latent variables
    if(std.lv) {
        # fix metric by fixing the variance of the latent variable
        lv.var.idx <- which(op == "~~" &
                            lhs %in% lv.names & lhs == rhs)
        ustart[lv.var.idx] <- 1.0
          free[lv.var.idx] <- 0L
    }
    if(auto.fix.first) {
        # fix metric by fixing the loading of the first indicator
        mm.idx <- which(op == "=~")
        first.idx <- mm.idx[which(!duplicated(lhs[mm.idx]))]
        ustart[first.idx] <- 1.0
          free[first.idx] <- 0L
    }

    # 2. fix residual variance of single indicators to zero
    if(auto.var && auto.fix.single) {
        mm.idx <- which(op == "=~")
        T <- table(lhs[mm.idx])
        if(any(T == 1L)) {
            # ok, we have a LV with only a single indicator
            lv.names.single <- names(T)[T == 1L]
            # get corresponding indicator if unique
            lhs.mm <- lhs[mm.idx]; rhs.mm <- rhs[mm.idx]
            single.ind <- rhs.mm[which(lhs.mm %in% lv.names.single &
                                       !(duplicated(rhs.mm) |
                                         duplicated(rhs.mm, fromLast=TRUE)))]
            # is the indicator unique?
            if(length(single.ind)) {
                var.idx <- which(op == "~~" & lhs %in% single.ind
                                            & rhs %in% single.ind
                                            & lhs == rhs
                                            & user == 0L)
                ustart[var.idx] <- 0.0
                  free[var.idx] <- 0L
            }
        }
    }

    # 3. orthogonal=TRUE?
    if(orthogonal) {
        # FIXME: only lv.x.idx for now
        lv.cov.idx <- which(op == "~~" &
                            lhs %in% lv.names &
                            lhs != rhs &
                            user == 0L)
        ustart[lv.cov.idx] <- 0.0
          free[lv.cov.idx] <- 0L
    }

    # 4. intercepts
    if(meanstructure) {
        if(categorical) {
            # zero intercepts/means ordinal variables
                   ov.int.idx <- which(op == "~1" &
                                       lhs %in% ov.names.ord &
                                       user == 0L)
            ustart[ov.int.idx] <- 0.0
              free[ov.int.idx] <- 0L
        }
        if(int.ov.free == FALSE) {
            # zero intercepts/means observed variables
                   ov.int.idx <- which(op == "~1" &
                                       lhs %in% ov.names &
                                       user == 0L)
            ustart[ov.int.idx] <- 0.0
              free[ov.int.idx] <- 0L
        }
        if(int.lv.free == FALSE) {
            # zero intercepts/means latent variables
                   lv.int.idx <- which(op == "~1" &
                                       lhs %in% lv.names &
                                       user == 0L)
            ustart[lv.int.idx] <- 0.0
              free[lv.int.idx] <- 0L
        }
    }

    # 5. handle exogenous `fixed.x' covariates
    if(length(ov.names.x) > 0 && fixed.x) {
        # 1. variances/covariances
               exo.idx  <- which(op == "~~" &
                                 rhs %in% ov.names.x &
                                 user == 0L)
        ustart[exo.idx] <- as.numeric(NA) # should be overriden later!
          free[exo.idx] <- 0L
           exo[exo.idx] <- 1L

        # 2. intercepts
               exo.int.idx  <- which(op == "~1" &
                                     lhs %in% ov.names.x &
                                     user == 0L)
        ustart[exo.int.idx] <- as.numeric(NA) # should be overriden later!
          free[exo.int.idx] <- 0L
           exo[exo.int.idx] <- 1L
    }

    # 5b. residual variances of ordinal variables?
    if(length(ov.names.ord) > 0L) {
        ord.idx <- which(lhs %in% ov.names.ord &
                         op == "~~" &
                         lhs == rhs)
        ustart[ord.idx] <- 1L ## FIXME!! or 0?? (0 breaks ex3.12)
          free[ord.idx] <- 0L
    }

    # 5c latent response scales of ordinal variables?
    if(length(ov.names.ord) > 0L) {
        delta.idx <- which(op == "~*~")
        ustart[delta.idx] <- 1.0
          free[delta.idx] <- 0L
    }

    # group proportions (group 1L)
    if(group.w.free) {
        group.idx <- which(lhs == "group" & op == "%")
        #if(ngroups > 1L) {
              free[ group.idx ] <- 1L
            ustart[ group.idx ] <- as.numeric(NA)
        #} else {
        #      free[ group.idx ] <- 0L
        #    ustart[ group.idx ] <- 0.0 # last group
        #}
    }

    # 6. multiple groups?
    group <- rep(1L, length(lhs))
    if(ngroups > 1) {
        group   <- rep(1:ngroups, each=length(lhs))
        user    <- rep(user,    times=ngroups)
        lhs     <- rep(lhs,     times=ngroups)
        op      <- rep(op,      times=ngroups)
        rhs     <- rep(rhs,     times=ngroups)
        free    <- rep(free,    times=ngroups)
        ustart  <- rep(ustart,  times=ngroups)
        mod.idx <- rep(mod.idx, times=ngroups)
        label   <- rep(label,   times=ngroups)
        exo     <- rep(exo,     times=ngroups)

        # specific changes per group
        for(g in 2:ngroups) {
            # label
            # label[group == g] <- paste(label[group == 1], ".g", g, sep="")

            # free/fix intercepts
            if(meanstructure) {
                int.idx  <- which(op == "~1" &
                                  lhs %in% lv.names &
                                  user == 0L &
                                  group == g)
                if(int.lv.free == FALSE && g > 1 &&
                   ("intercepts" %in% group.equal ||
                    "thresholds" %in% group.equal) &&
                   !("means" %in% group.equal) ) {
                      free[ int.idx ] <- 1L
                    ustart[ int.idx ] <- as.numeric(NA)
                }
            }

            # latent response scaling
            if(parameterization == "delta") {
                if(any(op == "~*~" & group == g) &&
                   ("thresholds" %in% group.equal)) {
                    delta.idx <- which(op == "~*~" & group == g)
                      free[ delta.idx ] <- 1L
                    ustart[ delta.idx ] <- as.numeric(NA)
                }
            } else if(parameterization == "theta") {
                if(any(op == "~*~" & group == g) &&
                   ("thresholds" %in% group.equal)) {
                    var.ord.idx <- which(op == "~~" & group == g &
                                         lhs %in% ov.names.ord & lhs == rhs)
                      free[ var.ord.idx ] <- 1L
                    ustart[ var.ord.idx ] <- as.numeric(NA)
                }
            } else {
                stop("lavaan ERROR: parameterization ", parameterization,
                     " not supported")
            }

            # group proportions
            if(group.w.free) {
                group.idx <- which(lhs == "group" & op == "%" & group == g)
                #if(g == ngroups) {
                #      free[ group.idx ] <- 0L
                #    ustart[ group.idx ] <- 0.0 # last group
                #} else {
                      free[ group.idx ] <- 1L
                    ustart[ group.idx ] <- as.numeric(NA)
                #}
            }
        } # g
    } # ngroups

    # construct LIST
    #LIST  <- data.frame(
    LIST   <- list(     id          = seq_along(lhs),
                        lhs         = lhs,
                        op          = op,
                        rhs         = rhs,
                        user        = user,
                        group       = group,
                        mod.idx     = mod.idx,
                        free        = free,
                        ustart      = ustart,
                        exo         = exo,
                        label       = label    ,
                        # IF we add these, also change
                        eq.id       = rep(0L,  length(lhs)),
                        unco        = rep(0L,  length(lhs))
                   )
    #                   stringsAsFactors=FALSE)

    LIST
}

lav_partable_matrixrep <- function(partable, target = NULL,
                                       representation = "LISREL") {

    if(is.null(target)) {
        target <- partable
    } 

    if(representation == "LISREL") {
        REP <- representation.LISREL(partable = partable, target = target,
                                     extra = FALSE)
    } else {
        stop("only LISREL representation has been implemented")
    }

    if(is.data.frame(target)) {
        target <- cbind(target, as.data.frame(REP, stringsAsFactors = FALSE))
    } else {
        target$mat <- REP$mat
        target$row <- REP$row
        target$col <- REP$col
    }

    target
}
    
