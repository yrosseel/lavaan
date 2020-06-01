# constructor for the ltavParTable model description
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
# - 26 April 2016: handle multiple 'blocks' (levels, classes, groups, ...)
# - 24 March 2019: handle efa sets
# - 23 May   2020: support for random slopes

lavaanify <- lavParTable <- function(

                      model            = NULL,
                      meanstructure    = FALSE,
                      int.ov.free      = FALSE,
                      int.lv.free      = FALSE,
                      orthogonal       = FALSE,
                      orthogonal.y     = FALSE,
                      orthogonal.x     = FALSE,
                      orthogonal.efa   = FALSE,
                      std.lv           = FALSE,
                      effect.coding    = "",
                      conditional.x    = FALSE,
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
                      auto.efa         = FALSE,

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
    if(!is.null(constraints) && any(nchar(constraints) > 0L)) {
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

    # bogus varTable? (if data.type == "none")
    if(!is.null(varTable)) {
        if(all(varTable$nobs == 0)) {
            varTable <- NULL
        }
    }

    # check for wrongly specified variances/covariances/intercepts
    # of exogenous variables in model syntax (if fixed.x=TRUE)
    if(fixed.x && warn) { # we ignore the groups here!
        # we only call this function for the warning message
        tmp <- lav_partable_vnames(FLAT, "ov.x", warn = TRUE); rm(tmp)
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
            auto.efa        = TRUE
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
            auto.efa        = TRUE
        }
    }

    # check for meanstructure
    if(any(FLAT$op == "~1")) {
        meanstructure <- TRUE
    }

    # check for block identifiers in the syntax (op = ":")
    n.block.flat <- length(which(FLAT$op == ":"))

    # for each non-empty `block' in n.block.flat, produce a USER
    if(n.block.flat > 0L) {

        # what are the block lhs labels?
        BLOCKS <- tolower(FLAT$lhs[FLAT$op == ":"])
        BLOCK.lhs <- unique(BLOCKS)

        # block op == ":" indices
        BLOCK.op.idx <- which(FLAT$op == ":")

        # check for wrong spelled 'group' lhs
        if(length(grep("group", BLOCK.lhs)) > 1L) {
            warning("lavaan WARNING: ambiguous block identifiers for group:",
                    "\n\t\t  ", paste(BLOCK.lhs[grep("group", BLOCK.lhs)],
                                      collapse = ", "))
        }

        # no empty :rhs fields allowed!
        if( any( nchar(FLAT$rhs[BLOCK.op.idx]) == 0L ) ) {
            empty.idx <- nchar(FLAT$rhs[BLOCK.op.idx]) == 0L
            txt <- paste(FLAT$lhs[BLOCK.op.idx][empty.idx], ":")
            stop("lavaan ERROR: syntax contains block identifiers with ",
                 "missing numbers/labels:\n\t\t", txt)
        }

        # check for 'group' (needed?)
        if("group" %in% BLOCK.lhs) {
            # how many group blocks?
            group.block.idx <- FLAT$op == ":" & FLAT$lhs == "group"
            n.group.flat <- length( unique(FLAT$rhs[group.block.idx]) )

            if(n.group.flat > 0L && n.group.flat != ngroups) {
                stop("lavaan ERROR: syntax defines ", n.group.flat, " groups; ",
                     "data (or argument ngroups) suggests ", ngroups, " groups")
            }
        }

        # split the FLAT data.frame per `block', create LIST
        # for each `block', and rbind them together, adding block columns
        FLAT <- as.data.frame(FLAT, stringsAsFactors = FALSE)
        BLOCK.op.idx <- c(BLOCK.op.idx, nrow(FLAT) + 1L)
        BLOCK.rhs <- rep("0", length(BLOCK.lhs))
        block.id <- 0L

        for(block in seq_len(n.block.flat)) {

            # fill BLOC.rhs value
            block.lhs <- FLAT$lhs[BLOCK.op.idx[block]]
            block.rhs <- FLAT$rhs[BLOCK.op.idx[block]]
            BLOCK.rhs[ which(block.lhs == BLOCK.lhs) ] <- block.rhs

            # another block identifier?
            if(BLOCK.op.idx[block+1] - BLOCK.op.idx[block] == 1L) {
                next
            }
            block.id <- block.id + 1L


            FLAT.block <- FLAT[(BLOCK.op.idx[block]+1L):(BLOCK.op.idx[block+1]-1L),]
            # rm 'block' column (if any) in FLAT.block
            FLAT.block$block <- NULL

            # new in 0.6-7: check for random slopes, add them here
            if(block.lhs == "level" &&
               block > 1L && # FIXME: multigroup,multilevel
               !is.null(FLAT$rs) &&
               any(nchar(FLAT$rs) > 0L)) {
                lv.names.rs <- unique(FLAT$rs[nchar(FLAT$rs) > 0L])
                for(i in 1:length(lv.names.rs)) {
                    # add phantom latent variable
                    TMP <- FLAT.block[1,]
                    TMP$lhs <- lv.names.rs[i]; TMP$op <- "=~"
                    TMP$rhs     <- lv.names.rs[i]
                    TMP$mod.idx <- max(FLAT$mod.idx) + i
                    TMP$fixed   <- "0"
                    TMP$start   <- ""; TMP$lower   <- ""; TMP$upper   <- ""
                    TMP$label   <- ""; TMP$prior   <- ""; TMP$efa     <- ""
                    TMP$rs      <- lv.names.rs[i]
                    FLAT.block <- rbind(FLAT.block, TMP, deparse.level = 0L)
                    MOD <- c(MOD, list(list(fixed = 0)))
                }
            }

            LIST.block <- lav_partable_flat(FLAT.block, blocks = BLOCK.lhs,
                block.id = block.id,
                meanstructure = meanstructure,
                int.ov.free = int.ov.free, int.lv.free = int.lv.free,
                orthogonal = orthogonal, orthogonal.y = orthogonal.y,
                orthogonal.x = orthogonal.x, orthogonal.efa = orthogonal.efa,
                std.lv = std.lv,
                conditional.x = conditional.x, fixed.x = fixed.x,
                parameterization = parameterization,
                auto.fix.first = auto.fix.first,
                auto.fix.single = auto.fix.single,
                auto.var = auto.var, auto.cov.lv.x = auto.cov.lv.x,
                auto.cov.y = auto.cov.y, auto.th = auto.th,
                auto.delta = auto.delta, auto.efa = auto.efa,
                varTable = varTable, group.equal = NULL,
                group.w.free = group.w.free, ngroups = 1L)
            LIST.block <- as.data.frame(LIST.block, stringsAsFactors = FALSE)

            # add block columns with current values in BLOCK.rhs
            for(b in seq_len(length(BLOCK.lhs))) {
                block.lhs <- BLOCK.lhs[b]
                block.rhs <- BLOCK.rhs[b]
                LIST.block[block.lhs] <- rep(block.rhs, length(LIST.block$lhs))
            }

            if(!exists("LIST")) {
                LIST <- LIST.block
            } else {
                LIST.block$id <- LIST.block$id + max(LIST$id)
                LIST <- rbind(LIST, LIST.block)
            }
        }
        LIST <- as.list(LIST)

        # convert block columns to integers if possible
        for(b in seq_len(length(BLOCK.lhs))) {
            block.lhs <- BLOCK.lhs[b]
            block.rhs <- BLOCK.rhs[b]
            tmp <- try(scan(text = LIST[[block.lhs]], what = integer(),
                       quiet = TRUE), silent = TRUE)
            if(inherits(tmp, "integer")) {
                 LIST[[block.lhs]] <- tmp
            }
        }

    } else {
        LIST <- lav_partable_flat(FLAT, blocks = "group",
            meanstructure = meanstructure,
            int.ov.free = int.ov.free, int.lv.free = int.lv.free,
            orthogonal = orthogonal, orthogonal.y = orthogonal.y,
            orthogonal.x = orthogonal.x, orthogonal.efa = orthogonal.efa,
            std.lv = std.lv,
            conditional.x = conditional.x, fixed.x = fixed.x,
            parameterization = parameterization,
            auto.fix.first = auto.fix.first, auto.fix.single = auto.fix.single,
            auto.var = auto.var, auto.cov.lv.x = auto.cov.lv.x,
            auto.cov.y = auto.cov.y, auto.th = auto.th,
            auto.delta = auto.delta, auto.efa = auto.efa,
            varTable = varTable, group.equal = group.equal,
            group.w.free = group.w.free,
            ngroups = ngroups)
    }
    if(debug) {
        cat("[lavaan DEBUG]: parameter LIST without MODIFIERS:\n")
        print( as.data.frame(LIST, stringsAsFactors=FALSE) )
    }

    # handle multilevel-specific constraints
    multilevel <- FALSE
    if(!is.null(LIST$level)) {
        nlevels <- lav_partable_nlevels(LIST)
        if(nlevels > 1L) {
            multilevel <- TRUE
        }
    }
    if(multilevel && any(LIST$op == "~1")) {
        # fix ov intercepts for all within ov that also appear at level 2
        # FIXME: not tested with > 2 levels
        ov.names <- lav_partable_vnames(LIST, "ov") ## all names
        level.values <- lav_partable_level_values(LIST)
        other.names <- LIST$lhs[ LIST$op == "~1" &
                                 LIST$level %in% level.values[-1L] &
                                 LIST$lhs %in% ov.names]
        fix.names.idx <- which(LIST$op == "~1" &
                               LIST$level %in% level.values[1L] &
                               LIST$lhs %in% other.names)
        if(length(fix.names.idx) > 0L) {
            LIST$free[fix.names.idx] <- 0L
            LIST$ustart[fix.names.idx] <- 0
        }
    }
    if(multilevel && any(LIST$op == "|")) {
        # fix ALL thresholds at level 1
        level.values <- lav_partable_level_values(LIST)
        th.idx <- which(LIST$op == "|" &
                        LIST$level %in% level.values[1L])
        LIST$free[th.idx] <- 0L
        LIST$ustart[th.idx] <- 0

       # fix ALL scaling parmaters at higher levels
       scale.idx <- which(LIST$op == "~*~" &
                          LIST$level %in% level.values[-1L])
       LIST$free[scale.idx] <- 0L
       LIST$ustart[scale.idx] <- 1
    }

    # apply user-specified modifiers
    if(length(MOD)) {
        for(el in 1:length(MOD)) {
            idx <- which(LIST$mod.idx == el) # for each group

            # 0.5-21: check if idx exists
            # perhaps the corresponding element was duplicated, and removed
            if(length(idx) == 0L) {
                next
            }

            MOD.fixed <- MOD[[el]]$fixed
            MOD.start <- MOD[[el]]$start
            MOD.lower <- MOD[[el]]$lower
            MOD.upper <- MOD[[el]]$upper
            MOD.label <- MOD[[el]]$label
            MOD.prior <- MOD[[el]]$prior
            MOD.efa   <- MOD[[el]]$efa
            MOD.rs    <- MOD[[el]]$rs

            # check for single argument if multiple groups
            if(ngroups > 1L && length(idx) > 1L) {
                # Ok, this is not very consistent:
                # A) here we force same behavior across groups
                if(length(MOD.fixed) == 1L) MOD.fixed <- rep(MOD.fixed, ngroups)
                if(length(MOD.start) == 1L) MOD.start <- rep(MOD.start, ngroups)
                if(length(MOD.lower) == 1L) MOD.lower <- rep(MOD.lower, ngroups)
                if(length(MOD.upper) == 1L) MOD.upper <- rep(MOD.upper, ngroups)
                if(length(MOD.prior) == 1L) MOD.prior <- rep(MOD.prior, ngroups)
                if(length(MOD.efa)   == 1L) MOD.efa   <- rep(MOD.efa,   ngroups)
                if(length(MOD.rs)    == 1L) MOD.rs    <- rep(MOD.rs,    ngroups)
                # B) here we do NOT! otherwise, it would imply an equality
                #                    constraint...
                #    except if group.equal="loadings"!
                if(length(MOD.label) == 1L) {
                    if("loadings" %in% group.equal ||
                       "composite.loadings" %in% group.equal) {
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
                (!is.null(MOD.lower) && nidx != length(MOD.lower)) ||
                (!is.null(MOD.upper) && nidx != length(MOD.upper)) ||
                (!is.null(MOD.prior) && nidx != length(MOD.prior)) ||
                (!is.null(MOD.efa)   && nidx != length(MOD.efa))   ||
                (!is.null(MOD.rs)    && nidx != length(MOD.rs))    ||
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
            if(!is.null(MOD.efa)) {
                # do we already have a `efa' column? if not, create one
                if(is.null(LIST$efa)) {
                    LIST$efa <- character( length(LIST$lhs) )
                }
                LIST$efa[idx] <- MOD.efa
            }
            if(!is.null(MOD.rs)) {
                # do we already have a `rs' column? if not, create one
                if(is.null(LIST$rs)) {
                    LIST$rs <- character( length(LIST$lhs) )
                }
                LIST$rs[idx] <- MOD.rs

                LIST$free[  idx] <- 0L
                LIST$ustart[idx] <- as.numeric(NA) # 
            }
            if(!is.null(MOD.lower)) {
                # do we already have a `lower' column? if not, create one
                if(is.null(LIST$lower)) {
                    LIST$lower <- rep(-Inf, length(LIST$lhs) )
                }
                LIST$lower[idx] <- as.numeric(MOD.lower)
            }
            if(!is.null(MOD.upper)) {
                # do we already have a `upper' column? if not, create one
                if(is.null(LIST$upper)) {
                    LIST$upper <- rep(Inf, length(LIST$lhs) )
                }
                LIST$upper[idx] <- as.numeric(MOD.upper)
            }
            if(!is.null(MOD.label)) {
                LIST$label[idx] <- MOD.label
            }
        }
    }
    # remove mod.idx column
    LIST$mod.idx <- NULL

    # if lower/upper values were added, fix non-free values to ustart values
    # new in 0.6-6
    if(!is.null(LIST$lower)) {
        fixed.idx <- which(LIST$free == 0L)
        if(length(fixed.idx) > 0L) {
            LIST$lower[fixed.idx] <- LIST$ustart[fixed.idx]
        }
    }
    if(!is.null(LIST$upper)) {
        fixed.idx <- which(LIST$free == 0L)
        if(length(fixed.idx) > 0L) {
            LIST$upper[fixed.idx] <- LIST$ustart[fixed.idx]
        }
    }

    # if rs column is present, add RS.names to ALL rows where they are used
    if(!is.null(LIST$rs)) {
        RS.names <- unique(LIST$rs[ nchar(LIST$rs) > 0L ])
        for(i in seq_len(length(RS.names))) {
            lhs.idx <- which(LIST$lhs == RS.names[i] &
                             LIST$op == "=~")
            LIST$rs[lhs.idx] <- RS.names[i]
        }
    }

    if(debug) {
        cat("[lavaan DEBUG]: parameter LIST with MODIFIERS:\n")
        print( as.data.frame(LIST, stringsAsFactors=FALSE) )
    }

    # get 'virtual' parameter labels
    if(n.block.flat > 1L) {
        blocks <- BLOCK.lhs
    } else {
        blocks <- "group"
    }
    LABEL <- lav_partable_labels(partable = LIST,
                                 blocks = blocks,
                                 group.equal = group.equal,
                                 group.partial = group.partial)

    if(debug) {
        cat("[lavaan DEBUG]: parameter LIST with LABELS:\n")
        tmp <- LIST; tmp$LABEL <- LABEL
        print( as.data.frame(tmp, stringsAsFactors=FALSE) )
    }

    # handle EFA equality constraints
    # YR 14 Jan 2020: 0.6-6 does no longer impose 'explicit' constraints
    #                 if we only need to fix a parameter to 0/1
    # Note: we should also check if they are really needed:
    #       eg., if all the factor-loadings of the 'second' set (time/group)
    #       are constrained to be equal to the factor-loadings of the first
    #       set, no further constraints are needed
    if(auto.efa && !is.null(LIST$efa)) {
        # for each set, for each block
        nblocks <- lav_partable_nblocks(LIST)
        set.names <- lav_partable_efa_values(LIST)
        nsets <- length(set.names)

        for(b in seq_len(nblocks)) {
            for(s in seq_len(nsets)) {
                # lv's for this block/set
                lv.nam.efa <- unique(LIST$lhs[LIST$op == "=~" &
                                              LIST$block == b &
                                              LIST$efa == set.names[s]])
                if(length(lv.nam.efa) == 1L) {
                    # nothing to do (warn?)
                    next
                }

                # equality constraints on ALL factor loadings in this set?
                # --> no constraints are needed
                #if(b > 1L || s > 1L) {
                #    # collect plabels for this set, if any
                #    set.idx <- which(LIST$op == "=~" &
                #                     LIST$block == b &
                #                     LIST$lhs %in% lv.nam.efa)
                #    plabels.set <- LIST$plabel[ set.idx ]
                #
                #    # check simple == equalities
                #    # TODO: what about more complex equality constraints?
                #    eq.idx <- which(LIST$op == "==")
                #    con.labels <- c(LIST$lhs[eq.idx],
                #                    LIST$rhs[eq.idx])
                #
                #    # TODO: what about exceptions?
                #    if(all(plabels.set %in% con.labels)) {
                #        next
                #    }
                #}

                # 1. echelon pattern
                nfac <- length(lv.nam.efa)
                for(f in seq_len(nfac)) {
                    if(f == 1L) {
                        next
                    }
                    nzero <- (f - 1L)
                    ind.idx <- which(LIST$op == "=~" &
                                     LIST$block == b &
                                     LIST$lhs %in% lv.nam.efa[f])
                    if(length(ind.idx) < nzero) {
                        stop("lavaan ERROR: efa factor ", lv.nam.efa[f],
                             " has not enough indicators for echelon pattern")
                    }

                    # fix to zero
                    LIST$free[  ind.idx[seq_len(nzero)]] <- 0L
                    LIST$ustart[ind.idx[seq_len(nzero)]] <- 0
                    LIST$user[  ind.idx[seq_len(nzero)]] <- 7L

                }

                # 2. covariances constrained to zero (only if oblique rotation)
                if(!orthogonal.efa) {
                    # skip if user == 1 (user-override!)
                    cov.idx <- which(LIST$op == "~~" &
                                     LIST$block == b &
                                     LIST$user == 0L &
                                     LIST$lhs %in% lv.nam.efa &
                                     LIST$rhs %in% lv.nam.efa &
                                     LIST$lhs != LIST$rhs)

                    # fix to zero
                    LIST$free[  cov.idx] <- 0L
                    LIST$ustart[cov.idx] <- 0
                    LIST$user[  cov.idx] <- 7L
                }

            } # sets
        } # blocks
    } # auto.efa

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

            # new in 0.6-6: make sure lower/upper constraints are equal too
            if(!is.null(LIST$lower) &&
               length(unique(LIST$lower[all.idx])) > 0L) {
                non.inf <- which(is.finite(LIST$lower[all.idx]))
                if(length(non.inf) > 0L) {
                    smallest.val <- min(LIST$lower[all.idx][non.inf])
                    LIST$lower[all.idx] <- smallest.val
                }
            }
            if(!is.null(LIST$upper) &&
                length(unique(LIST$upper[all.idx])) > 0L) {
                non.inf <- which(is.finite(LIST$upper[all.idx]))
                if(length(non.inf) > 0L) {
                    largest.val <- max(LIST$upper[all.idx][non.inf])
                    LIST$upper[all.idx] <- largest.val
                }
            }

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
                #if(any(nchar(LIST$label[ref.idx])  > 0)) {
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



    # handle constraints (if any) (NOT per group, but overall - 0.4-11)
    if(length(CON) > 0L) {
        nCon <- length(CON)
        IDX <- length(LIST$id) + seq_len(nCon)
        # grow LIST with length(CON) extra rows
        LIST <- lapply(LIST, function(x) {
                    if(is.character(x)) {
                        c(x, rep("", nCon))
                    } else {
                        c(x, rep(NA, nCon))
                    } })

        # fill in some columns
        LIST$id[IDX]    <- IDX
        LIST$lhs[IDX]   <- unlist(lapply(CON, "[[", "lhs"))
        LIST$op[IDX]    <- unlist(lapply(CON, "[[",  "op"))
        LIST$rhs[IDX]   <- unlist(lapply(CON, "[[", "rhs"))
        LIST$user[IDX]  <- unlist(lapply(CON, "[[", "user"))

        # zero is nicer?
        LIST$free[IDX] <- rep(0L, nCon)
        LIST$exo[IDX]  <- rep(0L, nCon)
        LIST$block[IDX] <- rep(0L, nCon)

        if(!is.null(LIST$group)) {
            if(is.character(LIST$group)) {
                LIST$group[IDX] <- rep("", nCon)
            } else {
                LIST$group[IDX] <- rep(0L, nCon)
            }
        }
        if(!is.null(LIST$level)) {
            if(is.character(LIST$level)) {
                LIST$level[IDX] <- rep("", nCon)
            } else {
                LIST$level[IDX] <- rep(0L, nCon)
            }
        }
        if(!is.null(LIST$class)) {
            if(is.character(LIST$class)) {
                LIST$class[IDX] <- rep("", nCon)
            } else {
                LIST$class[IDX] <- rep(0L, nCon)
            }
        }
    }

    # put lhs of := elements in label column
    def.idx <- which(LIST$op == ":=")
    LIST$label[def.idx] <- LIST$lhs[def.idx]


    # handle effect.coding related equality constraints
    if(is.logical(effect.coding) && effect.coding) {
        effect.coding <- c("loadings", "intercepts")
    } else if(!is.character(effect.coding)) {
        stop("lavaan ERROR: effect.coding argument must be a character string")
    }
    if(any(c("loadings", "intercepts") %in% effect.coding)) {
        TMP <- list()
        # for each block
        nblocks <- lav_partable_nblocks(LIST)
        for(b in seq_len(nblocks)) {
            # lv's for this block/set
            lv.names <- unique(LIST$lhs[ LIST$op == "=~" &
                                         LIST$block == b ])

            if(length(lv.names) == 0L) {
                next
            }

            int.plabel <- character(0L)
            for(lv in lv.names) {

                # ind.names
                ind.names <- LIST$rhs[ LIST$op == "=~" &
                                       LIST$block == b &
                                       LIST$lhs == lv ]

                if("loadings" %in% effect.coding) {
                    # factor loadings indicators of this lv
                    loadings.idx <- which( LIST$op == "=~" &
                                           LIST$block == b &
                                           LIST$rhs %in% ind.names &
                                           LIST$lhs == lv )

                    # all free?
                    if(length(loadings.idx) > 0L &&
                       all(LIST$free[ loadings.idx ] > 0L)) {
                        # add eq constraint
                        plabel <- LIST$plabel[ loadings.idx ]

                        # Note:  we write them as
                        # .p1. == 3 - .p2. - .p3.
                        # instead of
                        # 3 ==  .p1.+.p2.+.p3.
                        # as this makes it easier to translate things to
                        # JAGS/stan

                        LHS <- plabel[1]
                        if(length(loadings.idx) > 1L) {
                            RHS <- paste(length(loadings.idx), "-",
                                         paste(plabel[-1], collapse = "-"),
                                         sep = "")
                        } else {
                            RHS <- length(loadings.idx)
                        }

                        TMP$lhs    <- c(TMP$lhs,    LHS)
                        TMP$op     <- c(TMP$op,     "==")
                        TMP$rhs    <- c(TMP$rhs,    RHS)
                        TMP$block  <- c(TMP$block,  0L)
                        TMP$user   <- c(TMP$user,   2L)
                        TMP$ustart <- c(TMP$ustart, as.numeric(NA))
                    }
                } # loadings

                if("intercepts" %in% effect.coding) {
                    # intercepts for indicators of this lv
                    intercepts.idx <- which( LIST$op == "~1" &
                                             LIST$block == b &
                                             LIST$lhs %in% ind.names)

                    # all free?
                    if(length(intercepts.idx) > 0L &&
                       all(LIST$free[ intercepts.idx ] > 0L)) {
                        # 1) add eq constraint
                        plabel <- LIST$plabel[ intercepts.idx ]

                        LHS <- plabel[1]
                        if(length(intercepts.idx) > 1L) {
                            RHS <- paste("0-",
                                         paste(plabel[-1], collapse = "-"),
                                         sep = "")
                        } else {
                            RHS <- 0L
                        }

                        TMP$lhs    <- c(TMP$lhs,    LHS)
                        TMP$op     <- c(TMP$op,     "==")
                        TMP$rhs    <- c(TMP$rhs,    RHS)
                        TMP$block  <- c(TMP$block,  0L)
                        TMP$user   <- c(TMP$user,   2L)
                        TMP$ustart <- c(TMP$ustart, as.numeric(NA))

                        # 2) release latent mean
                        lv.int.idx <- which(  LIST$op == "~1" &
                                              LIST$block == b &
                                              LIST$lhs == lv )
                        # free only if automatically added
                        if(length(lv.int.idx) > 0L &&
                           LIST$user[ lv.int.idx ] == 0L) {
                            LIST$free[ lv.int.idx ] <- 1L
                        }
                    }
                } # intercepts

            } # lv
        } # blocks

        LIST <- lav_partable_merge(LIST, TMP)
    }

    # mg.lv.variances
    if(ngroups > 1L && "mg.lv.variances" %in% effect.coding) {

        TMP <- list()

        # do not include 'EFA' lv's
        if(!is.null(LIST$efa)) {
            lv.names <- unique(LIST$lhs[ LIST$op == "=~" &
                                         !nchar(LIST$efa) > 0L ])
        } else {
            lv.names <- unique(LIST$lhs[ LIST$op == "=~" ])
        }
        group.values <- lav_partable_group_values(LIST)

        for(lv in lv.names) {
            # factor variances
            lv.var.idx <- which( LIST$op  == "~~" &
                                 LIST$lhs == lv &
                                 LIST$rhs == LIST$lhs &
                                 LIST$lhs == lv )

            # all free (but the first?)
            if(length(lv.var.idx) > 0L &&
               all(LIST$free[ lv.var.idx ][-1] > 0L)) {
                # 1) add eq constraint
                plabel <- LIST$plabel[ lv.var.idx ]

                LHS <- plabel[1]
                if(length(lv.var.idx) > 1L) {
                    RHS <- paste(length(lv.var.idx), "-",
                                 paste(plabel[-1], collapse = "-"),
                                 sep = "")
                } else {
                    RHS <- length(lv.var.idx)
                }

                TMP$lhs    <- c(TMP$lhs,    LHS)
                TMP$op     <- c(TMP$op,     "==")
                TMP$rhs    <- c(TMP$rhs,    RHS)
                TMP$block  <- c(TMP$block,  0L)
                TMP$user   <- c(TMP$user,   2L)
                TMP$ustart <- c(TMP$ustart, as.numeric(NA))

                # 2) free lv variances first group
                lv.var.g1.idx <- which( LIST$op  == "~~" &
                                        LIST$group == group.values[1] &
                                        LIST$lhs == lv &
                                        LIST$rhs == LIST$lhs &
                                        LIST$lhs == lv )
                # free only if automatically added
                if(length(lv.var.g1.idx) > 0L &&
                   LIST$user[ lv.var.g1.idx ] == 0L) {
                    LIST$free[ lv.var.g1.idx ] <- 1L
                }
            }
        } # lv

        LIST <- lav_partable_merge(LIST, TMP)
    }

    # mg.lv.efa.variances
    if(ngroups > 1L && "mg.lv.efa.variances" %in% effect.coding) {

        TMP <- list()

        # only 'EFA' lv's
        if(!is.null(LIST$efa)) {
            lv.names <- unique(LIST$lhs[ LIST$op == "=~" &
                                         nchar(LIST$efa) > 0L ])
        } else {
            lv.names <- character(0L)
        }
        group.values <- lav_partable_group_values(LIST)

        for(lv in lv.names) {
            # factor variances
            lv.var.idx <- which( LIST$op  == "~~" &
                                 LIST$lhs == lv &
                                 LIST$rhs == LIST$lhs &
                                 LIST$lhs == lv )

            # all free (but the first?)
            if(length(lv.var.idx) > 0L &&
               all(LIST$free[ lv.var.idx ][-1] > 0L)) {
                # 1) add eq constraint
                plabel <- LIST$plabel[ lv.var.idx ]

                LHS <- plabel[1]
                if(length(lv.var.idx) > 1L) {
                    RHS <- paste(length(lv.var.idx), "-",
                                 paste(plabel[-1], collapse = "-"),
                                 sep = "")
                } else {
                    RHS <- length(lv.var.idx)
                }

                TMP$lhs    <- c(TMP$lhs,    LHS)
                TMP$op     <- c(TMP$op,     "==")
                TMP$rhs    <- c(TMP$rhs,    RHS)
                TMP$block  <- c(TMP$block,  0L)
                TMP$user   <- c(TMP$user,   2L)
                TMP$ustart <- c(TMP$ustart, as.numeric(NA))

                # 2) free lv variances first group
                lv.var.g1.idx <- which( LIST$op  == "~~" &
                                        LIST$group == group.values[1] &
                                        LIST$lhs == lv &
                                        LIST$rhs == LIST$lhs &
                                        LIST$lhs == lv )
                # free only if automatically added
                if(length(lv.var.g1.idx) > 0L &&
                   LIST$user[ lv.var.g1.idx ] == 0L) {
                    LIST$free[ lv.var.g1.idx ] <- 1L
                }
            }
        } # lv

        LIST <- lav_partable_merge(LIST, TMP)
    }


    # count free parameters
    idx.free <- which(LIST$free > 0)
    LIST$free[idx.free] <- seq_along(idx.free)
    # backwards compatibility...
    if(!is.null(LIST$unco)) {
         LIST$unco[idx.free] <- seq_along(idx.free)
    }

    if(debug) {
        cat("[lavaan DEBUG] lavParTable\n")
        print( as.data.frame(LIST) )
    }

    # data.frame?
    if(as.data.frame.) {
        LIST <- as.data.frame(LIST, stringsAsFactors = FALSE)
    }

    LIST
}

