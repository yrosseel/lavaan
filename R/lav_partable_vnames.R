# lav_partable_names
#
# YR. 29 june 2013
#  - as separate file; used to be in utils-user.R
#  - lav_partable_names (aka 'vnames') allows multiple options in 'type'
#    returning them all as a list (or just a vector if only 1 type is needed)

# public version
lavNames <- function(object, type = "ov", ...) {

    if(inherits(object, "lavaan") || inherits(object, "lavaanList")) {
         partable <- object@ParTable
    } else if(class(object) == "list" ||
              inherits(object, "data.frame")) {
        partable <- object
    } else if(class(object) == "character") {
        # just a model string?
        partable <- lavParseModelString(object)
    }

    lav_partable_vnames(partable, type = type, ...)
}

# alias for backwards compatibility
lavaanNames <- lavNames

# return variable names in a partable
# - the 'type' argument determines the status of the variable: observed,
#   latent, endo/exo/...; default = "ov", but most used is type = "all"
# - the 'group' argument either selects a single group (if group is an integer)
#   or returns a list per group
# - the 'level' argument either selects a single level (if level is an integer)
#   or returns a list per level
# - the 'block' argument either selects a single block (if block is an integer)
#   or returns a list per block
lav_partable_vnames <- function(partable, type = NULL, ...,
                                warn = FALSE, ov.x.fatal = FALSE) {

    # check for empy table
    if(length(partable$lhs) == 0) return(character(0L))

    # dotdotdot
    dotdotdot <- list(...)

    type.list <- c("ov",          # observed variables (ov)
                   "ov.x",        # (pure) exogenous observed variables
                   "ov.nox",      # non-exogenous observed variables
                   "ov.model",    # modeled observed variables (joint vs cond)
                   "ov.y",        # (pure) endogenous variables (dependent only)
                   "ov.num",      # numeric observed variables
                   "ov.ord",      # ordinal observed variables
                   "ov.ind",      # observed indicators of latent variables
                   "ov.orphan",   # lonely observed intercepts/variances
                   "ov.interaction", # interaction terms (with colon)
                   "th",          # thresholds ordinal only
                   "th.mean",     # thresholds ordinal + numeric variables

                   "lv",          # latent variables
                   "lv.regular",  # latent variables (defined by =~ only)
                   "lv.formative",# latent variables (defined by <~ only)
                   "lv.x",        # (pure) exogenous variables
                   "lv.y",        # (pure) endogenous variables
                   "lv.nox",      # non-exogenous latent variables
                   "lv.nonnormal",# latent variables with non-normal indicators
                   "lv.interaction", # interaction terms

                   "eqs.y",       # y's in regression
                   "eqs.x"        # x's in regression
                  )

    # sanity check
    stopifnot(is.list(partable),
              !missing(type),
              type %in% c(type.list, "all"))

    if(length(type) == 1L && type == "all") {
        type <- type.list
    }

    # ALWAYS need `block' column -- create one if missing
    if(is.null(partable$block)) {
        partable$block <- rep(1L, length(partable$lhs))
    }

    # nblocks -- block column is integer only
    nblocks <- lav_partable_nblocks(partable)

    # per default, use full partable
    block.select <- lav_partable_block_values(partable)

    # check for ... selection argument(s)
    ndotdotdot <- length(dotdotdot)
    if(ndotdotdot > 0L) {
        dot.names <- names(dotdotdot)
        block.select <- rep(TRUE, length(partable$lhs))
        for(dot in seq_len(ndotdotdot)) {
            # selection variable?
            block.var <- dot.names[dot]
            block.val <- dotdotdot[[block.var]]
            # do we have this 'block.var' in partable?
            if(is.null(partable[[block.var]])) {

                # for historical reasons, treat "group = 1" special
                if(block.var == "group" && block.val == 1L) {
                    partable$group <- rep(1L, length(partable$lhs))
                    # remove block == 0
                    idx <- which(partable$block == 0L)
                    if(length(idx) > 0L) {
                        partable$group[idx] <- 0L
                    }
                    block.select <- ( block.select &
                                  partable[[block.var]] %in% block.val )
                } else {
                    stop("lavaan ERROR: selection variable `",
                         block.var, " not found in the parameter table.")
                }

            } else {
                if(!all(block.val %in% partable[[block.var]])) {
                    stop("lavaan ERROR: ", block.var ,
                        " column does not contain value `", block.val, "'")
                }
                block.select <- ( block.select &
                                  partable[[block.var]] %in% block.val )
            }
        } # dot
        block.select <- unique(partable$block[block.select])

        if(length(block.select) == 0L) {
            warnings("lavaan WARNING: no blocks selected.")
        }
    }



    # output: list per block
    OUT                <- vector("list", length = nblocks)
    OUT$ov             <- vector("list", length = nblocks)
    OUT$ov.x           <- vector("list", length = nblocks)
    OUT$ov.nox         <- vector("list", length = nblocks)
    OUT$ov.model       <- vector("list", length = nblocks)
    OUT$ov.y           <- vector("list", length = nblocks)
    OUT$ov.num         <- vector("list", length = nblocks)
    OUT$ov.ord         <- vector("list", length = nblocks)
    OUT$ov.ind         <- vector("list", length = nblocks)
    OUT$ov.orphan      <- vector("list", length = nblocks)
    OUT$ov.interaction <- vector("list", length = nblocks)
    OUT$th             <- vector("list", length = nblocks)
    OUT$th.mean        <- vector("list", length = nblocks)

    OUT$lv             <- vector("list", length = nblocks)
    OUT$lv.regular     <- vector("list", length = nblocks)
    OUT$lv.formative   <- vector("list", length = nblocks)
    OUT$lv.x           <- vector("list", length = nblocks)
    OUT$lv.y           <- vector("list", length = nblocks)
    OUT$lv.nox         <- vector("list", length = nblocks)
    OUT$lv.nonnormal   <- vector("list", length = nblocks)
    OUT$lv.interaction <- vector("list", length = nblocks)

    OUT$eqs.y          <- vector("list", length = nblocks)
    OUT$eqs.x          <- vector("list", length = nblocks)

    for(b in block.select) {

        # always compute lv.names
        lv.names <- unique( partable$lhs[ partable$block == b  &
                                          (partable$op == "=~" |
                                           partable$op == "<~")  ] )

        # determine lv interactions
        int.names <- unique(partable$rhs[ partable$block == b  &
                                              grepl(":", partable$rhs) ] )
        n.int <- length(int.names)
        if(n.int > 0L) {
            ok.idx <- logical(n.int)
            for(iv in seq_len(n.int)) {
                NAMES <- strsplit(int.names[iv], ":", fixed = TRUE)[[1L]]

                # three scenario's:
                # - both variables are latent (ok)
                # - both variables are observed (ignore)
                # - only one latent (warn??) -> upgrade observed to latent
                # thus if at least one is in lv.names, we treat it as a
                # latent interaction
                if(sum(NAMES %in% lv.names) > 0L) {
                    ok.idx[iv] <- TRUE
                }
            }
            lv.interaction <- int.names[ok.idx]
            lv.names <- c(lv.names, lv.interaction)
        } else {
            lv.interaction <- character(0L)
        }

        # store lv
        if("lv" %in% type) {
            OUT$lv[[b]] <- lv.names
        }

        # regular latent variables ONLY (ie defined by =~ only)
        if("lv.regular" %in% type) {
            out <- unique( partable$lhs[ partable$block == b &
                                         partable$op == "=~"   ] )
            OUT$lv.regular[[b]] <- out
        }

        # interaction terms involving latent variables (only)
        if("lv.interaction" %in% type) {
            OUT$lv.interaction[[b]] <- lv.interaction
        }

        # formative latent variables ONLY (ie defined by <~ only)
        if("lv.formative" %in% type) {
            out <- unique( partable$lhs[ partable$block == b &
                                         partable$op == "<~"   ] )
            OUT$lv.formative[[b]] <- out
        }

        # eqs.y
        if(!(length(type) == 1L &&
           type %in% c("lv", "lv.regular", "lv.nonnormal"))) {
            eqs.y <- unique( partable$lhs[ partable$block == b  &
                                           partable$op == "~"     ] )
        }

        # store eqs.y
        if("eqs.y" %in% type) {
            OUT$eqs.y[[b]] <- eqs.y
        }

        # eqs.x
        if(!(length(type) == 1L &&
           type %in% c("lv", "lv.regular", "lv.nonnormal","lv.x"))) {
            eqs.x <- unique( partable$rhs[ partable$block == b  &
                                           (partable$op == "~"  |
                                            partable$op == "<~")  ] )
        }

        # store eqs.x
        if("eqs.x" %in% type) {
            OUT$eqs.x[[b]] <- eqs.x
        }

        # v.ind -- indicators of latent variables
        if(!(length(type) == 1L &&
           type %in% c("lv", "lv.regular", "lv.nonnormal"))) {
            v.ind <- unique( partable$rhs[ partable$block == b  &
                                           partable$op == "=~"    ] )
        }

        # ov.*
        if(!(length(type) == 1L &&
             type %in% c("lv", "lv.regular", "lv.nonnormal", "lv.x","lv.y"))) {
            # 1. indicators, which are not latent variables themselves
            ov.ind <- v.ind[ !v.ind %in% lv.names ]
            # 2. dependent ov's
            ov.y <- eqs.y[ !eqs.y %in% c(lv.names, ov.ind) ]
            # 3. independent ov's
            ov.x <- eqs.x[ !eqs.x %in% c(lv.names, ov.ind, ov.y) ]
        }

        # observed variables
        # easy approach would be: everything that is not in lv.names,
        # but the main purpose here is to 'order' the observed variables
        # according to 'type' (indicators, ov.y, ov.x, orphans)
        if(!(length(type) == 1L &&
             type %in% c("lv", "lv.regular", "lv.nonnormal", "lv.x","lv.y"))) {

            # 4. orphaned covariances
            ov.cov <- c(partable$lhs[ partable$block == b &
                                      partable$op == "~~" &
                                     !partable$lhs %in% lv.names ],
                        partable$rhs[ partable$block == b &
                                      partable$op == "~~" &
                                     !partable$rhs %in% lv.names ])
            # 5. orphaned intercepts/thresholds
            ov.int <- partable$lhs[ partable$block == b &
                                    (partable$op == "~1" |
                                     partable$op == "|") &
                                    !partable$lhs %in% lv.names ]

            ov.tmp <- c(ov.ind, ov.y, ov.x)
            ov.extra <- unique(c(ov.cov, ov.int)) # must be in this order!
                                                  # so that
                                                  # lav_partable_independence
                                                  # retains the same order
            ov.names <- c(ov.tmp, ov.extra[ !ov.extra %in% ov.tmp ])
        }

        # store ov?
        if("ov" %in% type) {
            OUT$ov[[b]] <- ov.names
        }

        if("ov.ind" %in% type) {
            OUT$ov.ind[[b]] <- ov.ind
        }

        if("ov.interaction" %in% type) {
            ov.int.names <- ov.names[ grepl(":", ov.names) ]
            n.int <- length(ov.int.names)
            if(n.int > 0L) {

                ov.names.noint <- ov.names[!ov.names %in% ov.int.names]

                ok.idx <- logical(n.int)
                for(iv in seq_len(n.int)) {
                    NAMES <- strsplit(ov.int.names[iv], ":", fixed = TRUE)[[1L]]

                    # two scenario's:
                    # - both variables are in ov.names.noint (ok)
                    # - at least one variables is NOT in ov.names.noint (ignore)
                    if(all(NAMES %in% ov.names.noint)) {
                        ok.idx[iv] <- TRUE
                    }
                }
                ov.interaction <- ov.int.names[ok.idx]
            } else {
                ov.interaction <- character(0L)
            }

            OUT$ov.interaction[[b]] <- ov.interaction
        }


        # exogenous `x' covariates
        if(any(type %in% c("ov.x","ov.nox", "ov.model",
                           "th.mean","lv.nonnormal"))) {
            # correction: is any of these ov.names.x mentioned as a variance,
            #             covariance, or intercept?
            # this should trigger a warning in lavaanify()
            if(is.null(partable$user)) { # FLAT!
                partable$user <-  rep(1L, length(partable$lhs))
            }
            vars <- c( partable$lhs[ partable$block == b  &
                                     partable$op == "~1"  &
                                     partable$user == 1     ],
                       partable$lhs[ partable$block == b  &
                                     partable$op == "~~"  &
                                     partable$user == 1     ],
                       partable$rhs[ partable$block == b  &
                                     partable$op == "~~"  &
                                     partable$user == 1     ] )
            idx.no.x <- which(ov.x %in% vars)
            if(length(idx.no.x)) {
                if(ov.x.fatal) {
                   stop("lavaan ERROR: model syntax contains variance/covariance/intercept formulas\n  involving (an) exogenous variable(s): [",
                            paste(ov.x[idx.no.x], collapse=" "),
                            "];\n  Please remove them and try again.")
                }
                if(warn) {
                    warning("lavaan WARNING: model syntax contains variance/covariance/intercept formulas\n  involving (an) exogenous variable(s): [",
                            paste(ov.x[idx.no.x], collapse=" "),
                            "];\n  Please use fixed.x=FALSE or leave them alone")
                }
                ov.x <- ov.x[-idx.no.x]
            }
            ov.tmp.x <- ov.x

            # extra
            if(!is.null(partable$exo)) {
                ov.cov <- c(partable$lhs[ partable$block == b &
                                          partable$op == "~~" &
                                          partable$exo == 1L],
                            partable$rhs[ partable$block == b &
                                          partable$op == "~~" &
                                          partable$exo == 1L])
                ov.int <- partable$lhs[ partable$block == b &
                                        partable$op == "~1" &
                                        partable$exo == 1L ]
                ov.extra <- unique(c(ov.cov, ov.int))
                ov.tmp.x <- c(ov.tmp.x, ov.extra[ !ov.extra %in% ov.tmp.x ])
            }

            ov.names.x <- ov.tmp.x
        }

        # store ov.x?
        if("ov.x" %in% type) {
            OUT$ov.x[[b]] <- ov.names.x
        }

        # story ov.orphan?
        if("ov.orphan" %in% type) {
            OUT$ov.orphan[[b]] <- ov.extra
        }

        # ov's withouth ov.x
        if(any(type %in% c("ov.nox", "ov.model",
                           "th.mean", "lv.nonnormal"))) {
            ov.names.nox <- ov.names[! ov.names %in% ov.names.x ]
        }

        # store ov.nox
        if("ov.nox" %in% type) {
            OUT$ov.nox[[b]] <- ov.names.nox
        }

        # store ov.model
        if("ov.model" %in% type) {
            # if no conditional.x, this is just ov
            # else, this is ov.nox
            if(any( partable$block == b & partable$op == "~" &
                                          partable$exo == 1L )) {
                OUT$ov.model[[b]] <- ov.names.nox
            } else {
                OUT$ov.model[[b]] <- ov.names
            }
        }

        # ov's strictly ordered
        if(any(type %in% c("ov.ord", "th", "th.mean",
                           "ov.num", "lv.nonnormal"))) {
            tmp <- unique(partable$lhs[ partable$block == b &
                                        partable$op == "|" ])
            ord.names <- ov.names[ ov.names %in% tmp ]
        }

        if("ov.ord" %in% type) {
            OUT$ov.ord[[b]] <- ord.names
        }

        # ov's strictly numeric
        if(any(type %in% c("ov.num", "lv.nonnormal"))) {
            ov.num <- ov.names[! ov.names %in% ord.names ]
        }

        if("ov.num" %in% type) {
            OUT$ov.num[[b]] <- ov.num
        }

        # nonnormal lv's
        if("lv.nonnormal" %in% type) {
            # regular lv's
            lv.reg <- unique( partable$lhs[ partable$block == b &
                                            partable$op == "=~"   ] )
            if(length(lv.reg) > 0L) {
                out <- unlist( lapply(lv.reg, function(x) {
                    # get indicators for this lv
                    tmp.ind <- unique( partable$rhs[ partable$block == b &
                                                     partable$op == "=~" &
                                                     partable$lhs == x     ] )
                    if(!all(tmp.ind %in% ov.num)) {
                        return(x)
                    } else {
                        return(character(0))
                    }
                    }) )
                OUT$lv.nonnormal[[b]] <- out
            } else {
                OUT$lv.nonnormal[[b]] <- character(0)
            }
        }

        if(any(c("th","th.mean") %in% type)) {
            TH.lhs <- partable$lhs[ partable$block == b &
                                    partable$op == "|" ]
            TH.rhs <- partable$rhs[ partable$block == b &
                                    partable$op == "|" ]
        }

        # threshold
        if("th" %in% type) {
            if(length(ord.names) > 0L) {
                # return in the right order (following ord.names!)
                out <- unlist(lapply(ord.names, function(x) {
                                  idx <- which(x == TH.lhs)
                                  TH <- unique(paste(TH.lhs[idx], "|",
                                                     TH.rhs[idx], sep=""))
                                  # make sure the th's are in increasing order
                                  # sort(TH)
                                  # NO!, don't do that; t10 will be before t2
                                  # fixed in 0.6-1 (bug report from Myrsini)
                             }))
            } else {
                out <- character(0L)
            }
            OUT$th[[b]] <- out
        }

        # thresholds and mean/intercepts of numeric variables
        if("th.mean" %in% type) {
            # if fixed.x -> use ov.names.nox
            # else -> use ov.names
            if(is.null(partable$exo) || all(partable$exo == 0L)) {
                OV.NAMES <- ov.names
            } else {
                OV.NAMES <- ov.names.nox
            }
            if(length(OV.NAMES) > 0L) {
                # return in the right order (following ov.names.nox!)
                out <- unlist(lapply(OV.NAMES, function(x) {
                              if(x %in% ord.names) {
                                  idx <- which(x == TH.lhs)
                                  TH <- unique(paste(TH.lhs[idx], "|",
                                                     TH.rhs[idx], sep=""))
                                  # make sure the th's are in increasing order
                                  #sort(TH)
                              } else {
                                  x
                              }
                         }))
            } else {
                out <- character(0L)
            }
            OUT$th.mean[[b]] <- out
        }


        # exogenous lv's
        if(any(c("lv.x","lv.nox") %in% type)) {
            tmp <- lv.names[ !lv.names %in% c(v.ind, eqs.y) ]
            lv.names.x <- lv.names[ lv.names %in% tmp ]
        }

        if("lv.x" %in% type) {
            OUT$lv.x[[b]] <- lv.names.x
        }

        # dependent ov (but not also indicator or x)
        if("ov.y" %in% type) {
            tmp <- eqs.y[ !eqs.y %in% c(v.ind, eqs.x, lv.names) ]
            OUT$ov.y[[b]] <- ov.names[ ov.names %in% tmp ]
        }

        # dependent lv (but not also indicator or x)
        if("lv.y" %in% type) {
            tmp <- eqs.y[ !eqs.y %in% c(v.ind, eqs.x) &
                           eqs.y %in% lv.names ]
            OUT$lv.y[[b]] <- lv.names[ lv.names %in% tmp ]
        }

        # non-exogenous latent variables
        if("lv.nox" %in% type) {
            OUT$lv.nox[[b]] <- lv.names[! lv.names %in% lv.names.x ]
        }

    }

    # to mimic old behaviour, if length(type) == 1L
    if(length(type) == 1L) {
        OUT <- OUT[[type]]
        # to mimic old behaviour, if specific block is requested
        if(ndotdotdot == 0L) {
            OUT <- unique(unlist(OUT))
        } else if(length(block.select) == 1L) {
            OUT <- OUT[[block.select]]
        } else {
            OUT <- OUT[block.select]
        }
    } else {
        OUT <- OUT[type]
    }

    OUT
}

# alias for backward compatibility
vnames <- lav_partable_vnames
