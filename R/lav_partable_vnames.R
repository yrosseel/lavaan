# lav_partable_names 
#
# YR. 29 june 2013 
#  - as separate file; used to be in utils-user.R
#  - lav_partable_names (aka 'vnames') allows multiple options in 'type'
#    returning them all as a list (or just a vector if only 1 type is needed)

# public version
lavNames <- function(object, type = "ov", group = NULL) {

    if(inherits(object, "lavaan")) {
         partable <- object@ParTable
    } else if(class(object) == "list" ||
              inherits(object, "data.frame")) {
        partable <- object
    }

    lav_partable_vnames(partable, type = type, group = group)
}

# alias for backwards compatibility
lavaanNames <- lavNames

# return variable names in a partable
# - the 'type' argument determines the status of the variable: observed, 
#   latent, endo/exo/...; default = "ov", but most used is type = "all"
# - the 'group' argument either selects a single group (if group is an integer)
#   or returns a list per group
lav_partable_vnames <- function(partable, type = NULL, group = NULL, 
                                warn = FALSE, ov.x.fatal = FALSE) {

    # check for empy table
    if(length(partable$lhs) == 0) return(character(0L))

    type.list <- c("ov",          # observed variables (ov)
                   "ov.x",        # (pure) exogenous observed variables
                   "ov.nox",      # non-exogenous observed variables
                   "ov.model",    # modeled observed variables (joint vs cond)
                   "ov.y",        # (pure) endogenous variables (dependent only)
                   "ov.num",      # numeric observed variables
                   "ov.ord",      # ordinal observed variables
                   "ov.ind",      # observed indicators of latent variables
                   "ov.orphan",   # lonely observed intercepts/variances
                   "th",          # thresholds ordinal only
                   "th.mean",     # thresholds ordinal + numeric variables

                   "lv",          # latent variables
                   "lv.regular",  # latent variables (defined by =~ only)
                   "lv.formative",# latent variables (defined by <~ only)
                   "lv.x",        # (pure) exogenous variables
                   "lv.y",        # (pure) endogenous variables
                   "lv.nox",      # non-exogenous latent variables
                   "lv.nonnormal",# latent variables with non-normal indicators
     
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

    # if `group' is missing in partable, just add group=1L 
    if(is.null(partable$group)) {
        partable$group <- rep(1L, length(partable$lhs))
    }
    ngroups <- max(partable$group)

    # handle group argument
    group.orig <- group
    if(is.numeric(group)) {
        group <- as.integer(group)
        stopifnot(all(group %in% partable$group))
    } else if(is.null(group) || group == "list") {
        group <- seq_len(ngroups)
    }

    # output: list per group
    OUT              <- vector("list", length=ngroups)
    OUT$ov           <- vector("list", length=ngroups) 
    OUT$ov.x         <- vector("list", length=ngroups)
    OUT$ov.nox       <- vector("list", length=ngroups)
    OUT$ov.model     <- vector("list", length=ngroups)
    OUT$ov.y         <- vector("list", length=ngroups)
    OUT$ov.num       <- vector("list", length=ngroups)
    OUT$ov.ord       <- vector("list", length=ngroups)
    OUT$ov.ind       <- vector("list", length=ngroups)
    OUT$ov.orphan    <- vector("list", length=ngroups)
    OUT$th           <- vector("list", length=ngroups)
    OUT$th.mean      <- vector("list", length=ngroups)

    OUT$lv           <- vector("list", length=ngroups)
    OUT$lv.regular   <- vector("list", length=ngroups)
    OUT$lv.formative <- vector("list", length=ngroups)
    OUT$lv.x         <- vector("list", length=ngroups)
    OUT$lv.y         <- vector("list", length=ngroups)
    OUT$lv.nox       <- vector("list", length=ngroups)
    OUT$lv.nonnormal <- vector("list", length=ngroups)

    OUT$eqs.y        <- vector("list", length=ngroups)
    OUT$eqs.x        <- vector("list", length=ngroups)

    for(g in group) {

        # always compute lv.names
        lv.names <- unique( partable$lhs[ partable$group == g  &
                                          (partable$op == "=~" | 
                                           partable$op == "<~")  ] )

        # store lv
        if("lv" %in% type) {
            OUT$lv[[g]] <- lv.names
        } 

        # regular latent variables ONLY (ie defined by =~ only)
        if("lv.regular" %in% type) {
            out <- unique( partable$lhs[ partable$group == g &
                                         partable$op == "=~"   ] )
            OUT$lv.regular[[g]] <- out
        }

        # formative latent variables ONLY (ie defined by <~ only)
        if("lv.formative" %in% type) {
            out <- unique( partable$lhs[ partable$group == g &
                                         partable$op == "<~"   ] )
            OUT$lv.formative[[g]] <- out
        }

        # eqs.y
        if(!(length(type) == 1L && 
           type %in% c("lv", "lv.regular", "lv.nonnormal"))) {
            eqs.y <- unique( partable$lhs[ partable$group == g  &
                                           partable$op == "~"     ] )
        }

        # store eqs.y
        if("eqs.y" %in% type) {
            OUT$eqs.y[[g]] <- eqs.y
        }
       
        # eqs.x
        if(!(length(type) == 1L && 
           type %in% c("lv", "lv.regular", "lv.nonnormal","lv.x"))) {
            eqs.x <- unique( partable$rhs[ partable$group == g  &
                                           (partable$op == "~"  |
                                            partable$op == "<~")  ] )
        }

        # store eqs.x
        if("eqs.x" %in% type) {
            OUT$eqs.x[[g]] <- eqs.x
        }

        # v.ind -- indicators of latent variables
        if(!(length(type) == 1L && 
           type %in% c("lv", "lv.regular", "lv.nonnormal"))) {
            v.ind <- unique( partable$rhs[ partable$group == g  &
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
            ov.cov <- c(partable$lhs[ partable$group == g &
                                      partable$op == "~~" &
                                     !partable$lhs %in% lv.names ], 
                        partable$rhs[ partable$group == g &
                                      partable$op == "~~" &
                                     !partable$rhs %in% lv.names ])
            # 5. orphaned intercepts/thresholds
            ov.int <- partable$lhs[ partable$group == g &
                                    (partable$op == "~1" | 
                                     partable$op == "|") &
                                    !partable$lhs %in% lv.names ]

            ov.tmp <- c(ov.ind, ov.y, ov.x)
            ov.extra <- unique(c(ov.cov, ov.int))
            ov.names <- c(ov.tmp, ov.extra[ !ov.extra %in% ov.tmp ])
        }

        # store ov?
        if("ov" %in% type) {
            OUT$ov[[g]] <- ov.names
        }

        if("ov.ind" %in% type) {
            OUT$ov.ind[[g]] <- ov.ind
        }

        # exogenous `x' covariates
        if(any(type %in% c("ov.x","ov.nox","ov.num", "ov.model",
                           "th.mean","lv.nonnormal"))) {
            # correction: is any of these ov.names.x mentioned as a variance,
            #             covariance, or intercept? 
            # this should trigger a warning in lavaanify()
            if(is.null(partable$user)) { # FLAT!
                partable$user <-  rep(1L, length(partable$lhs))
            }
            vars <- c( partable$lhs[ partable$group == g  &
                                     partable$op == "~1"  & 
                                     partable$user == 1     ],
                       partable$lhs[ partable$group == g  &
                                     partable$op == "~~"  & 
                                     partable$user == 1     ],
                       partable$rhs[ partable$group == g  &
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
                ov.cov <- c(partable$lhs[ partable$group == g &
                                          partable$op == "~~" & 
                                          partable$exo == 1L],
                            partable$rhs[ partable$group == g &
                                          partable$op == "~~" & 
                                          partable$exo == 1L])
                ov.int <- partable$lhs[ partable$group == g &
                                        partable$op == "~1" & 
                                        partable$exo == 1L ]
                ov.extra <- unique(c(ov.cov, ov.int))
                ov.tmp.x <- c(ov.tmp.x, ov.extra[ !ov.extra %in% ov.tmp.x ])
            }

            ov.names.x <- ov.tmp.x
        }
  
        # store ov.x?
        if("ov.x" %in% type) {
            OUT$ov.x[[g]] <- ov.names.x
        }

        # story ov.orphan?
        if("ov.orphan" %in% type) {
            OUT$ov.orphan[[g]] <- ov.extra   
        }

        # ov's withouth ov.x
        if(any(type %in% c("ov.nox", "ov.num", "ov.model",
                           "th.mean", "lv.nonnormal"))) {
            ov.names.nox <- ov.names[! ov.names %in% ov.names.x ]
        }

        # store ov.nox
        if("ov.nox" %in% type) {
            OUT$ov.nox[[g]] <- ov.names.nox
        }

        # store ov.model
        if("ov.model" %in% type) {
            # if no conditional.x, this is just ov
            # else, this is ov.nox
            if(any( partable$group == g & partable$op == "~" &
                                          partable$exo == 1L )) {
                OUT$ov.model[[g]] <- ov.names.nox
            } else {
                OUT$ov.model[[g]] <- ov.names
            }
        }

        # ov's strictly ordered
        if(any(type %in% c("ov.ord", "th", "th.mean", 
                           "ov.num", "lv.nonnormal"))) {
            tmp <- unique(partable$lhs[ partable$group == g &
                                        partable$op == "|" ])
            ord.names <- ov.names[ ov.names %in% tmp ]
        }

        if("ov.ord" %in% type) {
            OUT$ov.ord[[g]] <- ord.names
        }

        # ov's strictly numeric (but no x)
        if(any(type %in% c("ov.num", "lv.nonnormal"))) {
            ov.num <- ov.names.nox[! ov.names.nox %in% ord.names ]
        }

        if("ov.num" %in% type) {
            OUT$ov.num[[g]] <- ov.num
        }

        # nonnormal lv's
        if("lv.nonnormal" %in% type) {
            # regular lv's
            lv.reg <- unique( partable$lhs[ partable$group == g &
                                            partable$op == "=~"   ] )
            out <- unlist( lapply(lv.reg, function(x) {
                # get indicators for this lv
                tmp.ind <- unique( partable$rhs[ partable$group == g &
                                                 partable$op == "=~" &
                                                 partable$lhs == x     ] )
                if(!all(tmp.ind %in% ov.num)) {
                    return(x)
                } else {
                    return(character(0))
                }
                }) )
            OUT$lv.nonnormal[[g]] <- out
        }

        if(any(c("th","th.mean") %in% type)) {
            TH.lhs <- partable$lhs[ partable$group == g &
                                    partable$op == "|" ]
            TH.rhs <- partable$rhs[ partable$group == g &
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
                                  sort(TH)
                             }))
            } else {
                out <- character(0L)
            }
            OUT$th[[g]] <- out
        }

        # thresholds and mean/intercepts of numeric variables
        if("th.mean" %in% type) {
            if(length(ov.names.nox) > 0L) {
                # return in the right order (following ov.names.nox!)
                out <- unlist(lapply(ov.names.nox, function(x) {
                              if(x %in% ord.names) {
                                  idx <- which(x == TH.lhs)
                                  TH <- unique(paste(TH.lhs[idx], "|",
                                                     TH.rhs[idx], sep=""))
                                  # make sure the th's are in increasing order
                                  sort(TH)
                              } else {
                                  x
                              }
                         }))
            } else {
                out <- character(0L)
            }
            OUT$th.mean[[g]] <- out
        }


        # exogenous lv's
        if(any(c("lv.x","lv.nox") %in% type)) {
            tmp <- lv.names[ !lv.names %in% c(v.ind, eqs.y) ]
            lv.names.x <- lv.names[ lv.names %in% tmp ]
        }

        if("lv.x" %in% type) {
            OUT$lv.x[[g]] <- lv.names.x
        }
 
        # dependent ov (but not also indicator or x)
        if("ov.y" %in% type) {
            tmp <- eqs.y[ !eqs.y %in% c(v.ind, eqs.x, lv.names) ]
            OUT$ov.y[[g]] <- ov.names[ ov.names %in% tmp ]
        }

        # dependent lv (but not also indicator or x)
        if("lv.y" %in% type) {
            tmp <- eqs.y[ !eqs.y %in% c(v.ind, eqs.x) &
                           eqs.y %in% lv.names ]
            OUT$lv.y[[g]] <- lv.names[ lv.names %in% tmp ]
        }

        # non-exogenous latent variables
        if("lv.nox" %in% type) {
            OUT$lv.nox[[g]] <- lv.names[! lv.names %in% lv.names.x ]
        }

    }

    # to mimic old behaviour, if length(type) == 1L
    if(length(type) == 1L) {
        OUT <- OUT[[type]]
        # to mimic old behaviour, if specific group is requested
        if(is.null(group.orig)) {
            OUT <- unique(unlist(OUT))
        } else if(is.numeric(group.orig) && length(group.orig) == 1L) {
            if(length(group.orig) == 1L) {
                OUT <- OUT[[group.orig]]
            } else {
                OUT <- OUT[group.orig]
            }
        }
    } else {
        OUT <- OUT[type]
    }

    OUT
}

# alias for backward compatibility
vnames <- lav_partable_vnames
