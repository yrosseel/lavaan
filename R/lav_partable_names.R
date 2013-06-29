# lav_partable_names 
#
# YR. 29 june 2013 (as separate file; used to be in utils-user.R)

# public version
lavNames <- lavaanNames <- function(object, type="ov", group=NULL) {

    if(class(object) == "lavaan") {
         partable <- object@ParTable
    } else if(class(object) == "list" ||
              class(object) == "data.frame") {
        partable <- object
    }

    lav_partable_names(partable, type=type, group=group)
}

# return variable names in a partable
# - the 'type' argument determines the status of the variable (observed, latent,#   endo/exo/...
# - the 'group' argument either selects a single group (if group is an integer)
#   or returns a list per group if omitted
lav_partable_names <- vnames <- function(partable, type = NULL, group = NULL, 
                                         warn = FALSE, ov.x.fatal = FALSE) {

    # sanity check
    stopifnot(is.list(partable), !missing(type),
              type %in% c("lv",   "ov", "lv.regular",
                          "lv.x", "ov.x", "ov.num",
                          "ov.ord", "th", "th.mean",
                          "lv.y", "ov.y",
                          "ov.nox"))

    # if `group' is missing in partable, just add group=1L 
    if(is.null(partable$group)) {
        partable$group <- rep(1L, length(partable$lhs))
    }
    ngroups <- max(partable$group)

    # handle group argument
    group.orig <- group
    if(is.numeric(group)) {
        group <- as.integer(group)
        stopifnot(group %in% partable$group)
    } else if(is.null(group) || group == "list") {
        group <- seq_len(ngroups)
    }

    # output: list per group
    OUT <- vector("list", length=ngroups)
    for(g in group) {

        # prepare some common sets
        if(type %in% c("lv", "ov", "ov.x", "lv.x", "ov.y", "lv.y")) {
            lv.names <- unique( partable$lhs[ partable$group == g  &
                                              (partable$op == "=~" | 
                                               partable$op == "<~")  ] )
        }

        # indicators
        if(type %in% c("ov", "ov.x","lv.x","ov.y","lv.y")) {
            # indicators of latent variables
            v.ind <- unique( partable$rhs[ partable$group == g  &
                                           partable$op == "=~"    ] )
        }
   
        # eqs.y
        if(type %in% c("ov","ov.x","lv.x","ov.y","lv.y")) {
            eqs.y <- unique( partable$lhs[ partable$group == g  &
                                           partable$op == "~"     ] )
        }
       
        # eqs.x
        if(type %in% c("ov","ov.x","ov.y","lv.y")) {
            eqs.x <- unique( partable$rhs[ partable$group == g  &
                                           (partable$op == "~"  |
                                            partable$op == "<~")  ] )
        }




        # regular latent variables: lhs =~ and formative variables lhs <~
        if(type == "lv") {
            out <- lv.names
        } else 

        # regular latent variables ONLY (ie defined by =~ only)
        if(type == "lv.regular") {
            out <- unique( partable$lhs[ partable$group == g &
                                         partable$op == "=~"   ] )
        } else
    
        # observed variables 
        # easy approach would be: everything that is not in lv.names,
        # but the main purpose here is to 'order' the observed variables
        # according to 'type' (indicators, ov.y, ov.x, orphans)
        if(type == "ov") {
            # 1. indicators, which are not latent variables themselves
            ov.ind <- v.ind[ !v.ind %in% lv.names ]
            # 2. dependent ov's
            ov.y <- eqs.y[ !eqs.y %in% c(lv.names, ov.ind) ]
            # 3. independent ov's
            ov.x <- eqs.x[ !eqs.x %in% c(lv.names, ov.ind, ov.y) ]
            # order is important:
            out <- c(ov.ind, ov.y, ov.x)

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
            extra <- unique(c(ov.cov, ov.int))
            extra.idx <- which(!extra %in% out)
            out <- c(out, extra[extra.idx])
        } else

        # exogenous `x' covariates
        if(type == "ov.x") {
            # 1. indicators, which are not latent variables themselves
            ov.ind <- v.ind[ !v.ind %in% lv.names ]
            # 2. dependent ov's
            ov.y <- eqs.y[ !eqs.y %in% c(lv.names, ov.ind) ]
            # 3. independent ov's
            ov.x <- eqs.x[ !eqs.x %in% c(lv.names, ov.ind, ov.y) ]

            # correction: is any of these ov.names.x mentioned as a variance,
            #             covariance, or intercept? 
            # this should trigger a warning in lavaanify()
            if(is.null(partable$user)) { # FLAT!
                vars <- c( partable$lhs[ partable$group == g  &
                                         partable$op == "~1"    ],
                           partable$lhs[ partable$group == g  &
                                         partable$op == "~~"    ],
                           partable$rhs[ partable$group == g  &
                                         partable$op == "~~"    ])
            } else {
                vars <- c( partable$lhs[ partable$group == g  &
                                         partable$op == "~1"  & 
                                         partable$user == 1     ],
                           partable$lhs[ partable$group == g  &
                                         partable$op == "~~"  & 
                                         partable$user == 1     ],
                           partable$rhs[ partable$group == g  &
                                         partable$op == "~~"  & 
                                         partable$user == 1     ] )
            }
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
            out <- ov.x

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
                extra <- unique(c(ov.cov, ov.int))
                extra.idx <- which(!extra %in% out)
                out <- c(out, extra[extra.idx])
            }
        } else

    # ov's withouth ov.x
        if(type == "ov.nox") {
            out <- lav_partable_names(partable, "ov", group = g)
            ov.names.x <- lav_partable_names(partable, "ov.x", group = g)
            idx <- which(out %in% ov.names.x)
            if(length(idx)) out <- out[-idx]
        } else

        # ov's strictly ordered
        if(type == "ov.ord") {
            ov.names <- lav_partable_names(partable, "ov", group = g)
            tmp <- unique(partable$lhs[ partable$group == g &
                                        partable$op == "|" ])
            out <- ov.names[ov.names %in% tmp]
        } else

        # ov's strictly numeric (but no x)
        if(type == "ov.num") {
            out <- lav_partable_names(partable, "ov.nox", group = g)
            ord.names <- unique(partable$lhs[ partable$group == g &
                                              partable$op == "|" ])
            idx <- which(out %in% ord.names)
            if(length(idx)) out <- out[-idx]
        } else

        # threshold
        if(type == "th") {
            ## FIXME!! do some elegantly!
            ord.names <- lav_partable_names(partable, "ov.ord", group = g)
            if(length(ord.names) > 0L) {
                lhs <- partable$lhs[ partable$group == g &
                                     partable$op == "|" ]
                rhs <- partable$rhs[ partable$group == g &
                                     partable$op == "|" ]
                TH <- unique(paste(lhs, "|", rhs, sep=""))
                # return in the right order
                out <- unlist(lapply(ord.names, 
                              function(x) paste(x, "|t", 1:length(grep(paste("^",x,"\\|",sep=""),TH)), sep="")))
            } else {
                out <- character(0L)
            }
        } else

        # thresholds and mean/intercepts of numeric variables
        if(type == "th.mean") {
            ## FIXME!! do some elegantly!
            ord.names <- lav_partable_names(partable, "ov.ord", group = g)
            all.names <- lav_partable_names(partable, "ov.nox", group = g)
            lhs <- partable$lhs[ partable$group == g &
                                 partable$op == "|" ]
            rhs <- partable$rhs[ partable$group == g &
                                 partable$op == "|" ]
            TH <- unique(paste(lhs, "|", rhs, sep=""))
            # return in the right order
            out <- unlist(lapply(all.names,
                          function(x) {
                          if(x %in% ord.names) {
                               paste(x, "|t", 
                                   1:length(grep(paste("^",x,"\\|",sep=""),TH)), sep="")
                          } else {
                              x
                          }
                          }))
        } else


        # exogenous lv's
        if(type == "lv.x") {
            tmp <- lv.names[ !lv.names %in% c(v.ind, eqs.y) ]
            # make sure order is the same as lv.names
            out <- lv.names[ which(lv.names %in% tmp) ]
        } else
 
        # dependent ov (but not also indicator or x)
        if(type == "ov.y") {
            ov.names <- lav_partable_names(partable, "ov", group = g)
            tmp <- eqs.y[ !eqs.y %in% c(v.ind, eqs.x, lv.names) ]
            # make sure order is the same as ov.names
            out <- ov.names[ which(ov.names %in% tmp) ]
        } else

        # dependent lv (but not also indicator or x)
        if(type == "lv.y") {
            tmp <- eqs.y[ !eqs.y %in% c(v.ind, eqs.x) &
                           eqs.y %in% lv.names ]
            # make sure order is the same as lv.names
            out <- lv.names[ which(lv.names %in% tmp) ]
        }

        OUT[[g]] <- out
    }

    # to mimic old behaviour
    if(is.null(group.orig)) {
        out <- unique(unlist(OUT))
    } else if(is.numeric(group.orig)) {
        out <- OUT[[group.orig]]
    } else { # group == "list"
        out <- OUT
    }

    out
}


# return 'attributes' of a lavaan partable -- generate a new set if necessary
lav_partable_attributes <- function(partable, pa=NULL) {

    if(is.null(pa)) {
        # attached to partable?
        pa <- attributes(partable)
        if(!is.null(pa$ov.names)) { 
            return(pa)
        } else {
            pa <- list()
        }
    }

    # ov.names
    pa$ov.names <- lav_partable_names(partable, type="ov", group="list")

    # lv.names
    pa$lv.names <- lav_partable_names(partable, type="lv", group="list")

    pa
}
