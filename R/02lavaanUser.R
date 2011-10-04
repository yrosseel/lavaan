# constructor for the lavaanUser model description
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


lavaanify <- function(model.syntax    = NULL, 
                      meanstructure   = FALSE,
                      int.ov.free     = FALSE,
                      int.lv.free     = FALSE,
                      orthogonal      = FALSE, 
                      std.lv          = FALSE,
                      fixed.x         = TRUE,
                      constraints     = NULL,

                      auto.fix.first  = FALSE,
                      auto.fix.single = FALSE,
                      auto.var        = FALSE,
                      auto.cov.lv.x   = FALSE,
                      auto.cov.y      = FALSE,

                      ngroups         = 1L,
                      group.equal     = NULL,
                      group.partial   = NULL,
                      debug           = FALSE,
                      warn            = TRUE,
                      
                      as.data.frame.   = TRUE) {


    # parse the model syntax and flatten the user-specified model
    # return a data.frame, where each line is a model element (rhs, op, lhs)
    FLAT <- flatten.model.syntax(model.syntax=model.syntax, warn=warn)
    # user-specified *modifiers* are returned as an attribute
    MOD  <- attr(FLAT, "modifiers"); attr(FLAT, "modifiers") <- NULL
    # user-specified *constraints* are returned as an attribute
    CON  <- attr(FLAT, "constraints"); attr(FLAT, "constraints") <- NULL

    # extra constraints?
    if(!is.null(constraints) && nchar(constraints) > 0L) {
        FLAT2 <- flatten.model.syntax(model.syntax=constraints, warn=warn)
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

    # check for `empty' model
    if(!any(c(FLAT$op == "=~", FLAT$op == "~"))) {
        # empty model; abort or warning??
        #stop("lavaan ERROR: model does not contain a measurement model",
        #     "or at least one regression formula")
        warning("lavaan WARNING: model does not contain a measurement model",
              "or at least one regression formula")
    }

    # check for meanstructure
    if(any(FLAT$op == "~1")) meanstructure <- TRUE

    # extract `names' of various types of variables:
    lv.names     <- vnames(FLAT, type="lv")     # latent variables
    ov.names     <- vnames(FLAT, type="ov")     # observed variables
    ov.names.x   <- vnames(FLAT, type="ov.x")   # exogenous x covariates 
    ov.names.nox <- vnames(FLAT, type="ov.nox")
    lv.names.x   <- vnames(FLAT, type="lv.x")   # exogenous lv
    ov.names.y   <- vnames(FLAT, type="ov.y")   # dependent ov
    lv.names.y   <- vnames(FLAT, type="lv.y")   # dependent lv
    #lvov.names.y <- c(ov.names.y, lv.names.y)
    lvov.names.y <- c(lv.names.y, ov.names.y)


    #### prepare LIST
    ####    - LIST should only contain *non-zero* parameters (fixed or free)
    ####    - includes user-specified elements, but also `default' elements:
    ####     1) (residual) variances and covariances
    ####     2) ov/lv intercepts


    ### DEFAULT elements: parameters that are typically not specified by
    ###                   users, but should typically be considered, 
    ###                   either free or fixed
    lhs <- rhs <- character(0)

    # 1. default (residual) variances and covariances

    # a) (residual) VARIANCES (all ov's, except exo and regular lv's)
    if(auto.var) {
        lhs <- c(lhs, ov.names.nox, lv.names)
        rhs <- c(rhs, ov.names.nox, lv.names)
    }

    # b) `independent` latent variable COVARIANCES (lv.names.x)
    if(auto.cov.lv.x && length(lv.names.x) > 1L) {
        tmp <- combn(lv.names.x, 2)
        lhs <- c(lhs, tmp[1,]) # to fill upper.tri
        rhs <- c(rhs, tmp[2,])
    }

    # c) `dependent` latent variables COVARIANCES (lv.y.idx + ov.y.lv.idx)
    if(auto.cov.y && length(lvov.names.y) > 1L) {
        tmp <- combn(lvov.names.y, 2L)
        lhs <- c(lhs, tmp[1,]) # to fill upper.tri
        rhs <- c(rhs, tmp[2,])
    }

    # d) exogenous x covariates: VARIANCES + COVARIANCES
    if((nx <- length(ov.names.x)) > 0L) {
        idx <- lower.tri(matrix(0, nx, nx), diag=TRUE)
        lhs <- c(lhs, rep(ov.names.x,  each=nx)[idx]) # fill upper.tri
        rhs <- c(rhs, rep(ov.names.x, times=nx)[idx])
    }
    op <- rep("~~", length(lhs))

    # 2. INTERCEPTS
    if(meanstructure) {
        int.lhs <- c(ov.names, lv.names)
        lhs <- c(lhs, int.lhs)
        rhs <- c(rhs, rep("",   length(int.lhs)))
        op  <- c(op,  rep("~1", length(int.lhs)))
    }

    DEFAULT <- data.frame(lhs=lhs, op=op, rhs=rhs, 
                          mod.idx=rep(0L, length(lhs)),
                          stringsAsFactors=FALSE)

    if(debug) {
        cat("[lavaan DEBUG]: DEFAULT (auto-added) parameters:\n")
        print(DEFAULT)
    }

    # 4. USER: user-specified elements
    lhs     <- FLAT$lhs
     op     <- FLAT$op
    rhs     <- FLAT$rhs
    mod.idx <- FLAT$mod.idx

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
        #    cat("lavaan WARNING: lhs and rhs are not BOTH in ov.names or lv.names.\n")
        #    cat("lhs = ", lhs[i], "\n")
        #    cat("rhs = ", rhs[i], "\n")
        #    cat("ov.names = "); print(ov.names)
        #    cat("lv.names = "); print(lv.names)
        #    stop("lavaan ERROR: problem in model syntax")
        }
    }

    USER <- data.frame(lhs=lhs, op=op, rhs=rhs, mod.idx=mod.idx,
                       stringsAsFactors=FALSE)

    if(debug) {
        cat("[lavaan DEBUG]: USER (user-specified) parameters:\n")
        print(USER)
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
                   "intercepts" %in% group.equal &&
                   !("means" %in% group.equal) ) {
                      free[ int.idx ] <- 1L
                    ustart[ int.idx ] <- as.numeric(NA)
                }
            }
        } # g
    } # ngroups


    # construct LIST
    #LIST  <- data.frame(
    LIST   <- list(     id          = 1:length(lhs),
                        lhs         = lhs, 
                        op          = op, 
                        rhs         = rhs, 
                        user        = user,
                        group       = group,
                        mod.idx     = mod.idx,
                        free        = free, 
                        ustart      = ustart,
                        fixed.x     = exo,
                        label       = label,
                        equal       = rep("",  length(lhs)),
                        eq.id       = rep(0L,  length(lhs)),
                        free.uncon  = rep(0L,  length(lhs))
                   )
    #                   stringsAsFactors=FALSE)
 
    # apply user-specified modifiers
    if(length(MOD)) {
        for(el in 1:length(MOD)) {
            idx <- which(LIST$mod.idx == el) # for each group
            MOD.fixed <- MOD[[el]]$fixed; MOD.start <- MOD[[el]]$start
            MOD.label <- MOD[[el]]$label; MOD.equal <- MOD[[el]]$equal

            # check for single argument if multiple groups
            if(ngroups > 1) {
                # Ok, this is not very consistent:
                # A) here we force same behavior across groups
                if(length(MOD.fixed) == 1) MOD.fixed <- rep(MOD.fixed, ngroups)
                if(length(MOD.start) == 1) MOD.start <- rep(MOD.start, ngroups)
                # B) here we do NOT!
                if(length(MOD.label) == 1) MOD.label <-
                    #c(MOD.label, paste(MOD.label, ".g", 2:ngroups, sep=""))
                    c(MOD.label, rep("", (ngroups-1L)) )
                if(length(MOD.equal) == 1) MOD.equal <- 
                    c(MOD.equal, rep("", (ngroups)) )
            }

            # check for wrong number of arguments if multiple groups
            if( (!is.null(MOD.fixed) && ngroups != length(MOD.fixed)) ||
                (!is.null(MOD.start) && ngroups != length(MOD.start)) ||
                (!is.null(MOD.label) && ngroups != length(MOD.label)) ||
                (!is.null(MOD.equal) && ngroups != length(MOD.equal)) ) {
                el.idx <- which(LIST$mod.idx == el)
                stop("lavaan ERROR: wrong number of arguments in modifier of ", LIST$lhs[el.idx], LIST$op[el.idx], LIST$rhs[el.idx])
            }

            # apply modifiers
            if(!is.null(MOD.fixed)) {
                na.idx <- which(is.na(MOD.fixed))
                not.na.idx <- which(!is.na(MOD.fixed))
                LIST$ustart[idx][not.na.idx] <- MOD.fixed[not.na.idx]
                LIST$free[  idx][not.na.idx] <- 0L
                LIST$free[  idx][na.idx]     <- 1L # eg factor loading
            }
            if(!is.null(MOD.start)) {
                LIST$ustart[idx] <- MOD.start
            }
            if(!is.null(MOD.label)) {
                LIST$label[idx] <- MOD.label
            }
            if(!is.null(MOD.equal)) {
                # NA or "" remain free
                not.na.idx <- which(!is.na(MOD.equal) & nchar(MOD.equal))
                LIST$equal[idx][not.na.idx] <- MOD.equal[not.na.idx]
                LIST$free[ idx][not.na.idx] <- 0L
            }
        }
    }
    # remove mod.idx column
    LIST$mod.idx <- NULL

    # check if CON contains *simple* equality constraints (eg b1 == b2)
    # FIXME!!!
    # b1 == b3
    # b2 == b3 does not work
    #if(length(CON) > 0L) {
    #    el.idx <- integer(0L)
    #    for(el in 1:length(CON)) {
    #        if(CON[[el]]$op == "==") {
    #            lhs.idx <- which(LIST$label == CON[[el]]$lhs)
    #            rhs.idx <- which(LIST$label == CON[[el]]$rhs)
    #            if(length(lhs.idx) && length(rhs.idx)) {
    #                # flag this constraint (to be removed)
    #                el.idx <- c(el.idx, el)
    #                # fill in equal and fix
    #                if(LIST$equal[lhs.idx] == "") {
    #                    LIST$equal[rhs.idx] <- LIST$label[lhs.idx]
    #                } else {
    #                    LIST$equal[rhs.idx] <- LIST$equal[lhs.idx]
    #                }
    #                LIST$free[ rhs.idx] <- 0L
    #            }
    #        }
    #    }
    #    if(length(el.idx) > 0L) CON <- CON[-el.idx]
    #}

    # group.equal and group.partial
    if(ngroups > 1 && length(group.equal) > 0) {
 
        # create `default' labels for all parameters
        LABEL <- getParameterLabels(LIST)

        # LOADINGS
        if("loadings" %in% group.equal) {
            g1.idx <- which(LIST$op == "=~" & # LIST$free > 0 &
                            LIST$group == 1)
            for(g in 2:ngroups) {
                idx <- which(LIST$op == "=~" & # LIST$free > 0 &
                             LIST$group == g)
                LIST$free[ idx] <- 0L
                LIST$equal[idx] <- LABEL[g1.idx]
            }
        }

        # INTERCEPTS (OV)
        if("intercepts" %in% group.equal) {
            g1.idx <- which(LIST$op == "~1"  & # LIST$free > 0 &
                            LIST$group == 1  & LIST$lhs %in% ov.names.nox)
            for(g in 2:ngroups) {
                idx <- which(LIST$op == "~1"  & # LIST$free > 0 &
                             LIST$group == g  & LIST$lhs %in% ov.names.nox)
                LIST$free[ idx] <- 0L
                LIST$equal[idx] <- LABEL[g1.idx]
            }
        }

        # MEANS (LV)
        if("means" %in% group.equal) {
            g1.idx <- which(LIST$op == "~1" & # LIST$free > 0 &
                            LIST$group == 1 &
                            LIST$lhs %in% lv.names)
            for(g in 2:ngroups) {
                idx <- which(LIST$op == "~1" & # LIST$free > 0 &
                             LIST$group == g &
                             LIST$lhs %in% lv.names)
                LIST$free[ idx] <- 0L
                LIST$equal[idx] <- LABEL[g1.idx]
            }
        }
     
        # REGRESSIONS
        if("regressions" %in% group.equal) {
            g1.idx <- which(LIST$op == "~" & # LIST$free > 0 &
                            LIST$group == 1)
            for(g in 2:ngroups) {
                idx <- which(LIST$op == "~" & # LIST$free > 0 &
                             LIST$group == g)
                LIST$free[ idx] <- 0L
                LIST$equal[idx] <- LABEL[g1.idx]
            }
        }

        # RESIDUAL variances (FIXME: OV ONLY!)
        if("residuals" %in% group.equal) {
            g1.idx <- which(LIST$op == "~~" & # LIST$free > 0 &
                            LIST$group == 1 &
                            LIST$lhs %in% ov.names.nox &
                            LIST$lhs == LIST$rhs)
            for(g in 2:ngroups) {
                idx <- which(LIST$op == "~~" & # LIST$free > 0 &
                             LIST$group == g &
                             LIST$lhs %in% ov.names.nox &
                             LIST$lhs == LIST$rhs)
                LIST$free[ idx] <- 0L
                LIST$equal[idx] <- LABEL[g1.idx]
            }
        }

        # RESIDUAL covariances (FIXME: OV ONLY!)
        if("residual.covariances" %in% group.equal) {
            g1.idx <- which(LIST$op == "~~" & # LIST$free > 0 &
                            LIST$group == 1 &
                            LIST$lhs %in% ov.names.nox &
                            LIST$lhs != LIST$rhs)
            for(g in 2:ngroups) {
                idx <- which(LIST$op == "~~" & # LIST$free > 0 &
                             LIST$group == g &
                             LIST$lhs %in% ov.names.nox &
                             LIST$lhs != LIST$rhs)
                LIST$free[ idx] <- 0L
                LIST$equal[idx] <- LABEL[g1.idx]
            }
        }

        # LV VARIANCES
        if("lv.variances" %in% group.equal) {
            g1.idx <- which(LIST$op == "~~" &
                            LIST$group == 1 &
                            LIST$lhs %in% lv.names &
                            LIST$lhs == LIST$rhs)
            for(g in 2:ngroups) {
                idx <- which(LIST$op == "~~" & # LIST$free > 0 &
                             LIST$group == g &
                             LIST$lhs %in% lv.names &
                             LIST$lhs == LIST$rhs)
                LIST$free[ idx] <- 0L
                LIST$equal[idx] <- LABEL[g1.idx]
            }
        }

        # LV COVARIANCES
        if("lv.covariances" %in% group.equal) {
            g1.idx <- which(LIST$op == "~~" & # LIST$free > 0 &
                            LIST$group == 1 &
                            LIST$lhs %in% lv.names &
                            LIST$lhs != LIST$rhs)
            for(g in 2:ngroups) {
                idx <- which(LIST$op == "~~" & # LIST$free > 0 &
                             LIST$group == g &
                             LIST$lhs %in% lv.names &
                             LIST$lhs != LIST$rhs)
                LIST$free[ idx] <- 0L
                LIST$equal[idx] <- LABEL[g1.idx]
            }
        }

        # PARTIAL: undo constraints for user-specified parameters
        if(length(group.partial) > 0) {
            # strip white space
            group.partial <- gsub("[[:space:]]+", "", group.partial)

            # allow labels without a group prefix
            # add group prefix for all groups > 1
            group.partial.orig <- group.partial
            group.partial <- character(0)
            for(v in group.partial.orig) {
                # prefix?
                if(grepl("*.g[0-9]$", v)) {
                    group.partial <- c(group.partial, v)
                } else {
                    # add all groups > 1
                    for(g in 2:ngroups) {
                        v2 <- paste(v, ".g", 2:ngroups, sep="")
                        group.partial <- c(group.partial, v2)
                    }
                }
            }

            # check if all group.partial names are found
            if(warn) {
                idx.not <- which(!group.partial %in% LABEL)
                if(length(idx.not) > 0L) 
                    warning("lavaan WARNING: some parameter names in group.partial are not found:", group.partial[idx.not])
            }   

            # free up this parameter
            for(g in 2:ngroups) {
                idx <- which(LABEL %in% group.partial)
                LIST$free[idx] <- 1L
                LIST$equal[idx] <- ""
            }

        }
    }

    # handle equality constraints
    # two types: 

    # 1. those with the same (non-empty) label (user-specified only!)
    idx.eq.label <- which(nchar(LIST$label) & duplicated(LIST$label))
    if(length(idx.eq.label) > 0) {
        for(idx in idx.eq.label) {
            eq.label <- LIST$label[idx]
            ref.idx <- which(LIST$label == eq.label)[1] # the first one only
            LIST$eq.id[idx] <- ref.idx
            LIST$free[idx] <- 0L  # fix!
        }
    }


    # 2. those explicitly marked with the equal() modifier
    idx.eq.label <- which(LIST$equal != "")
    if(length(idx.eq.label) > 0) {
        if(!exists("LABEL")) LABEL <- getParameterLabels(LIST)
        for(idx in idx.eq.label) {
            eq.label <- LIST$equal[idx]
            ref.idx <- which(LABEL == eq.label)[1] # the first one only
            if(length(ref.idx) == 0) {
                stop("equality label [", eq.label, "] not found")
            }
            LIST$eq.id[idx] <- ref.idx
        }
    }
    # remove label and equal columns
    #LIST$equal <- NULL

    # count free parameters
    idx.free <- which(LIST$free > 0)
    LIST$free[idx.free] <- 1:length(idx.free)

    # insert 'target/reference' id's in eq.id/eq.label columns
    # so the eq.id/eq.label columns reflect sets of equal parameters
    # the 'reference' parameter has free > 0; the 'equal to' pars free == 0
    idx <- which(LIST$eq.id > 0)
    id.idx <- unique(  LIST$eq.id[idx] )
    LIST$eq.id[id.idx] <- id.idx

    # 2. add free counter to this element
    idx.equal <- which(LIST$eq.id > 0)
    LIST$free[idx.equal] <- LIST$free[ LIST$eq.id[idx.equal] ]

    # 3. which parameters would be free without equality constraints?
    idx.free.uncon <- which(LIST$free > 0)
    LIST$free.uncon[idx.free.uncon] <- 1:length(idx.free.uncon)


    # handle constraints (if any) (NOT per group, but overall - 0.4-11)
    if(length(CON) > 0L) {
        #cat("DEBUG:\n"); print(CON)
        lhs = unlist(lapply(CON, "[[", "lhs"))
         op = unlist(lapply(CON, "[[",  "op"))
        rhs = unlist(lapply(CON, "[[", "rhs"))
        LIST$id         <- c(LIST$id,         max(LIST$id) + 1:length(lhs) )
        LIST$lhs        <- c(LIST$lhs,        lhs)
        LIST$op         <- c(LIST$op,         op)
        LIST$rhs        <- c(LIST$rhs,        rhs)
        LIST$user       <- c(LIST$user,       rep(1L, length(lhs)) )
        LIST$group      <- c(LIST$group,      rep(0L, each=length(CON)))
        LIST$free       <- c(LIST$free,       rep(0L, length(lhs)) )
        LIST$ustart     <- c(LIST$ustart,     rep(as.numeric(NA), length(lhs)))
        LIST$fixed.x    <- c(LIST$fixed.x,    rep(0L, length(lhs)) )
        LIST$label      <- c(LIST$label,      rep("",  length(lhs)) )
        LIST$equal      <- c(LIST$equal,      rep("",  length(lhs)) )
        LIST$eq.id      <- c(LIST$eq.id,      rep(0L,  length(lhs)) )
        LIST$free.uncon <- c(LIST$free.uncon, rep(0L,  length(lhs)) )
    }

    # put lhs of := elements in label column
    def.idx <- which(LIST$op == ":=")
    LIST$label[def.idx] <- LIST$lhs[def.idx]


    if(debug) { 
        cat("[lavaan DEBUG] lavaanUser\n")
        print( as.data.frame(LIST) )
    }

    # data.frame?
    if(as.data.frame.) LIST <- as.data.frame(LIST, stringsAsFactors = FALSE)

    LIST
}

flatten.model.syntax <- function(model.syntax='', warn=TRUE, debug=FALSE) {

    # check for empty syntax
    if(length(model.syntax) == 0) {
        stop("lavaan ERROR: empty model syntax")
    }

    # replace semicolons by newlines
    model.syntax <- gsub(";","\n", model.syntax)

    # break up in lines 
    model <- unlist( strsplit(model.syntax, "\n") )

    # remove comments starting with '#'
    model <- gsub("#.*","", model)

    # remove comments starting with '!'
    model <- gsub("!.*","", model)

    # strip all white space
    model <- gsub("[[:space:]]+", "", model)

    # remove empty lines
    idx <- which(nzchar(model))
    model <- model[idx]

    # check for multi-line formulas: they contain no "~" or "=" character
    # but before we do that, we remove all modifiers
    # to avoid confusion with for example equal("f1=~x1") statements
    model.simple <- gsub("\\(.*\\)\\*", "MODIFIER*", model)

    start.idx <- grep("[~=<>]", model.simple)
    end.idx <- c( start.idx[-1]-1, length(model) )
    model.orig    <- model
    model <- character( length(start.idx) )
    for(i in 1:length(start.idx)) {
        model[i] <- paste(model.orig[start.idx[i]:end.idx[i]], collapse="")
    }

    # ok, in all remaining lines, we should have a '~' operator
    # OR one of '=', '<' '>' outside the ""
    model.simple <- gsub("\\\".[^\\\"]*\\\"", "LABEL", model)
    idx.wrong <- which(!grepl("[~=<>]", model.simple))
    if(length(idx.wrong) > 0) {
        cat("lavaan: missing operator in formula(s):\n")
        print(model[idx.wrong])
        stop("lavaan ERROR: syntax error in lavaan model syntax")
    }
  
    if(debug) {
       cat("DEBUG: before gsub label():\n")
       cat(model, "\n\n")
    }

    # EXPERIMENTAL: auto-add label(" ") around literal labels
    # for example in   y ~ b1*x1 + b2*x2 + b3*x3
    # but only if the formula contains the "~" operator (excluding constraints)
    idx.formula <- which(grepl("[~]", model))
    # first, we replace NA* by as.numeric(NA)*
    model[idx.formula] <-
        gsub("(NA)\\*", "as.numeric(NA)\\*", model[idx.formula])
    # second, we replace ab* by label("ab")*
    model[idx.formula] <- 
        gsub("([[:alpha:]][[:alnum:]]*)\\*", "label(\"\\1\")\\*", model[idx.formula])

    if(debug) {
       cat("DEBUG: after gsub label():\n")
       cat(model, "\n\n")
    }



    # main operation: flatten formulas into single bivariate pieces
    # with a left-hand-side (lhs), an operator (eg "=~"), and a 
    # right-hand-side (rhs)
    # both lhs and rhs can have a modifier 
    # (but we ignore the lhs modifier for now)
    FLAT.lhs         <- character(0)
    #FLAT.lhs.mod    <- character(0)
    FLAT.op          <- character(0)
    FLAT.rhs         <- character(0)
    FLAT.rhs.mod.idx <- integer(0)
    FLAT.fixed       <- character(0)  # only for display purposes! 
    FLAT.start       <- character(0)  # only for display purposes!
    FLAT.label       <- character(0)  # only for display purposes!
    FLAT.equal       <- character(0)  # only for display purposes!
    FLAT.idx <- 0L
    MOD.idx  <- 0L
    CON.idx  <- 0L
    MOD <- vector("list", length=0L)
    CON <- vector("list", length=0L)
    for(i in 1:length(model)) {
        x <- model[i]

        # 1. which operator is used?
        line.simple <- gsub("\\\".[^\\\"]*\\\"", "LABEL", x)
        # 1.a "=~" operator?
        if(grepl("=~", line.simple, fixed=TRUE)) {
            op <- "=~"
        # 1b. "~~" operator?
        } else if(grepl("~~", line.simple, fixed=TRUE)) {
            op <- "~~"
        # 1c. "~" operator?
        } else if(grepl("~", line.simple, fixed=TRUE)) {
            op <- "~"           
        # 1d, "==" operator?
        } else if(grepl("==", line.simple, fixed=TRUE)) {
            op <- "=="  
        # 1e, "<" operator?
        } else if(grepl("<", line.simple, fixed=TRUE)) {
            op <- "<"
        # 1f, ">" operator?
        } else if(grepl(">", line.simple, fixed=TRUE)) {
            op <- ">"
        # 1g, ":=" operator?
        } else if(grepl(":=", line.simple, fixed=TRUE)) {
            op <- ":="
        } else {
            stop("unknown operator in ", model[i])
        }

        # 2. split by operator (only the *first* occurence!)
        # check first if equal/label modifier has been used on the LEFT!
        if(substr(x,1,5) == "equal") stop("equal modifier can not be used on the left-hand side of the operator")
        if(substr(x,1,5) == "label") stop("label modifier can not be used on the left-hand side of the operator")
        op.idx <- regexpr(op, x)
        lhs <- substr(x, 1, op.idx-1)
        rhs <- substr(x, op.idx+attr(op.idx, "match.length"), nchar(x))

        # 2b. if operator is "==" or "<" or ">" or ":=", put it in CON
        if(op == "==" || op == "<" || op == ">" || op == ":=") {
            CON.idx <- CON.idx + 1L
            CON[[CON.idx]] <- list(op=op, lhs=lhs, rhs=rhs)
            next
        }

        # 3. parse left hand
        #    lhs modifiers will be ignored for now
        lhs.formula <- as.formula(paste("~",lhs))
        out <- parse.rhs(rhs=lhs.formula[[2]])
        lhs.names <- names(out)

        # 4. parse rhs (as rhs of a single-sided formula)
        rhs.formula <- as.formula(paste("~",rhs))
        out <- parse.rhs(rhs=rhs.formula[[2]])

        # for each lhs element
        for(l in 1:length(lhs.names)) {

            # for each rhs element
            for(j in 1:length(out)) {

                # catch intercepts
                if((op == "~" || op == "~1") && names(out)[j] == "intercept") {
                    op <- "~1"
                    rhs.name <- ""
                } else {
                    rhs.name <- names(out)[j]
                }

                # check if we not already have this combination
                # 1. asymmetric (=~, ~, ~1)
                idx <- which(FLAT.lhs == lhs.names[l] &
                             FLAT.op  == op &
                             FLAT.rhs == rhs.name)
                if(length(idx) > 0) {
                    stop("lavaan ERROR: duplicate model element in: ", model[i])
                }
                # 2. symmetric (~~)
                idx <- which(FLAT.lhs == rhs.name &
                             FLAT.op  == "~~" &
                             FLAT.rhs == lhs.names[l])
                if(length(idx) > 0) {
                    stop("lavaan ERROR: duplicate model element in: ", model[i])
                }

                FLAT.idx <- FLAT.idx + 1
                FLAT.lhs[FLAT.idx] <- lhs.names[l]
                FLAT.op[ FLAT.idx] <- op
                FLAT.rhs[FLAT.idx] <- rhs.name
                FLAT.fixed[FLAT.idx] <- ""
                FLAT.start[FLAT.idx] <- ""
                FLAT.label[FLAT.idx] <- ""
                FLAT.equal[FLAT.idx] <- ""

                mod <- list()
                rhs.mod <- 0L
                if(length(out[[j]]$fixed) > 0) {
                    mod$fixed <- out[[j]]$fixed
                    FLAT.fixed[FLAT.idx] <- paste(mod$fixed, collapse=";")
                    rhs.mod <- 1L
                }
                if(length(out[[j]]$start) > 0) {
                    mod$start <- out[[j]]$start
                    FLAT.start[FLAT.idx] <- paste(mod$start, collapse=";")
                    rhs.mod <- 1L
                }
                if(length(out[[j]]$label) > 0) {
                    mod$label <- out[[j]]$label
                    FLAT.label[FLAT.idx] <- paste(mod$label, collapse=";")
                    rhs.mod <- 1L
                }
                if(length(out[[j]]$equal) > 0) {
                    mod$equal <- out[[j]]$equal
                    FLAT.equal[FLAT.idx] <- paste(mod$equal, collapse=";")
                    rhs.mod <- 1L
    
                    # if there is also a 'fixed', remove it with a warning
                    if(length(mod$fixed) > 0) {
                        if(warn) {
                            warning("lavaan WARNING: if an equality constraint is placed on a parameter,\n",
                                    "                you cannot also fix the value of that parameter:\n\n", 
                                    "                ", rhs)
                        } 
                        mod$fixed <- NULL
                        FLAT.fixed[FLAT.idx] <- ""
                    }
                }
                if(op == "~1" && rhs == "0") {
                    mod$fixed <- 0
                    FLAT.fixed[FLAT.idx] <- paste(mod$fixed, collapse=";")
                    rhs.mod <- 1L
                }
                if(op == "=~" && rhs == "0") {
                    mod$fixed <- 0
                    FLAT.rhs[FLAT.idx] <- FLAT.lhs[FLAT.idx]
                    FLAT.fixed[FLAT.idx] <- paste(mod$fixed, collapse=";")
                    rhs.mod <- 1L
                }

                FLAT.rhs.mod.idx[FLAT.idx] <- rhs.mod

                if(rhs.mod > 0) {
                    MOD.idx <- MOD.idx + 1
                    MOD[[MOD.idx]] <- mod
                }
            } # rhs elements
        } # lhs elements
    } # model elements

    # enumerate modifier indices
    mod.idx <- which(FLAT.rhs.mod.idx > 0)
    FLAT.rhs.mod.idx[ mod.idx ] <- 1:length(mod.idx)

    FLAT <- data.frame(lhs=FLAT.lhs, op=FLAT.op, rhs=FLAT.rhs,
                       mod.idx=FLAT.rhs.mod.idx,
                       fixed=FLAT.fixed, start=FLAT.start,
                       label=FLAT.label, equal=FLAT.equal,
                       stringsAsFactors=FALSE)
    attr(FLAT, "modifiers") <- MOD
    attr(FLAT, "constraints") <- CON

    FLAT
}

parse.rhs <- function(rhs, debug=FALSE, warn=TRUE) {

    rhs.names <- all.vars(rhs, unique=FALSE)
    nels <- length(rhs.names)
    var.names <- all.vars(rhs, unique=TRUE)
    nvar <- length(var.names)

    if(nels == 0) {
        rhs.names <- "intercept"
        nels <- 1
        var.names <- "intercept"
        nvar <- 1
    }

    if(debug) {
        cat("[lavaan DEBUG] number of elements in the expression: ", nels, "\n")
        cat("[lavaan DEBUG] number of unique variables: ", nvar, "\n")
        cat("[lavaan DEBUG] names: ", var.names, "\n")
    }

    # default
    out <- vector("list", length=nvar); names(out) <- var.names

    # parse modifiers
    for(i in nels:1) {
        if(debug) {
            cat("[lavaan DEBUG] element = ", i," rhs = \n")
            print(as.list(rhs)); cat("\n")
        }
        if(i > 1) {
            operator <- rhs[[1]]
            if(operator != "+" ) {
                print(as.list(rhs))
                stop("lavaan ERROR: i'm confused: main operator is not '+' in ",rhs)
            }
            last.term <- rhs[[3]]
        } else {
            last.term <- rhs
        }
        ii <- which(var.names %in% rhs.names[i])
        if(debug) {
            cat("[lavaan DEBUG] last.term = ")
            print(as.list(last.term)); cat("\n")
        }
        if( length(last.term) == 3 ) {  # there is at least one modifier!
            if(last.term[[1]] != "*") { # and the last operator must be a '*'
                stop("lavaan ERROR: i'm confused: last term operator is not '*'\n in ", rhs)
            }
            # extract coefficient
            cof <- try( eval(last.term[[2]], envir=NULL, enclos=NULL),
                       silent=TRUE)
            if( is.numeric(cof) ) {
                out[[ii]]$fixed <- cof
            } else if( is.na(cof) ) {
                out[[ii]]$fixed <- as.numeric(NA)
            } else if( length(last.term[[2]]) > 2) { # a multiple modifiers!
                cat("lavaan WARNING: \n")
                cat("multiple modifiers not allowed in\n"); print(rhs)
            } else if( length(last.term[[2]]) == 2 ) { # a single modifier
                if( last.term[[2]][[1]] == "start" ) {
                    # extract number(s) inside start()
                    cof <- eval( last.term[[2]][[2]], envir=NULL, enclos=NULL )
                    out[[ii]]$start <- cof
                } else if( last.term[[2]][[1]] == "equal" ) {
                    # extract label(s) inside equal()
                    cof <- eval( last.term[[2]][[2]], envir=NULL, enclos=NULL )
                    out[[ii]]$equal <- cof
                } else if( last.term[[2]][[1]] == "label" ) {
                    label <- eval( last.term[[2]][[2]], envir=NULL,enclos=NULL )
                    out[[ii]]$label <- label
                } else {
                    stop("lavaan ERROR: unknown function name `",
                         last.term[[2]][[1]],
                         "' in ", rhs)
                }
            }
        } # modifier in this term

        if(i > 1) {
            rhs <- rhs[[2]]
        }

    } # nels

    out
}


