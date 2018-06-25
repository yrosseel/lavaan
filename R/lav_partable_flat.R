
lav_partable_flat <- function(FLAT = NULL,
                              blocks           = "group",
                              block.id         = NULL,
                              meanstructure    = FALSE,
                              int.ov.free      = FALSE,
                              int.lv.free      = FALSE,
                              orthogonal       = FALSE,
                              std.lv           = FALSE,
                              conditional.x    = FALSE,
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
        if(length(idx) > 0L) {
            ov.names.ord2 <- ov.names.ord2[-idx]
        }

        # remove those that do appear in the model syntax
        idx <- which(!ov.names.ord2 %in% ov.names)
        if(length(idx) > 0L) {
            ov.names.ord2 <- ov.names.ord2[-idx]
        }
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
                stop("lavaan ERROR: expected ", max(0,nth),
                     " threshold(s) for variable ",
                     sQuote(o), "; syntax contains ", nth.in.partable, "\n")
            }
        }
    }

    if(length(ov.names.ord) > 0L)
        categorical <- TRUE

    # std.lv = TRUE, group.equal includes "loadings": give warning
    if(ngroups > 1L && std.lv && "loadings" %in% group.equal) {
        # suggested by Michael Hallquist
        warning("lavaan WARNING: std.lv = TRUE forces all variances to be unity in all groups, despite group.equal = \"loadings\"")
    }

    lhs <- rhs <- character(0)

    # 1. THRESHOLDS (based on varTable)
    #    NOTE: - new in 0.5-18: ALWAYS include threshold parameters in partable,
    #            but only free them if auto.th = TRUE
    #          - only ov.names.ord2, because ov.names.ord1 are already in USER
    #            and we only need to add 'default' parameters here
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
    if((nx <- length(ov.names.x)) > 0L) {
        idx <- lower.tri(matrix(0, nx, nx), diag=TRUE)
        lhs <- c(lhs, rep(ov.names.x,  each=nx)[idx]) # fill upper.tri
        rhs <- c(rhs, rep(ov.names.x, times=nx)[idx])
    }

    # create 'op' (thresholds come first, then variances)
    op <- rep("~~", length(lhs)); op[seq_len(nth)] <- "|"

    # LATENT RESPONSE SCALES (DELTA)
    #    NOTE: - new in 0.5-19: ALWAYS include scaling parameters in partable,
    #            but only free them if auto.delta = TRUE (and parameterization
    #            is "delta"
    #if(auto.delta && auto.th && length(ov.names.ord) > 0L &&
    #   # length(lv.names) > 0L &&
    #   (ngroups > 1L || any(FLAT$op == "~*~") || parameterization == "theta")) {
    if(length(ov.names.ord) > 0L) {
        lhs <- c(lhs, ov.names.ord)
        rhs <- c(rhs, ov.names.ord)
         op <- c(op,  rep("~*~", length(ov.names.ord)))
    }

    # 3. INTERCEPTS
    if(meanstructure) {
        #if(conditional.x) {
        #    ov.int <- ov.names.nox
        #} else {
            ov.int <- ov.names
        #}
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
    USER <- data.frame(lhs=lhs, op=op, rhs=rhs, mod.idx=mod.idx,
                       stringsAsFactors=FALSE)

    # check for duplicated elements in USER
    TMP <- USER[,1:3]
    idx <- which(duplicated(TMP))
    if(length(idx) > 0L) {
        txt <- sapply(1:length(idx), function(i) {
            paste("    ", TMP[idx[i],"lhs"],
                          TMP[idx[i], "op"],
                          TMP[idx[i],"rhs"]) })
        warning("duplicated elements in model syntax have been ignored:\n",
                paste(txt, collapse = "\n"))
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
    # the LIST for a single group/block
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

    # 5. handle exogenous `x' covariates
    if(length(ov.names.x) > 0) {

        # 1. variances/covariances
        exo.var.idx  <- which(op == "~~" &
                          rhs %in% ov.names.x &
                          user == 0L)
        if(fixed.x) {
            ustart[exo.var.idx] <- as.numeric(NA) # should be overriden later!
              free[exo.var.idx] <- 0L
               exo[exo.var.idx] <- 1L
        } #else if(conditional.x) {
          #     exo[exo.var.idx] <- 1L
        #}

        # 2. intercepts
        exo.int.idx  <- which(op == "~1" &
                              lhs %in% ov.names.x &
                              user == 0L)
        if(fixed.x) {
            ustart[exo.int.idx] <- as.numeric(NA) # should be overriden later!
              free[exo.int.idx] <- 0L
               exo[exo.int.idx] <- 1L
        } #else if(conditional.x) {
          #     exo[exo.int.idx] <- 1L
        #}

        # 3. regressions ov + lv
        exo.reg.idx <- which(op == "~" &
                             lhs %in% c(lv.names, ov.names.nox) &
                             rhs %in% ov.names.x)
        if(conditional.x) {
            exo[exo.reg.idx] <- 1L
        }
    }

    # 5b. residual variances of ordinal variables?
    if(length(ov.names.ord) > 0L) {
        ord.idx <- which(lhs %in% ov.names.ord &
                         op == "~~" &
                         user == 0L & ## New in 0.6-1
                         lhs == rhs)
        ustart[ord.idx] <- 1L ## FIXME!! or 0?? (0 breaks ex3.12)
          free[ord.idx] <- 0L
    }

    # 5c latent response scales of ordinal variables?
    #    by default, all fixed to 1.0
    if(length(ov.names.ord) > 0L) {
        delta.idx <- which(op == "~*~" &
                           user == 0L) ## New in 0.6-1
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
            if(auto.delta && parameterization == "delta") {
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
    LIST <- list( id          = seq_along(lhs),
                  lhs         = lhs,
                  op          = op,
                  rhs         = rhs,
                  user        = user)

    # add block column (before group/level columns)
    if(!is.null(block.id)) {
        # only one block
        LIST$block <- rep(block.id, length(lhs))
    } else {
        # block is a combination of at least group, level, ...
        # for now, only group
        LIST$block <- group
    }

    # block columns (typically only group)
    for(block in blocks) {
        if(block == "group") {
            LIST[[ block ]] <- group
        } else {
            LIST[[block]] <- rep(0L, length(lhs))
        }
    }

    # other columns
    LIST2 <- list(mod.idx     = mod.idx,
                  free        = free,
                  ustart      = ustart,
                  exo         = exo,
                  label       = label)

    LIST <- c(LIST, LIST2)
}

