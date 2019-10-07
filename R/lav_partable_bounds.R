lav_partable_add_bounds <- function(partable       = NULL,
                                    lavpta         = NULL,
                                    lavh1          = NULL,
                                    lavdata        = NULL,
                                    lavsamplestats = NULL,
                                    lavoptions     = NULL) {

    # no support from effect.coding (for now)
    if(!is.null(lavoptions$effect.coding) &&
       nchar(lavoptions$effect.coding[1L]) > 0L) {
        warning("lavaan WARNING: automatic bounds not available (yet) if effect.coding is used")
       return(partable)
    }

    # check optim.bounds
    if(is.null(lavoptions$optim.bounds)) {
        # <0.6-6 version
        return(partable)
    } else {
        optim.bounds <- lavoptions$optim.bounds

        # check the elements
        if(is.null(optim.bounds$lower)) {
             optim.bounds$lower <- character(0L)
        } else {
             optim.bounds$lower <- as.character(optim.bounds$lower)
        }
        if(is.null(optim.bounds$upper)) {
             optim.bounds$upper <- character(0L)
        } else {
             optim.bounds$upper <- as.character(optim.bounds$upper)
        }

        if(is.null(optim.bounds$min.reliability.marker)) {
            optim.bounds$min.reliability.marker <- 0.0
        } else {
            if(optim.bounds$min.reliability.marker < 0 ||
               optim.bounds$min.reliability.marker > 1.0) {
                stop("lavaan ERROR: optim.bounds$min.reliability.marker ",
                     "is out of range: ", optim.bounds$min.reliability.marker)
            }
        }

        if(is.null(optim.bounds$factor)) {
            optim.bounds$factor <- 1.0
        }
    }

    # shortcut
    REL <- optim.bounds$min.reliability.marker

    # nothing to do
    if(length(optim.bounds$lower) == 0L &&
       length(optim.bounds$upper) == 0L) {
        return(partable)
    } else {
        # we compute ALL bounds, then we select what we need
        # (otherwise, we can not use the 'factor')

        if(!is.null(partable$lower)) {
            lower.user <- partable$lower
        } else {
            partable$lower <- lower.user <- rep(-Inf, length(partable$lhs))
        }
        if(!is.null(partable$upper)) {
            upper.user <- partable$upper
        } else {
            partable$upper <- upper.user <- rep(+Inf, length(partable$lhs))
        }

        # the 'automatic' bounds
        lower.auto <- rep(-Inf, length(partable$lhs))
        upper.auto <- rep(+Inf, length(partable$lhs))
    }

    # make sure we have lavpta
    if(is.null(lavpta)) {
        lavpta <- lav_partable_attributes(partable)
    }

    # check blocks
    if(is.null(partable$block)) {
        partable$block <- rep(1L, length(partable$lhs))
    }
    block.values <- lav_partable_block_values(partable)

    # check groups
    if(is.null(partable$group)) {
        partable$group <- rep(1L, length(partable$lhs))
    }
    group.values <- lav_partable_group_values(partable)
    ngroups <- length(group.values)

    # compute bounds per group ### TODO: add levels/classes/...
    b <- 0L
    for(g in seq_len(ngroups)) {

        # next block
        b <- b + 1L

        # for this block
        ov.names  <- lavpta$vnames$ov[[b]]
        lv.names  <- lavpta$vnames$lv[[b]]
        lv.marker <- lavpta$vnames$lv.marker[[b]]

        # OV.VAR for this group
        if(lavsamplestats@missing.flag && lavdata@nlevels == 1L) {
            OV.VAR <- diag(lavsamplestats@missing.h1[[g]]$sigma)
        } else {
            if(lavoptions$conditional.x) {
                OV.VAR <- diag(lavsamplestats@res.cov[[g]])
            } else {
                OV.VAR <- diag(lavsamplestats@cov[[g]])
            }
        }


        # we 'process' the parameters per 'type', so we can choose
        # to apply (or not) upper/lower bounds for each type separately

        ################################
        ## 1. (residual) ov variances ##
        ################################
        par.idx <- which(partable$group == group.values[g] &
                         partable$op == "~~" &
                         partable$lhs %in% ov.names &
                         partable$lhs == partable$rhs)

        if(length(par.idx) > 0L) {
            # lower == 0
            lower.auto[par.idx ] <- 0

            # upper == var(ov)
            var.idx <- match(partable$lhs[par.idx], ov.names)
            upper.auto[par.idx] <- OV.VAR[var.idx]

            # enlarge?
            bound.range <- upper.auto[par.idx] - lower.auto[par.idx]
            if( is.finite(optim.bounds$factor) ) {
                new.range <- bound.range * optim.bounds$factor
                diff <- abs(new.range - bound.range)
                lower.auto[par.idx] <- lower.auto[par.idx] - diff
                upper.auto[par.idx] <- upper.auto[par.idx] + diff
            }

            # requested?
            if("ov.var" %in% optim.bounds$lower) {
                partable$lower[par.idx] <- lower.auto[par.idx]
            }
            if("ov.var" %in% optim.bounds$upper) {
                partable$upper[par.idx] <- upper.auto[par.idx]
            }
        } # (res) ov variances

        ################################
        ## 2. (residual) lv variances ##
        ################################
        par.idx <- which(partable$group == group.values[g] &
                         partable$op == "~~" &
                         partable$lhs %in% lv.names &
                         partable$lhs == partable$rhs)

        if(length(par.idx) > 0L) {

            if(lavoptions$std.lv) {

                # lower == 0
                lower.auto[par.idx] <- 0.0

                # upper == 1.0
                upper.auto[par.idx] <- 1.0

            } else {

                for(lv in 1:length(par.idx)) {

                    this.idx <- par.idx[lv]
                    this.lv.name <- partable$lhs[ this.idx ]
                    this.lv.marker <- lv.marker[ match(this.lv.name, lv.names) ]

                    if(nchar(this.lv.marker) > 0L &&
                       this.lv.marker %in% ov.names) {
                        marker.var <- OV.VAR[ match(this.lv.marker, ov.names) ]

                        if(this.lv.name %in% lavpta$vnames$eqs.y[[b]]) {
                            # endogenous LV
                            lower.auto[this.idx] <- 0
                            upper.auto[this.idx] <- marker.var
                        } else {
                            # exogenous LV
                            LOWER <- marker.var - (1 - REL)*marker.var
                            lower.auto[this.idx] <- LOWER
                            upper.auto[this.idx] <- marker.var
                        }
                    } else {
                        # just in case
                        lower.auto[this.idx] <- 0.0
                        upper.auto[this.idx] <- max(OV.VAR) # check!!
                    }
                } # lv
            } # ULI

            # enlarge?
            bound.range <- upper.auto[par.idx] - lower.auto[par.idx]
            if( is.finite(optim.bounds$factor) ) {
                new.range <- bound.range * optim.bounds$factor
                diff <- abs(new.range - bound.range)
                lower.auto[par.idx] <- lower.auto[par.idx] - diff
                upper.auto[par.idx] <- upper.auto[par.idx] + diff
            }

            # requested?
            if("lv.var" %in% optim.bounds$lower) {
                partable$lower[par.idx] <- lower.auto[par.idx]
            }
            if("lv.var" %in% optim.bounds$upper) {
                partable$upper[par.idx] <- upper.auto[par.idx]
            }
        } # lv variances


        #############################################
        ## 3. factor loadings (ov indicators only) ##
        #############################################
        ov.ind.names <- lavpta$vnames$ov.ind[[b]]
        par.idx <- which(partable$group == group.values[g] &
                         partable$op == "=~" &
                         partable$lhs %in% lv.names &
                         partable$rhs %in% ov.ind.names)

        if(length(par.idx) > 0L) {
            if(lavoptions$std.lv) {

                var.idx <- match(partable$rhs[par.idx], ov.names)
                # lower == -1 * sqrt(var(ov))
                lower.auto[par.idx] <- -1 * sqrt(OV.VAR[var.idx])

                # upper == 1.0
                upper.auto[par.idx] <- +1 * sqrt(OV.VAR[var.idx])

            } else {

                lv.Markers <- lv.marker[partable$lhs[par.idx]]
                var.s1 <- OV.VAR[match(lv.marker[partable$lhs[par.idx]],
                                       ov.names)]
                tmp <- var.s1 - (1 - REL)*var.s1
                tmp[is.na(tmp)] <- 0 # just in case...

                var.all <- OV.VAR[ match(partable$rhs[par.idx], ov.names) ]
                lower.auto[par.idx] <- -1 * sqrt(var.all/tmp) # -Inf if REL==0
                upper.auto[par.idx] <- +1 * sqrt(var.all/tmp) # +Inf if REL==0
            }

            # enlarge?
            bound.range <- upper.auto[par.idx] - lower.auto[par.idx]
            if( is.finite(optim.bounds$factor) ) {
                new.range <- bound.range * optim.bounds$factor
                ok.idx <- is.finite(new.range)
                if(length(ok.idx) > 0L) {
                    diff <- abs(new.range[ok.idx] - bound.range[ok.idx])
                    lower.auto[par.idx][ok.idx] <-
                    lower.auto[par.idx][ok.idx] - diff
                    upper.auto[par.idx][ok.idx] <-
                    upper.auto[par.idx][ok.idx] + diff
                }
            }

            # requested?
            if("loadings" %in% optim.bounds$lower) {
                partable$lower[par.idx] <- lower.auto[par.idx]
            }
            if("loadings" %in% optim.bounds$upper) {
                partable$upper[par.idx] <- upper.auto[par.idx]
            }
        } # lambda


    } # g

    # overwrite with lower.user (except -Inf)
    not.inf.idx <- which(lower.user > -Inf)
    if(length(not.inf.idx) > 0L) {
        partable$lower[not.inf.idx] <- lower.user[not.inf.idx]
    }

    # overwrite with upper.user (except +Inf)
    not.inf.idx <- which(upper.user < +Inf)
    if(length(not.inf.idx) > 0L) {
        partable$upper[not.inf.idx] <- upper.user[not.inf.idx]
    }

    # non-free
    non.free.idx <- which(partable$free == 0L)
    if(length(non.free.idx) > 0L && !is.null(partable$start)) {
        partable$lower[non.free.idx] <- partable$start[non.free.idx]
        partable$upper[non.free.idx] <- partable$start[non.free.idx]
    }

    partable
}
