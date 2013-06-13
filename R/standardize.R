standardize.est.lv.x <- function(x, object, partable=NULL, cov.std=TRUE) {
    standardize.est.lv(object=object, partable=partable, est=x, cov.std=cov.std)
}

standardize.est.all.x <- function(x, object, partable=NULL, cov.std=TRUE) {
    standardize.est.all(object=object, partable=partable, est=x, est.std=NULL,
                        cov.std=cov.std)
}

standardize.est.all.nox.x <- function(x, object, partable=NULL, cov.std=TRUE) {
    standardize.est.all.nox(object=object, partable=partable, est=x, 
                            est.std=NULL, cov.std=cov.std)
}

unstandardize.est.ov.x <- function(x, object) {
    partable <- object@ParTable
    partable$ustart <- x
    unstandardize.est.ov(partable=partable, ov.var=object@SampleStats@var, 
                         cov.std=TRUE)
}

standardize.est.lv <- function(object, partable=NULL, est=NULL,
                               cov.std = TRUE) {

    if(is.null(partable)) partable <- object@ParTable
    if(is.null(est))   est <- object@Fit@est

    out <- est; N <- length(est)
    stopifnot(N == length(partable$lhs))

    GLIST <- object@Model@GLIST
    nmat <- object@Model@nmat

    # compute ETA
    LV.ETA <- computeVETA(object@Model, samplestats=object@SampleStats)
    
    for(g in 1:object@Data@ngroups) {

        ov.names <- vnames(object@ParTable, "ov", group=g) # not user, 
                                                       # which may be incomplete
        lv.names <- vnames(object@ParTable, "lv", group=g)
       
        # shortcut: no latents in this group, nothing to do
        if(length(lv.names) == 0L)
            next

        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST <- GLIST[ mm.in.group ]

        ETA2 <- diag(LV.ETA[[g]])
        ETA  <- sqrt(ETA2)

        # 1a. "=~" regular indicators
        idx <- which(partable$op == "=~" & !(partable$rhs %in% lv.names) & 
                     partable$group == g)
        out[idx] <- out[idx] * ETA[ match(partable$lhs[idx], lv.names) ]

        # 1b. "=~" regular higher-order lv indicators
        idx <- which(partable$op == "=~" & !(partable$rhs %in% ov.names) &
                     partable$group == g)
        out[idx] <- ( out[idx] * ETA[ match(partable$lhs[idx], lv.names) ]
                               / ETA[ match(partable$rhs[idx], lv.names) ] )

        # 1c. "=~" indicators that are both in ov and lv
        #idx <- which(partable$op == "=~" & partable$rhs %in% ov.names
        #                             & partable$rhs %in% lv.names &
        #             partable$group == g)

        # 2. "~" regressions (and "<~")
        idx <- which((partable$op == "~" | partable$op == "<~") & 
                     partable$lhs %in% lv.names &
                     partable$group == g)
        out[idx] <- out[idx] / ETA[ match(partable$lhs[idx], lv.names) ] 

        idx <- which((partable$op == "~" | partable$op == "<~") & 
                     partable$rhs %in% lv.names &
                     partable$group == g)
        out[idx] <- out[idx] * ETA[ match(partable$rhs[idx], lv.names) ]

        # 3a. "~~" ov
        #idx <- which(partable$op == "~~" & !(partable$lhs %in% lv.names) & 
        #             partable$group == g)

        # 3b. "~~" lv
        # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances 
        #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
        #            were i.var and j.var where diagonal elements of ETA
        #
        #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
        #            elements are the 'PSI' diagonal elements!!

        # variances
        rv.idx <- which(partable$op == "~~" & partable$rhs %in% lv.names &
                        partable$lhs == partable$rhs &
                        partable$group == g)
        out[rv.idx] <- ( out[rv.idx] / ETA[ match(partable$lhs[rv.idx], lv.names) ]
                                     / ETA[ match(partable$rhs[rv.idx], lv.names) ] )

        # covariances
        idx <- which(partable$op == "~~" & partable$rhs %in% lv.names &
                     partable$lhs != partable$rhs &
                     partable$group == g)
        if(cov.std == FALSE) {
            out[idx] <- ( out[idx] / ETA[ match(partable$lhs[idx], lv.names) ]
                                   / ETA[ match(partable$rhs[idx], lv.names) ] )
        } else {
            RV   <- sqrt(est[rv.idx])
            rv.names <- partable$lhs[rv.idx]
            out[idx] <- ( out[idx] / RV[ match(partable$lhs[idx], rv.names) ]
                                   / RV[ match(partable$rhs[idx], rv.names) ] )
        }

        # 4a. "~1" ov
        #idx <- which(partable$op == "~1" & !(partable$lhs %in% lv.names) &
        #             partable$group == g)

        # 4b. "~1" lv
        idx <- which(partable$op == "~1" & partable$lhs %in% lv.names &
                     partable$group == g)
        out[idx] <- out[idx] / ETA[ match(partable$lhs[idx], lv.names) ]
    }

    # 5a ":="
    idx <- which(partable$op == ":=")
    if(length(idx) > 0L) {
        x <- out[ partable$free & !duplicated(partable$free) ]
        out[idx] <- object@Model@def.function(x)
    }

    # 5b "=="
    idx <- which(partable$op == "==")
    if(length(idx) > 0L) {
        x <- out[ partable$free & !duplicated(partable$free) ]
        out[idx] <- object@Model@ceq.function(x)
    }

    # 5c. "<" or ">"
    idx <- which((partable$op == "<" | partable$op == ">"))
    if(length(idx) > 0L) {
        x <- out[ partable$free & !duplicated(partable$free) ]
        out[idx] <- object@Model@cin.function(x)
    }

    out
}

standardize.est.all <- function(object, partable=NULL, est=NULL, est.std=NULL,
                                cov.std = TRUE) {

    if(is.null(partable)) partable <- object@ParTable
    if(is.null(est))   est <- object@Fit@est
    if(is.null(est.std)) {
        est.std <- standardize.est.lv(object, partable=partable, est=est)
    } 

    out <- est.std; N <- length(est.std)
    stopifnot(N == length(partable$lhs))

    GLIST <- object@Model@GLIST

    VY <- computeVY(object@Model, samplestats=object@SampleStats)

    for(g in 1:object@Data@ngroups) {

        ov.names <- vnames(object@ParTable, "ov", group=g) # not user
        lv.names <- vnames(object@ParTable, "lv", group=g)

        OV  <- sqrt(VY[[g]])

        if(object@Model@categorical) {
            # extend OV with ov.names.x
            ov.names.x <- vnames(object@ParTable, "ov.x", group=g)
            ov.names <- c(ov.names, ov.names.x)
            OV <- c(OV, sqrt(diag(object@SampleStats@cov.x[[g]])))
        }

        # 1a. "=~" regular indicators
        idx <- which(partable$op == "=~" & !(partable$rhs %in% lv.names) &
                     partable$group == g)
        out[idx] <- out[idx] / OV[ match(partable$rhs[idx], ov.names) ]

        # 1b. "=~" regular higher-order lv indicators

        # 1c. "=~" indicators that are both in ov and lv
        #idx <- which(partable$op == "=~" & partable$rhs %in% ov.names
        #                             & partable$rhs %in% lv.names &
        #             partable$group == g)

        # 2. "~" regressions (and "<~")
        idx <- which((partable$op == "~" | partable$op == "<~") & 
                     partable$lhs %in% ov.names &
                     partable$group == g)
        out[idx] <- out[idx] / OV[ match(partable$lhs[idx], ov.names) ]

        idx <- which((partable$op == "~" | partable$op == "<~") & 
                     partable$rhs %in% ov.names &
                     partable$group == g)
        out[idx] <- out[idx] * OV[ match(partable$rhs[idx], ov.names) ]

        # 3a. "~~" ov
        # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances 
        #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
        #            were i.var and j.var where diagonal elements of OV
        #
        #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
        #            elements are the 'THETA' diagonal elements!!

        # variances
        rv.idx <- which(partable$op == "~~" & !(partable$lhs %in% lv.names) & 
                        partable$lhs == partable$rhs &
                        partable$group == g)
        out[rv.idx] <- ( out[rv.idx] / OV[ match(partable$lhs[rv.idx], ov.names) ]
                                     / OV[ match(partable$rhs[rv.idx], ov.names) ] )

        # covariances
        idx <- which(partable$op == "~~" & !(partable$lhs %in% lv.names) &
                     partable$lhs != partable$rhs &
                     partable$group == g)
        if(length(idx) > 0L) {
            if(cov.std == FALSE) {
                out[idx] <- ( out[idx] / OV[ match(partable$lhs[idx], ov.names) ]
                                       / OV[ match(partable$rhs[idx], ov.names) ] )
            } else {
                RV   <- sqrt(est[rv.idx])
                rv.names <- partable$lhs[rv.idx]
                out[idx] <- ( out[idx] / RV[ match(partable$lhs[idx], rv.names) ]
                                       / RV[ match(partable$rhs[idx], rv.names) ] )
            }
        }

        # 3b. "~~" lv
        #idx <- which(partable$op == "~~" & partable$rhs %in% lv.names &
        #             partable$group == g)

        # 4a. "~1" ov
        idx <- which(partable$op == "~1" & !(partable$lhs %in% lv.names) &
                     partable$group == g)
        out[idx] <- out[idx] / OV[ match(partable$lhs[idx], ov.names) ]

        # 4b. "~1" lv
        #idx <- which(partable$op == "~1" & partable$lhs %in% lv.names &
        #             partable$group == g)

        # 4c. "|" thresholds
        idx <- which(partable$op == "|" & !(partable$lhs %in% lv.names) &
                     partable$group == g)
        out[idx] <- out[idx] / OV[ match(partable$lhs[idx], ov.names) ]

        # 4d. "~*~" scales
        idx <- which(partable$op == "~*~" & !(partable$lhs %in% lv.names) &
                     partable$group == g)
        out[idx] <- 1.0
    }

    # 5a ":="
    idx <- which(partable$op == ":=")
    if(length(idx) > 0L) {
        x <- out[ partable$free & !duplicated(partable$free) ]
        out[idx] <- object@Model@def.function(x)
    }

    # 5b "=="
    idx <- which(partable$op == "==")
    if(length(idx) > 0L) {
        x <- out[ partable$free & !duplicated(partable$free) ]
        out[idx] <- object@Model@ceq.function(x)
    }

    # 5c. "<" or ">"
    idx <- which((partable$op == "<" | partable$op == ">"))
    if(length(idx) > 0L) {
        x <- out[ partable$free & !duplicated(partable$free) ]
        out[idx] <- object@Model@cin.function(x)
    }

    out
}


standardize.est.all.nox <- function(object, partable=NULL, est=NULL, 
                                    est.std=NULL, cov.std = TRUE) {

    if(is.null(partable)) partable <- object@ParTable
    if(is.null(est))   est <- object@Fit@est
    if(is.null(est.std)) {
        est.std <- standardize.est.lv(object, partable=partable, est=est)
    } 

    out <- est.std; N <- length(est.std)
    stopifnot(N == length(partable$lhs))

    GLIST <- object@Model@GLIST

    VY <- computeVY(object@Model, samplestats=object@SampleStats)

    for(g in 1:object@Data@ngroups) {

        ov.names     <- vnames(object@ParTable, "ov",     group=g) # not user
        ov.names.x   <- vnames(object@ParTable, "ov.x",   group=g)
        ov.names.nox <- vnames(object@ParTable, "ov.nox", group=g)
        lv.names     <- vnames(object@ParTable, "lv",     group=g)

        OV  <- sqrt(VY[[g]])

        if(object@Model@categorical) {
            # extend OV with ov.names.x
            ov.names.x <- vnames(object@ParTable, "ov.x", group=g)
            ov.names <- c(ov.names, ov.names.x)
            OV <- c(OV, sqrt(diag(object@SampleStats@cov.x[[g]])))
        }

        # 1a. "=~" regular indicators
        idx <- which(partable$op == "=~" & !(partable$rhs %in% lv.names) &
                     partable$group == g)
        out[idx] <- out[idx] / OV[ match(partable$rhs[idx], ov.names) ]

        # 1b. "=~" regular higher-order lv indicators

        # 1c. "=~" indicators that are both in ov and lv
        #idx <- which(partable$op == "=~" & partable$rhs %in% ov.names
        #                             & partable$rhs %in% lv.names &
        #             partable$group == g)

        # 2. "~" regressions (and "<~")
        idx <- which((partable$op == "~" | partable$op == "<~") & 
                     partable$lhs %in% ov.names &
                     partable$group == g)
        out[idx] <- out[idx] / OV[ match(partable$lhs[idx], ov.names) ]

        idx <- which((partable$op == "~" | partable$op == "<~") & 
                     partable$rhs %in% ov.names.nox &
                     partable$group == g)
        out[idx] <- out[idx] * OV[ match(partable$rhs[idx], ov.names.nox) ]

        # 3a. "~~" ov
        # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances 
        #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
        #            were i.var and j.var where diagonal elements of OV
        #
        #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
        #            elements are the 'THETA' diagonal elements!!

        # variances
        rv.idx <- which(partable$op == "~~" & !(partable$lhs %in% lv.names) & 
                        !(partable$lhs %in% ov.names.x) &
                        partable$lhs == partable$rhs &
                        partable$group == g)
        out[rv.idx] <- ( out[rv.idx] / OV[ match(partable$lhs[rv.idx], ov.names) ]
                                     / OV[ match(partable$rhs[rv.idx], ov.names) ] )

        # covariances
        idx <- which(partable$op == "~~" & !(partable$lhs %in% lv.names) &
                     !(partable$lhs %in% ov.names.x) &
                     !(partable$rhs %in% ov.names.x) &
                     partable$lhs != partable$rhs &
                     partable$group == g)
        if(length(idx) > 0L) {
            if(cov.std == FALSE) {
                out[idx] <- ( out[idx] / OV[ match(partable$lhs[idx], ov.names) ]
                                       / OV[ match(partable$rhs[idx], ov.names) ] )
            } else {
                RV   <- sqrt(est[rv.idx])
                rv.names <- partable$lhs[rv.idx]
                out[idx] <- ( out[idx] / RV[ match(partable$lhs[idx], rv.names) ]
                                       / RV[ match(partable$rhs[idx], rv.names) ] )
            }
        }

        # 3b. "~~" lv
        #idx <- which(partable$op == "~~" & partable$rhs %in% lv.names &
        #             partable$group == g)

        # 4a. "~1" ov
        idx <- which(partable$op == "~1" & !(partable$lhs %in% lv.names) &
                     !(partable$lhs %in% ov.names.x) &
                     partable$group == g)
        out[idx] <- out[idx] / OV[ match(partable$lhs[idx], ov.names) ]

        # 4b. "~1" lv
        #idx <- which(partable$op == "~1" & partable$lhs %in% lv.names &
        #             partable$group == g)

        # 4c. "|" thresholds
        idx <- which(partable$op == "|" & !(partable$lhs %in% lv.names) &
                     partable$group == g)
        out[idx] <- out[idx] / OV[ match(partable$lhs[idx], ov.names) ]

        # 4d. "~*~" scales
        idx <- which(partable$op == "~*~" & !(partable$lhs %in% lv.names) &
                     partable$group == g)
        out[idx] <- 1.0
    }

    # 5a ":="
    idx <- which(partable$op == ":=")
    if(length(idx) > 0L) {
        x <- out[ partable$free & !duplicated(partable$free) ]
        out[idx] <- object@Model@def.function(x)
    }

    # 5b "=="
    idx <- which(partable$op == "==")
    if(length(idx) > 0L) {
        x <- out[ partable$free & !duplicated(partable$free) ]
        out[idx] <- object@Model@ceq.function(x)
    }

    # 5c. "<" or ">"
    idx <- which((partable$op == "<" | partable$op == ">"))
    if(length(idx) > 0L) {
        x <- out[ partable$free & !duplicated(partable$free) ]
        out[idx] <- object@Model@cin.function(x)
    }

    out
}

unstandardize.est.ov <- function(partable, ov.var=NULL, cov.std=TRUE) {

    # check if ustart is missing; if so, look for est
    if(is.null(partable$ustart)) 
        partable$ustart <- partable$est
  
    # check if group is missing
    if(is.null(partable$group)) 
        partable$group <- rep(1L, length(partable$ustart))

    stopifnot(!any(is.na(partable$ustart)))
    est <- out <- partable$ustart
    N <- length(est)
    ngroups <- max(partable$group)

    # if ov.var is NOT a list, make a list
    if(!is.list(ov.var)) {
        tmp <- ov.var
        ov.var <- vector("list", length=ngroups)
        ov.var[1:ngroups] <- list(tmp)
    }

    for(g in 1:ngroups) {

        ov.names <- vnames(partable, "ov", group=g) # not user
        lv.names <- vnames(partable, "lv", group=g)

        OV  <- sqrt(ov.var[[g]])

        # 1a. "=~" regular indicators
        idx <- which(partable$op == "=~" & !(partable$rhs %in% lv.names) &
                     partable$group == g)
        out[idx] <- out[idx] * OV[ match(partable$rhs[idx], ov.names) ]

        # 1b. "=~" regular higher-order lv indicators

        # 1c. "=~" indicators that are both in ov and lv
        #idx <- which(partable$op == "=~" & partable$rhs %in% ov.names
        #                             & partable$rhs %in% lv.names &
        #             partable$group == g)

        # 2. "~" regressions (and "<~")
        idx <- which((partable$op == "~" | partable$op == "<~") & 
                     partable$lhs %in% ov.names &
                     partable$group == g)
        out[idx] <- out[idx] * OV[ match(partable$lhs[idx], ov.names) ]

        idx <- which((partable$op == "~" | partable$op == "<~") & 
                     partable$rhs %in% ov.names &
                     partable$group == g)
        out[idx] <- out[idx] / OV[ match(partable$rhs[idx], ov.names) ]

        # 3a. "~~" ov
        # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances 
        #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
        #            were i.var and j.var where diagonal elements of OV
        #
        #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
        #            elements are the 'THETA' diagonal elements!!

        # variances
        rv.idx <- which(partable$op == "~~" & !(partable$lhs %in% lv.names) & 
                        partable$lhs == partable$rhs &
                        partable$group == g)
        out[rv.idx] <- ( out[rv.idx] * OV[ match(partable$lhs[rv.idx], ov.names) ]
                                     * OV[ match(partable$rhs[rv.idx], ov.names) ] )

        # covariances
        idx <- which(partable$op == "~~" & !(partable$lhs %in% lv.names) &
                     partable$lhs != partable$rhs &
                     partable$group == g)
        if(length(idx) > 0L) {
            if(cov.std == FALSE) {
                out[idx] <- ( out[idx] * OV[ match(partable$lhs[idx], ov.names) ]
                                       * OV[ match(partable$rhs[idx], ov.names) ] )
            } else {
                # RV   <- sqrt(est[rv.idx])
                RV   <- sqrt(out[rv.idx])
                rv.names <- partable$lhs[rv.idx]
                out[idx] <- ( out[idx] * RV[ match(partable$lhs[idx], rv.names) ]
                                       * RV[ match(partable$rhs[idx], rv.names) ] )
            }
        }

        # 3b. "~~" lv
        #idx <- which(partable$op == "~~" & partable$rhs %in% lv.names &
        #             partable$group == g)

        # 4a. "~1" ov
        idx <- which(partable$op == "~1" & !(partable$lhs %in% lv.names) &
                     partable$group == g)
        out[idx] <- out[idx] * OV[ match(partable$lhs[idx], ov.names) ]

        # 4b. "~1" lv
        #idx <- which(partable$op == "~1" & partable$lhs %in% lv.names &
        #             partable$group == g)

    }

    # 5a ":="
    # 5b "=="
    # 5c. "<" or ">"

    out
}

