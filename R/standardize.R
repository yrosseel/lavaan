standardize.est.lv <- function(object, partable=NULL, est=NULL,
                               cov.std = TRUE) {

    if(is.null(partable)) partable <- object@ParTable
    if(is.null(est))   est <- object@Fit@est

    out <- est; N <- length(est)
    stopifnot(N == length(partable$lhs))

    GLIST <- object@Model@GLIST
    nmat <- object@Model@nmat
    
    for(g in 1:object@Sample@ngroups) {

        ov.names <- vnames(object@ParTable, "ov", group=g) # not user, 
                                                       # which may be incomplete
        lv.names <- vnames(object@ParTable, "lv", group=g)

        # which mm belong to group g?
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST <- GLIST[ mm.in.group ]

        if(object@Model@representation == "LISREL") {
            # replace negative variances by 0
            idx <- which(diag( MLIST$theta) < 0); diag(MLIST$theta)[idx] <- 0
            # bug reported by Daniel Oberski:
            # als er een negatieve psi is klopt de gestandaardiseerde waarde 
            # niet omdat je hem op 0 zet. Dat hoeft echter niet want een i
            # negatieve psi impliceert niet per se een negatieve variantie i
            # van de latente variabele. 
            #idx <- which(diag( MLIST$psi  ) < 0); diag(MLIST$psi)[idx] <- 0

            if("beta" %in% names(MLIST)) {
                tmp <- -1.0 * MLIST$beta; diag(tmp) <- 1.0
                ibeta.inv <- solve(tmp)
                LV <- (ibeta.inv %*% MLIST$psi %*% t(ibeta.inv))
            } else {
                LV <- MLIST$psi
            }
        }

        ETA2 <- diag(LV)
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

        # 5a ":="
        idx <- which(partable$op == ":=" & partable$group == g)
        if(length(idx) > 0L) {
            x <- out[ partable$free & !duplicated(partable$free) ]
            out[idx] <- object@Model@def.function(x)
        }

        # 5b "=="
        idx <- which(partable$op == "==" & partable$group == g)
        if(length(idx) > 0L) {
            x <- out[ partable$free & !duplicated(partable$free) ]
            out[idx] <- object@Model@ceq.function(x)
        }

        # 5c. "<" or ">"
        idx <- which((partable$op == "<" | partable$op == ">") & partable$group == g)
        if(length(idx) > 0L) {
            x <- out[ partable$free & !duplicated(partable$free) ]
            out[idx] <- object@Model@cin.function(x)
        }

        
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

    Sigma.hat <- object@Fit@Sigma.hat

    for(g in 1:object@Sample@ngroups) {

        ov.names <- vnames(object@ParTable, "ov", group=g) # not user
        lv.names <- vnames(object@ParTable, "lv", group=g)

        OV2 <- diag(Sigma.hat[[g]])
        OV  <- sqrt(OV2)

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

        # 5a ":="
        idx <- which(partable$op == ":=" & partable$group == g)
        if(length(idx) > 0L) {
            x <- out[ partable$free & !duplicated(partable$free) ]
            out[idx] <- object@Model@def.function(x)
        }

        # 5b "=="
        idx <- which(partable$op == "==" & partable$group == g)
        if(length(idx) > 0L) {
            x <- out[ partable$free & !duplicated(partable$free) ]
            out[idx] <- object@Model@ceq.function(x)
        }

        # 5c. "<" or ">"
        idx <- which((partable$op == "<" | partable$op == ">") & partable$group == g)
        if(length(idx) > 0L) {
            x <- out[ partable$free & !duplicated(partable$free) ]
            out[idx] <- object@Model@cin.function(x)
        }
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

    Sigma.hat <- object@Fit@Sigma.hat

    for(g in 1:object@Sample@ngroups) {

        ov.names     <- vnames(object@ParTable, "ov",     group=g) # not user
        ov.names.x   <- vnames(object@ParTable, "ov.x",   group=g)
        ov.names.nox <- vnames(object@ParTable, "ov.nox", group=g)
        lv.names     <- vnames(object@ParTable, "lv",     group=g)

        OV2 <- diag(Sigma.hat[[g]])
        OV  <- sqrt(OV2)

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

        # 5a ":="
        idx <- which(partable$op == ":=" & partable$group == g)
        if(length(idx) > 0L) {
            x <- out[ partable$free & !duplicated(partable$free) ]
            out[idx] <- object@Model@def.function(x)
        }

        # 5b "=="
        idx <- which(partable$op == "==" & partable$group == g)
        if(length(idx) > 0L) {
            x <- out[ partable$free & !duplicated(partable$free) ]
            out[idx] <- object@Model@ceq.function(x)
        }

        # 5c. "<" or ">"
        idx <- which((partable$op == "<" | partable$op == ">") & partable$group == g)
        if(length(idx) > 0L) {
            x <- out[ partable$free & !duplicated(partable$free) ]
            out[idx] <- object@Model@cin.function(x)
        }
    }

    out
}
