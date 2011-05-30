standardize.est.lv <- function(object, user=NULL, est=NULL,
                               cov.std = TRUE) {

    if(is.null(user)) user <- object@User
    if(is.null(est))   est <- object@Fit@est

    out <- est; N <- length(est)
    stopifnot(N == length(user$lhs))

    ov.names <- vnames(object@User, "ov") # not user, which may be incomplete
    lv.names <- vnames(object@User, "lv")
    GLIST <- object@Model@GLIST
    nmat <- object@Model@nmat
    
    for(g in 1:object@Sample@ngroups) {

        # which mm belong to group g?
        mm.in.group <- nmat * (g - 1L) + 1:nmat
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
        idx <- which(user$op == "=~" & !(user$rhs %in% lv.names) & 
                     user$group == g)
        out[idx] <- out[idx] * ETA[ match(user$lhs[idx], lv.names) ]

        # 1b. "=~" regular higher-order lv indicators
        idx <- which(user$op == "=~" & !(user$rhs %in% ov.names) &
                     user$group == g)
        out[idx] <- ( out[idx] * ETA[ match(user$lhs[idx], lv.names) ]
                               / ETA[ match(user$rhs[idx], lv.names) ] )

        # 1c. "=~" indicators that are both in ov and lv
        #idx <- which(user$op == "=~" & user$rhs %in% ov.names
        #                             & user$rhs %in% lv.names &
        #             user$group == g)

        # 2. "~" regressions
        idx <- which(user$op == "~" & user$lhs %in% lv.names &
                     user$group == g)
        out[idx] <- out[idx] / ETA[ match(user$lhs[idx], lv.names) ] 

        idx <- which(user$op == "~" & user$rhs %in% lv.names &
                     user$group == g)
        out[idx] <- out[idx] * ETA[ match(user$rhs[idx], lv.names) ]

        # 3a. "~~" ov
        #idx <- which(user$op == "~~" & !(user$lhs %in% lv.names) & 
        #             user$group == g)

        # 3b. "~~" lv
        # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances 
        #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
        #            were i.var and j.var where diagonal elements of ETA
        #
        #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
        #            elements are the 'PSI' diagonal elements!!

        # variances
        rv.idx <- which(user$op == "~~" & user$rhs %in% lv.names &
                        user$lhs == user$rhs &
                        user$group == g)
        out[rv.idx] <- ( out[rv.idx] / ETA[ match(user$lhs[rv.idx], lv.names) ]
                                     / ETA[ match(user$rhs[rv.idx], lv.names) ] )

        # covariances
        idx <- which(user$op == "~~" & user$rhs %in% lv.names &
                     user$lhs != user$rhs &
                     user$group == g)
        if(cov.std == FALSE) {
            out[idx] <- ( out[idx] / ETA[ match(user$lhs[idx], lv.names) ]
                                   / ETA[ match(user$rhs[idx], lv.names) ] )
        } else {
            RV   <- sqrt(est[rv.idx])
            rv.names <- user$lhs[rv.idx]
            out[idx] <- ( out[idx] / RV[ match(user$lhs[idx], rv.names) ]
                                   / RV[ match(user$rhs[idx], rv.names) ] )
        }

        # 4a. "~1" ov
        #idx <- which(user$op == "~1" & !(user$lhs %in% lv.names) &
        #             user$group == g)

        # 4b. "~1" lv
        idx <- which(user$op == "~1" & user$lhs %in% lv.names &
                     user$group == g)
        out[idx] <- out[idx] / ETA[ match(user$lhs[idx], lv.names) ]

    }

    out
}

standardize.est.all <- function(object, user=NULL, est=NULL, est.std=NULL,
                                cov.std = TRUE) {

    if(is.null(user)) user <- object@User
    if(is.null(est))   est <- object@Fit@est
    if(is.null(est.std)) {
        est.std <- standardize.est.lv(object, user=user, est=est)
    } 

    out <- est.std; N <- length(est.std)
    stopifnot(N == length(user$lhs))

    ov.names <- vnames(object@User, "ov") # not user, which may be incomplet
    lv.names <- vnames(object@User, "lv")
    GLIST <- object@Model@GLIST
    nmat <- object@Model@nmat

    Sigma.hat <- computeSigmaHat(object@Model)

    for(g in 1:object@Sample@ngroups) {

        OV2 <- diag(Sigma.hat[[g]])
        OV  <- sqrt(OV2)

        # 1a. "=~" regular indicators
        idx <- which(user$op == "=~" & !(user$rhs %in% lv.names) &
                     user$group == g)
        out[idx] <- out[idx] / OV[ match(user$rhs[idx], ov.names) ]

        # 1b. "=~" regular higher-order lv indicators

        # 1c. "=~" indicators that are both in ov and lv
        #idx <- which(user$op == "=~" & user$rhs %in% ov.names
        #                             & user$rhs %in% lv.names &
        #             user$group == g)

        # 2. "~" regressions
        idx <- which(user$op == "~" & user$lhs %in% ov.names &
                     user$group == g)
        out[idx] <- out[idx] / OV[ match(user$lhs[idx], ov.names) ]

        idx <- which(user$op == "~" & user$rhs %in% ov.names &
                     user$group == g)
        out[idx] <- out[idx] * OV[ match(user$rhs[idx], ov.names) ]

        # 3a. "~~" ov
        # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances 
        #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
        #            were i.var and j.var where diagonal elements of OV
        #
        #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
        #            elements are the 'THETA' diagonal elements!!

        # variances
        rv.idx <- which(user$op == "~~" & !(user$lhs %in% lv.names) & 
                        user$lhs == user$rhs &
                        user$group == g)
        out[rv.idx] <- ( out[rv.idx] / OV[ match(user$lhs[rv.idx], ov.names) ]
                                     / OV[ match(user$rhs[rv.idx], ov.names) ] )

        # covariances
        idx <- which(user$op == "~~" & !(user$lhs %in% lv.names) &
                     user$lhs != user$rhs &
                     user$group == g)
        if(cov.std == FALSE) {
            out[idx] <- ( out[idx] / OV[ match(user$lhs[idx], ov.names) ]
                                   / OV[ match(user$rhs[idx], ov.names) ] )
        } else {
            RV   <- sqrt(est[rv.idx])
            rv.names <- user$lhs[rv.idx]
            out[idx] <- ( out[idx] / RV[ match(user$lhs[idx], rv.names) ]
                                   / RV[ match(user$rhs[idx], rv.names) ] )
        }

        # 3b. "~~" lv
        #idx <- which(user$op == "~~" & user$rhs %in% lv.names &
        #             user$group == g)

        # 4a. "~1" ov
        idx <- which(user$op == "~1" & !(user$lhs %in% lv.names) &
                     user$group == g)
        out[idx] <- out[idx] / OV[ match(user$lhs[idx], ov.names) ]

        # 4b. "~1" lv
        #idx <- which(user$op == "~1" & user$lhs %in% lv.names &
        #             user$group == g)

    }

    out
}
