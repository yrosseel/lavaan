# univariate modification indices
#

modindices <- function(object, 
                       standardized = TRUE, 

                       # power statistics?
                       power = FALSE, 
                       delta = 0.1, 
                       alpha = 0.05, 
                       high.power = 0.75,

                       # customize output
                       sort. = FALSE, 
                       minimum.value = 0.0, 
                       maximum.number = nrow(LIST),
                       na.remove = FALSE, 
                       op = NULL) {

    # check if model has converged
    if(object@Fit@npar > 0L && !object@Fit@converged) {
        warning("lavaan WARNINGS: model did not converge")
    }

    # sanity check
    if(power) standardized <- TRUE

    # create extended parameter LIST, including model parameters
    partable <- object@ParTable
    LIST <- lav_partable_merge(lav_partable_full(object@ParTable, free = TRUE),
                               partable[c("lhs","op","rhs","group","free")], 
                               remove = TRUE, warn = FALSE)

    # add matrix representation
    if(object@Model@representation == "LISREL") {
        REP <- representation.LISREL(partable = object@ParTable, target = LIST,
                                     extra = FALSE)
    } else {
        stop("only LISREL representation has been implemented")
    }
    LIST <- cbind(LIST, as.data.frame(REP, stringsAsFactors = FALSE))

    # here we remove `non-existing' parameters (depends on the matrix
    # representation (eg in LISREL rep, there is no ~~ between lv and ov)
    idx <- which( nchar(LIST$mat) > 0L &
                  !is.na(LIST$row) & LIST$row > 0L &
                  !is.na(LIST$col) & LIST$col > 0L )
    LIST <- LIST[idx,]

    # check for duplicated matrix elements
    idx <- which( duplicated(LIST[,c("group", "mat","row","col")], 
                  fromLast=TRUE) )
    # FIXME!!! what to do here?
    # remove....
    if(length(idx) > 0L) LIST <- LIST[-idx,]

    # here we could/should remove elements that will produce NA anyways...
    # eg - first indicator of factors
    #    - regressions that are already free covariances
    # TODO?? (now, we just compute them, get NA, and remove them)

    # master index: *all* elements that we will feed to computeDelta
    LIST$id <- 1:nrow(LIST)

    # compute Delta for remaining parameters
    ngroups  <- object@Model@ngroups
    mmNumber <- object@Model@nmat
    mmNames  <- names(object@Model@GLIST)
    m.el.idx <- x.el.idx <- vector("list", length=length(object@Model@GLIST))
    offset <- 0L
    for(g in 1:ngroups) {
        for(mm in 1:mmNumber[g]) {
            # offset in GLIST
            offset <- offset + 1L

            # select elements for this matrix
            idx <- which(LIST$group == g & LIST$mat == mmNames[offset])

            tmp <- matrix(0L, nrow=nrow(object@Model@GLIST[[mm]]),
                              ncol=ncol(object@Model@GLIST[[mm]]))

            # assign id values
            tmp[ cbind(LIST$row[idx], LIST$col[idx]) ] <- LIST$id[idx]
            if(object@Model@isSymmetric[mm]) {
                # everything is in upper tri, but we only extract lower tri
                tmp <- t(tmp)
            }
            m.el.idx[[offset]] <-     which(tmp > 0)
            x.el.idx[[offset]] <- tmp[which(tmp > 0)]
        }
    }
    Delta <- computeDelta(lavmodel = object@Model, 
                          m.el.idx. = m.el.idx, x.el.idx. = x.el.idx)

    # compute information matrix
    E <- computeExpectedInformation(lavmodel       = object@Model, 
                                    lavsamplestats = object@SampleStats,
                                    estimator      = object@Options$estimator,
                                    Delta          = Delta)
    Q <- (1/object@SampleStats@ntotal) * E

    # list!
    DX <- lav_model_gradient(lavmodel       = object@Model, 
                             GLIST          = NULL, 
                             lavsamplestats = object@SampleStats,
                             lavdata        = object@Data,
                             lavcache       = object@Cache,
                             type           = "allofthem", 
                             estimator      = object@Options$estimator,
                             group.weight   = TRUE,
                             constraints    = TRUE, #### FIXME???
                             Delta          = Delta,
                             m.el.idx       = m.el.idx,
                             x.el.idx       = x.el.idx)

    # flatten list DX to a vector dx
    dx <- numeric(0)
    for(mm in 1:length(object@Model@GLIST)) {
        dx[ x.el.idx[[mm]] ] <- DX[[mm]][ m.el.idx[[mm]] ] 
    }
                             
    # Saris, Satorra & Sorbom 1987
    # partition Q into Q_11, Q_22 and Q_12/Q_21
    # which elements of Q correspond with 'free' and 'nonfree' parameters?
    #all.idx      <- LIST$id
    #free.idx     <- LIST$id[LIST$free > 0L & !duplicated(LIST$free)]
    #eq.idx       <- LIST$id[LIST$eq.id > 0L]
    #eq.id        <- LIST$eq.id[ eq.idx ]

    # NOTE: since 0.5-18, we do not make a distinction anymore between
    #       equality constrained parameters (using labels), and general
    #       linear constraints (using ==)
    #       As a result, the 'free' parameters are not necessarily a clean
    #       subset of the total set of parameters
    #
    #       But we can compute what would happen if we release a constraint,
    #       one at a time
    #       FIXME!!!


    free.idx    <- which(LIST$free  > 0L)
    nonfree.idx <- which(LIST$free == 0L)
    eq.idx <- integer(0L)

    # partition Q
    Q11 <- Q[nonfree.idx, nonfree.idx]
    Q12 <- Q[nonfree.idx,    free.idx]
    Q21 <- Q[free.idx,    nonfree.idx]
    Q22 <- Q[free.idx,       free.idx]

    # take care of constraints (if any)
    # handle constraints
    if(nrow(object@Model@con.jac) > 0L) {
        H <- object@Model@con.jac
        inactive.idx <- attr(H, "inactive.idx")
        lambda <- object@Model@con.lambda # lagrangean coefs
        if(length(inactive.idx) > 0L) {
            H <- H[-inactive.idx,,drop=FALSE]
            lambda <- lambda[-inactive.idx]
        }
        # if length(eq.idx) > 0L, remove them from the columns of H
        #if(length(eq.idx) > 0L) {
        #    free.again.idx <- LIST$id[LIST$free > 0L & !duplicated(LIST$free)]
        #    idx <- LIST$free[ free.again.idx[ free.again.idx %in% eq.idx ] ]
        #    H <- H[,-idx]
        #}
        if(nrow(H) > 0L) {
            H0 <- matrix(0,nrow(H),nrow(H))
            H10 <- matrix(0, ncol(Q22), nrow(H))
            DL <- 2*diag(lambda, nrow(H), nrow(H))
            E3 <- rbind( cbind(     Q22,  H10, t(H)),
                         cbind(t(H10),     DL,  H0),
                         cbind(     H,     H0,  H0)  )
            Q22.inv <- MASS::ginv(E3)[1:ncol(Q22), 1:ncol(Q22)]
            # FIXME: better include inactive + slacks??
        } else {
            Q22.inv <- solve(Q22)
        }
    } else {
        Q22.inv <- solve(Q22)
    }

    #print( dim(Q22.inv) )
    #print( Q22.inv[1:5, 1:5] )
    # NOTE: we could perhaps use Q22.inv <- vcov(fit) * N^2?

    V <- Q11 - Q12 %*% Q22.inv %*% Q21
    V.diag <- diag(V)
    # dirty hack: catch very small or negative values in diag(V)
    # this is needed eg when parameters are not identified if freed-up;
    idx <- which(V.diag < 1.0e-15); V.diag[idx] <- as.numeric(NA)

    # create and fill in mi
    mi <- numeric( length(dx) )
    mi[nonfree.idx] <- dx[nonfree.idx]^2 / V.diag

    # take care of equality constraints
    # Sorbom 1989 equations 12 and 13
    # surely this code needs some serious optimization
    # but it seems to work for now
    # we should use the K matrix somewhere...
    if(FALSE) {
        for(i in 1:length(eq.idx)) {
            # index in all.idx (only mi.fixed elements)
            theta2.idx <- eq.idx[i]
            thetac.idx <- eq.idx[eq.id == eq.id[i] & eq.idx != eq.idx[i]]
            #f.idx <- free.idx[-which(free.idx == theta2.idx)]
            f.idx <- free.idx
            Q22 <- Q[f.idx, f.idx, drop=FALSE]

            ivec <- as.matrix(rep(1,length(thetac.idx)))
            C11 <- Q[thetac.idx, thetac.idx, drop=FALSE]
            C21 <- Q[theta2.idx, thetac.idx, drop=FALSE]; C12 <- t(C21)
            C22 <- Q[theta2.idx, theta2.idx, drop=FALSE] 
            D1  <- Q[f.idx, thetac.idx, drop=FALSE]; 
            d2  <- Q[f.idx, theta2.idx, drop=FALSE]

            c11 <- t(ivec) %*% C11 %*% ivec
            c21 <- C21 %*% ivec; c12 <- t(c21)
            c22 <- C22

            C <- rbind( cbind(c11,c12), cbind(c21, c22) )
        
            D <- cbind(D1 %*% ivec, d2)

   
            e12 <- (D1 %*% ivec) + d2; e21 <- t(e12)
            e22 <- (t(ivec) %*% C11 %*% ivec) + (t(ivec) %*% C12) + (C21 %*% ivec) + C22
            E.star <- rbind( cbind(Q22, e12), cbind(e21, e22) )

            E.star.inv <- inv.chol(E.star)
            n.free <- length(f.idx); n.ref <- length(theta2.idx)
            E11 <- E.star.inv[1:n.free, 1:n.free]
            E21 <- t(E.star.inv[n.free + 1:n.ref, 1:n.free]); E12 <- t(E21)
            E22 <- as.numeric(E.star.inv[n.free + 1:n.ref, n.free + 1:n.ref])
            E.inv <- E11 - (E12 %*% E21)/E22

            H <- C - (t(D) %*% E.inv %*% D)
            h11 <- H[1,1]; h12 <- H[1,2]; h21 <- H[2,1]; h22 <- H[2,2]

            mi[theta2.idx] <- dx[theta2.idx]^2 * (h11 + 2*h21 + h22) / (h11*h22 - h12^2)
        }
    }

    # EPC
    d <- (-1 * object@SampleStats@ntotal) * dx
    # needed?
    d[which(abs(d) < 1e-15)] <- 1.0
    epc <- mi/d

    # FIXME: epc for equality constraints must be adapted!?
    # we need a reference for this!!!!
    #if(object@Model@eq.constraints) {
    #    for(i in 1:length(eq.idx)) {
    #        # this is just a temporary solution
    #        # until we get it right
    #        neq <- length(eq.idx[eq.id == eq.id[i]])
    #        epc[eq.idx[i]] <- epc[eq.idx[i]] / neq * (neq - 1)
    #    }
    #}

    LIST$mi <- mi
    if(length(object@Fit@test) > 1L) {
        LIST$mi.scaled <- mi / object@Fit@test[[2]]$scaling.factor
    }
    LIST$epc <- epc

    # remove some rows
    #idx <- which(LIST$free > 0L & !duplicated(LIST$free) & !LIST$eq.id > 0L)
    #LIST <- LIST[-idx,]
    if(length(free.idx) > 0L) {
        LIST <- LIST[-free.idx,]
    }

    # standardize?
    if(standardized) {
        LIST$sepc.lv  <- standardize.est.lv(object, partable=LIST, est=LIST$epc,
                                            cov.std=FALSE)
        LIST$sepc.all <- standardize.est.all(object, partable=LIST, est=LIST$epc,
                                             est.std=LIST$sepc.lv,
                                             cov.std=FALSE)
        LIST$sepc.nox <- standardize.est.all.nox(object,partable=LIST,est=LIST$epc,
                                             est.std=LIST$sepc.lv,
                                             cov.std=FALSE)
    }

    # power?
    if(power) {
        LIST$delta <- delta
        # FIXME: this is using epc in unstandardized metric
        #        this would be much more useful in standardized metric
        #        we need a standardize.est.all.reverse function...
        LIST$ncp <- (LIST$mi / LIST$epc^2) * (delta)^2
        LIST$power <- 1 - pchisq(qchisq((1.0 - alpha), df=1), 
                                 df=1, ncp=LIST$ncp)
        LIST$decision <- character( length(LIST$power) )

        # four possibilities (Table 6 in Saris, Satorra, van der Veld, 2009)
        mi.significant <- ifelse( 1 - pchisq(LIST$mi, df=1) < alpha,
                                  TRUE, FALSE )
        high.power <- LIST$power > high.power

        LIST$decision[ which(mi.significant &  high.power) ] <- "epc"
        LIST$decision[ which(mi.significant & !high.power) ] <- "***"
        LIST$decision[ which(!mi.significant & !high.power) ] <- "(i)"
    }


    # remove some columns
    LIST$free <- LIST$eq.id <- LIST$unco <- NULL
    LIST$mat <- LIST$row <- LIST$col <- LIST$id <- NULL
    #if(power) {
    #    LIST$epc <- NULL
    #    LIST$sepc.lv <- NULL
    #}

    class(LIST) <- c("lavaan.data.frame", "data.frame")

    # add eq constraints (if any)
    eq.idx <- which(partable$op == "==")
    if(length(eq.idx) > 0L) {
        N <- length(eq.idx)
        TMP <- data.frame(lhs = partable$lhs[eq.idx],
                           op = partable$op[eq.idx],
                          rhs = partable$rhs[eq.idx],
                          group = partable$group[eq.idx],
                          mi = rep(as.numeric(NA), N),
                          epc = rep(as.numeric(NA), N),
                          sepc.lv = rep(as.numeric(NA), N),
                          sepc.all = rep(as.numeric(NA), N),
                          sepc.nox = rep(as.numeric(NA), N) )
        if(!standardized) {
            TMP$sepc.lv <- TMP$sepc.all <- TMP$sepc.nox <- NULL
        }
        if(!is.null(LIST$mi.scaled)) {
            TMP$mi.scaled <- rep(as.numeric(NA), N)
        }
        LIST <- rbind(LIST, TMP)
    }

    # remove 'existing' parameters
    #TMP1 <- data.frame(lhs = partable$lhs,
    #                   op = partable$op,
    #                   rhs = partable$rhs,
    #                   group = partable$group)
    #TMP2 <- data.frame(lhs = LIST$lhs, 
    #                   op = LIST$op,
    #                   rhs = LIST$rhs,
    #                   group = LIST$group)
    #idx <- which(duplicated(rbind(TMP2, TMP1), fromLast = TRUE))
    #if(length(idx) > 0L) {
    #    LIST <- LIST[-idx,]
    #}
    if(ngroups == 1) LIST$group <- NULL
    
    # sort?
    if(sort.) {
        LIST <- LIST[order(LIST$mi, decreasing = TRUE),]
    }
    if(minimum.value > 0.0) {
        LIST <- LIST[!is.na(LIST$mi) & LIST$mi > minimum.value,]  
    }
    if(maximum.number < nrow(LIST)) {
        LIST <- LIST[seq_len(maximum.number),]
    }
    if(na.remove) {
        idx <- which(is.na(LIST$mi))
        if(length(idx) > 0) {
            LIST <- LIST[-idx,]
        }
    }
    if(!is.null(op)) {
        idx <- LIST$op %in% op
        if(length(idx) > 0) {
            LIST <- LIST[idx,]
        }
    }

    LIST
}

# aliases
modificationIndices <- modificationindices <- modindices
