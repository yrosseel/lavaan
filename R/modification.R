modificationIndices <- modificationindices <- modindices <- function(object, 
    standardized=TRUE, power=FALSE, delta=0.1, alpha=0.05, high.power=0.75) {

    if(power) standardized <- TRUE

    # get LIST parameter list
    LIST <- getUserListFull(object@User)
    LIST$free <- 0L; LIST$eq.id <- 0L; LIST$unco <- 0L
    LIST <- as.data.frame(LIST)

    # fill in USER information
    user <- object@User; N <- length(user$lhs)
    for(i in 1:N) {
        LIST.idx <- which(LIST$lhs == user$lhs[i] &
                          LIST$op  == user$op[i] &
                          LIST$rhs == user$rhs[i] &
                          LIST$group == user$group[i])
        LIST$free[LIST.idx] <- user$free[i]
        LIST$eq.id[LIST.idx] <- user$eq.id[i]
        LIST$unco[LIST.idx] <- user$unco[i]
    }

    # add matrix representation
    if(object@Model@representation == "LISREL") {
        REP <- representation.LISREL(user=object@User, target=LIST,
                                     extra=FALSE)
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

    # here we should remove elements that will produce NA anyways...
    # eg - first indicator of factors
    #    - regressions that are already free covariances
    # TODO
 

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
    Delta <- computeDelta(object@Model, m.el.idx=m.el.idx, x.el.idx=x.el.idx)

    # compute information matrix
    E <- computeExpectedInformation(object@Model, sample=object@Sample,
                                    estimator=object@Options$estimator,
                                    Delta=Delta)
    Q <- (1/object@Sample@ntotal) * E

    # list!
    DX <- computeGradient(object@Model, GLIST=NULL, sample=object@Sample,
                          type="allofthem", 
                          estimator=object@Options$estimator,
                          group.weight=TRUE)

    # flatten list DX to a vector dx
    dx <- numeric(0)
    for(mm in 1:length(object@Model@GLIST)) {
        dx[ x.el.idx[[mm]] ] <- DX[[mm]][ m.el.idx[[mm]] ] 
    }
                             
    # Saris, Satorra & Sorbom 1987
    # partition Q into Q_11, Q_22 and Q_12/Q_21
    # which elements of Q correspond with 'free' and 'nonfree' parameters?
    all.idx      <- LIST$id
    free.idx     <- LIST$id[LIST$free > 0L & !duplicated(LIST$free)]
    eq.idx       <- LIST$id[LIST$eq.id > 0L]
    eq.id        <- LIST$eq.id[ eq.idx ]

    # nonfree
    if(length(free.idx) > 0L) nonfree.idx <- all.idx[ -free.idx ]
    # if eqs.idx, remove them from free.idx
    if(length(eq.idx) > 0L) free.idx <- free.idx[-which(free.idx %in% eq.idx)]
 
    # check: nonfree.idx, free.idx and eq.idx should NOT overlap!

    # partition Q
    Q11 <- Q[nonfree.idx, nonfree.idx]
    Q12 <- Q[nonfree.idx,    free.idx]
    Q21 <- Q[free.idx,    nonfree.idx]
    Q22 <- Q[free.idx,       free.idx]; Q22.inv <- solve(Q22)

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
    if(object@Model@eq.constraints) {
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
    d <- (-1 * object@Sample@ntotal) * dx
    # needed?
    d[which(abs(d) < 1e-15)] <- 1.0
    epc <- mi/d

    # FIXME: epc for equality constraints must be adapted!?
    # we need a reference for this!!!!
    if(object@Model@eq.constraints) {
        for(i in 1:length(eq.idx)) {
            # this is just a temporary solution
            # until we get it right
            neq <- length(eq.idx[eq.id == eq.id[i]])
            epc[eq.idx[i]] <- epc[eq.idx[i]] / neq * (neq - 1)
        }
    }

    LIST$mi <- mi
    LIST$epc <- epc

    # remove some rows
    #idx <- which(LIST$free > 0L & !duplicated(LIST$free) & !LIST$eq.id > 0L)
    #LIST <- LIST[-idx,]

    # standardize?
    if(standardized) {
        LIST$sepc.lv  <- standardize.est.lv(object, user=LIST, est=LIST$epc,
                                            cov.std=FALSE)
        LIST$sepc.all <- standardize.est.all(object, user=LIST, est=LIST$epc,
                                             est.std=LIST$sepc.lv,
                                             cov.std=FALSE)
        LIST$sepc.nox <- standardize.est.all.nox(object,user=LIST,est=LIST$epc,
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
    if(power) {
        LIST$epc <- NULL
        LIST$sepc.lv <- NULL
    }
    if(ngroups == 1) LIST$group <- NULL

    class(LIST) <- c("lavaan.data.frame", "data.frame")

    LIST
}


