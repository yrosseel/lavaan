# lavPredict() contains a collection of `predict' methods
# the unifying theme is that they all rely on the factor scores
#
# fs: factor scores
# mu: expectation of the response (eg mean or probability=1)
# fy: conditional density of y (given lv's and exo) (under independence!)
#
# first version: YR 11 June 2013

# overload standard R function `predict'
setMethod("predict", "lavaan",
function(object, newdata=NULL) {
    lavPredict(object = object, newdata = newdata, type="FS", method="EBM")
})

# main function
lavPredict <- function(object, newdata=NULL, type="FS", method="EBM",
                       se.fit=FALSE,
                       label=TRUE) {

    stopifnot(inherits(object, "lavaan"))
    type <- tolower(type)
    if(type %in% c("factor", "factor.score", "factorscore"))
        type <- "fs"

    # need full data set supplied
    if(is.null(newdata)) {
        # use internal copy:
        if(is.null(object@Data@X[[1]])) {
            stop("no local copy of data; FIXME!")
        } else {
            data.obs <- object@Data@X
        }
        eXo <- object@Data@eXo
    } else {
        OV <- object@Data@ov
        newData <- lavData(data        = newdata,
                           group       = object@Data@group,
                           group.label = object@Data@group.label,
                           ov.names    = object@Data@ov.names,
                           ordered     = OV$name[ OV$type == "ordered" ],
                           ov.names.x  = object@Data@ov.names.x,
                           std.ov      = object@Data@std.ov,
                           missing     = object@Data@missing,
                           # warn      = FALSE,
                           allow.single.case = TRUE)
        data.obs <- newData@X
        eXo <- newData@eXo
    }

    # normal case?
    NORMAL <- all(object@Data@ov$type == "numeric")

    if(type == "fs") {
        if(NORMAL && method == "EBM") {
            out <- lav_predict_eta_normal(object = object, 
                data.obs = data.obs, eXo = eXo, label = label)
        } else {
        out <- lav_predict_eta_ebm(object = object, 
                data.obs = data.obs, eXo = eXo, label = label)
        }
    } else if(type == "mu") {
        out <- lav_predict_mu(object = object, data.obs = data.obs,
                              eXo = eXo, method = method, label = label)
    } else if(type == "fy") {
        out <- lav_predict_fy(object = object, data.obs = data.obs,
                              eXo = eXo, label = label)
    } else {
        stop("lavaan ERROR: type must be one of: FS MU FY")
    }

    # lavaan.matrix
    out <- lapply(X, "class<-", c("lavaan.matrix", "matrix"))

    if(object@Data@ngroups == 1L) {
        out <- out[[1L]]
    } else {
        out
    }

    out
}


## factor scores - EBM
lav_predict_eta_ebm <- function(object = NULL, data.obs = NULL, 
                                eXo = NULL, label = FALSE) {

    if(is.null(data.obs)) {
        data.obs <- object@Data@X
    }
    if(is.null(eXo)) {
        eXo <- object@Data@eXo
    }

    G <- object@Data@ngroups
    nmat <- object@Model@nmat
    FS <- vector("list", length=G)
    VETA <- computeVETA(object=object@Model, samplestats=object@SampleStats)
    VETA.inv <- lapply(VETA, solve)
    EETA <- computeEETA(object=object@Model, eXo=eXo)
    nfac <- ncol(VETA[[1]])
    TH <- computeTH(object@Model)
    th.idx <- object@Model@th.idx

    f.eta.i <- function(x, y.i, x.i, mu.i, g) {
        # conditional density of y, given eta.i(=x)
        log.fy <- lav_predict_fy_eta.i(object = object, y.i = y.i, x.i = x.i,
                                       eta.i = x, group = g, 
                                       TH = TH[[g]], th.idx = th.idx[[g]],
                                       log = TRUE)
        tmp <- as.numeric(0.5 * t(x - mu.i) %*% VETA.inv[[g]] %*% (x - mu.i))
        out <- tmp - sum(log.fy)
        #print(out)
        out
    }

    for(g in 1:G) {
        FS[[g]] <- matrix(0, nrow(data.obs[[g]]), nfac)
        # if no eXo, only one mu.i per group
        if(is.null(eXo[[g]])) {
            mu.i <- as.numeric(EETA[[g]])
            x.i <- NULL
        }

        # casewise for now
        N <- nrow(data.obs[[g]])
        for(i in 1:N) {
            if(!is.null(eXo[[g]])) {
                mu.i <- EETA[[g]][i,]
                 x.i <- eXo[[g]][i,]
            }
            y.i <- data.obs[[g]][i,]
            
            # find best values for eta.i
            out <- nlminb(start=numeric(nfac), objective=f.eta.i,
                            gradient=NULL, # for now
                            y.i=y.i, x.i=x.i, mu.i=mu.i, g=g)
            if(out$convergence == 0L) {
                eta.i <- out$par
            } else {
                eta.i <- rep(as.numeric(NA), nfac)
            }

            FS[[g]][i,] <- eta.i
        }

        if(label) {
            colnames(FS[[g]]) <- vnames(object@ParTable, type="lv", group=g)
        }

    }

    FS
}

## factor scores - normal case
lav_predict_eta_normal <- function(object = NULL, data.obs = NULL, 
                                   eXo = eXo, label = FALSE) {

    if(is.null(data.obs)) {
        data.obs <- object@Data@X
    }
    if(is.null(eXo)) {
        eXo <- object@Data@eXo
    }

    G <- object@Data@ngroups
    nmat <- object@Model@nmat
    FS <- vector("list", length=G)

    for(g in 1:G) {
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST     <- object@Model@GLIST[ mm.in.group ]

        nfac <- ncol(MLIST$lambda)

        LAMBDA <- MLIST$lambda
        PSI    <- MLIST$psi
        THETA  <- MLIST$theta
        BETA   <- MLIST$beta
        NU     <- MLIST$nu
        ALPHA  <- MLIST$alpha

        # nu, alpha? if not, set to colMeans/zero respectively
        if( is.null(ALPHA) ) {
            NU <- object@SampleStats@mean[[g]]
            ALPHA <- matrix(0, nfac, 1)
        }

        # beta?
        if( !is.null(BETA) ) {
            tmp <- -1.0 * BETA; diag(tmp) <- 1.0
            IB.inv <- solve(tmp)
            E.eta <- IB.inv %*% ALPHA
            V.eta <- IB.inv %*% PSI %*% t(IB.inv)
        } else {
            E.eta <- ALPHA
            V.eta <- PSI
        }

        # factor score coefficient matrix
        FSC = ( V.eta %*% t(LAMBDA) %*%
              solve( LAMBDA %*% V.eta %*% t(LAMBDA) + THETA ) )

        N <- nrow(data.obs[[g]])
        nvar <- ncol(object@SampleStats@cov[[g]])
        tmp1 <- matrix(NU, N, nvar, byrow=TRUE)
        tmp2 <- matrix(LAMBDA %*% E.eta, N, nvar, byrow=TRUE)
        tmp3 <- matrix(E.eta, N, nfac, byrow=TRUE)

        FS[[g]] <-
            (tmp3 + t(FSC %*% t(data.obs[[g]]-tmp1-tmp2)))[,1:nfac,drop=FALSE]

        if(label) {
            colnames(FS[[g]]) <- vnames(object@ParTable, type="lv", group=g)
        }

    }

    FS
}


# expectation of the response, conditional on the predicted latent
# variable scores
lav_predict_mu <- function(object = NULL, data.obs = NULL, eXo = NULL,
                           method = "EBM", label = FALSE) {

    # measurement part
    # y*_i = nu + lambda eta_i + K x_i + epsilon_i
    # 
    # where eta_i = predict(fit) = factor scores

    if(is.null(data.obs)) {
        data.obs <- object@Data@X
    }
    if(is.null(eXo)) {
        eXo <- object@Data@eXo
    }
    G <- object@Data@ngroups
    nmat <- object@Model@nmat
    MU <- vector("list", length=G)

    # normal case?
    NORMAL <- all(object@Data@ov$type == "numeric")

    # we need the factor scores (per groups)
    if(NORMAL && method == "EBM") {
        ETA <- lav_predict_eta_normal(object = object,
            data.obs = data.obs, eXo = eXo, label = label)
    } else if(method == "EBM") {
        ETA <- lav_predict_eta_ebm(object = object,
            data.obs = data.obs, eXo = eXo, label = label)
    } else {
        stop("lavaan ERROR: method ", method, " not (yet) supported of factor score prediction")
    }

    for(g in 1:G) {
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST     <- object@Model@GLIST[ mm.in.group ]

            NU <- MLIST$nu
        LAMBDA <- MLIST$lambda
         GAMMA <- MLIST$gamma
         DELTA <- MLIST$delta

        # remove dummy's from LAMBDA
        r.idx <- object@Model@ov.dummy.row.idx[[g]]
        c.idx <- object@Model@ov.dummy.col.idx[[g]]
        if(!is.null(c.idx)) {
            LAMBDA <- LAMBDA[,-c.idx,drop=FALSE]
        }

        # nu? if not, set to sample means
        if( is.null(NU) ) {
            NU <- object@SampleStats@mean[[g]]
        }

        # measurement model
        MU[[g]] <- sweep(ETA[[g]] %*% t(LAMBDA), MARGIN=2, NU, "+")

        # K + eXo?
        # note: K elements are in Gamma
        if(!is.null(eXo[[g]]) && !is.null(GAMMA)) {
            MU[[g]][,r.idx] <- 
                MU[[g]][,r.idx] + (eXo[[g]] %*% t(GAMMA[c.idx,,drop=FALSE]))
        }

        # delta?
        if(!is.null(DELTA)) {
            MU[[g]] <- sweep(MU[[g]], MARGIN=2, DELTA, "*")
        }

        if(label) {
            colnames(MU[[g]]) <- vnames(object@ParTable, type="ov", group=g)
        }

    }
    
    MU
}

# expectation of the response, conditional on the latent variables
# for given values 'eta' of a single observation (i) (eta.i)
lav_predict_mu_eta.i <- function(object = NULL, eta.i = NULL, x.i = NULL, 
                                 group = 1L) {

    # measurement part
    # y*_i = nu + lambda eta_i + K x_i + epsilon_i
    # 
    # where eta_i = predict(fit) = factor scores

    g <- group
    nmat <- object@Model@nmat
    nvar <- object@Model@nvar[g]
    mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
    MLIST     <- object@Model@GLIST[ mm.in.group ]
        NU <- MLIST$nu
    LAMBDA <- MLIST$lambda
     GAMMA <- MLIST$gamma
     DELTA <- MLIST$delta

    # remove dummy's from LAMBDA
    r.idx <- object@Model@ov.dummy.row.idx[[g]]
    c.idx <- object@Model@ov.dummy.col.idx[[g]]
    if(!is.null(c.idx)) {
        LAMBDA <- LAMBDA[,-c.idx,drop=FALSE]
    }

    # nu? if not, set to sample means
    if( is.null(NU) ) {
        NU <- object@SampleStats@mean[[g]]
    }

    # measurement model
    MU <- as.numeric(NU + LAMBDA %*% eta.i)

    # K + eXo?
    # note: K elements are in Gamma
    if(!is.null(x.i) && !is.null(GAMMA)) {
        MU.x <- numeric(nvar)
        MU.x[r.idx] <- (GAMMA[c.idx,,drop=FALSE] %*% x.i)
        MU <- MU + MU.x
    }

    # Delta?
    if(!is.null(DELTA)) {
        MU <- DELTA * MU
    }
    
    MU
}

# conditional density y -- assuming independence!!
# f(y_i | eta_i, x_i)
lav_predict_fy <- function(object = NULL, data.obs = NULL,
                           label = FALSE) {

    if(is.null(data.obs)) {
        data.obs <- object@Data@X
    }
    G <- object@Data@ngroups
    nmat <- object@Model@nmat
    FY <- vector("list", length=G)

    # we need the MUs (per group)
    MU <- lav_predict_mu(object = object, data.obs = data.obs)

    # all normal?
    NORMAL <- all(object@Data@ov$type == "numeric")

    for(g in 1:G) {
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST     <- object@Model@GLIST[ mm.in.group ]

        THETA <- MLIST$theta
        theta <- diag(THETA)

        if(NORMAL) {
            tmp <-  (data.obs[[g]] - MU[[g]])^2 
            tmp1 <- sweep(tmp, MARGIN=2, theta, "/")
            tmp2 <- exp( -0.5 * tmp1 )
            tmp3 <- sweep(tmp2, MARGIN=2, sqrt(2*pi*theta), "/")
            FY[[g]] <- tmp3
        } else {
            stop("not ready yet")
        }

        if(label) {
            colnames(FY[[g]]) <- vnames(object@ParTable, type="ov", group=g)
        }

    }

    FY
}

# conditional density y -- assuming independence!!
# f(y_i | eta_i, x_i)
# but for a single observation y_i (and x_i), for given values of eta_i
lav_predict_fy_eta.i <- function(object = NULL, y.i = NULL, x.i = NULL,
                                 eta.i = NULL, group = 1L, 
                                 TH = NULL, th.idx = NULL, log = TRUE) {

    g <- group
    nvar <- object@Model@nvar[g]
    nmat <- object@Model@nmat

    # we need the MU
    MU <- lav_predict_mu_eta.i(object = object, eta.i = eta.i, x.i = x.i,
                               group = group)

    # all normal?
    NORMAL <- all(object@Data@ov$type == "numeric")

    mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
    MLIST     <- object@Model@GLIST[ mm.in.group ]

    THETA <- MLIST$theta
    theta <- sqrt(diag(THETA))

    if(NORMAL) {
        FY <- dnorm(y.i, mean=MU, sd=theta, log = log)
    } else {
        FY <- numeric(nvar)
        for(v in seq_len(nvar)) {
            if(object@Data@ov$type[v] == "numeric") {
                ### FIXME!!! we can do all numeric vars at once!!
                FY[v] <- dnorm(y.i[v], mean=MU[v], sd=theta[v], log = log)
            } else if(object@Data@ov$type[v] == "ordered") {
                th.y <- TH[ th.idx == v ]; TH.Y <-  c(-Inf, th.y, Inf)
                k <- y.i[v]
                p1 <- pnorm( (TH.Y[ k + 1 ] - MU[v])/theta[v] )
                p2 <- pnorm( (TH.Y[ k     ] - MU[v])/theta[v] )
                prob <- (p1 - p2)
                if(log) {
                    FY[v] <- log(prob)
                } else {
                    FY[v] <- prob
                }
            } else {
                stop("lavaan ERROR: unknown type: ", 
                      object@Data@ov$type[v], " for variable", 
                      object@Data@ov$name[v])
            }
        }
    }

    FY
}

