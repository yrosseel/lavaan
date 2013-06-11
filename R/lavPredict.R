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
    lavPredict(object = object, newdata = newdata, type="FS")
})

# main function
lavPredict <- function(object, newdata=NULL, type="FS", se.fit=FALSE,
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
    }

    if(type == "fs") {
        out <- lav_predict_eta(object = object, data.obs = data.obs,
                               label = label)
    } else if(type == "mu") {
        out <- lav_predict_mu(object = object, data.obs = data.obs,
                              label = label)
    } else if(type == "fy") {
        out <- lav_predict_fy(object = object, data.obs = data.obs,
                              label = label)
    } else {
        stop("lavaan ERROR: type must be one of: FS MU FY")
    }

    if(object@Data@ngroups == 1L) {
        out <- out[[1L]]
    } else {
        out
    }

    out
}


## factor scores
## FIXME: 1) only works for LISREL representation for now!!!
##        2) if multiple group, return single data.frame with extra group col
##        3) only for continuous responses for now
lav_predict_eta <- function(object = NULL, data.obs = NULL, label = FALSE) {

    #if(object@SampleStats@missing.flag) {
    #    stop("FIXME: predict does not work with missing data (yet)!")
    #}
    if(object@Model@categorical) {
        stop("lavaan ERROR: predict does not work (yet) for categorical data")
    }

    G <- object@Data@ngroups
    nmat <- object@Model@nmat
    FS <- vector("list", length=G)

    for(g in 1:G) {
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST     <- object@Model@GLIST[ mm.in.group ]

        NFAC <- ncol(MLIST$lambda)

        LAMBDA <- MLIST$lambda
        PSI    <- MLIST$psi
        THETA  <- MLIST$theta
        BETA   <- MLIST$beta
        NU     <- MLIST$nu
        ALPHA  <- MLIST$alpha

        # nu, alpha? if not, set to colMeans/zero respectively
        if( is.null(ALPHA) ) {
            NU <- object@SampleStats@mean[[g]]
            ALPHA <- matrix(0, NFAC, 1)
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
        C = ( V.eta %*% t(LAMBDA) %*%
              solve( LAMBDA %*% V.eta %*% t(LAMBDA) + THETA ) )

        N <- nrow(data.obs[[g]])
        nvar <- ncol(object@SampleStats@cov[[g]])
        tmp1 <- matrix(NU, N, nvar, byrow=TRUE)
        tmp2 <- matrix(LAMBDA %*% E.eta, N, nvar, byrow=TRUE)
        tmp3 <- matrix(E.eta, N, NFAC, byrow=TRUE)

        FS[[g]] <-
            (tmp3 + t(C %*% t(data.obs[[g]]-tmp1-tmp2)))[,1:NFAC,drop=FALSE]

        if(label) {
            colnames(FS[[g]]) <- vnames(object@ParTable, type="lv", group=g)
        }

        class(FS[[g]]) <- c("lavaan.matrix", "matrix")
    }

    FS
}


# expectation of the response, conditional on the latent variables
lav_predict_mu <- function(object = NULL, data.obs = NULL,
                           label = FALSE) {

    # measurement part
    # y*_i = nu + lambda eta_i + K x_i + epsilon_i
    # 
    # where eta_i = predict(fit) = factor scores

    if(is.null(data.obs)) {
        data.obs <- object@Data@X
    }
    G <- object@Data@ngroups
    nmat <- object@Model@nmat
    MU <- vector("list", length=G)

    # we need the factor scores (per groups)
    ETA <- lav_predict_eta(object = object, data.obs = data.obs)

    for(g in 1:G) {
        mm.in.group <- 1:nmat[g] + cumsum(c(0,nmat))[g]
        MLIST     <- object@Model@GLIST[ mm.in.group ]

            NU <- MLIST$nu
        LAMBDA <- MLIST$lambda

        # nu? if not, set to sample means
        if( is.null(NU) ) {
            NU <- object@SampleStats@mean[[g]]
        }

        # FIXME!!! K and exogenous variables!

        # measurement model
        MU[[g]] <- sweep(ETA[[g]] %*% t(LAMBDA), MARGIN=2, NU, "+")

        if(label) {
            colnames(MU[[g]]) <- vnames(object@ParTable, type="ov", group=g)
        }

        class(MU[[g]]) <- c("lavaan.matrix", "matrix")

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

        class(FY[[g]]) <- c("lavaan.matrix", "matrix")

    }

    FY
}

