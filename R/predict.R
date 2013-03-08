
## currently, the 'predict' function is used only to 
## compute the 'factor scores' of the latent variables
## using the Emperical Bayes approach (aka the `regression approach')
## 

## FIXME: 1) only works for LISREL representation for now!!!
##        2) if multiple group, return single data.frame with extra group col
##        3) no support for missing data...
##        4) much more testing!

setMethod("predict", "lavaan",
function(object, newdata=NULL, ...) {

    #if(object@SampleStats@missing.flag) {
    #    stop("FIXME: predict does not work with missing data (yet)!")
    #}
    if(object@Model@categorical) {
        stop("lavaan ERROR: predict does not work (yet) for categorical data")
    }

    G <- object@Data@ngroups
    nmat <- object@Model@nmat
    FS <- vector("list", length=G)

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
                           missing     = object@Data@missing)
        data.obs <- newData@X
    }

    for(g in 1:G) {

        lv.names <- vnames(object@ParTable, type="lv", group=g)

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

        C = ( V.eta %*% t(LAMBDA) %*% 
              solve( LAMBDA %*% V.eta %*% t(LAMBDA) + THETA ) )

        N <- nrow(data.obs[[g]])
        nvar <- ncol(object@SampleStats@cov[[g]])
        tmp1 <- matrix(NU, N, nvar, byrow=TRUE)
        tmp2 <- matrix(LAMBDA %*% E.eta, N, nvar, byrow=TRUE)
        tmp3 <- matrix(E.eta, N, NFAC, byrow=TRUE)

        nfac <- length(lv.names)
        FS[[g]] <- 
            (tmp3 + t(C %*% t(data.obs[[g]]-tmp1-tmp2)))[,1:nfac,drop=FALSE]
        colnames(FS[[g]]) <- lv.names

        class(FS[[g]]) <- c("lavaan.matrix", "matrix")
    }

    if(G == 1) {
        FS <- FS[[1]]
    }

    FS 
})

