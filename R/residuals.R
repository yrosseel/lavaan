
setMethod("residuals", "lavaan",
function(object, type="raw", labels=TRUE) {
 
    # check estimator
    if(type %in% c("normalized", "standardized") && 
       object@Options$estimator != "ML") {
        stop("standardized and normalized residuals only availabe if estimator = ML (or MLF, MLR, MLM\n")
    }
    # NOTE: for some reason, Mplus does not compute the normalized/standardized
    # residuals if estimator = MLM !!!
 

    # check type
    if(!type %in% c("raw", "cor", "normalized", "standardized")) {
        stop("type must be one of \"raw\", \"cor\", \"normalized\" or \"standardized\"")
    }

    # check for 0 parameters if type == standardized
    if(type == "standardized" &&
       object@Fit@npar == 0) {
        stop("lavaan ERROR: can not compute standardized residuals if there are no free parameters in the model")
    }

    G <- object@Sample@ngroups
    meanstructure <- object@Model@meanstructure
    ov.names <- object@Data@ov.names

    # if type == standardized, we need VarCov and Delta
    if(type == "standardized") {
        # fixed.x idx?
        x.idx <- integer(0)
        if(object@Options$fixed.x) {
            x.idx <- match(vnames(object@User, "ov.x", group=1L),
                           object@Data@ov.names[[1L]]) ### FIXME!!!! will not
                                                       ### work for different
        }                                              ### models in groups

        if(length(x.idx) > 0L) {
            # we need to:
            # 1) to `augment' VarCov and Delta with the fixed.x  elements
            # 2) set cov between free and fixed.x elements in VarCov to zero

            # create 'augmented' User object (as if fixed.x=FALSE was used)
            augUser <- object@User
            idx <- which(augUser$exo > 0L)
            augUser$exo[       idx ] <- 0L
            augUser$free[      idx ] <- max(augUser$free) + 1:length(idx) 
            augUser$unco[idx ] <- max(augUser$unco) + 1:length(idx) 
            augModel <- Model(user           = augUser,
                              start          = object@Fit@est,
                              representation = object@Options$representation,
                              debug          = object@Options$debug)
            VarCov <- estimateVCOV(augModel, sample = object@Sample,
                                   options = object@Options)
            # set cov between free and fixed.x elements to zero
            ###
            ### FIXME: should we not do this on the information level,
            ###        *before* we compute VarCov?
            ###
            fixed.x.idx <- max(object@User$free) + 1:length(idx) 
            free.idx    <- 1:max(object@User$free)
            VarCov[free.idx, fixed.x.idx] <- 0.0
            VarCov[fixed.x.idx, free.idx] <- 0.0

            Delta  <- computeDelta(augModel)
        } else {
            VarCov <- estimateVCOV(object@Model, sample = object@Sample,
                                   options = object@Options)
            Delta  <- computeDelta(object@Model)
        }   
    }

    R <- vector("list", length=G)
    for(g in 1:G) {
        # sample moments
        if(!object@Sample@missing[[g]]$flag) {
            S <- object@Sample@cov[[g]]
            M <- object@Sample@mean[[g]]
        } else {
            S <- object@Sample@missing[[g]]$sigma
            M <- object@Sample@missing[[g]]$mu
        }
        if(!meanstructure) {
            M <- numeric( length(M) )
        }
        nvar <- ncol(S)

        # raw residuals (for this group
        if(type == "cor") {
            R[[g]]$cov  <- cov2cor(S) - cov2cor(object@Fit@Sigma.hat[[g]])
        } else {
            R[[g]]$cov  <- S - object@Fit@Sigma.hat[[g]]
        }
        R[[g]]$mean <- M - object@Fit@Mu.hat[[g]]
        if(labels) {
            rownames(R[[g]]$cov) <- colnames(R[[g]]$cov) <- ov.names[[g]]
        }

        if(type == "normalized" || type == "standardized") {
         
            # compute normalized residuals
            N <- object@Sample@nobs[[g]]; nvar <- length(R[[g]]$mean)
            idx.mean <- 1:nvar

            if(object@Options$se == "standard" ||
               object@Options$se == "none") {
                dS <- diag(S)
                Var.mean <- Var.sample.mean <- dS / N 
                Var.cov  <- Var.sample.cov  <- (tcrossprod(dS) + S^2) / N
                # this is identical to solve(A1)/N for complete data!!
            } else if(object@Options$se == "robust.mlr" ||
                      object@Options$se == "robust.mlm") {
                A1 <- compute.A1.sample(sample=object@Sample, group=g, 
                                        meanstructure=meanstructure,
                                        information=object@Options$information)
                B1 <- compute.B1.sample(sample=object@Sample, group=g,
                                        meanstructure=meanstructure)
                Info <- (solve(A1) %*% B1 %*% solve(A1)) / N
                Var.mean <- Var.sample.mean <- diag(Info)[idx.mean]
                Var.cov  <- Var.sample.cov  <- vech.reverse(diag(Info)[-idx.mean])
            } else if(object@Options$se == "first.order") {
                B1 <- compute.B1.sample(sample=object@Sample, group=g,
                                       meanstructure=meanstructure)
                Info <- solve(B1) / N
                Var.mean <- Var.sample.mean <- diag(Info)[idx.mean]
                Var.cov  <- Var.sample.cov  <- vech.reverse(diag(Info)[-idx.mean])
            }
        }

        if(type == "standardized") {

            Var.model <- diag(Delta[[g]] %*% VarCov %*% t(Delta[[g]]))
 
            if(meanstructure) {
                Var.model.mean <- Var.model[idx.mean]
                Var.model.cov  <- vech.reverse(Var.model[-idx.mean])
            } else {
                Var.model.mean <- rep(0, nvar)
                Var.model.cov  <- vech.reverse(Var.model)
            }

            Var.mean <- (Var.sample.mean - Var.model.mean)
            Var.cov  <- (Var.sample.cov  - Var.model.cov )

            # not for fixed x covariates
            if(length(x.idx) > 0L) {
                Var.mean[x.idx] <- 1.0
                Var.cov[x.idx,x.idx] <- 1.0
            }

            # avoid negative variances
            Var.mean[which(Var.mean < 0)] <- NA
            Var.cov[ which(Var.cov  < 0)] <- NA
        }

        # normalize/standardize
        if(type == "normalized" || type == "standardized") {
            # avoid small number (< 1.0e-15) to be divided
            # by another small number and get bigger...
            # FIXME!!!
            tol <- 1.0e-5
            R[[g]]$mean[ which(abs(R[[g]]$mean) < tol)] <- 0.0
            R[[g]]$cov[ which(abs(R[[g]]$cov) < tol)] <- 0.0
            
            R[[g]]$mean <- R[[g]]$mean / sqrt( Var.mean )
            R[[g]]$cov  <- R[[g]]$cov  / sqrt( Var.cov  )
        }

        # prepare for pretty printing
        R[[g]]$mean <- as.numeric(R[[g]]$mean)
        if(labels) names(R[[g]]$mean) <- ov.names[[g]]
        class(R[[g]]$mean) <- c("lavaan.vector", "numeric")
        class(R[[g]]$cov) <- c("lavaan.matrix.symmetric", "matrix")
    }

    # replace 'cov' by 'cor' if type == "cor"
    if(type == "cor") {
        R <- lapply(R, "names<-", c("cor", "mean") )
    }


    if(G == 1) {
        R <- R[[1]]
    } else {
        names(R) <- unlist(object@Sample@group.label)
    }

    R
})

setMethod("resid", "lavaan",
function(object, type="raw") {
    residuals(object, type=type)
})

