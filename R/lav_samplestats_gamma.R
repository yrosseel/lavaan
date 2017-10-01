# YR 21 March 2015
# new approach to compute 'Gamma': the asymptotic variance matrix of 
#                                  sqrt{N} times the
#                                  observed sample statistics (means + varcov)
#
# Gamma = N x ACOV[ ybar, vech(S) ]
#       = NACOV[ ybar, vech(S) ]
#
# - one single function for mean + cov
# - handle 'fixed.x' exogenous covariates
# - YR 3 Dec 2015: allow for conditional.x = TRUE

# generic public function
# input for lavGamma can be lavobject, lavdata, data.frame, or matrix
lavGamma <- function(object, group = NULL, missing = "listwise",
                     ov.names.x = NULL, fixed.x = FALSE, conditional.x = FALSE,
                     meanstructure = FALSE, slopestructure = FALSE,
                     Mplus.WLS = FALSE, add.labels) {
    
    if(inherits(object, "lavaan")) {
        lavdata <- object@Data
        if(missing(missing)) {
            missing <- object@Options$missing
            if(missing != "listwise") {
                ### FIXME!!!!!
                return(NULL) # for now!
            }
        } else {
            missing <- "listwise"
        }
    } else if(inherits(object, "lavData")) {
        lavdata <- object
    } else if(inherits(object, "data.frame") ||
              inherits(object, "matrix")) {
        NAMES <- names(object)
        if(!is.null(NAMES) && !is.null(group)) {
            NAMES <- NAMES[- match(group, NAMES)]
        }
        lavdata <- lavData(data = object, group = group,
                           ov.names = NAMES, ordered = NULL,
                           ov.names.x = ov.names.x, 
                           lavoptions = list(warn = FALSE,
                                             missing = missing))
    } else {
        stop("lavaan ERROR: lavGamma can not handle objects of class ",
             paste(class(object), collapse= " "))
    }

    # extract data
    Y <- lavdata@X

    # x-covariates?
    x.idx <- lapply(seq_len(lavdata@ngroups),
                    function(g) match(lavdata@ov.names.x[[g]],
                                      lavdata@ov.names[[g]]) )

    OUT <- lapply(seq_len(lavdata@ngroups),
              function(g) lav_samplestats_Gamma(Y              = Y[[g]],
                                                x.idx          = x.idx[[g]],
                                                fixed.x        = fixed.x,
                                                conditional.x  = conditional.x,
                                                meanstructure  = meanstructure,
                                                slopestructure = slopestructure,
                                                Mplus.WLS      = Mplus.WLS))
  
    OUT
}

# NOTE:
#  - three types:
#       1) plain         (conditional.x = FALSE, fixed.x = FALSE)
#       2) fixed.x       (conditional.x = FALSE, fixed.x = TRUE)
#       3) conditional.x (conditional.x = TRUE)
#  - if conditional.x = TRUE, we ignore fixed.x (can be TRUE or FALSE)

# NORMAL-THEORY
lav_samplestats_Gamma_NT <- function(Y              = NULL,
                                     wt             = NULL,
                                     COV            = NULL,
                                     MEAN           = NULL,
                                     rescale        = FALSE,
                                     x.idx          = integer(0L),
                                     fixed.x        = FALSE,
                                     conditional.x  = FALSE,
                                     meanstructure  = FALSE,
                                     slopestructure = FALSE) {

    # check arguments
    if(length(x.idx) == 0L) {
        conditional.x <- FALSE
        fixed.x       <- FALSE
    }

    if(is.null(COV)) {
        stopifnot(!is.null(Y))
    
        # coerce to matrix
        Y <- unname(as.matrix(Y)); N <- nrow(Y)
        if(is.null(wt)) {
            COV <- cov(Y)
        } else {
            out <- stats::cov.wt(Y, wt = wt, method = "ML")
            COV <- out$cov
        }
    }

    if(rescale && is.null(wt)) {
        COV <- COV * (N-1) / N # ML version
    }

    if(conditional.x && length(x.idx) > 0L && is.null(MEAN) &&
       (meanstructure || slopestructure)) {
       stopifnot(!is.null(Y))
       if(is.null(wt)) {
           MEAN <- colMeans(Y)
       } else {
           MEAN <- out$center
       }
    }

    # rename
    S <- COV
    M <- MEAN

    # unconditional
    if(!conditional.x) {

        # unconditional - stochastic x
        if(!fixed.x) {
            Gamma <- 2*lav_matrix_duplication_ginv_pre_post(S %x% S)
            if(meanstructure) {
                Gamma <- lav_matrix_bdiag(S, Gamma)
            }

        # unconditional - fixed x
        } else {
            # handle fixed.x = TRUE

            # cov(Y|X) = A - B C^{-1} B' 
            # where A = cov(Y), B = cov(Y,X), C = cov(X)
            A <- S[-x.idx, -x.idx, drop=FALSE]
            B <- S[-x.idx,  x.idx, drop=FALSE]
            C <- S[ x.idx,  x.idx, drop=FALSE]
            YbarX <- A - B %*% solve(C, t(B))

            # reinsert YbarX in Y+X (residual) covariance matrix
            YbarX.aug <- matrix(0, nrow = NROW(S), ncol = NCOL(S))
            YbarX.aug[ -x.idx, -x.idx ] <- YbarX

            # take difference
            R <- S - YbarX.aug

            Gamma.S <- 2*lav_matrix_duplication_ginv_pre_post(S %x% S)
            Gamma.R <- 2*lav_matrix_duplication_ginv_pre_post(R %x% R)
            Gamma <- Gamma.S - Gamma.R

            if(meanstructure) {
                Gamma <- lav_matrix_bdiag(YbarX.aug, Gamma)
            }
        }

    } else {
        # conditional.x 

        # 4 possibilities:
        # - no meanstructure, no slopes
        # -    meanstructure, no slopes
        # - no meanstructure, slopes
        # -    meanstructure, slopes

        # regress Y on X, and compute covariance of residuals 'R'
        A <- S[-x.idx, -x.idx, drop=FALSE]
        B <- S[-x.idx,  x.idx, drop=FALSE]
        C <- S[ x.idx,  x.idx, drop=FALSE]
        Cov.YbarX <- A - B %*% solve(C) %*% t(B)
        Gamma <- 2*lav_matrix_duplication_ginv_pre_post(Cov.YbarX %x% Cov.YbarX)

        if(meanstructure || slopestructure) {
            MY <- M[-x.idx]; MX <- M[x.idx]
            C3 <- rbind(c(1,MX),
                        cbind(MX, C + tcrossprod(MX)))
            #B3 <- cbind(MY, B + tcrossprod(MY,MX))
        }

        if(meanstructure) {
            if(slopestructure) {
                A11 <- solve(C3) %x% Cov.YbarX
            } else {
                A11 <- solve(C3)[1, 1, drop=FALSE] %x% Cov.YbarX
            }
        } else {
            if(slopestructure) {
                A11 <- solve(C3)[-1, -1, drop=FALSE] %x% Cov.YbarX
            } else {
                A11 <- matrix(0,0,0)
            }
        }

        Gamma <- lav_matrix_bdiag(A11, Gamma)
    }
   
    Gamma
}

# NOTE:
#  - three types:
#       1) plain         (conditional.x = FALSE, fixed.x = FALSE)
#       2) fixed.x       (conditional.x = FALSE, fixed.x = TRUE)
#       3) conditional.x (conditional.x = TRUE)
#  - if conditional.x = TRUE, we ignore fixed.x (can be TRUE or FALSE)

# ADF THEORY
lav_samplestats_Gamma <- function(Y, 
                                  x.idx          = integer(0L),
                                  fixed.x        = FALSE,
                                  conditional.x  = FALSE,
                                  meanstructure  = FALSE,
                                  slopestructure = FALSE, 
                                  Mplus.WLS      = FALSE,
                                  add.attributes = FALSE) {
    # coerce to matrix
    Y <- unname(as.matrix(Y)); N <- nrow(Y)

    # check arguments
    if(length(x.idx) == 0L) {
        conditional.x <- FALSE
        fixed.x <- FALSE
    }
    if(Mplus.WLS) {
        stopifnot(!conditional.x, !fixed.x)
    }

    if(!conditional.x) {
        # center only, so we can use crossprod instead of cov
        Yc <- base::scale(Y, center = TRUE, scale = FALSE)
        p <- ncol(Y)

        # create Z where the rows_i contain the following elements:
        #  - intercepts (if meanstructure is TRUE)
        #  - vech(Yc_i' %*% Yc_i) where Yc_i are the residuals
        idx1 <- lav_matrix_vech_col_idx(p)
        idx2 <- lav_matrix_vech_row_idx(p)
        if(meanstructure) {
            Z <- cbind(Yc, Yc[,idx1] * Yc[,idx2])
        } else {
            Z <- Yc[,idx1] * Yc[,idx2]
        }

        # handle fixed.x = TRUE
        if(fixed.x) {
            YX <- Yc
            # here, we do not need the intercepts, data is centered
            QR <- qr(Yc[, x.idx, drop = FALSE])
            RES <- qr.resid(QR, Yc[,-x.idx, drop = FALSE])
            # substract residuals from original Yc's
            YX[, -x.idx] <- Yc[, -x.idx, drop = FALSE] - RES
            if(meanstructure) {
                Z2 <- cbind(YX, YX[,idx1] * YX[,idx2])
            } else {
                Z2 <- YX[,idx1] * YX[,idx2]
            }
            # substract Z2 from original Z
            Z <- Z - Z2
        }

        #Gamma = (N-1)/N * cov(Z, use = "pairwise")
        # we center so we can use crossprod instead of cov
        Zc <- base::scale(Z, center = TRUE, scale = FALSE)

        # note: centering is the same as substracting lav_matrix_vech(S),
        # where S is the sample covariance matrix (divided by N)

        if(anyNA(Zc)) {
            Gamma <- lav_matrix_crossprod(Zc) / N
        } else {
            Gamma <- base::crossprod(Zc) / N
        }


    } else {
        # conditional.x 

        # 4 possibilities:
        # - no meanstructure, no slopes
        # -    meanstructure, no slopes
        # - no meanstructure, slopes
        # -    meanstructure, slopes

        # regress Y on X, and compute residuals
        X    <- cbind(1, Y[,  x.idx, drop = FALSE])
        QR   <- qr(X)
        RES  <- qr.resid(QR, Y[, -x.idx, drop = FALSE])
        p    <- ncol(RES)

        idx1 <- lav_matrix_vech_col_idx(p)
        idx2 <- lav_matrix_vech_row_idx(p)

        if(meanstructure || slopestructure) {
            XtX.inv <- unname(solve(crossprod(X)))
            Xi <- (X %*% XtX.inv) * N ## FIXME, shorter way?
            ncX <- NCOL(X); ncY <- NCOL(RES)
        }

        if(meanstructure) {
            if(slopestructure) {
                 Xi.idx <- rep(seq_len(ncX), each  = ncY)
                Res.idx <- rep(seq_len(ncY), times = ncX)
                Z <- cbind( Xi[, Xi.idx, drop = FALSE] *
                           RES[,Res.idx, drop = FALSE],
                           RES[,   idx1, drop = FALSE] * 
                           RES[,   idx2, drop = FALSE] )
            } else {
                Xi.idx <- rep(1L, each = ncY)
                Z <- cbind( Xi[, Xi.idx ,drop = FALSE] *
                           RES,
                           RES[,   idx1, drop = FALSE] * 
                           RES[,   idx2, drop = FALSE] )
            }
        } else {
            if(slopestructure) {
                 Xi.idx <- rep(seq_len(ncX), each  = ncY)
                 Xi.idx <- Xi.idx[ -seq_len(ncY) ]
                Res.idx <- rep(seq_len(ncY), times = (ncX - 1L))
                Z <- cbind( Xi[, Xi.idx, drop = FALSE] *
                           RES[,Res.idx, drop = FALSE],
                           RES[,   idx1, drop = FALSE] * 
                           RES[,   idx2, drop = FALSE] )
            } else {
                Z <- RES[,idx1,drop = FALSE] * RES[,idx2,drop = FALSE]
            }
        }

        Zc <- base::scale(Z, center = TRUE, scale = FALSE)
        if(anyNA(Zc)) {
            Gamma <- lav_matrix_crossprod(Zc) / N
        } else {
            Gamma <- base::crossprod(Zc) / N
        }
    }
   

    # only to mimic Mplus when estimator = "WLS"
    if(Mplus.WLS && !fixed.x && !conditional.x) {
        # adjust G_22 (the varcov part)
        S <- cov(Y, use = "pairwise")
        w <- lav_matrix_vech(S)
        w.biased <- (N-1)/N * w
        diff <- outer(w,w) - outer(w.biased, w.biased)
        if(meanstructure) {
            Gamma[-seq_len(p), -seq_len(p)] <-
                Gamma[-seq_len(p), -seq_len(p), drop = FALSE] - diff
        } else {
            Gamma <- Gamma - diff
        }

        if(meanstructure) {
            # adjust G_12/G_21 (third-order)
            # strange rescaling?
            N1 <- (N - 1) / N
            Gamma[seq_len(p),-seq_len(p)] <- Gamma[seq_len(p),-seq_len(p)] * N1
            Gamma[-seq_len(p),seq_len(p)] <- Gamma[-seq_len(p),seq_len(p)] * N1
        }
    }

    Gamma
}

