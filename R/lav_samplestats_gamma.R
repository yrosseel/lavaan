# YR 21 March 2015
# new approach to compute 'Gamma': the asymptotic variance matrix of the
#                                  observed sample statistics (means + varcov)
# - one single function for mean + cov
# - handle 'fixed.x' exogenous covariates

# generic public function
# input for lavGamma can be lavobject, lavdata, data.frame, or matrix
lavGamma <- function(object, group = NULL, missing = "listwise",
                     ov.names.x = NULL, meanstructure = FALSE,
                     Mplus.WLS = FALSE) {
    
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
                           ov.names.x = ov.names.x, warn = FALSE,
                           missing = missing)
    } else {
        stop("lavaan ERROR: lavGamma can not handle objects of class ",
             paste(class(object), collapse= " "))
    }

    # extract data
    Y <- lavdata@X

    # fixed.x?
    x.idx <- lapply(seq_len(lavdata@ngroups),
                    function(g) match(lavdata@ov.names.x[[g]],
                                      lavdata@ov.names[[g]]) )

    OUT <- lapply(seq_len(lavdata@ngroups),
              function(g) lav_samplestats_Gamma(Y = Y[[g]],
                                                x.idx = x.idx[[g]],
                                                meanstructure = meanstructure,
                                                Mplus.WLS = Mplus.WLS))
  
    OUT
}

lav_samplestats_Gamma <- function(Y, x.idx = integer(0L),
                                  meanstructure = FALSE,
                                  Mplus.WLS = FALSE) {
    # coerce to matrix
    Y <- as.matrix(Y)
    N <- nrow(Y); p <- ncol(Y)

    # center
    Yc <- base::scale(Y, center = TRUE, scale = FALSE)

    # create Z where the rows_i contain the following elements:
    #  - Yc_i (if meanstructure if TRUE)
    #  - vech(Yc_i' %*% Yc_i)
    idx1 <- lav_matrix_vech_col_idx(p)
    idx2 <- lav_matrix_vech_row_idx(p)
    if(meanstructure) {
        Z <- cbind(Yc, Yc[,idx1] * Yc[,idx2])
    } else {
        Z <- Yc[,idx1] * Yc[,idx2]
    }

    if(length(x.idx) > 0L) {
        YX <- Yc
        # regress Y on X, and compute residuals
        # note: since YC is centered, we do not need an intercept
        RES <- qr.resid(qr(cbind(Yc[, x.idx, drop = FALSE])),
                                 Yc[,-x.idx, drop = FALSE])
        # substract residuals from original Y's
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
    Zc <- base::scale(Z, center = TRUE, scale = FALSE)
    if(anyNA(Zc)) {
        Gamma <- lav_matrix_crossprod(Zc) / N
    } else {
        Gamma <- base::crossprod(Zc) / N
    }

    # only to mimic Mplus when estimator = "WLS"
    if(Mplus.WLS) {
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

