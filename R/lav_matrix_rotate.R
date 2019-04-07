# collection of functions that deal with rotation matrices
# initial version YR 3 April 2019

# generate random orthogonal rotation matrix
lav_matrix_rotate_gen <- function(M = 10L) {

    # create random normal matrix
    tmp <- matrix(rnorm(M*M), nrow = M, ncol = M)

    # use QR decomposition
    qr.out <- qr(tmp)

    # extra 'Q' part
    out <- qr.Q(qr.out)

    out
}

# check if ROT is an orthogonal matrix if orthogonal = TRUE, or normal if
# orthogonal = FALSE
lav_matrix_rotate_check <- function(ROT = NULL, orthogonal = TRUE,
                                    tolerance = sqrt(.Machine$double.eps)) {

    # we assume ROT is a matrix
    M <- nrow(ROT)

    # crossprod
    RR <- crossprod(ROT)

    # target
    if(orthogonal) {
        # ROT^T %*% ROT = I
        target <- diag(M)
    } else {
        # diagonal should be 1
        target <- RR
        diag(target) <- 1
    }

    # compare for near-equality
    res <- all.equal(target = target, current = RR, tolerance = tolerance)

    # return TRUE or FALSE
    if(is.logical(res) && res) {
        out <- TRUE
    } else {
        out <- FALSE
    }

    out
}

# get weights vector needed to weight the rows using Kaiser normalization
lav_matrix_rotate_kaiser_weights <- function(A = NULL) {
    1/sqrt( rowSums(A * A) )
}

# get weights vector needed to weight the rows using Cureton & Mulaik (1975)
# standardization
# see also Browne (2001) page 128-129
lav_matrix_rotate_cm_weights <- function(A = NULL) {

    P <- nrow(A); M <- ncol(A)

    # first principal component of AA'
    A.eigen <- eigen(tcrossprod(A), symmetric = TRUE)
    a <- A.eigen$vectors[,1] * sqrt(A.eigen$values[1])

    Kaiser.weights <- 1/sqrt( rowSums(A * A) )
    a.star <- abs(a * Kaiser.weights) # always between 0 and 1

    m.sqrt.inv <- 1/sqrt(M)
    acos.m.sqrt.inv <- acos(m.sqrt.inv)

    delta <- numeric(P)
    delta[a.star < m.sqrt.inv] <- pi/2

    tmp <- (acos.m.sqrt.inv - acos(a.star))/(acos.m.sqrt.inv - delta) * (pi/2)

    # add constant (see Cureton & Mulaik, 1975, page 187)
    cm <- cos(tmp)*cos(tmp) + 0.001

    # final weights = weighted by Kaiser weights
    cm * Kaiser.weights
}

# high-level function to rotate a single matrix 'A'
lav_matrix_rotate <- function(A           = NULL,      # original matrix
                              orthogonal  = FALSE,     # default is oblique
                              method      = "quartimin",
                              method.args = list(),    # rotation method args
                              ROT         = NULL,      # initial rotation matrix
                              ROT.check   = TRUE,      # check if init ROT is ok
                              rstarts     = 0L,        # number of random starts
                              row.weights = "default", # row weighting
                              std.ov      = FALSE,     # rescale ov
                              ov.var      = NULL,      # ov variances
                              verbose     = FALSE,     # show iterations
                              algorithm   = "gpa",     # not used for now
                              frob.tol    = 0.00001,   # stopping tol gpa 
                              max.iter    = 1000L) {   # max gpa iterations

    # check A
    if(!inherits(A, "matrix")) {
        stop("lavaan ERROR: A does not seem to a matrix")
    }

    P <- nrow(A); M <- ncol(A)
    if(M < 2L) { # single dimension
        res <- list(LAMBDA = A, PHI = matrix(1, 1, 1), ROT = matrix(1, 1, 1),
                orthogonal = orthogonal, method = "none",
                method.args = list(), row.weights = "none",
                algorithm = "none", iter = 0L, converged = TRUE,
                method.value = 0)
        return(res)
    }

    # check ROT
    if(!is.null(ROT) && ROT.check) {
        if(!inherits(ROT, "matrix")) {
            stop("lavaan ERROR: ROT does not seem to a matrix")
        }
        if(nrow(ROT) != M) {
            stop("lavaan ERROR: nrow(ROT) = ", nrow(ROT),
                 " does not equal ncol(A) = ", M)
        }
        if(nrow(ROT) != ncol(ROT)) {
            stop("lavaan ERROR: nrow(ROT) = ", nrow(ROT),
                 " does not equal ncol(ROT) = ", ncol(ROT))
        }
        # rotation matrix? ROT^T %*% ROT = I
        RR <- crossprod(ROT)
        if(!lav_matrix_rotate_check(ROT, orthogonal = orthogonal)) {
            stop("lavaan ERROR: ROT does not look like a rotation matrix")
        }
    }

    # determine method function name
    method <- tolower(method)
    if(method %in% c("crawfer", "crawford.ferguson", "crawford-ferguson",
                     "crawfordferguson")) {
        method <- "cf"
    }
    if(method %in% c("varimax", "quartimax", "orthomax", "cf", "oblimin",
                     "quartimin", "geomin")) {
        method.fname <- paste("lav_matrix_rotate_", method, sep = "")
    } else if(method %in% c("cf-quartimax", "cf-varimax", "cf-equamax",
                            "cf-parsimax", "cf-facparsim")) {
        method.fname <- "lav_matrix_rotate_cf"
        method.args$gamma <- switch(method,
            "cf-quartimax" = 0,
            "cf-varimax"   = 1 / P,
            "cf-equamax"   = M / (2 * P),
            "cf-parsimax"  = (M - 1) / (P + M - 2),
            "cf-facparsim" = 1)
    }

    # check if rotation method exists
    check <- try(get(method.fname), silent = TRUE)
    if(inherits(check, "try-error")) {
        stop("lavaan ERROR: unknown rotation method: ", method.fname)
    }

    # set orthogonal option
    if(missing(orthogonal)) {
        # the default is oblique, except for varimax
        if(method %in% c("varimax")) {
            orthogonal <- TRUE
        } else {
            orthogonal <- FALSE
        }
    }

    # set row.weights
    row.weights <- tolower(row.weights)
    if(row.weights == "default") {
        # the default is "none", except for varimax
        if(method %in% c("varimax")) {
            row.weights <- "kaiser"
        } else {
            row.weights <- "none"
        }
    }

    # set algorithm
    algorithm <- "gpa" # for now



    # 1. compute row weigths

    # 1.a cov -> cor?
    if(std.ov) {
        A <- A * 1/sqrt(ov.var)
    }

    if(row.weights == "none") {
        weights <- rep(1.0, P)
    } else if(row.weights == "kaiser") {
        weights <- lav_matrix_rotate_kaiser_weights(A)
    } else if(row.weights == "cureton-mulaik") {
         weights <- lav_matrix_rotate_cm_weights(A)
    }
    A <- A * weights


    # 2. rotate 

    # multiple random starts?
    if(rstarts > 0L) {
        # TODO
    } else {
        # initial rotation matrix
        ROT <- diag(M)
  
        # Gradient Projection Algorithm
        if(algorithm == "gpa") {
            ROT <- lav_matrix_rotate_gpa(A = A, orthogonal = orthogonal,
                                         ROT = ROT,
                                         method.fname = method.fname,
                                         method.args  = method.args,
                                         verbose = verbose,
                                         frob.tol = frob.tol,
                                         max.iter = max.iter)
        }
        info <- attr(ROT, "info"); attr(ROT, "info") <- NULL
    }

    # final rotation
    if(orthogonal) {
        LAMBDA <- A %*% ROT
        PHI    <- diag(ncol(LAMBDA)) # correlation matrix == I
    } else {
        LAMBDA <- t(solve(ROT, t(A)))
        PHI    <- crossprod(ROT)     # correction matrix
    }

    # 3. undo row weighting
    LAMBDA <- LAMBDA / weights
    # 3.b undo cov -> cor
    if(std.ov) {
        LAMBDA <- LAMBDA * sqrt(ov.var)
    }

    # 4. reorder/reflect columns
    # 4.a reorder using sums-squares of the columns
    L2 <- LAMBDA * LAMBDA
    order.idx <- order(colSums(L2), decreasing = TRUE)
    LAMBDA <- LAMBDA[, order.idx, drop = FALSE]

    # 4.b reflect so that column sum is always positive
    SUM <- colSums(LAMBDA)
    neg.idx <- which(SUM < 0)
    if(length(neg.idx) > 0L) {
        LAMBDA[, neg.idx] <- -1 * LAMBDA[, neg.idx, drop = FALSE]
    }

    # 5. return results as a list
    res <- list(LAMBDA = LAMBDA, PHI = PHI, ROT = ROT,
                orthogonal = orthogonal, method = method,
                method.args = method.args, row.weights = row.weights)
    res <- c(res, info)

    res
}


# Gradient Projection Algorithm (Jennrich 2001, 2002)
# 
# - this is a translation of the SAS PROC IML code presented in the Appendix
#   of Bernaards & Jennrich (2005)
# - as the orthogonal and oblique algorithm are so similar, they are
#   combined in a single function
# - the default is oblique rotation
#
lav_matrix_rotate_gpa <- function(A            = NULL,   # original matrix
                                  orthogonal   = FALSE,  # default is oblique
                                  ROT          = NULL,   # initial rotation
                                  method.fname = NULL,   # criterion function
                                  method.args  = list(), # optional method args
                                  verbose      = FALSE,
                                  frob.tol     = 0.00001,
                                  max.iter     = 1000L) {

    # number of columns
    M <- ncol(A)

    # transpose of A (not needed for orthogonal)
    At <- t(A)

    # check ROT
    if(is.null(ROT)) {
        ROT <- diag(M)
    } 

    # set initial value of alpha to 1 (is 0.5 not more logical?)
    alpha <- 1

    # initial rotation
    if(orthogonal) {
        LAMBDA <- A %*% ROT
    } else {
        LAMBDA <- t(solve(ROT, At))
    }

    # using the current LAMBDA, evaluate the user-specified
    # rotation criteron; return Q (the criterion) and its gradient Gq
    Q <- do.call(method.fname, 
                 c(list(LAMBDA = LAMBDA), method.args, list(grad = TRUE)))
    Gq <- attr(Q, "grad"); attr(Q, "grad") <- NULL
    Q.current <- Q

    # compute gradient GRAD of f() at ROT from the gradient Gq of Q at LAMBDA
    # in a manner appropiate for orthogonal or oblique rotation
    if(orthogonal) {
        GRAD <- crossprod(A, Gq)
    } else {
        GRAD <- -1 * solve(t(ROT), crossprod(Gq, LAMBDA))
    }

    # start iterations
    converged <- FALSE
    for(iter in seq_len(max.iter + 1L)) {

        # compute projection Gp of GRAD onto the linear manifold tangent at
        # ROT to the manifold of orthogonal or normal (for oblique) matrices
        #
        # this projection is zero if and only if ROT is a stationary point of
        # f() restricted to the orthogonal/normal matrices
        if(orthogonal) {
            M <- crossprod(ROT, GRAD)
            SYMM <- (M + t(M))/2
            Gp <- GRAD - (ROT %*% SYMM)
        } else {
            Gp <- GRAD - t( t(ROT) * colSums(ROT * GRAD) )
        }

        # check Frobenius norm of Gp
        frob <- sqrt( sum(Gp * Gp) )

        # if verbose, print
        if(verbose) {
            cat("iter = ", sprintf("%4d", iter-1),
                " Q = ", sprintf("%12.7f", Q.current),
                " frob.log10 = ", sprintf("%12.7f", log10(frob)),
                " alpha = ", sprintf("%4.2f", alpha), "\n")
        }

        if(frob < frob.tol) {
            converged <- TRUE
            break
        }

        # update
        alpha <- 2*alpha
        for(i in seq_len(100)) { # make option?

            # step in the negative projected gradient direction
            # (note, the original algorithm in Jennrich 2001 used G, not Gp)
            X <- ROT - alpha * Gp

            if(orthogonal) {
                # use SVD to compute the projection ROTt of X onto the manifold
                # of orthogonal matrices
                svd.out <- svd(X)
                U <- svd.out$u
                V <- svd.out$v
                ROTt <- U %*% t(V)
            } else {
                # compute the projection ROTt of X onto the manifold
                # of normal matrices
                v <- 1/sqrt(apply(X^2, 2, sum))
                ROTt <- X %*% diag(v)
            }

            # rotate again
            if(orthogonal) {
                LAMBDA <- A %*% ROTt
            } else {
                LAMBDA <- t(solve(ROTt, At))
            }

            # evaluate criterion
            Q.new <- do.call(method.fname, c(list(LAMBDA = LAMBDA), 
                                             method.args, list(grad = TRUE)))
            Gq <- attr(Q.new, "grad"); attr(Q.new, "grad") <- NULL

            # check stopping criterion
            if(Q.new < Q.current - 0.5*frob*frob*alpha) {
                break
            } else {
                alpha <- alpha/2
            }

            if(i == 100) {
                warning("lavaan WARNING: half-stepping failed in GPA\n")
            }
        }

        # update
        ROT <- ROTt
        Q.current <- Q.new

        if(orthogonal) {
            GRAD <- crossprod(A, Gq)
        } else {
            GRAD <- -1 * solve(t(ROT), crossprod(Gq, LAMBDA))
        }
    }

    # final rotation
    if(orthogonal) {
        LAMBDA <- A %*% ROT
        PHI    <- diag(ncol(LAMBDA)) # correlation matrix == I
    } else {
        LAMBDA <- t(solve(ROT, At))
        PHI    <- crossprod(ROT)     # correction matrix
    }

    # algorithm information 
    info <- list(algorithm = "gpa",
                 iter = iter - 1L,
                 converged = converged,
                 method.value = Q.current)

    attr(ROT, "info") <- info

    ROT
}

