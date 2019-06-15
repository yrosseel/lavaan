# rotation algorithms
#
# YR  3 April 2019 -- gradient projection algorithm
# YR 21 April 2019 -- pairwise rotation algorithm

# main function to rotate a single matrix 'A'
lav_matrix_rotate <- function(A           = NULL,      # original matrix
                              orthogonal  = FALSE,     # default is oblique
                              method      = "geomin",  # default rot method
                              method.args = list( geomin.epsilon = 0.01,
                                                  orthomax.gamma = 1,
                                                  cf.gamma       = 0,
                                                  oblimin.gamma  = 0,
                                                  target      = matrix(0,0,0),
                                                  target.mask = matrix(0,0,0) ),
                              init.ROT    = NULL,      # initial rotation matrix
                              init.ROT.check = TRUE,   # check if init ROT is ok
                              rstarts     = 100L,      # number of random starts
                              row.weights = "default", # row weighting
                              std.ov      = FALSE,     # rescale ov
                              ov.var      = NULL,      # ov variances
                              warn        = TRUE,      # show warnings?
                              verbose     = FALSE,     # show iterations
                              algorithm   = "gpa",     # not used for now
                              reflect     = TRUE,      # refect sign
                              order.lv.by = "index",   # how to order the lv's
                              gpa.tol     = 0.00001,   # stopping tol gpa
                              tol         = 1e-07,     # stopping tol others
                              keep.rep    = FALSE,     # store replications
                              max.iter    = 10000L) {  # max gpa iterations

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

    # check init.ROT
    if(!is.null(init.ROT) && init.ROT.check) {
        if(!inherits(init.ROT, "matrix")) {
            stop("lavaan ERROR: init.ROT does not seem to a matrix")
        }
        if(nrow(init.ROT) != M) {
            stop("lavaan ERROR: nrow(init.ROT) = ", nrow(init.ROT),
                 " does not equal ncol(A) = ", M)
        }
        if(nrow(init.ROT) != ncol(init.ROT)) {
            stop("lavaan ERROR: nrow(init.ROT) = ", nrow(init.ROT),
                 " does not equal ncol(init.ROT) = ", ncol(init.ROT))
        }
        # rotation matrix? init.ROT^T %*% init.ROT = I
        RR <- crossprod(init.ROT)
        if(!lav_matrix_rotate_check(init.ROT, orthogonal = orthogonal)) {
            stop("lavaan ERROR: init.ROT does not look like a rotation matrix")
        }
    }

    # determine method function name
    method <- tolower(method)
    if(method %in% c("crawfer", "crawford.ferguson", "crawford-ferguson",
                     "crawfordferguson")) {
        method <- "cf"
    }
    if(method %in% c("varimax", "quartimax", "orthomax", "cf", "oblimin",
                     "quartimin", "geomin", "entropy", "mccammon", "infomax",
                     "tandem1", "tandem2",
                     "oblimax", "bentler", "simplimax", "target", "pst")) {
        method.fname <- paste("lav_matrix_rotate_", method, sep = "")
    } else if(method %in% c("cf-quartimax", "cf-varimax", "cf-equamax",
                            "cf-parsimax", "cf-facparsim")) {
        method.fname <- "lav_matrix_rotate_cf"
        method.args$cf.gamma <- switch(method,
            "cf-quartimax" = 0,
            "cf-varimax"   = 1 / P,
            "cf-equamax"   = M / (2 * P),
            "cf-parsimax"  = (M - 1) / (P + M - 2),
            "cf-facparsim" = 1)
    } else {
        method.fname <- method
    }

    # check if rotation method exists
    check <- try(get(method.fname), silent = TRUE)
    if(inherits(check, "try-error")) {
        stop("lavaan ERROR: unknown rotation method: ", method.fname)
    }

    # if target, check target matrix
    if(method == "target" || method == "pst") {
        target <- method.args$target
        if(is.null(target) || !is.matrix(target)) {
            stop("lavaan ERROR: target matrix is NULL, or not a matrix")
        }
        # check dimension of target/A
        if(nrow(target) != nrow(A)) {
            stop("lavaan ERROR: nrow(target) != nrow(A)")
        }
        if(ncol(target) != ncol(A)) {
            stop("lavaan ERROR: ncol(target) != ncol(A)")
        }
    }
    if(method == "pst") {
        target.mask <- method.args$target.mask
        if(is.null(target.mask) || !is.matrix(target.mask)) {
            stop("lavaan ERROR: target.mask matrix is NULL, or not a matrix")
        }
        # check dimension of target.mask/A
        if(nrow(target.mask) != nrow(A)) {
            stop("lavaan ERROR: nrow(target.mask) != nrow(A)")
        }
        if(ncol(target.mask) != ncol(A)) {
            stop("lavaan ERROR: ncol(target.mask) != ncol(A)")
        }
    }
    # if NAs, force method to be 'pst' and create target.mask
    if(method == "target" && anyNA(target)) {
        method <- "pst"
        method.fname <- "lav_matrix_rotate_pst"
        target.mask <- matrix(1, nrow = nrow(target), ncol = ncol(target))
        target.mask[ is.na(target) ] <- 0
        method.args$target.mask <- target.mask
    }


    # set orthogonal option
    if(missing(orthogonal)) {
        # the default is oblique, except for varimax, entropy and mccammon
        if(method %in% c("varimax", "entropy", "mccammon",
                         "tandem1", "tandem2") ) {
            orthogonal <- TRUE
        } else {
            orthogonal <- FALSE
        }
    } else {
        if(!orthogonal && method %in% c("varimax", "entropy", "mccammon",
                                        "tandem1", "tandem2")) {
            warning("lavaan WARNING: rotation method ", dQuote(method), 
                    " may not work with oblique rotation.")
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

    # check algorithm
    algorithm <- tolower(algorithm)
    if(algorithm %in% c("gpa", "pairwise")) {
        # nothing to do
    } else {
        stop("lavaan ERROR: algorithm must be gpa or pairwise")
    }



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
    } else {
        stop("lavaan ERROR: row.weights can be none, kaiser or cureton-mulaik")
    }
    A <- A * weights


    # 2. rotate

    # multiple random starts?
    if(rstarts > 0L) {
        REP <- sapply(seq_len(rstarts), function(rep) {
            # random start (always orthogonal)
            init.ROT <- lav_matrix_rotate_gen(M = M, orthogonal = TRUE)

            if(verbose) {
                cat("\n")
                cat("rstart = ", sprintf("%4d", rep), " start:\n")
            }


            # choose rotation algorithm
            if(algorithm == "gpa") {
                ROT <- lav_matrix_rotate_gpa(A = A, orthogonal = orthogonal,
                                             init.ROT = init.ROT,
                                             method.fname = method.fname,
                                             method.args  = method.args,
                                             warn = warn,
                                             verbose = verbose,
                                             gpa.tol = gpa.tol,
                                             max.iter = max.iter)
                info <- attr(ROT, "info"); attr(ROT, "info") <- NULL
                res <- c(info$method.value, lav_matrix_vec(ROT))

            } else if(algorithm == "pairwise") {
                ROT <- lav_matrix_rotate_pairwise(A = A,
                                             orthogonal = orthogonal,
                                             init.ROT = init.ROT,
                                             method.fname = method.fname,
                                             method.args  = method.args,
                                             warn = warn,
                                             verbose = verbose,
                                             tol = tol,
                                             max.iter = max.iter)
                info <- attr(ROT, "info"); attr(ROT, "info") <- NULL
                res <- c(info$method.value, lav_matrix_vec(ROT))

            }

            if(verbose) {
                cat("rstart = ", sprintf("%4d", rep),
                    " end; current crit = ", sprintf("%17.15f", res[1]), "\n")
            }
            res
        })
        best.idx <- which.min(REP[1,])
        ROT <- matrix(REP[-1, best.idx], nrow = M, ncol = M)
        if(keep.rep) {
            info <- list(method.value = REP[1, best.idx], REP = REP)
        } else {
            info <- list(method.value = REP[1, best.idx])
        }
    } else {
        # initial rotation matrix
        if(is.null(init.ROT)) {
            init.ROT <- diag(M)
        }

        # Gradient Projection Algorithm
        if(algorithm == "gpa") {
            ROT <- lav_matrix_rotate_gpa(A = A, orthogonal = orthogonal,
                                         init.ROT = init.ROT,
                                         method.fname = method.fname,
                                         method.args  = method.args,
                                         warn = warn,
                                         verbose = verbose,
                                         gpa.tol = gpa.tol,
                                         max.iter = max.iter)
        } else if(algorithm == "pairwise") {
            ROT <- lav_matrix_rotate_pairwise(A = A,
                                          orthogonal = orthogonal,
                                          init.ROT = init.ROT,
                                          method.fname = method.fname,
                                          method.args  = method.args,
                                          warn = warn,
                                          verbose = verbose,
                                          tol = tol,
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

    # 4.a reflect so that column sum is always positive
    if(reflect) {
        SUM <- colSums(LAMBDA)
        neg.idx <- which(SUM < 0)
        if(length(neg.idx) > 0L) {
            LAMBDA[, neg.idx] <- -1 * LAMBDA[, neg.idx, drop = FALSE]
               ROT[, neg.idx] <- -1 *    ROT[, neg.idx, drop = FALSE]
            if(!orthogonal) {
                # recompute PHI
                PHI <- crossprod(ROT)
            }
        }
    }

    # 4.b reorder the columns
    if(order.lv.by == "sumofsquares") {
        L2 <- LAMBDA * LAMBDA
        order.idx <- order(colSums(L2), decreasing = TRUE)
     } else if(order.lv.by == "index") {
        # reorder using Asparouhov & Muthen 2009 criterion (see Appendix D)
        max.loading <- apply(abs(LAMBDA), 2, max)
        # 1: per factor, number of the loadings that are at least 0.8 of the
        #    highest loading of the factor
        # 2: mean of the index numbers
        average.index <- sapply(seq_len(ncol(LAMBDA)), function(i)
                             mean(which(abs(LAMBDA[,i]) >= 0.8 * max.loading[i])))
        # order of the factors
        order.idx <- base::order(average.index)
    } else if(order.lv.by == "none") {
        order.idx <- seq_len(ncol(LAMBDA))
    } else {
        stop("lavaan ERROR: order must be index, sumofsquares or none")
    }

    # do the same in PHI
    LAMBDA <- LAMBDA[, order.idx, drop = FALSE]
    PHI    <- PHI[order.idx, order.idx, drop = FALSE]

    # 6. return results as a list
    res <- list(LAMBDA = LAMBDA, PHI = PHI, ROT = ROT, order.idx = order.idx,
                orthogonal = orthogonal, method = method,
                method.args = method.args, row.weights = row.weights)

    # add method info
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
                                  init.ROT     = NULL,   # initial rotation
                                  method.fname = NULL,   # criterion function
                                  method.args  = list(), # optional method args
                                  warn         = TRUE,
                                  verbose      = FALSE,
                                  gpa.tol      = 0.00001,
                                  max.iter     = 10000L) {

    # number of columns
    M <- ncol(A)

    # transpose of A (not needed for orthogonal)
    At <- t(A)

    # check init.ROT
    if(is.null(init.ROT)) {
        ROT <- diag(M)
    } else {
        ROT <- init.ROT
    }

    # set initial value of alpha to 1
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
        GRAD <- -1 * solve(t(init.ROT), crossprod(Gq, LAMBDA))
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
            MM <- crossprod(ROT, GRAD)
            SYMM <- (MM + t(MM))/2
            Gp <- GRAD - (ROT %*% SYMM)
        } else {
            Gp <- GRAD - t( t(ROT) * colSums(ROT * GRAD) )
        }

        # check Frobenius norm of Gp
        frob <- sqrt( sum(Gp * Gp) )

        # if verbose, print
        if(verbose) {
            cat("iter = ", sprintf("%4d", iter-1),
                " Q = ", sprintf("%9.7f", Q.current),
                " frob.log10 = ", sprintf("%10.7f", log10(frob)),
                " alpha = ", sprintf("%9.7f", alpha), "\n")
        }

        if(frob < gpa.tol) {
            converged <- TRUE
            break
        }

        # update
        alpha <- 2*alpha
        for(i in seq_len(1000)) { # make option?

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

            if(warn && i == 1000) {
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

    } # iter

    # warn if no convergence
    if(!converged && warn) {
        warning("lavaan WARNING: ",
                "GP rotation algorithm did not converge after ",
                max.iter, " iterations")
    }

    # algorithm information
    info <- list(algorithm = "gpa",
                 iter = iter - 1L,
                 converged = converged,
                 method.value = Q.current)

    attr(ROT, "info") <- info

    ROT
}


# pairwise rotation algorithm with direct line search
#
# based on Kaiser's (1959) algorithm and Jennrich and Sampson (1966) algorithm
# but to make it generic, a line search is used
#
# - orthogonal: rotate one pair of columns (=plane) at a time
# - oblique: rotate 1 factor in one pair of columns (=plane) at a time
#            note: in the oblique case, (1,2) is not the same as (2,1)
# - BUT use optimize() to find the optimal angle (for each plane)
#   (see Browne, 2001, page 130)
# - repeat until the changes in the f() criterion are below tol
#

lav_matrix_rotate_pairwise <- function(A            = NULL,   # original matrix
                                       orthogonal   = FALSE,
                                       init.ROT     = NULL,
                                       method.fname = NULL,   # crit function
                                       method.args  = list(), # method args
                                       warn         = TRUE,
                                       verbose      = FALSE,
                                       tol          = 1e-8,
                                       max.iter     = 1000L) {

    # number of columns
    M <- ncol(A)

    # initial LAMBDA + PHI
    if(is.null(init.ROT)) {
        LAMBDA <- A
        if(!orthogonal) {
            PHI <- diag(M)
        }
    } else {
        if(orthogonal) {
            LAMBDA <- A %*% init.ROT
        } else {
            LAMBDA <- t(solve(init.ROT, t(A)))
            PHI <- crossprod(init.ROT)
        }
    }

    # using the current LAMBDA, evaluate the user-specified
    # rotation criteron; return Q (the criterion) only
    Q.current <- do.call(method.fname, c(list(LAMBDA = LAMBDA),
                                         method.args, list(grad = FALSE)))

    # if verbose, print
    if(verbose) {
        cat("iter = ",    sprintf("%4d", 0),
            " Q = ",      sprintf("%13.11f", Q.current), "\n")
    }

    # plane combinations
    if(orthogonal) {
        PLANE <- utils::combn(M, 2)
    } else {
        tmp <- utils::combn(M, 2)
        PLANE <- cbind(tmp, tmp[c(2,1),,drop = FALSE])
    }


    # define objective function -- orthogonal
    objf_orth <- function(theta = 0, A = NULL, col1 = 0L, col2 = 0L) {

        # construct ROT
        ROT <- diag(M)
        ROT[col1, col1] <- base::cos(theta)
        ROT[col1, col2] <- base::sin(theta)
        ROT[col2, col1] <- -1 * base::sin(theta)
        ROT[col2, col2] <- base::cos(theta)

        # rotate
        LAMBDA <- A %*% ROT

        # evaluate criterion
        Q <- do.call(method.fname, c(list(LAMBDA = LAMBDA),
                     method.args, list(grad = FALSE)))
        Q
    }

    # define objective function -- oblique
    objf_obliq <- function(delta = 0, A = NULL, col1 = 0L, col2 = 0L,
                           phi12 = 0) {

        # construct ROT
        ROT <- diag(M)

        # gamma
        gamma2 <- 1 + (2 * delta * phi12) + (delta * delta)

        ROT[col1, col1] <- sqrt(abs(gamma2))
        ROT[col1, col2] <- -1 * delta
        ROT[col2, col1] <- 0
        ROT[col2, col2] <- 1

        # rotate
        LAMBDA <- A %*% ROT

        # evaluate criterion
        Q <- do.call(method.fname, c(list(LAMBDA = LAMBDA),
                     method.args, list(grad = FALSE)))
        Q
    }

    # start iterations
    converged <- FALSE
    Q.old <- Q.current
    for(iter in seq_len(max.iter)) {

        # rotate - one cycle
        for(pl in seq_len(ncol(PLANE))) {

            # choose plane
            col1 <- PLANE[1, pl]
            col2 <- PLANE[2, pl]

            # optimize
            if(orthogonal) {
                out <- optimize(f = objf_orth, interval = c(-pi/4, +pi/4),
                                A = LAMBDA, col1 = col1, col2 = col2,
                                maximum = FALSE, tol = .Machine$double.eps^0.25)
                # best rotation - for this plane
                theta <- out$minimum

                # construct ROT
                ROT <- diag(M)
                ROT[col1, col1] <- base::cos(theta)
                ROT[col1, col2] <- base::sin(theta)
                ROT[col2, col1] <- -1 * base::sin(theta)
                ROT[col2, col2] <- base::cos(theta)

            } else {
                phi12 <- PHI[col1, col2]
                out <- optimize(f = objf_obliq, interval = c(-1, +1),
                                A = LAMBDA, col1 = col1, col2 = col2,
                                phi12 = phi12,
                                maximum = FALSE, tol = .Machine$double.eps^0.25)

                # best rotation - for this plane
                delta <- out$minimum

                # construct ROT
                ROT <- diag(M)

                # gamma
                gamma2 <- 1 + (2 * delta * phi12) + (delta * delta)
                gamma <- sqrt(abs(gamma2))

                ROT[col1, col1] <- gamma
                ROT[col1, col2] <- -1 * delta
                ROT[col2, col1] <- 0
                ROT[col2, col2] <- 1
            }

            # rotate
            LAMBDA <- LAMBDA %*% ROT

            if(!orthogonal) {
                # rotate PHI
                PHI[col1, ] <- (1/gamma)*PHI[col1,] + (delta/gamma)*PHI[col2,]
                PHI[, col1] <- PHI[col1, ]
                PHI[col1, col1] <- 1
            }

        } # all planes

        # check for convergence
        Q.current <- do.call(method.fname, c(list(LAMBDA = LAMBDA),
                             method.args, list(grad = FALSE)))

        # absolute change in Q
        diff <- abs(Q.old - Q.current)

        # if verbose, print
        if(verbose) {
            cat("iter = ",    sprintf("%4d", iter),
                " Q = ",      sprintf("%13.11f", Q.current),
               " change = ", sprintf("%13.11f", diff), "\n")
        }

        if(diff < tol) {
            converged <- TRUE
            break
        } else {
            Q.old <- Q.current
        }

    } # iter

    # warn if no convergence
    if(!converged && warn) {
        warning("lavaan WARNING: ",
                "pairwise rotation algorithm did not converge after ",
                max.iter, " iterations")
    }

    # compute final rotation matrix
    if(orthogonal) {
        ROT <- solve(crossprod(A), crossprod(A, LAMBDA))
    } else {
        # to be compatible with GPa
        ROTt.inv <- solve(crossprod(A), crossprod(A, LAMBDA))
        ROT <- solve(t(ROTt.inv))
    }

    # algorithm information
    info <- list(algorithm = "pairwise",
                 iter = iter,
                 converged = converged,
                 method.value = Q.current)

    attr(ROT, "info") <- info

    ROT
}

