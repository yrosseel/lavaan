# efa related functions
# YR - April 2019
#
# the lav_model_efa_rotate_x() function was based on a script orginally
# written by Florian Scharf (Muenster University, Germany)

# pre-rotate: reflect and/or re-order so the 'unrotated' solution is
# in sync with the rotated solution
lav_model_efa_rotate_pre <- function(lavmodel = NULL, lavoptions = NULL) {

    if(lavmodel@nefa == 0L || lavoptions$rotation == "none") {
        return(lavmodel)
    }

    # rotation options
    ropts <- lavoptions$rotation.args

    # GLIST
    GLIST <- lavmodel@GLIST

    # for now, rotate per group (not per block)
    for(g in seq_len(lavmodel@ngroups)) {

        # select model matrices for this group
        mm.in.group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0,lavmodel@nmat))[g]
        MLIST <- GLIST[ mm.in.group]
        Hg <- diag( ncol(MLIST$lambda) )

        # reconstruct full LAMBDA (in case of dummy ov's)
        LAMBDA.g <- computeLAMBDA.LISREL(MLIST = MLIST,
                        ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
                        ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
                        ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
                        ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]],
                        remove.dummy.lv = TRUE)

        # fill in optimal rotation for each set
        for(set in seq_len(lavmodel@nefa)) {

            # which ov/lv's are involved in this set?
            ov.idx <- lavmodel@ov.efa.idx[[g]][[set]]
            lv.idx <- lavmodel@lv.efa.idx[[g]][[set]]

            # empty set?
            if(length(ov.idx) == 0L) {
                next
            }

            # just 1 factor?
            if(length(lv.idx) < 2L) {
                next
            }

            # unrotated 'A' for this set
            A <- LAMBDA.g[ov.idx, lv.idx, drop = FALSE]

            # dummy rotation matrix
            ROT <- diag(length(lv.idx))

            # reflect so that column sum is always positive
            if(ropts$reflect) {
                SUM <- colSums(A)
                neg.idx <- which(SUM < 0)
                if(length(neg.idx) > 0L) {
                     diag(ROT)[neg.idx] <- -1
                }
            }

            # reorder the columns
            if(ropts$order.lv.by == "sumofsquares") {
                L2 <- A * A
                order.idx <- base::order(colSums(L2), decreasing = TRUE)
            } else if(ropts$order.lv.by == "index") {
                # reorder using Asparouhov & Muthen 2009 criterion (see Appendix D)
                max.loading <- apply(abs(A), 2, max)
                # 1: per factor, number of the loadings that are at least 0.8 of the
                #    highest loading of the factor
                # 2: mean of the index numbers
                average.index <- sapply(seq_len(ncol(A)), function(i)
                               mean(which(abs(A[,i]) >= 0.8 * max.loading[i])))
                # order of the factors
                order.idx <- base::order(average.index)
            } else if(ropts$order.lv.by == "none") {
                order.idx <- seq_len(ncol(A))
            } else {
                stop("lavaan ERROR: order must be index, sumofsquares or none")
            }
            ROT <- ROT[, order.idx, drop = FALSE]

            # fill in optimal rotation for this set
            Hg[lv.idx, lv.idx] <- ROT

        } # set

        # 1. lambda
        MLIST$lambda <- t(solve(Hg, t(MLIST$lambda)))

        # 3. beta
        if(!is.null(MLIST$beta)) {
            MLIST$beta <- t(Hg) %*% t(solve(Hg, t(MLIST$beta)))
        }

        # 4. alpha
        if(!is.null(MLIST$alpha)) {
            MLIST$alpha <- t(Hg) %*% MLIST$alpha
        }

        # store rotated matrices in GLIST
        GLIST[ mm.in.group ] <- MLIST

    } # group

    # store 'pre-rotated' GLIST
    lavmodel@GLIST <- GLIST

    # return updated lavmodel
    lavmodel
}

# rotate solution
lav_model_efa_rotate <- function(lavmodel = NULL, x.orig = NULL,
                                 lavoptions = NULL) {

    if(lavmodel@nefa == 0L || lavoptions$rotation == "none") {
        return(lavmodel)
    }

    # extract unrotated parameters from lavmodel
    if(is.null(x.orig)) {
        x.orig <- lav_model_get_parameters(lavmodel, type = "free",
                                           extra = FALSE)
    }

    # rotate
    x.rot <- lav_model_efa_rotate_x(x = x.orig, lavmodel = lavmodel,
                                    lavoptions = lavoptions, extra = TRUE,
                                    type = "user")
    extra <- attr(x.rot, "extra"); attr(x.rot, "extra") <- NULL

    # store full rotation matrix (per group)
    lavmodel@H <- extra$H
    lavmodel@lv.order <- extra$lv.order
    lavmodel@GLIST <- extra$GLIST

    # store ceq.efa.JAC (to use for bordered information)
    if(lavoptions$se == "none" || lavoptions$se == "bootstrap" ||
       lavoptions$rotation.se == "delta") {
        JAC <- matrix(0, 0, length(x.rot))
    } else {
        JAC <- numDeriv::jacobian(func = lav_model_efa_rotate_border_x,
                                  x = x.rot, lavmodel = lavmodel,
                                  lavoptions = lavoptions, type = "user")
    }
    lavmodel@ceq.efa.JAC <- JAC

    # place rotated parameters back in lavmodel
    # lavmodel <- lav_model_set_parameters(lavmodel = lavmodel, x = x.rot)

    # return updated lavmodel
    lavmodel
}


# lower-level function, needed for numDeriv
lav_model_efa_rotate_x <- function(x, lavmodel = NULL, lavoptions = NULL,
                                   init.rot = NULL, extra = FALSE,
                                   type = "free") {

    # cat("in = \n"); print(x); cat("\n")

    # extract rotation options from lavoptions
    method <- lavoptions$rotation
    if(method == "none") {
        return(x)
    }
    ropts <- lavoptions$rotation.args

    # place parameters into model matrices
    lavmodel.orig <- lav_model_set_parameters(lavmodel, x = x)

    # GLIST
    GLIST <- lavmodel.orig@GLIST

    # H per group
    H <- vector("list", lavmodel@ngroups)
    ORDER <- vector("list", lavmodel@ngroups)

    # for now, rotate per group (not per block)
    for(g in seq_len(lavmodel@ngroups)) {

        # select model matrices for this group
        mm.in.group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0,lavmodel@nmat))[g]
        MLIST <- GLIST[ mm.in.group ]

        # general rotation matrix (all latent variables)
        H[[g]] <- Hg <- diag( ncol(MLIST$lambda) )
        lv.order <- seq_len( ncol(MLIST$lambda) )

        # reconstruct full LAMBDA (in case of dummy ov's)
        LAMBDA.g <- computeLAMBDA.LISREL(MLIST = MLIST,
                        ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
                        ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
                        ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
                        ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]],
                        remove.dummy.lv = TRUE)
        # reconstruct full THETA (in case of dummy ov's)
        THETA.g <- computeTHETA.LISREL(MLIST = MLIST,
                        ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
                        ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
                        ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
                        ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]])

        # fill in optimal rotation for each set
        for(set in seq_len(lavmodel@nefa)) {

            # which ov/lv's are involved in this set?
            ov.idx <- lavmodel@ov.efa.idx[[g]][[set]]
            lv.idx <- lavmodel@lv.efa.idx[[g]][[set]]

            # empty set?
            if(length(ov.idx) == 0L) {
                next
            }

            # just 1 factor?
            if(length(lv.idx) < 2L) {
                next
            }

            # unrotated 'A' for this set
            A <- LAMBDA.g[ov.idx, lv.idx, drop = FALSE]

            # std.ov? we use diagonal of Sigma for this set of ov's only
            if(ropts$std.ov) {
                THETA <- THETA.g[ov.idx, ov.idx, drop = FALSE]
                Sigma <- tcrossprod(A) + THETA
                this.ov.var <- diag(Sigma)
            } else {
                this.ov.var <- NULL
            }

            # init.rot?
            if(!is.null(init.rot) && lavoptions$rotation.args$jac.init.rot) {
                init.ROT <- init.rot[[g]][lv.idx, lv.idx, drop = FALSE]
                rstarts <- 0
            } else {
                init.ROT <- NULL
                rstarts <- ropts$rstarts
            }

            # rotate
            res <- lav_matrix_rotate(A           = A,
                                     orthogonal  = ropts$orthogonal,
                                     method      = method,
                                     method.args = list(
                                         geomin.epsilon = ropts$geomin.epsilon,
                                         orthomax.gamma = ropts$orthomax.gamma,
                                         cf.gamma       = ropts$orthomax.gamma,
                                         oblimin.gamma  = ropts$oblimin.gamma,
                                         target         = ropts$target,
                                         target.mask    = ropts$target.mask),
                                     init.ROT    = init.ROT,
                                     init.ROT.check = FALSE,
                                     rstarts     = rstarts,
                                     row.weights = ropts$row.weights,
                                     std.ov      = ropts$std.ov,
                                     ov.var      = this.ov.var,
                                     verbose     = ropts$verbose,
                                     warn        = ropts$warn,
                                     algorithm   = ropts$algorithm,
                                     reflect     = ropts$reflect,
                                     order.lv.by = ropts$order.lv.by,
                                     gpa.tol     = ropts$gpa.tol,
                                     tol         = ropts$tol,
                                     max.iter    = ropts$max.iter)

            # extract rotation matrix (note, in Asp & Muthen, 2009; this is H')
            # note: as of 0.6-6, order.idx has already been applied to ROT,
            #       so no need to reorder rows/columns after rotation
            H.efa <- res$ROT

            # fill in optimal rotation for this set
            Hg[lv.idx, lv.idx] <- H.efa

            # keep track of possible re-orderings
            lv.order[ lv.idx ] <- lv.idx[res$order.idx]
        } # set

        # rotate all the SEM parametersa

        # 1. lambda
        MLIST$lambda <- t(solve(Hg, t(MLIST$lambda)))
        # reorder
        #MLIST$lambda <- MLIST$lambda[, lv.order, drop = FALSE]

        # 2. psi (note: eq 22 Asp & Muthen, 2009: transpose reversed)
        MLIST$psi <- t(Hg) %*% MLIST$psi %*% Hg
        # reorder rows
        #MLIST$psi <- MLIST$psi[lv.order, , drop = FALSE]
        # reorder cols
        #MLIST$psi <- MLIST$psi[, lv.order, drop = FALSE]

        # 3. beta
        if(!is.null(MLIST$beta)) {
            MLIST$beta <- t(Hg) %*% t(solve(Hg, t(MLIST$beta)))
            # reorder rows
            #MLIST$beta <- MLIST$beta[lv.order, , drop = FALSE]
            # reorder cols
            #MLIST$beta <- MLIST$beta[, lv.order, drop = FALSE]
        }

        # 4. alpha
        if(!is.null(MLIST$alpha)) {
            MLIST$alpha <- t(Hg) %*% MLIST$alpha
            # reorder rows
            #MLIST$alpha <- MLIST$alpha[lv.order, , drop = FALSE]
        }

        # no need for rotation: nu, theta

        # store rotated matrices in GLIST
        GLIST[ mm.in.group ] <- MLIST

        # store rotation matrix + lv.order
        H[[g]] <- Hg
        ORDER[[g]] <- lv.order

    } # group

    # extract all rotated parameter estimates
    x.rot <- lav_model_get_parameters(lavmodel, GLIST = GLIST, type = type)

    # extra?
    if(extra) {
        attr(x.rot, "extra") <- list(GLIST = GLIST, H = H, lv.order = ORDER)
    }

    # cat("out = \n"); print(x.rot); cat("\n")

    # return rotated parameter estimates as a vector
    x.rot
}


# lower-level function, needed for numDeriv
lav_model_efa_rotate_border_x <- function(x, lavmodel = NULL,
                                             lavoptions = NULL) {

    # extract rotation options from lavoptions
    method <- lavoptions$rotation
    ropts <- lavoptions$rotation.args
    method.args <- list( geomin.epsilon = ropts$geomin.epsilon,
                         orthomax.gamma = ropts$orthomax.gamma,
                         cf.gamma       = ropts$orthomax.gamma,
                         oblimin.gamma  = ropts$oblimin.gamma,
                         target         = ropts$target,
                         target.mask    = ropts$target.mask )

    # place parameters into model matrices
    lavmodel <- lav_model_set_parameters(lavmodel, x = x)

    # GLIST
    GLIST <- lavmodel@GLIST

    # res
    res <- numeric(0L)

    # per group (not per block)
    for(g in seq_len(lavmodel@ngroups)) {

        # select model matrices for this group
        mm.in.group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0,lavmodel@nmat))[g]
        MLIST <- GLIST[ mm.in.group ]

        # reconstruct full LAMBDA (in case of dummy ov's)
        LAMBDA.g <- computeLAMBDA.LISREL(MLIST = MLIST,
                        ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
                        ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
                        ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
                        ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]],
                        remove.dummy.lv = TRUE)
        # reconstruct full THETA (in case of dummy ov's)
        THETA.g <- computeTHETA.LISREL(MLIST = MLIST,
                        ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
                        ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
                        ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
                        ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]])

        # for each set
        for(set in seq_len(lavmodel@nefa)) {

            # which ov/lv's are involved in this set?
            ov.idx <- lavmodel@ov.efa.idx[[g]][[set]]
            lv.idx <- lavmodel@lv.efa.idx[[g]][[set]]

            # empty set?
            if(length(ov.idx) == 0L) {
                next
            }

            # just 1 factor?
            if(length(lv.idx) < 2L) {
                next
            }

            A <- LAMBDA.g[ov.idx, lv.idx, drop = FALSE]
            P <- nrow(A); M <- ncol(A)

            # for oblique, we also need PSI
            if(!ropts$orthogonal) {
                PSI <- MLIST$psi[lv.idx, lv.idx, drop = FALSE]
            }

            # std.ov? we use diagonal of Sigma for this set of ov's only
            if(ropts$std.ov) {
                THETA <- THETA.g[ov.idx, ov.idx, drop = FALSE]
                Sigma <- tcrossprod(A) + THETA
                this.ov.var <- diag(Sigma)
            } else {
                this.ov.var <- rep(1, P)
            }

            # choose method
            method <- tolower(method)
            if(method %in% c("varimax", "quartimax", "orthomax", "cf",
                             "oblimin", "quartimin", "geomin",
                             "entropy", "mccammon", "infomax",
                             "tandem1", "tandem2",
                             "oblimax", "bentler", "simplimax", "target",
                             "pst")) {
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

            # 1. compute row weigths

            # 1.a cov -> cor?
            if(ropts$std.ov) {
                A <- A * 1/sqrt(this.ov.var)
            }

            if(ropts$row.weights == "none") {
                weights <- rep(1.0, P)
            } else if(ropts$row.weights == "kaiser") {
                weights <- lav_matrix_rotate_kaiser_weights(A)
            } else if(ropts$row.weights == "cureton-mulaik") {
                weights <- lav_matrix_rotate_cm_weights(A)
            } else {
                stop("lavaan ERROR: row.weights can be none, kaiser or cureton-mulaik")
            }
            A <- A * weights

            # evaluate rotation criterion, extract GRAD
            Q <- do.call(method.fname,
                 c(list(LAMBDA = A), method.args, list(grad = TRUE)))
            Gq <- attr(Q, "grad"); attr(Q, "grad") <- NULL

            # compute 'Z'
            Z <- crossprod(A, Gq)

            # compute constraints
            if(ropts$orthogonal) {
                # the constraint: Z == diagonal
                # or in other words, the non-diagonal elements of
                # Z - t(Z) are all zero
                tmp <- Z - t(Z)
                this.res <- lav_matrix_vech(tmp, diagonal = FALSE)
            } else {
                PSI.z <- PSI * diag(Z) # rescale rows only
                tmp <- Z - PSI.z
                out1 <- lav_matrix_vech( tmp, diagonal = FALSE)
                out2 <- lav_matrix_vechu(tmp, diagonal = FALSE)
                this.res <- c(out1, out2)
            }

            res <- c(res, this.res)

        } # set

    } # group

    # return constraint vector
    res
}

