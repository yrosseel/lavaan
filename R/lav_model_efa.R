# efa related functions
# YR - April 2019
#
# the lav_model_efa_rotate_x() function was based on a script orginally
# written by Florian Scharf (Muenster University, Germany)

# rotate solution
lav_model_efa_rotate <- function(lavmodel = NULL, lavoptions = NULL) {
  if (lavmodel@nefa == 0L || lavoptions$rotation == "none") {
    return(lavmodel)
  }

  # extract unrotated parameters from lavmodel
  x.orig <- lav_model_get_parameters(lavmodel, type = "free", extra = FALSE)

  # rotate, extract information from 'extra' attribute
  tmp <- lav_model_efa_rotate_x(
    x = x.orig, lavmodel = lavmodel,
    lavoptions = lavoptions, extra = TRUE
  )
  extra <- attr(tmp, "extra")
  attr(tmp, "extra") <- NULL

  # store full rotation matrix (per group)
  # lavmodel@H <- extra$H
  # lavmodel@lv.order <- extra$lv.order
  # lavmodel@GLIST <- extra$GLIST

  # return updated lavmodel
  # lavmodel

  out <- list(GLIST = extra$GLIST, H = extra$H, lv.order = extra$lv.order)

  out
}


# lower-level function, needed for numDeriv
lav_model_efa_rotate_x <- function(x, lavmodel = NULL, lavoptions = NULL,
                                   init.rot = NULL, extra = FALSE,
                                   type = "free") {
  # extract rotation options from lavoptions
  method <- lavoptions$rotation
  if (method == "none") {
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
  for (g in seq_len(lavmodel@ngroups)) {
    # select model matrices for this group
    mm.in.group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0, lavmodel@nmat))[g]
    MLIST <- GLIST[mm.in.group]

    # general rotation matrix (all latent variables)
    H[[g]] <- Hg <- diag(ncol(MLIST$lambda))
    lv.order <- seq_len(ncol(MLIST$lambda))

    # reconstruct full LAMBDA (in case of dummy ov's)
    LAMBDA.g <- computeLAMBDA.LISREL(
      MLIST = MLIST,
      ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
      ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
      ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
      ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]],
      remove.dummy.lv = TRUE
    )
    # reconstruct full THETA (in case of dummy ov's)
    THETA.g <- computeTHETA.LISREL(
      MLIST = MLIST,
      ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
      ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
      ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
      ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]]
    )

    # fill in optimal rotation for each set
    for (set in seq_len(lavmodel@nefa)) {
      # which ov/lv's are involved in this set?
      ov.idx <- lavmodel@ov.efa.idx[[g]][[set]]
      lv.idx <- lavmodel@lv.efa.idx[[g]][[set]]

      # empty set?
      if (length(ov.idx) == 0L) {
        next
      }

      # just 1 factor?
      if (length(lv.idx) < 2L) {
        # new in 0.6-18: reflect if needed
        if (lavoptions$rotation.args$reflect) {
          tmp <- LAMBDA.g[ov.idx, lv.idx, drop = TRUE]
          if (sum(tmp) < 0) {
            MLIST$lambda[ov.idx, 1] <- -1 * MLIST$lambda[ov.idx, 1]
          }
        }
        next
      }

      # unrotated 'A' for this set
      A <- LAMBDA.g[ov.idx, lv.idx, drop = FALSE]

      # std.ov? we use diagonal of Sigma for this set of ov's only
      if (ropts$std.ov) {
        THETA <- THETA.g[ov.idx, ov.idx, drop = FALSE]
        Sigma <- tcrossprod(A) + THETA
        this.ov.var <- diag(Sigma)
      } else {
        this.ov.var <- NULL
      }

      # init.rot?
      if (!is.null(init.rot) && lavoptions$rotation.args$jac.init.rot) {
        init.ROT <- init.rot[[g]][lv.idx, lv.idx, drop = FALSE]
        rstarts <- 0
      } else {
        init.ROT <- NULL
        rstarts <- ropts$rstarts
      }

      # rotate this set
      res <- lav_matrix_rotate(
        A = A,
        orthogonal = ropts$orthogonal,
        method = method,
        method.args = list(
          geomin.epsilon = ropts$geomin.epsilon,
          orthomax.gamma = ropts$orthomax.gamma,
          cf.gamma       = ropts$orthomax.gamma,
          oblimin.gamma  = ropts$oblimin.gamma,
          promax.kappa   = ropts$promax.kappa,
          target         = ropts$target,
          target.mask    = ropts$target.mask
        ),
        init.ROT = init.ROT,
        init.ROT.check = FALSE,
        rstarts = rstarts,
        row.weights = ropts$row.weights,
        std.ov = ropts$std.ov,
        ov.var = this.ov.var,
        verbose = ropts$verbose,
        warn = ropts$warn,
        algorithm = ropts$algorithm,
        reflect = ropts$reflect,
        order.lv.by = ropts$order.lv.by,
        gpa.tol = ropts$gpa.tol,
        tol = ropts$tol,
        max.iter = ropts$max.iter
      )

      # extract rotation matrix (note, in Asp & Muthen, 2009; this is H')
      # note: as of 0.6-6, order.idx has already been applied to ROT,
      #       so no need to reorder rows/columns after rotation
      H.efa <- res$ROT

      # fill in optimal rotation for this set
      Hg[lv.idx, lv.idx] <- H.efa

      # keep track of possible re-orderings
      lv.order[lv.idx] <- lv.idx[res$order.idx]
    } # set

    # rotate all the SEM parametersa

    # 1. lambda
    MLIST$lambda <- t(solve(Hg, t(MLIST$lambda)))

    # 2. psi (note: eq 22 Asp & Muthen, 2009: transpose reversed)
    MLIST$psi <- t(Hg) %*% MLIST$psi %*% Hg

    # 3. beta
    if (!is.null(MLIST$beta)) {
      MLIST$beta <- t(Hg) %*% t(solve(Hg, t(MLIST$beta)))
    }

    # 4. alpha
    if (!is.null(MLIST$alpha)) {
      MLIST$alpha <- t(Hg) %*% MLIST$alpha
    }

    # no need for rotation: nu, theta

    # store rotated matrices in GLIST
    GLIST[mm.in.group] <- MLIST

    # store rotation matrix + lv.order
    H[[g]] <- Hg
    ORDER[[g]] <- lv.order
  } # group

  # extract all rotated parameter estimates
  x.rot <- lav_model_get_parameters(lavmodel, GLIST = GLIST, type = type)

  # extra?
  if (extra) {
    attr(x.rot, "extra") <- list(GLIST = GLIST, H = H, lv.order = ORDER)
  }

  # return rotated parameter estimates as a vector
  x.rot
}


# lower-level function, needed for numDeriv
lav_model_efa_rotate_border_x <- function(x, lavmodel = NULL,
                                          lavoptions = NULL,
                                          lavpartable = NULL) {
  # extract rotation options from lavoptions
  method <- lavoptions$rotation
  ropts <- lavoptions$rotation.args
  method.args <- list(
    geomin.epsilon = ropts$geomin.epsilon,
    orthomax.gamma = ropts$orthomax.gamma,
    cf.gamma = ropts$orthomax.gamma,
    oblimin.gamma = ropts$oblimin.gamma,
    promax.kappa = ropts$oblimin.kappa,
    target = ropts$target,
    target.mask = ropts$target.mask
  )

  # place parameters into model matrices
  lavmodel <- lav_model_set_parameters(lavmodel, x = x)

  # GLIST
  GLIST <- lavmodel@GLIST

  # res
  res <- numeric(0L)

  # per group (not per block)
  for (g in seq_len(lavmodel@ngroups)) {
    # select model matrices for this group
    mm.in.group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0, lavmodel@nmat))[g]
    MLIST <- GLIST[mm.in.group]

    # reconstruct full LAMBDA (in case of dummy ov's)
    LAMBDA.g <- computeLAMBDA.LISREL(
      MLIST = MLIST,
      ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
      ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
      ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
      ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]],
      remove.dummy.lv = TRUE
    )
    # reconstruct full THETA (in case of dummy ov's)
    THETA.g <- computeTHETA.LISREL(
      MLIST = MLIST,
      ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
      ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
      ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
      ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]]
    )

    # setnames
    set.names <- lav_partable_efa_values(lavpartable)

    # for each set
    for (set in seq_len(lavmodel@nefa)) {
      # check if we have any user=7 elements in this set
      # if not, skip constraints
      ind.idx <- which(lavpartable$op == "=~" &
        lavpartable$group == g &
        lavpartable$efa == set.names[set])
      if (!any(lavpartable$user[ind.idx] == 7L)) {
        next
      }

      # which ov/lv's are involved in this set?
      ov.idx <- lavmodel@ov.efa.idx[[g]][[set]]
      lv.idx <- lavmodel@lv.efa.idx[[g]][[set]]

      # empty set?
      if (length(ov.idx) == 0L) {
        next
      }

      # just 1 factor?
      if (length(lv.idx) < 2L) {
        next
      }

      A <- LAMBDA.g[ov.idx, lv.idx, drop = FALSE]
      P <- nrow(A)
      M <- ncol(A)

      # for oblique, we also need PSI
      if (!ropts$orthogonal) {
        PSI <- MLIST$psi[lv.idx, lv.idx, drop = FALSE]
      }

      # std.ov? we use diagonal of Sigma for this set of ov's only
      if (ropts$std.ov) {
        THETA <- THETA.g[ov.idx, ov.idx, drop = FALSE]
        Sigma <- tcrossprod(A) + THETA
        this.ov.var <- diag(Sigma)
      } else {
        this.ov.var <- rep(1, P)
      }

      # choose method
      method <- tolower(method)
      if (method %in% c(
        "cf-quartimax", "cf-varimax", "cf-equamax",
        "cf-parsimax", "cf-facparsim"
      )) {
        method.fname <- "lav_matrix_rotate_cf"
        method.args$cf.gamma <- switch(method,
          "cf-quartimax" = 0,
          "cf-varimax"   = 1 / P,
          "cf-equamax"   = M / (2 * P),
          "cf-parsimax"  = (M - 1) / (P + M - 2),
          "cf-facparsim" = 1
        )
      } else {
        method.fname <- paste("lav_matrix_rotate_", method, sep = "")
      }

      # check if rotation method exists
      check <- try(get(method.fname), silent = TRUE)
      if (inherits(check, "try-error")) {
        lav_msg_stop(gettextf("unknown rotation method: %s", method.fname))
      }

      # 1. compute row weigths

      # 1.a cov -> cor?
      if (ropts$std.ov) {
        A <- A * 1 / sqrt(this.ov.var)
      }

      if (ropts$row.weights == "none") {
        weights <- rep(1.0, P)
      } else if (ropts$row.weights == "kaiser") {
        weights <- lav_matrix_rotate_kaiser_weights(A)
      } else if (ropts$row.weights == "cureton-mulaik") {
        weights <- lav_matrix_rotate_cm_weights(A)
      } else {
        lav_msg_stop(gettextf("row.weights can be",
          lav_msg_view(c("none", "kaiser", "cureton-mulaik"), "or")))
      }
      A <- A * weights

      # evaluate rotation criterion, extract GRAD
      Q <- do.call(
        method.fname,
        c(list(LAMBDA = A), method.args, list(grad = TRUE))
      )
      Gq <- attr(Q, "grad")
      attr(Q, "grad") <- NULL

      # compute 'Z'
      Z <- crossprod(A, Gq)

      # compute constraints
      if (ropts$orthogonal) {
        # the constraint: Z == diagonal
        # or in other words, the non-diagonal elements of
        # Z - t(Z) are all zero
        tmp <- Z - t(Z)
        this.res <- lav_matrix_vech(tmp, diagonal = FALSE)
      } else {
        PSI.z <- PSI * diag(Z) # rescale rows only
        tmp <- Z - PSI.z
        out1 <- lav_matrix_vech(tmp, diagonal = FALSE)
        out2 <- lav_matrix_vechu(tmp, diagonal = FALSE)
        this.res <- c(out1, out2)
      }

      res <- c(res, this.res)
    } # set
  } # group

  # return constraint vector
  res
}
