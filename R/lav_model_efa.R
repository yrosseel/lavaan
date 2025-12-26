# efa related functions
# YR - April 2019
#
# the lav_model_efa_rotate_x() function was based on a script orginally
# written by Florian Scharf (Muenster University, Germany)

# rotate solution
lav_model_efa_rotate <- function(lavmodel = NULL, lavoptions = NULL) {

  # sanity check
  if (lavmodel@nefa == 0L || lavoptions$rotation == "none") {
    return(lavmodel)
  }

  # save warn and verbos settings
  current.warn <- lav_warn()
  current.verbose <- lav_verbose()

  # extract unrotated parameters from lavmodel
  x.orig <- lav_model_get_parameters(lavmodel, type = "free", extra = FALSE)

  # rotate, extract information from 'extra' attribute
  tmp <- lav_model_efa_rotate_x(
    x = x.orig, lavmodel = lavmodel,
    lavoptions = lavoptions, extra = TRUE
  )
  extra <- attr(tmp, "extra")
  attr(tmp, "extra") <- NULL

  out <- list(GLIST = extra$GLIST, H = extra$H, lv.order = extra$lv.order)

  # restore warn/verbose settings
  lav_warn(current.warn)
  lav_verbose(current.verbose)

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

  # new in 0.6-22 -- three options:
  # 1. rotate per group (default)
  # 2. group.equal = "loadings": use first lambda matrix only
  # 3. rotate all groups together + agreement (eg De Roover & Vermunt, 2019)
  mg.group.equal.flag <- mg.agreement.flag <- FALSE
  if (lavmodel@ngroups > 1) {
    if("loadings" %in% lavoptions$group.equal) {
      mg.group.equal.flag <- TRUE
    } else if (ropts$mg.agreement) {
      mg.agreement.flag <- TRUE
    }
  }

  # option 1 + 2 (single group or multipgroup + no agreement)
  if (!mg.agreement.flag) {
    for (g in seq_len(lavmodel@ngroups)) {
      # select model matrices for this group
      mm.in.group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0, lavmodel@nmat))[g]
      MLIST <- GLIST[mm.in.group]

      # general rotation matrix (all latent variables)
      H[[g]] <- Hg <- diag(ncol(MLIST$lambda))
      lv.order <- seq_len(ncol(MLIST$lambda))

      # reconstruct full LAMBDA (in case of dummy ov's)
      LAMBDA.g <- lav_lisrel_lambda(
        MLIST = MLIST,
        ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
        ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
        ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
        ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]],
        remove.dummy.lv = TRUE
      )
      # reconstruct full THETA (in case of dummy ov's)
      THETA.g <- lav_lisrel_theta(
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
          if (mg.group.equal.flag) {
            # compute mean THETA across groups
            theta.idx <- which(names(lavmodel@GLIST) == "theta")
            THETA <- Reduce("+", lavmodel@GLIST[theta.idx]) / lavmodel@ngroups
          } else {
            THETA <- THETA.g[ov.idx, ov.idx, drop = FALSE]
          }
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
        # set warn and verbose to ropts-values
        current.warn <- lav_warn()
        current.verbose <- lav_verbose()
        if (lav_warn(ropts$warn))
          on.exit(lav_warn(current.warn), TRUE)
        if (lav_verbose(ropts$verbose))
          on.exit(lav_verbose(current.verbose), TRUE)
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
          algorithm = ropts$algorithm,
          reflect = ropts$reflect,
          order.lv.by = ropts$order.lv.by,
          gpa.tol = ropts$gpa.tol,
          tol = ropts$tol,
          max.iter = ropts$max.iter,
  		group = g
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

      # rotate all the SEM parameters

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

      # check for group.equal
      if (mg.group.equal.flag) {
        break # only rotation once
      }
    } # group

    # if group.equal = "loadings", take care of the other groups
    if (mg.group.equal.flag) {
      MLIST1 <- MLIST # first group
      for (g in 2:lavmodel@ngroups) {
        mm.in.group <- ( seq_len(lavmodel@nmat[g]) +
                         cumsum(c(0, lavmodel@nmat))[g] )
        MLIST <- GLIST[mm.in.group]
        MLIST$lambda <- MLIST1$lambda
        MLIST$psi <- t(Hg) %*% MLIST$psi %*% Hg
        if (!is.null(MLIST$beta)) {
          MLIST$beta <- t(Hg) %*% t(solve(Hg, t(MLIST$beta)))
        }
        if (!is.null(MLIST$alpha)) {
          MLIST$alpha <- t(Hg) %*% MLIST$alpha
        }
        GLIST[mm.in.group] <- MLIST
        H[[g]] <- Hg
        ORDER[[g]] <- lv.order
      }
    }

  # option 3: rotation + agreement
  } else {
    lambdaList <- vector("list", length = lavmodel@ngroups)
    for (g in seq_len(lavmodel@ngroups)) {
      # select model matrices for this group
      mm.in.group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0, lavmodel@nmat))[g]
      MLIST <- GLIST[mm.in.group]

      # reconstruct full LAMBDA (in case of dummy ov's)
      LAMBDA.g <- lav_lisrel_lambda(
        MLIST = MLIST,
        ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
        ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
        ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
        ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]],
        remove.dummy.lv = TRUE
      )

      lambdaList[[g]] <- LAMBDA.g

      # prepare H and ORDER
      H[[g]] <- diag(ncol(LAMBDA.g))
      ORDER[[g]] <- seq_len(ncol(LAMBDA.g))
    } # group

    # fill in optimal rotation for each set
    for (set in seq_len(lavmodel@nefa)) {
      # which ov/lv's are involved in this set?
      ov.idx <- lavmodel@ov.efa.idx[[1]][[set]] # assumed equal across groups!
      lv.idx <- lavmodel@lv.efa.idx[[1]][[set]] # assumed equal across groups!

      # empty set?
      if (length(ov.idx) == 0L) {
        next
      }

      # just 1 factor?
#       if (length(lv.idx) < 2L) {
#         # new in 0.6-18: reflect if needed
#         if (lavoptions$rotation.args$reflect) {
#           tmp <- LAMBDA.g[ov.idx, lv.idx, drop = TRUE]
#           if (sum(tmp) < 0) {
#             MLIST$lambda[ov.idx, 1] <- -1 * MLIST$lambda[ov.idx, 1]
#           }
#         }
#         next
#       }

      # unrotated 'A' for this set
      Alist <- lapply(lambdaList,
        function (x) { x[ov.idx, lv.idx, drop = FALSE] })
      # std.ov? we use diagonal of Sigma for this set of ov's only
#      if (ropts$std.ov) {
#        if (mg.group.equal.flag) {
#          # compute mean THETA across groups
#          theta.idx <- which(names(lavmodel@GLIST) == "theta")
#          THETA <- Reduce("+", lavmodel@GLIST[theta.idx]) / lavmodel@ngroups
#        } else {
#          THETA <- THETA.g[ov.idx, ov.idx, drop = FALSE]
#        }
#        Sigma <- tcrossprod(A) + THETA
#        this.ov.var <- diag(Sigma)
#      } else {
        this.ov.var <- NULL
#      }

      # init.rot?
      if (!is.null(init.rot) && lavoptions$rotation.args$jac.init.rot) {
        init.rotList <- vector("list", length = lavmodel@ngroups)
        init.rotList[[g]] <- init.rot[[g]][lv.idx, lv.idx, drop = FALSE]
        rstarts <- 0
      } else {
        init.rotList <- NULL
        rstarts <- ropts$rstarts
      }

      # set warn and verbose to ropts-values
      current.warn <- lav_warn()
      current.verbose <- lav_verbose()
      if (lav_warn(ropts$warn))
        on.exit(lav_warn(current.warn), TRUE)
      if (lav_verbose(ropts$verbose))
        on.exit(lav_verbose(current.verbose), TRUE)
      if (lav_verbose()) {
        cat("\n")
      }

      # rotate this set
      res <- lav_matrix_rotate_mg(
        Alist = Alist,
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
        init.rotList = init.rotList,
        init.ROT.check = FALSE,
        rstarts = rstarts,
        #row.weights = ropts$row.weights,
        #std.ov = ropts$std.ov,
        #ov.var = this.ov.var,
        #algorithm = ropts$algorithm,
        reflect = ropts$reflect,
        order.lv.by = ropts$order.lv.by,
        #gpa.tol = ropts$gpa.tol,
        tol = ropts$tol,
        max.iter = ropts$max.iter
      )
      for (g in seq_len(lavmodel@ngroups)) {
        # fill in optimal rotation for this set, for this group
        H[[g]][lv.idx, lv.idx] <- res$rotList[[g]]

        # keep track of possible re-orderings
        ORDER[[g]][lv.idx] <- lv.idx[res$orderList[[g]]]
      }
    } # set

    # rotate all the SEM parameters
    for (g in seq_len(lavmodel@ngroups)) {
      mm.in.group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0, lavmodel@nmat))[g]
      MLIST <- GLIST[mm.in.group]
      Hg <- H[[g]]

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

      GLIST[mm.in.group] <- MLIST
    }

    # rescale lambda again, so that sum-of-squares per columns is the same
    # for all groups
    lambdaList <- GLIST[names(GLIST) == "lambda"]
    psiList    <- GLIST[names(GLIST) == "psi"]

    ss_list <- lapply(lambdaList, function(x) colSums(x^2))
    ss_ave <- Reduce("+", ss_list) / length(ss_list)
    scale_factor <- vector("list", length = ncol(lambdaList[[1]]))
    for(g in seq_len(lavmodel@ngroups)) {
      scale_factor[[g]] <- (sqrt(ss_list[[g]]) / sqrt(ss_ave))
      lambdaList[[g]] <- t( t(lambdaList[[g]]) / scale_factor[[g]] )
    }

    # rescale psiList accordingly
    for(g in seq_len(lavmodel@ngroups)) {
      psiList[[g]] <- ( diag(scale_factor[[g]]) %*% psiList[[g]] %*%
                        diag(scale_factor[[g]]) )
    }

    # put back in GLIST
    GLIST[names(GLIST) == "lambda"] <- lambdaList
    GLIST[names(GLIST) == "psi"]    <- psiList

  } # option 3

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

    # group-specific method.args
	this.method.args <- method.args

	# set group-specific target/target.mask (if needed)
    # if target, check target matrix
    if (method == "target.strict" || method == "pst") {
      target <- method.args$target
      if (is.list(target)) {
        this.method.args$target <- target[[g]]
      }
    }
    if (method == "pst") {
      target.mask <- method.args$target.mask
      if (is.list(target.mask)) {
        this.method.args$target.mask <- target.mask[[g]]
      }
    }

    # select model matrices for this group
    mm.in.group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0, lavmodel@nmat))[g]
    MLIST <- GLIST[mm.in.group]

    # reconstruct full LAMBDA (in case of dummy ov's)
    LAMBDA.g <- lav_lisrel_lambda(
      MLIST = MLIST,
      ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
      ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
      ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
      ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]],
      remove.dummy.lv = TRUE
    )
    # reconstruct full THETA (in case of dummy ov's)
    THETA.g <- lav_lisrel_theta(
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
        this.method.args$cf.gamma <- switch(method,
          "cf-quartimax" = 0,
          "cf-varimax"   = 1 / P,
          "cf-equamax"   = M / (2 * P),
          "cf-parsimax"  = (M - 1) / (P + M - 2),
          "cf-facparsim" = 1
        )
      } else if (method == "target.strict") {
        method.fname <- "lav_matrix_rotate_target"
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
        c(list(LAMBDA = A), this.method.args, list(grad = TRUE))
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
