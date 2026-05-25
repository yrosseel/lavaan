# efa related functions
# YR - April 2019
#
# the lav_model_efa_rotate_x() function was based on a script originally
# written by Florian Scharf (Muenster University, Germany)

# rotate solution
lav_model_efa_rotate <- function(lavmodel = NULL, lavoptions = NULL) {

  # sanity check
  if (lavmodel@nefa == 0L || lavoptions$rotation == "none") {
    return(lavmodel)
  }

  # save warn and verbose settings
  current_warn <- lav_warn()
  current_verbose <- lav_verbose()

  # extract unrotated parameters from lavmodel
  x_orig <- lav_model_get_parameters(lavmodel, type = "free", extra = FALSE)

  # rotate, extract information from 'extra' attribute
  tmp <- lav_model_efa_rotate_x(
    x = x_orig, lavmodel = lavmodel,
    lavoptions = lavoptions, extra = TRUE
  )
  extra <- attr(tmp, "extra")
  attr(tmp, "extra") <- NULL

  out <- list(GLIST = extra$GLIST, H = extra$H, lv.order = extra$lv.order)

  # restore warn/verbose settings
  lav_warn(current_warn)
  lav_verbose(current_verbose)

  out
}


# lower-level function, needed for numDeriv
lav_model_efa_rotate_x <- function(x, lavmodel = NULL, lavoptions = NULL,
                                   init_rot = NULL, extra = FALSE,
                                   type = "free") {
  # extract rotation options from lavoptions
  method <- lavoptions$rotation
  if (method == "none") {
    return(x)
  }
  ropts <- lavoptions$rotation.args

  # place parameters into model matrices
  lavmodel_orig <- lav_model_set_parameters(lavmodel, x = x)

  # GLIST
  glist <- lavmodel_orig@GLIST

  # H per group
  h <- vector("list", lavmodel@ngroups)
  order_1 <- vector("list", lavmodel@ngroups)

  # new in 0.6-22 -- three options:
  # 1. rotate per group (default)
  # 2. group.equal = "loadings": use first lambda matrix only
  # 3. rotate all groups together + agreement (eg De Roover & Vermunt, 2019)
  mg_group_equal_flag <- mg_agreement_flag <- FALSE
  if (lavmodel@ngroups > 1) {
    if ("loadings" %in% lavoptions$group.equal) {
      mg_group_equal_flag <- TRUE
    } else if (ropts$mg_agreement) {
      mg_agreement_flag <- TRUE
    }
  }

  # option 1 + 2 (single group or multipgroup + no agreement)
  if (!mg_agreement_flag) {
    for (g in seq_len(lavmodel@ngroups)) {
      # select model matrices for this group
      mm_in_group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0, lavmodel@nmat))[g]
      mlist <- glist[mm_in_group]

      # general rotation matrix (all latent variables)
      h[[g]] <- hg <- diag(ncol(mlist$lambda))
      lv_order <- seq_len(ncol(mlist$lambda))

      # reconstruct full mm_lambda (in case of dummy ov's)
      lambda_g <- lav_lisrel_lambda(
        mlist = mlist,
        ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[g]],
        ov_x_dummy_ov_idx = lavmodel@ov.x.dummy.ov.idx[[g]],
        ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[g]],
        ov_x_dummy_lv_idx = lavmodel@ov.x.dummy.lv.idx[[g]],
        remove_dummy_lv = TRUE
      )
      # reconstruct full THETA (in case of dummy ov's)
      theta_g <- lav_lisrel_theta(
        mlist = mlist,
        ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[g]],
        ov_x_dummy_ov_idx = lavmodel@ov.x.dummy.ov.idx[[g]],
        ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[g]],
        ov_x_dummy_lv_idx = lavmodel@ov.x.dummy.lv.idx[[g]]
      )

      # fill in optimal rotation for each set
      for (set in seq_len(lavmodel@nefa)) {
        # which ov/lv's are involved in this set?
        ov_idx <- lavmodel@ov.efa.idx[[g]][[set]]
        lv_idx <- lavmodel@lv.efa.idx[[g]][[set]]

        # empty set?
        if (length(ov_idx) == 0L) {
          next
        }

        # just 1 factor?
        if (length(lv_idx) < 2L) {
          # new in 0.6-18: reflect if needed
          if (lavoptions$rotation.args$reflect) {
            tmp <- lambda_g[ov_idx, lv_idx, drop = TRUE]
            if (sum(tmp) < 0) {
              mlist$lambda[ov_idx, 1] <- -1 * mlist$lambda[ov_idx, 1]
            }
          }
          next
        }

        # unrotated 'A' for this set
        a <- lambda_g[ov_idx, lv_idx, drop = FALSE]

        # std.ov? we use diagonal of Sigma for this set of ov's only
        if (ropts$std_ov) {
          if (mg_group_equal_flag) {
            # compute mean THETA across groups
            theta_idx <- which(names(lavmodel@GLIST) == "theta")
            mm_theta <- Reduce("+", lavmodel@GLIST[theta_idx]) /
                                                    lavmodel@ngroups
          } else {
            mm_theta <- theta_g[ov_idx, ov_idx, drop = FALSE]
          }
          sigma_1 <- tcrossprod(a) + mm_theta
          this_ov_var <- diag(sigma_1)
        } else {
          this_ov_var <- NULL
        }

        # init.rot?
        if (!is.null(init_rot) && lavoptions$rotation.args$jac_init_rot) {
          init_rot_1 <- init_rot[[g]][lv_idx, lv_idx, drop = FALSE]
          rstarts <- 0
        } else {
          init_rot_1 <- NULL
          rstarts <- ropts$rstarts
        }
        # set warn and verbose to ropts-values
        current_warn <- lav_warn()
        current_verbose <- lav_verbose()
        if (lav_warn(ropts$warn))
          on.exit(lav_warn(current_warn), TRUE)
        if (lav_verbose(ropts$verbose))
          on.exit(lav_verbose(current_verbose), TRUE)
        # rotate this set
        res <- lav_mat_rotate(
          a = a,
          orthogonal = ropts$orthogonal,
          method = method,
          method_args = list(
            geomin_epsilon = ropts$geomin_epsilon,
            orthomax_gamma = ropts$orthomax_gamma,
            cf_gamma       = ropts$orthomax_gamma,
            oblimin_gamma  = ropts$oblimin_gamma,
            promax_kappa   = ropts$promax_kappa,
            target         = ropts$target,
            target_mask    = ropts$target_mask
          ),
          init_rot = init_rot_1,
          init_rot_check = FALSE,
          rstarts = rstarts,
          row_weights = ropts$row_weights,
          std_ov = ropts$std_ov,
          ov_var = this_ov_var,
          algorithm = ropts$algorithm,
          reflect = ropts$reflect,
          order_lv_by = ropts$order_lv_by,
          gpa_tol = ropts$gpa_tol,
          tol = ropts$tol,
          max_iter = ropts$max_iter,
          group = g
        )

        # extract rotation matrix (note, in Asp & Muthen, 2009; this is H')
        # note: as of 0.6-6, order.idx has already been applied to ROT,
        #       so no need to reorder rows/columns after rotation
        h_efa <- res$ROT

        # fill in optimal rotation for this set
        hg[lv_idx, lv_idx] <- h_efa

        # keep track of possible re-orderings
        lv_order[lv_idx] <- lv_idx[res$order.idx]
      } # set

      # rotate all the SEM parameters

      # 1. lambda
      mlist$lambda <- t(solve(hg, t(mlist$lambda)))

      # 2. psi (note: eq 22 Asp & Muthen, 2009: transpose reversed)
      mlist$psi <- t(hg) %*% mlist$psi %*% hg

      # 3. beta
      if (!is.null(mlist$beta)) {
        mlist$beta <- t(hg) %*% t(solve(hg, t(mlist$beta)))
      }

      # 4. alpha
      if (!is.null(mlist$alpha)) {
        mlist$alpha <- t(hg) %*% mlist$alpha
      }

      # no need for rotation: nu, theta

      # store rotated matrices in GLIST
      glist[mm_in_group] <- mlist

      # store rotation matrix + lv.order
      h[[g]] <- hg
      order_1[[g]] <- lv_order

      # check for group.equal
      if (mg_group_equal_flag) {
        break # only rotation once
      }
    } # group

    # if group.equal = "loadings", take care of the other groups
    if (mg_group_equal_flag) {
      mlist1 <- mlist # first group
      for (g in 2:lavmodel@ngroups) {
        mm_in_group <- (seq_len(lavmodel@nmat[g]) +
                         cumsum(c(0, lavmodel@nmat))[g])
        mlist <- glist[mm_in_group]
        mlist$lambda <- mlist1$lambda
        mlist$psi <- t(hg) %*% mlist$psi %*% hg
        if (!is.null(mlist$beta)) {
          mlist$beta <- t(hg) %*% t(solve(hg, t(mlist$beta)))
        }
        if (!is.null(mlist$alpha)) {
          mlist$alpha <- t(hg) %*% mlist$alpha
        }
        glist[mm_in_group] <- mlist
        h[[g]] <- hg
        order_1[[g]] <- lv_order
      }
    }

  # option 3: rotation + agreement
  } else {
    lambda_list <- vector("list", length = lavmodel@ngroups)
    for (g in seq_len(lavmodel@ngroups)) {
      # select model matrices for this group
      mm_in_group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0, lavmodel@nmat))[g]
      mlist <- glist[mm_in_group]

      # reconstruct full mm_lambda (in case of dummy ov's)
      lambda_g <- lav_lisrel_lambda(
        mlist = mlist,
        ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[g]],
        ov_x_dummy_ov_idx = lavmodel@ov.x.dummy.ov.idx[[g]],
        ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[g]],
        ov_x_dummy_lv_idx = lavmodel@ov.x.dummy.lv.idx[[g]],
        remove_dummy_lv = TRUE
      )

      lambda_list[[g]] <- lambda_g

      # prepare H and ORDER
      h[[g]] <- diag(ncol(lambda_g))
      order_1[[g]] <- seq_len(ncol(lambda_g))
    } # group

    # fill in optimal rotation for each set
    for (set in seq_len(lavmodel@nefa)) {
      # which ov/lv's are involved in this set?
      ov_idx <- lavmodel@ov.efa.idx[[1]][[set]] # assumed equal across groups!
      lv_idx <- lavmodel@lv.efa.idx[[1]][[set]] # assumed equal across groups!

      # empty set?
      if (length(ov_idx) == 0L) {
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
      alist_1 <- lapply(lambda_list,
        function(x) x[ov_idx, lv_idx, drop = FALSE])
      # std.ov? we use diagonal of Sigma for this set of ov's only
#      if (ropts$std_ov) {
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
        this_ov_var <- NULL
#      }

      # init.rot?
      if (!is.null(init_rot) && lavoptions$rotation.args$jac_init_rot) {
        init_rot_list <- vector("list", length = lavmodel@ngroups)
        init_rot_list[[g]] <- init_rot[[g]][lv_idx, lv_idx, drop = FALSE]
        rstarts <- 0
      } else {
        init_rot_list <- NULL
        rstarts <- ropts$rstarts
      }

      # set warn and verbose to ropts-values
      current_warn <- lav_warn()
      current_verbose <- lav_verbose()
      if (lav_warn(ropts$warn))
        on.exit(lav_warn(current_warn), TRUE)
      if (lav_verbose(ropts$verbose))
        on.exit(lav_verbose(current_verbose), TRUE)
      if (lav_verbose()) {
        cat("\n")
      }

      # rotate this set
      res <- lav_mat_rotate_mg(
        a_list = alist_1,
        orthogonal = ropts$orthogonal,
        method = method,
        method_args = list(
          geomin_epsilon = ropts$geomin_epsilon,
          orthomax_gamma = ropts$orthomax_gamma,
          cf_gamma       = ropts$orthomax_gamma,
          oblimin_gamma  = ropts$oblimin_gamma,
          promax_kappa   = ropts$promax_kappa,
          target         = ropts$target,
          target_mask    = ropts$target_mask
        ),
        init_rot_list = init_rot_list,
        init_rot_check = FALSE,
        rstarts = rstarts,
        #row_weights = ropts$row_weights,
        #std.ov = ropts$std_ov,
        #ov.var = this.ov.var,
        mg_algorithm = ropts$mg_agreement_algorithm,
        mg_agreement_method = ropts$mg_agreement_method,
        mg_agreement_weight = ropts$mg_agreement_weight,
        reflect = ropts$reflect,
        order_lv_by = ropts$order_lv_by,
        #gpa.tol = ropts$gpa_tol,
        tol = ropts$tol,
        max_iter = ropts$max_iter
      )
      for (g in seq_len(lavmodel@ngroups)) {
        # fill in optimal rotation for this set, for this group
        h[[g]][lv_idx, lv_idx] <- res$rotList[[g]]

        # keep track of possible re-orderings
        order_1[[g]][lv_idx] <- lv_idx[res$orderList[[g]]]
      }
    } # set

    # rotate all the SEM parameters
    for (g in seq_len(lavmodel@ngroups)) {
      mm_in_group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0, lavmodel@nmat))[g]
      mlist <- glist[mm_in_group]
      hg <- h[[g]]

      # 1. lambda
      mlist$lambda <- t(solve(hg, t(mlist$lambda)))

      # 2. psi (note: eq 22 Asp & Muthen, 2009: transpose reversed)
      mlist$psi <- t(hg) %*% mlist$psi %*% hg

      # 3. beta
      if (!is.null(mlist$beta)) {
        mlist$beta <- t(hg) %*% t(solve(hg, t(mlist$beta)))
      }

      # 4. alpha
      if (!is.null(mlist$alpha)) {
        mlist$alpha <- t(hg) %*% mlist$alpha
      }
      # no need for rotation: nu, theta

      glist[mm_in_group] <- mlist
    }

    # rescale lambda again, so that sum-of-squares per columns is the same
    # for all groups
    lambda_list <- glist[names(glist) == "lambda"]
    psi_list    <- glist[names(glist) == "psi"]

    ss_list <- lapply(lambda_list, function(x) colSums(x^2))
    ss_ave <- Reduce("+", ss_list) / length(ss_list)
    scale_factor <- vector("list", length = ncol(lambda_list[[1]]))
    for (g in seq_len(lavmodel@ngroups)) {
      scale_factor[[g]] <- (sqrt(ss_list[[g]]) / sqrt(ss_ave))
      lambda_list[[g]] <- t(t(lambda_list[[g]]) / scale_factor[[g]])
    }

    # rescale psiList accordingly
    for (g in seq_len(lavmodel@ngroups)) {
      psi_list[[g]] <- (diag(scale_factor[[g]]) %*% psi_list[[g]] %*%
                        diag(scale_factor[[g]]))
    }

    # put back in GLIST
    glist[names(glist) == "lambda"] <- lambda_list
    glist[names(glist) == "psi"]    <- psi_list

  } # option 3

  # extract all rotated parameter estimates
  x_rot <- lav_model_get_parameters(lavmodel, GLIST = glist, type = type)

  # extra?
  if (extra) {
    attr(x_rot, "extra") <- list(GLIST = glist, H = h, lv.order = order_1)
  }

  # return rotated parameter estimates as a vector
  x_rot
}


# lower-level function, needed for numDeriv
lav_model_efa_rotate_border_x <- function(x, lavmodel = NULL,
                                          lavoptions = NULL,
                                          lavpartable = NULL) {
  # extract rotation options from lavoptions
  method <- lavoptions$rotation
  ropts <- lavoptions$rotation.args
  method_args <- list(
    geomin_epsilon = ropts$geomin_epsilon,
    orthomax_gamma = ropts$orthomax_gamma,
    cf_gamma = ropts$orthomax_gamma,
    oblimin_gamma = ropts$oblimin_gamma,
    promax_kappa = ropts$oblimin_kappa,
    target = ropts$target,
    target_mask = ropts$target_mask
  )

  # place parameters into model matrices
  lavmodel <- lav_model_set_parameters(lavmodel, x = x)

  # GLIST
  glist <- lavmodel@GLIST

  # res
  res <- numeric(0L)

  # per group (not per block)
  for (g in seq_len(lavmodel@ngroups)) {

    # group-specific method.args
  this_method_args <- method_args

  # set group-specific target/target_mask (if needed)
    # if target, check target matrix
    if (method == "target.strict" || method == "pst") {
      target <- method_args$target
      if (is.list(target)) {
        this_method_args$target <- target[[g]]
      }
    }
    if (method == "pst") {
      target_mask <- method_args$target_mask
      if (is.list(target_mask)) {
        this_method_args$target_mask <- target_mask[[g]]
      }
    }

    # select model matrices for this group
    mm_in_group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0, lavmodel@nmat))[g]
    mlist <- glist[mm_in_group]

    # reconstruct full mm_lambda (in case of dummy ov's)
    lambda_g <- lav_lisrel_lambda(
      mlist = mlist,
      ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[g]],
      ov_x_dummy_ov_idx = lavmodel@ov.x.dummy.ov.idx[[g]],
      ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[g]],
      ov_x_dummy_lv_idx = lavmodel@ov.x.dummy.lv.idx[[g]],
      remove_dummy_lv = TRUE
    )
    # reconstruct full THETA (in case of dummy ov's)
    theta_g <- lav_lisrel_theta(
      mlist = mlist,
      ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[g]],
      ov_x_dummy_ov_idx = lavmodel@ov.x.dummy.ov.idx[[g]],
      ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[g]],
      ov_x_dummy_lv_idx = lavmodel@ov.x.dummy.lv.idx[[g]]
    )

    # setnames
    set_names <- lav_pt_efa_values(lavpartable)

    # for each set
    for (set in seq_len(lavmodel@nefa)) {
      # check if we have any user=7 elements in this set
      # if not, skip constraints
      ind_idx <- which(lavpartable$op == "=~" &
        lavpartable$group == g &
        lavpartable$efa == set_names[set])
      if (!any(lavpartable$user[ind_idx] == 7L)) {
        next
      }

      # which ov/lv's are involved in this set?
      ov_idx <- lavmodel@ov.efa.idx[[g]][[set]]
      lv_idx <- lavmodel@lv.efa.idx[[g]][[set]]

      # empty set?
      if (length(ov_idx) == 0L) {
        next
      }

      # just 1 factor?
      if (length(lv_idx) < 2L) {
        next
      }

      a <- lambda_g[ov_idx, lv_idx, drop = FALSE]
      p <- nrow(a)
      m <- ncol(a)

      # for oblique, we also need PSI
      if (!ropts$orthogonal) {
        mm_psi <- mlist$psi[lv_idx, lv_idx, drop = FALSE]
      }

      # std.ov? we use diagonal of Sigma for this set of ov's only
      if (ropts$std_ov) {
        mm_theta <- theta_g[ov_idx, ov_idx, drop = FALSE]
        sigma_1 <- tcrossprod(a) + mm_theta
        this_ov_var <- diag(sigma_1)
      } else {
        this_ov_var <- rep(1, p)
      }

      # choose method
      method <- tolower(method)
      if (method %in% c(
        "cf-quartimax", "cf-varimax", "cf-equamax",
        "cf-parsimax", "cf-facparsim"
      )) {
        method_fname <- "lav_mat_rotate_cf"
        this_method_args$cf_gamma <- switch(method,
          "cf-quartimax" = 0,
          "cf-varimax"   = 1 / p,
          "cf-equamax"   = m / (2 * p),
          "cf-parsimax"  = (m - 1) / (p + m - 2),
          "cf-facparsim" = 1
        )
      } else if (method == "target.strict") {
        method_fname <- "lav_mat_rotate_target"
      } else {
        method_fname <- paste("lav_mat_rotate_", method, sep = "")
      }

      # check if rotation method exists
      check <- try(get(method_fname), silent = TRUE)
      if (inherits(check, "try-error")) {
        lav_msg_stop(gettextf("unknown rotation method: %s", method_fname))
      }

      # 1. compute row weigths

      # 1.a cov -> cor?
      if (ropts$std_ov) {
        a <- a * 1 / sqrt(this_ov_var)
      }

      if (ropts$row_weights == "none") {
        weights <- rep(1.0, p)
      } else if (ropts$row_weights == "kaiser") {
        weights <- lav_mat_rotate_kaiser_weights(a)
      } else if (ropts$row_weights == "cureton-mulaik") {
        weights <- lav_mat_rotate_cm_weights(a)
      } else {
        lav_msg_stop(gettextf("row_weights can be",
          lav_msg_view(c("none", "kaiser", "cureton-mulaik"), "or")))
      }
      a <- a * weights

      # evaluate rotation criterion, extract GRAD
      q_1 <- do.call(
        method_fname,
        c(list(mm_lambda = a), this_method_args, list(grad = TRUE))
      )
      gq <- attr(q_1, "grad")
      attr(q_1, "grad") <- NULL

      # compute 'Z'
      z <- crossprod(a, gq)

      # compute constraints
      if (ropts$orthogonal) {
        # the constraint: Z == diagonal
        # or in other words, the non-diagonal elements of
        # Z - t(Z) are all zero
        tmp <- z - t(z)
        this_res <- lav_mat_vech(tmp, diagonal = FALSE)
      } else {
        psi_z <- mm_psi * diag(z) # rescale rows only
        tmp <- z - psi_z
        out1 <- lav_mat_vech(tmp, diagonal = FALSE)
        out2 <- lav_mat_vechu(tmp, diagonal = FALSE)
        this_res <- c(out1, out2)
      }

      res <- c(res, this_res)
    } # set
  } # group

  # return constraint vector
  res
}
