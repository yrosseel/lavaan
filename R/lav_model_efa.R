# efa related functions
# YR - April 2019
#
# the lav_model_efa_rotate_x() function was based on a script originally
# written by Florian Scharf (Muenster University, Germany)

# apply a (block-level) rotation matrix 'hg' to all factor-related model
# matrices of a single block: lambda, psi, beta and alpha. This is the common
# transform used both when rotating the solution (lav_model_efa_rotate_x()) and
# when aligning the rotated factors of a bootstrap sample to the reference
# solution (lav_efa_bootstrap_align_coef()). nu and theta need no rotation.
#
# The general formulas below hold for any invertible 'hg'. For an orthonormal
# signed-permutation 'hg' (as used by the bootstrap alignment) they reduce to
# lambda %*% hg, t(hg) %*% beta %*% hg, etc., since solve(t(hg)) == hg.
# translate a subset of the latent variables by t (a vector over all lv's,
# zero for the others): eta' = eta + t. The model-implied moments are
# unchanged when alpha' = alpha + (I - B) t and nu' = nu - lambda t. Used
# by the alignment of the intercepts (mg_agreement_intercepts)
lav_model_efa_apply_translation <- function(mlist = NULL, t_vec = NULL) {
  if (!is.null(mlist$alpha)) {
    shift <- t_vec
    if (!is.null(mlist$beta)) {
      shift <- shift - as.numeric(mlist$beta %*% t_vec)
    }
    mlist$alpha <- mlist$alpha + shift
  }
  if (!is.null(mlist$nu)) {
    mlist$nu <- mlist$nu - as.numeric(mlist$lambda %*% t_vec)
  }
  mlist
}

lav_model_efa_apply_hg <- function(mlist = NULL, hg = NULL) {
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

  mlist
}

# reconstruct the full block-g LAMBDA from a block's model matrices, filling in
# the dummy latent variables used for observed-variable predictors/outcomes (and
# then removing them). Used by the rotation routines below.
lav_model_efa_block_lambda <- function(lavmodel = NULL, mlist = NULL, g = 1L) {
  lav_lisrel_lambda(
    mlist = mlist,
    ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[g]],
    ov_x_dummy_ov_idx = lavmodel@ov.x.dummy.ov.idx[[g]],
    ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[g]],
    ov_x_dummy_lv_idx = lavmodel@ov.x.dummy.lv.idx[[g]],
    remove_dummy_lv = TRUE
  )
}

# select the target (or target_mask) matrix for a single efa/esem set: a
# *named* list is keyed by the efa block labels (issue #516), so we return the
# element matching 'set_name'; an unnamed list is interpreted per group and is
# returned as-is (the group argument of lav_mat_rotate() selects the element
# downstream); a plain matrix (used for all sets) is also returned as-is.
lav_model_efa_target_set <- function(x = NULL, set_name = NULL,
                                     what = "target") {
  if (!is.list(x)) {
    return(x)
  }
  x_names <- names(x)
  if (is.null(x_names) || !any(nzchar(x_names))) {
    # unnamed list: per-group target/target_mask
    return(x)
  }
  if (length(set_name) == 0L || is.na(set_name)) {
    lav_msg_stop(gettextf(
      "rotation_args$%s is a named list, but the efa block labels could not be found in the model.",
      what))
  }
  idx <- match(set_name, x_names)
  if (is.na(idx)) {
    lav_msg_stop(gettextf(
      "rotation_args$%1$s is a named list, but contains no element for efa block \"%2$s\".",
      what, set_name))
  }
  x[[idx]]
}

# early validation of *named* per-set target/target_mask lists (issue #516):
# all names must be efa block labels, and every efa block with at least two
# factors needs an entry. Called by lav_model() so a mismatch is caught
# before the model is estimated.
lav_model_efa_check_target <- function(ropts = NULL, lavpartable = NULL,
                                       efa_values = NULL) {
  for (what in c("target", "target_mask")) {
    x <- ropts[[what]]
    if (!is.list(x)) {
      next
    }
    x_names <- names(x)
    if (is.null(x_names) || !any(nzchar(x_names))) {
      # unnamed list: per-group target/target_mask
      next
    }
    unknown <- setdiff(x_names, efa_values)
    if (length(unknown) > 0L) {
      lav_msg_stop(gettextf(
        "some names of the rotation_args$%1$s list (%2$s) do not match any efa block label in the model (%3$s).",
        what, lav_msg_view(unknown, log_sep = "none"),
        lav_msg_view(efa_values, log_sep = "none")))
    }
    for (set in efa_values) {
      nfac <- length(unique(lavpartable$lhs[lavpartable$op == "=~" &
        lavpartable$efa == set]))
      if (nfac > 1L && !set %in% x_names) {
        lav_msg_stop(gettextf(
          "the rotation_args$%1$s list contains no element for efa block \"%2$s\".",
          what, set))
      }
    }
  }
  invisible(NULL)
}

# warn for regressions among the latent variables of the same efa block: the
# factors of a block are rotated jointly, and their relations are already
# determined by the rotation (the factor correlations), so such a
# regression coefficient is not identified (it stays at its starting value
# and the rotated solution shows a factor covariance instead). Called once
# from step 9 (lav_lavaan_step09_model()), before the model is estimated.
lav_model_efa_check_regressions <- function(lavpartable = NULL,
                                            efa_values = NULL) {
  reg_idx <- which(lavpartable$op == "~")
  if (length(reg_idx) == 0L) {
    return(invisible(NULL))
  }
  bad <- character(0L)
  for (set in efa_values) {
    set_lv <- unique(lavpartable$lhs[lavpartable$op == "=~" &
      lavpartable$efa == set])
    this_idx <- reg_idx[lavpartable$lhs[reg_idx] %in% set_lv &
      lavpartable$rhs[reg_idx] %in% set_lv]
    if (length(this_idx) > 0L) {
      bad <- c(bad, unique(paste(lavpartable$lhs[this_idx], "~",
        lavpartable$rhs[this_idx])))
    }
  }
  if (length(bad) > 0L) {
    lav_msg_warn(gettextf(
      "the model contains regressions among latent variables of the same
       efa block (%s). The relations among the factors of an efa block are
       determined by the rotation, and these regression coefficients are
       not identified; their estimates are not meaningful.",
      lav_msg_view(bad, log_sep = "none")))
  }
  invisible(NULL)
}

# the names of the latent variables of the full block-g LAMBDA (as built by
# lav_model_efa_block_lambda(): the dummy latent variables are removed), and
# of its rows (the observed variables)
lav_model_efa_block_lv_names <- function(lavmodel = NULL, mlist = NULL,
                                         g = 1L) {
  mm_in_group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0, lavmodel@nmat))[g]
  lambda_mm <- mm_in_group[which(names(mlist) == "lambda")[1L]]
  lv_names <- lavmodel@dimNames[[lambda_mm]][[2L]]
  dummy_idx <- c(lavmodel@ov.y.dummy.lv.idx[[g]],
                 lavmodel@ov.x.dummy.lv.idx[[g]])
  if (length(dummy_idx) > 0L) {
    lv_names <- lv_names[-dummy_idx]
  }
  lv_names
}
lav_model_efa_block_ov_names <- function(lavmodel = NULL, mlist = NULL,
                                         g = 1L) {
  mm_in_group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0, lavmodel@nmat))[g]
  lambda_mm <- mm_in_group[which(names(mlist) == "lambda")[1L]]
  lavmodel@dimNames[[lambda_mm]][[1L]]
}

# collect the rotation information of one efa set (one block) in a list,
# from the result of lav_mat_rotate() (single group) or lav_mat_rotate_mg()
# (multiple groups + agreement; 'group' selects the group). This is what
# lavInspect(fit, "rotation") returns (per block, per efa set).
#
# rot: the rotation matrix T of the set (rows: the unrotated factors,
#   columns: the rotated factors, in their final order), such that
#   LAMBDA = A %*% T (orthogonal) or A %*% solve(t(T)) (oblique)
# phi: the factor correlation matrix of the rotated factors (T'T)
# order: the column permutation (order_lv_by) applied to the rotated
#   solution: order[k] is the column of the rotated solution that was
#   placed in position k
# starts: the random starts (rotation.args$rstarts): the criterion value
#   at the end of each run, and the corresponding rotation matrix
#   (before reflection and reordering); a failed run has value Inf
lav_model_efa_rotation_info <- function(res = NULL, set_name = NULL,
                                        lv_names = NULL, ov_names = NULL,
                                        group = NULL) {
  m <- length(lv_names)
  mg_flag <- !is.null(group)
  if (mg_flag) {
    rot <- res$rotList[[group]]
    order_idx <- res$orderList[[group]]
    ginfo <- res$group_info[[group]]
    rep_mat <- ginfo$REP
  } else {
    rot <- res$ROT
    order_idx <- res$order.idx
    rep_mat <- res$REP
  }
  dimnames(rot) <- list(lv_names, lv_names)
  # the covariance matrix of the rotated factors (a correlation matrix,
  # except in the multiple-group agreement case, where the column
  # rescaling is folded into rot)
  phi <- crossprod(rot)
  dimnames(phi) <- list(lv_names, lv_names)
  starts <- NULL
  if (!is.null(rep_mat)) {
    # one column per start: value, iterations, converged, vec(rot)
    rep_mat <- as.matrix(rep_mat)
    starts <- list(
      value = rep_mat[1L, ],
      iterations = as.integer(rep_mat[2L, ]),
      converged = rep_mat[3L, ] == 1,
      rot = lapply(seq_len(ncol(rep_mat)), function(i) {
        matrix(rep_mat[-(1:3), i], nrow = m, ncol = m)
      })
    )
  }
  out <- list(
    efa = set_name,
    lv_names = lv_names,
    ov_names = ov_names,
    method = res$method,
    orthogonal = res$orthogonal,
    row_weights = res$row_weights,
    algorithm = res$algorithm,
    converged = res$converged,
    iterations = res$iter,
    criterion = res$method.value,
    rot = rot,
    phi = phi,
    order = order_idx,
    starts = starts
  )
  if (mg_flag) {
    # the joint (agreement) rotation started from a per-group rotation:
    # its criterion value, convergence and random starts are per group
    out$mg_agreement <- TRUE
    out$mg_agreement_scales <- res$mg_agreement_scales
    out$mg_agreement_intercepts <- res$mg_agreement_intercepts
    out$mg_agreement_identification <- res$mg_agreement_identification
    if (!is.null(res$meanList)) {
      means <- res$meanList[[group]]
      names(means) <- lv_names
      out$means <- means
    }
    out$group_criterion <- ginfo$method.value
    out$group_converged <- ginfo$converged
    out$group_iterations <- ginfo$iter
  }
  out
}

# reconstruct the full block-g THETA from a block's model matrices (filling in
# the dummy observed variables). Used by the rotation routines below.
lav_model_efa_block_theta <- function(lavmodel = NULL, mlist = NULL, g = 1L) {
  lav_lisrel_theta(
    mlist = mlist,
    ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[g]],
    ov_x_dummy_ov_idx = lavmodel@ov.x.dummy.ov.idx[[g]],
    ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[g]],
    ov_x_dummy_lv_idx = lavmodel@ov.x.dummy.lv.idx[[g]]
  )
}

# rotate solution
# nobs: the (effective) sample sizes per group; only needed for the
# multiple-group agreement rotation with mg_agreement_group_weights = "nobs"
lav_model_efa_rotate <- function(lavmodel = NULL, lavoptions = NULL,
                                 nobs = NULL) {

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
    lavoptions = lavoptions, extra = TRUE, nobs = nobs
  )
  extra <- attr(tmp, "extra")
  attr(tmp, "extra") <- NULL

  out <- list(GLIST = extra$GLIST, H = extra$H, lv.order = extra$lv.order,
              rotation = extra$rotation)

  # restore warn/verbose settings
  lav_warn(current_warn)
  lav_verbose(current_verbose)

  out
}


# lower-level function, needed for numDeriv
lav_model_efa_rotate_x <- function(x, lavmodel = NULL, lavoptions = NULL,
                                   init_rot = NULL, extra = FALSE,
                                   type = "free", nobs = NULL) {
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

  # H per block
  h <- vector("list", lavmodel@nblocks)
  order_1 <- vector("list", lavmodel@nblocks)
  # rotation information per block and per efa set (only if extra = TRUE;
  # stored in the @internal slot of the lavaan object, see lavInspect())
  rotation_info <- vector("list", lavmodel@nblocks)

  # the aligned factor means per group (intercept alignment), zero if none
  mean_1 <- vector("list", length = lavmodel@ngroups)

  # new in 0.7-1 -- three options:
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
    # iterate over blocks, not groups: each block (eg the within and the
    # between level of a twolevel model) is rotated independently. The
    # efa-related slots (ov.efa.idx, lv.efa.idx, nmat, the dummy-ov indices
    # and H) are all indexed per block; for single-level models block == group,
    # so multigroup behaviour is unchanged.
    for (g in seq_len(lavmodel@nblocks)) {
      # select model matrices for this group
      mm_in_group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0, lavmodel@nmat))[g]
      mlist <- glist[mm_in_group]

      # general rotation matrix (all latent variables)
      h[[g]] <- hg <- diag(ncol(mlist$lambda))
      lv_order <- seq_len(ncol(mlist$lambda))

      # reconstruct full mm_lambda + THETA (in case of dummy ov's)
      lambda_g <- lav_model_efa_block_lambda(lavmodel, mlist, g)
      theta_g <- lav_model_efa_block_theta(lavmodel, mlist, g)

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
            # compute mean THETA across groups, for this set of ov's only
            theta_list <- lapply(seq_len(lavmodel@ngroups), function(gg) {
              mm_in_gg <- (seq_len(lavmodel@nmat[gg]) +
                            cumsum(c(0, lavmodel@nmat))[gg])
              lav_model_efa_block_theta(lavmodel, glist[mm_in_gg], gg)
            })
            mm_theta <- (Reduce("+", theta_list) /
                          lavmodel@ngroups)[ov_idx, ov_idx, drop = FALSE]
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
        # per-set target/target_mask? (named list, keyed by efa block label)
        set_name <- names(lavmodel@lv.efa.idx[[g]])[set]
        this_target <- lav_model_efa_target_set(
          ropts$target, set_name, what = "target")
        this_target_mask <- lav_model_efa_target_set(
          ropts$target_mask, set_name, what = "target_mask")
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
            cf_gamma       = ropts$cf_gamma,
            oblimin_gamma  = ropts$oblimin_gamma,
            promax_kappa   = ropts$promax_kappa,
            simplimax_k    = ropts$simplimax_k,
            target         = this_target,
            target_mask    = this_target_mask
          ),
          init_rot = init_rot_1,
          init_rot_check = FALSE,
          rstarts = rstarts,
          row_weights = ropts$row_weights,
          std_ov = ropts$std_ov,
          ov_var = this_ov_var,
          algorithm = ropts$algorithm,
          pairwise_accelerate = ropts$pairwise_accelerate,
          reflect = ropts$reflect,
          order_lv_by = ropts$order_lv_by,
          gpa_tol = ropts$gpa_tol,
          gpa_algorithm = ropts$gpa_algorithm,
          gpa_fwindow = ropts$gpa_fwindow,
          tol = ropts$tol,
          keep_rep = extra,
          max_iter = ropts$max_iter,
          group = g
        )

        # collect rotation information for this set
        if (extra) {
          lv_names_g <- lav_model_efa_block_lv_names(lavmodel, mlist, g)
          ov_names_g <- lav_model_efa_block_ov_names(lavmodel, mlist, g)
          rotation_info[[g]][[set_name]] <- lav_model_efa_rotation_info(
            res = res, set_name = set_name,
            lv_names = lv_names_g[lv_idx], ov_names = ov_names_g[ov_idx]
          )
        }

        # extract rotation matrix (note, in Asp & Muthen, 2009; this is H')
        # note: as of 0.6-6, order.idx has already been applied to ROT,
        #       so no need to reorder rows/columns after rotation
        h_efa <- res$ROT

        # fill in optimal rotation for this set
        hg[lv_idx, lv_idx] <- h_efa

        # keep track of possible re-orderings
        lv_order[lv_idx] <- lv_idx[res$order.idx]
      } # set

      # rotate all the SEM parameters (lambda, psi, beta, alpha)
      mlist <- lav_model_efa_apply_hg(mlist, hg)

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
        # transform psi/beta/alpha; lambda is forced equal across groups
        mlist <- lav_model_efa_apply_hg(mlist, hg)
        mlist$lambda <- mlist1$lambda
        glist[mm_in_group] <- mlist
        h[[g]] <- hg
        order_1[[g]] <- lv_order
        rotation_info[[g]] <- rotation_info[[1L]]
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
      lambda_g <- lav_model_efa_block_lambda(lavmodel, mlist, g)

      lambda_list[[g]] <- lambda_g

      # prepare H, ORDER and the factor means
      h[[g]] <- diag(ncol(lambda_g))
      order_1[[g]] <- seq_len(ncol(lambda_g))
      mean_1[[g]] <- numeric(ncol(lambda_g))
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

      # note: single-factor sets and std.ov are not (yet) supported in the
      # multigroup rotation + agreement path; this_ov_var is always NULL here

      # unrotated 'A' for this set
      alist_1 <- lapply(lambda_list,
        function(x) x[ov_idx, lv_idx, drop = FALSE])
      this_ov_var <- NULL

      # mg_agreement_std_ov: the agreement objective (and the per-group
      # pre-rotation) on the loadings of the standardized indicators,
      # standardized over the entire population (the model-implied
      # variances and means of the groups pooled by sample size); the
      # rotation matrices, scales and means are the same in both metrics
      ov_var_list <- NULL
      sd_pooled <- NULL
      if (isTRUE(ropts$mg_agreement_std_ov)) {
        n_g <- unlist(nobs)
        if (length(n_g) != lavmodel@ngroups || any(!is.finite(n_g))) {
          lav_msg_stop(gettext(
            "mg_agreement_std_ov = TRUE requires the sample sizes of all
             groups."))
        }
        ov_var_list <- vector("list", length = lavmodel@ngroups)
        mu_list <- vector("list", length = lavmodel@ngroups)
        for (g in seq_len(lavmodel@ngroups)) {
          mm_in_group <- (seq_len(lavmodel@nmat[g]) +
                          cumsum(c(0, lavmodel@nmat))[g])
          mlist_g <- glist[mm_in_group]
          lambda_g <- lav_model_efa_block_lambda(lavmodel, mlist_g, g)
          theta_g <- lav_model_efa_block_theta(lavmodel, mlist_g, g)
          sigma_g <- lambda_g %*% mlist_g$psi %*% t(lambda_g) + theta_g
          ov_var_list[[g]] <- diag(sigma_g)[ov_idx]
          mu_list[[g]] <- if (is.null(mlist_g$nu)) {
            numeric(length(ov_idx))
          } else {
            as.numeric(mlist_g$nu[ov_idx, 1L])
          }
        }
        w_n <- n_g / sum(n_g)
        mu_pooled <- Reduce("+", Map("*", mu_list, w_n))
        var_pooled <- Reduce("+", Map(function(v, mu, w) {
          w * (v + (mu - mu_pooled)^2)
        }, ov_var_list, mu_list, w_n))
        sd_pooled <- sqrt(var_pooled)
        alist_1 <- lapply(alist_1, function(x) x / sd_pooled)
      }

      # intercept alignment: the intercepts of the indicators of this set
      nu_list <- NULL
      if (ropts$mg_agreement_intercepts != "none") {
        nu_list <- vector("list", length = lavmodel@ngroups)
        for (g in seq_len(lavmodel@ngroups)) {
          mm_in_group <- (seq_len(lavmodel@nmat[g]) +
                          cumsum(c(0, lavmodel@nmat))[g])
          mlist_g <- glist[mm_in_group]
          if (is.null(mlist_g$nu)) {
            lav_msg_stop(gettext(
              "the alignment of the intercepts (mg_agreement_intercepts)
               requires a mean structure (meanstructure = TRUE)."))
          }
          if (lavmodel@categorical || lavmodel@conditional.x) {
            lav_msg_stop(gettext(
              "the alignment of the intercepts (mg_agreement_intercepts) is
               not available for categorical data or conditional.x = TRUE."))
          }
          dummy_idx <- c(lavmodel@ov.y.dummy.ov.idx[[g]],
                         lavmodel@ov.x.dummy.ov.idx[[g]])
          if (length(intersect(ov_idx, dummy_idx)) > 0L) {
            lav_msg_stop(gettext(
              "the alignment of the intercepts (mg_agreement_intercepts) is
               not available when an indicator of an efa block is also a
               predictor or an outcome in a regression."))
          }
          nu_list[[g]] <- as.numeric(mlist_g$nu[ov_idx, 1L])
          if (!is.null(sd_pooled)) {
            nu_list[[g]] <- nu_list[[g]] / sd_pooled
          }
        }
      }

      # per-set target/target_mask? (named list, keyed by efa block label)
      set_name <- names(lavmodel@lv.efa.idx[[1]])[set]
      this_target <- lav_model_efa_target_set(
        ropts$target, set_name, what = "target")
      this_target_mask <- lav_model_efa_target_set(
        ropts$target_mask, set_name, what = "target_mask")

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

      # group weights for the agreement criterion: proportional to the
      # sample sizes (scaled to average one), to their square roots (so
      # that a pair of groups is weighted by sqrt(n_g n_h), as in the
      # alignment method), or none
      mg_group_weights <- NULL
      if (ropts$mg_agreement_group_weights %in% c("nobs", "sqrt.nobs")) {
        n_g <- unlist(nobs)
        if (length(n_g) != lavmodel@ngroups || any(!is.finite(n_g)) ||
            any(n_g <= 0)) {
          lav_msg_stop(gettextf(
            "mg_agreement_group_weights = %s requires the sample sizes of
             all groups.", dQuote(ropts$mg_agreement_group_weights)))
        }
        mg_group_weights <- n_g / mean(n_g)
        if (ropts$mg_agreement_group_weights == "sqrt.nobs") {
          mg_group_weights <- sqrt(mg_group_weights)
        }
      }

      # rotate this set
      res <- lav_mat_rotate_mg(
        a_list = alist_1,
        orthogonal = ropts$orthogonal,
        method = method,
        method_args = list(
          geomin_epsilon = ropts$geomin_epsilon,
          orthomax_gamma = ropts$orthomax_gamma,
          cf_gamma       = ropts$cf_gamma,
          oblimin_gamma  = ropts$oblimin_gamma,
          promax_kappa   = ropts$promax_kappa,
          simplimax_k    = ropts$simplimax_k,
          target         = this_target,
          target_mask    = this_target_mask
        ),
        init_rot_list = init_rot_list,
        init_rot_check = FALSE,
        rstarts = rstarts,
        #row_weights = ropts$row_weights,
        #std.ov = ropts$std_ov,
        #ov.var = this.ov.var,
        mg_algorithm = ropts$mg_agreement_algorithm,
        mg_agreement_method = ropts$mg_agreement_method,
        mg_agreement_epsilon = ropts$mg_agreement_epsilon,
        mg_agreement_weight = ropts$mg_agreement_weight,
        mg_agreement_group_weights = mg_group_weights,
        mg_agreement_scales = ropts$mg_agreement_scales,
        mg_agreement_rotate = ropts$mg_agreement_rotate,
        mg_agreement_intercepts = ropts$mg_agreement_intercepts,
        mg_agreement_identification = ropts$mg_agreement_identification,
        nu_list = nu_list,
        std_ov = isTRUE(ropts$mg_agreement_std_ov),
        ov_var_list = ov_var_list,
        pairwise_accelerate = ropts$pairwise_accelerate,
        reflect = ropts$reflect,
        order_lv_by = ropts$order_lv_by,
        gpa_tol = ropts$gpa_tol,
        tol = ropts$tol,
        keep_rep = extra,
        max_iter = ropts$max_iter
      )
      # rescale the rotated factors of this set so that the sum of squares
      # per column is the same for all groups (the agreement criterion is
      # evaluated on such rescaled loadings; the scale moves into psi).
      # The rescaling is folded into the rotation matrices, so that H is
      # the complete transformation and lambda, psi, beta and alpha stay
      # consistent (before 0.7-3, lambda and psi were rescaled afterwards,
      # beta and alpha were not, and the non-efa factors were rescaled too)
      # With mg_agreement_scales = "free", the scales were part of the
      # objective and are already in the rotation matrices
      if (ropts$mg_agreement_scales != "free") {
        ss_list <- lapply(res$lambdaList, function(x) colSums(x * x))
        ss_ave <- lav_mat_rotate_mg_ss_ref(ss_list,
                                           ropts$mg_agreement_identification)
        for (g in seq_len(lavmodel@ngroups)) {
          scale_g <- sqrt(ss_list[[g]] / ss_ave)
          res$rotList[[g]] <- t(t(res$rotList[[g]]) * scale_g)
        }
      }
      for (g in seq_len(lavmodel@ngroups)) {
        # fill in optimal rotation for this set, for this group
        h[[g]][lv_idx, lv_idx] <- res$rotList[[g]]

        # keep track of possible re-orderings
        order_1[[g]][lv_idx] <- lv_idx[res$orderList[[g]]]

        # the aligned factor means of this set (translation, applied
        # after the rotation below)
        if (!is.null(res$meanList)) {
          mean_1[[g]][lv_idx] <- res$meanList[[g]]
        }

        # collect rotation information for this set, for this group
        if (extra) {
          mm_in_group <- (seq_len(lavmodel@nmat[g]) +
                          cumsum(c(0, lavmodel@nmat))[g])
          mlist <- glist[mm_in_group]
          lv_names_g <- lav_model_efa_block_lv_names(lavmodel, mlist, g)
          ov_names_g <- lav_model_efa_block_ov_names(lavmodel, mlist, g)
          rotation_info[[g]][[set_name]] <- lav_model_efa_rotation_info(
            res = res, set_name = set_name,
            lv_names = lv_names_g[lv_idx], ov_names = ov_names_g[ov_idx],
            group = g
          )
        }
      }
    } # set

    # rotate all the SEM parameters
    for (g in seq_len(lavmodel@ngroups)) {
      mm_in_group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0, lavmodel@nmat))[g]
      mlist <- glist[mm_in_group]
      hg <- h[[g]]

      # rotate all the SEM parameters (lambda, psi, beta, alpha); the
      # rescaling to equal column sums of squares across groups is part
      # of hg (see above)
      mlist <- lav_model_efa_apply_hg(mlist, hg)

      # intercept alignment: translate the (rotated) efa factors by their
      # aligned means t: eta' = eta + t keeps the model-implied moments
      # when alpha' = alpha + (I - B) t and nu' = nu - lambda t
      t_g <- mean_1[[g]]
      if (any(t_g != 0)) {
        mlist <- lav_model_efa_apply_translation(mlist, t_g)
      }

      glist[mm_in_group] <- mlist
    }
  } # option 3

  # extract all rotated parameter estimates
  x_rot <- lav_model_get_parameters(lavmodel, glist = glist, type = type)

  # extra?
  if (extra) {
    attr(x_rot, "extra") <- list(GLIST = glist, H = h, lv.order = order_1,
                                 rotation = rotation_info)
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
    cf_gamma = ropts$cf_gamma,
    oblimin_gamma = ropts$oblimin_gamma,
    promax_kappa = ropts$promax_kappa,
    simplimax_k  = ropts$simplimax_k,
    target = ropts$target,
    target_mask = ropts$target_mask
  )

  # place parameters into model matrices
  lavmodel <- lav_model_set_parameters(lavmodel, x = x)

  # GLIST
  glist <- lavmodel@GLIST

  # res
  res <- numeric(0L)

  # per block (group x level): the rotation constraints must be generated for
  # every rotated block, otherwise the bordered information matrix is rank
  # deficient (and the rotated vcov is not positive definite). The efa-related
  # slots are indexed per block; for single-level models block == group.
  for (g in seq_len(lavmodel@nblocks)) {

    # group-specific method.args
  this_method_args <- method_args

  # set group-specific target/target_mask (if needed)
    # if target, check target matrix; an *unnamed* list is interpreted per
    # group; a *named* list is keyed by the efa block labels and is resolved
    # per set (in the set loop below)
    if (method == "target.strict" || method == "pst") {
      target <- method_args$target
      if (is.list(target) &&
          (is.null(names(target)) || !any(nzchar(names(target))))) {
        this_method_args$target <- target[[g]]
      }
    }
    if (method == "pst") {
      target_mask <- method_args$target_mask
      if (is.list(target_mask) &&
          (is.null(names(target_mask)) || !any(nzchar(names(target_mask))))) {
        this_method_args$target_mask <- target_mask[[g]]
      }
    }

    # select model matrices for this group
    mm_in_group <- seq_len(lavmodel@nmat[g]) + cumsum(c(0, lavmodel@nmat))[g]
    mlist <- glist[mm_in_group]

    # reconstruct full mm_lambda + THETA (in case of dummy ov's)
    lambda_g <- lav_model_efa_block_lambda(lavmodel, mlist, g)
    theta_g <- lav_model_efa_block_theta(lavmodel, mlist, g)

    # setnames
    set_names <- lav_pt_efa_values(lavpartable)

    # for each set
    for (set in seq_len(lavmodel@nefa)) {
      # check if we have any user=7 elements in this set
      # if not, skip constraints
      ind_idx <- which(lavpartable$op == "=~" &
        lavpartable$block == g &
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

      # per-set target/target_mask? (named list, keyed by efa block label)
      set_method_args <- this_method_args
      set_method_args$target <- lav_model_efa_target_set(
        this_method_args$target, set_names[set], what = "target")
      set_method_args$target_mask <- lav_model_efa_target_set(
        this_method_args$target_mask, set_names[set], what = "target_mask")

      # for oblique, we also need PSI
      if (!ropts$orthogonal) {
        mm_psi <- mlist$psi[lv_idx, lv_idx, drop = FALSE]
      }

      # std.ov? we use diagonal of Sigma for this set of ov's only
      # note: unlike in lav_model_efa_rotate(), 'a' is the *rotated* lambda
      # here, so for oblique rotation the implied variances involve PSI
      if (ropts$std_ov) {
        mm_theta <- theta_g[ov_idx, ov_idx, drop = FALSE]
        if (ropts$orthogonal) {
          sigma_1 <- tcrossprod(a) + mm_theta
        } else {
          sigma_1 <- a %*% mm_psi %*% t(a) + mm_theta
        }
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
        set_method_args$cf_gamma <- switch(method,
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
        c(list(mm_lambda = a), set_method_args, list(grad = TRUE))
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
  } # block

  # return constraint vector
  res
}
