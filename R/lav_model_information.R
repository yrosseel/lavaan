# here, we compute various versions of the `information' matrix
# NOTE:
# 1) we ALWAYS compute the UNIT information (not the total information)
#
# 2) by default, we ignore the constraints (we deal with this when we
#    take the inverse later on)

lav_model_info <- function(lavmodel = NULL,
                                  lavsamplestats = NULL,
                                  lavdata = NULL,
                                  lavimplied = NULL,
                                  lavh1 = NULL,
                                  delta = NULL,
                                  lavcache = NULL,
                                  lavoptions = NULL,
                                  extra = FALSE,
                                  augmented = FALSE,
                                  inverted = FALSE,
                                  use_ginv = FALSE) {

  information <- lavoptions$information[1] # ALWAYS used the first one
  # caller can control it

  # rotation?
  # if(!is.null(lavoptions$rotation) && lavoptions$rotation != "none") {
  #    use.ginv <- TRUE
  # }

  if (is.null(lavh1)) {
    lavh1 <- lav_h1_implied_logl(
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      lavoptions = lavoptions
    )
  }

  # compute information matrix
  if (information == "observed") {
    if (lavsamplestats@missing.flag || lavdata@nlevels > 1L) {
      group_weight <- FALSE
    } else {
      group_weight <- TRUE
    }
    m_e <- lav_model_info_observed(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, group_weight = group_weight,
      lavoptions = lavoptions, extra = extra,
      augmented = augmented, inverted = inverted, use_ginv = use_ginv
    )
  } else if (information == "expected") {
    m_e <- lav_model_info_expected(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions, extra = extra,
      augmented = augmented, inverted = inverted, use_ginv = use_ginv
    )
  } else if (information == "first.order") {
    m_e <- lav_model_info_firstorder(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions, # extra = extra,
      check_pd = FALSE,
      augmented = augmented, inverted = inverted, use_ginv = use_ginv
    )
  }

  # information, augmented information, or inverted information
  m_e
}

# sum a list of per-group information blocks
lav_model_info_group_sum <- function(info_group) {
  information <- info_group[[1]]
  if (length(info_group) > 1L) {
    for (g in 2:length(info_group)) {
      information <- information + info_group[[g]]
    }
  }
  information
}

# assemble the parameter-space information from the h1 information:
#
#     Info = sum_g fg * t(Delta_g) %*% A1_g %*% Delta_g
#
# with a diagonal shortcut when A1_g is stored as the diagonal vector of
# the weight matrix (the DWLS/ULS @WLS.VD slot); note that
# lav_model_info_expected_mlm() passes diagonal = FALSE, as its A1 is
# always a full matrix, whatever the estimator
lav_model_info_assemble <- function(lavmodel = NULL,
                                    lavsamplestats = NULL,
                                    delta = NULL,
                                    a1 = NULL,
                                    diagonal =
                                      lavmodel@estimator %in%
                                        c("DWLS", "ULS")) {
  info_group <- vector("list", length = lavsamplestats@ngroups)
  for (g in 1:lavsamplestats@ngroups) {
    # note LISREL documentation suggests (Ng - 1) instead of Ng...
    fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
    # compute information for this group
    if (diagonal) {
      # diagonal weight matrix
      delta2 <- sqrt(a1[[g]]) * delta[[g]]
      info_group[[g]] <- fg * crossprod(delta2)
    } else {
      # full weight matrix
      info_group[[g]] <- fg * (crossprod(delta[[g]], a1[[g]]) %*% delta[[g]])
    }
  }

  # assemble over groups
  lav_model_info_group_sum(info_group)
}

# fisher/expected information
#
# information = Delta' I1 Delta, where I1 is the unit information of
# the saturated model (evaluated either at the structured or unstructured
# estimates)
lav_model_info_expected <- function(lavmodel = NULL,
                                           lavsamplestats = NULL,
                                           lavdata = NULL,
                                           lavoptions = NULL,
                                           lavimplied = NULL,
                                           lavh1 = NULL,
                                           delta = NULL,
                                           lavcache = NULL,
                                           extra = FALSE,
                                           augmented = FALSE,
                                           inverted = FALSE,
                                           use_ginv = FALSE) {
  if (inverted) {
    augmented <- TRUE
  }

  # 1. Delta
  if (is.null(delta)) {
    delta <- lav_model_delta(lavmodel = lavmodel)
  }


  # two-level least-squares? (WLS.V/WLS.VD hold the stacked two-level
  # weight matrix; the standard Delta' A1 Delta path applies per group)
  wls_2l <- (lavdata@nlevels > 1L &&
    lavmodel@estimator %in% c("WLS", "DWLS", "ULS"))

  # GLS/ML: A1 is bdiag(W, 0.5 t(D) (W %x% W) D) with W the inverse of
  # the sample (GLS, ML-unstructured) or model-implied (ML-structured)
  # covariance matrix; when the caller does not need the full matrix
  # (extra = FALSE), stream A1 %*% Delta column-wise instead of forming
  # the (large) A1 -- O(nvar^3) per free parameter instead of O(nvar^4)
  gls_stream <- (lavmodel@estimator == "GLS" && !extra &&
    lavdata@nlevels == 1L &&
    is.null(lavsamplestats@WLS.V[[1]]) &&
    !lavmodel@correlation &&
    !lavmodel@conditional.x &&
    !lavmodel@group.w.free)
  ml_stream <- (lavmodel@estimator == "ML" && !extra &&
    lavdata@nlevels == 1L &&
    !lavsamplestats@missing.flag &&
    !lavmodel@correlation &&
    !lavmodel@conditional.x &&
    !lavmodel@group.w.free)

  # 2. H1 information (single level, or two-level least-squares)
  # for two-level ML there is no pstar-space A1 here (the information is
  # computed directly in parameter space below); a1 stays NULL and the
  # extra = TRUE attribute "WLS.V" is simply absent
  a1 <- NULL
  if (gls_stream || ml_stream) {
    # no A1 needed
  } else if (lavdata@nlevels == 1L || wls_2l) {
    a1 <- lav_model_h1_info_expected(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      lavimplied = lavimplied,
      lavh1 = lavh1,
      lavcache = lavcache
    )
  } else {
    # force conditional.x = FALSE
    lavimplied <- lav_model_implied_cond2uncond(lavimplied)
  }

  # 3. compute Information per group + assemble over groups
  if (gls_stream || ml_stream) {
    # ML: structured (default) evaluates W at the model-implied moments,
    # unstructured at the sample moments (mirroring the corresponding
    # branch in lav_model_h1_information); GLS is always sample-based
    ml_structured <- TRUE
    if (!is.null(lavoptions) &&
      !is.null(lavoptions$h1.information[1]) &&
      lavoptions$h1.information[1] == "unstructured") {
      ml_structured <- FALSE
    }
    if (ml_stream && ml_structured && length(lavimplied) == 0L) {
      lavimplied <- lav_model_implied(lavmodel = lavmodel)
    }
    info_group <- vector("list", length = lavsamplestats@ngroups)
    for (g in 1:lavsamplestats@ngroups) {
      fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
      v11_scale <- 1.0
      if (gls_stream) {
        w_inv <- lavsamplestats@icov[[g]]
        # mimic Mplus: V11 rescale (full-data input only, see
        # lav_model_h1_information)
        if (isTRUE(lavmodel@estimator.args$gls.v11.mplus) &&
          lavmodel@meanstructure && lavdata@data.type == "full") {
          v11_scale <- lavsamplestats@nobs[[g]] /
            (lavsamplestats@nobs[[g]] - 1)
        }
      } else if (ml_structured) {
        # same inversion route as lav_mvn_info_expected
        w_inv <- lav_mvn_sigma_inv(
          sigma_1 = lavimplied$cov[[g]],
          sinv_method = "eigen"
        )
      } else {
        w_inv <- lavsamplestats@icov[[g]]
      }
      wd <- lav_samp_wls_v_nt_prod(
        m_icov = w_inv,
        x = delta[[g]],
        meanstructure = lavmodel@meanstructure,
        fixed_x = lavmodel@fixed.x,
        x_idx = lavsamplestats@x.idx[[g]],
        v11_scale = v11_scale
      )
      info_group[[g]] <- fg * crossprod(delta[[g]], wd)
    }
    information <- lav_model_info_group_sum(info_group)
  } else if (lavdata@nlevels > 1L && !wls_2l) {
    # multilevel (ML): the information is computed directly in parameter
    # space (no pstar-space A1)
    info_group <- vector("list", length = lavsamplestats@ngroups)
    for (g in 1:lavsamplestats@ngroups) {
      fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal

      # here, we assume only 2 levels, at [[1]] and [[2]]
      if (lavoptions$h1.information[1] == "structured") {
        implied <- lavimplied
      } else {
        implied <- lavh1$implied
      }
      sigma_w <- implied$cov[[(g - 1) * 2 + 1]]
      mu_w <- implied$mean[[(g - 1) * 2 + 1]]
      sigma_b <- implied$cov[[(g - 1) * 2 + 2]]
      mu_b <- implied$mean[[(g - 1) * 2 + 2]]

      info_g <-
        lav_mvn_cl_info_expected_delta(
          lp = lavdata@Lp[[g]],
          delta = delta[[g]],
          mu_w = mu_w,
          sigma_w = sigma_w,
          mu_b = mu_b,
          sigma_b = sigma_b,
          sinv_method = "eigen"
        )
      info_group[[g]] <- fg * info_g
    } # g
    information <- lav_model_info_group_sum(info_group)
  } else {
    information <- lav_model_info_assemble(
      lavmodel = lavmodel, lavsamplestats = lavsamplestats,
      delta = delta, a1 = a1
    )
  }

  # 5. augmented information?
  if (augmented) {
    information <-
      lav_model_info_augment_invert(
        lavmodel = lavmodel,
        information = information,
        inverted = inverted,
        use_ginv = use_ginv
      )
  }

  if (extra) {
    attr(information, "Delta") <- delta
    attr(information, "WLS.V") <- a1 # unweighted
  }

  # possibly augmented/inverted
  information
}

# only for Mplus MLM
lav_model_info_expected_mlm <- function(lavmodel = NULL,
                                               lavsamplestats = NULL,
                                               m_delta = NULL,
                                               extra = FALSE,
                                               augmented = FALSE,
                                               inverted = FALSE,
                                               use_ginv = FALSE) {
  if (inverted) {
    augmented <- TRUE
  }

  if (is.null(m_delta)) {
    m_delta <- lav_model_delta(lavmodel = lavmodel)
  }

  # compute A1
  a1 <- vector("list", length = lavsamplestats@ngroups)
  if (lavmodel@group.w.free) {
    gw <- unlist(lav_model_gw(lavmodel = lavmodel))
  }
  for (g in 1:lavsamplestats@ngroups) {
    a1[[g]] <- lav_mvn_h1_info_expected(
      sample_cov     = lavsamplestats@cov[[g]],
      sample_cov_inv = lavsamplestats@icov[[g]],
      x_idx          = lavsamplestats@x.idx[[g]]
    )
    # the same as GLS... (except for the N/N-1 scaling)
    if (lavmodel@group.w.free) {
      # unweight!!
      a <- exp(gw[g]) / lavsamplestats@nobs[[g]]
      # a <- exp(gw[g]) * lavsamplestats@ntotal / lavsamplestats@nobs[[g]]
      a1[[g]] <- lav_mat_bdiag(matrix(a, 1, 1), a1[[g]])
    }
  }

  # compute Information per group + assemble over groups
  information <- lav_model_info_assemble(
    lavmodel = lavmodel, lavsamplestats = lavsamplestats,
    delta = m_delta, a1 = a1, diagonal = FALSE
  )

  # augmented information?
  if (augmented) {
    information <-
      lav_model_info_augment_invert(
        lavmodel = lavmodel,
        information = information,
        inverted = inverted,
        use_ginv = use_ginv
      )
  }

  if (extra) {
    attr(information, "Delta") <- m_delta
    attr(information, "WLS.V") <- a1 # unweighted
  }

  information
}


lav_model_info_observed <- function(lavmodel = NULL,
                                           lavsamplestats = NULL,
                                           lavdata = NULL,
                                           lavimplied = NULL,
                                           lavh1 = NULL,
                                           lavcache = NULL,
                                           lavoptions = NULL,
                                           extra = FALSE,
                                           group_weight = TRUE,
                                           augmented = FALSE,
                                           inverted = FALSE,
                                           use_ginv = FALSE) {
  if (inverted) {
    augmented <- TRUE
  }

  # observed.information:
  #     - "hessian": second derivative of objective function
  #     - "h1": observed information matrix of saturated (h1) model,
  #             pre- and post-multiplied by the jacobian of the model
  #             parameters (Delta), usually evaluated at the structured
  #             sample statistics (but this depends on the h1.information
  #             option)
  if (!is.null(lavoptions) &&
    !is.null(lavoptions$observed.information[1]) &&
    lavoptions$observed.information[1] == "h1") {
    observed_information <- "h1"
  } else {
    observed_information <- "hessian"
  }

  # HESSIAN based
  if (observed_information == "hessian") {
    hessian <- lav_model_hessian(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      lavcache = lavcache,
      group_weight = group_weight,
      ceq_simple = FALSE
    )

    # NOTE! What is the relationship between the Hessian of the objective
    # function, and the `information' matrix (unit or total)

    # 1. in lavaan, we ALWAYS minimize, so the Hessian is already pos def
    # 2. currently, all estimators give unit information, except MML and PML
    #    so, no need to divide by N
    information <- hessian

    # divide by 'N' for MML and PML
    if (lavmodel@estimator == "PML" || lavmodel@estimator == "MML") {
      information <- information / lavsamplestats@ntotal
      # HJ: Does this need to be divided by sum of weights instead?
    }

    # if multilevel, we should divide by 'J', the number of clusters
    if (lavdata@nlevels > 1L) {
      nc <- 0
      for (g in 1:lavsamplestats@ngroups) {
        nc <- nc + lavdata@Lp[[g]]$nclusters[[2]]
      }
      information <- information * lavsamplestats@ntotal / nc
    }
  }

  # using 'observed h1 information'
  # we need DELTA and 'WLS.V' (=A1)

  if (observed_information == "h1" || extra) {
    # 1. Delta
    delta <- lav_model_delta(lavmodel = lavmodel)

    # 2. H1 information

    a1 <- lav_model_h1_info_observed(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      lavimplied = lavimplied,
      lavh1 = lavh1,
      lavcache = lavcache
    )
  }

  if (observed_information == "h1") {
    # compute Information per group + assemble over groups
    information <- lav_model_info_assemble(
      lavmodel = lavmodel, lavsamplestats = lavsamplestats,
      delta = delta, a1 = a1
    )
  }

  # augmented information?
  if (augmented) {
    information <-
      lav_model_info_augment_invert(
        lavmodel = lavmodel,
        information = information,
        inverted = inverted,
        use_ginv = use_ginv
      )
  }

  if (extra) {
    attr(information, "Delta") <- delta
    attr(information, "WLS.V") <- a1
  }

  information
}

# outer product of the case-wise scores (gradients)
# HJ 18/10/23: Need to divide sum of crossproduct of individual log-likelihoods
# by sum of weights rather than sample size.
lav_model_info_firstorder <- function(lavmodel = NULL,
                                             lavsamplestats = NULL,
                                             lavdata = NULL,
                                             lavimplied = NULL,
                                             lavh1 = NULL,
                                             lavcache = NULL,
                                             lavoptions = NULL,
                                             check_pd = FALSE,
                                             extra = FALSE,
                                             augmented = FALSE,
                                             inverted = FALSE,
                                             use_ginv = FALSE) {
  if (!lavmodel@estimator %in% c("ML", "PML")) {
    lav_msg_stop(gettext(
      "information = \"first.order\" not available for estimator"),
      sQuote(lavmodel@estimator))
  }

  if (inverted) {
    augmented <- TRUE
  }

  # random slopes (rv() modifier): the h1 route is not available;
  # compute the first-order information directly from (numerically
  # obtained) per-cluster scores in model-parameter space; the 'unit'
  # convention for two-level models is per *cluster*: B0 =
  # crossprod(scores) / nclusters (matching the hessian-based observed
  # information, which is also expressed per cluster)
  if (length(lavmodel@rv.ov) > 0L || length(lavmodel@rv.lv) > 0L) {
    g <- 1L # single group only (enforced elsewhere)
    rs <- lavcache[[g]]$rs
    if (is.null(rs)) {
      rs_info <- lav_mvn_cl_rs_info(lavmodel = lavmodel,
                                    lavdata = lavdata)
      rs_stats <- lav_mvn_cl_rs_stats(
        y1 = lavdata@X[[g]], lp = lavdata@Lp[[g]], rs_info = rs_info
      )
      rs <- list(info = rs_info, stats = rs_stats)
    }
    sc <- lav_mvn_cl_rs_scores(lavmodel = lavmodel, rs = rs)
    nclusters <- lavdata@Lp[[g]]$nclusters[[2]]
    information <- crossprod(sc) / nclusters
    b0_group <- list(information)

    if (augmented) {
      information <-
        lav_model_info_augment_invert(
          lavmodel = lavmodel,
          information = information,
          check_pd = check_pd,
          inverted = inverted,
          use_ginv = use_ginv
        )
    }
    if (extra) {
      attr(information, "B0.group") <- b0_group
    }
    return(information)
  }

  b0_group <- vector("list", lavsamplestats@ngroups)

  # 1. Delta
  delta <- lav_model_delta(lavmodel = lavmodel)

  # 2. H1 information
  b1 <- lav_model_h1_info_firstorder(
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    lavimplied = lavimplied,
    lavh1 = lavh1,
    lavcache = lavcache
  )

  # 3. compute Information per group
  info_group <- vector("list", length = lavsamplestats@ngroups)
  for (g in 1:lavsamplestats@ngroups) {
    # unweighted (needed in lav_test?)
    b0_group[[g]] <- t(delta[[g]]) %*% b1[[g]] %*% delta[[g]]

    # >>>>>>>> HJ/MK PML CODE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    wt <- lavdata@weights[[g]]
    if (is.null(wt)) {
      fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
    } else {
      totalwt <- sum(unlist(lavdata@weights))
      fg <- sum(wt) / totalwt
    }

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # compute information for this group
    if (lavmodel@estimator == "PML") {
      # (fixed July 2026, resolving the "WHAT IS THE ROLE OF fg?" note):
      # for PML, b0_group is a TOTAL (summed-over-cases) information --
      # the crossprod of the casewise pairwise scores -- so the correct
      # aggregate is the plain sum, divided by N further below. The extra
      # fg weight double-counted the group sizes and, combined with the
      # matching error in the hessian-based information, inflated all
      # multigroup PML standard errors by sqrt(N/n_g). Single-group fits
      # are unaffected (fg = 1).
      info_group[[g]] <- b0_group[[g]]
    } else {
      info_group[[g]] <- fg * b0_group[[g]]
    }
  }

  # 4. assemble over groups
  information <- lav_model_info_group_sum(info_group)

  # NOTE: for MML and PML, we get 'total' information (instead of unit) divide
  # by 'N' for MML and PML. For weighted sample, use the sum of weights
  # instead of sample size

  # >>>>>>>> HJ/MK PML CODE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  if (lavmodel@estimator == "PML" || lavmodel@estimator == "MML") {
    if (length(lavdata@sampling.weights) == 0) {
      the_n <- lavsamplestats@ntotal
    } else {
      the_n <- sum(unlist(lavdata@weights))
    }
    information <- information / the_n
    for (g in 1:lavsamplestats@ngroups) {
      b0_group[[g]] <- b0_group[[g]] / the_n
    }
  }
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # augmented information?
  if (augmented) {
    information <-
      lav_model_info_augment_invert(
        lavmodel = lavmodel,
        information = information,
        check_pd = check_pd,
        inverted = inverted,
        use_ginv = use_ginv
      )
  }

  if (extra) {
    attr(information, "B0.group") <- b0_group
    attr(information, "Delta") <- delta
    attr(information, "WLS.V") <- b1
  }

  information
}


# create augmented information matrix (if needed), and take the inverse
# (if inverted = TRUE), returning only the [1:npar, 1:npar] elements
#
# rm.idx is used by lav_sam_step2_se; it used when the structural model
# contains more parameters than the joint model; therefore, information
# will be 'too small', and we need to remove some columns in H
#
# the constrained inverse is computed via the null-space identity: the
# [1:npar, 1:npar] block of the inverse of the bordered (KKT) matrix
# equals Z (Z' I Z)^-1 Z', with Z an orthonormal basis of the null space
# of the active constraint Jacobian. This needs one solve of order
# npar - rank(h) instead of the SVD of the (npar + 2m) x (npar + 2m)
# bordered matrix, and avoids the MASS::ginv() tolerance altogether. When
# Z'IZ is singular (the model is not identified even on the constrained
# tangent space), we fall back to the bordered matrix + Moore-Penrose
# route, preserving the legacy behavior in that regime. The full bordered
# matrix is still built when inverted = FALSE (lavInspect
# "augmented.information").
lav_model_info_augment_invert <- function(lavmodel = NULL,
                                                 information = NULL,
                                                 inverted = FALSE,
                                                 check_pd = FALSE,
                                                 use_ginv = FALSE,
                                                 rm_idx = integer(0L)) {
  npar <- nrow(information)

  # collect the ACTIVE constraint rows (row-normalized), if any
  h <- NULL
  lambda <- numeric(0L)
  border_type <- "none"
  if (nrow(lavmodel@con.jac) > 0L) {
    border_type <- "general"
    h <- lavmodel@con.jac
    # read the attribute BEFORE any subsetting: `[` drops attributes, so
    # reading it after h[, -rm_idx] silently kept the inactive inequality
    # rows (eg synthesized bounds) and imposed them as active equalities
    # whenever rm_idx was non-empty (the SAM step-2 path)
    inactive_idx <- attr(h, "inactive.idx")
    lambda <- lavmodel@con.lambda # lagrangean coefs
    if (length(inactive_idx) > 0L) {
      h <- h[-inactive_idx, , drop = FALSE]
      lambda <- lambda[-inactive_idx]
    }
    if (length(rm_idx) > 0L) {
      h <- h[, -rm_idx, drop = FALSE]
    }
    # normalize the rows of h: the null space of h (and therefore the
    # constrained inverse) is invariant to row scaling, but crossprod(h)
    # puts the SQUARED row norms into the spectrum of the bordered
    # matrix, and MASS::ginv() uses a tolerance RELATIVE to the largest
    # singular value: a badly scaled constraint (eg 1000*a == 1000*b)
    # raised the cutoff enough to zero genuine directions of the
    # information, silently collapsing the standard errors. Rows that
    # lost all their entries to rm_idx (constraints involving only
    # removed parameters) are dropped.
    row_norm <- sqrt(rowSums(h * h))
    nonzero_idx <- which(row_norm > 0)
    if (length(nonzero_idx) < nrow(h)) {
      h <- h[nonzero_idx, , drop = FALSE]
      lambda <- lambda[nonzero_idx]
      row_norm <- row_norm[nonzero_idx]
    }
    if (nrow(h) > 0L) {
      h <- h / row_norm
      # the multiplier of a rescaled constraint rescales inversely
      lambda <- lambda * row_norm
    }
  } else if (lavmodel@ceq.simple.only) {
    border_type <- "simple"
    h <- t(lav_mat_ortho_complement(lavmodel@ceq.simple.K))
    if (length(rm_idx) > 0L) {
      h <- h[, -rm_idx, drop = FALSE]
      # rows may have lost (some of) their entries; drop empty rows and
      # re-normalize the rest (see the note above)
      row_norm <- sqrt(rowSums(h * h))
      nonzero_idx <- which(row_norm > 0)
      h <- h[nonzero_idx, , drop = FALSE] / row_norm[nonzero_idx]
    }
    lambda <- numeric(nrow(h))
  }
  is_augmented <- !is.null(h) && nrow(h) > 0L

  # Jacobi (diagonal) preconditioning (new in 0.7-2): badly scaled
  # variables induce a purely DIAGONAL ill-conditioning of the
  # information matrix (the condition number grows with the fourth
  # power of the sd ratios), which can break solve() although the
  # underlying problem is perfectly well posed. Conjugate the whole
  # computation by P = diag(1/sqrt(diag(I))): I_p = P I P has unit
  # diagonal, the constraint jacobian transforms as h_p = h P (the
  # restricted inverse transforms covariantly), and the result maps
  # back as V = P V_p P. Only activated for extreme diagonal spreads,
  # so well-scaled problems take the exact same path as before.
  d_pre <- sqrt(abs(diag(information)))
  d_pre[!is.finite(d_pre) | d_pre < .Machine$double.eps] <- 1
  precondition <- (max(d_pre) / min(d_pre)) > 1e6
  if (precondition) {
    info_p <- information / tcrossprod(d_pre)
    h_p <- NULL
    if (is_augmented) {
      h_p <- t(t(h) / d_pre)
    }
  }

  # build the full bordered (KKT) matrix; only needed when the caller
  # asks for the un-inverted augmented matrix, or as the Moore-Penrose
  # fallback of the inverted path below
  build_border <- function() {
    h0 <- matrix(0, nrow(h), nrow(h))
    # adding crossprod(h) leaves the [1:npar, 1:npar] block of the
    # inverse unchanged (h x = 0 implies (info + h'h) x = info x), but
    # makes the bordered matrix nonsingular when the information is
    # singular only along the constraint-normal directions
    info <- information + crossprod(h)
    if (border_type == "general") {
      # note: the middle 2*diag(lambda) block is fully decoupled (zero
      # off-diagonal blocks on both sides), so it cannot affect the
      # [1:npar, 1:npar] slice of the inverse; it is kept for the shape
      # of the full augmented matrix (lavInspect "augmented.information")
      h10 <- matrix(0, npar, nrow(h))
      dl <- 2 * diag(lambda, nrow(h), nrow(h))
      rbind(
        cbind(info, h10, t(h)),
        cbind(t(h10), dl, h0),
        cbind(h, h0, h0)
      )
    } else {
      rbind(
        cbind(info, t(h)),
        cbind(h, h0)
      )
    }
  }

  if (check_pd) {
    # for a constrained model, the relevant identification check is the
    # information restricted to the constrained tangent space (a bordered
    # matrix itself is a saddle matrix and always indefinite)
    # (under preconditioning, check the conjugated matrix: the
    # eigenvalue SIGNS are invariant under the congruence)
    if (is_augmented) {
      if (precondition) {
        z <- lav_mat_ortho_complement(t(h_p))
        m_check <- crossprod(z, info_p %*% z)
      } else {
        z <- lav_mat_ortho_complement(t(h))
        m_check <- crossprod(z, information %*% z)
      }
    } else if (precondition) {
      m_check <- info_p
    } else {
      m_check <- information
    }
    eigvals <- eigen(m_check,
      symmetric = TRUE,
      only.values = TRUE
    )$values
    if (any(eigvals < -1 * .Machine$double.eps^(3 / 4))) {
      lav_msg_warn(gettext(
        "information matrix is not positive definite;
        the model may not be identified"))
    }
  }

  if (!inverted) {
    if (is_augmented) {
      return(build_border())
    }
    return(information)
  }

  # inverted
  if (!is_augmented) {
    info_1 <- if (precondition) info_p else information
    if (use_ginv) {
      # note: default tol in MASS::ginv is sqrt(.Machine$double.eps)
      #       which seems a bit too conservative
      #       from 0.5-20, we changed this to .Machine$double.eps^(3/4)
      out <- try(
        MASS::ginv(info_1,
          tol = .Machine$double.eps^(3 / 4)
        ),
        silent = TRUE
      )
    } else {
      out <- try(solve(info_1), silent = TRUE)
    }
    if (precondition && !inherits(out, "try-error")) {
      out <- out / tcrossprod(d_pre)
    }
    return(out)
  }

  # constrained inverse: null-space route
  out <- try(
    {
      if (precondition) {
        z <- lav_mat_ortho_complement(t(h_p))
        ziz <- crossprod(z, info_p %*% z)
        v <- z %*% solve(ziz, t(z))
        v <- v / tcrossprod(d_pre)
      } else {
        z <- lav_mat_ortho_complement(t(h))
        ziz <- crossprod(z, information %*% z)
        v <- z %*% solve(ziz, t(z))
      }
      # enforce exact symmetry
      (v + t(v)) / 2
    },
    silent = TRUE
  )
  if (!inherits(out, "try-error")) {
    return(out)
  }

  # fallback: bordered matrix + Moore-Penrose
  try(
    MASS::ginv(build_border(),
      tol = .Machine$double.eps^(3 / 4)
    )[1:npar,
      1:npar,
      drop = FALSE
    ],
    silent = TRUE
  )
}
