# here, we compute various versions of the `information' matrix
# NOTE:
# 1) we ALWAYS compute the UNIT information (not the total information)
#
# 2) by default, we ignore the constraints (we deal with this when we
#    take the inverse later on)

lav_model_information <- function(lavmodel = NULL,
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
    m_e <- lav_model_information_observed(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, group_weight = group_weight,
      lavoptions = lavoptions, extra = extra,
      augmented = augmented, inverted = inverted, use_ginv = use_ginv
    )
  } else if (information == "expected") {
    m_e <- lav_model_information_expected(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions, extra = extra,
      augmented = augmented, inverted = inverted, use_ginv = use_ginv
    )
  } else if (information == "first.order") {
    m_e <- lav_model_information_firstorder(
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

# fisher/expected information
#
# information = Delta' I1 Delta, where I1 is the unit information of
# the saturated model (evaluated either at the structured or unstructured
# estimates)
lav_model_information_expected <- function(lavmodel = NULL,
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


  # 2. H1 information (single level)
  if (lavdata@nlevels == 1L) {
    a1 <- lav_model_h1_information_expected(
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

  # 3. compute Information per group
  info_group <- vector("list", length = lavsamplestats@ngroups)
  for (g in 1:lavsamplestats@ngroups) {
    # note LISREL documentation suggests (Ng - 1) instead of Ng...
    fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal

    # multilevel
    if (lavdata@nlevels > 1L) {
      # here, we assume only 2 levels, at [[1]] and [[2]]
      if (lavoptions$h1.information[1] == "structured") {
        sigma_w <- lavimplied$cov[[(g - 1) * 2 + 1]]
        mu_w <- lavimplied$mean[[(g - 1) * 2 + 1]]
        sigma_b <- lavimplied$cov[[(g - 1) * 2 + 2]]
        mu_b <- lavimplied$mean[[(g - 1) * 2 + 2]]
      } else {
        sigma_w <- lavh1$implied$cov[[(g - 1) * 2 + 1]]
        mu_w <- lavh1$implied$mean[[(g - 1) * 2 + 1]]
        sigma_b <- lavh1$implied$cov[[(g - 1) * 2 + 2]]
        mu_b <- lavh1$implied$mean[[(g - 1) * 2 + 2]]
      }
      lp <- lavdata@Lp[[g]]

      info_g <-
        lav_mvnorm_cluster_information_expected_delta(
          lp = lp,
          delta = delta[[g]],
          mu_w = mu_w,
          sigma_w = sigma_w,
          mu_b = mu_b,
          sigma_b = sigma_b,
          sinv_method = "eigen"
        )
      info_group[[g]] <- fg * info_g
    } else {
      # compute information for this group
      if (lavmodel@estimator %in% c("DWLS", "ULS")) {
        # diagonal weight matrix
        delta2 <- sqrt(a1[[g]]) * delta[[g]]
        info_group[[g]] <- fg * crossprod(delta2)
      } else {
        # full weight matrix
        # if (lav_use_lavaanC()) {
        # # (i) use of m_crossprod with sparse matrix on the left:
        # # Info.group[[g]] <- fg * lavaanC::m_crossprod(Delta[[g]],
        # #                     lavaanC::m_prod(A1[[g]], Delta[[g]], "R"), "L")
        # #
        # # (ii) use of m_prod on transposed sparse first matrix,
        # #                                                 faster than (i):
        #   Info.group[[g]] <- fg * lavaanC::m_prod(t(Delta[[g]]),
        #                       lavaanC::m_prod(A1[[g]], Delta[[g]], "R"), "L")
        # } else {
          info_group[[g]] <-
            fg * lav_matrix_delta_A_delta(delta[[g]], a1[[g]])
        # }
      }
    }
  } # g

  # 4. assemble over groups
  information <- info_group[[1]]
  if (lavsamplestats@ngroups > 1) {
    for (g in 2:lavsamplestats@ngroups) {
      information <- information + info_group[[g]]
    }
  }

  # 5. augmented information?
  if (augmented) {
    information <-
      lav_model_information_augment_invert(
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
lav_model_information_expected_mlm <- function(lavmodel = NULL,
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
    a1[[g]] <- lav_mvnorm_h1_information_expected(
      sample_cov     = lavsamplestats@cov[[g]],
      sample_cov_inv = lavsamplestats@icov[[g]],
      x_idx          = lavsamplestats@x.idx[[g]]
    )
    # the same as GLS... (except for the N/N-1 scaling)
    if (lavmodel@group.w.free) {
      # unweight!!
      a <- exp(gw[g]) / lavsamplestats@nobs[[g]]
      # a <- exp(gw[g]) * lavsamplestats@ntotal / lavsamplestats@nobs[[g]]
      a1[[g]] <- lav_matrix_bdiag(matrix(a, 1, 1), a1[[g]])
    }
  }

  # compute Information per group
  info_group <- vector("list", length = lavsamplestats@ngroups)
  for (g in 1:lavsamplestats@ngroups) {
    fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
    # compute information for this group
    info_group[[g]] <- fg * (t(m_delta[[g]]) %*% a1[[g]] %*% m_delta[[g]])
  }

  # assemble over groups
  information <- info_group[[1]]
  if (lavsamplestats@ngroups > 1) {
    for (g in 2:lavsamplestats@ngroups) {
      information <- information + info_group[[g]]
    }
  }

  # augmented information?
  if (augmented) {
    information <-
      lav_model_information_augment_invert(
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


lav_model_information_observed <- function(lavmodel = NULL,
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

    a1 <- lav_model_h1_information_observed(
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
    # compute Information per group
    info_group <- vector("list", length = lavsamplestats@ngroups)
    for (g in 1:lavsamplestats@ngroups) {
      fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
      # compute information for this group
      if (lavmodel@estimator %in% c("DWLS", "ULS")) {
        # diagonal weight matrix
        delta2 <- sqrt(a1[[g]]) * delta[[g]]
        info_group[[g]] <- fg * crossprod(delta2)
      } else {
        # full weight matrix
        info_group[[g]] <-
          fg * lav_matrix_delta_A_delta(delta[[g]], a1[[g]])
      }
    }

    # assemble over groups
    information <- info_group[[1]]
    if (lavsamplestats@ngroups > 1) {
      for (g in 2:lavsamplestats@ngroups) {
        information <- information + info_group[[g]]
      }
    }
  }

  # augmented information?
  if (augmented) {
    information <-
      lav_model_information_augment_invert(
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
lav_model_information_firstorder <- function(lavmodel = NULL,           # nolint
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

  b0_group <- vector("list", lavsamplestats@ngroups)

  # 1. Delta
  delta <- lav_model_delta(lavmodel = lavmodel)

  # 2. H1 information
  b1 <- lav_model_h1_information_firstorder(
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

    # NOTE: UNSURE ABOUT THIS PART. WHAT IS THE ROLE OF fg?

    wt <- lavdata@weights[[g]]
    if (is.null(wt)) {
      fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
    } else {
      totalwt <- sum(unlist(lavdata@weights))
      fg <- sum(wt) / totalwt
    }

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # compute information for this group
    info_group[[g]] <- fg * b0_group[[g]]
  }

  # 4. assemble over groups
  information <- info_group[[1]]
  if (lavsamplestats@ngroups > 1) {
    for (g in 2:lavsamplestats@ngroups) {
      information <- information + info_group[[g]]
    }
  }

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
      lav_model_information_augment_invert(
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
lav_model_information_augment_invert <- function(lavmodel = NULL,   # nolint
                                                 information = NULL,
                                                 inverted = FALSE,
                                                 check_pd = FALSE,
                                                 use_ginv = FALSE,
                                                 rm_idx = integer(0L)) {
  npar <- nrow(information)
  is_augmented <- FALSE

  # handle constraints
  if (nrow(lavmodel@con.jac) > 0L) {
    h <- lavmodel@con.jac
    if (length(rm_idx) > 0L) {
      h <- h[, -rm_idx, drop = FALSE]
    }
    inactive_idx <- attr(h, "inactive.idx")
    lambda <- lavmodel@con.lambda # lagrangean coefs
    if (length(inactive_idx) > 0L) {
      h <- h[-inactive_idx, , drop = FALSE]
      lambda <- lambda[-inactive_idx]
    }
    if (nrow(h) > 0L) {
      is_augmented <- TRUE
      h0 <- matrix(0, nrow(h), nrow(h))
      h10 <- matrix(0, ncol(information), nrow(h))
      dl <- 2 * diag(lambda, nrow(h), nrow(h))
      # FIXME: better include inactive + slacks??
      # INFO <- information
      # or
      info <- information + crossprod(h)
      e3 <- rbind(
        cbind(info, h10, t(h)),
        cbind(t(h10), dl, h0),
        cbind(h, h0, h0)
      )
      information <- e3
    }
  } else if (lavmodel@ceq.simple.only) {
    h <- t(lav_matrix_orthogonal_complement(lavmodel@ceq.simple.K))
    if (length(rm_idx) > 0L) {
      h <- h[, -rm_idx, drop = FALSE]
    }
    if (nrow(h) > 0L) {
      is_augmented <- TRUE
      h0 <- matrix(0, nrow(h), nrow(h))
      h10 <- matrix(0, ncol(information), nrow(h))
      info <- information + crossprod(h)
      e2 <- rbind(
        cbind(info, t(h)),
        cbind(h, h0)
      )
      information <- e2
    }
  }

  if (check_pd) {
    eigvals <- eigen(information,
      symmetric = TRUE,
      only.values = TRUE
    )$values
    if (any(eigvals < -1 * .Machine$double.eps^(3 / 4))) {
      lav_msg_warn(gettext(
        "information matrix is not positive definite;
        the model may not be identified"))
    }
  }

  if (inverted) {
    if (is_augmented) {
      # note: default tol in MASS::ginv is sqrt(.Machine$double.eps)
      #       which seems a bit too conservative
      #       from 0.5-20, we changed this to .Machine$double.eps^(3/4)
      information <-
        try(
          MASS::ginv(information,
            tol = .Machine$double.eps^(3 / 4)
          )[1:npar,
            1:npar,
            drop = FALSE
          ],
          silent = TRUE
        )
    } else {
      if (use_ginv) {
        information <- try(
          MASS::ginv(information,
            tol = .Machine$double.eps^(3 / 4)
          ),
          silent = TRUE
        )
      } else {
        information <- try(solve(information), silent = TRUE)
      }
    }
  }

  # augmented/inverted information
  information
}
