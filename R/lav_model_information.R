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
                                  Delta = NULL,
                                  lavcache = NULL,
                                  lavoptions = NULL,
                                  extra = FALSE,
                                  augmented = FALSE,
                                  inverted = FALSE,
                                  use.ginv = FALSE) {
  if (.hasSlot(lavmodel, "estimator")) {
    estimator <- lavmodel@estimator
  } else {
    estmator <- lavoptions$estimator
  }
  information <- lavoptions$information[1] # ALWAYS used the first one
  # called can control it

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
      group.weight <- FALSE
    } else {
      group.weight <- TRUE
    }
    E <- lav_model_information_observed(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, group.weight = group.weight,
      lavoptions = lavoptions, extra = extra,
      augmented = augmented, inverted = inverted, use.ginv = use.ginv
    )
  } else if (information == "expected") {
    E <- lav_model_information_expected(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions, extra = extra,
      augmented = augmented, inverted = inverted, use.ginv = use.ginv
    )
  } else if (information == "first.order") {
    E <- lav_model_information_firstorder(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = lavoptions, # extra = extra,
      check.pd = FALSE,
      augmented = augmented, inverted = inverted, use.ginv = use.ginv
    )
  }

  # information, augmented information, or inverted information
  E
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
                                           Delta = NULL,
                                           lavcache = NULL,
                                           extra = FALSE,
                                           augmented = FALSE,
                                           inverted = FALSE,
                                           use.ginv = FALSE) {
  if (inverted) {
    augmented <- TRUE
  }

  # 1. Delta
  if (is.null(Delta)) {
    Delta <- computeDelta(lavmodel = lavmodel)
  }


  # 2. H1 information (single level)
  if (lavdata@nlevels == 1L) {
    A1 <- lav_model_h1_information_expected(
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
  Info.group <- vector("list", length = lavsamplestats@ngroups)
  for (g in 1:lavsamplestats@ngroups) {
    # note LISREL documentation suggests (Ng - 1) instead of Ng...
    fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal

    # multilevel
    if (lavdata@nlevels > 1L) {
      # here, we assume only 2 levels, at [[1]] and [[2]]
      if (lavoptions$h1.information[1] == "structured") {
        Sigma.W <- lavimplied$cov[[(g - 1) * 2 + 1]]
        Mu.W <- lavimplied$mean[[(g - 1) * 2 + 1]]
        Sigma.B <- lavimplied$cov[[(g - 1) * 2 + 2]]
        Mu.B <- lavimplied$mean[[(g - 1) * 2 + 2]]
      } else {
        Sigma.W <- lavh1$implied$cov[[(g - 1) * 2 + 1]]
        Mu.W <- lavh1$implied$mean[[(g - 1) * 2 + 1]]
        Sigma.B <- lavh1$implied$cov[[(g - 1) * 2 + 2]]
        Mu.B <- lavh1$implied$mean[[(g - 1) * 2 + 2]]
      }
      Lp <- lavdata@Lp[[g]]

      Info.g <-
        lav_mvnorm_cluster_information_expected_delta(
          Lp = Lp,
          Delta = Delta[[g]],
          Mu.W = Mu.W,
          Sigma.W = Sigma.W,
          Mu.B = Mu.B,
          Sigma.B = Sigma.B,
          Sinv.method = "eigen"
        )
      Info.group[[g]] <- fg * Info.g
    } else {
      # compute information for this group
      if (lavmodel@estimator %in% c("DWLS", "ULS")) {
        # diagonal weight matrix
        Delta2 <- sqrt(A1[[g]]) * Delta[[g]]
        Info.group[[g]] <- fg * crossprod(Delta2)
      } else {
        # full weight matrix
        Info.group[[g]] <-
          fg * (crossprod(Delta[[g]], A1[[g]]) %*% Delta[[g]])
      }
    }
  } # g

  # 4. assemble over groups
  Information <- Info.group[[1]]
  if (lavsamplestats@ngroups > 1) {
    for (g in 2:lavsamplestats@ngroups) {
      Information <- Information + Info.group[[g]]
    }
  }

  # 5. augmented information?
  if (augmented) {
    Information <-
      lav_model_information_augment_invert(
        lavmodel = lavmodel,
        information = Information,
        inverted = inverted,
        use.ginv = use.ginv
      )
  }

  if (extra) {
    attr(Information, "Delta") <- Delta
    attr(Information, "WLS.V") <- A1 # unweighted
  }

  # possibly augmented/inverted
  Information
}

# only for Mplus MLM
lav_model_information_expected_MLM <- function(lavmodel = NULL,
                                               lavsamplestats = NULL,
                                               Delta = NULL,
                                               extra = FALSE,
                                               augmented = FALSE,
                                               inverted = FALSE,
                                               use.ginv = FALSE) {
  if (inverted) {
    augmented <- TRUE
  }

  if (is.null(Delta)) {
    Delta <- computeDelta(lavmodel = lavmodel)
  }

  # compute A1
  A1 <- vector("list", length = lavsamplestats@ngroups)
  if (lavmodel@group.w.free) {
    GW <- unlist(computeGW(lavmodel = lavmodel))
  }
  for (g in 1:lavsamplestats@ngroups) {
    A1[[g]] <- lav_mvnorm_h1_information_expected(
      sample.cov     = lavsamplestats@cov[[g]],
      sample.cov.inv = lavsamplestats@icov[[g]],
      x.idx          = lavsamplestats@x.idx[[g]]
    )
    # the same as GLS... (except for the N/N-1 scaling)
    if (lavmodel@group.w.free) {
      # unweight!!
      a <- exp(GW[g]) / lavsamplestats@nobs[[g]]
      # a <- exp(GW[g]) * lavsamplestats@ntotal / lavsamplestats@nobs[[g]]
      A1[[g]] <- lav_matrix_bdiag(matrix(a, 1, 1), A1[[g]])
    }
  }

  # compute Information per group
  Info.group <- vector("list", length = lavsamplestats@ngroups)
  for (g in 1:lavsamplestats@ngroups) {
    fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
    # compute information for this group
    Info.group[[g]] <- fg * (t(Delta[[g]]) %*% A1[[g]] %*% Delta[[g]])
  }

  # assemble over groups
  Information <- Info.group[[1]]
  if (lavsamplestats@ngroups > 1) {
    for (g in 2:lavsamplestats@ngroups) {
      Information <- Information + Info.group[[g]]
    }
  }

  # augmented information?
  if (augmented) {
    Information <-
      lav_model_information_augment_invert(
        lavmodel = lavmodel,
        information = Information,
        inverted = inverted,
        use.ginv = use.ginv
      )
  }

  if (extra) {
    attr(Information, "Delta") <- Delta
    attr(Information, "WLS.V") <- A1 # unweighted
  }

  Information
}


lav_model_information_observed <- function(lavmodel = NULL,
                                           lavsamplestats = NULL,
                                           lavdata = NULL,
                                           lavimplied = NULL,
                                           lavh1 = NULL,
                                           lavcache = NULL,
                                           lavoptions = NULL,
                                           extra = FALSE,
                                           group.weight = TRUE,
                                           augmented = FALSE,
                                           inverted = FALSE,
                                           use.ginv = FALSE) {
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
    observed.information <- "h1"
  } else {
    observed.information <- "hessian"
  }

  # HESSIAN based
  if (observed.information == "hessian") {
    Hessian <- lav_model_hessian(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      lavcache = lavcache,
      group.weight = group.weight,
      ceq.simple = FALSE
    )

    # NOTE! What is the relationship between the Hessian of the objective
    # function, and the `information' matrix (unit or total)

    # 1. in lavaan, we ALWAYS minimize, so the Hessian is already pos def
    # 2. currently, all estimators give unit information, except MML and PML
    #    so, no need to divide by N
    Information <- Hessian

    # divide by 'N' for MML and PML
    if (lavmodel@estimator == "PML" || lavmodel@estimator == "MML") {
      Information <- Information / lavsamplestats@ntotal
      # HJ: Does this need to be divided by sum of weights instead?
    }

    # if multilevel, we should divide by 'J', the number of clusters
    if (lavdata@nlevels > 1L) {
      NC <- 0
      for (g in 1:lavsamplestats@ngroups) {
        NC <- NC + lavdata@Lp[[g]]$nclusters[[2]]
      }
      Information <- Information * lavsamplestats@ntotal / NC
    }
  }

  # using 'observed h1 information'
  # we need DELTA and 'WLS.V' (=A1)

  if (observed.information == "h1" || extra) {
    # 1. Delta
    Delta <- computeDelta(lavmodel = lavmodel)

    # 2. H1 information

    A1 <- lav_model_h1_information_observed(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavoptions = lavoptions,
      lavimplied = lavimplied,
      lavh1 = lavh1,
      lavcache = lavcache
    )
  }

  if (observed.information == "h1") {
    # compute Information per group
    Info.group <- vector("list", length = lavsamplestats@ngroups)
    for (g in 1:lavsamplestats@ngroups) {
      fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
      # compute information for this group
      if (lavmodel@estimator %in% c("DWLS", "ULS")) {
        # diagonal weight matrix
        Delta2 <- sqrt(A1[[g]]) * Delta[[g]]
        Info.group[[g]] <- fg * crossprod(Delta2)
      } else {
        # full weight matrix
        Info.group[[g]] <-
          fg * (crossprod(Delta[[g]], A1[[g]]) %*% Delta[[g]])
      }
    }

    # assemble over groups
    Information <- Info.group[[1]]
    if (lavsamplestats@ngroups > 1) {
      for (g in 2:lavsamplestats@ngroups) {
        Information <- Information + Info.group[[g]]
      }
    }
  }

  # augmented information?
  if (augmented) {
    Information <-
      lav_model_information_augment_invert(
        lavmodel = lavmodel,
        information = Information,
        inverted = inverted,
        use.ginv = use.ginv
      )
  }

  if (extra) {
    attr(Information, "Delta") <- Delta
    attr(Information, "WLS.V") <- A1
  }

  Information
}

# outer product of the case-wise scores (gradients)
# HJ 18/10/23: Need to divide sum of crossproduct of individual log-likelihoods
# by sum of weights rather than sample size.
lav_model_information_firstorder <- function(lavmodel = NULL,
                                             lavsamplestats = NULL,
                                             lavdata = NULL,
                                             lavimplied = NULL,
                                             lavh1 = NULL,
                                             lavcache = NULL,
                                             lavoptions = NULL,
                                             check.pd = FALSE,
                                             extra = FALSE,
                                             augmented = FALSE,
                                             inverted = FALSE,
                                             use.ginv = FALSE) {
  if (!lavmodel@estimator %in% c("ML", "PML")) {
    lav_msg_stop(gettext(
      "information = \"first.order\" not available for estimator"),
      sQuote(lavmodel@estimator))
  }

  if (inverted) {
    augmented <- TRUE
  }

  B0.group <- vector("list", lavsamplestats@ngroups)

  # 1. Delta
  Delta <- computeDelta(lavmodel = lavmodel)

  # 2. H1 information
  B1 <- lav_model_h1_information_firstorder(
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    lavimplied = lavimplied,
    lavh1 = lavh1,
    lavcache = lavcache
  )

  # 3. compute Information per group
  Info.group <- vector("list", length = lavsamplestats@ngroups)
  for (g in 1:lavsamplestats@ngroups) {
    # unweighted (needed in lav_test?)
    B0.group[[g]] <- t(Delta[[g]]) %*% B1[[g]] %*% Delta[[g]]

    # >>>>>>>> HJ/MK PML CODE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # NOTE: UNSURE ABOUT THIS PART. WHAT IS THE ROLE OF fg?

    if (.hasSlot(lavdata, "weights")) {
      wt <- lavdata@weights[[g]]
    } else { # pre-0.6 object
      wt <- NULL
    }
    if (is.null(wt)) {
      fg <- lavsamplestats@nobs[[g]] / lavsamplestats@ntotal
    } else {
      totalwt <- sum(unlist(lavdata@weights))
      fg <- sum(wt) / totalwt
    }

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # compute information for this group
    Info.group[[g]] <- fg * B0.group[[g]]
  }

  # 4. assemble over groups
  Information <- Info.group[[1]]
  if (lavsamplestats@ngroups > 1) {
    for (g in 2:lavsamplestats@ngroups) {
      Information <- Information + Info.group[[g]]
    }
  }

  # NOTE: for MML and PML, we get 'total' information (instead of unit) divide
  # by 'N' for MML and PML. For weighted sample, use the sum of weights
  # instead of sample size

  # >>>>>>>> HJ/MK PML CODE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  if (lavmodel@estimator == "PML" || lavmodel@estimator == "MML") {
    if (!.hasSlot(lavdata, "sampling.weights") ||
      length(lavdata@sampling.weights) == 0) {
      the_N <- lavsamplestats@ntotal
    } else {
      the_N <- sum(unlist(lavdata@weights))
    }
    Information <- Information / the_N
    for (g in 1:lavsamplestats@ngroups) {
      B0.group[[g]] <- B0.group[[g]] / the_N
    }
  }
  # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  # augmented information?
  if (augmented) {
    Information <-
      lav_model_information_augment_invert(
        lavmodel = lavmodel,
        information = Information,
        check.pd = check.pd,
        inverted = inverted,
        use.ginv = use.ginv
      )
  }

  if (extra) {
    attr(Information, "B0.group") <- B0.group
    attr(Information, "Delta") <- Delta
    attr(Information, "WLS.V") <- B1
  }

  Information
}


# create augmented information matrix (if needed), and take the inverse
# (if inverted = TRUE), returning only the [1:npar, 1:npar] elements
#
# rm.idx is used by lav_sam_step2_se; it used when the structural model
# contains more parameters than the joint model; therefore, information
# will be 'too small', and we need to remove some columns in H
lav_model_information_augment_invert <- function(lavmodel = NULL,
                                                 information = NULL,
                                                 inverted = FALSE,
                                                 check.pd = FALSE,
                                                 use.ginv = FALSE,
                                                 rm.idx = integer(0L)) {
  npar <- nrow(information)
  is.augmented <- FALSE

  # handle constraints
  if (nrow(lavmodel@con.jac) > 0L) {
    H <- lavmodel@con.jac
    if (length(rm.idx) > 0L) {
      H <- H[, -rm.idx, drop = FALSE]
    }
    inactive.idx <- attr(H, "inactive.idx")
    lambda <- lavmodel@con.lambda # lagrangean coefs
    if (length(inactive.idx) > 0L) {
      H <- H[-inactive.idx, , drop = FALSE]
      lambda <- lambda[-inactive.idx]
    }
    if (nrow(H) > 0L) {
      is.augmented <- TRUE
      H0 <- matrix(0, nrow(H), nrow(H))
      H10 <- matrix(0, ncol(information), nrow(H))
      DL <- 2 * diag(lambda, nrow(H), nrow(H))
      # FIXME: better include inactive + slacks??
      # INFO <- information
      # or
      INFO <- information + crossprod(H)
      E3 <- rbind(
        cbind(INFO, H10, t(H)),
        cbind(t(H10), DL, H0),
        cbind(H, H0, H0)
      )
      information <- E3
    }
  } else if (.hasSlot(lavmodel, "ceq.simple.only") &&
    lavmodel@ceq.simple.only) {
    H <- t(lav_matrix_orthogonal_complement(lavmodel@ceq.simple.K))
    if (length(rm.idx) > 0L) {
      H <- H[, -rm.idx, drop = FALSE]
    }
    if (nrow(H) > 0L) {
      is.augmented <- TRUE
      H0 <- matrix(0, nrow(H), nrow(H))
      H10 <- matrix(0, ncol(information), nrow(H))
      INFO <- information + crossprod(H)
      E2 <- rbind(
        cbind(INFO, t(H)),
        cbind(H, H0)
      )
      information <- E2
    }
  }

  if (check.pd) {
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
    if (is.augmented) {
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
      if (use.ginv) {
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

lav_model_information_expected_2l <- function(lavmodel = NULL,
                                              lavsamplestats = NULL,
                                              lavdata = NULL,
                                              lavoptions = NULL,
                                              lavimplied = NULL,
                                              lavh1 = NULL,
                                              g = 1L) {
  # see Yuan & Bentler (2002), p.549 top line
  # I.j = nj. Delta.mu' sigma.j.inv +
  #       Delta.sigma.j' W.j Delta.sigma.j +
  #       (nj-1) Delta.sigma.w' W.w Delta.sigma.w
  #
  # where
  # - sigma.j = sigma.w + n.j * sigma.b
  # - W.w = 1/2 * D'(sigma.w.inv %x% sigma.w.inv) D
  # - W.j = 1/2 * D'(sigma.j.inv %x% sigma.j.inv) D
}
