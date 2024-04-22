# contributed by Ed Merkle (17 Jan 2013)
# WLS version contributed by Franz Classe (March 2024)
# (adapted for inclusion in lavaan by YR)


# YR 12 Feb 2013: small changes to match the results of lav_model_gradient
#                 in the multiple group case
# YR 30 May 2014: handle 1-variable case (fixing apply in lines 56, 62, 108)
# YR 05 Nov 2015: add remove.duplicated = TRUE, to cope with strucchange in
#                 case of simple equality constraints
# YR 19 Nov 2015: if constraints have been used, compute case-wise Lagrange
#                 multipliers, and define the scores as: SC + (t(R) lambda)
# YR 05 Feb 2016: catch conditional.x = TRUE: no support (for now), until
#                 we can use the generic 0.6 infrastructure for scores,
#                 including the missing-values case
# YR 16 Feb 2016: adapt to changed @Mp slot elements; add remove.empty.cases=
#                 argument
# YR 12 Mar 2024: make lintr (more) happy; include WLS code from Franz Classe
#                 move ML-specific code to lav_scores_ml() function

estfun.lavaan <- lavScores <- function(object, scaling = FALSE, # nolint
                                       ignore.constraints = FALSE,
                                       remove.duplicated = TRUE,
                                       remove.empty.cases = TRUE) {
  stopifnot(inherits(object, "lavaan"))

  # what if estimator is not ML or WLS?
  # avoid hard error (using stop); throw a warning, and return an empty matrix
  if (!object@Options$estimator %in% c("ML", "WLS")) {
    lav_msg_warn(gettext("scores only availlabe if estimator is ML"))
    return(matrix(0, 0, 0))
  }

  # check if conditional.x = TRUE
  if (object@Model@conditional.x) {
    lav_msg_stop(gettext("scores not available (yet) if conditional.x = TRUE"))
  }

  # shortcuts
  lavdata <- object@Data
  lavmodel <- object@Model
  lavsamplestats <- object@SampleStats
  lavoptions <- object@Options

  ## number variables/sample size
  # ntab <- unlist(lavsamplestats@nobs)
  ## change in 0.5-17: we keep the 'empty cases'
  ##                   and 'fill' in the scores at their 'case.idx'
  ##                   later, we remove the 'empty rows'
  # ntot <- max( object@Data@case.idx[[ object@Data@ngroups ]] )
  ntab <- unlist(lavdata@norig)
  ntot <- sum(ntab)
  npar <- lav_object_inspect_npar(object)

  if (object@Options$estimator == "ML") {
    moments <- fitted(object)
    score_matrix <- lav_scores_ml(
      ntab = ntab, ntot = ntot, npar = npar,
      moments = moments, lavdata = lavdata, lavsamplestats = lavsamplestats,
      lavmodel = lavmodel, lavoptions = lavoptions, scaling = scaling
    )
  } else if (object@Options$estimator == "WLS") {
    # check if ALL observed variables are ordered
    ov.names <- unlist(lavdata@ov.names)
    ov.idx <- which(lavdata@ov$name %in% ov.names)
    if (!all(lavdata@ov$type[ov.idx] == "ordered")) {
      lav_msg_stop(gettext(
        "WLS scores only available if all observed variables are ordered."))
    }

    # compute WLS scores
    score_matrix <- lav_scores_wls(
      ntab = ntab, ntot = ntot, npar = npar,
      lavdata = lavdata, lavsamplestats = lavsamplestats,
      lavmodel = lavmodel, lavoptions = lavoptions
    )
  } else {
    # should not happen
    lav_msg_fixme("this should not happen")
  }

  # handle empty rows
  if (remove.empty.cases) {
    # empty.idx <- which( apply(score_matrix, 1L,
    #                        function(x) sum(is.na(x))) == ncol(score_matrix) )
    empty.idx <- unlist(lapply(lavdata@Mp, "[[", "empty.idx"))
    if (length(empty.idx) > 0L) {
      score_matrix <- score_matrix[-empty.idx, , drop = FALSE]
    }
  }

  # provide column names
  colnames(score_matrix) <- names(lav_object_inspect_coef(object,
    type = "free", add.labels = TRUE
  ))

  # handle general constraints, so that the sum of the columns equals zero
  if (!ignore.constraints &&
    sum(
      lavmodel@ceq.linear.idx, lavmodel@ceq.nonlinear.idx,
      lavmodel@cin.linear.idx, lavmodel@cin.nonlinear.idx
    ) > 0) {
    r_matrix <- object@Model@con.jac[, ]
    pre <- lav_constraints_lambda_pre(object)
    # LAMBDA <- -1 * t(pre %*% t(score_matrix))
    # RLAMBDA <- t(t(r_matrix) %*% t(LAMBDA))
    score_matrix <- score_matrix - t(t(r_matrix) %*% pre %*% t(score_matrix))
  }

  # handle simple equality constraints
  if (remove.duplicated && lavmodel@eq.constraints) {
    simple.flag <- lav_constraints_check_simple(lavmodel)
    if (simple.flag) {
      k_matrix <- lav_constraints_R2K(lavmodel)
      score_matrix <- score_matrix %*% k_matrix
    } else {
      lav_msg_warn(gettext(
        "remove.duplicated is TRUE, but equality constraints do not appear
        to be simple; returning full scores"))
    }
  }

  score_matrix
}

lav_scores_ml <- function(ntab = 0L,
                          ntot = 0L,
                          npar = 0L,
                          moments = NULL,
                          lavdata = NULL,
                          lavsamplestats = NULL,
                          lavmodel = NULL,
                          lavoptions = NULL,
                          scaling = FALSE) {
  score_matrix <- matrix(NA, ntot, npar)

  # Delta matrix
  Delta <- computeDelta(lavmodel = lavmodel)

  # rename moments
  moments.groups <- moments

  for (g in 1:lavsamplestats@ngroups) {
    if (lavsamplestats@ngroups > 1) {
      moments <- moments.groups[[g]]
    }
    sigma_hat <- moments$cov

    if (lavoptions$likelihood == "wishart") {
      nobs1 <- lavsamplestats@nobs[[g]] / (lavsamplestats@nobs[[g]] - 1)
    } else {
      nobs1 <- 1
    }

    if (!lavsamplestats@missing.flag) { # complete data
      # if(lavmodel@meanstructure) { # mean structure
      nvar <- ncol(lavsamplestats@cov[[g]])
      mu_hat <- moments$mean
      X <- lavdata@X[[g]]
      sigma_inv <- inv.chol(sigma_hat, logdet = FALSE)
      group.w <- (unlist(lavsamplestats@nobs) / lavsamplestats@ntotal)

      J <- matrix(1, 1L, ntab[g]) ## FIXME: needed? better maybe rowSums/colSums?
      J2 <- matrix(1, nvar, nvar)
      diag(J2) <- 0.5

      if (lavmodel@meanstructure) {
        ## scores_h1 (H1 = saturated model)
        mean.diff <- t(t(X) - mu_hat %*% J)
        dx_mu <- -1 * mean.diff %*% sigma_inv
        dx_sigma <- t(matrix(apply(
          mean.diff, 1L,
          function(x) {
            lav_matrix_vech(-J2 *
              (sigma_inv %*% (tcrossprod(x) * nobs1 - sigma_hat) %*% sigma_inv))
          }
        ), ncol = nrow(mean.diff)))

        scores_h1 <- cbind(dx_mu, dx_sigma)
      } else {
        mean.diff <- t(t(X) - lavsamplestats@mean[[g]] %*% J)
        dx_sigma <- t(matrix(apply(
          mean.diff, 1L,
          function(x) {
            lav_matrix_vech(-J2 *
              (sigma_inv %*% (tcrossprod(x) * nobs1 - sigma_hat) %*% sigma_inv))
          }
        ), ncol = nrow(mean.diff)))
        scores_h1 <- dx_sigma
      }
      ## FIXME? Seems like we would need group.w even in the
      ##        complete-data case:
      ## if(scaling){
      ##  scores_h1 <- group.w[g] * scores_h1
      ## }

      # } else {
      #  ## no mean structure
      #  stop("Score calculation with no mean structure is not implemented.")
      # }
    } else { # incomplete data
      nsub <- ntab[g]
      M <- lavsamplestats@missing[[g]]
      Mp <- lavdata@Mp[[g]]
      # pat.idx <- match(MP1$id, MP1$order)
      group.w <- (unlist(lavsamplestats@nobs) / lavsamplestats@ntotal)

      mu_hat <- moments$mean
      nvar <- ncol(lavsamplestats@cov[[g]])
      score.sigma <- matrix(0, nsub, nvar * (nvar + 1) / 2)
      score.mu <- matrix(0, nsub, nvar)

      for (p in seq_along(length(M))) {
        ## Data
        # X <- M[[p]][["X"]]
        case.idx <- Mp$case.idx[[p]]
        var.idx <- M[[p]][["var.idx"]]
        X <- lavdata@X[[g]][case.idx, var.idx, drop = FALSE]
        nobs <- M[[p]][["freq"]]
        ## Which unique entries of covariance matrix are estimated?
        ## (Used to keep track of scores in score.sigma)
        var.idx.mat <- tcrossprod(var.idx)
        sigma.idx <-
          which(var.idx.mat[lower.tri(var.idx.mat, diag = TRUE)] == 1)

        J <- matrix(1, 1L, nobs) # [var.idx]
        J2 <- matrix(1, nvar, nvar)[var.idx, var.idx, drop = FALSE]
        diag(J2) <- 0.5
        sigma_inv <- inv.chol(sigma_hat[var.idx, var.idx, drop = FALSE],
          logdet = FALSE
        )
        Mu <- mu_hat[var.idx]
        mean.diff <- t(t(X) - Mu %*% J)

        ## Scores for missing pattern p within group g
        score.mu[case.idx, var.idx] <- -1 * mean.diff %*% sigma_inv
        score.sigma[case.idx, sigma.idx] <- t(matrix(apply(
          mean.diff, 1L,
          function(x) {
            lav_matrix_vech(-J2 *
              (sigma_inv %*% (tcrossprod(x) -
                sigma_hat[var.idx, var.idx, drop = FALSE]) %*% sigma_inv))
          }
        ), ncol = nrow(mean.diff)))
      }

      scores_h1 <- cbind(score.mu, score.sigma)
      if (scaling) {
        scores_h1 <- group.w[g] * scores_h1
      }
    } # missing

    # if(lavmodel@eq.constraints) {
    #    Delta <- Delta %*% lavmodel@eq.constraints.K
    #    #x <- as.numeric(lavmodel@eq.constraints.K %*% x) +
    #    #                lavmodel@eq.constraints.k0
    # }
    wi <- lavdata@case.idx[[g]]
    score_matrix[wi, ] <- -scores_h1 %*% Delta[[g]]
    if (scaling) {
      score_matrix[wi, ] <- (-1 / ntot) * score_matrix[wi, ]
    }
  } # g

  score_matrix
}

# this function is based on code originally written by Franz Classe (Munich)
lav_scores_wls <- function(ntab = 0L,
                           ntot = 0L,
                           npar = 0L,
                           lavdata = NULL,
                           lavsamplestats = NULL,
                           lavmodel = NULL,
                           lavoptions = NULL) {
  # internal function
  doDummySingleVar <- function(X, lv, ntot, num) {
    Xd <- matrix(NA, nrow = ntot, ncol = lv[num] - 1)
    x <- X[, num]
    minx <- min(x)
    categ <- minx - 1
    v <- 1
    while (categ < lv[num] - 1) {
      categ <- categ + 1
      Xd[, v] <- ifelse(x > categ, 1, 0)
      v <- v + 1
    }

    Xd
  }

  # containere for scores
  score_matrix <- matrix(NA, ntot, npar)

  # Delta matrix
  Delta <- computeDelta(lavmodel = lavmodel)

  # shortcuts
  lv <- lavdata@ov[["nlev"]]

  for (g in 1:lavsamplestats@ngroups) {
    nvar <- ncol(lavsamplestats@cov[[g]])
    X <- lavdata@X[[g]]

    # convert categorical data to dummy variables
    # FIXME: skip continuous variables
    Xd <- do.call(
      cbind,
      lapply(
        1:nvar,
        function(i) doDummySingleVar(X, lv, ntot, i)
      )
    )

    # e1
    musd <- colMeans(Xd)
    e1 <- t(t(Xd) - musd)

    # e2
    mus <- colMeans(X)
    y_minus_mu <- t(apply(X, 1L, function(x) x - mus))
    s_vech <- t(apply(y_minus_mu, 1L, function(i) {
      lavaan::lav_matrix_vech(tcrossprod(i), diagonal = FALSE)
    })) # s=c( (y1-mu1)(y2-mu2)....
    sigma <- colMeans(s_vech)
    e2 <- t(apply(s_vech, 1L, function(x) x - sigma))

    # e
    e <- cbind(e1, e2)

    # weight matrix
    W <- lavsamplestats@WLS.V[[g]]

    # combine matrices
    wi <- lavdata@case.idx[[g]]
    score_matrix[wi, ] <- t(t(Delta[[g]]) %*% W %*% t(e))
  } # g

  score_matrix
}
