# find ALL (eligible) model-implied instrumental variables (no pruning yet)
#
# two algorithms:
#   1) using 'total effects' of errors/disturbances, based on Bollen & Curran
#      (2004) [algorithm = "bc2004"]
#   2) treating errors/disturbances as latent variables (as in the
#      comprehensive RAM model); similar to miivs() in MIIVsem package
#      [algorithm = "miivsem"] (the default)
#
# YR 29 Dec 2025 - first version

lav_model_find_iv <- function(lavobject = NULL, lavmodel = NULL,
                              lavpta = NULL, algorithm = "miivsem",
                              output = "list", drop.list.single.group = FALSE) {
  # check output
  output <- tolower(output)
  stopifnot(output %in% c("list", "table"))

  # lavobject or components?
  if (!is.null(lavobject)) {
    stopifnot(inherits(lavobject, "lavaan"))
    lavpta <- lavobject@pta
    lavmodel <- lavobject@Model
  }

  # sanity checks
  nblocks <- lavpta$nblocks
  for (b in seq_len(nblocks)) {
    lv.idx <- lavpta$vidx$lv.regular[[b]]
    lv.marker <- lavpta$vnames$lv.marker[[b]]
    if (length(lv.idx) > 0L) {
      # do have 'clear' marker/scaling indicators?
      empty.idx <- which(nchar(lv.marker) == 0L)
      if (length(empty.idx) > 0L) {
        tmp_string <- paste(names(lv.marker)[empty.idx], collapse = " ")
        lav_msg_stop(gettextf("no clear marker/scaling indicator found for
                               factor(s): %s", tmp_string))
      }
      # marker/scaling indicator cannot be a dependent variable in an ~ equation
      badmarker.idx <- which(lv.marker %in% lavpta$vnames$eqs.y[[b]])
      if (length(badmarker.idx) > 0L) {
        tmp_string <- paste(lv.marker[badmarker.idx], collapse = " ")
        lav_msg_stop(gettextf("marker/scaling indicator cannot be a dependent
                               variable in a regression: %s", tmp_string))
      }
    }
  }

  algorithm <- tolower(algorithm)
  if (algorithm == "bc2004") {
    iv_list <- lav_model_find_iv_bc2004(lavmodel = lavmodel, lavpta = lavpta)
  } else if (algorithm == "miivsem") {
    iv_list <- lav_model_find_iv_miivsem(lavmodel = lavmodel, lavpta = lavpta)
  }

  if (output == "table") {
    table <- vector("list", length = nblocks)
    for (b in seq_len(nblocks)) {
      eqs <- iv_list[[b]]
      lhs <- sapply(eqs, "[[", "lhs")
      rhs <- sapply(lapply(eqs, "[[", "rhs"), paste, collapse = ", ")
      lhs_new <- sapply(eqs, "[[", "lhs_new")
      rhs_new <- sapply(lapply(eqs, "[[", "rhs_new"), paste, collapse = ", ")
      miiv <- sapply(lapply(eqs, "[[", "miiv"), paste, collapse = ", ")
      table[[b]] <- data.frame(
        lhs = lhs, rhs = rhs,
        lhs_new = lhs_new, rhs_new = rhs_new, miiv = miiv
      )
      class(table[[b]]) <- c("lavaan.data.frame", "data.frame")
    }
    out <- table
  } else {
    out <- iv_list
  }

  if (nblocks == 1L && drop.list.single.group) {
    out <- out[[1]]
  }

  out
}


# algorithm 1:
# loosely based on Bollen & Bauer (2004), but with support for higher-order
# factors
lav_model_find_iv_bc2004 <- function(lavmodel = NULL, lavpta = NULL) {
  # check representation
  if (lavmodel@representation != "LISREL") {
    lav_msg_stop(gettext(
      "this function only works with LISREL representation",
      " (for now)"
    ))
  }

  # create model matrices with 'user/partable' entries
  glist <- lavmodel@GLIST
  for (mm in seq_along(glist)) {
    dimnames(glist[[mm]]) <- lavmodel@dimNames[[mm]]
    m.idx <- lavmodel@m.user.idx[[mm]]
    x.idx <- lavmodel@x.user.idx[[mm]]
    glist[[mm]][, ] <- 0.0
    glist[[mm]][m.idx] <- x.idx
  }

  # number of blocks
  nblocks <- lavpta$nblocks
  lambda.idx <- which(names(glist) == "lambda")
  beta.idx <- which(names(glist) == "beta")
  psi.idx <- which(names(glist) == "psi")
  theta.idx <- which(names(glist) == "theta")

  # repeat for every block
  iv_list <- vector("list", length = nblocks)
  for (b in seq_len(nblocks)) {
    # extract information to create Bollen & Bauer (2004) matrices
    ov.names <- lavpta$vnames$ov[[b]]
    nvar <- length(ov.names)
    lv.names <- lavpta$vnames$lv.regular[[b]]
    lv.marker <- lavpta$vnames$lv.marker[[b]]
    lv.marker.orig <- lavpta$vnames$lv.marker[[b]]
    lv.ho <- lavpta$vnames$lv.ho[[b]]


    lv.idx <- lavpta$vidx$lv.regular[[b]]
    lv.x.idx <- lavpta$vidx$lv.x[[b]]
    lv.marker.idx <- lavpta$vidx$lv.marker[[b]]

    # keep track of higher order factors
    lv.ho.idx <- which(is.na(lv.marker.idx)) # same as lavpta$vidx$lv.ho
    if (length(lv.ho.idx) > 0L) {
      lv_in_marker <- lv.marker[lv.marker %in% lv.names]
      while (length(lv_in_marker) > 0L) {
        new_markers <- unname(lv.marker[lv_in_marker])
        lv.marker[match(lv_in_marker, lv.marker)] <- new_markers
        lv_in_marker <- lv.marker[lv.marker %in% lv.names]
      }
      lv.marker.idx <- match(lv.marker, ov.names)
    }

    # model matrices for this block
    lambda <- glist[[lambda.idx[b]]]
    beta <- glist[[beta.idx[b]]]
    theta <- glist[[theta.idx[b]]]
    psi <- glist[[psi.idx[b]]]

    # binary model matrices: nonzero = 1, zero = 0
    if (length(lv.idx) == 0L) {
      lambda_bin <- diag(nrow = ncol(lambda))
    } else {
      lambda_bin <- (lambda != 0) * 1L
    }
    if (!is.null(beta)) {
      beta_bin <- (beta != 0) * 1L
      ibinv <- solve(diag(nrow(beta_bin)) - beta_bin)
      ibinv_bin <- (ibinv != 0) * 1L
    } else {
      beta_bin <- matrix(0, nrow = 0L, ncol = ncol(lambda))
      ibinv_bin <- diag(nrow = ncol(lambda_bin))
    }

    # construct pred
    if (length(lv.idx) > 0L) {
      reg_bin <- rbind(lambda_bin, beta_bin)
      pred_orig <- rbind(lambda, beta)
    } else {
      reg_bin <- beta_bin
      pred_orig <- beta
    }
    pred.idx <- numeric(ncol(reg_bin))
    if (length(lv.marker.idx) > 0L) {
      pred.idx[lv.idx] <- lv.marker.idx
    }
    if (length(lavmodel@ov.x.dummy.ov.idx[[b]]) > 0L) {
      pred.idx[lavmodel@ov.x.dummy.lv.idx[[b]]] <-
        lavmodel@ov.x.dummy.ov.idx[[b]]
    }
    if (length(lavmodel@ov.y.dummy.ov.idx[[b]]) > 0L) {
      pred.idx[lavmodel@ov.y.dummy.lv.idx[[b]]] <-
        lavmodel@ov.y.dummy.ov.idx[[b]]
    }
    tmp <- t(t(reg_bin) * pred.idx)
    # remove scaling '1' for markers (if any)
    if (length(lv.marker.idx) > 0L) {
      row.idx <- match(lv.marker.orig, rownames(tmp))
      col.idx <- match(names(lv.marker.orig), colnames(tmp))
      tmp[cbind(row.idx, col.idx)] <- 0L
    }
    # keep only the 'ov' part
    pred <- tmp[seq_along(ov.names), , drop = FALSE]
    # replace marker rows, with beta entries of their corresponding lv's
    if (length(lv.idx) > 0L && !is.null(beta)) {
      tmp_lv <- tmp[length(ov.names) + lv.idx, , drop = FALSE]
      orig_lv <- pred_orig[length(ov.names) + lv.idx, , drop = FALSE]
      # tmp[lv.marker.idx, ] <- tmp[length(ov.names) + lv.idx, ]
      # careful, lv.marker.idx may contain the same element multiple times!
      zerolv.idx <- which(apply(tmp_lv, 1, function(x) all(x == 0)))
      if (length(zerolv.idx) > 0L) {
        tmp_lv <- tmp_lv[-zerolv.idx, , drop = FALSE]
        orig_lv <- orig_lv[-zerolv.idx, , drop = FALSE]
      }
      r.idx <- match(lv.marker[rownames(tmp_lv)], rownames(pred))
      pred[r.idx, ] <- tmp_lv
      pred_orig[r.idx, ] <- orig_lv
      rownames(pred_orig)[r.idx] <- rownames(tmp_lv)
      pred_orig <- pred_orig[1:nvar, , drop = FALSE]
    }
    # non-zero entries of 'pred' contain the marker idx of the latent
    # predictors

    # construct comp (composite error term)
    # - for each dependent: the corresponding element in Theta/Psi (y_res)
    # - for each (latent) predictor: the corresponding *marker* element in
    #   Theta (x_res)
    # - for the markers: the disturbance (psi) of the corresponding
    #   lv's (included in x_res)
    # - for higher-order factors: the disturbance (psi) of their 'marker'
    #   factor
    y_res <- diag(theta)
    if (length(lv.marker.idx) > 0L) {
      x_res <- pred
      # x_res[pred != 0] <- y_res[as.vector(pred[pred != 0])]
      x_res[pred != 0] <- y_res[pred]
      x_res[cbind(lv.marker.idx, seq_along(lv.idx))] <- diag(psi)[lv.idx]
      if (length(lv.ho.idx) > 0L) {
        # higher-order factors: what are their indicators?
        # add diag(psi) value in this row
        x_res <- t(apply(x_res, 1, function(x) {
          res <- x
          hof <- x[lv.ho.idx]
          hof_nonzero <- hof[hof != 0]
          if (length(hof_nonzero) > 0L) {
            # keep going until we only have first-order factors (or ovs)
            target <- unname(lv.marker.orig[names(hof_nonzero)])
            ho_in_target <- target[target %in% lv.ho]
            while (length(ho_in_target) > 0L) {
              new <- unname(lv.marker.orig[ho_in_target])
              target <- c(target, new)
              ho_in_target <- new[new %in% lv.ho]
            }

            if (length(target)) {
              res[target] <- diag(psi)[target]
            }
          }
          res
        }))
      }
      # comp <- cbind(y_res, x_res)[1:nvar, ]
      comp <- cbind(y_res, x_res)
    } else {
      # no latent variables
      comp <- cbind(y_res)
    }

    # construct total: total effects of errors/disturbances on
    #                  EACH observed variable in the model
    # - direct effect of 'epsilon': diag(theta)
    # - total effect of 'disturbance': Lambda %*% (I - Beta)^{-1}
    lv_res <- diag(psi)
    lv_res[lv.x.idx] <- 0 # remove 'total/exo' variances from psi
    y_res <- diag(theta)
    t_res <- t(t(lambda_bin %*% ibinv_bin) * lv_res)
    total <- cbind(y_res, t_res)

    # construct piv (potential ivs):
    # initial set of potential instruments for each equation, by selecting
    # observed variables that are unaffected by the disturbances or
    # uniquenesses in that equation
    piv <- matrix(1:nvar, nrow = nvar, ncol = nvar, byrow = TRUE)
    colnames(piv) <- rownames(piv) <- ov.names
    for (i in 1:nvar) {
      res_id <- unique(comp[i, comp[i, ] != 0])
      tmp <- apply(total, 2, function(x) !x %in% res_id)
      rm.idx <- which(!apply(tmp, 1, all))
      piv[i, rm.idx] <- 0
    }
    # diagonal must be zero (not needed, but just to be safe)
    diag(piv) <- 0
    # ov.y -> zero column
    if (length(lavpta$vidx$ov.y[[b]]) > 0L) {
      piv[, lavpta$vidx$ov.y[[b]]] <- 0
    }
    # ov.x -> zero row
    if (length(lavpta$vidx$ov.x[[b]]) > 0L) {
      piv[lavpta$vidx$ov.x[[b]], ] <- 0
    }


    # construct iv
    # Starting from piv, remove any potential ivs that are affected by a
    # disturbance that correlates with any disturbance in the composite
    # disturbance term
    iv <- piv
    for (i in 1:nvar) {
      for (p in 1:nvar) {
        # do check for this iv, if false, set to zero
        comp.idx <- which(y_res %in% comp[i, ])
        total.idx <- which(y_res %in% total[piv[i, p], ])
        cov_values <- theta[as.matrix(expand.grid(comp.idx, total.idx))]
        if (any(cov_values != 0)) {
          iv[i, p] <- 0
        }
      }
    }

    # remove 'empty' rows
    empty.idx <- which(apply(pred, 1L, function(x) all(x == 0)) &
                         rownames(pred) %in% ov.names)
    if (length(empty.idx) > 0L) {
      iv <- iv[-empty.idx, , drop = FALSE]
      comp <- comp[-empty.idx, , drop = FALSE]
      pred <- pred[-empty.idx, , drop = FALSE]
      pred_orig <- pred_orig[-empty.idx, , drop = FALSE]
    }
    # replace colnames of pred by marker names
    idx <- match(names(lv.marker), colnames(pred))
    colnames(pred)[idx] <- lv.marker

    # prepare list
    eqs <- vector("list", length = nrow(iv))
    for (j in seq_along(eqs)) {
      cet <- colnames(comp)[which(comp[j, ] != 0)]
      x.idx <- which(pred[j, ] != 0)
      eqs[[j]] <- list(
        lhs_new = rownames(pred)[j],
        rhs_new = colnames(pred[j, x.idx, drop = FALSE]),
        lhs = rownames(pred_orig)[j],
        rhs = colnames(pred_orig[j, x.idx, drop = FALSE]),
        pt = unname(pred_orig[j, x.idx]),
        cet = cet,
        markers = unique(lv.marker[cet[-1]]),
        miiv = colnames(iv[j, iv[j, ] != 0, drop = FALSE])
      )
    }

    iv_list[[b]] <- eqs
  }

  iv_list
}

# algorithm 2: create a model where error/disturbance terms are latent
#              variables (a bit like the comprehensive RAM representation)
#              and compute 'big_sigma':
#              - to identify ov's that are uncorrelated with the errors in
#                the composite error term
#              - to identify ov's that are correlated with at least one
#                marker
#              the intersection give the (full list of) miivs
#
# the implementation below is loosely based on (but not identical to) the
# 'miivs()' function from the MIIVsem package
#
lav_model_find_iv_miivsem <- function(lavmodel = NULL, lavpta = NULL) {
  # check representation
  if (lavmodel@representation != "LISREL") {
    lav_msg_stop(gettext(
      "this function only works with LISREL representation",
      " (for now)"
    ))
  }

  # create model matrices with 'user/partable' entries
  glist <- lavmodel@GLIST
  for (mm in seq_along(glist)) {
    dimnames(glist[[mm]]) <- lavmodel@dimNames[[mm]]
    m.idx <- lavmodel@m.user.idx[[mm]]
    x.idx <- lavmodel@x.user.idx[[mm]]
    glist[[mm]][, ] <- 0.0
    glist[[mm]][m.idx] <- x.idx
  }
  lambda.idx <- which(names(glist) == "lambda")
  beta.idx <- which(names(glist) == "beta")
  psi.idx <- which(names(glist) == "psi")
  theta.idx <- which(names(glist) == "theta")

  # nblocks
  nblocks <- lavmodel@nblocks

  # repeat for every block
  iv_list <- vector("list", length = nblocks)
  for (b in seq_len(nblocks)) {
    lv.names <- lavpta$vnames$lv.regular[[b]]
    ov.names <- lavpta$vnames$ov[[b]]
    lv.marker <- lavpta$vnames$lv.marker[[b]]

    # model matrices for this block
    lambda <- glist[[lambda.idx[b]]]
    beta <- glist[[beta.idx[b]]]
    theta <- glist[[theta.idx[b]]]
    psi <- glist[[psi.idx[b]]]

    # 'exogenous' variables (in beta/psi) (including dummy ov's)
    if (is.null(beta)) {
      lv.x.idx <- lavpta$vidx$lv.x[[b]]
    } else {
      lv.x.idx <- which(apply(beta, 1L, function(x) all(x == 0)))
    }

    # lambda+beta
    lambda_beta <- pred_orig <- rbind(lambda, beta)
    lambda_beta[lambda_beta != 0] <- as.numeric(NA)
    # add markers
    row.idx <- match(lavpta$vnames$lv.marker[[b]], rownames(lambda_beta))
    col.idx <- match(names(lavpta$vnames$lv.marker[[b]]), colnames(lambda_beta))
    lambda_beta[cbind(row.idx, col.idx)] <- 1
    empty.idx <- which(apply(lambda_beta, 1L, function(x) all(x == 0)))
    if (length(empty.idx) > 0) {
      lambda_beta <- lambda_beta[-empty.idx, , drop = FALSE]
      pred_orig <- pred_orig[-empty.idx, , drop = FALSE]
    }

    # error + zero matrix, to create gamma and beta
    error_matrix <- diag(nrow = nrow(lambda_beta))
    rownames(error_matrix) <- rownames(lambda_beta)
    colnames(error_matrix) <- paste("e.", rownames(lambda_beta), sep = "")
    zero_matrix <- error_matrix * 0
    colnames(zero_matrix) <- rownames(lambda_beta)

    # construct gamma
    gamma <- cbind(
      lambda_beta[, lv.x.idx, drop = FALSE],
      error_matrix
    )
    gamma_orig <- cbind(
      pred_orig[, lv.x.idx, drop = FALSE],
      error_matrix
    )

    # construct beta
    beta <- beta_orig <- zero_matrix
    tmp <- lambda_beta[, -lv.x.idx, drop = FALSE]
    beta[, colnames(tmp)] <- tmp
    tmp <- pred_orig[, -lv.x.idx, drop = FALSE]
    beta_orig[, colnames(tmp)] <- tmp

    # construct Phi
    Phi <- lav_matrix_bdiag(
      psi[lv.x.idx, lv.x.idx, drop = FALSE],
      theta,
      psi[-lv.x.idx, -lv.x.idx, drop = FALSE]
    )
    empty.idx <- which(apply(Phi, 1L, function(x) all(x == 0)))
    if (length(empty.idx) > 0) {
      Phi <- Phi[-empty.idx, -empty.idx, drop = FALSE]
    }
    colnames(Phi) <- rownames(Phi) <- colnames(gamma)
    Phi[Phi != 0] <- as.numeric(NA)

    tmp <- crossprod(gamma)
    tmp[, ] <- 0
    Beta <- lav_matrix_bdiag(beta, tmp)
    I <- diag(nrow(Beta))
    diag(tmp) <- 1
    Gamma <- rbind(gamma, tmp)

    # to compute big_sigma, we can replace NA by 1
    Beta[is.na(Beta)] <- 1
    Gamma[is.na(Gamma)] <- 1
    Phi[is.na(Phi)] <- 1
    ib_inv <- solve(I - Beta)
    big_sigma <- ib_inv %*% Gamma %*% Phi %*% t(Gamma) %*% t(ib_inv)
    sigma_ov <- big_sigma[ov.names, , drop = FALSE] # only rows with ov.names

    # matrix of all regressions
    gamma_beta <- cbind(gamma, beta)
    gamma_beta_orig <- cbind(gamma_orig, beta_orig)
    e_names <- colnames(gamma_beta)[grep("e\\.", colnames(gamma_beta))]

    # select only those rows of gamma_beta that have at least one free (NA)
    # element
    free.idx <- which(apply(is.na(gamma_beta), 1, any))
    eqs_y_free <- unique(rownames(gamma_beta)[free.idx])
    eqs <- vector("list", length = length(eqs_y_free))
    for (j in seq_along(eqs)) {
      # lhs
      eq_lhs <- eq_lhs_orig <- eqs_y_free[j]

      # lhs + rhs
      zero.idx <- which(gamma_beta[eq_lhs, ] != 0) # fixed or error term
      na.idx <- which(is.na(gamma_beta[eq_lhs, ])) # free/NA
      eq_rhs <- eq_rhs_orig <- colnames(gamma_beta)[c(zero.idx, na.idx)]
      eq_all <- c(eq_lhs, eq_rhs)
      pt <- unname(gamma_beta_orig[eq_lhs, na.idx])

      # markers + composite error term (cet)
      cet <- paste("e.", eq_lhs, sep = "")

      # replace all latent variables by their markers
      markers <- character(0L)
      lv_in_eq <- eq_all[eq_all %in% lv.names]
      if (length(lv_in_eq) > 0L) {
        markers <- unname(lv.marker[lv_in_eq])
        # replace in eq_all
        eq_all[match(lv_in_eq, eq_all)] <- markers
        # add marker errors to composite error term
        cet <- c(cet, paste("e.", markers, sep = ""))
      }

      # any markers that are still latent? (eg higher-order factors)
      if (any(eq_all %in% lv.names)) {
        lv_in_eq <- eq_all[eq_all %in% lv.names]
        while (length(lv_in_eq) > 0L) {
          new_markers <- unname(lv.marker[lv_in_eq])
          eq_all[match(lv_in_eq, eq_all)] <- new_markers
          markers[match(lv_in_eq, markers)] <- new_markers
          cet <- c(cet, paste("e.", new_markers, sep = ""))
          lv_in_eq <- eq_all[eq_all %in% lv.names]
        }
      }
      markers <- unique(markers)
      eq_lhs <- eq_all[1]
      eq_rhs <- eq_all[-1]

      # observed variables that are uncorrelated with all the components of
      # the composite error term
      ov_uncorrelated_with_cet <-
        ov.names[apply(sigma_ov[, cet, drop = FALSE] == 0, 1L, all)]

      # observed variables that have at least one non-zero implied correlation
      # with any of the marker indicators (vector i).
      ov_correlated_with_markers <-
        ov.names[apply(sigma_ov[, markers, drop = FALSE] != 0, 1, all)]

      # valid instruments
      miivs <- intersect(ov_uncorrelated_with_cet, ov_correlated_with_markers)

      eqs[[j]] <- list(
        lhs_new = eq_lhs,
        rhs_new = setdiff(eq_rhs, e_names),
        lhs = eq_lhs_orig,
        rhs = setdiff(eq_rhs_orig, e_names),
        pt = pt,
        cet = cet,
        markers = markers,
        miiv = miivs
      )
    }

    iv_list[[b]] <- eqs
  } # nblocks


  iv_list
}
