# find ALL (eligible) *model-implied* instrumental variables (no pruning yet)
#
# three algorithms:
#   1) using 'total effects' of errors/disturbances, based on Bollen & Curran
#      (2004) [algorithm = "bc2004"]
#   2) treating errors/disturbances as latent variables (as in the
#      comprehensive RAM model); similar to miivs() in MIIVsem package
#      [algorithm = "miivsem"]
#   3) using an explicity expression (due to Albert Maydeu-Olivares) for
#      cov(u,y) [algorithm = "covuy"] (the default)
#
# YR 29 Dec 2025 - first version (bc2004 + miivsem)
# YR 24 Jan 2025 - covuy algorithm (including higher-order factors)

lav_model_find_iv <- function(lavobject = NULL, lavmodel = NULL,
                              lavpta = NULL, algorithm = "covuy",
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
  } else {
    iv_list <- lav_model_find_iv_covuy(lavmodel = lavmodel, lavpta = lavpta)
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
  if (lavmodel@meanstructure) {
    nu.idx <- which(names(glist) == "nu")
    alpha.idx <- which(names(glist) == "alpha")
  }

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

    # dummy y?
    if (length(lavmodel@ov.y.dummy.ov.idx[[b]]) > 0L) {
      lv.idx <- c(lv.idx, lavmodel@ov.y.dummy.lv.idx[[b]])
      lv.marker <- c(lv.marker.orig, ov.names[lavmodel@ov.y.dummy.ov.idx[[b]]])
      names(lv.marker) <- c(
        names(lv.marker.orig),
        ov.names[lavmodel@ov.y.dummy.ov.idx[[b]]]
      )
      lv.marker.idx <- c(lv.marker.idx, lavmodel@ov.y.dummy.ov.idx[[b]])
    }

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
    if (lavmodel@meanstructure) {
      nu <- glist[[nu.idx[b]]]
      alpha <- glist[[alpha.idx[b]]]
    }

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
      if (lavmodel@meanstructure) {
        int_orig <- rbind(nu, alpha)
      }
    } else {
      reg_bin <- beta_bin
      pred_orig <- beta
      if (lavmodel@meanstructure) {
        int_orig <- alpha
      }
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
    # dummy y + regular lv's
    if (length(lv.idx) > 0L && !is.null(beta)) {
      tmp_lv <- tmp[length(ov.names) + lv.idx, , drop = FALSE]
      orig_lv <- pred_orig[length(ov.names) + lv.idx, , drop = FALSE]
      if (lavmodel@meanstructure) {
        int_lv <- int_orig[length(ov.names) + lv.idx, , drop = FALSE]
      }

      # tmp[lv.marker.idx, ] <- tmp[length(ov.names) + lv.idx, ]
      # careful, lv.marker.idx may contain the same element multiple times!
      zerolv.idx <- which(apply(tmp_lv, 1, function(x) all(x == 0)))
      if (length(zerolv.idx) > 0L) {
        tmp_lv <- tmp_lv[-zerolv.idx, , drop = FALSE]
        orig_lv <- orig_lv[-zerolv.idx, , drop = FALSE]
        if (lavmodel@meanstructure) {
          int_lv <- int_lv[-zerolv.idx, , drop = FALSE]
        }
      }
      r.idx <- match(lv.marker[rownames(tmp_lv)], rownames(pred))
      pred[r.idx, ] <- tmp_lv
      pred_orig[r.idx, ] <- orig_lv
      rownames(pred_orig)[r.idx] <- rownames(tmp_lv)
      pred_orig <- pred_orig[1:nvar, , drop = FALSE]
      if (lavmodel@meanstructure) {
        int_orig[r.idx, ] <- int_lv
        rownames(int_orig)[r.idx] <- rownames(int_lv)
        int_orig <- int_orig[1:nvar, , drop = FALSE]
      }
    }
    # non-zero entries of 'pred' contain the marker idx of the latent
    # predictors (and perhaps the dummy ov indices)

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
      # check if we need instruments at all
      rhs <- pred[i, pred[i, ] != 0, drop = FALSE]
      if (length(rhs) > 0L && all(comp[i, colnames(rhs)] == 0)) {
        # no instruments needed
        iv[i, -rhs] <- 0
      } else {
        for (p in 1:nvar) {
          # do check for this iv, if false, set to zero
          comp.idx <- which(y_res %in% comp[i, ])
          total.idx <- which(y_res %in% total[piv[i, p], ])
          cov_values <- theta[as.matrix(expand.grid(comp.idx, total.idx))]
          if (any(cov_values != 0)) {
            iv[i, p] <- 0
          }
        }
      } # instruments needed
    }

    # remove 'empty' rows
    empty.idx <- which(apply(pred, 1L, function(x) all(x == 0)) &
      rownames(pred) %in% ov.names)
    if (length(empty.idx) > 0L) {
      iv <- iv[-empty.idx, , drop = FALSE]
      comp <- comp[-empty.idx, , drop = FALSE]
      pred <- pred[-empty.idx, , drop = FALSE]
      pred_orig <- pred_orig[-empty.idx, , drop = FALSE]
      if (lavmodel@meanstructure) {
        int_orig <- int_orig[-empty.idx, , drop = FALSE]
      }
    }
    # replace colnames of pred by marker names
    idx <- match(names(lv.marker), colnames(pred))
    colnames(pred)[idx] <- lv.marker

    # prepare list
    eqs <- vector("list", length = nrow(iv))
    for (j in seq_along(eqs)) {
      cet <- colnames(comp)[which(comp[j, ] != 0)]
      x.idx <- which(pred[j, ] != 0)
      ptint <- integer(0L)
      if (lavmodel@meanstructure) {
        ptint <- unname(int_orig[j, 1])
      }
      eqs[[j]] <- list(
        lhs_new = rownames(pred)[j],
        rhs_new = colnames(pred[j, x.idx, drop = FALSE]),
        lhs = rownames(pred_orig)[j],
        rhs = colnames(pred_orig[j, x.idx, drop = FALSE]),
        pt = unname(pred_orig[j, x.idx]),
        ptint = ptint,
        cet = cet,
        markers = unique(lv.marker[cet[-1]]),
        miiv = colnames(iv[j, iv[j, ] != 0, drop = FALSE])
      )
    }

    # reorder so that ov lhs come first, then lv lhs
    lhs <- sapply(eqs, "[[", "lhs")
    ov.idx <- which(lhs %in% ov.names)
    lv.idx <- seq_along(lhs)[-ov.idx]
    eqs <- eqs[c(ov.idx, lv.idx)]

    iv_list[[b]] <- eqs
  } # blocks

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
  if (lavmodel@meanstructure) {
    nu.idx <- which(names(glist) == "nu")
    alpha.idx <- which(names(glist) == "alpha")
  }

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
    if (lavmodel@meanstructure) {
      nu <- glist[[nu.idx[b]]]
      alpha <- glist[[alpha.idx[b]]]
    }


    # 'exogenous' variables (in beta/psi) (including dummy ov's)
    if (is.null(beta)) {
      lv.x.idx <- lavpta$vidx$lv.x[[b]]
    } else {
      lv.x.idx <- which(apply(beta, 1L, function(x) all(x == 0)))
    }

    # lambda+beta
    lambda_beta <- pred_orig <- rbind(lambda, beta)
    if (lavmodel@meanstructure) {
      int_orig <- rbind(nu, alpha)
    }
    lambda_beta[lambda_beta != 0] <- as.numeric(NA)
    # add markers
    row.idx <- match(lavpta$vnames$lv.marker[[b]], rownames(lambda_beta))
    col.idx <- match(names(lavpta$vnames$lv.marker[[b]]), colnames(lambda_beta))
    lambda_beta[cbind(row.idx, col.idx)] <- 1
    empty.idx <- which(apply(lambda_beta, 1L, function(x) all(x == 0)))
    if (length(empty.idx) > 0) {
      lambda_beta <- lambda_beta[-empty.idx, , drop = FALSE]
      pred_orig <- pred_orig[-empty.idx, , drop = FALSE]
      if (lavmodel@meanstructure) {
        int_orig <- int_orig[-empty.idx, , drop = FALSE]
      }
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
      ptint <- integer(0L)
      if (lavmodel@meanstructure) {
        ptint <- unname(int_orig[eq_lhs, 1L])
      }

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
        ptint = ptint,
        cet = cet,
        markers = markers,
        miiv = miivs
      )
    }

    iv_list[[b]] <- eqs
  } # nblocks


  iv_list
}

# algorithm 3:
# - use random data to create sigma (used to check if predictors are correlated
#   with instruments)
# - use an explicit expression to compute cov(u,y), where u_j is the disturbance
#   term of an equation (used to verify that instruments are not correlated
#   with this disturbance term)
lav_model_find_iv_covuy <- function(lavmodel = NULL, lavpta = NULL) {
  # check representation
  if (lavmodel@representation != "LISREL") {
    lav_msg_stop(gettext(
      "this function only works with LISREL representation",
      " (for now)"
    ))
  }

  # create model matrices with 1) 'user/partable' entries, 2) random entries
  glist <- glistr <- lavmodel@GLIST
  glist_names <- names(glist)
  for (mm in seq_along(glist)) {
    dimnames(glist[[mm]]) <- lavmodel@dimNames[[mm]]
    m.idx <- lavmodel@m.user.idx[[mm]]
    x.idx <- lavmodel@x.user.idx[[mm]]
    # partable entries
    glist[[mm]][, ] <- 0.0
    glist[[mm]][m.idx] <- x.idx
    # random entries (overriding user specified!)
    f.idx <- lavmodel@m.free.idx[[mm]]
    glistr[[mm]][f.idx] <- (runif(length(f.idx), min = 0.1, max = 0.5) *
      sign(runif(length(f.idx), min = -1, max = 1)))
    if (glist_names[mm] %in% c("theta", "psi")) {
      tmp <- glistr[[mm]]
      # make sure diagonal is large enough (but only for nonzero values)
      zerodiag.idx <- which(diag(tmp) == 0)
      diag(tmp) <- diag(tmp) + 1.5
      diag(tmp)[zerodiag.idx] <- 0
      # make symmetric
      glistr[[mm]] <- (tmp + t(tmp)) / 2
    }
  }


  # generate random sigma per block
  sigma <- lav_model_sigma(lavmodel, GLIST = glistr, extra = FALSE)
  sigma_aug <- lav_model_cov_both(lavmodel, GLIST = glistr)

  # number of blocks
  nblocks <- lavpta$nblocks
  lambda.idx <- which(glist_names == "lambda")
  beta.idx <- which(glist_names == "beta")
  if (lavmodel@meanstructure) {
    nu.idx <- which(names(glist) == "nu")
    alpha.idx <- which(names(glist) == "alpha")
  }


  # repeat for every block
  iv_list <- vector("list", length = nblocks)
  for (b in seq_len(nblocks)) {
    # extract information to create Bollen & Bauer (2004) matrices
    ov.names <- lavmodel@dimNames[[lambda.idx[b]]][[1]]
    lv.names <- lavmodel@dimNames[[lambda.idx[b]]][[2]]
    both.names <- c(ov.names, lv.names)
    nvar <- length(ov.names)
    nfac <- length(lv.names)

    # model matrices for this block
    lambda <- glist[[lambda.idx[b]]]
    beta <- glist[[beta.idx[b]]]
    if (lavmodel@meanstructure) {
      nu <- glist[[nu.idx[b]]]
      alpha <- glist[[alpha.idx[b]]]
    }

    # sigma + sigma_aug
    this_sigma <- sigma[[b]]
    rownames(this_sigma) <- colnames(this_sigma) <- ov.names
    this_sigma_aug <- sigma_aug[[b]]
    rownames(this_sigma_aug) <- colnames(this_sigma_aug) <- both.names

    lambdar <- glistr[[lambda.idx[b]]]
    if (is.null(glistr[[beta.idx[b]]])) {
      betar <- diag(1, nrow = nfac)
    } else {
      betar <- glistr[[beta.idx[b]]]
    }

    # all markers (including higher-order)
    lv.marker <- lv.names
    names(lv.marker) <- lv.names
    lv.marker[match(names(lavpta$vnames$lv.marker[[b]]), lv.names)] <-
      lavpta$vnames$lv.marker[[b]]
    lv.marker.idx <- match(lv.marker, both.names)

    # construct cov_u_y
    pi_mat <- rbind(lambdar, betar)
    rownames(pi_mat) <- both.names
    colnames(pi_mat) <- lv.names
    pi_mat[lv.marker.idx, ] <- betar
    pi_zero_mat <- matrix(0, nrow = nrow(pi_mat), ncol = nrow(pi_mat))
    pi_zero_mat[, lv.marker.idx] <- pi_mat
    cov_u_y_aug <- this_sigma_aug - pi_zero_mat %*% this_sigma_aug
    # zap small elements to be exactly zero
    cov_u_y_aug[abs(cov_u_y_aug) < 1e-07] <- 0.0

    if (length(lavpta$vnames$lv.ho[[b]]) > 0L) {
      # 'add' lv marker rows to ov marker rows
      ov.marker <- lv.marker
      lv.names.nox <- lv.names[!lv.names %in% ov.names]
      lv_in_marker <- lv_beta <- ov.marker[ov.marker %in% lv.names.nox]
      while (length(lv_in_marker) > 0L) {
        new_markers <- unname(ov.marker[lv_in_marker])
        ov.marker[match(lv_in_marker, ov.marker)] <- new_markers
        lv_in_marker <- ov.marker[ov.marker %in% lv.names.nox]
      }
      ho.idx <- match(names(lv_beta), both.names)
      target.idx <- match(ov.marker[names(lv_beta)], both.names)
      cov_u_y_aug[target.idx, ] <-
        cov_u_y_aug[target.idx, ] + cov_u_y_aug[ho.idx, ]
    } else {
      ov.marker <- lv.marker
    }
    cov_u_y <- cov_u_y_aug[seq_len(nvar), seq_len(nvar)]

    # construct pred
    pred <- rbind(lambda, beta)
    if (lavmodel@meanstructure) {
      int_orig <- rbind(nu, alpha)
    }
    # remove scaling '1' for markers (if any)
    if (length(lv.marker.idx) > 0L) {
      pred[lv.marker.idx, ] <- 0
    }
    pred_orig <- pred

    # replace col/rownames of pred by marker names
    idx <- match(names(ov.marker), colnames(pred))
    idx <- idx[!is.na(idx)]
    colnames(pred)[idx] <- ov.marker
    idx <- match(names(ov.marker), rownames(pred))
    idx <- idx[!is.na(idx)]
    rownames(pred)[idx] <- ov.marker

    # remove 'empty' rows
    empty.idx <- which(apply(pred, 1L, function(x) all(x == 0)))
    if (length(empty.idx) > 0L) {
      pred <- pred[-empty.idx, , drop = FALSE]
      pred_orig <- pred_orig[-empty.idx, , drop = FALSE]
      if (lavmodel@meanstructure) {
        int_orig <- int_orig[-empty.idx, , drop = FALSE]
      }
    }

    # prepare list
    eqs <- vector("list", length = nrow(pred))
    for (j in seq_along(eqs)) {
      lhs <- rownames(pred)[j]
      x.idx <- which(pred[j, ] != 0)
      rhs <- colnames(pred)[x.idx]
      ptint <- integer(0L)
      if (lavmodel@meanstructure) {
        ptint <- unname(int_orig[j, 1])
      }

      # instruments needed?
      iv_flag <- TRUE
      if (all(cov_u_y[lhs, rhs] == 0)) {
        iv_flag <- FALSE
        miiv <- rhs
      } else {
        # correlated with at least one predictor
        correlated_with_pred <- apply(
          this_sigma[, rhs, drop = FALSE], 1L,
          function(x) any(x != 0)
        )
        # instruments should be uncorrelated with the  term of the eq (u_j)
        uncorrelated_with_u <- cov_u_y[lhs, ] == 0
        miiv <- ov.names[correlated_with_pred & uncorrelated_with_u]
      }

      eqs[[j]] <- list(
        lhs_new = rownames(pred)[j],
        rhs_new = colnames(pred[j, x.idx, drop = FALSE]),
        lhs = rownames(pred_orig)[j],
        rhs = colnames(pred_orig[j, x.idx, drop = FALSE]),
        pt = unname(pred_orig[j, x.idx]),
        ptint = ptint,
        iv_flag = iv_flag,
        miiv = miiv
      )
    }

    iv_list[[b]] <- eqs
  }

  iv_list
}
