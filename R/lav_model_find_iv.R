# find ALL (eligible) *model-implied* instrumental variables (no pruning yet)
#
# three algorithms:
#   1) using 'total effects' of errors/disturbances, based on Bollen & Bauer
#      (2004) [algorithm = "bb2004"]
#   2) treating errors/disturbances as latent variables (as in the
#      comprehensive RAM model); similar to miivs() in MIIVsem package
#      [algorithm = "miivsem"]
#   3) using an explicit expression (due to Albert Maydeu-Olivares) for
#      cov(u,y) [algorithm = "covuy"] (the default)
#
# YR 29 Dec 2025 - first version (bb2004 + miivsem)
# YR 24 Jan 2025 - covuy algorithm (including higher-order factors)

lav_model_find_iv <- function(lavobject = NULL, lavmodel = NULL,
                              lavpartable = NULL, algorithm = "covuy",
                              output = "list", drop_list_single_group = FALSE) {
  # check output
  output <- tolower(output)
  stopifnot(output %in% c("list", "table"))

  # lavobject or components?
  if (!is.null(lavobject)) {
    stopifnot(inherits(lavobject, "lavaan"))
    lavpartable <- lavobject@ParTable
    lavpta <- lavobject@pta
    lavmodel <- lavobject@Model
  } else {
    # get lavpta
    lavpta <- lav_pt_attributes(lavpartable)
  }

  # sanity checks
  nblocks <- lavpta$nblocks
  for (b in seq_len(nblocks)) {
    lv_idx <- lavpta$vidx$lv.regular[[b]]
    lv_marker <- lavpta$vnames$lv.marker[[b]]
    if (length(lv_idx) > 0L) {
      # do have 'clear' marker/scaling indicators?
      empty_idx <- which(nchar(lv_marker) == 0L)
      if (length(empty_idx) > 0L) {
        tmp_string <- paste(names(lv_marker)[empty_idx], collapse = " ")
        lav_msg_stop(gettextf("no clear marker/scaling indicator found for
                               factor(s): %s", tmp_string))
      }
      # marker/scaling indicator cannot be a dependent variable in an ~ equation
      badmarker_idx <- which(lv_marker %in% lavpta$vnames$eqs.y[[b]])
      if (length(badmarker_idx) > 0L) {
        tmp_string <- paste(lv_marker[badmarker_idx], collapse = " ")
        lav_msg_stop(gettextf("marker/scaling indicator cannot be a dependent
                               variable in a regression: %s", tmp_string))
      }
    }
  }

  # first find all model-implied instruments
  algorithm <- tolower(algorithm)
  if (algorithm == "bb2004") {
    iv_list <- lav_model_find_iv_bb2004(lavmodel = lavmodel, lavpta = lavpta)
  } else if (algorithm == "miivsem") {
    iv_list <- lav_model_find_iv_miivsem(lavmodel = lavmodel, lavpta = lavpta)
  } else {
    iv_list <- lav_model_find_iv_covuy(lavmodel = lavmodel, lavpta = lavpta)
  }

  # check for user-specified instruments
  if (any(lavpartable$op == "|~")) {
    for (b in seq_len(nblocks)) {
      iv_list[[b]] <- lapply(iv_list[[b]], function(eq) {
        lhs <- eq$lhs[1] # we assume a single lhs
        iv_idx <- which(lavpartable$op == "|~" & lavpartable$lhs == lhs)
        if (length(iv_idx)) {
          # override iv
          eq$iv <- lavpartable$rhs[iv_idx]
          eq$iv_type <- "user"
        }
        eq
      })
    } # blocks
  }

  if (output == "table") {
    table <- vector("list", length = nblocks)
    for (b in seq_len(nblocks)) {
      eqs <- iv_list[[b]]
      lhs <- sapply(eqs, "[[", "lhs")
      rhs <- sapply(lapply(eqs, "[[", "rhs"), paste, collapse = " + ")
      lhs_new <- sapply(eqs, "[[", "lhs_new")
      rhs_new <- sapply(lapply(eqs, "[[", "rhs_new"), paste, collapse = " + ")
      iv <- sapply(lapply(eqs, "[[", "iv"), paste, collapse = ", ")
      type <- sapply(eqs, "[[", "iv_type")
      table[[b]] <- data.frame(
        lhs = lhs, rhs = rhs,
        lhs_new = lhs_new, rhs_new = rhs_new, type = type, iv = iv
      )
      class(table[[b]]) <- c("lavaan.data.frame", "data.frame")
    }
    out <- table
  } else {
    out <- iv_list
  }

  if (nblocks == 1L && drop_list_single_group) {
    out <- out[[1]]
  }

  out
}


# algorithm 1:
# loosely based on Bollen & Bauer (2004), but with support for higher-order
# factors
lav_model_find_iv_bb2004 <- function(lavmodel = NULL, lavpta = NULL) {
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
    m_idx <- lavmodel@m.user.idx[[mm]]
    x_idx <- lavmodel@x.user.idx[[mm]]
    glist[[mm]][, ] <- 0.0
    glist[[mm]][m_idx] <- x_idx
  }

  # number of blocks
  nblocks <- lavpta$nblocks
  lambda_idx <- which(names(glist) == "lambda")
  beta_idx <- which(names(glist) == "beta")
  psi_idx <- which(names(glist) == "psi")
  theta_idx <- which(names(glist) == "theta")
  if (lavmodel@meanstructure) {
    nu_idx <- which(names(glist) == "nu")
    alpha_idx <- which(names(glist) == "alpha")
  }

  # repeat for every block
  iv_list <- vector("list", length = nblocks)
  for (b in seq_len(nblocks)) {
    # extract information to create Bollen & Bauer (2004) matrices
    ov_names <- lavpta$vnames$ov[[b]]
    nvar <- length(ov_names)
    lv_names <- lavpta$vnames$lv.regular[[b]]
    lv_marker <- lavpta$vnames$lv.marker[[b]]
    lv_marker_orig <- lavpta$vnames$lv.marker[[b]]
    lv_ho <- lavpta$vnames$lv.ho[[b]]

    lv_idx <- lavpta$vidx$lv.regular[[b]]
    lv_x_idx <- lavpta$vidx$lv.x[[b]]
    lv_marker_idx <- lavpta$vidx$lv.marker[[b]]

    # dummy y?
    if (length(lavmodel@ov.y.dummy.ov.idx[[b]]) > 0L) {
      lv_idx <- c(lv_idx, lavmodel@ov.y.dummy.lv.idx[[b]])
      lv_marker <- c(lv_marker_orig, ov_names[lavmodel@ov.y.dummy.ov.idx[[b]]])
      names(lv_marker) <- c(
        names(lv_marker_orig),
        ov_names[lavmodel@ov.y.dummy.ov.idx[[b]]]
      )
      lv_marker_idx <- c(lv_marker_idx, lavmodel@ov.y.dummy.ov.idx[[b]])
    }

    # keep track of higher order factors
    lv_ho_idx <- which(is.na(lv_marker_idx)) # same as lavpta$vidx$lv.ho
    if (length(lv_ho_idx) > 0L) {
      lv_in_marker <- lv_marker[lv_marker %in% lv_names]
      while (length(lv_in_marker) > 0L) {
        new_markers <- unname(lv_marker[lv_in_marker])
        lv_marker[match(lv_in_marker, lv_marker)] <- new_markers
        lv_in_marker <- lv_marker[lv_marker %in% lv_names]
      }
      lv_marker_idx <- match(lv_marker, ov_names)
    }

    # model matrices for this block
    lambda <- glist[[lambda_idx[b]]]
    beta <- glist[[beta_idx[b]]]
    theta <- glist[[theta_idx[b]]]
    psi <- glist[[psi_idx[b]]]
    if (lavmodel@meanstructure) {
      nu <- glist[[nu_idx[b]]]
      alpha <- glist[[alpha_idx[b]]]
    }
    if (is.null(beta)) {
      beta <- matrix(0, ncol(lambda), ncol(lambda))
      colnames(beta) <- rownames(beta) <- colnames(lambda)
    }

    # binary model matrices: nonzero = 1, zero = 0
    if (length(lv_idx) == 0L) {
      lambda_bin <- diag(nrow = ncol(lambda))
    } else {
      lambda_bin <- (lambda != 0) * 1L
    }
    beta_bin <- (beta != 0) * 1L
    ibinv <- solve(diag(nrow(beta_bin)) - beta_bin)
    ibinv_bin <- (ibinv != 0) * 1L

    # construct pred
    if (length(lv_idx) > 0L) {
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
    pred_idx <- numeric(ncol(reg_bin))
    if (length(lv_marker_idx) > 0L) {
      pred_idx[lv_idx] <- lv_marker_idx
    }
    if (length(lavmodel@ov.x.dummy.ov.idx[[b]]) > 0L) {
      pred_idx[lavmodel@ov.x.dummy.lv.idx[[b]]] <-
        lavmodel@ov.x.dummy.ov.idx[[b]]
    }
    if (length(lavmodel@ov.y.dummy.ov.idx[[b]]) > 0L) {
      pred_idx[lavmodel@ov.y.dummy.lv.idx[[b]]] <-
        lavmodel@ov.y.dummy.ov.idx[[b]]
    }
    tmp <- t(t(reg_bin) * pred_idx)
    # remove scaling '1' for markers (if any)
    if (length(lv_marker_idx) > 0L) {
      row_idx <- match(lv_marker_orig, rownames(tmp))
      col_idx <- match(names(lv_marker_orig), colnames(tmp))
      tmp[cbind(row_idx, col_idx)] <- 0L
    }
    # keep only the 'ov' part
    pred <- tmp[seq_along(ov_names), , drop = FALSE]
    # replace marker rows, with beta entries of their corresponding lv's
    # dummy y + regular lv's
    if (length(lv_idx) > 0L && !is.null(beta)) {
      tmp_lv <- tmp_lv2 <- tmp[length(ov_names) + lv_idx, , drop = FALSE]
      orig_lv <- pred_orig[length(ov_names) + lv_idx, , drop = FALSE]
      if (lavmodel@meanstructure) {
        int_lv <- int_orig[length(ov_names) + lv_idx, , drop = FALSE]
        tmp_lv2 <- cbind(int_lv, tmp_lv)
      }

      # tmp[lv.marker.idx, ] <- tmp[length(ov.names) + lv.idx, ]
      # careful, lv.marker.idx may contain the same element multiple times!
      zerolv_idx <- which(apply(tmp_lv2, 1, function(x) all(x == 0)))
      if (length(zerolv_idx) > 0L) {
        tmp_lv <- tmp_lv[-zerolv_idx, , drop = FALSE]
        orig_lv <- orig_lv[-zerolv_idx, , drop = FALSE]
        if (lavmodel@meanstructure) {
          int_lv <- int_lv[-zerolv_idx, , drop = FALSE]
        }
      }
      r_idx <- match(lv_marker[rownames(tmp_lv)], rownames(pred))
      pred[r_idx, ] <- tmp_lv
      pred_orig[r_idx, ] <- orig_lv
      rownames(pred_orig)[r_idx] <- rownames(tmp_lv)
      pred_orig <- pred_orig[1:nvar, , drop = FALSE]
      if (lavmodel@meanstructure) {
        int_orig[r_idx, ] <- int_lv
        # also exogenous ov's!
        lv_x_dummy_idx <- lavmodel@ov.x.dummy.lv.idx[[b]]
        if (length(lv_x_dummy_idx) > 0L) {
          int_orig[lavmodel@ov.x.dummy.ov.idx[[b]], 1] <-
                                                alpha[lv_x_dummy_idx, 1]
        }
        rownames(int_orig)[r_idx] <- rownames(int_lv)
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
    if (length(lv_marker_idx) > 0L) {
      x_res <- pred
      # x_res[pred != 0] <- y_res[as.vector(pred[pred != 0])]
      x_res[pred != 0] <- y_res[pred]
      x_res[cbind(lv_marker_idx, seq_along(lv_idx))] <- diag(psi)[lv_idx]
      if (length(lv_ho_idx) > 0L) {
        # higher-order factors: what are their indicators?
        # add diag(psi) value in this row
        x_res <- t(apply(x_res, 1, function(x) {
          res <- x
          hof <- x[lv_ho_idx]
          hof_nonzero <- hof[hof != 0]
          if (length(hof_nonzero) > 0L) {
            # keep going until we only have first-order factors (or ovs)
            target <- unname(lv_marker_orig[names(hof_nonzero)])
            ho_in_target <- target[target %in% lv_ho]
            while (length(ho_in_target) > 0L) {
              new <- unname(lv_marker_orig[ho_in_target])
              target <- c(target, new)
              ho_in_target <- new[new %in% lv_ho]
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
    # - total effect of 'disturbance': Lambda %*% (ii - Beta)^{-1}
    lv_res <- diag(psi)
    lv_res[lv_x_idx] <- 0 # remove 'total/exo' variances from psi
    y_res <- diag(theta)
    t_res <- t(t(lambda_bin %*% ibinv_bin) * lv_res)
    total <- cbind(y_res, t_res)

    # construct piv (potential ivs):
    # initial set of potential instruments for each equation, by selecting
    # observed variables that are unaffected by the disturbances or
    # uniquenesses in that equation
    piv <- matrix(1:nvar, nrow = nvar, ncol = nvar, byrow = TRUE)
    colnames(piv) <- rownames(piv) <- ov_names
    for (i in 1:nvar) {
      res_id <- unique(comp[i, comp[i, ] != 0])
      tmp <- apply(total, 2, function(x) !x %in% res_id)
      rm_idx <- which(!apply(tmp, 1, all))
      piv[i, rm_idx] <- 0
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
          comp_idx <- which(y_res %in% comp[i, ])
          total_idx <- which(y_res %in% total[piv[i, p], ])
          cov_values <- theta[as.matrix(expand.grid(comp_idx, total_idx))]
          if (any(cov_values != 0)) {
            iv[i, p] <- 0
          }
        }
      } # instruments needed
    }

    # add intercepts here
    if (lavmodel@meanstructure) {
      colnames(int_orig) <- "1"
      pred <- cbind(int_orig, pred)
      pred_orig <- cbind(int_orig, pred_orig)

      idx <- match(names(lv_marker), rownames(pred))
      rownames(pred)[idx] <- lv_marker
    }

    # remove 'empty' rows
    empty_idx <- which(apply(pred, 1L, function(x) all(x == 0)) &
      rownames(pred) %in% ov_names)
    if (length(empty_idx) > 0L) {
      iv <- iv[-empty_idx, , drop = FALSE]
      comp <- comp[-empty_idx, , drop = FALSE]
      pred <- pred[-empty_idx, , drop = FALSE]
      pred_orig <- pred_orig[-empty_idx, , drop = FALSE]
      if (lavmodel@meanstructure) {
        int_orig <- int_orig[-empty_idx, , drop = FALSE]
      }
    }
    # replace row/colnames of pred by marker names
    idx <- match(names(lv_marker), colnames(pred))
    colnames(pred)[idx] <- lv_marker

    # prepare list
    eqs <- vector("list", length = nrow(iv))
    for (j in seq_along(eqs)) {
      cet <- colnames(comp)[which(comp[j, ] != 0)]
      x_idx <- which(pred[j, ] != 0)
      ptint <- integer(0L)
      if (lavmodel@meanstructure) {
        x_idx <- x_idx[-1]
        ptint <- unname(pred_orig[j, 1])
      }
      if (lavmodel@meanstructure && length(x_idx) == 0L) {
        rhs_new <- "1"
        rhs <- "1"
        miiv <- "1"
      } else {
        rhs_new <- colnames(pred[j, x_idx, drop = FALSE])
        rhs <- colnames(pred_orig[j, x_idx, drop = FALSE])
        miiv <- colnames(iv[j, iv[j, ] != 0, drop = FALSE])
      }
      iv_type <- "miiv"
      if (identical(rhs_new, miiv)) {
        iv_type <- "ols"
      }
      eqs[[j]] <- list(
        lhs_new = rownames(pred)[j],
        rhs_new = rhs_new,
        lhs = rownames(pred_orig)[j],
        rhs = rhs,
        pt = unname(pred_orig[j, x_idx]),
        ptint = ptint,
        cet = cet,
        markers = unique(lv_marker[cet[-1]]),
        iv_type = iv_type,
        iv = miiv,
        miiv = miiv
      )
    }

    # reorder so that ov lhs come first, then lv lhs
    lhs <- sapply(eqs, "[[", "lhs")
    ov_idx <- match(ov_names[ov_names %in% lhs], lhs)
    lv_idx <- seq_along(lhs)[-ov_idx]
    eqs <- eqs[c(ov_idx, lv_idx)]

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
    m_idx <- lavmodel@m.user.idx[[mm]]
    x_idx <- lavmodel@x.user.idx[[mm]]
    glist[[mm]][, ] <- 0.0
    glist[[mm]][m_idx] <- x_idx
  }
  lambda_idx <- which(names(glist) == "lambda")
  beta_idx <- which(names(glist) == "beta")
  psi_idx <- which(names(glist) == "psi")
  theta_idx <- which(names(glist) == "theta")
  if (lavmodel@meanstructure) {
    nu_idx <- which(names(glist) == "nu")
    alpha_idx <- which(names(glist) == "alpha")
  }

  # nblocks
  nblocks <- lavmodel@nblocks

  # repeat for every block
  iv_list <- vector("list", length = nblocks)
  for (b in seq_len(nblocks)) {
    lv_names <- lavpta$vnames$lv.regular[[b]]
    ov_names <- lavpta$vnames$ov[[b]]
    lv_marker <- lavpta$vnames$lv.marker[[b]]

    # model matrices for this block
    lambda <- glist[[lambda_idx[b]]]
    beta <- glist[[beta_idx[b]]]
    theta <- glist[[theta_idx[b]]]
    psi <- glist[[psi_idx[b]]]
    if (lavmodel@meanstructure) {
      nu <- glist[[nu_idx[b]]]
      alpha <- glist[[alpha_idx[b]]]
    }


    # 'exogenous' variables (in beta/psi) (including dummy ov's)
    if (is.null(beta)) {
      lv_x_idx <- lavpta$vidx$lv.x[[b]]
    } else {
      lv_x_idx <- which(apply(beta, 1L, function(x) all(x == 0)))
    }

    # lambda+beta
    if (is.null(beta)) {
      beta <- matrix(0, ncol(lambda), ncol(lambda))
    }
    lambda_beta <- rbind(lambda, beta)
    pred_orig <- lambda_beta
    lambda_beta[lambda_beta != 0] <- as.numeric(NA)
    # add markers
    row_idx <- match(lavpta$vnames$lv.marker[[b]], rownames(lambda_beta))
    col_idx <- match(names(lavpta$vnames$lv.marker[[b]]), colnames(lambda_beta))
    lambda_beta[cbind(row_idx, col_idx)] <- 1
    pred <- lambda_beta
    pred_orig_full <- pred_orig
    empty_idx <- which(apply(lambda_beta, 1L, function(x) all(x == 0)))
    if (length(empty_idx) > 0) {
      lambda_beta <- lambda_beta[-empty_idx, , drop = FALSE]
      pred_orig <- pred_orig[-empty_idx, , drop = FALSE]
    }

    # error + zero matrix, to create gamma and beta
    error_matrix <- diag(nrow = nrow(lambda_beta))
    rownames(error_matrix) <- rownames(lambda_beta)
    colnames(error_matrix) <- paste("e.", rownames(lambda_beta), sep = "")
    zero_matrix <- error_matrix * 0
    colnames(zero_matrix) <- rownames(lambda_beta)

    # construct gamma
    gamma <- cbind(
      lambda_beta[, lv_x_idx, drop = FALSE],
      error_matrix
    )
    # gamma_orig <- cbind(
    #   pred_orig[, lv_x_idx, drop = FALSE],
    #   error_matrix
    # )

    # construct beta
    beta <- beta_orig <- zero_matrix
    tmp <- lambda_beta[, -lv_x_idx, drop = FALSE]
    beta[, colnames(tmp)] <- tmp
    tmp <- pred_orig[, -lv_x_idx, drop = FALSE]
    beta_orig[, colnames(tmp)] <- tmp

    # construct Phi
    phi <- lav_mat_bdiag(
      psi[lv_x_idx, lv_x_idx, drop = FALSE],
      theta,
      psi[-lv_x_idx, -lv_x_idx, drop = FALSE]
    )
    empty_idx <- which(apply(phi, 1L, function(x) all(x == 0)))
    if (length(empty_idx) > 0) {
      phi <- phi[-empty_idx, -empty_idx, drop = FALSE]
    }
    colnames(phi) <- rownames(phi) <- colnames(gamma)
    phi[phi != 0] <- as.numeric(NA)

    tmp <- crossprod(gamma)
    tmp[, ] <- 0
    beta_1 <- lav_mat_bdiag(beta, tmp)
    ii <- diag(nrow(beta_1))
    diag(tmp) <- 1
    gamma <- rbind(gamma, tmp)

    # to compute big_sigma, we can replace NA by 1
    beta_1[is.na(beta_1)] <- 1
    gamma[is.na(gamma)] <- 1
    phi[is.na(phi)] <- 1
    ib_inv <- solve(ii - beta_1)
    big_sigma <- ib_inv %*% gamma %*% phi %*% t(gamma) %*% t(ib_inv)
    sigma_ov <- big_sigma[ov_names, , drop = FALSE] # only rows with ov.names
    # e_names <- colnames(gamma)[grep("e\\.", colnames(gamma))]


    # matrix of all regressions (nothing has been removed yet, here)
    gamma_beta <- pred
    gamma_beta_orig <- pred_orig_full

    # add intercepts here
    if (lavmodel@meanstructure) {
      int <- rbind(nu, alpha)
      colnames(int) <- "1"
      int_orig <- int
      int[int != 0] <- as.numeric(NA)
      idx <- match(lv_marker, rownames(int))
      int[idx, 1] <- 1
      gamma_beta <- cbind(int, gamma_beta)
      gamma_beta_orig <- cbind(int_orig, gamma_beta_orig)
    }

    # select only those rows of gamma_beta that have at least one free (NA)
    # element
    free_idx <- which(apply(is.na(gamma_beta), 1, any))
    eqs_y_free <- unique(rownames(gamma_beta)[free_idx])
    eqs <- vector("list", length = length(eqs_y_free))
    for (j in seq_along(eqs)) {
      # lhs
      eq_lhs <- eq_lhs_orig <- eqs_y_free[j]

      # index in gamma_beta
      jj <- free_idx[j]

      # lhs + rhs
      zero_idx <- which(gamma_beta[jj, ] != 0) # fixed or error term
      na_idx <- which(is.na(gamma_beta[jj, ])) # free/NA
      eq_rhs <- eq_rhs_orig <- colnames(gamma_beta)[c(zero_idx, na_idx)]
      pt <- unname(gamma_beta_orig[jj, na_idx])
      ptint <- integer(0L)
      if (lavmodel@meanstructure) {
        ptint <- unname(int_orig[jj, 1L])
        eq_rhs <- eq_rhs[-1]
        eq_rhs_orig <- eq_rhs_orig[-1]
      }
      eq_all <- c(eq_lhs, eq_rhs)

      # markers + composite error term (cet)
      cet <- paste("e.", eq_lhs, sep = "")

      # replace all latent variables by their markers
      markers <- character(0L)
      lv_in_eq <- eq_all[eq_all %in% lv_names]
      if (length(lv_in_eq) > 0L) {
        markers <- unname(lv_marker[lv_in_eq])
        # replace in eq_all
        eq_all[match(lv_in_eq, eq_all)] <- markers
        # add marker errors to composite error term
        cet <- c(cet, paste("e.", markers, sep = ""))
      }

      # any markers that are still latent? (eg higher-order factors)
      if (any(eq_all %in% lv_names)) {
        lv_in_eq <- eq_all[eq_all %in% lv_names]
        while (length(lv_in_eq) > 0L) {
          new_markers <- unname(lv_marker[lv_in_eq])
          eq_all[match(lv_in_eq, eq_all)] <- new_markers
          markers[match(lv_in_eq, markers)] <- new_markers
          cet <- c(cet, paste("e.", new_markers, sep = ""))
          lv_in_eq <- eq_all[eq_all %in% lv_names]
        }
      }
      markers <- unique(markers)
      eq_lhs <- eq_all[1]
      eq_rhs <- eq_all[-1]

      if (length(eq_rhs) > 0L) {
        # observed variables that are uncorrelated with all the components of
        # the composite error term
        ov_uncorrelated_with_cet <-
          ov_names[apply(sigma_ov[, cet, drop = FALSE] == 0, 1L, all)]

        # observed variables that have at least one non-zero implied correlation
        # with any of the marker indicators (vector i).
        ov_correlated_with_markers <-
          ov_names[apply(sigma_ov[, markers, drop = FALSE] != 0, 1, all)]

        # valid instruments
        miivs <- intersect(ov_uncorrelated_with_cet, ov_correlated_with_markers)
      } else {
        eq_rhs <- "1"
        eq_rhs_orig <- "1"
        miivs <- "1"
      }

      iv_type <- "miiv"
      if (identical(eq_rhs, miivs)) {
        iv_type <- "ols"
      }

      eqs[[j]] <- list(
        lhs_new = eq_lhs,
        rhs_new = eq_rhs,
        lhs = eq_lhs_orig,
        rhs = eq_rhs_orig,
        pt = pt,
        ptint = ptint,
        cet = cet,
        markers = markers,
        iv_type = iv_type,
        iv = miivs,
        miiv = miivs
      )
    }

    # reorder so that ov lhs come first, then lv lhs
    lhs <- sapply(eqs, "[[", "lhs")
    ov_idx <- match(ov_names[ov_names %in% lhs], lhs)
    lv_idx <- seq_along(lhs)[-ov_idx]
    eqs <- eqs[c(ov_idx, lv_idx)]

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
    m_idx <- lavmodel@m.user.idx[[mm]]
    x_idx <- lavmodel@x.user.idx[[mm]]
    # partable entries
    glist[[mm]][, ] <- 0.0
    glist[[mm]][m_idx] <- x_idx
    # random entries (overriding user specified!)
    f_idx <- lavmodel@m.free.idx[[mm]]
    glistr[[mm]][f_idx] <- (runif(length(f_idx), min = 0.1, max = 0.5) *
      sign(runif(length(f_idx), min = -1, max = 1)))
    if (glist_names[mm] %in% c("theta", "psi")) {
      tmp <- glistr[[mm]]
      # make sure diagonal is large enough (but only for nonzero values)
      zerodiag_idx <- which(diag(tmp) == 0)
      diag(tmp) <- diag(tmp) + 1.5
      diag(tmp)[zerodiag_idx] <- 0
      # make symmetric
      glistr[[mm]] <- (tmp + t(tmp)) / 2
    }
  }

  # generate random sigma per block
  sigma <- lav_model_sigma(lavmodel, glist = glistr, extra = FALSE)
  sigma_aug <- lav_model_cov_both(lavmodel, glist = glistr)

  # number of blocks
  nblocks <- lavpta$nblocks
  lambda_idx <- which(glist_names == "lambda")
  beta_idx <- which(glist_names == "beta")
  if (lavmodel@meanstructure) {
    nu_idx <- which(names(glist) == "nu")
    alpha_idx <- which(names(glist) == "alpha")
  }


  # repeat for every block
  iv_list <- vector("list", length = nblocks)
  for (b in seq_len(nblocks)) {
    # extract information to create Bollen & Bauer (2004) matrices
    ov_names <- lavmodel@dimNames[[lambda_idx[b]]][[1]]
    lv_names <- lavmodel@dimNames[[lambda_idx[b]]][[2]]
    both_names <- c(ov_names, lv_names)
    nvar <- length(ov_names)
    nfac <- length(lv_names)

    # model matrices for this block
    lambda <- glist[[lambda_idx[b]]]
    beta <- glist[[beta_idx[b]]]
    if (is.null(beta)) {
      beta <- matrix(0, nfac, nfac)
    }
    if (lavmodel@meanstructure) {
      nu <- glist[[nu_idx[b]]]
      alpha <- glist[[alpha_idx[b]]]
    }

    # sigma + sigma_aug
    this_sigma <- sigma[[b]]
    rownames(this_sigma) <- colnames(this_sigma) <- ov_names
    this_sigma_aug <- sigma_aug[[b]]
    rownames(this_sigma_aug) <- colnames(this_sigma_aug) <- both_names

    lambdar <- glistr[[lambda_idx[b]]]
    if (is.null(glistr[[beta_idx[b]]])) {
      betar <- diag(1, nrow = nfac)
    } else {
      betar <- glistr[[beta_idx[b]]]
    }

    # all markers (including higher-order)
    lv_marker <- lv_names
    names(lv_marker) <- lv_names
    lv_marker[match(names(lavpta$vnames$lv.marker[[b]]), lv_names)] <-
      lavpta$vnames$lv.marker[[b]]
    lv_marker_idx <- match(lv_marker, both_names)

    # construct cov_u_y
    pi_mat <- rbind(lambdar, betar)
    rownames(pi_mat) <- both_names
    colnames(pi_mat) <- lv_names
    pi_mat[lv_marker_idx, ] <- betar
    pi_zero_mat <- matrix(0, nrow = nrow(pi_mat), ncol = nrow(pi_mat))
    pi_zero_mat[, lv_marker_idx] <- pi_mat
    cov_u_y_aug <- this_sigma_aug - pi_zero_mat %*% this_sigma_aug
    # zap small elements to be exactly zero
    cov_u_y_aug[abs(cov_u_y_aug) < 1e-07] <- 0.0

    if (length(lavpta$vnames$lv.ho[[b]]) > 0L) {
      # 'add' lv marker rows to ov marker rows
      ov_marker <- lv_marker
      lv_names_nox <- lv_names[!lv_names %in% ov_names]
      lv_in_marker <- lv_beta <- ov_marker[ov_marker %in% lv_names_nox]
      while (length(lv_in_marker) > 0L) {
        new_markers <- unname(ov_marker[lv_in_marker])
        ov_marker[match(lv_in_marker, ov_marker)] <- new_markers
        lv_in_marker <- ov_marker[ov_marker %in% lv_names_nox]
      }
      ho_idx <- match(names(lv_beta), both_names)
      target_idx <- match(ov_marker[names(lv_beta)], both_names)
      cov_u_y_aug[target_idx, ] <-
        cov_u_y_aug[target_idx, ] + cov_u_y_aug[ho_idx, ]
    } else {
      ov_marker <- lv_marker
    }
    cov_u_y <- cov_u_y_aug[seq_len(nvar), seq_len(nvar)]

    # construct pred
    pred <- rbind(lambda, beta)

    # if meanstructure, add intercepts
    if (lavmodel@meanstructure) {
      int <- rbind(nu, alpha)
      colnames(int) <- "1" # shorter than "intercept"
      pred <- cbind(int, pred)
    }

    pred_orig <- pred
    # replace col/rownames of pred by marker names
    idx <- match(names(ov_marker), colnames(pred))
    idx <- idx[!is.na(idx)]
    colnames(pred)[idx] <- ov_marker
    idx <- match(names(ov_marker), rownames(pred))
    idx <- idx[!is.na(idx)]
    rownames(pred)[idx] <- ov_marker

    # remove scaling '1' for markers (if any)
    if (length(lv_marker_idx) > 0L) {
      pred[lv_marker_idx, ] <- 0
    }

    # remove 'empty' rows
    empty_idx <- which(apply(pred, 1L, function(x) all(x == 0)))
    if (length(empty_idx) > 0L) {
      pred <- pred[-empty_idx, , drop = FALSE]
      pred_orig <- pred_orig[-empty_idx, , drop = FALSE]
    }

    # prepare list
    eqs <- vector("list", length = nrow(pred))
    for (j in seq_along(eqs)) {
      lhs <- rownames(pred)[j]
      x_idx <- which(pred[j, ] != 0)
      rhs <- colnames(pred)[x_idx]
      ptint <- integer(0L)
      if (lavmodel@meanstructure && pred[j, 1] != 0) {
        ptint <- unname(pred_orig[j, 1])
        rhs <- rhs[-1]
        x_idx <- x_idx[-1]
      }

      # instruments needed?
      iv_flag <- TRUE
      if (lavmodel@meanstructure && length(x_idx) == 0L) {
        iv_flag <- FALSE
        miiv <- "1"
      } else if (all(cov_u_y[lhs, rhs] == 0)) {
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
        miiv <- ov_names[correlated_with_pred & uncorrelated_with_u]
      }

      # prepare rhs, rhs_new, pt
      if (length(x_idx) == 0) {
        rhs_new <- "1"
        rhs <- "1"
        pt <- integer(0L)
      } else {
        rhs_new <- colnames(pred[j, x_idx, drop = FALSE])
        rhs <- colnames(pred_orig[j, x_idx, drop = FALSE])
        pt <- as.integer(unname(pred_orig[j, x_idx]))
      }

      iv_type <- "miiv"
      if (identical(rhs_new, miiv)) {
        iv_type <- "ols"
      }

      eqs[[j]] <- list(
        lhs_new = rownames(pred)[j],
        rhs_new = rhs_new,
        lhs = rownames(pred_orig)[j],
        rhs = rhs,
        pt = pt,
        ptint = ptint,
        iv_flag = iv_flag,
        iv_type = iv_type,
        iv = miiv,
        miiv = miiv
      )
    }

    # reorder so that ov lhs come first, then lv lhs
    lhs <- sapply(eqs, "[[", "lhs")
    ov_idx <- match(ov_names[ov_names %in% lhs], lhs)
    lv_idx <- seq_along(lhs)[-ov_idx]
    eqs <- eqs[c(ov_idx, lv_idx)]

    iv_list[[b]] <- eqs
  }

  iv_list
}
