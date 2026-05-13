# and matrix-representation specific functions:
# - lav_model_sigma
# - lav_model_mu
# - derivative.F

# initial version: YR 2011-01-21: LISREL stuff
# updates:          YR 2011-12-01: group specific extraction
#                   YR 2012-05-17: thresholds
#                   YR 2021-10-04: rename representation.LISREL -> lav_lisrel

lav_lisrel <- function(lavpartable = NULL,
                       target = NULL,
                       extra = FALSE,
                       allow_composites = TRUE,
                       remove_nonexisting = TRUE) {
  # prepare target list
  if (is.null(target)) target <- lavpartable

  stopifnot(!is.null(target$block))

  # prepare output
  n <- length(target$lhs)
  tmp_mat <- character(n)
  tmp_row <- integer(n)
  tmp_col <- integer(n)

  # global settings
  meanstructure <- any(lavpartable$op == "~1")
  # categorical <- any(lavpartable$op == "|")
  composites <- any(lavpartable$op == "<~") && allow_composites
  # group_w_free <- any(lavpartable$lhs == "group" & lavpartable$op == "%")

  # gamma? only if conditional.x
  if (any(lavpartable$op %in% c("~", "<~") & lavpartable$exo == 1L) &&
    !composites) {
    gamma <- TRUE
  } else {
    gamma <- FALSE
  }

  # number of blocks
  nblocks <- lav_partable_nblocks(lavpartable)

  # multilevel?
  nlevels <- lav_partable_nlevels(lavpartable)
  ngroups <- lav_partable_ngroups(lavpartable)

  ov_dummy_names_nox <- vector("list", nblocks)
  ov_dummy_names_x <- vector("list", nblocks)
  if (extra) {
    rep_mm_names <- vector("list", nblocks)
    rep_mm_number <- vector("list", nblocks)
    rep_mm_rows <- vector("list", nblocks)
    rep_mm_cols <- vector("list", nblocks)
    rep_mm_dim_names <- vector("list", nblocks)
    rep_mm_symmetric <- vector("list", nblocks)
  }

  for (g in 1:nblocks) {
    # info from user model per block
    if (gamma) {
      ov_names <- lav_partable_vnames(lavpartable, "ov.nox", block = g)
    } else {
      ov_names <- lav_partable_vnames(lavpartable, "ov", block = g)
    }
    nvar <- length(ov_names)
    lv_names <- lav_partable_vnames(lavpartable, "lv", block = g)
    nfac <- length(lv_names)
    ov_th <- lav_partable_vnames(lavpartable, "th", block = g)
    nth <- length(ov_th)
    ov_names_x <- lav_partable_vnames(lavpartable, "ov.x", block = g)
    nexo <- length(ov_names_x)
    ov_names_nox <- lav_partable_vnames(lavpartable, "ov.nox", block = g)

    # in this representation, we need to create 'phantom/dummy' latent
    # variables for all `x' and `y' variables not in lv.names
    # (only y if conditional.x = TRUE)

    # regression dummys
    if (gamma) {
      tmp_names <-
        unique(lavpartable$lhs[(lavpartable$op == "~" |
          lavpartable$op == "<~") &
          lavpartable$block == g])
      # new in 0.6-12: fix for multilevel + conditional.x: split ov.x
      # are removed from ov.x
      if (nlevels > 1L) {
        if (ngroups == 1L) {
          other_block_names <- lav_partable_vnames(lavpartable, "ov",
            block = seq_len(nblocks)[-g]
          )
        } else {
          # TEST ME
          this_group <- ceiling(g / nlevels)
          blocks_within_group <- (this_group - 1L) * nlevels + seq_len(nlevels)
          other_block_names <- lav_partable_vnames(lavpartable,
            "ov",
            block = blocks_within_group[-g]
          )
        }
        if (length(ov_names_x) > 0L) {
          idx <- which(ov_names_x %in% other_block_names)
          if (length(idx) > 0L) {
            tmp_names <- unique(c(tmp_names, ov_names_x[idx]))
            ov_names_nox <- unique(c(ov_names_nox, ov_names_x[idx]))
            ov_names_x <- ov_names_x[-idx]
            nexo <- length(ov_names_x)
            ov_names <- ov_names_nox
            nvar <- length(ov_names)
          }
        }
      }
    } else {
      if (composites) {
        tmp_names <-
          unique(c(
            lavpartable$lhs[(lavpartable$op == "~") &
              lavpartable$block == g],
            lavpartable$rhs[(lavpartable$op == "~") &
              lavpartable$block == g]
          ))
      } else {
        # old behavior < 0.6-20
        tmp_names <-
          unique(c(
            lavpartable$lhs[(lavpartable$op == "~" |
              lavpartable$op == "<~") &
              lavpartable$block == g],
            lavpartable$rhs[(lavpartable$op == "~" |
              lavpartable$op == "<~") &
              lavpartable$block == g]
          ))
      }
    }
    dummy_names1 <- tmp_names[!tmp_names %in% lv_names]
    # covariances involving dummys
    dummy_cov_idx <- which(lavpartable$op == "~~" & lavpartable$block == g &
      (lavpartable$lhs %in% dummy_names1 |
        lavpartable$rhs %in% dummy_names1))
    # new in 0.5-21: also include covariances involving these covariances...
    dummy_cov_idx1 <- which(lavpartable$op == "~~" & lavpartable$block == g &
      (lavpartable$lhs %in% lavpartable$lhs[dummy_cov_idx] |
        lavpartable$rhs %in% lavpartable$rhs[dummy_cov_idx]))
    dummy_cov_idx <- unique(c(dummy_cov_idx, dummy_cov_idx1))

    dummy_names2 <- unique(c(
      lavpartable$lhs[dummy_cov_idx],
      lavpartable$rhs[dummy_cov_idx]
    ))


    # new in 0.6-7: ~~ between latent and observed
    dummy_cov_ov_lv_idx1 <- which(lavpartable$op == "~~" &
      lavpartable$block == g &
      lavpartable$lhs %in% ov_names &
      lavpartable$rhs %in% lv_names)
    dummy_cov_ov_lv_idx2 <- which(lavpartable$op == "~~" &
      lavpartable$block == g &
      lavpartable$lhs %in% lv_names &
      lavpartable$rhs %in% ov_names)
    dummy_names3 <- unique(c(
      lavpartable$lhs[dummy_cov_ov_lv_idx1],
      lavpartable$rhs[dummy_cov_ov_lv_idx2]
    ))

    # new in 0.6-10: ~~ between observed and observed, but not in ~
    # dummy_orphan_idx <- which(lavpartable$op == "~~" &
    #   lavpartable$block == g &
    #   lavpartable$lhs %in% ov_names &
    #   lavpartable$rhs %in% ov_names &
    #   (!lavpartable$lhs %in% c(
    #     dummy_names1,
    #     dummy_names2
    #   ) |
    #     !lavpartable$rhs %in% c(
    #       dummy_names1,
    #       dummy_names2
    #     )))

    # collect all dummy variables
    dummy_names <- unique(c(dummy_names1, dummy_names2, dummy_names3))


    if (length(dummy_names)) {
      # make sure order is the same as ov.names
      ov_dummy_names_nox[[g]] <-
        ov_names_nox[ov_names_nox %in% dummy_names]
      ov_dummy_names_x[[g]] <-
        ov_names_x[ov_names_x %in% dummy_names]

      # combine them, make sure order is identical to ov.names
      tmp <- ov_names[ov_names %in% dummy_names]

      # same for ov.names.x (if they are not in ov.names) (conditional.x)
      if (length(ov_names_x) > 0L) {
        tmp_x <- ov_names_x[ov_names_x %in% dummy_names]
        tmp <- unique(c(tmp, tmp_x))
      }

      # extend lv.names
      lv_names <- c(lv_names, tmp)
      nfac <- length(lv_names)

      # add 'dummy' =~ entries
      # dummy_mat <- rep("lambda", length(dummy_names))
    } else {
      ov_dummy_names_nox[[g]] <- character(0)
      ov_dummy_names_x[[g]] <- character(0)
    }

    # 1a. "=~" regular indicators
    idx <- which(target$block == g &
      target$op == "=~" & !(target$rhs %in% lv_names))
    tmp_mat[idx] <- "lambda"
    tmp_row[idx] <- match(target$rhs[idx], ov_names)
    tmp_col[idx] <- match(target$lhs[idx], lv_names)

    # 1b. "=~" regular higher-order lv indicators
    idx <- which(target$block == g &
      target$op == "=~" & !(target$rhs %in% ov_names))
    tmp_mat[idx] <- "beta"
    tmp_row[idx] <- match(target$rhs[idx], lv_names)
    tmp_col[idx] <- match(target$lhs[idx], lv_names)

    # 1c. "=~" indicators that are both in ov and lv
    idx <- which(target$block == g &
      target$op == "=~" & target$rhs %in% ov_names &
      target$rhs %in% lv_names)
    tmp_mat[idx] <- "beta"
    tmp_row[idx] <- match(target$rhs[idx], lv_names)
    tmp_col[idx] <- match(target$lhs[idx], lv_names)

    # 1d. "<~" indicators
    if (composites) {
      idx <- which(target$block == g &
        target$op == "<~" & !(target$rhs %in% lv_names))
      tmp_mat[idx] <- "wmat"
      tmp_row[idx] <- match(target$rhs[idx], ov_names)
      tmp_col[idx] <- match(target$lhs[idx], lv_names)
    }

    # 2. "~" regressions
    if (gamma) {
      # gamma
      if (composites) {
        idx <- which(target$rhs %in% ov_names_x &
          target$block == g & target$op == "~")
      } else {
        idx <- which(target$rhs %in% ov_names_x &
          target$block == g & (target$op == "~" |
          target$op == "<~"))
      }
      tmp_mat[idx] <- "gamma"
      tmp_row[idx] <- match(target$lhs[idx], lv_names)
      tmp_col[idx] <- match(target$rhs[idx], ov_names_x)

      # beta
      if (composites) {
        idx <- which(!target$rhs %in% ov_names_x &
          target$block == g & target$op == "~")
      } else {
        idx <- which(!target$rhs %in% ov_names_x &
          target$block == g & (target$op == "~" |
          target$op == "<~"))
      }
      tmp_mat[idx] <- "beta"
      tmp_row[idx] <- match(target$lhs[idx], lv_names)
      tmp_col[idx] <- match(target$rhs[idx], lv_names)
    } else {
      if (composites) {
        idx <- which(target$block == g & target$op == "~")
      } else {
        idx <- which(target$block == g & (target$op == "~" |
          target$op == "<~"))
      }
      tmp_mat[idx] <- "beta"
      tmp_row[idx] <- match(target$lhs[idx], lv_names)
      tmp_col[idx] <- match(target$rhs[idx], lv_names)
    }

    # 3a. "~~" ov
    idx <- which(target$block == g &
      target$op == "~~" & !(target$lhs %in% lv_names))
    tmp_mat[idx] <- "theta"
    tmp_row[idx] <- match(target$lhs[idx], ov_names)
    tmp_col[idx] <- match(target$rhs[idx], ov_names)

    # 3aa. "~~" ov.x
    if (gamma) {
      idx <- which(target$block == g &
        target$op == "~~" & (target$lhs %in% ov_names_x))
      tmp_mat[idx] <- "cov.x"
      tmp_row[idx] <- match(target$lhs[idx], ov_names_x)
      tmp_col[idx] <- match(target$rhs[idx], ov_names_x)
    }

    # 3b. "~~" lv
    idx <- which(target$block == g &
      target$op == "~~" & target$rhs %in% lv_names)
    tmp_mat[idx] <- "psi"
    tmp_row[idx] <- match(target$lhs[idx], lv_names)
    tmp_col[idx] <- match(target$rhs[idx], lv_names)

    # 4a. "~1" ov
    idx <- which(target$block == g &
      target$op == "~1" & !(target$lhs %in% lv_names))
    tmp_mat[idx] <- "nu"
    tmp_row[idx] <- match(target$lhs[idx], ov_names)
    tmp_col[idx] <- 1L

    # 4aa, "~1" ov.x
    if (gamma) {
      idx <- which(target$block == g &
        target$op == "~1" & (target$lhs %in% ov_names_x))
      tmp_mat[idx] <- "mean.x"
      tmp_row[idx] <- match(target$lhs[idx], ov_names_x)
      tmp_col[idx] <- 1L
    }

    # 4b. "~1" lv
    idx <- which(target$block == g &
      target$op == "~1" & target$lhs %in% lv_names)
    tmp_mat[idx] <- "alpha"
    tmp_row[idx] <- match(target$lhs[idx], lv_names)
    tmp_col[idx] <- 1L

    # 5. "|" th
    label <- paste(target$lhs, target$op, target$rhs, sep = "")
    idx <- which(target$block == g &
      target$op == "|" & label %in% ov_th)
    th <- paste(target$lhs[idx], "|", target$rhs[idx], sep = "")
    tmp_mat[idx] <- "tau"
    tmp_row[idx] <- match(th, ov_th)
    tmp_col[idx] <- 1L

    # 6. "~*~" scales
    idx <- which(target$block == g &
      target$op == "~*~")
    tmp_mat[idx] <- "delta"
    tmp_row[idx] <- match(target$lhs[idx], ov_names)
    tmp_col[idx] <- 1L

    # new 0.5-12: catch lower-elements in theta/psi
    idx_lower <- which(tmp_mat %in% c("theta", "psi") & tmp_row > tmp_col)
    if (length(idx_lower) > 0L) {
      tmp <- tmp_row[idx_lower]
      tmp_row[idx_lower] <- tmp_col[idx_lower]
      tmp_col[idx_lower] <- tmp
    }

    # new 0.5-16: group weights
    idx <- which(target$block == g & target$lhs == "group" &
      target$op == "%")
    tmp_mat[idx] <- "gw"
    tmp_row[idx] <- 1L
    tmp_col[idx] <- 1L

    # new in 0.6-22: instruments
    idx <- which(target$block == g & target$lhs == "group" &
      target$op == "|~")
    tmp_mat[idx] <- "miiv"
    tmp_row[idx] <- 0L
    tmp_col[idx] <- 0L

    if (extra) {
      # mRows
      mm_rows <- list(
        tau = nth,
        delta = nvar,
        nu = nvar,
        lambda = nvar,
        wmat = nvar,
        theta = nvar,
        alpha = nfac,
        beta = nfac,
        gamma = nfac,
        cov.x = nexo,
        mean.x = nexo,
        gw = 1L,
        psi = nfac
      )

      # mCols
      mm_cols <- list(
        tau = 1L,
        delta = 1L,
        nu = 1L,
        lambda = nfac,
        wmat = nfac,
        theta = nvar,
        alpha = 1L,
        beta = nfac,
        gamma = nexo,
        cov.x = nexo,
        mean.x = 1L,
        gw = 1L,
        psi = nfac
      )

      # dimNames for LISREL model matrices
      mm_dim_names <- list(
        tau = list(ov_th, "threshold"),
        delta = list(ov_names, "scales"),
        nu = list(ov_names, "intercept"),
        lambda = list(ov_names, lv_names),
        wmat = list(ov_names, lv_names),
        theta = list(ov_names, ov_names),
        alpha = list(lv_names, "intercept"),
        beta = list(lv_names, lv_names),
        gamma = list(lv_names, ov_names_x),
        cov.x = list(ov_names_x, ov_names_x),
        mean.x = list(ov_names_x, "intercepts"),
        gw = list("group", "weight"),
        psi = list(lv_names, lv_names)
      )

      # isSymmetric
      mm_symmetric <- list(
        tau = FALSE,
        delta = FALSE,
        nu = FALSE,
        lambda = FALSE,
        wmat = FALSE,
        theta = TRUE,
        alpha = FALSE,
        beta = FALSE,
        gamma = FALSE,
        cov.x = TRUE,
        mean.x = FALSE,
        gw = FALSE,
        psi = TRUE
      )

      # which mm's do we need? (always include lambda, theta and psi)
      # new: 0.6 this block only!!
      idx_1 <- which(target$block == g)
      if ("wmat" %in% tmp_mat[idx_1]) {
        mm_names <- c("lambda", "wmat", "theta", "psi")
      } else {
        mm_names <- c("lambda", "theta", "psi")
      }

      if ("beta" %in% tmp_mat[idx_1]) {
        mm_names <- c(mm_names, "beta")
      }
      if (meanstructure) {
        mm_names <- c(mm_names, "nu", "alpha")
      }
      if ("tau" %in% tmp_mat[idx_1]) {
        mm_names <- c(mm_names, "tau")
      }
      if ("delta" %in% tmp_mat[idx_1]) {
        mm_names <- c(mm_names, "delta")
      }
      if ("gamma" %in% tmp_mat[idx_1]) {
        mm_names <- c(mm_names, "gamma")
      }
      if ("gw" %in% tmp_mat[idx_1]) {
        mm_names <- c(mm_names, "gw")
      }
      if ("cov.x" %in% tmp_mat[idx_1]) {
        mm_names <- c(mm_names, "cov.x")
      }
      if ("mean.x" %in% tmp_mat[idx_1]) {
        mm_names <- c(mm_names, "mean.x")
      }

      rep_mm_names[[g]] <- mm_names
      rep_mm_number[[g]] <- length(mm_names)
      rep_mm_rows[[g]] <- unlist(mm_rows[mm_names])
      rep_mm_cols[[g]] <- unlist(mm_cols[mm_names])
      rep_mm_dim_names[[g]] <- mm_dim_names[mm_names]
      rep_mm_symmetric[[g]] <- unlist(mm_symmetric[mm_names])
    } # extra
  } # nblocks

  rep_1 <- list(
    mat = tmp_mat,
    row = tmp_row,
    col = tmp_col
  )

  # remove non-existing (NAs)?
  # here we remove `non-existing' parameters; this depends on the matrix
  # representation (eg in LISREL rep, there is no ~~ between lv and ov)
  # if(remove.nonexisting) {
  #    idx <- which( nchar(REP$mat) > 0L &
  #                 !is.na(REP$row) & REP$row > 0L &
  #                 !is.na(REP$col) & REP$col > 0L )
  #   # but keep ==, :=, etc.
  #   idx <- c(idx, which(lavpartable$op %in% c("==", ":=", "<", ">")))
  #   REP$mat <- REP$mat[idx]
  #   REP$row <- REP$row[idx]
  #   REP$col <- REP$col[idx]
  #

  # always add 'ov.dummy.*.names' attributes
  attr(rep_1, "ov.dummy.names.nox") <- ov_dummy_names_nox
  attr(rep_1, "ov.dummy.names.x") <- ov_dummy_names_x

  if (extra) {
    attr(rep_1, "mmNames") <- rep_mm_names
    attr(rep_1, "mmNumber") <- rep_mm_number
    attr(rep_1, "mmRows") <- rep_mm_rows
    attr(rep_1, "mmCols") <- rep_mm_cols
    attr(rep_1, "mmDimNames") <- rep_mm_dim_names
    attr(rep_1, "mmSymmetric") <- rep_mm_symmetric
  }

  rep_1
}


# ETA:
# 1) EETA
# 2) EETAx
# 3) VETA
# 4) VETAx

# 1) EETA
# compute E(ETA): expected value of latent variables (marginal over x)
# - if no eXo (and GAMMA):
#     E(ETA) = (I-B)^-1 ALPHA
# - if eXo and GAMMA:
#     E(ETA) = (I-B)^-1 ALPHA + (I-B)^-1 GAMMA mean.x
lav_lisrel_eeta <- function(mlist = NULL, mean_x = NULL,
                            sample_mean = NULL,
                            ov_y_dummy_ov_idx = NULL,
                            ov_x_dummy_ov_idx = NULL,
                            ov_y_dummy_lv_idx = NULL,
                            ov_x_dummy_lv_idx = NULL) {
  mm_beta <- mlist$beta
  mm_gamma <- mlist$gamma

  # ALPHA? (reconstruct, but no 'fix')
  mm_alpha <- lav_lisrel_alpha0(
    mlist = mlist, sample_mean = sample_mean,
    ov_y_dummy_ov_idx = ov_y_dummy_ov_idx,
    ov_x_dummy_ov_idx = ov_x_dummy_ov_idx,
    ov_y_dummy_lv_idx = ov_y_dummy_lv_idx,
    ov_x_dummy_lv_idx = ov_x_dummy_lv_idx
  )

  # BETA?
  if (!is.null(mm_beta)) {
    ib_inv <- lav_lisrel_ibinv(mlist = mlist)
    # GAMMA?
    if (!is.null(mm_gamma)) {
      eeta <- as.vector(ib_inv %*% mm_alpha + ib_inv %*% mm_gamma %*% mean_x)
    } else {
      eeta <- as.vector(ib_inv %*% mm_alpha)
    }
  } else {
    # GAMMA?
    if (!is.null(mm_gamma)) {
      eeta <- as.vector(mm_alpha + mm_gamma %*% mean_x)
    } else {
      eeta <- as.vector(mm_alpha)
    }
  }

  eeta
}

# 2) EETAx
# compute E(ETA|x_i): conditional expected value of latent variable,
#                     given specific value of x_i
# - if no eXo (and GAMMA):
#     E(ETA) = (I-B)^-1 ALPHA
#     we return a matrix of size [nobs x nfac] replicating E(ETA)
# - if eXo and GAMMA:
#     E(ETA|x_i) = (I-B)^-1 ALPHA + (I-B)^-1 GAMMA x_i
#     we return  a matrix of size [nobs x nfac]
#
lav_lisrel_eetax <- function(mlist = NULL, exo = NULL, n = nrow(exo),
                             sample_mean = NULL,
                             ov_y_dummy_ov_idx = NULL,
                             ov_x_dummy_ov_idx = NULL,
                             ov_y_dummy_lv_idx = NULL,
                             ov_x_dummy_lv_idx = NULL) {
  mm_lambda <- mlist$lambda
  mm_beta <- mlist$beta
  mm_gamma <- mlist$gamma
  nfac <- ncol(mm_lambda)
  # if eXo, N must be nrow(eXo)
  if (!is.null(exo)) {
    n <- nrow(exo)
  }

  # ALPHA?
  mm_alpha <- lav_lisrel_alpha0(
    mlist = mlist, sample_mean = sample_mean,
    ov_y_dummy_ov_idx = ov_y_dummy_ov_idx,
    ov_x_dummy_ov_idx = ov_x_dummy_ov_idx,
    ov_y_dummy_lv_idx = ov_y_dummy_lv_idx,
    ov_x_dummy_lv_idx = ov_x_dummy_lv_idx
  )

  # construct [nobs x nfac] matrix (repeating ALPHA)
  eeta <- matrix(mm_alpha, n, nfac, byrow = TRUE)

  # put back eXo values if dummy
  if (length(ov_x_dummy_lv_idx) > 0L) {
    eeta[, ov_x_dummy_lv_idx] <- exo
  }

  # BETA?
  if (!is.null(mm_beta)) {
    ib_inv <- lav_lisrel_ibinv(mlist = mlist)
    eeta <- eeta %*% t(ib_inv)
  }

  # GAMMA?
  if (!is.null(mm_gamma)) {
    if (!is.null(mm_beta)) {
      eeta <- eeta + exo %*% t(ib_inv %*% mm_gamma)
    } else {
      eeta <- eeta + exo %*% t(mm_gamma)
    }
  }

  eeta
}

# 3) VETA
# compute V(ETA): variances/covariances of latent variables
# - if no eXo (and GAMMA)
#     V(ETA) = (I-B)^-1 PSI (I-B)^-T
# - if eXo and GAMMA: (cfr lisrel submodel 3a with ksi=x)
#     V(ETA) = (I-B)^-1 [ GAMMA  cov.x t(GAMMA) + PSI] (I-B)^-T
lav_lisrel_veta <- function(mlist = NULL) {
  # mm_lambda <- mlist$lambda
  # nvar <- nrow(mm_lambda)
  mm_psi <- mlist$psi
  # mm_theta <- mlist$theta
  mm_beta <- mlist$beta
  mm_gamma <- mlist$gamma

  if (!is.null(mm_gamma)) {
    cov_x <- mlist$cov.x
    # we treat 'x' as 'ksi' in the LISREL model; cov.x is PHI
    mm_psi <- tcrossprod(mm_gamma %*% cov_x, mm_gamma) + mm_psi
  }

  # beta?
  if (is.null(mm_beta)) {
    veta <- mm_psi
  } else {
    ib_inv <- lav_lisrel_ibinv(mlist = mlist)
    veta <- tcrossprod(ib_inv %*% mm_psi, ib_inv)
  }

  veta
}

# 4) VETAx
# compute V(ETA|x_i): variances/covariances of latent variables
#     V(ETA) = (I-B)^-1 PSI (I-B)^-T  + remove dummies
lav_lisrel_vetax <- function(mlist = NULL, lv_dummy_idx = NULL) {
  mm_psi <- mlist$psi
  mm_beta <- mlist$beta

  # beta?
  if (is.null(mm_beta)) {
    veta <- mm_psi
  } else {
    ib_inv <- lav_lisrel_ibinv(mlist = mlist)
    veta <- tcrossprod(ib_inv %*% mm_psi, ib_inv)
  }

  # remove dummy lv?
  if (!is.null(lv_dummy_idx)) {
    veta <- veta[-lv_dummy_idx, -lv_dummy_idx, drop = FALSE]
  }

  veta
}


# Y
# 1) EY
# 2) EYx
# 3) EYetax
# 4) VY
# 5) VYx
# 6) VYetax

# 1) EY
# compute E(Y): expected value of observed
# E(Y) = NU + LAMBDA %*% E(eta)
#      = NU + LAMBDA %*% (IB.inv %*% ALPHA) # no exo, no GAMMA
#      = NU + LAMBDA %*% (IB.inv %*% ALPHA + IB.inv %*% GAMMA %*% mean.x) # eXo
# if DELTA -> E(Y) = delta * E(Y)
#
# this is similar to lav_model_mu but:
# - we ALWAYS compute NU+ALPHA, even if meanstructure=FALSE
# - never used if GAMMA, since we then have categorical variables, and the
#   'part 1' structure contains the (thresholds +) intercepts, not
#   the means
lav_lisrel_ey <- function(mlist = NULL, mean_x = NULL, sample_mean = NULL,
                          ov_y_dummy_ov_idx = NULL,
                          ov_x_dummy_ov_idx = NULL,
                          ov_y_dummy_lv_idx = NULL,
                          ov_x_dummy_lv_idx = NULL, delta = TRUE) {
  mm_lambda <- mlist$lambda

  # get NU, but do not 'fix'
  mm_nu <- lav_lisrel_nu0(
    mlist = mlist, sample_mean = sample_mean,
    ov_y_dummy_ov_idx = ov_y_dummy_ov_idx,
    ov_x_dummy_ov_idx = ov_x_dummy_ov_idx,
    ov_y_dummy_lv_idx = ov_y_dummy_lv_idx,
    ov_x_dummy_lv_idx = ov_x_dummy_lv_idx
  )

  # compute E(ETA)
  eeta <- lav_lisrel_eeta(
    mlist = mlist, sample_mean = sample_mean,
    mean_x = mean_x,
    ov_y_dummy_ov_idx = ov_y_dummy_ov_idx,
    ov_x_dummy_ov_idx = ov_x_dummy_ov_idx,
    ov_y_dummy_lv_idx = ov_y_dummy_lv_idx,
    ov_x_dummy_lv_idx = ov_x_dummy_lv_idx
  )

  # EY
  ey <- as.vector(mm_nu) + as.vector(mm_lambda %*% eeta)

  # if delta, scale
  if (delta && !is.null(mlist$delta)) {
    ey <- ey * as.vector(mlist$delta)
  }

  ey
}

# 3) EYetax
# compute E(Y|eta_i,x_i): conditional expected value of observed variable
#                         given specific value of eta_i AND x_i
#
# E(y*_i|eta_i, x_i) = NU + LAMBDA eta_i + KAPPA x_i
#
# where eta_i = predict(fit) = factor scores OR specific values for eta_i
# (as in GH integration)
#
# if nexo = 0, and eta_i is single row, YHAT is the same for each observation
# in this case, we return a single row, unless Nobs > 1L, in which case
# we return Nobs identical rows
#
# NOTE: we assume that any effect of x_i on eta_i has already been taken
#       care off

# categorical version
lav_lisrel_eyetax <- function(mlist = NULL,
                              exo = NULL,
                              eta = NULL,
                              n = nrow(exo),
                              sample_mean = NULL,
                              ov_y_dummy_ov_idx = NULL,
                              ov_x_dummy_ov_idx = NULL,
                              ov_y_dummy_lv_idx = NULL,
                              ov_x_dummy_lv_idx = NULL,
                              delta = TRUE) {
  mm_lambda <- mlist$lambda
  mm_beta <- mlist$beta
  if (!is.null(exo)) {
    n <- nrow(exo)
  } else if (!is.null(n)) {
    # nothing to do
  } else {
    n <- 1L
  }

  # create ETA matrix
  if (nrow(eta) == 1L) {
    eta <- matrix(eta, n, ncol(eta), byrow = TRUE)
  }

  # always augment ETA with 'dummy values' (0 for ov.y, eXo for ov.x)
  # ndummy <- length(c(ov.y.dummy.lv.idx, ov.x.dummy.lv.idx))
  # if(ndummy > 0L) {
  #    ETA2 <- cbind(ETA, matrix(0, N, ndummy))
  # } else {
  eta2 <- eta
  # }

  # only if we have dummy ov.y, we need to compute the 'yhat' values
  # beforehand
  if (length(ov_y_dummy_lv_idx) > 0L) {
    # insert eXo values
    if (length(ov_x_dummy_lv_idx) > 0L) {
      eta2[, ov_x_dummy_lv_idx] <- exo
    }
    # zero ov.y values
    if (length(ov_y_dummy_lv_idx) > 0L) {
      eta2[, ov_y_dummy_lv_idx] <- 0
    }

    # ALPHA? (reconstruct, but no 'fix')
    mm_alpha <- lav_lisrel_alpha0(
      mlist = mlist, sample_mean = sample_mean,
      ov_y_dummy_ov_idx = ov_y_dummy_ov_idx,
      ov_x_dummy_ov_idx = ov_x_dummy_ov_idx,
      ov_y_dummy_lv_idx = ov_y_dummy_lv_idx,
      ov_x_dummy_lv_idx = ov_x_dummy_lv_idx
    )
    # BETA?
    if (!is.null(mm_beta)) {
      eta2 <- sweep(tcrossprod(eta2, mm_beta), 2L, STATS = mm_alpha, FUN = "+")
    } else {
      eta2 <- sweep(eta2, 2L, STATS = mm_alpha, FUN = "+")
    }

    # put back eXo values
    if (length(ov_x_dummy_lv_idx) > 0L) {
      eta2[, ov_x_dummy_lv_idx] <- exo
    }

    # put back ETA values for the 'real' latent variables
    dummy_idx <- c(ov_x_dummy_lv_idx, ov_y_dummy_lv_idx)
    if (length(dummy_idx) > 0L) {
      lv_regular_idx <- seq_len(min(dummy_idx) - 1L)
      eta2[, lv_regular_idx] <- eta[, lv_regular_idx, drop = FALSE]
    }
  }

  # get NU, but do not 'fix'
  mm_nu <- lav_lisrel_nu0(
    mlist = mlist,
    sample_mean = sample_mean,
    ov_y_dummy_ov_idx = ov_y_dummy_ov_idx,
    ov_x_dummy_ov_idx = ov_x_dummy_ov_idx,
    ov_y_dummy_lv_idx = ov_y_dummy_lv_idx,
    ov_x_dummy_lv_idx = ov_x_dummy_lv_idx
  )

  # EYetax
  eyetax <- sweep(tcrossprod(eta2, mm_lambda), 2L, STATS = mm_nu, FUN = "+")

  # if delta, scale
  if (delta && !is.null(mlist$delta)) {
    eyetax <- sweep(eyetax, 2L, STATS = mlist$delta, FUN = "*")
  }

  eyetax
}

# unconditional version
lav_lisrel_eyetax3 <- function(mlist = NULL,
                               eta = NULL,
                               sample_mean = NULL,
                               mean_x = NULL,
                               ov_y_dummy_ov_idx = NULL,
                               ov_x_dummy_ov_idx = NULL,
                               ov_y_dummy_lv_idx = NULL,
                               ov_x_dummy_lv_idx = NULL,
                               delta = TRUE) {
  mm_lambda <- mlist$lambda

  # special case: empty lambda
  if (ncol(mm_lambda) == 0L) {
    return(matrix(sample_mean,
      nrow(eta), length(sample_mean),
      byrow = TRUE
    ))
  }

  # lv idx
  dummy_idx <- c(ov_y_dummy_lv_idx, ov_x_dummy_lv_idx)
  if (length(dummy_idx) > 0L) {
    nondummy_idx <- seq_len(min(dummy_idx) - 1L)
  } else {
    nondummy_idx <- seq_len(ncol(mlist$lambda))
  }

  # beta?
  if (is.null(mlist$beta) || length(ov_y_dummy_lv_idx) == 0L ||
    length(nondummy_idx) == 0L) {
    lambda__ib_inv <- mm_lambda
  } else {
    # only keep those columns of BETA that correspond to the
    # the `regular' latent variables
    # (ie. ignore the structural part altogether)
    mlist2 <- mlist
    mlist2$beta[, dummy_idx] <- 0
    ib_inv <- lav_lisrel_ibinv(mlist = mlist2)
    lambda__ib_inv <- mm_lambda %*% ib_inv
  }

  # compute model-implied means
  ey <- lav_lisrel_ey(
    mlist = mlist, mean_x = mean_x,
    sample_mean = sample_mean,
    ov_y_dummy_ov_idx = ov_y_dummy_ov_idx,
    ov_x_dummy_ov_idx = ov_x_dummy_ov_idx,
    ov_y_dummy_lv_idx = ov_y_dummy_lv_idx,
    ov_x_dummy_lv_idx = ov_x_dummy_lv_idx
  )

  eeta <- lav_lisrel_eeta(
    mlist = mlist, sample_mean = sample_mean,
    mean_x = mean_x,
    ov_y_dummy_ov_idx = ov_y_dummy_ov_idx,
    ov_x_dummy_ov_idx = ov_x_dummy_ov_idx,
    ov_y_dummy_lv_idx = ov_y_dummy_lv_idx,
    ov_x_dummy_lv_idx = ov_x_dummy_lv_idx
  )

  # center regular lv only
  eta[, nondummy_idx] <- sweep(eta[, nondummy_idx, drop = FALSE], 2L,
    STATS = eeta[nondummy_idx], FUN = "-"
  )

  # project from lv to ov, if we have any lv
  if (length(nondummy_idx) > 0) {
    eyetax <- sweep(
      tcrossprod(
        eta[, nondummy_idx, drop = FALSE],
        lambda__ib_inv[, nondummy_idx, drop = FALSE]
      ),
      2L,
      STATS = ey, FUN = "+"
    )
  } else {
    eyetax <- eta
  }

  # put back eXo variables
  if (length(ov_x_dummy_lv_idx) > 0L) {
    eyetax[, ov_x_dummy_ov_idx] <- eta[, ov_x_dummy_lv_idx, drop = FALSE]
  }

  # if delta, scale
  if (delta && !is.null(mlist$delta)) {
    eyetax <- sweep(eyetax, 2L, STATS = mlist$delta, FUN = "*")
  }

  eyetax
}

# 4) VY
# compute the *un*conditional variance/covariance of y: V(Y) or V(Y*)
# 'unconditional' model-implied (co)variances
#  - same as Sigma.hat if all Y are continuous
#  - diagonal is 1.0 (or delta^2) if categorical
#  - if also Gamma, cov.x is used (only if conditional.x)
#    only in THIS case, VY is different from diag(VYx)
#
# V(Y) = LAMBDA V(ETA) t(LAMBDA) + THETA
lav_lisrel_vy <- function(mlist = NULL) {
  mm_lambda <- mlist$lambda
  mm_theta <- mlist$theta

  veta <- lav_lisrel_veta(mlist = mlist)
  vy <- tcrossprod(mm_lambda %*% veta, mm_lambda) + mm_theta
  vy
}

# 5) VYx
# compute V(Y*|x_i) == model-implied covariance matrix
# this equals V(Y*) if no (explicit) eXo no GAMMA
#
# in >0.6-20: special treatment for composites
#
lav_lisrel_sigma <- function(mlist = NULL, delta = TRUE) {
  mm_lambda <- mlist$lambda
  nvar <- nrow(mm_lambda)
  mm_psi <- mlist$psi
  mm_theta <- mlist$theta
  mm_beta <- mlist$beta
  mm_wmat <- mlist$wmat

  # standard: no composites
  if (is.null(mm_wmat)) {
    # beta?
    if (is.null(mm_beta)) {
      lambda__ib_inv <- mm_lambda
    } else {
      ib_inv <- lav_lisrel_ibinv(mlist = mlist)
      lambda__ib_inv <- mm_lambda %*% ib_inv
    }
    # compute V(Y*|x_i)
    vyx <- tcrossprod(lambda__ib_inv %*% mm_psi, lambda__ib_inv) + mm_theta

    # composites, or mix of composites and latent variables
  } else {
    # - first join LAMBDA and WMAT
    # - create 'T' matrix: - identity for regular lv's,
    #                      - THETA block-diagonal for composites
    # - create C_0: VETA, but zero diagonal elements for composites
    cov_idx <- which(apply(
      mm_lambda, 1L,
      function(x) sum(x == 0) == ncol(mm_lambda)
    ))
    clv_idx <- which(apply(
      mm_lambda, 2L,
      function(x) sum(x == 0) == nrow(mm_lambda)
    ))
    # regular latent variables
    # rlv_idx <- seq_len(ncol(mm_lambda))[-clv_idx]

    # combine LAMBDA and WMAT
    lw <- mm_lambda + mm_wmat

    tmat <- diag(nrow(mm_lambda))
    tmat[cov_idx, cov_idx] <- mm_theta[cov_idx, cov_idx]
    wtw <- t(lw[, clv_idx, drop = FALSE]) %*% tmat %*%
                               lw[, clv_idx, drop = FALSE]
    wtw_inv <- solve(wtw)
    wtw_inv_1 <- diag(ncol(mm_lambda))
    wtw_inv_1[clv_idx, clv_idx] <- wtw_inv

    if (is.null(mm_beta)) {
      ib_inv <- diag(nrow(mm_psi))
    } else {
      ib_inv <- lav_lisrel_ibinv(mlist = mlist)
    }
    veta <- ib_inv %*% mm_psi %*% t(ib_inv)
    c0 <- veta
    diag(c0)[clv_idx] <- 0

    vyx <- tmat %*% lw %*% wtw_inv_1 %*% c0 %*%
                       t(wtw_inv_1) %*% t(lw) %*% tmat + mm_theta
  }

  # if delta, scale
  if (delta && !is.null(mlist$delta)) {
    mm_delta <- diag(mlist$delta[, 1L], nrow = nvar, ncol = nvar)
    vyx <- mm_delta %*% vyx %*% mm_delta
  }

  vyx
}

lav_lisrel_implied_fast <- function(mlist = NULL, th_idx = NULL,
                                    need_sigma = FALSE,
                                    need_mu = FALSE,
                                    need_th = FALSE,
                                    need_pi = FALSE,
                                    delta = TRUE) {
  mm_lambda <- mlist$lambda
  nvar <- nrow(mm_lambda)
  mm_beta <- mlist$beta
  mm_wmat <- mlist$wmat

  out <- list()

  # sigma, mu, thresholds, and slopes all use LAMBDA %*% (I - BETA)^-1
  # in the common LISREL path. Compute it once for this block.
  need_lambda_ib_inv <- (need_sigma && is.null(mm_wmat)) ||
    need_mu || need_th || need_pi
  ib_inv <- NULL
  lambda__ib_inv <- NULL
  if (need_lambda_ib_inv) {
    if (is.null(mm_beta)) {
      lambda__ib_inv <- mm_lambda
    } else {
      ib_inv <- lav_lisrel_ibinv(mlist = mlist)
      lambda__ib_inv <- mm_lambda %*% ib_inv
    }
  }

  if (need_sigma) {
    mm_psi <- mlist$psi
    mm_theta <- mlist$theta

    # standard: no composites
    if (is.null(mm_wmat)) {
      vyx <- tcrossprod(lambda__ib_inv %*% mm_psi, lambda__ib_inv) + mm_theta

      # composites, or mix of composites and latent variables
    } else {
      cov_idx <- which(apply(
        mm_lambda, 1L,
        function(x) sum(x == 0) == ncol(mm_lambda)
      ))
      clv_idx <- which(apply(
        mm_lambda, 2L,
        function(x) sum(x == 0) == nrow(mm_lambda)
      ))

      # combine LAMBDA and WMAT
      lw <- mm_lambda + mm_wmat

      tmat <- diag(nrow(mm_lambda))
      tmat[cov_idx, cov_idx] <- mm_theta[cov_idx, cov_idx]
      wtw <- t(lw[, clv_idx, drop = FALSE]) %*% tmat %*%
                                 lw[, clv_idx, drop = FALSE]
      wtw_inv <- solve(wtw)
      wtw_inv_1 <- diag(ncol(mm_lambda))
      wtw_inv_1[clv_idx, clv_idx] <- wtw_inv

      if (is.null(mm_beta)) {
        ib_inv_sigma <- diag(nrow(mm_psi))
      } else {
        if (is.null(ib_inv)) {
          ib_inv <- lav_lisrel_ibinv(mlist = mlist)
        }
        ib_inv_sigma <- ib_inv
      }
      veta <- ib_inv_sigma %*% mm_psi %*% t(ib_inv_sigma)
      c0 <- veta
      diag(c0)[clv_idx] <- 0

      vyx <- tmat %*% lw %*% wtw_inv_1 %*% c0 %*%
                         t(wtw_inv_1) %*% t(lw) %*% tmat + mm_theta
    }

    # if delta, scale
    if (delta && !is.null(mlist$delta)) {
      mm_delta <- diag(mlist$delta[, 1L], nrow = nvar, ncol = nvar)
      vyx <- mm_delta %*% vyx %*% mm_delta
    }

    out$sigma <- vyx
  }

  if (need_mu) {
    mm_nu <- mlist$nu
    mm_alpha <- mlist$alpha

    # shortcut
    if (is.null(mm_alpha) || is.null(mm_nu)) {
      out$mu <- matrix(0, nrow(mm_lambda), 1L)
    } else {
      out$mu <- mm_nu + lambda__ib_inv %*% mm_alpha
    }
  }

  if (need_th) {
    mm_tau <- mlist$tau
    nth <- nrow(mm_tau)

    # missing alpha
    if (is.null(mlist$alpha)) {
      mm_alpha <- matrix(0, ncol(mm_lambda), 1L)
    } else {
      mm_alpha <- mlist$alpha
    }

    # missing nu
    if (is.null(mlist$nu)) {
      mm_nu <- matrix(0, nvar, 1L)
    } else {
      mm_nu <- mlist$nu
    }

    if (is.null(th_idx)) {
      th_idx <- seq_len(nth)
      nlev <- rep(1L, nvar)
      k_nu <- diag(nvar)
    } else {
      nlev <- tabulate(th_idx, nbins = nvar)
      nlev[nlev == 0L] <- 1L
      k_nu <- matrix(0, sum(nlev), nvar)
      k_nu[cbind(seq_len(sum(nlev)), rep(seq_len(nvar), times = nlev))] <- 1.0
    }

    # shortcut
    if (is.null(mm_tau)) {
      out$th <- matrix(0, length(th_idx), 1L)
    } else {
      pi0 <- mm_nu + lambda__ib_inv %*% mm_alpha

      # interleave th's with zeros where we have numeric variables
      th <- numeric(length(th_idx))
      th[th_idx > 0L] <- mm_tau[, 1L]

      th_1 <- th - (k_nu %*% pi0)

      # if delta, scale
      if (delta && !is.null(mlist$delta)) {
        delta_diag <- mlist$delta[, 1L]
        delta_star_diag <- rep(delta_diag, times = nlev)
        th_1 <- th_1 * delta_star_diag
      }

      out$th <- as.vector(th_1)
    }
  }

  if (need_pi) {
    mm_gamma <- mlist$gamma

    # shortcut
    if (is.null(mm_gamma)) {
      out$pi <- matrix(0, nrow(mm_lambda), 0L)
    } else {
      pi0 <- lambda__ib_inv %*% mm_gamma

      # if delta, scale
      if (delta && !is.null(mlist$delta)) {
        delta_diag <- mlist$delta[, 1L]
        pi0 <- pi0 * delta_diag
      }

      out$pi <- pi0
    }
  }

  out
}


### compute model-implied sample statistics
#
# 1) MuHat (similar to EY, but continuous only)
# 2) TH
# 3) PI
# 4) SigmaHat == VYx

# compute MuHat for a single block/group; only for the continuous case (no eXo)
#
# this is a special case of E(Y) where
# - we have no (explicit) eXogenous variables
# - only continuous
lav_lisrel_mu <- function(mlist = NULL) {
  mm_nu <- mlist$nu
  mm_alpha <- mlist$alpha
  mm_lambda <- mlist$lambda
  mm_beta <- mlist$beta

  # shortcut
  if (is.null(mm_alpha) || is.null(mm_nu)) {
    return(matrix(0, nrow(mm_lambda), 1L))
  }

  # beta?
  if (is.null(mm_beta)) {
    lambda__ib_inv <- mm_lambda
  } else {
    ib_inv <- lav_lisrel_ibinv(mlist = mlist)
    lambda__ib_inv <- mm_lambda %*% ib_inv
  }

  # compute Mu Hat
  mu_hat <- mm_nu + lambda__ib_inv %*% mm_alpha

  mu_hat
}

# compute TH for a single block/group
lav_lisrel_th <- function(mlist = NULL, th_idx = NULL, delta = TRUE) {
  mm_lambda <- mlist$lambda
  nvar <- nrow(mm_lambda)
  nfac <- ncol(mm_lambda)
  mm_beta <- mlist$beta
  mm_tau <- mlist$tau
  nth <- nrow(mm_tau)

  # missing alpha
  if (is.null(mlist$alpha)) {
    mm_alpha <- matrix(0, nfac, 1L)
  } else {
    mm_alpha <- mlist$alpha
  }

  # missing nu
  if (is.null(mlist$nu)) {
    mm_nu <- matrix(0, nvar, 1L)
  } else {
    mm_nu <- mlist$nu
  }

  if (is.null(th_idx)) {
    th_idx <- seq_len(nth)
    nlev <- rep(1L, nvar)
    k_nu <- diag(nvar)
  } else {
    nlev <- tabulate(th_idx, nbins = nvar)
    nlev[nlev == 0L] <- 1L
    k_nu <- matrix(0, sum(nlev), nvar)
    k_nu[cbind(seq_len(sum(nlev)), rep(seq_len(nvar), times = nlev))] <- 1.0
  }

  # shortcut
  if (is.null(mm_tau)) {
    return(matrix(0, length(th_idx), 1L))
  }

  # beta?
  if (is.null(mm_beta)) {
    lambda__ib_inv <- mm_lambda
  } else {
    ib_inv <- lav_lisrel_ibinv(mlist = mlist)
    lambda__ib_inv <- mm_lambda %*% ib_inv
  }

  # compute pi0
  pi0 <- mm_nu + lambda__ib_inv %*% mm_alpha

  # interleave th's with zeros where we have numeric variables
  th <- numeric(length(th_idx))
  th[th_idx > 0L] <- mm_tau[, 1L]

  # compute TH
  th_1 <- th - (k_nu %*% pi0)

  # if delta, scale
  if (delta && !is.null(mlist$delta)) {
    delta_diag <- mlist$delta[, 1L]
    delta_star_diag <- rep(delta_diag, times = nlev)
    th_1 <- th_1 * delta_star_diag
  }

  as.vector(th_1)
}

# compute PI for a single block/group
lav_lisrel_pi <- function(mlist = NULL, delta = TRUE) {
  mm_lambda <- mlist$lambda
  mm_beta <- mlist$beta
  mm_gamma <- mlist$gamma

  # shortcut
  if (is.null(mm_gamma)) {
    return(matrix(0, nrow(mm_lambda), 0L))
  }

  # beta?
  if (is.null(mm_beta)) {
    lambda__ib_inv <- mm_lambda
  } else {
    ib_inv <- lav_lisrel_ibinv(mlist = mlist)
    lambda__ib_inv <- mm_lambda %*% ib_inv
  }

  # compute PI
  pi0 <- lambda__ib_inv %*% mm_gamma

  # if delta, scale
  if (delta && !is.null(mlist$delta)) {
    delta_diag <- mlist$delta[, 1L]
    pi0 <- pi0 * delta_diag
  }

  pi0
}

lav_lisrel_lambda <- function(mlist = NULL,
                              ov_y_dummy_ov_idx = NULL,
                              ov_x_dummy_ov_idx = NULL,
                              ov_y_dummy_lv_idx = NULL,
                              ov_x_dummy_lv_idx = NULL,
                              remove_dummy_lv = FALSE) {
  lv_dummy_idx <- c(ov_y_dummy_lv_idx, ov_x_dummy_lv_idx)

  # fix LAMBDA
  mm_lambda <- mlist$lambda
  if (length(ov_y_dummy_ov_idx) > 0L && !is.null(mlist$beta)) {
    mm_lambda[ov_y_dummy_ov_idx, ] <- mlist$beta[ov_y_dummy_lv_idx, ]
  }

  # remove dummy lv?
  if (remove_dummy_lv && length(lv_dummy_idx) > 0L) {
    mm_lambda <- mm_lambda[, -lv_dummy_idx, drop = FALSE]
  }

  mm_lambda
}

lav_lisrel_theta <- function(mlist = NULL,
                             ov_y_dummy_ov_idx = NULL,
                             ov_x_dummy_ov_idx = NULL,
                             ov_y_dummy_lv_idx = NULL,
                             ov_x_dummy_lv_idx = NULL) {
  ov_dummy_idx <- c(ov_y_dummy_ov_idx, ov_x_dummy_ov_idx)
  lv_dummy_idx <- c(ov_y_dummy_lv_idx, ov_x_dummy_lv_idx)

  # fix THETA
  mm_theta <- mlist$theta
  if (length(ov_dummy_idx) > 0L) {
    mm_theta[ov_dummy_idx, ov_dummy_idx] <-
      mlist$psi[lv_dummy_idx, lv_dummy_idx]
  }

  mm_theta
}

# compute (I - BETA)^{-1}
# new in 0.6-22: check structure of BETA
#  1. BETA absent / all-zero        -> identity
#  1b.BETA is complex               -> general solve()
#  2. BETA strictly lower triangular-> forwardsolve
#  3. BETA strictly upper triangular-> backsolve
#  4. BETA is a DAG in any order    -> Kahn topological sort, permute to
#                                      lower tri, forwardsolve, unpermute
#  5. BETA has directed cycles      -> general solve()
lav_lisrel_ibinv <- function(mlist = NULL) {
  mm_beta <- mlist$beta
  nr   <- nrow(mlist$psi)

  # case 1: no BETA, or BETA is identically zero
  if (is.null(mm_beta) || all(mm_beta == 0)) {
    return(diag(nr))
  }

  # forwardsolve/backsolve do not support complex values; solve() does
  if (is.complex(mm_beta)) {
    tmp <- -mm_beta
    diag(tmp) <- 1
    return(solve(tmp))
  }

  # case 2: strictly lower triangular
  if (all(mm_beta[upper.tri(mm_beta)] == 0)) {
    return(forwardsolve(diag(nr) - mm_beta, diag(nr)))
  }

  # case 3: strictly upper triangular
  if (all(mm_beta[lower.tri(mm_beta)] == 0)) {
    return(backsolve(diag(nr) - mm_beta, diag(nr)))
  }

  # case 4: recursive/DAG
  # BETA[i,j] != 0  <=>  directed edge j -> i
  indegree <- rowSums(mm_beta != 0)
  result <- integer(nr)
  k <- 0L
  queue <- which(indegree == 0)
  while (length(queue) > 0L) {
    v <- queue[1L]
    queue <- queue[-1L]
    k <- k + 1L
    result[k] <- v
    for (u in which(mm_beta[, v] != 0)) {
      indegree[u] <- indegree[u] - 1L
      if (indegree[u] == 0L) queue <- c(queue, u)
    }
  }
  if (k == nr) {
    # recursive model: permute B to strictly lower triangular, solve, unpermute
    ib_inv <- forwardsolve(diag(nr) - mm_beta[result, result], diag(nr))
    inv_order <- order(result)
    return(ib_inv[inv_order, inv_order])
  }

  # case 5: non-recursive model
  tmp <- -mm_beta
  diag(tmp) <- 1
  solve(tmp)
}

# only if ALPHA=NULL but we need it anyway
# we 'reconstruct' ALPHA here (including dummy entries), no fixing
#
# without any dummy variables, this is just the zero vector
# but if we have dummy variables, we need to fill in their values
#
#
lav_lisrel_alpha0 <- function(mlist = NULL, sample_mean = NULL,
                              ov_y_dummy_ov_idx = NULL,
                              ov_x_dummy_ov_idx = NULL,
                              ov_y_dummy_lv_idx = NULL,
                              ov_x_dummy_lv_idx = NULL) {
  if (!is.null(mlist$alpha)) {
    return(mlist$alpha)
  }

  mm_lambda <- mlist$lambda
  nfac <- ncol(mm_lambda)

  ov_dummy_idx <- c(ov_y_dummy_ov_idx, ov_x_dummy_ov_idx)
  lv_dummy_idx <- c(ov_y_dummy_lv_idx, ov_x_dummy_lv_idx)

  if (length(ov_dummy_idx) > 0L) {
    mm_alpha <- matrix(0, nfac, 1L)
    # Note: instead of sample.mean, we need 'intercepts'
    # sample.mean = NU + LAMBDA..IB.inv %*% ALPHA
    # so,
    # solve(LAMBDA..IB.inv) %*% (sample.mean - NU) = ALPHA
    # where
    # - LAMBDA..IB.inv only contains 'dummy' variables, and is square
    # - NU elements are not needed (since not in ov.dummy.idx)
    ib_inv <- lav_lisrel_ibinv(mlist = mlist)
    lambda__ib_inv <- mm_lambda %*% ib_inv
    lambda__ib_inv_dummy <- lambda__ib_inv[ov_dummy_idx, lv_dummy_idx]
    mm_alpha[lv_dummy_idx] <-
      solve(lambda__ib_inv_dummy, sample_mean[ov_dummy_idx])
  } else {
    mm_alpha <- matrix(0, nfac, 1L)
  }

  mm_alpha
}

# only if NU=NULL but we need it anyway
#
# since we have no meanstructure, we can assume NU is unrestricted
# and contains either:
#     1) the sample means (if not eXo)
#     2) the intercepts, if we have exogenous covariates
#        since sample.mean = NU + LAMBDA %*% E(eta)
#        we have NU = sample.mean - LAMBDA %*% E(eta)
lav_lisrel_nu0 <- function(mlist = NULL, sample_mean = NULL,
                           ov_y_dummy_ov_idx = NULL,
                           ov_x_dummy_ov_idx = NULL,
                           ov_y_dummy_lv_idx = NULL,
                           ov_x_dummy_lv_idx = NULL) {
  if (!is.null(mlist$nu)) {
    return(mlist$nu)
  }

  # if nexo > 0, subtract lambda %*% EETA
  if (length(ov_x_dummy_ov_idx) > 0L) {
    eeta <- lav_lisrel_eeta(mlist,
      mean_x = NULL,
      sample_mean = sample_mean,
      ov_y_dummy_ov_idx = ov_y_dummy_ov_idx,
      ov_x_dummy_ov_idx = ov_x_dummy_ov_idx,
      ov_y_dummy_lv_idx = ov_y_dummy_lv_idx,
      ov_x_dummy_lv_idx = ov_x_dummy_lv_idx
    )

    # 'regress' NU on X
    mm_nu <- sample_mean - mlist$lambda %*% eeta

    # just to make sure we have exact zeroes for all dummies
    mm_nu[c(ov_y_dummy_ov_idx, ov_x_dummy_ov_idx)] <- 0
  } else {
    # unrestricted mean
    mm_nu <- sample_mean
  }

  mm_nu
}

# set (total/residual) variances of composites
# and while we at it, also set intercepts of composites
lav_lisrel_comp_set_intresvar <- function(mlist = NULL,
                                          tol = .Machine$double.eps,
                                          debug = FALSE) {
  mm_lambda <- mlist$lambda
  mm_beta <- mlist$beta
  mm_psi <- mlist$psi
  mm_wmat <- mlist$wmat
  # mm_theta <- mlist$theta

  # std.lv or not?
  marker_idx <- lav_utils_get_marker(mlist$wmat)
  std_lv <- FALSE
  if (all(is.na(marker_idx))) {
    std_lv <- TRUE
  }

  # housekeeping
  ovc_idx <- which(apply(
    mm_lambda, 1L,
    function(x) sum(x == 0) == ncol(mm_lambda)
  ))
  lvc_idx <- which(apply(
    mm_lambda, 2L,
    function(x) sum(x == 0) == nrow(mm_lambda)
  ))
  lvc_flag <- logical(nrow(mm_psi))
  lvc_flag[lvc_idx] <- TRUE
  tmat <- diag(nrow(mm_lambda))
  tmat[ovc_idx, ovc_idx] <- mlist$theta[ovc_idx, ovc_idx]

  if (std_lv) {
    target_psi <- rep(1, ncol(mm_wmat))
  } else {
    # total variances composites
    target_psi <- diag(t(mm_wmat) %*% tmat %*% mm_wmat)
  }
  # fill in PSI element for non-composites
  target_psi[!lvc_flag] <- diag(mm_psi)[!lvc_flag]

  # initial values (including exogenous variances)
  diag(mm_psi) <- target_psi

  # no regressions
  if (is.null(mm_beta)) {
    # store PSI
    mlist$psi <- mm_psi

    # fix intercept (if needed)
    if (!is.null(mlist$alpha)) {
      tmp <- t(mm_wmat) %*% mlist$nu
      mlist$alpha[lvc_idx, 1L] <- tmp[lvc_idx, 1L]
    }

    return(mlist)
  }

  # set (residual) variances in psi
  abs_beta <- abs(mlist$beta)
  # x_idx <- which(apply(abs_beta, 1L, sum) == 0 & lvc_flag)
  y_idx <- which(apply(abs_beta, 1L, sum) != 0 & lvc_flag)

  # compute IB.inv
  ib_inv <- lav_lisrel_ibinv(mlist)

  # check if BETA is acyclic
  if (!lav_graph_is_acyclic(mm_beta)) {
    # damn, we have a cyclic model; use nlminb()
    mm_psi <- lav_mlist_target_psi(
      ib_inv = ib_inv, mm_psi = mm_psi,
      target_psi = target_psi, y_idx = y_idx
    )
  } else {
    # for an acyclic model, we should be able to find the
    # residual analytically; simply by computing the model-based
    # total variances of the RHS of each regression, and set the
    # residual of the y variable so that it is exactly equal to unity
    # (yes, this will result in negative residual variances if needed)

    # ideally, we first sort the variables 'in topological order'; then
    # we need only one run; but here we are somewhat lazy, and we
    # use a few runs, each time setting more variables right

    # get ancestors list for each node/variable
    ancestors <- lav_graph_get_ancestors(mm_beta)

    # for each y variable, compute IB.inv %*% psi %*% t(IB.inv), without y
    ny <- length(y_idx)
    max_rep <- ny * 4
    for (rep in seq_len(max_rep)) {
      if (debug) {
        cat("rep = ", rep, "\n")
      }
      # check current diagonal
      current_diag <- diag(ib_inv %*% mm_psi %*% t(ib_inv))
      if (debug) {
        cat("target.psi   = ", target_psi[y_idx], "\n")
        cat("current.diag = ", current_diag[y_idx], "\n")
      }
      if (all(abs(current_diag[y_idx] - target_psi[y_idx]) < tol)) {
        # we are done, bail out
        break
      }
      for (i in seq_len(ny)) {
        this_y_idx <- y_idx[i]
        this_x_idx <- ancestors[[this_y_idx]]
        ib_inv_y <- ib_inv[this_y_idx, this_x_idx, drop = FALSE]
        psi_x <- mm_psi[this_x_idx, this_x_idx, drop = FALSE]
        var_y <- drop(ib_inv_y %*% psi_x %*% t(ib_inv_y))
        mm_psi[this_y_idx, this_y_idx] <- target_psi[this_y_idx] - var_y
      }
    }

    # final check?
    current_diag <- diag(ib_inv %*% mm_psi %*% t(ib_inv))
    if (debug) {
      cat("final current.diag = ", current_diag[y_idx], "\n")
    }
    # don't be too strict here
    if (any(abs(current_diag[y_idx] - target_psi[y_idx]) > sqrt(tol))) {
      # as a last resort, use optimization
      mm_psi <- lav_mlist_target_psi(
        ib_inv = ib_inv, mm_psi = mm_psi,
        target_psi = target_psi, y_idx = y_idx
      )
    }
  } # acyclic

  # store PSI
  mlist$psi <- mm_psi

  # fix composite mean (if needed)
  if (!is.null(mlist$alpha)) {
    tmp <- t(mm_wmat) %*% mlist$nu # total means
    # step 1: set exogenous alpha values
    mlist$alpha[lvc_idx, 1L] <- tmp[lvc_idx, 1L]
    # step 2: set endogenous alpha values
    exo_only <- mlist$alpha
    exo_only[y_idx] <- 0
    # compensate alpha so that total mean stays the same
    tmp2 <- ib_inv %*% exo_only
    mlist$alpha[y_idx, 1L] <- mlist$alpha[y_idx, 1L] - tmp2[y_idx, 1L]
  }

  mlist
}

# if DELTA parameterization, compute residual elements (in theta, or psi)
# - typically for (endogenous) observed *categorical* variables only
# - but could also be all (endogenous) observed variables, if correlation = TRUE
# - or all (endogenous) latent and observed variables, if ov.only = FALSE
#
# new version YR 29 Oct 2024: try harder for psi elements (but this only works
#                             for acyclic models)
#             YR 01 Nov 2024: for non-acyclic models: use optimization
lav_lisrel_residual_variances <- function(mlist = NULL,
                                          num_idx = NULL,
                                          ov_y_dummy_ov_idx = NULL,
                                          ov_y_dummy_lv_idx = NULL,
                                          ov_only = TRUE,
                                          tol = .Machine$double.eps,
                                          debug = FALSE) {
  mm_beta <- mlist$beta
  mm_psi <- mlist$psi
  if (is.null(mlist$delta)) {
    delta <- rep(1, nrow(mlist$lambda))
  } else {
    delta <- mlist$delta
  }

  # remove num.idx from ov.y.dummy.*
  if (length(num_idx) > 0L && length(ov_y_dummy_ov_idx) > 0L) {
    n_idx <- which(ov_y_dummy_ov_idx %in% num_idx)
    if (length(n_idx) > 0L) {
      ov_y_dummy_ov_idx <- ov_y_dummy_ov_idx[-n_idx]
      ov_y_dummy_lv_idx <- ov_y_dummy_lv_idx[-n_idx]
    }
  }

  # if delta, the target may not be unity, but DELTA^(-2)
  target_all <- 1 / (delta * delta) # often the unit vector
  target_psi <- rep(1, nrow(mm_psi))
  if (length(ov_y_dummy_ov_idx) > 0L) {
    target_psi[ov_y_dummy_lv_idx] <- target_all[ov_y_dummy_ov_idx]
  }

  # phase 1: set (residual) variances in psi
  if (!is.null(mm_beta) && length(ov_y_dummy_ov_idx) > 0L) {
    abs_beta <- abs(mlist$beta)
    x_idx <- which(apply(abs_beta, 1L, sum) == 0)
    if (ov_only) {
      y_idx <- ov_y_dummy_lv_idx
      # remove x.idx elements (fixed.x = FALSE)
      if (any(y_idx %in% x_idx)) {
        y_idx <- y_idx[-which(y_idx %in% x_idx)]
      }
    } else {
      y_idx <- which(apply(abs_beta, 1L, sum) != 0)
    }

    nr <- nrow(mm_beta)
    ib <- -mm_beta
    ib[lav_matrix_diag_idx(nr)] <- 1
    ib_inv <- solve(ib)


    # check if BETA is acyclic
    if (!lav_graph_is_acyclic(mm_beta)) {
      # damn, we have a cyclic model; use nlminb()
      mm_psi <- lav_mlist_target_psi(
        ib_inv = ib_inv, mm_psi = mm_psi,
        target_psi = target_psi, y_idx = y_idx
      )
    } else {
      # for an acyclic model, we should be able to find the
      # residual analytically; simply by computing the model-based
      # total variances of the RHS of each regression, and set the
      # residual of the y variable so that it is exactly equal to unity
      # (yes, this will result in negative residual variances if needed)

      # ideally, we first sort the variables 'in topological order'; then
      # we need only one run; but here we are somewhat lazy, and we
      # use a few runs, each time setting more variables right

      # get ancestors list for each node/variable
      ancestors <- lav_graph_get_ancestors(mm_beta)

      # for each y variable, compute IB.inv %*% psi %*% t(IB.inv), without y
      ny <- length(y_idx)
      max_rep <- ny * 4
      for (rep in seq_len(max_rep)) {
        if (debug) {
          cat("rep = ", rep, "\n")
        }
        # check current diagonal
        current_diag <- diag(ib_inv %*% mm_psi %*% t(ib_inv))
        if (debug) {
          cat("target.psi   = ", target_psi[y_idx], "\n")
          cat("current.diag = ", current_diag[y_idx], "\n")
        }
        if (all(abs(current_diag[y_idx] - target_psi[y_idx]) < tol)) {
          # we are done, bail out
          break
        }
        for (i in seq_len(ny)) {
          this_y_idx <- y_idx[i]
          this_x_idx <- ancestors[[this_y_idx]]
          ib_inv_y <- ib_inv[this_y_idx, this_x_idx, drop = FALSE]
          psi_x <- mm_psi[this_x_idx, this_x_idx, drop = FALSE]
          var_y <- drop(ib_inv_y %*% psi_x %*% t(ib_inv_y))
          mm_psi[this_y_idx, this_y_idx] <- target_psi[this_y_idx] - var_y
        }
      }

      # final check?
      current_diag <- diag(ib_inv %*% mm_psi %*% t(ib_inv))
      if (debug) {
        cat("final current.diag = ", current_diag[y_idx], "\n")
      }
      # don't be too strict here
      if (any(abs(current_diag[y_idx] - target_psi[y_idx]) > sqrt(tol))) {
        # as a last resort, use optimization
        mm_psi <- lav_mlist_target_psi(
          ib_inv = ib_inv, mm_psi = mm_psi,
          target_psi = target_psi, y_idx = y_idx
        )
      }
    } # acyclic
  } # phase 1

  # store PSI
  mlist$psi <- mm_psi

  # phase 2: set residual variances in theta

  # force non-numeric theta elements to be zero
  if (length(num_idx) > 0L) {
    diag(mlist$theta)[-num_idx] <- 0.0
  } else {
    diag(mlist$theta) <- 0.0
  }

  sigma_hat <- lav_lisrel_sigma(mlist = mlist, delta = FALSE)
  diag_sigma <- diag(sigma_hat)
  # theta = DELTA^(-2) - diag( LAMBDA (I-B)^-1 PSI (I-B)^-T t(LAMBDA) )
  theta_diag <- target_all - diag_sigma
  not_idx <- unique(c(num_idx, ov_y_dummy_ov_idx))
  if (length(not_idx) > 0L) {
    diag(mlist$theta)[-not_idx] <- theta_diag[-not_idx]
  } else {
    diag(mlist$theta) <- theta_diag
  }

  mlist
}

# if THETA parameterization, compute delta elements
# of observed categorical variables, as a function of other model parameters
lav_lisrel_delta <- function(mlist = NULL, num_idx = NULL) {
  sigma_hat <- lav_lisrel_sigma(mlist = mlist, delta = FALSE)
  diag_sigma <- diag(sigma_hat)

  # (1/delta^2) = diag( LAMBDA (I-B)^-1 PSI (I-B)^-T t(LAMBDA) ) + THETA
  # tmp <- diag.Sigma + THETA
  tmp <- diag_sigma
  tmp[tmp < 0] <- as.numeric(NA)
  mlist$delta[, 1L] <- sqrt(1 / tmp)

  # numeric delta's stay 1.0
  if (length(num_idx) > 0L) {
    mlist$delta[num_idx] <- 1.0
  }

  mlist
}

# compute Sigma/ETA: variances/covariances of BOTH observed and latent variables
lav_lisrel_cov_both <- function(mlist = NULL, delta = TRUE) {
  mm_lambda <- mlist$lambda
  nvar <- nrow(mm_lambda)
  mm_psi <- mlist$psi
  nlat <- nrow(mm_psi)
  mm_theta <- mlist$theta
  mm_beta <- mlist$beta

  # 'extend' matrices
  lambda2 <- rbind(mm_lambda, diag(nlat))
  theta2 <- lav_matrix_bdiag(mm_theta, matrix(0, nlat, nlat))


  # beta?
  if (is.null(mm_beta)) {
    lambda__ib_inv <- lambda2
  } else {
    ib_inv <- lav_lisrel_ibinv(mlist = mlist)
    lambda__ib_inv <- lambda2 %*% ib_inv
  }

  # compute augment COV matrix
  cov_1 <- tcrossprod(lambda__ib_inv %*% mm_psi, lambda__ib_inv) + theta2

  # if delta, scale
  if (delta && !is.null(mlist$delta)) {
    mm_delta <- diag(mlist$delta[, 1L], nrow = nvar, ncol = nvar)
    cov_1[seq_len(nvar), seq_len(nvar)] <-
      mm_delta %*% cov_1[seq_len(nvar), seq_len(nvar)] %*% mm_delta
  }


  # if GAMMA, also x part
  mm_gamma <- mlist$gamma
  if (!is.null(mm_gamma)) {
    cov_x <- mlist$cov.x
    if (is.null(mm_beta)) {
      sx <- tcrossprod(mm_gamma %*% cov_x, mm_gamma)
    } else {
      ib_inv__gamma <- ib_inv %*% mm_gamma
      sx <- tcrossprod(ib_inv__gamma %*% cov_x, ib_inv__gamma)
    }
    cov_1[(nvar + 1):(nvar + nlat), (nvar + 1):(nvar + nlat)] <-
      cov_1[(nvar + 1):(nvar + nlat), (nvar + 1):(nvar + nlat)] + sx
  }

  cov_1
}


# derivative of the objective function
lav_lisrel_df_dmlist <- function(mlist = NULL, omega = NULL, omega_mu = NULL) {
  mm_lambda <- mlist$lambda
  mm_psi <- mlist$psi
  mm_beta <- mlist$beta
  mm_alpha <- mlist$alpha
  # mm_wmat <- mlist$wmat

  lambda_deriv <- NULL
  beta_deriv <- NULL
  theta_deriv <- NULL
  psi_deriv <- NULL
  nu_deriv <- NULL
  alpha_deriv <- NULL
  group_w_deriv <- NULL
  wmat_deriv <- NULL

  # beta?
  if (is.null(mm_beta)) {
    lambda__ib_inv <- mm_lambda
  } else {
    ib_inv <- lav_lisrel_ibinv(mlist = mlist)
    lambda__ib_inv <- mm_lambda %*% ib_inv
  }

  # meanstructure?
  meanstructure <- FALSE
  if (!is.null(omega_mu)) meanstructure <- TRUE

  # group weight?
  group_w_free <- FALSE
  if (!is.null(mlist$gw)) group_w_free <- TRUE

  # pre-compute some values
  t_lambda__ib_inv <- t(lambda__ib_inv)
  if (!is.null(mm_beta)) {
    omega__lambda__ib_inv__psi__t_ib_inv <-                 # nolint
      (omega %*% lambda__ib_inv %*% mm_psi %*% t(ib_inv))
  } else {
    omega__lambda <- omega %*% mm_lambda
  }

  # 1. LAMBDA
  if (!is.null(mm_beta)) {
    if (meanstructure) {
      lambda_deriv <- -1.0 * (omega_mu %*% t(mm_alpha) %*% t(ib_inv) +
        omega__lambda__ib_inv__psi__t_ib_inv)
    } else {
      lambda_deriv <- -1.0 * omega__lambda__ib_inv__psi__t_ib_inv
    }
  } else {
    # no BETA
    if (meanstructure) {
      lambda_deriv <- -1.0 * (omega_mu %*% t(mm_alpha) +
        omega__lambda %*% mm_psi)
    } else {
      lambda_deriv <- -1.0 * (omega__lambda %*% mm_psi)
    }
  }

  # 2. BETA
  if (!is.null(mm_beta)) {
    if (meanstructure) {
      beta_deriv <- -1.0 * ((t(ib_inv) %*%
        (t(mm_lambda) %*% omega_mu %*% t(mm_alpha)) %*%
        t(ib_inv)) +
        (t_lambda__ib_inv %*%
          omega__lambda__ib_inv__psi__t_ib_inv))
    } else {
      beta_deriv <- -1.0 * (t_lambda__ib_inv %*%
        omega__lambda__ib_inv__psi__t_ib_inv)
    }
  }

  # 3. PSI
  psi_deriv <- -1.0 * (t_lambda__ib_inv %*% omega %*% lambda__ib_inv)
  diag(psi_deriv) <- 0.5 * diag(psi_deriv)

  # 4. THETA
  theta_deriv <- -1.0 * omega
  diag(theta_deriv) <- 0.5 * diag(theta_deriv)

  if (meanstructure) {
    # 5. NU
    nu_deriv <- -1.0 * omega_mu

    # 6. ALPHA
    alpha_deriv <- -1.0 * t(t(omega_mu) %*% lambda__ib_inv)
  }

  if (group_w_free) {
    group_w_deriv <- 0.0
  }

  list(
    lambda = lambda_deriv,
    wmat = wmat_deriv,
    beta = beta_deriv,
    theta = theta_deriv,
    psi = psi_deriv,
    nu = nu_deriv,
    alpha = alpha_deriv,
    gw = group_w_deriv
  )
}

# dSigma/dx -- per model matrix
lav_lisrel_dsigma_dx <- function(mlist = NULL,
                                 m = "lambda",
                                 # all model matrix elements, or only a few?
                                 # NOTE: for symmetric matrices,
                                 # we assume that the have full size
                                 # (nvar*nvar) (but already correct for
                                 # symmetry)
                                 idx = seq_along(mlist[[m]]),
                                 vech = TRUE,
                                 delta = TRUE) {
  mm_lambda <- mlist$lambda
  nvar <- nrow(mm_lambda)
  nfac <- ncol(mm_lambda)
  mm_psi <- mlist$psi
  mm_wmat <- mlist$wmat

  # for composites (vec version)
  compute_sigma <- function(x, mm = "wmat", mlist = NULL) {
    mlist_1 <- mlist
    if (mm %in% c("psi", "theta")) {
      mlist_1[[mm]] <- lav_matrix_vech_reverse(x)
    } else {
      mlist_1[[mm]][, ] <- x
    }
    lav_matrix_vec(lav_lisrel_sigma(mlist_1))
  }

  composites <- FALSE
  if (!is.null(mm_wmat)) {
    composites <- TRUE
  }

  # only lower.tri part of sigma (not same order as elimination matrix?)
  v_idx <- lav_matrix_vech_idx(nvar)
  pstar <- nvar * (nvar + 1) / 2

  # shortcut for gamma, nu, alpha, tau,.... : empty matrix
  if (m == "nu" || m == "alpha" || m == "tau" || m == "gamma" ||
    m == "gw" || m == "cov.x" || m == "mean.x") {
    return(matrix(0.0, nrow = pstar, ncol = length(idx)))
  }

  # Delta?
  delta_flag <- FALSE
  if (delta && !is.null(mlist$delta)) {
    mm_delta <- mlist$delta
    delta_flag <- TRUE
  } else if (m == "delta") { # modindices?
    return(matrix(0.0, nrow = pstar, ncol = length(idx)))
  }

  # beta?
  if (!is.null(mlist$ibeta.inv)) {
    ib_inv <- mlist$ibeta.inv
  } else {
    ib_inv <- lav_lisrel_ibinv(mlist = mlist)
  }

  # pre
  # if(m == "lambda" || m == "beta")
  #    IK <- diag(nvar*nvar) + lav_matrix_commutation(nvar, nvar)
  if (m == "lambda" || m == "beta") {
    l1 <- mm_lambda %*% ib_inv %*% mm_psi %*% t(ib_inv)
  }
  if (m == "beta" || m == "psi") {
    lambda__ib_inv <- mm_lambda %*% ib_inv
  }

  # here we go:
  if (m == "lambda") {
    kol_idx <- matrix(1:(nvar * nfac), nvar, nfac, byrow = TRUE)[idx]
    dx <- (l1 %x% diag(nvar))[, idx, drop = FALSE] +
      (diag(nvar) %x% l1)[, kol_idx, drop = FALSE]
  } else if (m == "beta") {
    if (composites) {
      dx <- lav_func_jacobian_complex(
        func = compute_sigma,
        x = lav_matrix_vec(mlist$beta),
        mm = "beta", mlist = mlist
      )
      dx <- dx[, idx, drop = FALSE]
    } else {
      kol_idx <- matrix(1:(nfac * nfac), nfac, nfac, byrow = TRUE)[idx]
      dx <- (l1 %x% lambda__ib_inv)[, idx, drop = FALSE] +
        (lambda__ib_inv %x% l1)[, kol_idx, drop = FALSE]
      # this is not really needed (because we select idx=m.el.idx)
      # but just in case we need all elements of beta...
      dx[, which(idx %in% lav_matrix_diag_idx(nfac))] <- 0.0
    }
  } else if (m == "psi") {
    if (composites) {
      tmp <- lav_func_jacobian_complex(
        func = compute_sigma,
        x = lav_matrix_vech(mlist$psi),
        mm = "psi", mlist = mlist
      )
      dx <- matrix(0, nrow = nrow(tmp), ncol = length(mm_psi))
      dx[, lav_matrix_vech_idx(nrow(mm_psi))] <- tmp
      dx[, lav_matrix_vechu_idx(nrow(mm_psi), diagonal = FALSE)] <-
        dx[, lav_matrix_vech_idx(nrow(mm_psi), diagonal = FALSE), drop = FALSE]
      dx <- dx[, idx, drop = FALSE]
    } else {
      dx <- (lambda__ib_inv %x% lambda__ib_inv)
      # symmetry correction, but keeping all duplicated elements
      # since we depend on idx=m.el.idx
      lower_idx <- lav_matrix_vech_idx(nfac, diagonal = FALSE)
      upper_idx <- lav_matrix_vechru_idx(nfac, diagonal = FALSE)
      offdiag_sum <- dx[, lower_idx] + dx[, upper_idx]
      dx[, c(lower_idx, upper_idx)] <- cbind(offdiag_sum, offdiag_sum)
      dx <- dx[, idx, drop = FALSE]
    }
  } else if (m == "theta") {
    # DX <- diag(nvar*nvar) # very sparse...
    dx <- matrix(0, nvar * nvar, length(idx))
    dx[cbind(idx, seq_along(idx))] <- 1
    # symmetry correction not needed, since all off-diagonal elements
    # are zero?
  } else if (m == "delta") {
    omega <- lav_lisrel_sigma(mlist, delta = FALSE)
    dd <- diag(mm_delta[, 1], nvar, nvar)
    dd_omega <- (dd %*% omega)
    a <- dd_omega %x% diag(nvar)
    m_b <- diag(nvar) %x% dd_omega
    dx <- a[, lav_matrix_diag_idx(nvar), drop = FALSE] +
      m_b[, lav_matrix_diag_idx(nvar), drop = FALSE]
    dx <- dx[, idx, drop = FALSE]
  } else if (m == "wmat") {
    # just a dummy to get us going
    dx <- lav_func_jacobian_complex(
      func = compute_sigma,
      x = lav_matrix_vec(mm_wmat),
      mm = "wmat", mlist = mlist
    )
    dx <- dx[, idx, drop = FALSE]

    # KOL.idx <- matrix(1:(nvar * nfac), nvar, nfac, byrow = TRUE)[idx]
    # VETA <- IB.inv %*% PSI %*% t(IB.inv)
    # C0 <- VETA; diag(C0) <- 0
    # cov.idx <- which(apply(LAMBDA, 1L,
    #                        function(x) sum(x == 0) == ncol(LAMBDA)))
    # clv.idx <- which(apply(LAMBDA, 2L,
    #                        function(x) sum(x == 0) == nrow(LAMBDA)))
    # Tmat <- diag(nrow(LAMBDA))
    # Tmat[cov.idx, cov.idx] <- MLIST$theta[cov.idx, cov.idx]
    # L1 <- Tmat %*% WMAT %*% C0
    # DX <- (L1 %x% diag(nvar))[, idx, drop = FALSE] +
    #   (diag(nvar) %x% L1)[, KOL.idx, drop = FALSE]
    # DX <- DX * nfac
  } else {
    lav_msg_stop(gettext("wrong model matrix name:"), m)
  }

  if (delta_flag && !m == "delta") {
    dx <- dx * as.vector(mm_delta %x% mm_delta)
  }

  # vech?
  if (vech) {
    dx <- dx[v_idx, , drop = FALSE]
  }

  dx
}

# dMu/dx -- per model matrix
lav_lisrel_dmu_dx <- function(mlist = NULL,
                              m = "alpha",
                              # all model matrix elements, or only a few?
                              idx = seq_along(mlist[[m]])) {
  mm_lambda <- mlist$lambda
  nvar <- nrow(mm_lambda)
  nfac <- ncol(mm_lambda)
  # mm_wmat <- mlist$wmat

  # shortcut for empty matrices
  if (m == "gamma" || m == "psi" || m == "theta" || m == "tau" ||
    m == "delta" || m == "gw" || m == "cov.x" || m == "mean.x") {
    return(matrix(0.0, nrow = nvar, ncol = length(idx)))
  }

  # missing alpha
  if (is.null(mlist$alpha)) {
    mm_alpha <- matrix(0, nfac, 1L)
  } else {
    mm_alpha <- mlist$alpha
  }


  # beta?
  if (!is.null(mlist$ibeta.inv)) {
    ib_inv <- mlist$ibeta.inv
  } else {
    ib_inv <- lav_lisrel_ibinv(mlist = mlist)
  }

  if (m == "nu") {
    dx <- diag(nvar)
  } else if (m == "lambda") {
    dx <- t(ib_inv %*% mm_alpha) %x% diag(nvar)
  } else if (m == "wmat") {
    # dummy, just to get us going
    dx <- t(ib_inv %*% mm_alpha) %x% diag(nvar)
  } else if (m == "beta") {
    dx <- t(ib_inv %*% mm_alpha) %x% (mm_lambda %*% ib_inv)
    # this is not really needed (because we select idx=m.el.idx)
    dx[, lav_matrix_diag_idx(nfac)] <- 0.0
  } else if (m == "alpha") {
    dx <- mm_lambda %*% ib_inv
  } else {
    lav_msg_stop(gettext("wrong model matrix names:"), m)
  }

  dx <- dx[, idx, drop = FALSE]
  dx
}

# dTh/dx -- per model matrix
lav_lisrel_dth_dx <- function(mlist = NULL,
                              m = "tau",
                              # all model matrix elements, or only a few?
                              idx = seq_along(mlist[[m]]),
                              th_idx = NULL,
                              delta = TRUE) {
  mm_lambda <- mlist$lambda
  nvar <- nrow(mm_lambda)
  nfac <- ncol(mm_lambda)
  mm_tau <- mlist$tau
  nth <- nrow(mm_tau)

  # missing alpha
  if (is.null(mlist$alpha)) {
    mm_alpha <- matrix(0, nfac, 1L)
  } else {
    mm_alpha <- mlist$alpha
  }

  # missing nu
  if (is.null(mlist$nu)) {
    mm_nu <- matrix(0, nvar, 1L)
  } else {
    mm_nu <- mlist$nu
  }

  # Delta?
  delta_flag <- FALSE
  if (delta && !is.null(mlist$delta)) {
    mm_delta <- mlist$delta
    delta_flag <- TRUE
  }

  if (is.null(th_idx)) {
    th_idx <- seq_len(nth)
    nlev <- rep(1L, nvar)
    k_nu <- diag(nvar)
  } else {
    nlev <- tabulate(th_idx, nbins = nvar)
    nlev[nlev == 0L] <- 1L
    k_nu <- matrix(0, sum(nlev), nvar)
    k_nu[cbind(seq_len(sum(nlev)), rep(seq_len(nvar), times = nlev))] <- 1.0
  }

  # shortcut for empty matrices
  if (m == "gamma" || m == "psi" || m == "theta" || m == "gw" ||
    m == "cov.x" || m == "mean.x") {
    return(matrix(0.0, nrow = length(th_idx), ncol = length(idx)))
  }

  # beta?
  if (!is.null(mlist$ibeta.inv)) {
    ib_inv <- mlist$ibeta.inv
  } else {
    ib_inv <- lav_lisrel_ibinv(mlist = mlist)
  }

  if (m == "tau") {
    dx <- matrix(0, nrow = length(th_idx), ncol = nth)
    dx[th_idx > 0L, ] <- diag(nth)
    if (delta_flag) {
      dx <- dx * as.vector(k_nu %*% mm_delta)
    }
  } else if (m == "nu") {
    dx <- (-1) * k_nu
    if (delta_flag) {
      dx <- dx * as.vector(k_nu %*% mm_delta)
    }
  } else if (m == "lambda") {
    dx <- (-1) * t(ib_inv %*% mm_alpha) %x% diag(nvar)
    dx <- k_nu %*% dx
    if (delta_flag) {
      dx <- dx * as.vector(k_nu %*% mm_delta)
    }
  } else if (m == "beta") {
    dx <- (-1) * t(ib_inv %*% mm_alpha) %x% (mm_lambda %*% ib_inv)
    # this is not really needed (because we select idx=m.el.idx)
    dx[, lav_matrix_diag_idx(nfac)] <- 0.0
    dx <- k_nu %*% dx
    if (delta_flag) {
      dx <- dx * as.vector(k_nu %*% mm_delta)
    }
  } else if (m == "alpha") {
    dx <- (-1) * mm_lambda %*% ib_inv
    dx <- k_nu %*% dx
    if (delta_flag) {
      dx <- dx * as.vector(k_nu %*% mm_delta)
    }
  } else if (m == "delta") {
    dx1 <- matrix(0, nrow = length(th_idx), ncol = 1)
    dx1[th_idx > 0L, ] <- mm_tau
    dx2 <- mm_nu + mm_lambda %*% ib_inv %*% mm_alpha
    dx2 <- k_nu %*% dx2
    dx <- k_nu * as.vector(dx1 - dx2)
  } else {
    lav_msg_stop(gettext("wrong model matrix names:"), m)
  }

  dx <- dx[, idx, drop = FALSE]
  dx
}

# dPi/dx -- per model matrix
lav_lisrel_dpi_dx <- function(mlist = NULL,
                              m = "lambda",
                              # all model matrix elements, or only a few?
                              idx = seq_along(mlist[[m]])) {
  mm_lambda <- mlist$lambda
  nvar <- nrow(mm_lambda)
  nfac <- ncol(mm_lambda)
  mm_gamma <- mlist$gamma
  nexo <- ncol(mm_gamma)

  # Delta?
  delta_flag <- FALSE
  if (!is.null(mlist$delta)) {
    delta_diag <- mlist$delta[, 1L]
    delta_flag <- TRUE
  }

  # shortcut for empty matrices
  if (m == "tau" || m == "nu" || m == "alpha" || m == "psi" ||
    m == "theta" || m == "gw" || m == "cov.x" || m == "mean.x") {
    return(matrix(0.0, nrow = nvar * nexo, ncol = length(idx)))
  }

  # beta?
  if (!is.null(mlist$ibeta.inv)) {
    ib_inv <- mlist$ibeta.inv
  } else {
    ib_inv <- lav_lisrel_ibinv(mlist = mlist)
  }

  if (m == "lambda") {
    dx <- t(ib_inv %*% mm_gamma) %x% diag(nvar)
    if (delta_flag) {
      dx <- dx * delta_diag
    }
  } else if (m == "beta") {
    dx <- t(ib_inv %*% mm_gamma) %x% (mm_lambda %*% ib_inv)
    # this is not really needed (because we select idx=m.el.idx)
    dx[, lav_matrix_diag_idx(nfac)] <- 0.0
    if (delta_flag) {
      dx <- dx * delta_diag
    }
  } else if (m == "gamma") {
    dx <- diag(nexo) %x% (mm_lambda %*% ib_inv)
    if (delta_flag) {
      dx <- dx * delta_diag
    }
  } else if (m == "delta") {
    pre <- rep(1, nexo) %x% diag(nvar)
    dx <- pre * as.vector(mm_lambda %*% ib_inv %*% mm_gamma)
  } else {
    lav_msg_stop(gettext("wrong model matrix names:"), m)
  }

  dx <- dx[, idx, drop = FALSE]
  dx
}

# dGW/dx -- per model matrix
lav_lisrel_dgw_dx <- function(mlist = NULL,
                              m = "gw",
                              # all model matrix elements, or only a few?
                              idx = seq_along(mlist[[m]])) {
  # shortcut for empty matrices
  if (m != "gw") {
    return(matrix(0.0, nrow = 1L, ncol = length(idx)))
  } else {
    # m == "gw"
    dx <- matrix(1.0, 1, 1)
  }

  dx <- dx[, idx, drop = FALSE]
  dx
}

# dlambda/dx -- per model matrix
lav_lisrel_dlambda_dx <- function(mlist = NULL,
                                  m = "lambda",
                                  # all model matrix elements, or only a few?
                                  idx = seq_along(mlist[[m]])) {
  mm_lambda <- mlist$lambda

  # shortcut for empty matrices
  if (m != "lambda") {
    return(matrix(0.0, nrow = length(mm_lambda), ncol = length(idx)))
  } else {
    # m == "lambda"
    dx <- diag(1, nrow = length(mm_lambda), ncol = length(mm_lambda))
  }

  dx <- dx[, idx, drop = FALSE]
  dx
}

# dpsi/dx -- per model matrix - FIXME!!!!!
lav_lisrel_dpsi_dx <- function(mlist = NULL,
                               m = "psi",
                               # all model matrix elements, or only a few?
                               idx = seq_along(mlist[[m]])) {
  mm_psi <- mlist$psi
  nfac <- nrow(mm_psi)
  v_idx <- lav_matrix_vech_idx(nfac)

  # shortcut for empty matrices
  if (m != "psi") {
    dx <- matrix(0.0, nrow = length(mm_psi), ncol = length(idx))
    return(dx[v_idx, , drop = FALSE])
  } else {
    # m == "psi"
    dx <- diag(1, nrow = length(mm_psi), ncol = length(mm_psi))
  }

  dx <- dx[v_idx, idx, drop = FALSE]
  dx
}

# dtheta/dx -- per model matrix
lav_lisrel_dtheta_dx <- function(mlist = NULL,
                                 m = "theta",
                                 # all model matrix elements, or only a few?
                                 idx = seq_along(mlist[[m]])) {
  mm_theta <- mlist$theta
  nvar <- nrow(mm_theta)
  v_idx <- lav_matrix_vech_idx(nvar)

  # shortcut for empty matrices
  if (m != "theta") {
    dx <- matrix(0.0, nrow = length(mm_theta), ncol = length(idx))
    return(dx[v_idx, , drop = FALSE])
  } else {
    # m == "theta"
    dx <- diag(1, nrow = length(mm_theta), ncol = length(mm_theta))
  }

  dx <- dx[v_idx, idx, drop = FALSE]
  dx
}


# dbeta/dx -- per model matrix
lav_lisrel_dbeta_dx <- function(mlist = NULL,
                                m = "beta",
                                # all model matrix elements, or only a few?
                                idx = seq_along(mlist[[m]])) {
  mm_beta <- mlist$beta

  # shortcut for empty matrices
  if (m != "beta") {
    return(matrix(0.0, nrow = length(mm_beta), ncol = length(idx)))
  } else {
    # m == "beta"
    dx <- diag(1, nrow = length(mm_beta), ncol = length(mm_beta))
  }

  dx <- dx[, idx, drop = FALSE]
  dx
}

# dgamma/dx -- per model matrix
lav_lisrel_dgamma_dx <- function(mlist = NULL,
                                 m = "gamma",
                                 # all model matrix elements, or only a few?
                                 idx = seq_along(mlist[[m]])) {
  mm_gamma <- mlist$gamma

  # shortcut for empty matrices
  if (m != "gamma") {
    return(matrix(0.0, nrow = length(mm_gamma), ncol = length(idx)))
  } else {
    # m == "gamma"
    dx <- diag(1, nrow = length(mm_gamma), ncol = length(mm_gamma))
  }

  dx <- dx[, idx, drop = FALSE]
  dx
}

# dnu/dx -- per model matrix
lav_lisrel_dnu_dx <- function(mlist = NULL,
                              m = "nu",
                              # all model matrix elements, or only a few?
                              idx = seq_along(mlist[[m]])) {
  mm_nu <- mlist$nu

  # shortcut for empty matrices
  if (m != "nu") {
    return(matrix(0.0, nrow = length(mm_nu), ncol = length(idx)))
  } else {
    # m == "nu"
    dx <- diag(1, nrow = length(mm_nu), ncol = length(mm_nu))
  }

  dx <- dx[, idx, drop = FALSE]
  dx
}

# dtau/dx -- per model matrix
lav_lisrel_dtau_dx <- function(mlist = NULL,
                               m = "tau",
                               # all model matrix elements, or only a few?
                               idx = seq_along(mlist[[m]])) {
  mm_tau <- mlist$tau

  # shortcut for empty matrices
  if (m != "tau") {
    return(matrix(0.0, nrow = length(mm_tau), ncol = length(idx)))
  } else {
    # m == "tau"
    dx <- diag(1, nrow = length(mm_tau), ncol = length(mm_tau))
  }

  dx <- dx[, idx, drop = FALSE]
  dx
}


# dalpha/dx -- per model matrix
lav_lisrel_dalpha_dx <- function(mlist = NULL,
                                 m = "alpha",
                                 # all model matrix elements, or only a few?
                                 idx = seq_along(mlist[[m]])) {
  mm_alpha <- mlist$alpha

  # shortcut for empty matrices
  if (m != "alpha") {
    return(matrix(0.0, nrow = length(mm_alpha), ncol = length(idx)))
  } else {
    # m == "alpha"
    dx <- diag(1, nrow = length(mm_alpha), ncol = length(mm_alpha))
  }

  dx <- dx[, idx, drop = FALSE]
  dx
}

# MLIST = NULL; meanstructure=TRUE; th=TRUE; delta=TRUE; pi=TRUE; gw=FALSE
# lav_matrix_vech_idx <- lavaan:::lav_matrix_vech_idx;
# lav_matrix_vechru_idx <- lavaan:::lav_matrix_vechru_idx
# vec <- lavaan:::vec;
# lav_func_jacobian_complex <- lavaan:::lav_func_jacobian_complex
# lav_lisrel_sigma <- lavaan:::lav_lisrel_sigma
# lav_lisrel_delta <- lavaan:::lav_lisrel_delta
lav_lisrel_test_derivatives <- function(mlist = NULL,
                                        nvar = NULL, nfac = NULL, nexo = NULL,
                                        th_idx = NULL, num_idx = NULL,
                                        meanstructure = TRUE,
                                        th = TRUE, delta = TRUE, pi = TRUE,
                                        gw = FALSE, theta = FALSE) {
  if (is.null(mlist)) {
    # create artificial matrices, compare 'numerical' vs 'analytical'
    # derivatives
    # nvar <- 12; nfac <- 3; nexo <- 4 # this combination is special?
    if (is.null(nvar)) {
      nvar <- 20
    }
    if (is.null(nfac)) {
      nfac <- 6
    }
    if (is.null(nexo)) {
      nexo <- 5
    }
    if (is.null(num_idx)) {
      num_idx <- sort(sample(seq_len(nvar), ceiling(nvar / 2)))
    }
    if (is.null(th_idx)) {
      th_idx <- integer(0L)
      for (i in seq_len(nvar)) {
        if (i %in% num_idx) {
          th_idx <- c(th_idx, 0)
        } else {
          th_idx <- c(th_idx, rep(i, sample(c(1, 1, 2, 6), 1L)))
        }
      }
    }
    nth <- sum(th_idx > 0L)

    mlist <- list()
    mlist$lambda <- matrix(0, nvar, nfac)
    mlist$beta <- matrix(0, nfac, nfac)
    mlist$theta <- matrix(0, nvar, nvar)
    mlist$psi <- matrix(0, nfac, nfac)
    if (meanstructure) {
      mlist$alpha <- matrix(0, nfac, 1L)
      mlist$nu <- matrix(0, nvar, 1L)
    }
    if (th) mlist$tau <- matrix(0, nth, 1L)
    if (delta) mlist$delta <- matrix(0, nvar, 1L)
    mlist$gamma <- matrix(0, nfac, nexo)
    if (gw) mlist$gw <- matrix(0, 1L, 1L)

    # feed random numbers
    mlist <- lapply(mlist, function(x) {
      x[, ] <- rnorm(length(x))
      x
    })
    # fix
    diag(mlist$beta) <- 0.0
    diag(mlist$theta) <- diag(mlist$theta) * diag(mlist$theta) * 10
    diag(mlist$psi) <- diag(mlist$psi) * diag(mlist$psi) * 10
    mlist$psi[lav_matrix_vechru_idx(nfac)] <-
      mlist$psi[lav_matrix_vech_idx(nfac)]
    mlist$theta[lav_matrix_vechru_idx(nvar)] <-
      mlist$theta[lav_matrix_vech_idx(nvar)]
    if (delta) mlist$delta[, ] <- abs(mlist$delta) * 10
  } else {
    nvar <- nrow(mlist$lambda)
  }

  compute_sigma <- function(x, mm = "lambda", mlist = NULL) {
    mlist_1 <- mlist
    if (mm %in% c("psi", "theta")) {
      mlist_1[[mm]] <- lav_matrix_vech_reverse(x)
    } else {
      mlist_1[[mm]][, ] <- x
    }
    if (theta) {
      mlist_1 <- lav_lisrel_delta(mlist = mlist_1, num_idx = num_idx)
    }
    lav_matrix_vech(lav_lisrel_sigma(mlist_1))
  }

  compute_mu <- function(x, mm = "lambda", mlist = NULL) {
    mlist_1 <- mlist
    if (mm %in% c("psi", "theta")) {
      mlist_1[[mm]] <- lav_matrix_vech_reverse(x)
    } else {
      mlist_1[[mm]][, ] <- x
    }
    if (theta) {
      mlist_1 <- lav_lisrel_delta(mlist = mlist_1, num_idx = num_idx)
    }
    lav_lisrel_mu(mlist_1)
  }

  compute_th2 <- function(x, mm = "tau", mlist = NULL, th_idx) {
    mlist_1 <- mlist
    if (mm %in% c("psi", "theta")) {
      mlist_1[[mm]] <- lav_matrix_vech_reverse(x)
    } else {
      mlist_1[[mm]][, ] <- x
    }
    if (theta) {
      mlist_1 <- lav_lisrel_delta(mlist = mlist_1, num_idx = num_idx)
    }
    lav_lisrel_th(mlist_1, th_idx = th_idx)
  }

  compute_pi <- function(x, mm = "lambda", mlist = NULL) {
    mlist_1 <- mlist
    if (mm %in% c("psi", "theta")) {
      mlist_1[[mm]] <- lav_matrix_vech_reverse(x)
    } else {
      mlist_1[[mm]][, ] <- x
    }
    if (theta) {
      mlist_1 <- lav_lisrel_delta(mlist = mlist_1, num_idx = num_idx)
    }
    lav_lisrel_pi(mlist_1)
  }

  compute_gw <- function(x, mm = "gw", mlist = NULL) {
    mlist_1 <- mlist
    if (mm %in% c("psi", "theta")) {
      mlist_1[[mm]] <- lav_matrix_vech_reverse(x)
    } else {
      mlist_1[[mm]][, ] <- x
    }
    if (theta) {
      mlist_1 <- lav_lisrel_delta(mlist = mlist_1, num_idx = num_idx)
    }
    mlist_1$gw[1, 1]
  }

  # if theta, set MLIST$delta
  if (theta) {
    mlist <- lav_lisrel_delta(mlist = mlist, num_idx = num_idx)
  }

  for (mm in names(mlist)) {
    if (mm %in% c("psi", "theta")) {
      x <- lav_matrix_vech(mlist[[mm]])
    } else {
      x <- lav_matrix_vec(mlist[[mm]])
    }
    if (mm == "delta" && theta) next
    if (lav_debug()) {
      cat("### mm = ", mm, "\n")
    }

    # 1. sigma
    dx1 <- lav_func_jacobian_complex(func = compute_sigma,
                                      x = x, mm = mm, mlist = mlist)
    dx2 <- lav_lisrel_dsigma_dx(
      mlist = mlist, m = mm, idx = seq_along(mlist[[mm]]),
      delta = !theta
    )
    if (mm %in% c("psi", "theta")) {
      # remove duplicated columns of symmetric matrices
      idx <- lav_matrix_vechru_idx(sqrt(ncol(dx2)), diagonal = FALSE)
      if (length(idx) > 0L) dx2 <- dx2[, -idx]
    }
    if (theta) {
      sigma_hat <- lav_lisrel_sigma(mlist = mlist, delta = FALSE)
      r <- lav_deriv_cov2cor(sigma_hat, num.idx = num_idx)

      dx3 <- dx2
      dx2 <- r %*% dx2
    }
    if (lav_debug()) {
      cat("[SIGMA] mm = ", sprintf("%-8s:", mm), "DX1 (numerical):\n")
      print(zapsmall(dx1))
      cat("\n")
      cat("[SIGMA] mm = ", sprintf("%-8s:", mm), "DX2 (analytical):\n")
      print(dx2)
      cat("\n")
      if (theta) {
        cat("[SIGMA] mm = ", sprintf("%-8s:", mm), "DX3 (analytical):\n")
        print(dx3)
        cat("\n")
      }
    }
    cat(
      "[SIGMA] mm = ", sprintf("%-8s:", mm), "sum delta = ",
      sprintf("%12.9f", sum(dx1 - dx2)), "  max delta = ",
      sprintf("%12.9f", max(dx1 - dx2)), "\n"
    )

    # 2. mu
    dx1 <- lav_func_jacobian_complex(func = compute_mu,
                                    x = x, mm = mm, mlist = mlist)
    dx2 <- lav_lisrel_dmu_dx(
      mlist = mlist,
      m = mm, idx = seq_along(mlist[[mm]])
    )
    if (mm %in% c("psi", "theta")) {
      # remove duplicated columns of symmetric matrices
      idx <- lav_matrix_vechru_idx(sqrt(ncol(dx2)), diagonal = FALSE)
      if (length(idx) > 0L) dx2 <- dx2[, -idx]
    }
    cat(
      "[MU   ] mm = ", sprintf("%-8s:", mm), "sum delta = ",
      sprintf("%12.9f", sum(dx1 - dx2)), "  max delta = ",
      sprintf("%12.9f", max(dx1 - dx2)), "\n"
    )
    if (lav_debug()) {
      cat("[MU   ] mm = ", sprintf("%-8s:", mm), "DX1 (numerical):\n")
      print(zapsmall(dx1))
      cat("\n")
      cat("[MU   ] mm = ", sprintf("%-8s:", mm), "DX2 (analytical):\n")
      print(dx2)
      cat("\n")
    }

    # 3. th
    if (th) {
      dx1 <- lav_func_jacobian_complex(
        func = compute_th2, x = x, mm = mm, mlist = mlist, th_idx = th_idx
      )
      dx2 <- lav_lisrel_dth_dx(
        mlist = mlist, m = mm, idx = seq_along(mlist[[mm]]),
        th_idx = th_idx,
        delta = TRUE
      )
      if (theta) {
        # 1. compute dDelta.dx
        dx_sigma <-
          lav_lisrel_dsigma_dx(
            m = mm, idx = seq_along(mlist[[mm]]),
            mlist = mlist, delta = !theta
          )
        var_idx <- which(!lav_matrix_vech_idx(nvar) %in%
          lav_matrix_vech_idx(nvar, diagonal = FALSE))
        sigma_hat <- lav_lisrel_sigma(mlist = mlist, delta = FALSE)
        dsigma <- diag(sigma_hat)
        # dy/ddsigma = -0.5/(ddsigma*sqrt(ddsigma))
        d_delta_dx <- dx_sigma[var_idx, ] * -0.5 / (dsigma * sqrt(dsigma))

        # 2. compute dth.dDelta
        dth_d_delta <-
          lav_lisrel_dth_dx(
            mlist = mlist,
            m = "delta",
            idx = seq_along(mlist[["delta"]]),
            th_idx = th_idx
          )

        # 3. add dth.dDelta %*% dDelta.dx
        no_num_idx <- which(th_idx > 0)
        dx2[no_num_idx, ] <- dx2[no_num_idx, , drop = FALSE] +
          (dth_d_delta %*% d_delta_dx)[no_num_idx, , drop = FALSE]
        # DX2 <- DX2 + dth.dDelta %*% dDelta.dx
      }
      if (mm %in% c("psi", "theta")) {
        # remove duplicated columns of symmetric matrices
        idx <- lav_matrix_vechru_idx(sqrt(ncol(dx2)), diagonal = FALSE)
        if (length(idx) > 0L) dx2 <- dx2[, -idx]
      }
      cat(
        "[TH   ] mm = ", sprintf("%-8s:", mm), "sum delta = ",
        sprintf("%12.9f", sum(dx1 - dx2)), "  max delta = ",
        sprintf("%12.9f", max(dx1 - dx2)), "\n"
      )
      if (lav_debug()) {
        cat("[TH   ] mm = ", sprintf("%-8s:", mm), "DX1 (numerical):\n")
        print(zapsmall(dx1))
        cat("\n")
        cat("[TH   ] mm = ", sprintf("%-8s:", mm), "DX2 (analytical):\n")
        print(dx2)
        cat("\n")
      }
    }

    # 4. pi
    if (pi) {
      dx1 <- lav_func_jacobian_complex(func = compute_pi,
                                         x = x, mm = mm, mlist = mlist)
      dx2 <- lav_lisrel_dpi_dx(
        mlist = mlist,
        m = mm, idx = seq_along(mlist[[mm]])
      )
      if (mm %in% c("psi", "theta")) {
        # remove duplicated columns of symmetric matrices
        idx <- lav_matrix_vechru_idx(sqrt(ncol(dx2)), diagonal = FALSE)
        if (length(idx) > 0L) dx2 <- dx2[, -idx]
      }
      if (theta) {
        # 1. compute dDelta.dx
        dx_sigma <-
          lav_lisrel_dsigma_dx(
            mlist = mlist, m = mm, idx = seq_along(mlist[[mm]]),
            delta = !theta
          )
        if (mm %in% c("psi", "theta")) {
          # remove duplicated columns of symmetric matrices
          idx <- lav_matrix_vechru_idx(sqrt(ncol(dx_sigma)), diagonal = FALSE)
          if (length(idx) > 0L) dx_sigma <- dx_sigma[, -idx]
        }
        var_idx <- which(!lav_matrix_vech_idx(nvar) %in%
          lav_matrix_vech_idx(nvar, diagonal = FALSE))
        sigma_hat <- lav_lisrel_sigma(mlist = mlist, delta = FALSE)
        dsigma <- diag(sigma_hat)
        # dy/ddsigma = -0.5/(ddsigma*sqrt(ddsigma))
        d_delta_dx <- dx_sigma[var_idx, ] * -0.5 / (dsigma * sqrt(dsigma))

        # 2. compute dpi.dDelta
        dpi_d_delta <-
          lav_lisrel_dpi_dx(
            mlist = mlist,
            m = "delta",
            idx = seq_along(mlist[["delta"]])
          )

        # 3. add dpi.dDelta %*% dDelta.dx
        no_num_idx <- which(!seq.int(1L, nvar) %in% num_idx)
        no_num_idx <- rep(seq.int(0, nexo - 1) * nvar,
          each = length(no_num_idx)
        ) + no_num_idx
        dx2[no_num_idx, ] <- dx2[no_num_idx, , drop = FALSE] +
          (dpi_d_delta %*% d_delta_dx)[no_num_idx, , drop = FALSE]
      }
      cat(
        "[PI   ] mm = ", sprintf("%-8s:", mm), "sum delta = ",
        sprintf("%12.9f", sum(dx1 - dx2)), "  max delta = ",
        sprintf("%12.9f", max(dx1 - dx2)), "\n"
      )
      if (lav_debug()) {
        cat("[PI   ] mm = ", sprintf("%-8s:", mm), "DX1 (numerical):\n")
        print(zapsmall(dx1))
        cat("\n")
        cat("[PI   ] mm = ", sprintf("%-8s:", mm), "DX2 (analytical):\n")
        print(dx2)
        cat("\n")
      }
    }

    # 5. gw
    if (gw) {
      dx1 <- lav_func_jacobian_complex(func = compute_gw,
                                          x = x, mm = mm, mlist = mlist)
      dx2 <- lav_lisrel_dgw_dx(
        mlist = mlist,
        m = mm, idx = seq_along(mlist[[mm]])
      )
      if (mm %in% c("psi", "theta")) {
        # remove duplicated columns of symmetric matrices
        idx <- lav_matrix_vechru_idx(sqrt(ncol(dx2)), diagonal = FALSE)
        if (length(idx) > 0L) dx2 <- dx2[, -idx]
      }
      cat(
        "[GW   ] mm = ", sprintf("%-8s:", mm), "sum delta = ",
        sprintf("%12.9f", sum(dx1 - dx2)), "  max delta = ",
        sprintf("%12.9f", max(dx1 - dx2)), "\n"
      )
      if (lav_debug()) {
        cat("[GW   ] mm = ", sprintf("%-8s:", mm), "DX1 (numerical):\n")
        print(dx1)
        cat("\n\n")
        cat("[GW   ] mm = ", sprintf("%-8s:", mm), "DX2 (analytical):\n")
        print(dx2)
        cat("\n\n")
      }
    }
  }

  mlist$th.idx <- th_idx
  mlist$num.idx <- num_idx

  mlist
}

# check for marker indicators:
#   - if std.lv = FALSE: a single '1' per factor, everything else zero
#   - if std.lv = TRUE: a single non-zero value per factor, everything else zero
lav_utils_get_marker <- function(mm_lambda = NULL, std_lv = FALSE) {
  mm_lambda <- as.matrix(mm_lambda)
  # nvar <- nrow(mm_lambda)
  nfac <- ncol(mm_lambda)

  # round values
  mm_lambda <- round(mm_lambda, 3L)

  marker_idx <- numeric(nfac)
  for (f in seq_len(nfac)) {
    if (std_lv) {
      marker_idx[f] <- which(rowSums(cbind(
        mm_lambda[, f] != 0,
        mm_lambda[, -f] == 0
      )) == nfac)[1]
    } else {
      marker_idx[f] <- which(rowSums(cbind(
        mm_lambda[, f] == 1,
        mm_lambda[, -f] == 0
      )) == nfac)[1]
    }
  }

  marker_idx
}

# find the residual variances (diagonal element of PSI) in such a way
# so that the diagonal elements of IB.inv %*% PSI %*% t(IB.inv) are
# equal to the elements of target.psi (usually the 1 vector)
#
# YR 01 Nov 2024: initial version; no bounds for now... (so we may end up
#                 with negative variances)
lav_mlist_target_psi <- function(mm_beta = NULL, ib_inv = NULL, mm_psi = NULL,
                                 target_psi = NULL, y_idx = NULL) {
  nr <- nrow(mm_psi)

  # IB.inv (if not given)
  if (is.null(ib_inv)) {
    ib <- -mm_beta
    ib[lav_matrix_diag_idx(nr)] <- 1
    ib_inv <- solve(ib)
  }

  # target.psi
  if (is.null(target_psi)) {
    target_psi <- rep(1, nr)
  }

  # y.idx
  if (is.null(y_idx)) {
    y_idx <- seq_len(nr)
  }

  # cast the problem as a nonlinear optimization problem
  obj <- function(x) {
    # cat("x = ", x, "\n")
    # x are the diagonal elements of PSI
    this_psi <- mm_psi
    diag(this_psi)[y_idx] <- x
    veta <- ib_inv %*% this_psi %*% t(ib_inv)
    current_diag <- diag(veta)
    # ratio or difference?
    diff <- target_psi[y_idx] - current_diag[y_idx]
    out <- sum(diff * diff) # least squares
    out
  }

  veta <- ib_inv %*% mm_psi %*% t(ib_inv)
  x_start <- diag(veta)[y_idx]
  out <- nlminb(start = x_start, objective = obj)

  # return updated PSI matrix
  diag(mm_psi)[y_idx] <- out$par
  mm_psi
}

# second-order derivatives of Sigma wrt elements of Lambda/Beta/Psi/Theta
# (elementwise)
#
# work in progress
lav_lisrel_d2sigma_lambda_psi <- function(mlist = NULL,
                                          i = 1L, j = 1L, # lambda
                                          k = 1L, l = 1L) { # psi
  # helper functions
  # elementary matrix: p x q matrix with 1 in position (i,j), 0 elsewhere
  elem_mat <- function(i, j, p, q) {
    m_e <- matrix(0, p, q)
    m_e[i, j] <- 1
    m_e
  }

  # elementary symmetric matrix for vech index (k,l) of a q x q symmetric matrix
  elem_mat_sym <- function(k, l, q) {
    m_e <- matrix(0, q, q)
    m_e[k, l] <- 1
    if (k != l) m_e[l, k] <- 1
    m_e
  }

  mm_lambda <- mlist$lambda
  nvar <- nrow(mm_lambda) # p
  nfac <- ncol(mm_lambda) # q
  ib_inv <- mlist$IB.inv
  if (is.null(ib_inv)) {
    if (is.null(mlist$beta)) {
      ib_inv <- diag(1, nrow = nfac)
    } else {
      ib_inv <- lav_lisrel_ibinv(mlist = mlist)
    }
  }

  lib_inv <- mm_lambda %*% ib_inv

  d_lambda <- elem_mat(i, j, nvar, nfac) # p x q,  elementary matrix for Lambda
  d_psi <- elem_mat_sym(k, l, nfac) # q x q,  elementary symm matrix for Psi

  x <- d_lambda %*% ib_inv %*% d_psi %*% t(lib_inv) # p x p
  x + t(x) # p x p symmetric matrix
}

lav_lisrel_d2sigma_beta_psi <- function(mlist = NULL,
                                        i = 1L, j = 1L, # beta
                                        k = 1L, l = 1L) { # psi
  # helper functions
  # elementary matrix: p x q matrix with 1 in position (i,j), 0 elsewhere
  elem_mat <- function(i, j, p, q) {
    m_e_e <- matrix(0, p, q)
    m_e_e[i, j] <- 1
    m_e_e
  }

  # elementary symmetric matrix for vech index (k,l) of a q x q symmetric matrix
  elem_mat_sym <- function(k, l, q) {
    m_e_e <- matrix(0, q, q)
    m_e_e[k, l] <- 1
    if (k != l) m_e_e[l, k] <- 1
    m_e_e
  }

  mm_lambda <- mlist$lambda
  # nvar <- nrow(mm_lambda) # p
  nfac <- ncol(mm_lambda) # q
  ib_inv <- mlist$IB.inv
  if (is.null(ib_inv)) {
    if (is.null(mlist$beta)) {
      ib_inv <- diag(1, nrow = nfac)
    } else {
      ib_inv <- lav_lisrel_ibinv(mlist = mlist)
    }
  }
  lib_inv <- mm_lambda %*% ib_inv

  d_beta <- elem_mat(i, j, nfac, nfac) # q x q, elementary matrix for Beta
  d_psi <- elem_mat_sym(k, l, nfac) # q x q, elementary symm matrix for Psi

  x <- lib_inv %*% d_beta %*% ib_inv %*% d_psi %*% t(lib_inv) # p x p
  x + t(x) # p x p symmetric matrix
}


# Jacobian of the model-implied moments wrt the free parameters
# (single block, LISREL representation).
#
# Inputs:
#   MLIST         : named list of model matrices for the block (lambda, psi,
#                   theta, [beta], [delta], [wmat], [nu], [alpha], [tau])
#   m.free.idx    : list of length(MLIST); m.free.idx[[mm]] holds the
#                   matrix-element indices of the free elements of MLIST[[mm]]
#   x.free.idx    : list of length(MLIST); x.free.idx[[mm]] holds the
#                   parameter-vector positions of those free elements
#   nx.free       : total number of free parameters (column dimension of out)
#   meanstructure : logical
#   categorical   : logical
#   correlation   : logical (continuous correlation structure)
#   num.idx       : integer vector of numeric-variable indices (block)
#   th.idx        : integer vector of threshold indices (block)
#   group.w.free  : logical; if TRUE, prepend a row for the free group weight
#
# Output: a matrix with nx.free columns, rows are the model-implied moments
# in the order used elsewhere in lavaan:
# [group weight (if group.w.free) | thresholds | means] then [vech(Sigma)].
lav_lisrel_dimplied_dx <- function(mlist           = NULL,
                                   m_free_idx      = NULL,
                                   x_free_idx      = NULL,
                                   nx_free         = NULL,
                                   meanstructure   = FALSE,
                                   categorical     = FALSE,
                                   correlation     = FALSE,
                                   conditional_x   = FALSE,
                                   num_idx         = integer(0L),
                                   th_idx          = integer(0L),
                                   group_w_free    = FALSE,
                                   parameterization = "delta") {

  delta_flag <- beta_flag <- FALSE
  if (!is.null(mlist$delta) && any(mlist$delta[, 1] != 1)) {
    delta_flag <- TRUE
  }
  if (!is.null(mlist$beta)) {
    beta_flag <- TRUE
  }

  # model matrices in this block
  mnames <- names(mlist)

  mm_lambda_idx <- which(mnames == "lambda")
  x_lambda_idx <- x_free_idx[[mm_lambda_idx]]
  m_lambda_idx <- m_free_idx[[mm_lambda_idx]]
  n_lam <- length(m_lambda_idx)

  mm_psi_idx <- which(mnames == "psi")
  tmp <- x_free_idx[[mm_psi_idx]]
  x_psi_idx <- x_psi_idx <- tmp[!duplicated(tmp)]
  m_psi_idx <- m_free_idx[[mm_psi_idx]][!duplicated(tmp)]
  n_psi <- length(m_psi_idx)


  mm_theta_idx <- which(mnames == "theta")
  tmp <- x_free_idx[[mm_theta_idx]]
  x_theta_idx <- tmp[!duplicated(tmp)]
  m_theta_idx <- m_free_idx[[mm_theta_idx]][!duplicated(tmp)]
  n_the <- length(m_theta_idx)

  n_bet <- 0L
  x_beta_idx <- integer(0L)
  if (beta_flag) {
    mm_beta_idx <- which(mnames == "beta")
    x_beta_idx <- x_free_idx[[mm_beta_idx]]
    m_beta_idx <- m_free_idx[[mm_beta_idx]]
    n_bet <- length(m_beta_idx)
  }

  n_del <- 0L
  x_delta_idx <- integer(0L)
  if (delta_flag) {
    mm_delta_idx <- which(mnames == "delta")
    x_delta_idx <- x_free_idx[[mm_delta_idx]]
    m_delta_idx <- m_free_idx[[mm_delta_idx]]
    n_del <- length(m_delta_idx)
  }

  # composite: wmat
  n_wmat <- 0L
  x_wmat_idx <- integer(0L)
  m_wmat_idx <- integer(0L)
  wmat_flag <- !is.null(mlist$wmat)
  if (wmat_flag) {
    mm_wmat_idx <- which(mnames == "wmat")
    if (length(mm_wmat_idx) > 0L) {
      x_wmat_idx <- x_free_idx[[mm_wmat_idx]]
      m_wmat_idx <- m_free_idx[[mm_wmat_idx]]
      n_wmat <- length(m_wmat_idx)
    }
  }

  # in case !meanstructure or !categorical
  n_nu  <- 0L
  m_nu_idx    <- integer(0L)
  x_nu_idx    <- integer(0L)
  n_alp <- 0L
  m_alpha_idx <- integer(0L)
  x_alpha_idx <- integer(0L)
  n_gam <- 0L
  m_gamma_idx <- integer(0L)
  x_gamma_idx <- integer(0L)


  if (meanstructure) {
    mm_nu_idx <- which(mnames == "nu")
    x_nu_idx <- x_free_idx[[mm_nu_idx]]
    m_nu_idx <- m_free_idx[[mm_nu_idx]]
    n_nu <- length(m_nu_idx)

    mm_alpha_idx <- which(mnames == "alpha")
    x_alpha_idx <- x_free_idx[[mm_alpha_idx]]
    m_alpha_idx <- m_free_idx[[mm_alpha_idx]]
    n_alp <- length(m_alpha_idx)
  }

  if (conditional_x) {
    mm_gamma_idx <- which(mnames == "gamma")
    if (length(mm_gamma_idx) > 0L) {
      x_gamma_idx <- x_free_idx[[mm_gamma_idx]]
      m_gamma_idx <- m_free_idx[[mm_gamma_idx]]
      n_gam <- length(m_gamma_idx)
    }
  }

  if (categorical) {
    mm_th_idx <- which(mnames == "tau")
    x_th_idx <- x_free_idx[[mm_th_idx]]
    m_th_idx <- m_free_idx[[mm_th_idx]]
    n_th <- length(m_th_idx)
  }

  nvar <- nrow(mlist$lambda)
  nfac <- ncol(mlist$lambda)
  pstar <- nvar * (nvar + 1L) / 2L

  # precompute
  if (beta_flag) {
    a_1 <- lav_lisrel_ibinv(mlist)
    m <- mlist$lambda %*% a_1
    g <- a_1 %*% mlist$psi %*% t(a_1)
    lg <- mlist$lambda %*% g
  } else {
    m <- mlist$lambda
    g <- mlist$psi
    lg <- mlist$lambda %*% mlist$psi
  }

  if (delta_flag) {
    delta  <- as.vector(mlist$delta)
  }


  # vech structure for Sigma
  r_s <- lav_matrix_vech_row_idx(nvar)
  c_s <- lav_matrix_vech_col_idx(nvar)
  sigma_lut <- lav_matrix_vech_reverse(seq_len(pstar))

  if (delta_flag) {
    # scaling weights: delta[r] * delta[s] for each vech position
    delta_weight <- delta[r_s] * delta[c_s]
  }

  # helper: vec index -> (row, col)
  vec2rc <- function(idx, nr) {
    cbind(row = (idx - 1L) %% nr + 1L,
          col = (idx - 1L) %/% nr + 1L)
  }

  # prepare output matrix for sigma
  n_free <- n_lam + n_wmat + n_bet + n_psi + n_the + n_del
  jac_sigma <- matrix(0, pstar, n_free)

  col <- 1L


  if (!wmat_flag) {
    # Lambda [i,j]: dS0[r,s] = LG[s,j]*I(r==i) + LG[r,j]*I(s==i)
    if (n_lam > 0L) {
      rc  <- vec2rc(m_lambda_idx, nvar)
      i_v <- rc[, 1L]
      j_v <- rc[, 2L]

      t1 <- lg[c_s, j_v, drop = FALSE] * outer(r_s, i_v, `==`)
      t2 <- lg[r_s, j_v, drop = FALSE] * outer(c_s, i_v, `==`)

      dx <- t1 + t2
      if (delta_flag) {
        dx <- dx * delta_weight
      }
      jac_sigma[, col:(col + n_lam - 1L)] <- dx
      col <- col + n_lam
    }

    # Beta [i,j]: dS0[r,s] = M[r,i]*LG[s,j] + LG[r,j]*M[s,i]
    if (n_bet > 0L) {
      rc  <- vec2rc(m_beta_idx, nfac)
      i_v <- rc[, 1L]
      j_v <- rc[, 2L]

      t1 <- m[r_s,  i_v, drop = FALSE] * lg[c_s, j_v, drop = FALSE]
      t2 <- lg[r_s, j_v, drop = FALSE] * m[c_s,  i_v, drop = FALSE]

      dx <- t1 + t2
      if (delta_flag) {
        dx <- dx * delta_weight
      }
      jac_sigma[, col:(col + n_bet - 1L)] <- dx
      col <- col + n_bet
    }
   # Psi [k,l] symmetric: dS0[r,s] = M[r,k]*M[s,l] + M[r,l]*M[s,k]
    # diagonal (k==l): halve
    if (n_psi > 0L) {
      rc  <- vec2rc(m_psi_idx, nfac)
      k_v <- pmax(rc[, 1L], rc[, 2L])
      l_v <- pmin(rc[, 1L], rc[, 2L])

      t1 <- m[r_s, k_v, drop = FALSE] * m[c_s, l_v, drop = FALSE]
      t2 <- m[r_s, l_v, drop = FALSE] * m[c_s, k_v, drop = FALSE]
      dx <- t1 + t2

      diag_mask <- (k_v == l_v)
      if (any(diag_mask)) {
        dx[, diag_mask] <- dx[, diag_mask] * 0.5
      }

      if (delta_flag) {
        dx <- dx * delta_weight
      }
      jac_sigma[, col:(col + n_psi - 1L)] <- dx
      col <- col + n_psi
    }

    # Theta [k,l] symmetric: dS = Delta %*% dTheta %*% Delta
    #  -> unit vector at vech pos (k,l)
    if (n_the > 0L) {
      rc  <- vec2rc(m_theta_idx, nvar)
      k_v <- pmax(rc[, 1L], rc[, 2L])
      l_v <- pmin(rc[, 1L], rc[, 2L])

      jt <- matrix(0, pstar, n_the)
      vech_pos <- sigma_lut[cbind(k_v, l_v)]
      if (delta_flag) {
        jt[cbind(vech_pos, seq_len(n_the))] <- delta_weight[vech_pos]
      } else {
        jt[cbind(vech_pos, seq_len(n_the))] <- 1
      }

      jac_sigma[, col:(col + n_the - 1L)] <- jt
      col <- col + n_the
    }

    # Delta[k] (diagonal only):
    # dS[r,s] = I(r==k)*Sigma0[k,s]*delta[s] + delta[r]*Sigma0[r,k]*I(s==k)
    if (n_del > 0L) {
      sigma0 <- m %*% mlist$psi %*% t(m) + mlist$theta
      k_v <- m_delta_idx
      t1 <- sigma0[c_s, k_v, drop = FALSE] * delta[c_s] * outer(r_s, k_v, `==`)
      t2 <- delta[r_s] * sigma0[r_s, k_v, drop = FALSE] * outer(c_s, k_v, `==`)
      jac_sigma[, col:(col + n_del - 1L)] <- t1 + t2
    }

  } else {
    # composite ov (cov) and composite lv (clv)
    cov_idx <- which(apply(
      mlist$lambda, 1L,
      function(x) sum(x == 0) == ncol(mlist$lambda)
    ))
    clv_idx <- which(apply(
      mlist$lambda, 2L,
      function(x) sum(x == 0) == nrow(mlist$lambda)
    ))

    lw  <- mlist$lambda + mlist$wmat
    tmat <- diag(nvar)
    tmat[cov_idx, cov_idx] <- mlist$theta[cov_idx, cov_idx]
    wtw <- t(lw[, clv_idx, drop = FALSE]) %*% tmat %*%
             lw[, clv_idx, drop = FALSE]
    wtw_inv <- solve(wtw)
    wtw_inv_1 <- diag(nfac)
    wtw_inv_1[clv_idx, clv_idx] <- wtw_inv

    if (!beta_flag) {
      a_1 <- diag(nfac)
    }

    # C0: V(eta) with composite diagonal entries zeroed
    c0 <- g
    diag(c0)[clv_idx] <- 0

    # precompute
    lambdafull <- tmat %*% lw %*% wtw_inv_1
    lg_c <- lambdafull %*% c0 %*% wtw_inv_1
    h_c <- lambdafull %*% a_1
    f_c <- lambdafull %*% g
    sig0s <- lambdafull %*% c0 %*% t(lambdafull)
    lamc <- tmat %*% lw[, clv_idx, drop = FALSE] %*% wtw_inv
    # dLambdafull/dW[i,j] for i in ovc, j in clv expands to
    #   wtw.inv[j_loc,c_loc]*(Tmat - LAMc %*% wtw %*% t(LAMc))[r,i]
    #   - LAMc[i,c_loc]*LAMc[r,j_loc].
    # The "Pmat1" and "Pmat2" multipliers below feed into the wmat-block
    # T1+T2-T3-T4 expansion. Earlier versions used Tmat - LAMc %*% t(LAMc)
    # and LAMc %*% wtw.inv respectively, which are correct only when
    # wtw == I; the generally correct forms use the wtw factor explicitly:
    pmat1 <- tmat - lamc %*% wtw %*% t(lamc)
    pmat2 <- lamc

    # Beta/Psi correction matrices
    xrs_xss <- lambdafull[r_s, clv_idx, drop = FALSE] *
              lambdafull[c_s, clv_idx, drop = FALSE]
    a_sub  <- a_1[clv_idx, , drop = FALSE]
    g_sub  <- g[clv_idx, , drop = FALSE]
    # lambda
    if (n_lam > 0L) {
      rc  <- vec2rc(m_lambda_idx, nvar)
      i_v <- rc[, 1L]
      j_v <- rc[, 2L]

      t1 <- lg_c[c_s, j_v, drop = FALSE] * outer(r_s, i_v, `==`)
      t2 <- lg_c[r_s, j_v, drop = FALSE] * outer(c_s, i_v, `==`)

      dx <- t1 + t2
      if (delta_flag) {
        dx <- dx * delta_weight
      }
      jac_sigma[, col:(col + n_lam - 1L)] <- dx
      col <- col + n_lam
    }

    # wmat
    if (n_wmat > 0L) {
      rc <- vec2rc(m_wmat_idx, nvar)
      i_v <- rc[, 1L]
      j_v <- rc[, 2L]
      j_loc <- match(j_v, clv_idx)

      # path 1: xi_i terms
      t1 <- pmat1[r_s, i_v, drop = FALSE] * lg_c[c_s, j_v, drop = FALSE]
      t2 <- pmat1[c_s, i_v, drop = FALSE] * lg_c[r_s, j_v, drop = FALSE]

      # path 2: phi_j * sigma_i correction
      t3 <- pmat2[r_s, j_loc, drop = FALSE] * sig0s[c_s, i_v, drop = FALSE]
      t4 <- pmat2[c_s, j_loc, drop = FALSE] * sig0s[r_s, i_v, drop = FALSE]

      dx <- t1 + t2 - t3 - t4
      if (delta_flag) {
        dx <- dx * delta_weight
      }
      jac_sigma[, col:(col + n_wmat - 1L)] <- dx
      col <- col + n_wmat
    }

    # beta
    if (n_bet > 0L) {
      rc  <- vec2rc(m_beta_idx, nfac)
      i_v <- rc[, 1L]
      j_v <- rc[, 2L]

      t1 <- h_c[r_s, i_v, drop = FALSE] * f_c[c_s, j_v, drop = FALSE]
      t2 <- f_c[r_s, j_v, drop = FALSE] * h_c[c_s, i_v, drop = FALSE]

      ag     <- a_sub[, i_v, drop = FALSE] * g_sub[, j_v, drop = FALSE]
      r_corr <- 2 * xrs_xss %*% ag

      dx <- t1 + t2 - r_corr
      if (delta_flag) {
        dx <- dx * delta_weight
      }
      jac_sigma[, col:(col + n_bet - 1L)] <- dx
      col <- col + n_bet
    }

    # psi
    if (n_psi > 0L) {
      rc <- vec2rc(m_psi_idx, nfac)
      k_v <- pmax(rc[, 1L], rc[, 2L])
      l_v <- pmin(rc[, 1L], rc[, 2L])

      t1 <- h_c[r_s, k_v, drop = FALSE] * h_c[c_s, l_v, drop = FALSE]
      t2 <- h_c[r_s, l_v, drop = FALSE] * h_c[c_s, k_v, drop = FALSE]

      a_kl <- a_sub[, k_v, drop = FALSE] * a_sub[, l_v, drop = FALSE]
      r_corr <- 2 * xrs_xss %*% a_kl

      dx <- t1 + t2 - r_corr

      diag_mask <- (k_v == l_v)
      if (any(diag_mask)) {
        dx[, diag_mask] <- dx[, diag_mask] * 0.5
      }

      if (delta_flag) {
        dx <- dx * delta_weight
      }
      jac_sigma[, col:(col + n_psi - 1L)] <- dx
      col <- col + n_psi
    }

    # theta
    if (n_the > 0L) {
      rc  <- vec2rc(m_theta_idx, nvar)
      k_v <- pmax(rc[, 1L], rc[, 2L])
      l_v <- pmin(rc[, 1L], rc[, 2L])

      jt <- matrix(0, pstar, n_the)
      vech_pos <- sigma_lut[cbind(k_v, l_v)]
      if (delta_flag) {
        jt[cbind(vech_pos, seq_len(n_the))] <- delta_weight[vech_pos]
      } else {
        jt[cbind(vech_pos, seq_len(n_the))] <- 1
      }

      jac_sigma[, col:(col + n_the - 1L)] <- jt
      col <- col + n_the
    }

    # delta
    if (n_del > 0L) {
      sigma0_c <- sig0s + mlist$theta
      k_v      <- m_delta_idx
      t1 <- sigma0_c[c_s, k_v, drop = FALSE] * delta[c_s] *
                                              outer(r_s, k_v, `==`)
      t2 <- delta[r_s] * sigma0_c[r_s, k_v, drop = FALSE] *
                                              outer(c_s, k_v, `==`)
      jac_sigma[, col:(col + n_del - 1L)] <- t1 + t2
    }

  } # wmat.flag

  # theta parameterization (categorical / correlation only).
  #
  # In the "theta" approach (Muthen & Asparouhov, Mplus Web Note 4), the
  # diagonal of THETA is a free parameter and the diagonal scaling factor
  # Delta is *not* a free parameter; it is determined implicitly by
  #     Delta_j^{-2} = Sigma*_jj   (j ordinal),     Delta_j = 1   (j numeric).
  # The implied moments are
  #     Sigma_obs[r,s] = Sigma*[r,s] / (s_r * s_s),    s_j = sqrt(Sigma*_jj),
  #     TH_obs[t]     = (tau - nu - lambda * alpha) / s_{v(t)}.
  #
  # Up to this point, jac_sigma has been built using the prepopulated
  # MLIST$delta (= 1/s_j) values, i.e. it equals
  #     d(Delta * Sigma* * Delta) / dphi
  # treating Delta as a *constant*. The implicit dependence of Delta on
  # the other free parameters is added directly here, by writing
  #     d Sigma_obs / dphi
  #         =  (1/(s_r s_s)) * dSigma*_rs/dphi
  #          - 0.5 * Sigma_obs[r,s]/Sigma*_rr * dSigma*_rr/dphi   (r ordinal)
  #          - 0.5 * Sigma_obs[r,s]/Sigma*_ss * dSigma*_ss/dphi   (s ordinal)
  # The first term is already in jac_sigma (the Delta=const value). The two
  # correction terms are added below using jac_sigma diagonal rows
  # (which equal (1/Sigma*_jj) * dSigma*_jj/dphi for ordinal j and
  #  dSigma*_jj/dphi for numeric j -- both consistent with the formula
  #  above when combined with sigma_obs_v as the multiplier).
  jac_diag_theta <- NULL
  if (parameterization == "theta" && (categorical || correlation)) {
    sigma_star <- lav_lisrel_sigma(mlist = mlist, delta = FALSE)
    s_diag <- sqrt(diag(sigma_star))
    if (length(num_idx) > 0L) {
      s_diag[num_idx] <- 1.0
    }
    sigma_obs <- sigma_star / tcrossprod(s_diag)
    sigma_obs_v <- lav_matrix_vech(sigma_obs)

    diag_pos <- lav_matrix_diagh_idx(nvar)
    jac_diag_theta <- jac_sigma[diag_pos, , drop = FALSE]

    # If Delta itself is in the (compact) column block (e.g. the modindices
    # "extended" model where every model-matrix element is free), Sigma*
    # does not depend on Delta -- so dSigma*_jj/d(Delta_k) = 0. The chain
    # rule below uses jac_sigma's diagonal rows as a proxy for
    # (1/Sigma*_jj) * dSigma*_jj/dphi (correct for lambda/beta/psi/theta
    # columns), but those rows pick up the delta-block contribution of
    # jac_sigma[(r,r),:], which would spuriously feed back into the chain.
    # Zero those columns so the chain term at the delta columns is zero.
    if (n_del > 0L) {
      jac_diag_theta[,
        (ncol(jac_diag_theta) - n_del + 1L):ncol(jac_diag_theta)] <- 0
    }

    is_ord <- rep(TRUE, nvar)
    if (length(num_idx) > 0L) {
      is_ord[num_idx] <- FALSE
    }

    coef_r <- -0.5 * sigma_obs_v * is_ord[r_s]
    coef_s <- -0.5 * sigma_obs_v * is_ord[c_s]

    jac_sigma <- jac_sigma +
                 coef_r * jac_diag_theta[r_s, , drop = FALSE] +
                 coef_s * jac_diag_theta[c_s, , drop = FALSE]
  }

  # categorical: reorder
  if (categorical) {
    # reorder: first variances (of numeric), then covariances
    cov_idx <- lav_matrix_vech_idx(nvar)
    covd_idx <- lav_matrix_vech_idx(nvar, diagonal = FALSE)

    var_idx <- which(is.na(match(
      cov_idx,
      covd_idx
    )))[num_idx]
    cor_idx <- match(covd_idx, cov_idx)

    jac_sigma <- rbind(
      jac_sigma[var_idx, , drop = FALSE],
      jac_sigma[cor_idx, , drop = FALSE]
    )
  }

  # correlation structure
  if (!categorical && correlation) {
    rm_idx <- lav_matrix_diagh_idx(nvar)
    jac_sigma <- jac_sigma[-rm_idx, , drop = FALSE]
  }

  jac_th   <- NULL
  jac_beta <- NULL
  nth_full <- 0L
  nexo_g   <- 0L
  if (conditional_x && !is.null(mlist$gamma)) {
    nexo_g <- ncol(mlist$gamma)
  }
  if (!categorical) {
    if (conditional_x) {
      # conditional.x continuous case: implied moments are
      #   mu_i  = nu_i + (LAMBDA(I-B)^{-1} alpha)_i
      #   pi_{i,k} = (LAMBDA(I-B)^{-1} GAMMA)_{i,k}
      #   Sigma = LAMBDA(I-B)^{-1} PSI (I-B)^{-T} LAMBDA' + THETA
      # The Beta-style row layout is interleaved column-major of a
      # (nexo+1) x nvar matrix (intercept first, then nexo slopes),
      # matching what lav_mvreg_scores_* expects.
      if (beta_flag) {
        a <- drop(a_1 %*% mlist$alpha)
      } else {
        a <- if (!is.null(mlist$alpha)) as.numeric(mlist$alpha) else
             rep(0, nfac)
      }
      gamma_ib <- if (nexo_g > 0L) {
        if (beta_flag) a_1 %*% mlist$gamma else mlist$gamma
      } else {
        matrix(0, nfac, 0L)
      }

      n_beta_rows <- (nexo_g + 1L) * nvar
      mu_slots    <- (seq_len(nvar) - 1L) * (nexo_g + 1L) + 1L

      n_free_beta <- n_nu + n_lam + n_bet + n_alp + n_gam
      jac_beta <- matrix(0, n_beta_rows, n_free_beta)

      col <- 1L
      # nu[i]: only affects mu_i (intercept slot for variable i).
      if (n_nu > 0L) {
        jac_beta[cbind(mu_slots[m_nu_idx], col + seq_len(n_nu) - 1L)] <- 1
        col <- col + n_nu
      }

      # lambda[i,j]: dmu_i = a[j] * I(r==i);
      #                                dpi_{i,k} = gamma_IB[j,k] * I(r==i).
      if (n_lam > 0L) {
        rc <- vec2rc(m_lambda_idx, nvar)
        i_v <- rc[, 1L]
        j_v <- rc[, 2L]
        # mu-row: a[j_v] at row mu_slots[i_v]
        jac_beta[cbind(mu_slots[i_v], col + seq_len(n_lam) - 1L)] <- a[j_v]
        # pi-rows (k = 1..nexo)
        if (nexo_g > 0L) {
          for (k in seq_len(nexo_g)) {
            jac_beta[cbind(mu_slots[i_v] + k,
                           col + seq_len(n_lam) - 1L)] <- gamma_ib[j_v, k]
          }
        }
        col <- col + n_lam
      }

      # beta[i,j]: dmu_r = M[r,i]*a[j]; dpi_{r,k} = M[r,i]*gamma_IB[j,k].
      if (n_bet > 0L) {
        rc <- vec2rc(m_beta_idx, nfac)
        i_v <- rc[, 1L]
        j_v <- rc[, 2L]
        # mu rows
        jac_beta[mu_slots, col:(col + n_bet - 1L)] <-
          m[, i_v, drop = FALSE] * rep(a[j_v], each = nvar)
        # pi rows
        if (nexo_g > 0L) {
          for (k in seq_len(nexo_g)) {
            jac_beta[mu_slots + k, col:(col + n_bet - 1L)] <-
              m[, i_v, drop = FALSE] *
                rep(gamma_ib[j_v, k], each = nvar)
          }
        }
        col <- col + n_bet
      }

      # alpha[k]: only affects mu (dmu = M[,k]).
      if (n_alp > 0L) {
        jac_beta[mu_slots, col:(col + n_alp - 1L)] <-
          m[, m_alpha_idx, drop = FALSE]
        col <- col + n_alp
      }

      # gamma[i,k]: only affects pi_{r,k} -> M[r,i].
      if (n_gam > 0L) {
        rc <- vec2rc(m_gamma_idx, nfac)
        i_v <- rc[, 1L]
        k_v <- rc[, 2L]
        for (c2 in seq_len(n_gam)) {
          jac_beta[mu_slots + k_v[c2], col + c2 - 1L] <- m[, i_v[c2]]
        }
      }
    } else if (meanstructure) {
      # precompute: IB.inv * alpha
      if (beta_flag) {
        a <- drop(a_1 %*% mlist$alpha)
      } else {
        a <- mlist$alpha
      }

      n_free <- n_nu + n_lam + n_bet + n_alp
      jac_mean <- matrix(0, nvar, n_free)

      col <- 1L
      # nu[i]:  dmu = e_i
      if (n_nu > 0L) {
        jn <- matrix(0, nvar, n_nu)
        jn[cbind(m_nu_idx, seq_len(n_nu))] <- 1
        jac_mean[, col:(col + n_nu - 1L)] <- jn
        col <- col + n_nu
      }

      # Lambda[i,j]:  dmu[r] = a[j] * I(r == i)
      if (n_lam > 0L) {
        rc  <- vec2rc(m_lambda_idx, nvar)
        i_v <- rc[, 1L]
        j_v <- rc[, 2L]

        jl <- matrix(0, nvar, n_lam)
        jl[cbind(i_v, seq_len(n_lam))] <- a[j_v]
        jac_mean[, col:(col + n_lam - 1L)] <- jl
        col <- col + n_lam
      }

      # Beta[i,j]:  dmu[r] = M[r,i] * a[j]
      if (n_bet > 0L) {
        rc  <- vec2rc(m_beta_idx, nfac)
        i_v <- rc[, 1L]
        j_v <- rc[, 2L]

        jac_mean[, col:(col + n_bet - 1L)] <-
          m[, i_v, drop = FALSE] * rep(a[j_v], each = nvar)
        col <- col + n_bet
      }

      # alpha[k]:  dmu = M[,k]
      if (n_alp > 0L) {
        jac_mean[, col:(col + n_alp - 1L)] <- m[, m_alpha_idx, drop = FALSE]
      }
    } # meanstructure


  # categorical
  } else if (categorical) {

    th_idx_g <- th_idx
    nth_full  <- length(th_idx_g)

    # v_slot[t] = observed-variable index for TH element t
    nlev_v <- tabulate(th_idx_g, nbins = nvar)
    nlev_v[nlev_v == 0L] <- 1L
    v_slot <- rep(seq_len(nvar), times = nlev_v)

    # per-slot delta scale
    if (delta_flag) {
      delta_star <- delta[v_slot]
    } else {
      delta_star <- rep(1, nth_full)
    }

    # ord_slots[k] = position in TH vector that corresponds to tau row k
    ord_slots <- which(th_idx_g > 0L)

    # a_vec = IB.inv * alpha  (for lambda/beta derivatives)
    if (!is.null(mlist$alpha)) {
      a_vec <- if (beta_flag) {
        as.vector(a_1 %*% mlist$alpha)
      } else {
        as.vector(mlist$alpha)
      }
    } else {
      a_vec <- rep(0, nfac)
    }

    # nu_vec (intercepts; zero for purely ordinal models)
    nu_vec <- if (!is.null(mlist$nu)) as.vector(mlist$nu) else rep(0, nvar)

    # pi0 = nu + Lambda * IB.inv * alpha = nu + M * alpha
    if (!is.null(mlist$alpha)) {
      pi0 <- nu_vec + as.vector(m %*% as.vector(mlist$alpha))
    } else {
      pi0 <- nu_vec
    }

    # tau_full: tau values at ordinal TH slots, 0 at numeric TH slots
    tau_full <- numeric(nth_full)
    if (!is.null(mlist$tau) && n_th > 0L) {
      tau_full[ord_slots] <- as.vector(mlist$tau)
    }

    # unscaled TH: TH0[t] = tau_full[t] - pi0[v(t)]
    th0 <- tau_full - pi0[v_slot]

    # build jac_th  (nth_full rows, n_free_th columns)
    # column order: tau | delta | nu | lambda | beta | alpha
    n_free_th <- n_th + n_del + n_nu + n_lam + n_bet + n_alp
    jac_th <- matrix(0, nth_full, n_free_th)
    col_th <- 1L

    # tau[k]: delta_star[t] * I(t == ord_slots[k])
    if (n_th > 0L) {
      t_slots <- ord_slots[m_th_idx]
      jt <- matrix(0, nth_full, n_th)
      jt[cbind(t_slots, seq_len(n_th))] <- delta_star[t_slots]
      jac_th[, col_th:(col_th + n_th - 1L)] <- jt
      col_th <- col_th + n_th
    }

    # delta[k]: I(v_slot[t]==k) * TH0[t]
    if (n_del > 0L) {
      jac_th[, col_th:(col_th + n_del - 1L)] <-
        outer(v_slot, m_delta_idx, `==`) * th0
      col_th <- col_th + n_del
    }

    # nu[i]: -delta_star[t] * I(v_slot[t]==i)
    if (n_nu > 0L) {
      jac_th[, col_th:(col_th + n_nu - 1L)] <-
        -outer(v_slot, m_nu_idx, `==`) * delta_star
      col_th <- col_th + n_nu
    }

    # lambda[i,j]: -delta_star[t]*I(v_slot[t]==i)*a_vec[j]
    if (n_lam > 0L) {
      rc  <- vec2rc(m_lambda_idx, nvar)
      i_v <- rc[, 1L]
      j_v <- rc[, 2L]
      jac_th[, col_th:(col_th + n_lam - 1L)] <-
        -(outer(v_slot, i_v, `==`) * delta_star) *
         matrix(a_vec[j_v], nth_full, n_lam, byrow = TRUE)
      col_th <- col_th + n_lam
    }

    # beta[i,j]: -delta_star[t]*M[v_slot[t],i]*a_vec[j]
    if (n_bet > 0L) {
      rc  <- vec2rc(m_beta_idx, nfac)
      i_v <- rc[, 1L]
      j_v <- rc[, 2L]
      jac_th[, col_th:(col_th + n_bet - 1L)] <-
        -delta_star * m[v_slot, i_v, drop = FALSE] *
         matrix(a_vec[j_v], nth_full, n_bet, byrow = TRUE)
         matrix(a_vec[j_v], nth_full, n_bet, byrow = TRUE)
      col_th <- col_th + n_bet
    }

    # alpha[k]: -delta_star[t] * M[v_slot[t], k]
    if (n_alp > 0L) {
     jac_th[, col_th:(col_th + n_alp - 1L)] <-
       -delta_star * m[v_slot, m_alpha_idx, drop = FALSE]
    }
  }

  # categorical + conditional.x: build jac_pi (slope rows).
  #   pi_obs[r, k] = Delta[r] * (LAMBDA(I-B)^{-1} GAMMA)[r, k]
  # vec(PI) row layout is column-major: index (r, k) -> (k-1)*nvar + r.
  # Free params (in this column order): lambda | beta | gamma | delta.
  jac_pi <- NULL
  n_pi_rows <- 0L
  if (categorical && conditional_x && nexo_g > 0L) {
    n_pi_rows <- nvar * nexo_g

    gamma_ib <- if (beta_flag) a_1 %*% mlist$gamma else mlist$gamma
    delta_var <- if (delta_flag) delta else rep(1, nvar)

    n_free_pi <- n_lam + n_bet + n_gam + n_del
    jac_pi <- matrix(0, n_pi_rows, n_free_pi)
    col_pi <- 1L

    # lambda[i,j]: dpi[r,k] = Delta[r] * I(r==i) * gamma_IB[j, k]
    if (n_lam > 0L) {
      rc <- vec2rc(m_lambda_idx, nvar)
      i_v <- rc[, 1L]
      j_v <- rc[, 2L]
      for (k in seq_len(nexo_g)) {
        jac_pi[cbind((k - 1L) * nvar + i_v,
                     col_pi + seq_len(n_lam) - 1L)] <-
          delta_var[i_v] * gamma_ib[j_v, k]
      }
      col_pi <- col_pi + n_lam
    }

    # beta[i,j]: dpi[r,k] = Delta[r] * M[r,i] * gamma_IB[j, k]
    if (n_bet > 0L) {
      rc <- vec2rc(m_beta_idx, nfac)
      i_v <- rc[, 1L]
      j_v <- rc[, 2L]
      for (k in seq_len(nexo_g)) {
        rows_k <- ((k - 1L) * nvar + 1L):(k * nvar)
        jac_pi[rows_k, col_pi:(col_pi + n_bet - 1L)] <-
          delta_var * m[, i_v, drop = FALSE] *
            rep(gamma_ib[j_v, k], each = nvar)
      }
      col_pi <- col_pi + n_bet
    }

    # gamma[i,k']: dpi[r,k] = Delta[r] * M[r,i] * I(k==k')
    if (n_gam > 0L) {
      rc <- vec2rc(m_gamma_idx, nfac)
      i_v <- rc[, 1L]
      k_v <- rc[, 2L]
      for (c2 in seq_len(n_gam)) {
        rows_k <- ((k_v[c2] - 1L) * nvar + 1L):(k_v[c2] * nvar)
        jac_pi[rows_k, col_pi + c2 - 1L] <- delta_var * m[, i_v[c2]]
      }
      col_pi <- col_pi + n_gam
    }

    # delta[r] (delta-parameterization free):
    #                               dpi[r,k] = (M GAMMA)[r, k] (unscaled).
    if (n_del > 0L) {
      m_gamma <- m %*% mlist$gamma  # nvar x nexo
      for (c2 in seq_len(n_del)) {
        kd <- m_delta_idx[c2]
        for (k in seq_len(nexo_g)) {
          jac_pi[(k - 1L) * nvar + kd, col_pi + c2 - 1L] <- m_gamma[kd, k]
        }
      }
    }
  }

  # sigma
  out <- matrix(0, nrow = nrow(jac_sigma), ncol = nx_free)
  el_idx_sigma <- if (!wmat_flag) {
    c(x_lambda_idx, x_beta_idx, x_psi_idx, x_theta_idx, x_delta_idx)
  } else {
    c(x_lambda_idx, x_wmat_idx, x_beta_idx, x_psi_idx, x_theta_idx, x_delta_idx)
  }
  out[, el_idx_sigma] <- jac_sigma

  # meanstructure / conditional.x
  if (!categorical && conditional_x && !is.null(jac_beta)) {
    el_idx_beta <- c(x_nu_idx, x_lambda_idx, x_beta_idx, x_alpha_idx,
                     x_gamma_idx)
    outm <- matrix(0, nrow = nrow(jac_beta), ncol = nx_free)
    outm[, el_idx_beta] <- jac_beta
    out <- rbind(outm, out)
  } else if (!categorical && meanstructure) {
    # right order
    el_idx <- c(x_nu_idx, x_lambda_idx, x_beta_idx, x_alpha_idx)
    outm <- matrix(0, nrow = nrow(jac_mean), ncol = nx_free)
    outm[, el_idx] <- jac_mean
    out <- rbind(outm, out)
  } else if (categorical && !is.null(jac_th)) {
    el_idx_th <- c(
      if (n_th  > 0L) x_th_idx     else integer(0L),
      if (n_del > 0L) x_delta_idx  else integer(0L),
      if (n_nu  > 0L) x_nu_idx     else integer(0L),
      if (n_lam > 0L) x_lambda_idx else integer(0L),
      if (n_bet > 0L) x_beta_idx   else integer(0L),
      if (n_alp > 0L) x_alpha_idx  else integer(0L)
    )
    out_th <- matrix(0, nth_full, nx_free)
    out_th[, el_idx_th] <- jac_th

    # jac_diag_full only needed under theta parameterization (cached for
    # tau- and pi-side chain corrections).
    jac_diag_full <- NULL
    is_ord_v <- rep(TRUE, nvar)
    if (length(num_idx) > 0L) {
      is_ord_v[num_idx] <- FALSE
    }
    if (parameterization == "theta" && !is.null(jac_diag_theta)) {
      jac_diag_full <- matrix(0, nvar, nx_free)
      jac_diag_full[, el_idx_sigma] <- jac_diag_theta
    }

    # theta parameterization: chain rule for tau-side moments.
    #   TH_obs[t] = (tau - nu - lambda*alpha) / s_{v(t)},  v(t) ordinal,
    # was filled treating Delta as a constant. The implicit dependence
    # of Delta on the other free parameters adds
    #   d TH_obs[t] / dphi += -0.5 * TH_obs[t] / Sigma*_v(t),v(t)
    #                              * dSigma*_v(t),v(t) / dphi
    # which, in nx.free column space, is
    #   chain[t, :] = -0.5 * TH_obs[t] * jac_diag_full[v_slot[t], :]
    # Numeric phantom slots get no correction (Delta_num = 1 by construction,
    # so it does not depend on the free parameters).
    if (!is.null(jac_diag_full) && nth_full > 0L) {
      th_obs <- delta_star * th0 * is_ord_v[v_slot]
      out_th <- out_th +
                (-0.5 * th_obs) * jac_diag_full[v_slot, , drop = FALSE]
    }

    out <- rbind(out_th, out)

    # categorical + conditional.x: insert pi rows between th and sigma.
    if (!is.null(jac_pi)) {
      el_idx_pi <- c(
        if (n_lam > 0L) x_lambda_idx else integer(0L),
        if (n_bet > 0L) x_beta_idx   else integer(0L),
        if (n_gam > 0L) x_gamma_idx  else integer(0L),
        if (n_del > 0L) x_delta_idx  else integer(0L)
      )
      out_pi <- matrix(0, n_pi_rows, nx_free)
      out_pi[, el_idx_pi] <- jac_pi

      # theta parameterization: chain rule for pi-side moments.
      #   pi_obs[r, k] = Delta_r * unscaled_pi[r, k], with Delta_r = 1/s_r
      #   d pi_obs[r, k] / dphi += -0.5 * pi_obs[r, k] / Sigma*_rr * dSigma*_rr
      # only for ordinal r. Numeric variables have Delta_r = 1 (constant).
      if (!is.null(jac_diag_full)) {
        # construct vec(PI_obs) following the column-major row layout
        unscaled_pi <- m %*% mlist$gamma
        pi_obs_mat  <- if (delta_flag) unscaled_pi * delta else unscaled_pi
        pi_obs_v <- as.numeric(pi_obs_mat)              # nvar*nexo, vec
        r_pi     <- rep.int(seq_len(nvar), nexo_g)      # variable per row
        out_pi <- out_pi +
                  (-0.5 * pi_obs_v * is_ord_v[r_pi]) *
                    jac_diag_full[r_pi, , drop = FALSE]
      }

      # row order: out_th already on top; pi between th and sigma.
      # 'out' currently is rbind(out_th, sigma); reassemble:
      n_th_rows <- nrow(out_th)
      n_sig_rows <- nrow(out) - n_th_rows
      out <- rbind(out[seq_len(n_th_rows), , drop = FALSE],
                   out_pi,
                   out[n_th_rows + seq_len(n_sig_rows), , drop = FALSE])
    }
  }

  # composite chain-rule correction.
  # lav_lisrel_comp_set_intresvar() implicitly reparameterizes psi[clv,clv]
  # diagonal (and alpha[clv] when meanstructure is present) as a function of
  # the other free parameters (wmat / nu / theta / beta / lambda). The
  # analytical jac_sigma / jac_mean above were computed treating psi*/alpha*
  # as inputs; we add the chain term:
  #   chain[, j] = sum_{k in clv} dSigma/dpsi*[k,k] * dpsi*[k,k]/dphi_j
  #              + sum_{k in clv} dMu/dalpha*[k]   * dalpha*[k]/dphi_j
  # The closed-form dSigma/dpsi*[k,k] is the same expression used for free
  # psi[k,k] (T1+T2-R_corr halved), reused with k ranging over clv.idx;
  # dMu/dalpha*[k] = M[r,k] (zero for composite-OV rows by construction).
  # dpsi*/dphi and dalpha*/dphi are obtained by complex-step differentiation
  # through lav_lisrel_comp_set_intresvar() itself.
  if (wmat_flag && length(clv_idx) > 0L) {
    n_lvc <- length(clv_idx)

    # current free parameter vector
    x0_chain <- numeric(nx_free)
    for (mm in seq_along(mlist)) {
      if (length(m_free_idx[[mm]])) {
        x0_chain[x_free_idx[[mm]]] <- mlist[[mm]][m_free_idx[[mm]]]
      }
    }

    ml_tmpl    <- mlist
    m_free_cap <- m_free_idx
    x_free_cap <- x_free_idx
    clv_cap    <- clv_idx
    meanstr_cap <- meanstructure
    compute_lvc <- function(xx) {
      ml <- ml_tmpl
      for (mm in seq_along(ml)) {
        if (length(m_free_cap[[mm]])) {
          ml[[mm]][m_free_cap[[mm]]] <- xx[x_free_cap[[mm]]]
        }
      }
      ml <- lav_lisrel_comp_set_intresvar(ml)
      o <- ml$psi[cbind(clv_cap, clv_cap)]
      if (meanstr_cap && !is.null(ml$alpha)) {
        o <- c(o, ml$alpha[clv_cap, 1L])
      }
      o
    }
    j_lvc <- lav_func_jacobian_complex(func = compute_lvc, x = x0_chain)
    j_psi <- j_lvc[seq_len(n_lvc), , drop = FALSE]

    # dSigma/dpsi*[k,k] for k in clv.idx (closed form):
    #   H_c[r,k]*H_c[s,k] -
    #             sum_{k' in clv} Lambdafull[r,k']*Lambdafull[s,k']*A[k',k]^2
    h_c_lvc <- h_c[, clv_idx, drop = FALSE]
    m_psi   <- h_c_lvc[r_s, , drop = FALSE] *
                 h_c_lvc[c_s, , drop = FALSE]
    a_lvc_lvc <- a_1[clv_idx, clv_idx, drop = FALSE]
    m_psi   <- m_psi - xrs_xss %*% (a_lvc_lvc^2)

    sigma_chain <- m_psi %*% j_psi

    if (!categorical && meanstructure) {
      sigma_row_offset <- nvar
    } else if (categorical && !is.null(jac_th)) {
      sigma_row_offset <- nth_full
    } else {
      sigma_row_offset <- 0L
    }
    sig_rows <- sigma_row_offset + seq_len(nrow(jac_sigma))
    out[sig_rows, ] <- out[sig_rows, , drop = FALSE] + sigma_chain

    # dMu/dalpha*[k] for k in clv.idx: column M[, k]
    if (meanstructure && !categorical) {
      j_alpha  <- j_lvc[n_lvc + seq_len(n_lvc), , drop = FALSE]
      mu_chain <- m[, clv_idx, drop = FALSE] %*% j_alpha
      out[seq_len(nvar), ] <- out[seq_len(nvar), , drop = FALSE] + mu_chain
    }
  }

  # group weight: prepend a row with 1.0 at the gw column
  if (group_w_free) {
    mm_gw_idx <- which(mnames == "gw")
    if (length(mm_gw_idx) > 0L) {
      x_gw_idx <- x_free_idx[[mm_gw_idx]]
      out_gw <- matrix(0, nrow = 1L, ncol = nx_free)
      out_gw[1L, x_gw_idx] <- 1.0
      out <- rbind(out_gw, out)
    }
  }

  out
}
