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
                       allow.composites = TRUE,
                       remove.nonexisting = TRUE) {
  # prepare target list
  if (is.null(target)) target <- lavpartable

  stopifnot(!is.null(target$block))

  # prepare output
  N <- length(target$lhs)
  tmp.mat <- character(N)
  tmp.row <- integer(N)
  tmp.col <- integer(N)

  # global settings
  meanstructure <- any(lavpartable$op == "~1")
  categorical <- any(lavpartable$op == "|")
  composites <- any(lavpartable$op == "<~") && allow.composites
  group.w.free <- any(lavpartable$lhs == "group" & lavpartable$op == "%")

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

  ov.dummy.names.nox <- vector("list", nblocks)
  ov.dummy.names.x <- vector("list", nblocks)
  if (extra) {
    REP.mmNames <- vector("list", nblocks)
    REP.mmNumber <- vector("list", nblocks)
    REP.mmRows <- vector("list", nblocks)
    REP.mmCols <- vector("list", nblocks)
    REP.mmDimNames <- vector("list", nblocks)
    REP.mmSymmetric <- vector("list", nblocks)
  }

  for (g in 1:nblocks) {
    # info from user model per block
    if (gamma) {
      ov.names <- lav_partable_vnames(lavpartable, "ov.nox", block = g)
    } else {
      ov.names <- lav_partable_vnames(lavpartable, "ov", block = g)
    }
    nvar <- length(ov.names)
    lv.names <- lav_partable_vnames(lavpartable, "lv", block = g)
    nfac <- length(lv.names)
    ov.th <- lav_partable_vnames(lavpartable, "th", block = g)
    nth <- length(ov.th)
    ov.names.x <- lav_partable_vnames(lavpartable, "ov.x", block = g)
    nexo <- length(ov.names.x)
    ov.names.nox <- lav_partable_vnames(lavpartable, "ov.nox", block = g)

    # in this representation, we need to create 'phantom/dummy' latent
    # variables for all `x' and `y' variables not in lv.names
    # (only y if conditional.x = TRUE)

    # regression dummys
    if (gamma) {
      tmp.names <-
        unique(lavpartable$lhs[(lavpartable$op == "~" |
          lavpartable$op == "<~") &
          lavpartable$block == g])
      # new in 0.6-12: fix for multilevel + conditional.x: split ov.x
      # are removed from ov.x
      if (nlevels > 1L) {
        if (ngroups == 1L) {
          OTHER.BLOCK.NAMES <- lav_partable_vnames(lavpartable, "ov",
            block = seq_len(nblocks)[-g]
          )
        } else {
          # TEST ME
          this.group <- ceiling(g / nlevels)
          blocks.within.group <- (this.group - 1L) * nlevels + seq_len(nlevels)
          OTHER.BLOCK.NAMES <- lav_partable_vnames(lavpartable,
            "ov",
            block = blocks.within.group[-g]
          )
        }
        if (length(ov.names.x) > 0L) {
          idx <- which(ov.names.x %in% OTHER.BLOCK.NAMES)
          if (length(idx) > 0L) {
            tmp.names <- unique(c(tmp.names, ov.names.x[idx]))
            ov.names.nox <- unique(c(ov.names.nox, ov.names.x[idx]))
            ov.names.x <- ov.names.x[-idx]
            nexo <- length(ov.names.x)
            ov.names <- ov.names.nox
            nvar <- length(ov.names)
          }
        }
      }
    } else {
      if (composites) {
        tmp.names <-
          unique(c(
            lavpartable$lhs[(lavpartable$op == "~") &
              lavpartable$block == g],
            lavpartable$rhs[(lavpartable$op == "~") &
              lavpartable$block == g]
          ))
      } else {
        # old behavior < 0.6-20
        tmp.names <-
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
    dummy.names1 <- tmp.names[!tmp.names %in% lv.names]
    # covariances involving dummys
    dummy.cov.idx <- which(lavpartable$op == "~~" & lavpartable$block == g &
      (lavpartable$lhs %in% dummy.names1 |
        lavpartable$rhs %in% dummy.names1))
    # new in 0.5-21: also include covariances involving these covariances...
    dummy.cov.idx1 <- which(lavpartable$op == "~~" & lavpartable$block == g &
      (lavpartable$lhs %in% lavpartable$lhs[dummy.cov.idx] |
        lavpartable$rhs %in% lavpartable$rhs[dummy.cov.idx]))
    dummy.cov.idx <- unique(c(dummy.cov.idx, dummy.cov.idx1))

    dummy.names2 <- unique(c(
      lavpartable$lhs[dummy.cov.idx],
      lavpartable$rhs[dummy.cov.idx]
    ))


    # new in 0.6-7: ~~ between latent and observed
    dummy.cov.ov.lv.idx1 <- which(lavpartable$op == "~~" &
      lavpartable$block == g &
      lavpartable$lhs %in% ov.names &
      lavpartable$rhs %in% lv.names)
    dummy.cov.ov.lv.idx2 <- which(lavpartable$op == "~~" &
      lavpartable$block == g &
      lavpartable$lhs %in% lv.names &
      lavpartable$rhs %in% ov.names)
    dummy.names3 <- unique(c(
      lavpartable$lhs[dummy.cov.ov.lv.idx1],
      lavpartable$rhs[dummy.cov.ov.lv.idx2]
    ))

    # new in 0.6-10: ~~ between observed and observed, but not in ~
    dummy.orphan.idx <- which(lavpartable$op == "~~" &
      lavpartable$block == g &
      lavpartable$lhs %in% ov.names &
      lavpartable$rhs %in% ov.names &
      (!lavpartable$lhs %in% c(
        dummy.names1,
        dummy.names2
      ) |
        !lavpartable$rhs %in% c(
          dummy.names1,
          dummy.names2
        )))

    # collect all dummy variables
    dummy.names <- unique(c(dummy.names1, dummy.names2, dummy.names3))


    if (length(dummy.names)) {
      # make sure order is the same as ov.names
      ov.dummy.names.nox[[g]] <-
        ov.names.nox[ov.names.nox %in% dummy.names]
      ov.dummy.names.x[[g]] <-
        ov.names.x[ov.names.x %in% dummy.names]

      # combine them, make sure order is identical to ov.names
      tmp <- ov.names[ov.names %in% dummy.names]

      # same for ov.names.x (if they are not in ov.names) (conditional.x)
      if (length(ov.names.x) > 0L) {
        tmp.x <- ov.names.x[ov.names.x %in% dummy.names]
        tmp <- unique(c(tmp, tmp.x))
      }

      # extend lv.names
      lv.names <- c(lv.names, tmp)
      nfac <- length(lv.names)

      # add 'dummy' =~ entries
      dummy.mat <- rep("lambda", length(dummy.names))
    } else {
      ov.dummy.names.nox[[g]] <- character(0)
      ov.dummy.names.x[[g]] <- character(0)
    }

    # 1a. "=~" regular indicators
    idx <- which(target$block == g &
      target$op == "=~" & !(target$rhs %in% lv.names))
    tmp.mat[idx] <- "lambda"
    tmp.row[idx] <- match(target$rhs[idx], ov.names)
    tmp.col[idx] <- match(target$lhs[idx], lv.names)

    # 1b. "=~" regular higher-order lv indicators
    idx <- which(target$block == g &
      target$op == "=~" & !(target$rhs %in% ov.names))
    tmp.mat[idx] <- "beta"
    tmp.row[idx] <- match(target$rhs[idx], lv.names)
    tmp.col[idx] <- match(target$lhs[idx], lv.names)

    # 1c. "=~" indicators that are both in ov and lv
    idx <- which(target$block == g &
      target$op == "=~" & target$rhs %in% ov.names &
      target$rhs %in% lv.names)
    tmp.mat[idx] <- "beta"
    tmp.row[idx] <- match(target$rhs[idx], lv.names)
    tmp.col[idx] <- match(target$lhs[idx], lv.names)

    # 1d. "<~" indicators
    if (composites) {
      idx <- which(target$block == g &
        target$op == "<~" & !(target$rhs %in% lv.names))
      tmp.mat[idx] <- "wmat"
      tmp.row[idx] <- match(target$rhs[idx], ov.names)
      tmp.col[idx] <- match(target$lhs[idx], lv.names)
    }

    # 2. "~" regressions
    if (gamma) {
      # gamma
      if (composites) {
        idx <- which(target$rhs %in% ov.names.x &
          target$block == g & target$op == "~")
      } else {
        idx <- which(target$rhs %in% ov.names.x &
          target$block == g & (target$op == "~" |
          target$op == "<~"))
      }
      tmp.mat[idx] <- "gamma"
      tmp.row[idx] <- match(target$lhs[idx], lv.names)
      tmp.col[idx] <- match(target$rhs[idx], ov.names.x)

      # beta
      if (composites) {
        idx <- which(!target$rhs %in% ov.names.x &
          target$block == g & target$op == "~")
      } else {
        idx <- which(!target$rhs %in% ov.names.x &
          target$block == g & (target$op == "~" |
          target$op == "<~"))
      }
      tmp.mat[idx] <- "beta"
      tmp.row[idx] <- match(target$lhs[idx], lv.names)
      tmp.col[idx] <- match(target$rhs[idx], lv.names)
    } else {
      if (composites) {
        idx <- which(target$block == g & target$op == "~")
      } else {
        idx <- which(target$block == g & (target$op == "~" |
          target$op == "<~"))
      }
      tmp.mat[idx] <- "beta"
      tmp.row[idx] <- match(target$lhs[idx], lv.names)
      tmp.col[idx] <- match(target$rhs[idx], lv.names)
    }

    # 3a. "~~" ov
    idx <- which(target$block == g &
      target$op == "~~" & !(target$lhs %in% lv.names))
    tmp.mat[idx] <- "theta"
    tmp.row[idx] <- match(target$lhs[idx], ov.names)
    tmp.col[idx] <- match(target$rhs[idx], ov.names)

    # 3aa. "~~" ov.x
    if (gamma) {
      idx <- which(target$block == g &
        target$op == "~~" & (target$lhs %in% ov.names.x))
      tmp.mat[idx] <- "cov.x"
      tmp.row[idx] <- match(target$lhs[idx], ov.names.x)
      tmp.col[idx] <- match(target$rhs[idx], ov.names.x)
    }

    # 3b. "~~" lv
    idx <- which(target$block == g &
      target$op == "~~" & target$rhs %in% lv.names)
    tmp.mat[idx] <- "psi"
    tmp.row[idx] <- match(target$lhs[idx], lv.names)
    tmp.col[idx] <- match(target$rhs[idx], lv.names)

    # 4a. "~1" ov
    idx <- which(target$block == g &
      target$op == "~1" & !(target$lhs %in% lv.names))
    tmp.mat[idx] <- "nu"
    tmp.row[idx] <- match(target$lhs[idx], ov.names)
    tmp.col[idx] <- 1L

    # 4aa, "~1" ov.x
    if (gamma) {
      idx <- which(target$block == g &
        target$op == "~1" & (target$lhs %in% ov.names.x))
      tmp.mat[idx] <- "mean.x"
      tmp.row[idx] <- match(target$lhs[idx], ov.names.x)
      tmp.col[idx] <- 1L
    }

    # 4b. "~1" lv
    idx <- which(target$block == g &
      target$op == "~1" & target$lhs %in% lv.names)
    tmp.mat[idx] <- "alpha"
    tmp.row[idx] <- match(target$lhs[idx], lv.names)
    tmp.col[idx] <- 1L

    # 5. "|" th
    LABEL <- paste(target$lhs, target$op, target$rhs, sep = "")
    idx <- which(target$block == g &
      target$op == "|" & LABEL %in% ov.th)
    TH <- paste(target$lhs[idx], "|", target$rhs[idx], sep = "")
    tmp.mat[idx] <- "tau"
    tmp.row[idx] <- match(TH, ov.th)
    tmp.col[idx] <- 1L

    # 6. "~*~" scales
    idx <- which(target$block == g &
      target$op == "~*~")
    tmp.mat[idx] <- "delta"
    tmp.row[idx] <- match(target$lhs[idx], ov.names)
    tmp.col[idx] <- 1L

    # new 0.5-12: catch lower-elements in theta/psi
    idx.lower <- which(tmp.mat %in% c("theta", "psi") & tmp.row > tmp.col)
    if (length(idx.lower) > 0L) {
      tmp <- tmp.row[idx.lower]
      tmp.row[idx.lower] <- tmp.col[idx.lower]
      tmp.col[idx.lower] <- tmp
    }

    # new 0.5-16: group weights
    idx <- which(target$block == g & target$lhs == "group" &
      target$op == "%")
    tmp.mat[idx] <- "gw"
    tmp.row[idx] <- 1L
    tmp.col[idx] <- 1L

    # new in 0.6-22: instruments
    idx <- which(target$block == g & target$lhs == "group" &
      target$op == "|~")
    tmp.mat[idx] <- "miiv"
    tmp.row[idx] <- 0L
    tmp.col[idx] <- 0L

    if (extra) {
      # mRows
      mmRows <- list(
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
      mmCols <- list(
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
      mmDimNames <- list(
        tau = list(ov.th, "threshold"),
        delta = list(ov.names, "scales"),
        nu = list(ov.names, "intercept"),
        lambda = list(ov.names, lv.names),
        wmat = list(ov.names, lv.names),
        theta = list(ov.names, ov.names),
        alpha = list(lv.names, "intercept"),
        beta = list(lv.names, lv.names),
        gamma = list(lv.names, ov.names.x),
        cov.x = list(ov.names.x, ov.names.x),
        mean.x = list(ov.names.x, "intercepts"),
        gw = list("group", "weight"),
        psi = list(lv.names, lv.names)
      )

      # isSymmetric
      mmSymmetric <- list(
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
      IDX <- which(target$block == g)
      if ("wmat" %in% tmp.mat[IDX]) {
        mmNames <- c("lambda", "wmat", "theta", "psi")
      } else {
        mmNames <- c("lambda", "theta", "psi")
      }

      if ("beta" %in% tmp.mat[IDX]) {
        mmNames <- c(mmNames, "beta")
      }
      if (meanstructure) {
        mmNames <- c(mmNames, "nu", "alpha")
      }
      if ("tau" %in% tmp.mat[IDX]) {
        mmNames <- c(mmNames, "tau")
      }
      if ("delta" %in% tmp.mat[IDX]) {
        mmNames <- c(mmNames, "delta")
      }
      if ("gamma" %in% tmp.mat[IDX]) {
        mmNames <- c(mmNames, "gamma")
      }
      if ("gw" %in% tmp.mat[IDX]) {
        mmNames <- c(mmNames, "gw")
      }
      if ("cov.x" %in% tmp.mat[IDX]) {
        mmNames <- c(mmNames, "cov.x")
      }
      if ("mean.x" %in% tmp.mat[IDX]) {
        mmNames <- c(mmNames, "mean.x")
      }

      REP.mmNames[[g]] <- mmNames
      REP.mmNumber[[g]] <- length(mmNames)
      REP.mmRows[[g]] <- unlist(mmRows[mmNames])
      REP.mmCols[[g]] <- unlist(mmCols[mmNames])
      REP.mmDimNames[[g]] <- mmDimNames[mmNames]
      REP.mmSymmetric[[g]] <- unlist(mmSymmetric[mmNames])
    } # extra
  } # nblocks

  REP <- list(
    mat = tmp.mat,
    row = tmp.row,
    col = tmp.col
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
  attr(REP, "ov.dummy.names.nox") <- ov.dummy.names.nox
  attr(REP, "ov.dummy.names.x") <- ov.dummy.names.x

  if (extra) {
    attr(REP, "mmNames") <- REP.mmNames
    attr(REP, "mmNumber") <- REP.mmNumber
    attr(REP, "mmRows") <- REP.mmRows
    attr(REP, "mmCols") <- REP.mmCols
    attr(REP, "mmDimNames") <- REP.mmDimNames
    attr(REP, "mmSymmetric") <- REP.mmSymmetric
  }

  REP
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
lav_lisrel_eeta <- function(MLIST = NULL, mean.x = NULL,
                            sample.mean = NULL,
                            ov.y.dummy.ov.idx = NULL,
                            ov.x.dummy.ov.idx = NULL,
                            ov.y.dummy.lv.idx = NULL,
                            ov.x.dummy.lv.idx = NULL) {
  BETA <- MLIST$beta
  GAMMA <- MLIST$gamma

  # ALPHA? (reconstruct, but no 'fix')
  ALPHA <- lav_lisrel_alpha0(
    MLIST = MLIST, sample.mean = sample.mean,
    ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
    ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
    ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
    ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
  )

  # BETA?
  if (!is.null(BETA)) {
    IB.inv <- lav_lisrel_ibinv(MLIST = MLIST)
    # GAMMA?
    if (!is.null(GAMMA)) {
      eeta <- as.vector(IB.inv %*% ALPHA + IB.inv %*% GAMMA %*% mean.x)
    } else {
      eeta <- as.vector(IB.inv %*% ALPHA)
    }
  } else {
    # GAMMA?
    if (!is.null(GAMMA)) {
      eeta <- as.vector(ALPHA + GAMMA %*% mean.x)
    } else {
      eeta <- as.vector(ALPHA)
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
lav_lisrel_eetax <- function(MLIST = NULL, eXo = NULL, N = nrow(eXo),
                             sample.mean = NULL,
                             ov.y.dummy.ov.idx = NULL,
                             ov.x.dummy.ov.idx = NULL,
                             ov.y.dummy.lv.idx = NULL,
                             ov.x.dummy.lv.idx = NULL) {
  LAMBDA <- MLIST$lambda
  BETA <- MLIST$beta
  GAMMA <- MLIST$gamma
  nfac <- ncol(LAMBDA)
  # if eXo, N must be nrow(eXo)
  if (!is.null(eXo)) {
    N <- nrow(eXo)
  }

  # ALPHA?
  ALPHA <- lav_lisrel_alpha0(
    MLIST = MLIST, sample.mean = sample.mean,
    ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
    ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
    ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
    ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
  )

  # construct [nobs x nfac] matrix (repeating ALPHA)
  EETA <- matrix(ALPHA, N, nfac, byrow = TRUE)

  # put back eXo values if dummy
  if (length(ov.x.dummy.lv.idx) > 0L) {
    EETA[, ov.x.dummy.lv.idx] <- eXo
  }

  # BETA?
  if (!is.null(BETA)) {
    IB.inv <- lav_lisrel_ibinv(MLIST = MLIST)
    EETA <- EETA %*% t(IB.inv)
  }

  # GAMMA?
  if (!is.null(GAMMA)) {
    if (!is.null(BETA)) {
      EETA <- EETA + eXo %*% t(IB.inv %*% GAMMA)
    } else {
      EETA <- EETA + eXo %*% t(GAMMA)
    }
  }

  EETA
}

# 3) VETA
# compute V(ETA): variances/covariances of latent variables
# - if no eXo (and GAMMA)
#     V(ETA) = (I-B)^-1 PSI (I-B)^-T
# - if eXo and GAMMA: (cfr lisrel submodel 3a with ksi=x)
#     V(ETA) = (I-B)^-1 [ GAMMA  cov.x t(GAMMA) + PSI] (I-B)^-T
lav_lisrel_veta <- function(MLIST = NULL) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  PSI <- MLIST$psi
  THETA <- MLIST$theta
  BETA <- MLIST$beta
  GAMMA <- MLIST$gamma

  if (!is.null(GAMMA)) {
    COV.X <- MLIST$cov.x
    # we treat 'x' as 'ksi' in the LISREL model; cov.x is PHI
    PSI <- tcrossprod(GAMMA %*% COV.X, GAMMA) + PSI
  }

  # beta?
  if (is.null(BETA)) {
    VETA <- PSI
  } else {
    IB.inv <- lav_lisrel_ibinv(MLIST = MLIST)
    VETA <- tcrossprod(IB.inv %*% PSI, IB.inv)
  }

  VETA
}

# 4) VETAx
# compute V(ETA|x_i): variances/covariances of latent variables
#     V(ETA) = (I-B)^-1 PSI (I-B)^-T  + remove dummies
lav_lisrel_vetax <- function(MLIST = NULL, lv.dummy.idx = NULL) {
  PSI <- MLIST$psi
  BETA <- MLIST$beta

  # beta?
  if (is.null(BETA)) {
    VETA <- PSI
  } else {
    IB.inv <- lav_lisrel_ibinv(MLIST = MLIST)
    VETA <- tcrossprod(IB.inv %*% PSI, IB.inv)
  }

  # remove dummy lv?
  if (!is.null(lv.dummy.idx)) {
    VETA <- VETA[-lv.dummy.idx, -lv.dummy.idx, drop = FALSE]
  }

  VETA
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
lav_lisrel_ey <- function(MLIST = NULL, mean.x = NULL, sample.mean = NULL,
                          ov.y.dummy.ov.idx = NULL,
                          ov.x.dummy.ov.idx = NULL,
                          ov.y.dummy.lv.idx = NULL,
                          ov.x.dummy.lv.idx = NULL, delta = TRUE) {
  LAMBDA <- MLIST$lambda

  # get NU, but do not 'fix'
  NU <- lav_lisrel_nu0(
    MLIST = MLIST, sample.mean = sample.mean,
    ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
    ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
    ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
    ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
  )

  # compute E(ETA)
  EETA <- lav_lisrel_eeta(
    MLIST = MLIST, sample.mean = sample.mean,
    mean.x = mean.x,
    ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
    ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
    ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
    ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
  )

  # EY
  EY <- as.vector(NU) + as.vector(LAMBDA %*% EETA)

  # if delta, scale
  if (delta && !is.null(MLIST$delta)) {
    EY <- EY * as.vector(MLIST$delta)
  }

  EY
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
lav_lisrel_eyetax <- function(MLIST = NULL,
                              eXo = NULL,
                              ETA = NULL,
                              N = nrow(eXo),
                              sample.mean = NULL,
                              ov.y.dummy.ov.idx = NULL,
                              ov.x.dummy.ov.idx = NULL,
                              ov.y.dummy.lv.idx = NULL,
                              ov.x.dummy.lv.idx = NULL,
                              delta = TRUE) {
  LAMBDA <- MLIST$lambda
  BETA <- MLIST$beta
  if (!is.null(eXo)) {
    N <- nrow(eXo)
  } else if (!is.null(N)) {
    # nothing to do
  } else {
    N <- 1L
  }

  # create ETA matrix
  if (nrow(ETA) == 1L) {
    ETA <- matrix(ETA, N, ncol(ETA), byrow = TRUE)
  }

  # always augment ETA with 'dummy values' (0 for ov.y, eXo for ov.x)
  # ndummy <- length(c(ov.y.dummy.lv.idx, ov.x.dummy.lv.idx))
  # if(ndummy > 0L) {
  #    ETA2 <- cbind(ETA, matrix(0, N, ndummy))
  # } else {
  ETA2 <- ETA
  # }

  # only if we have dummy ov.y, we need to compute the 'yhat' values
  # beforehand
  if (length(ov.y.dummy.lv.idx) > 0L) {
    # insert eXo values
    if (length(ov.x.dummy.lv.idx) > 0L) {
      ETA2[, ov.x.dummy.lv.idx] <- eXo
    }
    # zero ov.y values
    if (length(ov.y.dummy.lv.idx) > 0L) {
      ETA2[, ov.y.dummy.lv.idx] <- 0
    }

    # ALPHA? (reconstruct, but no 'fix')
    ALPHA <- lav_lisrel_alpha0(
      MLIST = MLIST, sample.mean = sample.mean,
      ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
      ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
      ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
      ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
    )
    # BETA?
    if (!is.null(BETA)) {
      ETA2 <- sweep(tcrossprod(ETA2, BETA), 2L, STATS = ALPHA, FUN = "+")
    } else {
      ETA2 <- sweep(ETA2, 2L, STATS = ALPHA, FUN = "+")
    }

    # put back eXo values
    if (length(ov.x.dummy.lv.idx) > 0L) {
      ETA2[, ov.x.dummy.lv.idx] <- eXo
    }

    # put back ETA values for the 'real' latent variables
    dummy.idx <- c(ov.x.dummy.lv.idx, ov.y.dummy.lv.idx)
    if (length(dummy.idx) > 0L) {
      lv.regular.idx <- seq_len(min(dummy.idx) - 1L)
      ETA2[, lv.regular.idx] <- ETA[, lv.regular.idx, drop = FALSE]
    }
  }

  # get NU, but do not 'fix'
  NU <- lav_lisrel_nu0(
    MLIST = MLIST,
    sample.mean = sample.mean,
    ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
    ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
    ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
    ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
  )

  # EYetax
  EYetax <- sweep(tcrossprod(ETA2, LAMBDA), 2L, STATS = NU, FUN = "+")

  # if delta, scale
  if (delta && !is.null(MLIST$delta)) {
    EYetax <- sweep(EYetax, 2L, STATS = MLIST$delta, FUN = "*")
  }

  EYetax
}

# unconditional version
lav_lisrel_eyetax3 <- function(MLIST = NULL,
                               ETA = NULL,
                               sample.mean = NULL,
                               mean.x = NULL,
                               ov.y.dummy.ov.idx = NULL,
                               ov.x.dummy.ov.idx = NULL,
                               ov.y.dummy.lv.idx = NULL,
                               ov.x.dummy.lv.idx = NULL,
                               delta = TRUE) {
  LAMBDA <- MLIST$lambda

  # special case: empty lambda
  if (ncol(LAMBDA) == 0L) {
    return(matrix(sample.mean,
      nrow(ETA), length(sample.mean),
      byrow = TRUE
    ))
  }

  # lv idx
  dummy.idx <- c(ov.y.dummy.lv.idx, ov.x.dummy.lv.idx)
  if (length(dummy.idx) > 0L) {
    nondummy.idx <- seq_len(min(dummy.idx) - 1L)
  } else {
    nondummy.idx <- seq_len(ncol(MLIST$lambda))
  }

  # beta?
  if (is.null(MLIST$beta) || length(ov.y.dummy.lv.idx) == 0L ||
    length(nondummy.idx) == 0L) {
    LAMBDA..IB.inv <- LAMBDA
  } else {
    # only keep those columns of BETA that correspond to the
    # the `regular' latent variables
    # (ie. ignore the structural part altogether)
    MLIST2 <- MLIST
    MLIST2$beta[, dummy.idx] <- 0
    IB.inv <- lav_lisrel_ibinv(MLIST = MLIST2)
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv
  }

  # compute model-implied means
  EY <- lav_lisrel_ey(
    MLIST = MLIST, mean.x = mean.x,
    sample.mean = sample.mean,
    ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
    ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
    ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
    ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
  )

  EETA <- lav_lisrel_eeta(
    MLIST = MLIST, sample.mean = sample.mean,
    mean.x = mean.x,
    ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
    ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
    ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
    ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
  )

  # center regular lv only
  ETA[, nondummy.idx] <- sweep(ETA[, nondummy.idx, drop = FALSE], 2L,
    STATS = EETA[nondummy.idx], FUN = "-"
  )

  # project from lv to ov, if we have any lv
  if (length(nondummy.idx) > 0) {
    EYetax <- sweep(
      tcrossprod(
        ETA[, nondummy.idx, drop = FALSE],
        LAMBDA..IB.inv[, nondummy.idx, drop = FALSE]
      ),
      2L,
      STATS = EY, FUN = "+"
    )
  } else {
    EYetax <- ETA
  }

  # put back eXo variables
  if (length(ov.x.dummy.lv.idx) > 0L) {
    EYetax[, ov.x.dummy.ov.idx] <- ETA[, ov.x.dummy.lv.idx, drop = FALSE]
  }

  # if delta, scale
  if (delta && !is.null(MLIST$delta)) {
    EYetax <- sweep(EYetax, 2L, STATS = MLIST$delta, FUN = "*")
  }

  EYetax
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
lav_lisrel_vy <- function(MLIST = NULL) {
  LAMBDA <- MLIST$lambda
  THETA <- MLIST$theta

  VETA <- lav_lisrel_veta(MLIST = MLIST)
  VY <- tcrossprod(LAMBDA %*% VETA, LAMBDA) + THETA
  VY
}

# 5) VYx
# compute V(Y*|x_i) == model-implied covariance matrix
# this equals V(Y*) if no (explicit) eXo no GAMMA
#
# in >0.6-20: special treatment for composites
#
lav_lisrel_sigma <- function(MLIST = NULL, delta = TRUE) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  PSI <- MLIST$psi
  THETA <- MLIST$theta
  BETA <- MLIST$beta
  WMAT <- MLIST$wmat

  # standard: no composites
  if (is.null(WMAT)) {
    # beta?
    if (is.null(BETA)) {
      LAMBDA..IB.inv <- LAMBDA
    } else {
      IB.inv <- lav_lisrel_ibinv(MLIST = MLIST)
      LAMBDA..IB.inv <- LAMBDA %*% IB.inv
    }
    # compute V(Y*|x_i)
    VYx <- tcrossprod(LAMBDA..IB.inv %*% PSI, LAMBDA..IB.inv) + THETA

    # composites, or mix of composites and latent variables
  } else {
    # - first join LAMBDA and WMAT
    # - create 'T' matrix: - identity for regular lv's,
    #                      - THETA block-diagonal for composites
    # - create C_0: VETA, but zero diagonal elements for composites
    cov.idx <- which(apply(
      LAMBDA, 1L,
      function(x) sum(x == 0) == ncol(LAMBDA)
    ))
    clv.idx <- which(apply(
      LAMBDA, 2L,
      function(x) sum(x == 0) == nrow(LAMBDA)
    ))
    # regular latent variables
    rlv.idx <- seq_len(ncol(LAMBDA))[-clv.idx]

    # combine LAMBDA and WMAT
    LW <- LAMBDA + WMAT

    Tmat <- diag(nrow(LAMBDA))
    Tmat[cov.idx, cov.idx] <- THETA[cov.idx, cov.idx]
    wtw <- t(LW[, clv.idx, drop = FALSE]) %*% Tmat %*% LW[, clv.idx, drop = FALSE]
    wtw.inv <- solve(wtw)
    WTW.inv <- diag(ncol(LAMBDA))
    WTW.inv[clv.idx, clv.idx] <- wtw.inv

    if (is.null(BETA)) {
      IB.inv <- diag(nrow(PSI))
    } else {
      IB.inv <- lav_lisrel_ibinv(MLIST = MLIST)
    }
    VETA <- IB.inv %*% PSI %*% t(IB.inv)
    C0 <- VETA
    diag(C0)[clv.idx] <- 0

    VYx <- Tmat %*% LW %*% WTW.inv %*% C0 %*% t(WTW.inv) %*% t(LW) %*% Tmat + THETA
  }

  # if delta, scale
  if (delta && !is.null(MLIST$delta)) {
    DELTA <- diag(MLIST$delta[, 1L], nrow = nvar, ncol = nvar)
    VYx <- DELTA %*% VYx %*% DELTA
  }

  VYx
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
lav_lisrel_mu <- function(MLIST = NULL) {
  NU <- MLIST$nu
  ALPHA <- MLIST$alpha
  LAMBDA <- MLIST$lambda
  BETA <- MLIST$beta

  # shortcut
  if (is.null(ALPHA) || is.null(NU)) {
    return(matrix(0, nrow(LAMBDA), 1L))
  }

  # beta?
  if (is.null(BETA)) {
    LAMBDA..IB.inv <- LAMBDA
  } else {
    IB.inv <- lav_lisrel_ibinv(MLIST = MLIST)
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv
  }

  # compute Mu Hat
  Mu.hat <- NU + LAMBDA..IB.inv %*% ALPHA

  Mu.hat
}

# compute TH for a single block/group
lav_lisrel_th <- function(MLIST = NULL, th.idx = NULL, delta = TRUE) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)
  BETA <- MLIST$beta
  TAU <- MLIST$tau
  nth <- nrow(TAU)

  # missing alpha
  if (is.null(MLIST$alpha)) {
    ALPHA <- matrix(0, nfac, 1L)
  } else {
    ALPHA <- MLIST$alpha
  }

  # missing nu
  if (is.null(MLIST$nu)) {
    NU <- matrix(0, nvar, 1L)
  } else {
    NU <- MLIST$nu
  }

  if (is.null(th.idx)) {
    th.idx <- seq_len(nth)
    nlev <- rep(1L, nvar)
    K_nu <- diag(nvar)
  } else {
    nlev <- tabulate(th.idx, nbins = nvar)
    nlev[nlev == 0L] <- 1L
    K_nu <- matrix(0, sum(nlev), nvar)
    K_nu[cbind(seq_len(sum(nlev)), rep(seq_len(nvar), times = nlev))] <- 1.0
  }

  # shortcut
  if (is.null(TAU)) {
    return(matrix(0, length(th.idx), 1L))
  }

  # beta?
  if (is.null(BETA)) {
    LAMBDA..IB.inv <- LAMBDA
  } else {
    IB.inv <- lav_lisrel_ibinv(MLIST = MLIST)
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv
  }

  # compute pi0
  pi0 <- NU + LAMBDA..IB.inv %*% ALPHA

  # interleave th's with zeros where we have numeric variables
  th <- numeric(length(th.idx))
  th[th.idx > 0L] <- TAU[, 1L]

  # compute TH
  TH <- th - (K_nu %*% pi0)

  # if delta, scale
  if (delta && !is.null(MLIST$delta)) {
    DELTA.diag <- MLIST$delta[, 1L]
    DELTA.star.diag <- rep(DELTA.diag, times = nlev)
    TH <- TH * DELTA.star.diag
  }

  as.vector(TH)
}

# compute PI for a single block/group
lav_lisrel_pi <- function(MLIST = NULL, delta = TRUE) {
  LAMBDA <- MLIST$lambda
  BETA <- MLIST$beta
  GAMMA <- MLIST$gamma

  # shortcut
  if (is.null(GAMMA)) {
    return(matrix(0, nrow(LAMBDA), 0L))
  }

  # beta?
  if (is.null(BETA)) {
    LAMBDA..IB.inv <- LAMBDA
  } else {
    IB.inv <- lav_lisrel_ibinv(MLIST = MLIST)
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv
  }

  # compute PI
  PI <- LAMBDA..IB.inv %*% GAMMA

  # if delta, scale
  if (delta && !is.null(MLIST$delta)) {
    DELTA.diag <- MLIST$delta[, 1L]
    PI <- PI * DELTA.diag
  }

  PI
}

lav_lisrel_lambda <- function(MLIST = NULL,
                              ov.y.dummy.ov.idx = NULL,
                              ov.x.dummy.ov.idx = NULL,
                              ov.y.dummy.lv.idx = NULL,
                              ov.x.dummy.lv.idx = NULL,
                              remove.dummy.lv = FALSE) {
  lv.dummy.idx <- c(ov.y.dummy.lv.idx, ov.x.dummy.lv.idx)

  # fix LAMBDA
  LAMBDA <- MLIST$lambda
  if (length(ov.y.dummy.ov.idx) > 0L && !is.null(MLIST$beta)) {
    LAMBDA[ov.y.dummy.ov.idx, ] <- MLIST$beta[ov.y.dummy.lv.idx, ]
  }

  # remove dummy lv?
  if (remove.dummy.lv && length(lv.dummy.idx) > 0L) {
    LAMBDA <- LAMBDA[, -lv.dummy.idx, drop = FALSE]
  }

  LAMBDA
}

lav_lisrel_theta <- function(MLIST = NULL,
                             ov.y.dummy.ov.idx = NULL,
                             ov.x.dummy.ov.idx = NULL,
                             ov.y.dummy.lv.idx = NULL,
                             ov.x.dummy.lv.idx = NULL) {
  ov.dummy.idx <- c(ov.y.dummy.ov.idx, ov.x.dummy.ov.idx)
  lv.dummy.idx <- c(ov.y.dummy.lv.idx, ov.x.dummy.lv.idx)

  # fix THETA
  THETA <- MLIST$theta
  if (length(ov.dummy.idx) > 0L) {
    THETA[ov.dummy.idx, ov.dummy.idx] <-
      MLIST$psi[lv.dummy.idx, lv.dummy.idx]
  }

  THETA
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
lav_lisrel_ibinv <- function(MLIST = NULL) {
  BETA <- MLIST$beta
  nr   <- nrow(MLIST$psi)

  # case 1: no BETA, or BETA is identically zero
  if (is.null(BETA) || all(BETA == 0)) {
    return(diag(nr))
  }

  # forwardsolve/backsolve do not support complex values; solve() does
  if (is.complex(BETA)) {
    tmp <- -BETA; diag(tmp) <- 1; return(solve(tmp))
  }

  # case 2: strictly lower triangular
  if (all(BETA[upper.tri(BETA)] == 0)) {
    return(forwardsolve(diag(nr) - BETA, diag(nr)))
  }

  # case 3: strictly upper triangular
  if (all(BETA[lower.tri(BETA)] == 0)) {
    return(backsolve(diag(nr) - BETA, diag(nr)))
  }

  # case 4: recursive/DAG
  # BETA[i,j] != 0  <=>  directed edge j -> i
  indegree <- rowSums(BETA != 0)
  result <- integer(nr)
  k <- 0L
  queue <- which(indegree == 0)
  while (length(queue) > 0L) {
    v <- queue[1L]; queue <- queue[-1L]
    k <- k + 1L;   result[k] <- v
    for (u in which(BETA[, v] != 0)) {
      indegree[u] <- indegree[u] - 1L
      if (indegree[u] == 0L) queue <- c(queue, u)
    }
  }
  if (k == nr) {
    # recursive model: permute B to strictly lower triangular, solve, unpermute
    IB.inv <- forwardsolve(diag(nr) - BETA[result, result], diag(nr))
    inv.order <- order(result)
    return(IB.inv[inv.order, inv.order])
  }

  # case 5: non-recursive model
  tmp <- -BETA
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
lav_lisrel_alpha0 <- function(MLIST = NULL, sample.mean = NULL,
                              ov.y.dummy.ov.idx = NULL,
                              ov.x.dummy.ov.idx = NULL,
                              ov.y.dummy.lv.idx = NULL,
                              ov.x.dummy.lv.idx = NULL) {
  if (!is.null(MLIST$alpha)) {
    return(MLIST$alpha)
  }

  LAMBDA <- MLIST$lambda
  nfac <- ncol(LAMBDA)

  ov.dummy.idx <- c(ov.y.dummy.ov.idx, ov.x.dummy.ov.idx)
  lv.dummy.idx <- c(ov.y.dummy.lv.idx, ov.x.dummy.lv.idx)

  if (length(ov.dummy.idx) > 0L) {
    ALPHA <- matrix(0, nfac, 1L)
    # Note: instead of sample.mean, we need 'intercepts'
    # sample.mean = NU + LAMBDA..IB.inv %*% ALPHA
    # so,
    # solve(LAMBDA..IB.inv) %*% (sample.mean - NU) = ALPHA
    # where
    # - LAMBDA..IB.inv only contains 'dummy' variables, and is square
    # - NU elements are not needed (since not in ov.dummy.idx)
    IB.inv <- lav_lisrel_ibinv(MLIST = MLIST)
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv
    LAMBDA..IB.inv.dummy <- LAMBDA..IB.inv[ov.dummy.idx, lv.dummy.idx]
    ALPHA[lv.dummy.idx] <-
      solve(LAMBDA..IB.inv.dummy, sample.mean[ov.dummy.idx])
  } else {
    ALPHA <- matrix(0, nfac, 1L)
  }

  ALPHA
}

# only if NU=NULL but we need it anyway
#
# since we have no meanstructure, we can assume NU is unrestricted
# and contains either:
#     1) the sample means (if not eXo)
#     2) the intercepts, if we have exogenous covariates
#        since sample.mean = NU + LAMBDA %*% E(eta)
#        we have NU = sample.mean - LAMBDA %*% E(eta)
lav_lisrel_nu0 <- function(MLIST = NULL, sample.mean = NULL,
                           ov.y.dummy.ov.idx = NULL,
                           ov.x.dummy.ov.idx = NULL,
                           ov.y.dummy.lv.idx = NULL,
                           ov.x.dummy.lv.idx = NULL) {
  if (!is.null(MLIST$nu)) {
    return(MLIST$nu)
  }

  # if nexo > 0, subtract lambda %*% EETA
  if (length(ov.x.dummy.ov.idx) > 0L) {
    EETA <- lav_lisrel_eeta(MLIST,
      mean.x = NULL,
      sample.mean = sample.mean,
      ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
      ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
      ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
      ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
    )

    # 'regress' NU on X
    NU <- sample.mean - MLIST$lambda %*% EETA

    # just to make sure we have exact zeroes for all dummies
    NU[c(ov.y.dummy.ov.idx, ov.x.dummy.ov.idx)] <- 0
  } else {
    # unrestricted mean
    NU <- sample.mean
  }

  NU
}

# set (total/residual) variances of composites
# and while we at it, also set intercepts of composites
lav_lisrel_comp_set_intresvar <- function(MLIST = NULL,
                                          tol = .Machine$double.eps,
                                          debug = FALSE) {
  LAMBDA <- MLIST$lambda
  BETA <- MLIST$beta
  PSI <- MLIST$psi
  WMAT <- MLIST$wmat
  THETA <- MLIST$theta

  # std.lv or not?
  marker.idx <- lav_utils_get_marker(MLIST$wmat)
  std.lv <- FALSE
  if (all(is.na(marker.idx))) {
    std.lv <- TRUE
  }

  # housekeeping
  ovc.idx <- which(apply(
    LAMBDA, 1L,
    function(x) sum(x == 0) == ncol(LAMBDA)
  ))
  lvc.idx <- which(apply(
    LAMBDA, 2L,
    function(x) sum(x == 0) == nrow(LAMBDA)
  ))
  lvc.flag <- logical(nrow(PSI))
  lvc.flag[lvc.idx] <- TRUE
  Tmat <- diag(nrow(LAMBDA))
  Tmat[ovc.idx, ovc.idx] <- MLIST$theta[ovc.idx, ovc.idx]

  if (std.lv) {
    target.psi <- rep(1, ncol(WMAT))
  } else {
    # total variances composites
    target.psi <- diag(t(WMAT) %*% Tmat %*% WMAT)
  }
  # fill in PSI element for non-composites
  target.psi[!lvc.flag] <- diag(PSI)[!lvc.flag]

  # initial values (including exogenous variances)
  diag(PSI) <- target.psi

  # no regressions
  if (is.null(BETA)) {
    # store PSI
    MLIST$psi <- PSI

    # fix intercept (if needed)
    if (!is.null(MLIST$alpha)) {
      tmp <- t(WMAT) %*% MLIST$nu
      MLIST$alpha[lvc.idx, 1L] <- tmp[lvc.idx, 1L]
    }

    return(MLIST)
  }

  # set (residual) variances in psi
  abs.beta <- abs(MLIST$beta)
  x.idx <- which(apply(abs.beta, 1L, sum) == 0 & lvc.flag)
  y.idx <- which(apply(abs.beta, 1L, sum) != 0 & lvc.flag)

  # compute IB.inv
  IB.inv <- lav_lisrel_ibinv(MLIST)

  # check if BETA is acyclic
  if (!lav_graph_is_acyclic(BETA)) {
    # damn, we have a cyclic model; use nlminb()
    PSI <- lav_mlist_target_psi(
      IB.inv = IB.inv, PSI = PSI,
      target.psi = target.psi, y.idx = y.idx
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
    ancestors <- lav_graph_get_ancestors(BETA)

    # for each y variable, compute IB.inv %*% psi %*% t(IB.inv), without y
    ny <- length(y.idx)
    max_rep <- ny * 4
    for (rep in seq_len(max_rep)) {
      if (debug) {
        cat("rep = ", rep, "\n")
      }
      # check current diagonal
      current.diag <- diag(IB.inv %*% PSI %*% t(IB.inv))
      if (debug) {
        cat("target.psi   = ", target.psi[y.idx], "\n")
        cat("current.diag = ", current.diag[y.idx], "\n")
      }
      if (all(abs(current.diag[y.idx] - target.psi[y.idx]) < tol)) {
        # we are done, bail out
        break
      }
      for (i in seq_len(ny)) {
        this.y.idx <- y.idx[i]
        this.x.idx <- ancestors[[this.y.idx]]
        IB.inv.y <- IB.inv[this.y.idx, this.x.idx, drop = FALSE]
        PSI.x <- PSI[this.x.idx, this.x.idx, drop = FALSE]
        var.y <- drop(IB.inv.y %*% PSI.x %*% t(IB.inv.y))
        PSI[this.y.idx, this.y.idx] <- target.psi[this.y.idx] - var.y
      }
    }

    # final check?
    current.diag <- diag(IB.inv %*% PSI %*% t(IB.inv))
    if (debug) {
      cat("final current.diag = ", current.diag[y.idx], "\n")
    }
    # don't be too strict here
    if (any(abs(current.diag[y.idx] - target.psi[y.idx]) > sqrt(tol))) {
      # as a last resort, use optimization
      PSI <- lav_mlist_target_psi(
        IB.inv = IB.inv, PSI = PSI,
        target.psi = target.psi, y.idx = y.idx
      )
    }
  } # acyclic

  # store PSI
  MLIST$psi <- PSI

  # fix composite mean (if needed)
  if (!is.null(MLIST$alpha)) {
    tmp <- t(WMAT) %*% MLIST$nu # total means
    # step 1: set exogenous alpha values
    MLIST$alpha[lvc.idx, 1L] <- tmp[lvc.idx, 1L]
    # step 2: set endogenous alpha values
    exo_only <- MLIST$alpha; exo_only[y.idx] <- 0
    # compensate alpha so that total mean stays the same
    tmp2 <- IB.inv %*% exo_only
    MLIST$alpha[y.idx, 1L] <- MLIST$alpha[y.idx, 1L] - tmp2[y.idx, 1L]
  }

  MLIST
}

# if DELTA parameterization, compute residual elements (in theta, or psi)
# - typically for (endogenous) observed *categorical* variables only
# - but could also be all (endogenous) observed variables, if correlation = TRUE
# - or all (endogenous) latent and observed variables, if ov.only = FALSE
#
# new version YR 29 Oct 2024: try harder for psi elements (but this only works
#                             for acyclic models)
#             YR 01 Nov 2024: for non-acyclic models: use optimization
lav_lisrel_residual_variances <- function(MLIST = NULL,
                                          num.idx = NULL,
                                          ov.y.dummy.ov.idx = NULL,
                                          ov.y.dummy.lv.idx = NULL,
                                          ov.only = TRUE,
                                          tol = .Machine$double.eps,
                                          debug = FALSE) {
  BETA <- MLIST$beta
  PSI <- MLIST$psi
  if (is.null(MLIST$delta)) {
    delta <- rep(1, nrow(MLIST$lambda))
  } else {
    delta <- MLIST$delta
  }

  # remove num.idx from ov.y.dummy.*
  if (length(num.idx) > 0L && length(ov.y.dummy.ov.idx) > 0L) {
    n.idx <- which(ov.y.dummy.ov.idx %in% num.idx)
    if (length(n.idx) > 0L) {
      ov.y.dummy.ov.idx <- ov.y.dummy.ov.idx[-n.idx]
      ov.y.dummy.lv.idx <- ov.y.dummy.lv.idx[-n.idx]
    }
  }

  # if delta, the target may not be unity, but DELTA^(-2)
  target.all <- 1 / (delta * delta) # often the unit vector
  target.psi <- rep(1, nrow(PSI))
  if (length(ov.y.dummy.ov.idx) > 0L) {
    target.psi[ov.y.dummy.lv.idx] <- target.all[ov.y.dummy.ov.idx]
  }

  # phase 1: set (residual) variances in psi
  if (!is.null(BETA) && length(ov.y.dummy.ov.idx) > 0L) {
    abs.beta <- abs(MLIST$beta)
    x.idx <- which(apply(abs.beta, 1L, sum) == 0)
    if (ov.only) {
      y.idx <- ov.y.dummy.lv.idx
      # remove x.idx elements (fixed.x = FALSE)
      if (any(y.idx %in% x.idx)) {
        y.idx <- y.idx[-which(y.idx %in% x.idx)]
      }
    } else {
      y.idx <- which(apply(abs.beta, 1L, sum) != 0)
    }

    nr <- nrow(BETA)
    IB <- -BETA
    IB[lav_matrix_diag_idx(nr)] <- 1
    IB.inv <- solve(IB)


    # check if BETA is acyclic
    if (!lav_graph_is_acyclic(BETA)) {
      # damn, we have a cyclic model; use nlminb()
      PSI <- lav_mlist_target_psi(
        IB.inv = IB.inv, PSI = PSI,
        target.psi = target.psi, y.idx = y.idx
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
      ancestors <- lav_graph_get_ancestors(BETA)

      # for each y variable, compute IB.inv %*% psi %*% t(IB.inv), without y
      ny <- length(y.idx)
      max_rep <- ny * 4
      for (rep in seq_len(max_rep)) {
        if (debug) {
          cat("rep = ", rep, "\n")
        }
        # check current diagonal
        current.diag <- diag(IB.inv %*% PSI %*% t(IB.inv))
        if (debug) {
          cat("target.psi   = ", target.psi[y.idx], "\n")
          cat("current.diag = ", current.diag[y.idx], "\n")
        }
        if (all(abs(current.diag[y.idx] - target.psi[y.idx]) < tol)) {
          # we are done, bail out
          break
        }
        for (i in seq_len(ny)) {
          this.y.idx <- y.idx[i]
          this.x.idx <- ancestors[[this.y.idx]]
          IB.inv.y <- IB.inv[this.y.idx, this.x.idx, drop = FALSE]
          PSI.x <- PSI[this.x.idx, this.x.idx, drop = FALSE]
          var.y <- drop(IB.inv.y %*% PSI.x %*% t(IB.inv.y))
          PSI[this.y.idx, this.y.idx] <- target.psi[this.y.idx] - var.y
        }
      }

      # final check?
      current.diag <- diag(IB.inv %*% PSI %*% t(IB.inv))
      if (debug) {
        cat("final current.diag = ", current.diag[y.idx], "\n")
      }
      # don't be too strict here
      if (any(abs(current.diag[y.idx] - target.psi[y.idx]) > sqrt(tol))) {
        # as a last resort, use optimization
        PSI <- lav_mlist_target_psi(
          IB.inv = IB.inv, PSI = PSI,
          target.psi = target.psi, y.idx = y.idx
        )
      }
    } # acyclic
  } # phase 1

  # store PSI
  MLIST$psi <- PSI

  # phase 2: set residual variances in theta

  # force non-numeric theta elements to be zero
  if (length(num.idx) > 0L) {
    diag(MLIST$theta)[-num.idx] <- 0.0
  } else {
    diag(MLIST$theta) <- 0.0
  }

  Sigma.hat <- lav_lisrel_sigma(MLIST = MLIST, delta = FALSE)
  diag.Sigma <- diag(Sigma.hat)
  # theta = DELTA^(-2) - diag( LAMBDA (I-B)^-1 PSI (I-B)^-T t(LAMBDA) )
  theta.diag <- target.all - diag.Sigma
  not.idx <- unique(c(num.idx, ov.y.dummy.ov.idx))
  if (length(not.idx) > 0L) {
    diag(MLIST$theta)[-not.idx] <- theta.diag[-not.idx]
  } else {
    diag(MLIST$theta) <- theta.diag
  }

  MLIST
}

# if THETA parameterization, compute delta elements
# of observed categorical variables, as a function of other model parameters
lav_lisrel_delta <- function(MLIST = NULL, num.idx = NULL) {
  Sigma.hat <- lav_lisrel_sigma(MLIST = MLIST, delta = FALSE)
  diag.Sigma <- diag(Sigma.hat)

  # (1/delta^2) = diag( LAMBDA (I-B)^-1 PSI (I-B)^-T t(LAMBDA) ) + THETA
  # tmp <- diag.Sigma + THETA
  tmp <- diag.Sigma
  tmp[tmp < 0] <- as.numeric(NA)
  MLIST$delta[, 1L] <- sqrt(1 / tmp)

  # numeric delta's stay 1.0
  if (length(num.idx) > 0L) {
    MLIST$delta[num.idx] <- 1.0
  }

  MLIST
}

# compute Sigma/ETA: variances/covariances of BOTH observed and latent variables
lav_lisrel_cov_both <- function(MLIST = NULL, delta = TRUE) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  PSI <- MLIST$psi
  nlat <- nrow(PSI)
  THETA <- MLIST$theta
  BETA <- MLIST$beta

  # 'extend' matrices
  LAMBDA2 <- rbind(LAMBDA, diag(nlat))
  THETA2 <- lav_matrix_bdiag(THETA, matrix(0, nlat, nlat))


  # beta?
  if (is.null(BETA)) {
    LAMBDA..IB.inv <- LAMBDA2
  } else {
    IB.inv <- lav_lisrel_ibinv(MLIST = MLIST)
    LAMBDA..IB.inv <- LAMBDA2 %*% IB.inv
  }

  # compute augment COV matrix
  COV <- tcrossprod(LAMBDA..IB.inv %*% PSI, LAMBDA..IB.inv) + THETA2

  # if delta, scale
  if (delta && !is.null(MLIST$delta)) {
    DELTA <- diag(MLIST$delta[, 1L], nrow = nvar, ncol = nvar)
    COV[seq_len(nvar), seq_len(nvar)] <-
      DELTA %*% COV[seq_len(nvar), seq_len(nvar)] %*% DELTA
  }


  # if GAMMA, also x part
  GAMMA <- MLIST$gamma
  if (!is.null(GAMMA)) {
    COV.X <- MLIST$cov.x
    if (is.null(BETA)) {
      SX <- tcrossprod(GAMMA %*% COV.X, GAMMA)
    } else {
      IB.inv..GAMMA <- IB.inv %*% GAMMA
      SX <- tcrossprod(IB.inv..GAMMA %*% COV.X, IB.inv..GAMMA)
    }
    COV[(nvar + 1):(nvar + nlat), (nvar + 1):(nvar + nlat)] <-
      COV[(nvar + 1):(nvar + nlat), (nvar + 1):(nvar + nlat)] + SX
  }

  COV
}


# derivative of the objective function
lav_lisrel_df_dmlist <- function(MLIST = NULL, Omega = NULL, Omega.mu = NULL) {
  LAMBDA <- MLIST$lambda
  PSI <- MLIST$psi
  BETA <- MLIST$beta
  ALPHA <- MLIST$alpha
  WMAT <- MLIST$wmat

  LAMBDA.deriv <- NULL
  BETA.deriv <- NULL
  THETA.deriv <- NULL
  PSI.deriv <- NULL
  NU.deriv <- NULL
  ALPHA.deriv <- NULL
  GROUP.W.deriv <- NULL
  WMAT.deriv <- NULL

  # beta?
  if (is.null(BETA)) {
    LAMBDA..IB.inv <- LAMBDA
  } else {
    IB.inv <- lav_lisrel_ibinv(MLIST = MLIST)
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv
  }

  # meanstructure?
  meanstructure <- FALSE
  if (!is.null(Omega.mu)) meanstructure <- TRUE

  # group weight?
  group.w.free <- FALSE
  if (!is.null(MLIST$gw)) group.w.free <- TRUE

  # pre-compute some values
  tLAMBDA..IB.inv <- t(LAMBDA..IB.inv)
  if (!is.null(BETA)) {
    Omega..LAMBDA..IB.inv..PSI..tIB.inv <-
      (Omega %*% LAMBDA..IB.inv %*% PSI %*% t(IB.inv))
  } else {
    Omega..LAMBDA <- Omega %*% LAMBDA
  }

  # 1. LAMBDA
  if (!is.null(BETA)) {
    if (meanstructure) {
      LAMBDA.deriv <- -1.0 * (Omega.mu %*% t(ALPHA) %*% t(IB.inv) +
        Omega..LAMBDA..IB.inv..PSI..tIB.inv)
    } else {
      LAMBDA.deriv <- -1.0 * Omega..LAMBDA..IB.inv..PSI..tIB.inv
    }
  } else {
    # no BETA
    if (meanstructure) {
      LAMBDA.deriv <- -1.0 * (Omega.mu %*% t(ALPHA) +
        Omega..LAMBDA %*% PSI)
    } else {
      LAMBDA.deriv <- -1.0 * (Omega..LAMBDA %*% PSI)
    }
  }

  # 2. BETA
  if (!is.null(BETA)) {
    if (meanstructure) {
      BETA.deriv <- -1.0 * ((t(IB.inv) %*%
        (t(LAMBDA) %*% Omega.mu %*% t(ALPHA)) %*%
        t(IB.inv)) +
        (tLAMBDA..IB.inv %*%
          Omega..LAMBDA..IB.inv..PSI..tIB.inv))
    } else {
      BETA.deriv <- -1.0 * (tLAMBDA..IB.inv %*%
        Omega..LAMBDA..IB.inv..PSI..tIB.inv)
    }
  }

  # 3. PSI
  PSI.deriv <- -1.0 * (tLAMBDA..IB.inv %*% Omega %*% LAMBDA..IB.inv)
  diag(PSI.deriv) <- 0.5 * diag(PSI.deriv)

  # 4. THETA
  THETA.deriv <- -1.0 * Omega
  diag(THETA.deriv) <- 0.5 * diag(THETA.deriv)

  if (meanstructure) {
    # 5. NU
    NU.deriv <- -1.0 * Omega.mu

    # 6. ALPHA
    ALPHA.deriv <- -1.0 * t(t(Omega.mu) %*% LAMBDA..IB.inv)
  }

  if (group.w.free) {
    GROUP.W.deriv <- 0.0
  }

  list(
    lambda = LAMBDA.deriv,
    wmat = WMAT.deriv,
    beta = BETA.deriv,
    theta = THETA.deriv,
    psi = PSI.deriv,
    nu = NU.deriv,
    alpha = ALPHA.deriv,
    gw = GROUP.W.deriv
  )
}

# dSigma/dx -- per model matrix
lav_lisrel_dsigma_dx <- function(MLIST = NULL,
                                 m = "lambda",
                                 # all model matrix elements, or only a few?
                                 # NOTE: for symmetric matrices,
                                 # we assume that the have full size
                                 # (nvar*nvar) (but already correct for
                                 # symmetry)
                                 idx = seq_len(length(MLIST[[m]])),
                                 vech = TRUE,
                                 delta = TRUE) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)
  PSI <- MLIST$psi
  WMAT <- MLIST$wmat

  # for composites (vec version)
  compute.sigma <- function(x, mm = "wmat", MLIST = NULL) {
    mlist <- MLIST
    if (mm %in% c("psi", "theta")) {
      mlist[[mm]] <- lav_matrix_vech_reverse(x)
    } else {
      mlist[[mm]][, ] <- x
    }
    lav_matrix_vec(lav_lisrel_sigma(mlist))
  }

  composites <- FALSE
  if (!is.null(WMAT)) {
    composites <- TRUE
  }

  # only lower.tri part of sigma (not same order as elimination matrix?)
  v.idx <- lav_matrix_vech_idx(nvar)
  pstar <- nvar * (nvar + 1) / 2

  # shortcut for gamma, nu, alpha, tau,.... : empty matrix
  if (m == "nu" || m == "alpha" || m == "tau" || m == "gamma" ||
    m == "gw" || m == "cov.x" || m == "mean.x") {
    return(matrix(0.0, nrow = pstar, ncol = length(idx)))
  }

  # Delta?
  delta.flag <- FALSE
  if (delta && !is.null(MLIST$delta)) {
    DELTA <- MLIST$delta
    delta.flag <- TRUE
  } else if (m == "delta") { # modindices?
    return(matrix(0.0, nrow = pstar, ncol = length(idx)))
  }

  # beta?
  if (!is.null(MLIST$ibeta.inv)) {
    IB.inv <- MLIST$ibeta.inv
  } else {
    IB.inv <- lav_lisrel_ibinv(MLIST = MLIST)
  }

  # pre
  # if(m == "lambda" || m == "beta")
  #    IK <- diag(nvar*nvar) + lav_matrix_commutation(nvar, nvar)
  if (m == "lambda" || m == "beta") {
    L1 <- LAMBDA %*% IB.inv %*% PSI %*% t(IB.inv)
  }
  if (m == "beta" || m == "psi") {
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv
  }

  # here we go:
  if (m == "lambda") {
    KOL.idx <- matrix(1:(nvar * nfac), nvar, nfac, byrow = TRUE)[idx]
    DX <- (L1 %x% diag(nvar))[, idx, drop = FALSE] +
      (diag(nvar) %x% L1)[, KOL.idx, drop = FALSE]
  } else if (m == "beta") {
    if (composites) {
      DX <- lav_func_jacobian_complex(
        func = compute.sigma,
        x = lav_matrix_vec(MLIST$beta),
        mm = "beta", MLIST = MLIST
      )
      DX <- DX[, idx, drop = FALSE]
    } else {
      KOL.idx <- matrix(1:(nfac * nfac), nfac, nfac, byrow = TRUE)[idx]
      DX <- (L1 %x% LAMBDA..IB.inv)[, idx, drop = FALSE] +
        (LAMBDA..IB.inv %x% L1)[, KOL.idx, drop = FALSE]
      # this is not really needed (because we select idx=m.el.idx)
      # but just in case we need all elements of beta...
      DX[, which(idx %in% lav_matrix_diag_idx(nfac))] <- 0.0
    }
  } else if (m == "psi") {
    if (composites) {
      tmp <- lav_func_jacobian_complex(
        func = compute.sigma,
        x = lav_matrix_vech(MLIST$psi),
        mm = "psi", MLIST = MLIST
      )
      DX <- matrix(0, nrow = nrow(tmp), ncol = length(PSI))
      DX[, lav_matrix_vech_idx(nrow(PSI))] <- tmp
      DX[, lav_matrix_vechu_idx(nrow(PSI), diagonal = FALSE)] <-
        DX[, lav_matrix_vech_idx(nrow(PSI), diagonal = FALSE), drop = FALSE]
      DX <- DX[, idx, drop = FALSE]
    } else {
      DX <- (LAMBDA..IB.inv %x% LAMBDA..IB.inv)
      # symmetry correction, but keeping all duplicated elements
      # since we depend on idx=m.el.idx
      lower.idx <- lav_matrix_vech_idx(nfac, diagonal = FALSE)
      upper.idx <- lav_matrix_vechru_idx(nfac, diagonal = FALSE)
      offdiagSum <- DX[, lower.idx] + DX[, upper.idx]
      DX[, c(lower.idx, upper.idx)] <- cbind(offdiagSum, offdiagSum)
      DX <- DX[, idx, drop = FALSE]
    }
  } else if (m == "theta") {
    # DX <- diag(nvar*nvar) # very sparse...
    DX <- matrix(0, nvar * nvar, length(idx))
    DX[cbind(idx, seq_along(idx))] <- 1
    # symmetry correction not needed, since all off-diagonal elements
    # are zero?
  } else if (m == "delta") {
    Omega <- lav_lisrel_sigma(MLIST, delta = FALSE)
    DD <- diag(DELTA[, 1], nvar, nvar)
    DD.Omega <- (DD %*% Omega)
    A <- DD.Omega %x% diag(nvar)
    B <- diag(nvar) %x% DD.Omega
    DX <- A[, lav_matrix_diag_idx(nvar), drop = FALSE] +
      B[, lav_matrix_diag_idx(nvar), drop = FALSE]
    DX <- DX[, idx, drop = FALSE]
  } else if (m == "wmat") {
    # just a dummy to get us going
    DX <- lav_func_jacobian_complex(
      func = compute.sigma,
      x = lav_matrix_vec(WMAT),
      mm = "wmat", MLIST = MLIST
    )
    DX <- DX[, idx, drop = FALSE]

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

  if (delta.flag && !m == "delta") {
    DX <- DX * as.vector(DELTA %x% DELTA)
  }

  # vech?
  if (vech) {
    DX <- DX[v.idx, , drop = FALSE]
  }

  DX
}

# dMu/dx -- per model matrix
lav_lisrel_dmu_dx <- function(MLIST = NULL,
                              m = "alpha",
                              # all model matrix elements, or only a few?
                              idx = seq_len(length(MLIST[[m]]))) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)
  WMAT <- MLIST$wmat

  # shortcut for empty matrices
  if (m == "gamma" || m == "psi" || m == "theta" || m == "tau" ||
    m == "delta" || m == "gw" || m == "cov.x" || m == "mean.x") {
    return(matrix(0.0, nrow = nvar, ncol = length(idx)))
  }

  # missing alpha
  if (is.null(MLIST$alpha)) {
    ALPHA <- matrix(0, nfac, 1L)
  } else {
    ALPHA <- MLIST$alpha
  }


  # beta?
  if (!is.null(MLIST$ibeta.inv)) {
    IB.inv <- MLIST$ibeta.inv
  } else {
    IB.inv <- lav_lisrel_ibinv(MLIST = MLIST)
  }

  if (m == "nu") {
    DX <- diag(nvar)
  } else if (m == "lambda") {
    DX <- t(IB.inv %*% ALPHA) %x% diag(nvar)
  } else if (m == "wmat") {
    # dummy, just to get us going
    DX <- t(IB.inv %*% ALPHA) %x% diag(nvar)
  } else if (m == "beta") {
    DX <- t(IB.inv %*% ALPHA) %x% (LAMBDA %*% IB.inv)
    # this is not really needed (because we select idx=m.el.idx)
    DX[, lav_matrix_diag_idx(nfac)] <- 0.0
  } else if (m == "alpha") {
    DX <- LAMBDA %*% IB.inv
  } else {
    lav_msg_stop(gettext("wrong model matrix names:"), m)
  }

  DX <- DX[, idx, drop = FALSE]
  DX
}

# dTh/dx -- per model matrix
lav_lisrel_dth_dx <- function(MLIST = NULL,
                              m = "tau",
                              # all model matrix elements, or only a few?
                              idx = seq_len(length(MLIST[[m]])),
                              th.idx = NULL,
                              delta = TRUE) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)
  TAU <- MLIST$tau
  nth <- nrow(TAU)

  # missing alpha
  if (is.null(MLIST$alpha)) {
    ALPHA <- matrix(0, nfac, 1L)
  } else {
    ALPHA <- MLIST$alpha
  }

  # missing nu
  if (is.null(MLIST$nu)) {
    NU <- matrix(0, nvar, 1L)
  } else {
    NU <- MLIST$nu
  }

  # Delta?
  delta.flag <- FALSE
  if (delta && !is.null(MLIST$delta)) {
    DELTA <- MLIST$delta
    delta.flag <- TRUE
  }

  if (is.null(th.idx)) {
    th.idx <- seq_len(nth)
    nlev <- rep(1L, nvar)
    K_nu <- diag(nvar)
  } else {
    nlev <- tabulate(th.idx, nbins = nvar)
    nlev[nlev == 0L] <- 1L
    K_nu <- matrix(0, sum(nlev), nvar)
    K_nu[cbind(seq_len(sum(nlev)), rep(seq_len(nvar), times = nlev))] <- 1.0
  }

  # shortcut for empty matrices
  if (m == "gamma" || m == "psi" || m == "theta" || m == "gw" ||
    m == "cov.x" || m == "mean.x") {
    return(matrix(0.0, nrow = length(th.idx), ncol = length(idx)))
  }

  # beta?
  if (!is.null(MLIST$ibeta.inv)) {
    IB.inv <- MLIST$ibeta.inv
  } else {
    IB.inv <- lav_lisrel_ibinv(MLIST = MLIST)
  }

  if (m == "tau") {
    DX <- matrix(0, nrow = length(th.idx), ncol = nth)
    DX[th.idx > 0L, ] <- diag(nth)
    if (delta.flag) {
      DX <- DX * as.vector(K_nu %*% DELTA)
    }
  } else if (m == "nu") {
    DX <- (-1) * K_nu
    if (delta.flag) {
      DX <- DX * as.vector(K_nu %*% DELTA)
    }
  } else if (m == "lambda") {
    DX <- (-1) * t(IB.inv %*% ALPHA) %x% diag(nvar)
    DX <- K_nu %*% DX
    if (delta.flag) {
      DX <- DX * as.vector(K_nu %*% DELTA)
    }
  } else if (m == "beta") {
    DX <- (-1) * t(IB.inv %*% ALPHA) %x% (LAMBDA %*% IB.inv)
    # this is not really needed (because we select idx=m.el.idx)
    DX[, lav_matrix_diag_idx(nfac)] <- 0.0
    DX <- K_nu %*% DX
    if (delta.flag) {
      DX <- DX * as.vector(K_nu %*% DELTA)
    }
  } else if (m == "alpha") {
    DX <- (-1) * LAMBDA %*% IB.inv
    DX <- K_nu %*% DX
    if (delta.flag) {
      DX <- DX * as.vector(K_nu %*% DELTA)
    }
  } else if (m == "delta") {
    DX1 <- matrix(0, nrow = length(th.idx), ncol = 1)
    DX1[th.idx > 0L, ] <- TAU
    DX2 <- NU + LAMBDA %*% IB.inv %*% ALPHA
    DX2 <- K_nu %*% DX2
    DX <- K_nu * as.vector(DX1 - DX2)
  } else {
    lav_msg_stop(gettext("wrong model matrix names:"), m)
  }

  DX <- DX[, idx, drop = FALSE]
  DX
}

# dPi/dx -- per model matrix
lav_lisrel_dpi_dx <- function(MLIST = NULL,
                              m = "lambda",
                              # all model matrix elements, or only a few?
                              idx = seq_len(length(MLIST[[m]]))) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)
  GAMMA <- MLIST$gamma
  nexo <- ncol(GAMMA)

  # Delta?
  delta.flag <- FALSE
  if (!is.null(MLIST$delta)) {
    DELTA.diag <- MLIST$delta[, 1L]
    delta.flag <- TRUE
  }

  # shortcut for empty matrices
  if (m == "tau" || m == "nu" || m == "alpha" || m == "psi" ||
    m == "theta" || m == "gw" || m == "cov.x" || m == "mean.x") {
    return(matrix(0.0, nrow = nvar * nexo, ncol = length(idx)))
  }

  # beta?
  if (!is.null(MLIST$ibeta.inv)) {
    IB.inv <- MLIST$ibeta.inv
  } else {
    IB.inv <- lav_lisrel_ibinv(MLIST = MLIST)
  }

  if (m == "lambda") {
    DX <- t(IB.inv %*% GAMMA) %x% diag(nvar)
    if (delta.flag) {
      DX <- DX * DELTA.diag
    }
  } else if (m == "beta") {
    DX <- t(IB.inv %*% GAMMA) %x% (LAMBDA %*% IB.inv)
    # this is not really needed (because we select idx=m.el.idx)
    DX[, lav_matrix_diag_idx(nfac)] <- 0.0
    if (delta.flag) {
      DX <- DX * DELTA.diag
    }
  } else if (m == "gamma") {
    DX <- diag(nexo) %x% (LAMBDA %*% IB.inv)
    if (delta.flag) {
      DX <- DX * DELTA.diag
    }
  } else if (m == "delta") {
    PRE <- rep(1, nexo) %x% diag(nvar)
    DX <- PRE * as.vector(LAMBDA %*% IB.inv %*% GAMMA)
  } else {
    lav_msg_stop(gettext("wrong model matrix names:"), m)
  }

  DX <- DX[, idx, drop = FALSE]
  DX
}

# dGW/dx -- per model matrix
lav_lisrel_dgw_dx <- function(MLIST = NULL,
                              m = "gw",
                              # all model matrix elements, or only a few?
                              idx = seq_len(length(MLIST[[m]]))) {
  # shortcut for empty matrices
  if (m != "gw") {
    return(matrix(0.0, nrow = 1L, ncol = length(idx)))
  } else {
    # m == "gw"
    DX <- matrix(1.0, 1, 1)
  }

  DX <- DX[, idx, drop = FALSE]
  DX
}

# dlambda/dx -- per model matrix
lav_lisrel_dlambda_dx <- function(MLIST = NULL,
                                  m = "lambda",
                                  # all model matrix elements, or only a few?
                                  idx = seq_len(length(MLIST[[m]]))) {
  LAMBDA <- MLIST$lambda

  # shortcut for empty matrices
  if (m != "lambda") {
    return(matrix(0.0, nrow = length(LAMBDA), ncol = length(idx)))
  } else {
    # m == "lambda"
    DX <- diag(1, nrow = length(LAMBDA), ncol = length(LAMBDA))
  }

  DX <- DX[, idx, drop = FALSE]
  DX
}

# dpsi/dx -- per model matrix - FIXME!!!!!
lav_lisrel_dpsi_dx <- function(MLIST = NULL,
                               m = "psi",
                               # all model matrix elements, or only a few?
                               idx = seq_len(length(MLIST[[m]]))) {
  PSI <- MLIST$psi
  nfac <- nrow(PSI)
  v.idx <- lav_matrix_vech_idx(nfac)

  # shortcut for empty matrices
  if (m != "psi") {
    DX <- matrix(0.0, nrow = length(PSI), ncol = length(idx))
    return(DX[v.idx, , drop = FALSE])
  } else {
    # m == "psi"
    DX <- diag(1, nrow = length(PSI), ncol = length(PSI))
  }

  DX <- DX[v.idx, idx, drop = FALSE]
  DX
}

# dtheta/dx -- per model matrix
lav_lisrel_dtheta_dx <- function(MLIST = NULL,
                                 m = "theta",
                                 # all model matrix elements, or only a few?
                                 idx = seq_len(length(MLIST[[m]]))) {
  THETA <- MLIST$theta
  nvar <- nrow(THETA)
  v.idx <- lav_matrix_vech_idx(nvar)

  # shortcut for empty matrices
  if (m != "theta") {
    DX <- matrix(0.0, nrow = length(THETA), ncol = length(idx))
    return(DX[v.idx, , drop = FALSE])
  } else {
    # m == "theta"
    DX <- diag(1, nrow = length(THETA), ncol = length(THETA))
  }

  DX <- DX[v.idx, idx, drop = FALSE]
  DX
}


# dbeta/dx -- per model matrix
lav_lisrel_dbeta_dx <- function(MLIST = NULL,
                                m = "beta",
                                # all model matrix elements, or only a few?
                                idx = seq_len(length(MLIST[[m]]))) {
  BETA <- MLIST$beta

  # shortcut for empty matrices
  if (m != "beta") {
    return(matrix(0.0, nrow = length(BETA), ncol = length(idx)))
  } else {
    # m == "beta"
    DX <- diag(1, nrow = length(BETA), ncol = length(BETA))
  }

  DX <- DX[, idx, drop = FALSE]
  DX
}

# dgamma/dx -- per model matrix
lav_lisrel_dgamma_dx <- function(MLIST = NULL,
                                 m = "gamma",
                                 # all model matrix elements, or only a few?
                                 idx = seq_len(length(MLIST[[m]]))) {
  GAMMA <- MLIST$gamma

  # shortcut for empty matrices
  if (m != "gamma") {
    return(matrix(0.0, nrow = length(GAMMA), ncol = length(idx)))
  } else {
    # m == "gamma"
    DX <- diag(1, nrow = length(GAMMA), ncol = length(GAMMA))
  }

  DX <- DX[, idx, drop = FALSE]
  DX
}

# dnu/dx -- per model matrix
lav_lisrel_dnu_dx <- function(MLIST = NULL,
                              m = "nu",
                              # all model matrix elements, or only a few?
                              idx = seq_len(length(MLIST[[m]]))) {
  NU <- MLIST$nu

  # shortcut for empty matrices
  if (m != "nu") {
    return(matrix(0.0, nrow = length(NU), ncol = length(idx)))
  } else {
    # m == "nu"
    DX <- diag(1, nrow = length(NU), ncol = length(NU))
  }

  DX <- DX[, idx, drop = FALSE]
  DX
}

# dtau/dx -- per model matrix
lav_lisrel_dtau_dx <- function(MLIST = NULL,
                               m = "tau",
                               # all model matrix elements, or only a few?
                               idx = seq_len(length(MLIST[[m]]))) {
  TAU <- MLIST$tau

  # shortcut for empty matrices
  if (m != "tau") {
    return(matrix(0.0, nrow = length(TAU), ncol = length(idx)))
  } else {
    # m == "tau"
    DX <- diag(1, nrow = length(TAU), ncol = length(TAU))
  }

  DX <- DX[, idx, drop = FALSE]
  DX
}


# dalpha/dx -- per model matrix
lav_lisrel_dalpha_dx <- function(MLIST = NULL,
                                 m = "alpha",
                                 # all model matrix elements, or only a few?
                                 idx = seq_len(length(MLIST[[m]]))) {
  ALPHA <- MLIST$alpha

  # shortcut for empty matrices
  if (m != "alpha") {
    return(matrix(0.0, nrow = length(ALPHA), ncol = length(idx)))
  } else {
    # m == "alpha"
    DX <- diag(1, nrow = length(ALPHA), ncol = length(ALPHA))
  }

  DX <- DX[, idx, drop = FALSE]
  DX
}

# MLIST = NULL; meanstructure=TRUE; th=TRUE; delta=TRUE; pi=TRUE; gw=FALSE
# lav_matrix_vech_idx <- lavaan:::lav_matrix_vech_idx; lav_matrix_vechru_idx <- lavaan:::lav_matrix_vechru_idx
# vec <- lavaan:::vec; lav_func_jacobian_complex <- lavaan:::lav_func_jacobian_complex
# lav_lisrel_sigma <- lavaan:::lav_lisrel_sigma
# lav_lisrel_delta <- lavaan:::lav_lisrel_delta
lav_lisrel_test_derivatives <- function(MLIST = NULL,
                                        nvar = NULL, nfac = NULL, nexo = NULL,
                                        th.idx = NULL, num.idx = NULL,
                                        meanstructure = TRUE,
                                        th = TRUE, delta = TRUE, pi = TRUE,
                                        gw = FALSE, theta = FALSE) {
  if (is.null(MLIST)) {
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
    if (is.null(num.idx)) {
      num.idx <- sort(sample(seq_len(nvar), ceiling(nvar / 2)))
    }
    if (is.null(th.idx)) {
      th.idx <- integer(0L)
      for (i in seq_len(nvar)) {
        if (i %in% num.idx) {
          th.idx <- c(th.idx, 0)
        } else {
          th.idx <- c(th.idx, rep(i, sample(c(1, 1, 2, 6), 1L)))
        }
      }
    }
    nth <- sum(th.idx > 0L)

    MLIST <- list()
    MLIST$lambda <- matrix(0, nvar, nfac)
    MLIST$beta <- matrix(0, nfac, nfac)
    MLIST$theta <- matrix(0, nvar, nvar)
    MLIST$psi <- matrix(0, nfac, nfac)
    if (meanstructure) {
      MLIST$alpha <- matrix(0, nfac, 1L)
      MLIST$nu <- matrix(0, nvar, 1L)
    }
    if (th) MLIST$tau <- matrix(0, nth, 1L)
    if (delta) MLIST$delta <- matrix(0, nvar, 1L)
    MLIST$gamma <- matrix(0, nfac, nexo)
    if (gw) MLIST$gw <- matrix(0, 1L, 1L)

    # feed random numbers
    MLIST <- lapply(MLIST, function(x) {
      x[, ] <- rnorm(length(x))
      x
    })
    # fix
    diag(MLIST$beta) <- 0.0
    diag(MLIST$theta) <- diag(MLIST$theta) * diag(MLIST$theta) * 10
    diag(MLIST$psi) <- diag(MLIST$psi) * diag(MLIST$psi) * 10
    MLIST$psi[lav_matrix_vechru_idx(nfac)] <-
      MLIST$psi[lav_matrix_vech_idx(nfac)]
    MLIST$theta[lav_matrix_vechru_idx(nvar)] <-
      MLIST$theta[lav_matrix_vech_idx(nvar)]
    if (delta) MLIST$delta[, ] <- abs(MLIST$delta) * 10
  } else {
    nvar <- nrow(MLIST$lambda)
  }

  compute.sigma <- function(x, mm = "lambda", MLIST = NULL) {
    mlist <- MLIST
    if (mm %in% c("psi", "theta")) {
      mlist[[mm]] <- lav_matrix_vech_reverse(x)
    } else {
      mlist[[mm]][, ] <- x
    }
    if (theta) {
      mlist <- lav_lisrel_delta(MLIST = mlist, num.idx = num.idx)
    }
    lav_matrix_vech(lav_lisrel_sigma(mlist))
  }

  compute.mu <- function(x, mm = "lambda", MLIST = NULL) {
    mlist <- MLIST
    if (mm %in% c("psi", "theta")) {
      mlist[[mm]] <- lav_matrix_vech_reverse(x)
    } else {
      mlist[[mm]][, ] <- x
    }
    if (theta) {
      mlist <- lav_lisrel_delta(MLIST = mlist, num.idx = num.idx)
    }
    lav_lisrel_mu(mlist)
  }

  compute.th2 <- function(x, mm = "tau", MLIST = NULL, th.idx) {
    mlist <- MLIST
    if (mm %in% c("psi", "theta")) {
      mlist[[mm]] <- lav_matrix_vech_reverse(x)
    } else {
      mlist[[mm]][, ] <- x
    }
    if (theta) {
      mlist <- lav_lisrel_delta(MLIST = mlist, num.idx = num.idx)
    }
    lav_lisrel_th(mlist, th.idx = th.idx)
  }

  compute.pi <- function(x, mm = "lambda", MLIST = NULL) {
    mlist <- MLIST
    if (mm %in% c("psi", "theta")) {
      mlist[[mm]] <- lav_matrix_vech_reverse(x)
    } else {
      mlist[[mm]][, ] <- x
    }
    if (theta) {
      mlist <- lav_lisrel_delta(MLIST = mlist, num.idx = num.idx)
    }
    lav_lisrel_pi(mlist)
  }

  compute.gw <- function(x, mm = "gw", MLIST = NULL) {
    mlist <- MLIST
    if (mm %in% c("psi", "theta")) {
      mlist[[mm]] <- lav_matrix_vech_reverse(x)
    } else {
      mlist[[mm]][, ] <- x
    }
    if (theta) {
      mlist <- lav_lisrel_delta(MLIST = mlist, num.idx = num.idx)
    }
    mlist$gw[1, 1]
  }

  # if theta, set MLIST$delta
  if (theta) {
    MLIST <- lav_lisrel_delta(MLIST = MLIST, num.idx = num.idx)
  }

  for (mm in names(MLIST)) {
    if (mm %in% c("psi", "theta")) {
      x <- lav_matrix_vech(MLIST[[mm]])
    } else {
      x <- lav_matrix_vec(MLIST[[mm]])
    }
    if (mm == "delta" && theta) next
    if (lav_debug()) {
      cat("### mm = ", mm, "\n")
    }

    # 1. sigma
    DX1 <- lav_func_jacobian_complex(func = compute.sigma, x = x, mm = mm, MLIST = MLIST)
    DX2 <- lav_lisrel_dsigma_dx(
      MLIST = MLIST, m = mm, idx = seq_len(length(MLIST[[mm]])),
      delta = !theta
    )
    if (mm %in% c("psi", "theta")) {
      # remove duplicated columns of symmetric matrices
      idx <- lav_matrix_vechru_idx(sqrt(ncol(DX2)), diagonal = FALSE)
      if (length(idx) > 0L) DX2 <- DX2[, -idx]
    }
    if (theta) {
      sigma.hat <- lav_lisrel_sigma(MLIST = MLIST, delta = FALSE)
      R <- lav_deriv_cov2cor(sigma.hat, num.idx = num.idx)

      DX3 <- DX2
      DX2 <- R %*% DX2
    }
    if (lav_debug()) {
      cat("[SIGMA] mm = ", sprintf("%-8s:", mm), "DX1 (numerical):\n")
      print(zapsmall(DX1))
      cat("\n")
      cat("[SIGMA] mm = ", sprintf("%-8s:", mm), "DX2 (analytical):\n")
      print(DX2)
      cat("\n")
      if (theta) {
        cat("[SIGMA] mm = ", sprintf("%-8s:", mm), "DX3 (analytical):\n")
        print(DX3)
        cat("\n")
      }
    }
    cat(
      "[SIGMA] mm = ", sprintf("%-8s:", mm), "sum delta = ",
      sprintf("%12.9f", sum(DX1 - DX2)), "  max delta = ",
      sprintf("%12.9f", max(DX1 - DX2)), "\n"
    )

    # 2. mu
    DX1 <- lav_func_jacobian_complex(func = compute.mu, x = x, mm = mm, MLIST = MLIST)
    DX2 <- lav_lisrel_dmu_dx(
      MLIST = MLIST,
      m = mm, idx = seq_len(length(MLIST[[mm]]))
    )
    if (mm %in% c("psi", "theta")) {
      # remove duplicated columns of symmetric matrices
      idx <- lav_matrix_vechru_idx(sqrt(ncol(DX2)), diagonal = FALSE)
      if (length(idx) > 0L) DX2 <- DX2[, -idx]
    }
    cat(
      "[MU   ] mm = ", sprintf("%-8s:", mm), "sum delta = ",
      sprintf("%12.9f", sum(DX1 - DX2)), "  max delta = ",
      sprintf("%12.9f", max(DX1 - DX2)), "\n"
    )
    if (lav_debug()) {
      cat("[MU   ] mm = ", sprintf("%-8s:", mm), "DX1 (numerical):\n")
      print(zapsmall(DX1))
      cat("\n")
      cat("[MU   ] mm = ", sprintf("%-8s:", mm), "DX2 (analytical):\n")
      print(DX2)
      cat("\n")
    }

    # 3. th
    if (th) {
      DX1 <- lav_func_jacobian_complex(
        func = compute.th2, x = x, mm = mm, MLIST = MLIST,
        th.idx = th.idx
      )
      DX2 <- lav_lisrel_dth_dx(
        MLIST = MLIST, m = mm, idx = seq_len(length(MLIST[[mm]])),
        th.idx = th.idx,
        delta = TRUE
      )
      if (theta) {
        # 1. compute dDelta.dx
        dxSigma <-
          lav_lisrel_dsigma_dx(
            m = mm, idx = seq_len(length(MLIST[[mm]])),
            MLIST = MLIST, delta = !theta
          )
        var.idx <- which(!lav_matrix_vech_idx(nvar) %in%
          lav_matrix_vech_idx(nvar, diagonal = FALSE))
        sigma.hat <- lav_lisrel_sigma(MLIST = MLIST, delta = FALSE)
        dsigma <- diag(sigma.hat)
        # dy/ddsigma = -0.5/(ddsigma*sqrt(ddsigma))
        dDelta.dx <- dxSigma[var.idx, ] * -0.5 / (dsigma * sqrt(dsigma))

        # 2. compute dth.dDelta
        dth.dDelta <-
          lav_lisrel_dth_dx(
            MLIST = MLIST,
            m = "delta",
            idx = seq_len(length(MLIST[["delta"]])),
            th.idx = th.idx
          )

        # 3. add dth.dDelta %*% dDelta.dx
        no.num.idx <- which(th.idx > 0)
        DX2[no.num.idx, ] <- DX2[no.num.idx, , drop = FALSE] +
          (dth.dDelta %*% dDelta.dx)[no.num.idx, , drop = FALSE]
        # DX2 <- DX2 + dth.dDelta %*% dDelta.dx
      }
      if (mm %in% c("psi", "theta")) {
        # remove duplicated columns of symmetric matrices
        idx <- lav_matrix_vechru_idx(sqrt(ncol(DX2)), diagonal = FALSE)
        if (length(idx) > 0L) DX2 <- DX2[, -idx]
      }
      cat(
        "[TH   ] mm = ", sprintf("%-8s:", mm), "sum delta = ",
        sprintf("%12.9f", sum(DX1 - DX2)), "  max delta = ",
        sprintf("%12.9f", max(DX1 - DX2)), "\n"
      )
      if (lav_debug()) {
        cat("[TH   ] mm = ", sprintf("%-8s:", mm), "DX1 (numerical):\n")
        print(zapsmall(DX1))
        cat("\n")
        cat("[TH   ] mm = ", sprintf("%-8s:", mm), "DX2 (analytical):\n")
        print(DX2)
        cat("\n")
      }
    }

    # 4. pi
    if (pi) {
      DX1 <- lav_func_jacobian_complex(func = compute.pi, x = x, mm = mm, MLIST = MLIST)
      DX2 <- lav_lisrel_dpi_dx(
        MLIST = MLIST,
        m = mm, idx = seq_len(length(MLIST[[mm]]))
      )
      if (mm %in% c("psi", "theta")) {
        # remove duplicated columns of symmetric matrices
        idx <- lav_matrix_vechru_idx(sqrt(ncol(DX2)), diagonal = FALSE)
        if (length(idx) > 0L) DX2 <- DX2[, -idx]
      }
      if (theta) {
        # 1. compute dDelta.dx
        dxSigma <-
          lav_lisrel_dsigma_dx(
            MLIST = MLIST, m = mm, idx = seq_len(length(MLIST[[mm]])),
            delta = !theta
          )
        if (mm %in% c("psi", "theta")) {
          # remove duplicated columns of symmetric matrices
          idx <- lav_matrix_vechru_idx(sqrt(ncol(dxSigma)), diagonal = FALSE)
          if (length(idx) > 0L) dxSigma <- dxSigma[, -idx]
        }
        var.idx <- which(!lav_matrix_vech_idx(nvar) %in%
          lav_matrix_vech_idx(nvar, diagonal = FALSE))
        sigma.hat <- lav_lisrel_sigma(MLIST = MLIST, delta = FALSE)
        dsigma <- diag(sigma.hat)
        # dy/ddsigma = -0.5/(ddsigma*sqrt(ddsigma))
        dDelta.dx <- dxSigma[var.idx, ] * -0.5 / (dsigma * sqrt(dsigma))

        # 2. compute dpi.dDelta
        dpi.dDelta <-
          lav_lisrel_dpi_dx(
            MLIST = MLIST,
            m = "delta",
            idx = seq_len(length(MLIST[["delta"]]))
          )

        # 3. add dpi.dDelta %*% dDelta.dx
        no.num.idx <- which(!seq.int(1L, nvar) %in% num.idx)
        no.num.idx <- rep(seq.int(0, nexo - 1) * nvar,
          each = length(no.num.idx)
        ) + no.num.idx
        DX2[no.num.idx, ] <- DX2[no.num.idx, , drop = FALSE] +
          (dpi.dDelta %*% dDelta.dx)[no.num.idx, , drop = FALSE]
      }
      cat(
        "[PI   ] mm = ", sprintf("%-8s:", mm), "sum delta = ",
        sprintf("%12.9f", sum(DX1 - DX2)), "  max delta = ",
        sprintf("%12.9f", max(DX1 - DX2)), "\n"
      )
      if (lav_debug()) {
        cat("[PI   ] mm = ", sprintf("%-8s:", mm), "DX1 (numerical):\n")
        print(zapsmall(DX1))
        cat("\n")
        cat("[PI   ] mm = ", sprintf("%-8s:", mm), "DX2 (analytical):\n")
        print(DX2)
        cat("\n")
      }
    }

    # 5. gw
    if (gw) {
      DX1 <- lav_func_jacobian_complex(func = compute.gw, x = x, mm = mm, MLIST = MLIST)
      DX2 <- lav_lisrel_dgw_dx(
        MLIST = MLIST,
        m = mm, idx = seq_len(length(MLIST[[mm]]))
      )
      if (mm %in% c("psi", "theta")) {
        # remove duplicated columns of symmetric matrices
        idx <- lav_matrix_vechru_idx(sqrt(ncol(DX2)), diagonal = FALSE)
        if (length(idx) > 0L) DX2 <- DX2[, -idx]
      }
      cat(
        "[GW   ] mm = ", sprintf("%-8s:", mm), "sum delta = ",
        sprintf("%12.9f", sum(DX1 - DX2)), "  max delta = ",
        sprintf("%12.9f", max(DX1 - DX2)), "\n"
      )
      if (lav_debug()) {
        cat("[GW   ] mm = ", sprintf("%-8s:", mm), "DX1 (numerical):\n")
        print(DX1)
        cat("\n\n")
        cat("[GW   ] mm = ", sprintf("%-8s:", mm), "DX2 (analytical):\n")
        print(DX2)
        cat("\n\n")
      }
    }
  }

  MLIST$th.idx <- th.idx
  MLIST$num.idx <- num.idx

  MLIST
}

# check for marker indicators:
#   - if std.lv = FALSE: a single '1' per factor, everything else zero
#   - if std.lv = TRUE: a single non-zero value per factor, everything else zero
lav_utils_get_marker <- function(LAMBDA = NULL, std.lv = FALSE) {
  LAMBDA <- as.matrix(LAMBDA)
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)

  # round values
  LAMBDA <- round(LAMBDA, 3L)

  marker.idx <- numeric(nfac)
  for (f in seq_len(nfac)) {
    if (std.lv) {
      marker.idx[f] <- which(rowSums(cbind(
        LAMBDA[, f] != 0,
        LAMBDA[, -f] == 0
      )) == nfac)[1]
    } else {
      marker.idx[f] <- which(rowSums(cbind(
        LAMBDA[, f] == 1,
        LAMBDA[, -f] == 0
      )) == nfac)[1]
    }
  }

  marker.idx
}

# find the residual variances (diagonal element of PSI) in such a way
# so that the diagonal elements of IB.inv %*% PSI %*% t(IB.inv) are
# equal to the elements of target.psi (usually the 1 vector)
#
# YR 01 Nov 2024: initial version; no bounds for now... (so we may end up
#                 with negative variances)
lav_mlist_target_psi <- function(BETA = NULL, IB.inv = NULL, PSI = NULL,
                                 target.psi = NULL, y.idx = NULL) {
  nr <- nrow(PSI)

  # IB.inv (if not given)
  if (is.null(IB.inv)) {
    IB <- -BETA
    IB[lav_matrix_diag_idx(nr)] <- 1
    IB.inv <- solve(IB)
  }

  # target.psi
  if (is.null(target.psi)) {
    target.psi <- rep(1, nr)
  }

  # y.idx
  if (is.null(y.idx)) {
    y.idx <- seq_len(nr)
  }

  # cast the problem as a nonlinear optimization problem
  obj <- function(x) {
    # cat("x = ", x, "\n")
    # x are the diagonal elements of PSI
    this.PSI <- PSI
    diag(this.PSI)[y.idx] <- x
    VETA <- IB.inv %*% this.PSI %*% t(IB.inv)
    current.diag <- diag(VETA)
    # ratio or difference?
    diff <- target.psi[y.idx] - current.diag[y.idx]
    out <- sum(diff * diff) # least squares
    out
  }

  VETA <- IB.inv %*% PSI %*% t(IB.inv)
  x.start <- diag(VETA)[y.idx]
  out <- nlminb(start = x.start, objective = obj)

  # return updated PSI matrix
  diag(PSI)[y.idx] <- out$par
  PSI
}

# second-order derivatives of Sigma wrt elements of Lambda/Beta/Psi/Theta
# (elementwise)
#
# work in progress
lav_lisrel_d2sigma_lambda_psi <- function(MLIST = NULL,
                                          i = 1L, j = 1L, # lambda
                                          k = 1L, l = 1L) { # psi
  # helper functions
  # elementary matrix: p x q matrix with 1 in position (i,j), 0 elsewhere
  elem_mat <- function(i, j, p, q) {
    E <- matrix(0, p, q)
    E[i, j] <- 1
    E
  }

  # elementary symmetric matrix for vech index (k,l) of a q x q symmetric matrix
  elem_mat_sym <- function(k, l, q) {
    E <- matrix(0, q, q)
    E[k, l] <- 1
    if (k != l) E[l, k] <- 1
    E
  }

  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA) # p
  nfac <- ncol(LAMBDA) # q
  IB.inv <- MLIST$IB.inv
  if (is.null(IB.inv)) {
    if (is.null(MLIST$beta)) {
      IB.inv <- diag(1, nrow = nfac)
    } else {
      IB.inv <- lav_lisrel_ibinv(MLIST = MLIST)
    }
  }

  LIB.inv <- LAMBDA %*% IB.inv

  dLambda <- elem_mat(i, j, nvar, nfac) # p x q,  elementary matrix for Lambda
  dPsi <- elem_mat_sym(k, l, nfac) # q x q,  elementary symm matrix for Psi

  X <- dLambda %*% IB.inv %*% dPsi %*% t(LIB.inv) # p x p
  X + t(X) # p x p symmetric matrix
}

lav_lisrel_d2sigma_beta_psi <- function(MLIST = NULL,
                                        i = 1L, j = 1L, # beta
                                        k = 1L, l = 1L) { # psi
  # helper functions
  # elementary matrix: p x q matrix with 1 in position (i,j), 0 elsewhere
  elem_mat <- function(i, j, p, q) {
    E <- matrix(0, p, q)
    E[i, j] <- 1
    E
  }

  # elementary symmetric matrix for vech index (k,l) of a q x q symmetric matrix
  elem_mat_sym <- function(k, l, q) {
    E <- matrix(0, q, q)
    E[k, l] <- 1
    if (k != l) E[l, k] <- 1
    E
  }

  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA) # p
  nfac <- ncol(LAMBDA) # q
  IB.inv <- MLIST$IB.inv
  if (is.null(IB.inv)) {
    if (is.null(MLIST$beta)) {
      IB.inv <- diag(1, nrow = nfac)
    } else {
      IB.inv <- lav_lisrel_ibinv(MLIST = MLIST)
    }
  }
  LIB.inv <- LAMBDA %*% IB.inv

  dBeta <- elem_mat(i, j, nfac, nfac) # q x q, elementary matrix for Beta
  dPsi <- elem_mat_sym(k, l, nfac) # q x q, elementary symm matrix for Psi

  X <- LIB.inv %*% dBeta %*% IB.inv %*% dPsi %*% t(LIB.inv) # p x p
  X + t(X) # p x p symmetric matrix
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
lav_lisrel_dimplied_dx <- function(MLIST           = NULL,
                                   m.free.idx      = NULL,
                                   x.free.idx      = NULL,
                                   nx.free         = NULL,
                                   meanstructure   = FALSE,
                                   categorical     = FALSE,
                                   correlation     = FALSE,
                                   conditional.x   = FALSE,
                                   num.idx         = integer(0L),
                                   th.idx          = integer(0L),
                                   group.w.free    = FALSE,
                                   parameterization = "delta") {

  delta.flag <- beta.flag <- FALSE
  if (!is.null(MLIST$delta) && any(MLIST$delta[,1] != 1)) {
    delta.flag <- TRUE
  }
  if (!is.null(MLIST$beta)) {
    beta.flag <- TRUE
  }

  # model matrices in this block
  mnames <- names(MLIST)

  mm.lambda.idx <- which(mnames == "lambda")
  x.lambda.idx <- x.free.idx[[mm.lambda.idx]]
  m.lambda.idx <- m.free.idx[[mm.lambda.idx]]
  n_lam <- length(m.lambda.idx)

  mm.psi.idx <- which(mnames == "psi")
  tmp <- x.free.idx[[mm.psi.idx]]
  x.psi.idx <- x.psi.idx <- tmp[!duplicated(tmp)]
  m.psi.idx <- m.free.idx[[mm.psi.idx]][!duplicated(tmp)]
  n_psi <- length(m.psi.idx)


  mm.theta.idx <- which(mnames == "theta")
  tmp <- x.free.idx[[mm.theta.idx]]
  x.theta.idx <- tmp[!duplicated(tmp)]
  m.theta.idx <- m.free.idx[[mm.theta.idx]][!duplicated(tmp)]
  n_the <- length(m.theta.idx)

  n_bet <- 0L
  x.beta.idx <- integer(0L)
  if (beta.flag) {
    mm.beta.idx <- which(mnames == "beta")
    x.beta.idx <- x.free.idx[[mm.beta.idx]]
    m.beta.idx <- m.free.idx[[mm.beta.idx]]
    n_bet <- length(m.beta.idx)
  }

  n_del <- 0L
  x.delta.idx <- integer(0L)
  if (delta.flag) {
    mm.delta.idx <- which(mnames == "delta")
    x.delta.idx <- x.free.idx[[mm.delta.idx]]
    m.delta.idx <- m.free.idx[[mm.delta.idx]]
    n_del <- length(m.delta.idx)
  }

  # composite: wmat
  n_wmat <- 0L
  x.wmat.idx <- integer(0L)
  m.wmat.idx <- integer(0L)
  wmat.flag <- !is.null(MLIST$wmat)
  if (wmat.flag) {
    mm.wmat.idx <- which(mnames == "wmat")
    if (length(mm.wmat.idx) > 0L) {
      x.wmat.idx <- x.free.idx[[mm.wmat.idx]]
      m.wmat.idx <- m.free.idx[[mm.wmat.idx]]
      n_wmat <- length(m.wmat.idx)
    }
  }

  # in case !meanstructure or !categorical
  n_nu  <- 0L;  m.nu.idx    <- integer(0L);  x.nu.idx    <- integer(0L)
  n_alp <- 0L;  m.alpha.idx <- integer(0L);  x.alpha.idx <- integer(0L)
  n_gam <- 0L;  m.gamma.idx <- integer(0L);  x.gamma.idx <- integer(0L)


  if (meanstructure) {
    mm.nu.idx <- which(mnames == "nu")
    x.nu.idx <- x.free.idx[[mm.nu.idx]]
    m.nu.idx <- m.free.idx[[mm.nu.idx]]
    n_nu <- length(m.nu.idx)

    mm.alpha.idx <- which(mnames == "alpha")
    x.alpha.idx <- x.free.idx[[mm.alpha.idx]]
    m.alpha.idx <- m.free.idx[[mm.alpha.idx]]
    n_alp <- length(m.alpha.idx)
  }

  if (conditional.x) {
    mm.gamma.idx <- which(mnames == "gamma")
    if (length(mm.gamma.idx) > 0L) {
      x.gamma.idx <- x.free.idx[[mm.gamma.idx]]
      m.gamma.idx <- m.free.idx[[mm.gamma.idx]]
      n_gam <- length(m.gamma.idx)
    }
  }

  if (categorical) {
    mm.th.idx <- which(mnames == "tau")
    x.th.idx <- x.free.idx[[mm.th.idx]]
    m.th.idx <- m.free.idx[[mm.th.idx]]
    n_th <- length(m.th.idx)
  }

  nvar <- nrow(MLIST$lambda)
  nfac <- ncol(MLIST$lambda)
  pstar <- nvar * (nvar + 1L) / 2L

  # precompute
  if (beta.flag) {
    A <- lav_lisrel_ibinv(MLIST)
    M <- MLIST$lambda %*% A
    G <- A %*% MLIST$psi %*% t(A)
    LG <- MLIST$lambda %*% G
  } else {
    M <- MLIST$lambda
    G <- MLIST$psi
    LG <- MLIST$lambda %*% MLIST$psi
  }

  if (delta.flag) {
    delta  <- as.vector(MLIST$delta)
  }


  # vech structure for Sigma
  r_s <- lav_matrix_vech_row_idx(nvar)
  c_s <- lav_matrix_vech_col_idx(nvar)
  sigma_lut <- lav_matrix_vech_reverse(seq_len(pstar))

  if (delta.flag) {
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


  if (!wmat.flag) {
    # Lambda [i,j]: dS0[r,s] = LG[s,j]*I(r==i) + LG[r,j]*I(s==i)
    if (n_lam > 0L) {
      rc  <- vec2rc(m.lambda.idx, nvar);  i_v <- rc[, 1L];  j_v <- rc[, 2L]

      T1 <- LG[c_s, j_v, drop = FALSE] * outer(r_s, i_v, `==`)
      T2 <- LG[r_s, j_v, drop = FALSE] * outer(c_s, i_v, `==`)

      DX <- T1 + T2
      if (delta.flag) {
        DX <- DX * delta_weight
      }
      jac_sigma[, col:(col + n_lam - 1L)] <- DX
      col <- col + n_lam
    }

    # Beta [i,j]: dS0[r,s] = M[r,i]*LG[s,j] + LG[r,j]*M[s,i]
    if (n_bet > 0L) {
      rc  <- vec2rc(m.beta.idx, nfac);  i_v <- rc[, 1L];  j_v <- rc[, 2L]

      T1 <- M[r_s,  i_v, drop = FALSE] * LG[c_s, j_v, drop = FALSE]
      T2 <- LG[r_s, j_v, drop = FALSE] * M[c_s,  i_v, drop = FALSE]

      DX <- T1 + T2
      if (delta.flag) {
        DX <- DX * delta_weight
      }
      jac_sigma[, col:(col + n_bet - 1L)] <- DX
      col <- col + n_bet
    }
   # Psi [k,l] symmetric: dS0[r,s] = M[r,k]*M[s,l] + M[r,l]*M[s,k]
    # diagonal (k==l): halve
    if (n_psi > 0L) {
      rc  <- vec2rc(m.psi.idx, nfac)
      k_v <- pmax(rc[, 1L], rc[, 2L])
      l_v <- pmin(rc[, 1L], rc[, 2L])

      T1 <- M[r_s, k_v, drop = FALSE] * M[c_s, l_v, drop = FALSE]
      T2 <- M[r_s, l_v, drop = FALSE] * M[c_s, k_v, drop = FALSE]
      DX <- T1 + T2

      diag_mask <- (k_v == l_v)
      if (any(diag_mask)) {
        DX[, diag_mask] <- DX[, diag_mask] * 0.5
      }

      if (delta.flag) {
        DX <- DX * delta_weight
      }
      jac_sigma[, col:(col + n_psi - 1L)] <- DX
      col <- col + n_psi
    }

    # Theta [k,l] symmetric: dS = Delta %*% dTheta %*% Delta
    #  -> unit vector at vech pos (k,l)
    if (n_the > 0L) {
      rc  <- vec2rc(m.theta.idx, nvar)
      k_v <- pmax(rc[, 1L], rc[, 2L])
      l_v <- pmin(rc[, 1L], rc[, 2L])

      Jt <- matrix(0, pstar, n_the)
      vech_pos <- sigma_lut[cbind(k_v, l_v)]
      if (delta.flag) {
        Jt[cbind(vech_pos, seq_len(n_the))] <- delta_weight[vech_pos]
      } else {
        Jt[cbind(vech_pos, seq_len(n_the))] <- 1
      }

      jac_sigma[, col:(col + n_the - 1L)] <- Jt
      col <- col + n_the
    }

    # Delta[k] (diagonal only):
    # dS[r,s] = I(r==k)*Sigma0[k,s]*delta[s] + delta[r]*Sigma0[r,k]*I(s==k)
    if (n_del > 0L) {
      Sigma0 <- M %*% MLIST$psi %*% t(M) + MLIST$theta
      k_v <- m.delta.idx
      T1 <- Sigma0[c_s, k_v, drop = FALSE] * delta[c_s] * outer(r_s, k_v, `==`)
      T2 <- delta[r_s] * Sigma0[r_s, k_v, drop = FALSE] * outer(c_s, k_v, `==`)
      jac_sigma[, col:(col + n_del - 1L)] <- T1 + T2
    }

  } else {
    # composite ov (cov) and composite lv (clv)
    cov.idx <- which(apply(
      MLIST$lambda, 1L,
      function(x) sum(x == 0) == ncol(MLIST$lambda)
    ))
    clv.idx <- which(apply(
      MLIST$lambda, 2L,
      function(x) sum(x == 0) == nrow(MLIST$lambda)
    ))

    LW  <- MLIST$lambda + MLIST$wmat
    Tmat <- diag(nvar); Tmat[cov.idx, cov.idx] <- MLIST$theta[cov.idx, cov.idx]
    wtw <- t(LW[, clv.idx, drop = FALSE]) %*% Tmat %*%
             LW[, clv.idx, drop = FALSE]
    wtw.inv <- solve(wtw)
    WTW.inv <- diag(nfac)
    WTW.inv[clv.idx, clv.idx] <- wtw.inv

    if (!beta.flag) {
      A <- diag(nfac)
    }

    # C0: V(eta) with composite diagonal entries zeroed
    C0 <- G; diag(C0)[clv.idx] <- 0

    # precompute
    Lambdafull <- Tmat %*% LW %*% WTW.inv
    LG_c <- Lambdafull %*% C0 %*% WTW.inv
    H_c <- Lambdafull %*% A
    F_c <- Lambdafull %*% G
    Sig0s <- Lambdafull %*% C0 %*% t(Lambdafull)
    LAMc <- Tmat %*% LW[, clv.idx, drop = FALSE] %*% wtw.inv
    # dLambdafull/dW[i,j] for i in ovc, j in clv expands to
    #   wtw.inv[j_loc,c_loc]*(Tmat - LAMc %*% wtw %*% t(LAMc))[r,i]
    #   - LAMc[i,c_loc]*LAMc[r,j_loc].
    # The "Pmat1" and "Pmat2" multipliers below feed into the wmat-block
    # T1+T2-T3-T4 expansion. Earlier versions used Tmat - LAMc %*% t(LAMc)
    # and LAMc %*% wtw.inv respectively, which are correct only when
    # wtw == I; the generally correct forms use the wtw factor explicitly:
    Pmat1 <- Tmat - LAMc %*% wtw %*% t(LAMc)
    Pmat2 <- LAMc

    # Beta/Psi correction matrices
    XrsXss <- Lambdafull[r_s, clv.idx, drop = FALSE] *
              Lambdafull[c_s, clv.idx, drop = FALSE]
    A_sub  <- A[clv.idx, , drop = FALSE]
    G_sub  <- G[clv.idx, , drop = FALSE]
    # lambda
    if (n_lam > 0L) {
      rc  <- vec2rc(m.lambda.idx, nvar);  i_v <- rc[, 1L];  j_v <- rc[, 2L]

      T1 <- LG_c[c_s, j_v, drop = FALSE] * outer(r_s, i_v, `==`)
      T2 <- LG_c[r_s, j_v, drop = FALSE] * outer(c_s, i_v, `==`)

      DX <- T1 + T2
      if (delta.flag) {
        DX <- DX * delta_weight
      }
      jac_sigma[, col:(col + n_lam - 1L)] <- DX
      col <- col + n_lam
    }

    # wmat
    if (n_wmat > 0L) {
      rc <- vec2rc(m.wmat.idx, nvar)
      i_v <- rc[, 1L]
      j_v <- rc[, 2L]
      j_loc <- match(j_v, clv.idx)

      # path 1: xi_i terms
      T1 <- Pmat1[r_s, i_v, drop = FALSE] * LG_c[c_s, j_v, drop = FALSE]
      T2 <- Pmat1[c_s, i_v, drop = FALSE] * LG_c[r_s, j_v, drop = FALSE]

      # path 2: phi_j * sigma_i correction
      T3 <- Pmat2[r_s, j_loc, drop = FALSE] * Sig0s[c_s, i_v, drop = FALSE]
      T4 <- Pmat2[c_s, j_loc, drop = FALSE] * Sig0s[r_s, i_v, drop = FALSE]

      DX <- T1 + T2 - T3 - T4
      if (delta.flag) {
        DX <- DX * delta_weight
      }
      jac_sigma[, col:(col + n_wmat - 1L)] <- DX
      col <- col + n_wmat
    }

    # beta
    if (n_bet > 0L) {
      rc  <- vec2rc(m.beta.idx, nfac);  i_v <- rc[, 1L];  j_v <- rc[, 2L]

      T1 <- H_c[r_s, i_v, drop = FALSE] * F_c[c_s, j_v, drop = FALSE]
      T2 <- F_c[r_s, j_v, drop = FALSE] * H_c[c_s, i_v, drop = FALSE]

      AG     <- A_sub[, i_v, drop = FALSE] * G_sub[, j_v, drop = FALSE]
      R_corr <- 2 * XrsXss %*% AG

      DX <- T1 + T2 - R_corr
      if (delta.flag) {
        DX <- DX * delta_weight
      }
      jac_sigma[, col:(col + n_bet - 1L)] <- DX
      col <- col + n_bet
    }

    # psi
    if (n_psi > 0L) {
      rc <- vec2rc(m.psi.idx, nfac)
      k_v <- pmax(rc[, 1L], rc[, 2L])
      l_v <- pmin(rc[, 1L], rc[, 2L])

      T1 <- H_c[r_s, k_v, drop = FALSE] * H_c[c_s, l_v, drop = FALSE]
      T2 <- H_c[r_s, l_v, drop = FALSE] * H_c[c_s, k_v, drop = FALSE]

      A_kl <- A_sub[, k_v, drop = FALSE] * A_sub[, l_v, drop = FALSE]
      R_corr <- 2 * XrsXss %*% A_kl

      DX <- T1 + T2 - R_corr

      diag_mask <- (k_v == l_v)
      if (any(diag_mask)) {
        DX[, diag_mask] <- DX[, diag_mask] * 0.5
      }

      if (delta.flag) {
        DX <- DX * delta_weight
      }
      jac_sigma[, col:(col + n_psi - 1L)] <- DX
      col <- col + n_psi
    }

    # theta
    if (n_the > 0L) {
      rc  <- vec2rc(m.theta.idx, nvar)
      k_v <- pmax(rc[, 1L], rc[, 2L])
      l_v <- pmin(rc[, 1L], rc[, 2L])

      Jt <- matrix(0, pstar, n_the)
      vech_pos <- sigma_lut[cbind(k_v, l_v)]
      if (delta.flag) {
        Jt[cbind(vech_pos, seq_len(n_the))] <- delta_weight[vech_pos]
      } else {
        Jt[cbind(vech_pos, seq_len(n_the))] <- 1
      }

      jac_sigma[, col:(col + n_the - 1L)] <- Jt
      col <- col + n_the
    }

    # delta
    if (n_del > 0L) {
      Sigma0_c <- Sig0s + MLIST$theta
      k_v      <- m.delta.idx
      T1 <- Sigma0_c[c_s, k_v, drop = FALSE] * delta[c_s] * outer(r_s, k_v, `==`)
      T2 <- delta[r_s] * Sigma0_c[r_s, k_v, drop = FALSE] * outer(c_s, k_v, `==`)
      jac_sigma[, col:(col + n_del - 1L)] <- T1 + T2
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
    Sigma_star <- lav_lisrel_sigma(MLIST = MLIST, delta = FALSE)
    s_diag <- sqrt(diag(Sigma_star))
    if (length(num.idx) > 0L) {
      s_diag[num.idx] <- 1.0
    }
    Sigma_obs <- Sigma_star / tcrossprod(s_diag)
    sigma_obs_v <- lav_matrix_vech(Sigma_obs)

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
    if (length(num.idx) > 0L) {
      is_ord[num.idx] <- FALSE
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
    cov.idx <- lav_matrix_vech_idx(nvar)
    covd.idx <- lav_matrix_vech_idx(nvar, diagonal = FALSE)

    var.idx <- which(is.na(match(
      cov.idx,
      covd.idx
    )))[num.idx]
    cor.idx <- match(covd.idx, cov.idx)

    jac_sigma <- rbind(
      jac_sigma[var.idx, , drop = FALSE],
      jac_sigma[cor.idx, , drop = FALSE]
    )
  }

  # correlation structure
  if (!categorical && correlation) {
    rm.idx <- lav_matrix_diagh_idx(nvar)
    jac_sigma <- jac_sigma[-rm.idx, , drop = FALSE]
  }

  jac_th   <- NULL
  jac_beta <- NULL
  nth_full <- 0L
  nexo_g   <- 0L
  if (conditional.x && !is.null(MLIST$gamma)) {
    nexo_g <- ncol(MLIST$gamma)
  }
  if (!categorical) {
    if (conditional.x) {
      # conditional.x continuous case: implied moments are
      #   mu_i  = nu_i + (LAMBDA(I-B)^{-1} alpha)_i
      #   pi_{i,k} = (LAMBDA(I-B)^{-1} GAMMA)_{i,k}
      #   Sigma = LAMBDA(I-B)^{-1} PSI (I-B)^{-T} LAMBDA' + THETA
      # The Beta-style row layout is interleaved column-major of a
      # (nexo+1) x nvar matrix (intercept first, then nexo slopes),
      # matching what lav_mvreg_scores_* expects.
      if (beta.flag) {
        a <- drop(A %*% MLIST$alpha)
      } else {
        a <- if (!is.null(MLIST$alpha)) as.numeric(MLIST$alpha) else
             rep(0, nfac)
      }
      gamma_IB <- if (nexo_g > 0L) {
        if (beta.flag) A %*% MLIST$gamma else MLIST$gamma
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
        jac_beta[cbind(mu_slots[m.nu.idx], col + seq_len(n_nu) - 1L)] <- 1
        col <- col + n_nu
      }

      # lambda[i,j]: dmu_i = a[j] * I(r==i); dpi_{i,k} = gamma_IB[j,k] * I(r==i).
      if (n_lam > 0L) {
        rc <- vec2rc(m.lambda.idx, nvar); i_v <- rc[, 1L]; j_v <- rc[, 2L]
        # mu-row: a[j_v] at row mu_slots[i_v]
        jac_beta[cbind(mu_slots[i_v], col + seq_len(n_lam) - 1L)] <- a[j_v]
        # pi-rows (k = 1..nexo)
        if (nexo_g > 0L) {
          for (k in seq_len(nexo_g)) {
            jac_beta[cbind(mu_slots[i_v] + k,
                           col + seq_len(n_lam) - 1L)] <- gamma_IB[j_v, k]
          }
        }
        col <- col + n_lam
      }

      # beta[i,j]: dmu_r = M[r,i]*a[j]; dpi_{r,k} = M[r,i]*gamma_IB[j,k].
      if (n_bet > 0L) {
        rc <- vec2rc(m.beta.idx, nfac); i_v <- rc[, 1L]; j_v <- rc[, 2L]
        # mu rows
        jac_beta[mu_slots, col:(col + n_bet - 1L)] <-
          M[, i_v, drop = FALSE] * rep(a[j_v], each = nvar)
        # pi rows
        if (nexo_g > 0L) {
          for (k in seq_len(nexo_g)) {
            jac_beta[mu_slots + k, col:(col + n_bet - 1L)] <-
              M[, i_v, drop = FALSE] *
                rep(gamma_IB[j_v, k], each = nvar)
          }
        }
        col <- col + n_bet
      }

      # alpha[k]: only affects mu (dmu = M[,k]).
      if (n_alp > 0L) {
        jac_beta[mu_slots, col:(col + n_alp - 1L)] <-
          M[, m.alpha.idx, drop = FALSE]
        col <- col + n_alp
      }

      # gamma[i,k]: only affects pi_{r,k} -> M[r,i].
      if (n_gam > 0L) {
        rc <- vec2rc(m.gamma.idx, nfac); i_v <- rc[, 1L]; k_v <- rc[, 2L]
        for (c2 in seq_len(n_gam)) {
          jac_beta[mu_slots + k_v[c2], col + c2 - 1L] <- M[, i_v[c2]]
        }
      }
    } else if (meanstructure) {
      # precompute: IB.inv * alpha
      if (beta.flag) {
        a <- drop(A %*% MLIST$alpha)
      } else {
        a <- MLIST$alpha
      }

      n_free <- n_nu + n_lam + n_bet + n_alp
      jac_mean <- matrix(0, nvar, n_free)

      col <- 1L
      # nu[i]:  dmu = e_i
      if (n_nu > 0L) {
        Jn <- matrix(0, nvar, n_nu)
        Jn[cbind(m.nu.idx, seq_len(n_nu))] <- 1
        jac_mean[, col:(col + n_nu - 1L)] <- Jn
        col <- col + n_nu
      }

      # Lambda[i,j]:  dmu[r] = a[j] * I(r == i)
      if (n_lam > 0L) {
        rc  <- vec2rc(m.lambda.idx, nvar); i_v <- rc[, 1L];  j_v <- rc[, 2L]

        Jl <- matrix(0, nvar, n_lam)
        Jl[cbind(i_v, seq_len(n_lam))] <- a[j_v]
        jac_mean[, col:(col + n_lam - 1L)] <- Jl
        col <- col + n_lam
      }

      # Beta[i,j]:  dmu[r] = M[r,i] * a[j]
      if (n_bet > 0L) {
        rc  <- vec2rc(m.beta.idx, nfac); i_v <- rc[, 1L];  j_v <- rc[, 2L]

        jac_mean[, col:(col + n_bet - 1L)] <-
          M[, i_v, drop = FALSE] * rep(a[j_v], each = nvar)
        col <- col + n_bet
      }

      # alpha[k]:  dmu = M[,k]
      if (n_alp > 0L) {
        jac_mean[, col:(col + n_alp - 1L)] <- M[, m.alpha.idx, drop = FALSE]
      }
    } # meanstructure


  # categorical
  } else if (categorical) {

    th_idx_g <- th.idx
    nth_full  <- length(th_idx_g)

    # v_slot[t] = observed-variable index for TH element t
    nlev_v <- tabulate(th_idx_g, nbins = nvar)
    nlev_v[nlev_v == 0L] <- 1L
    v_slot <- rep(seq_len(nvar), times = nlev_v)

    # per-slot delta scale
    if (delta.flag) {
      delta_star <- delta[v_slot]
    } else {
      delta_star <- rep(1, nth_full)
    }

    # ord_slots[k] = position in TH vector that corresponds to tau row k
    ord_slots <- which(th_idx_g > 0L)

    # a_vec = IB.inv * alpha  (for lambda/beta derivatives)
    if (!is.null(MLIST$alpha)) {
      a_vec <- if (beta.flag) {
        as.vector(A %*% MLIST$alpha)
      } else {
        as.vector(MLIST$alpha)
      }
    } else {
      a_vec <- rep(0, nfac)
    }

    # nu_vec (intercepts; zero for purely ordinal models)
    nu_vec <- if (!is.null(MLIST$nu)) as.vector(MLIST$nu) else rep(0, nvar)

    # pi0 = nu + Lambda * IB.inv * alpha = nu + M * alpha
    if (!is.null(MLIST$alpha)) {
      pi0 <- nu_vec + as.vector(M %*% as.vector(MLIST$alpha))
    } else {
      pi0 <- nu_vec
    }

    # tau_full: tau values at ordinal TH slots, 0 at numeric TH slots
    tau_full <- numeric(nth_full)
    if (!is.null(MLIST$tau) && n_th > 0L) {
      tau_full[ord_slots] <- as.vector(MLIST$tau)
    }

    # unscaled TH: TH0[t] = tau_full[t] - pi0[v(t)]
    TH0 <- tau_full - pi0[v_slot]

    # build jac_th  (nth_full rows, n_free_th columns)
    # column order: tau | delta | nu | lambda | beta | alpha
    n_free_th <- n_th + n_del + n_nu + n_lam + n_bet + n_alp
    jac_th <- matrix(0, nth_full, n_free_th)
    col_th <- 1L

    # tau[k]: delta_star[t] * I(t == ord_slots[k])
    if (n_th > 0L) {
      t_slots <- ord_slots[m.th.idx]
      Jt <- matrix(0, nth_full, n_th)
      Jt[cbind(t_slots, seq_len(n_th))] <- delta_star[t_slots]
      jac_th[, col_th:(col_th + n_th - 1L)] <- Jt
      col_th <- col_th + n_th
    }

    # delta[k]: I(v_slot[t]==k) * TH0[t]
    if (n_del > 0L) {
      jac_th[, col_th:(col_th + n_del - 1L)] <-
        outer(v_slot, m.delta.idx, `==`) * TH0
      col_th <- col_th + n_del
    }

    # nu[i]: -delta_star[t] * I(v_slot[t]==i)
    if (n_nu > 0L) {
      jac_th[, col_th:(col_th + n_nu - 1L)] <-
        -outer(v_slot, m.nu.idx, `==`) * delta_star
      col_th <- col_th + n_nu
    }

    # lambda[i,j]: -delta_star[t]*I(v_slot[t]==i)*a_vec[j]
    if (n_lam > 0L) {
      rc  <- vec2rc(m.lambda.idx, nvar); i_v <- rc[, 1L]; j_v <- rc[, 2L]
      jac_th[, col_th:(col_th + n_lam - 1L)] <-
        -(outer(v_slot, i_v, `==`) * delta_star) *
         matrix(a_vec[j_v], nth_full, n_lam, byrow = TRUE)
      col_th <- col_th + n_lam
    }

    # beta[i,j]: -delta_star[t]*M[v_slot[t],i]*a_vec[j]
    if (n_bet > 0L) {
      rc  <- vec2rc(m.beta.idx, nfac); i_v <- rc[, 1L]; j_v <- rc[, 2L]
      jac_th[, col_th:(col_th + n_bet - 1L)] <-
        -delta_star * M[v_slot, i_v, drop = FALSE] *
         matrix(a_vec[j_v], nth_full, n_bet, byrow = TRUE)
         matrix(a_vec[j_v], nth_full, n_bet, byrow = TRUE)
      col_th <- col_th + n_bet
    }

    # alpha[k]: -delta_star[t] * M[v_slot[t], k]
    if (n_alp > 0L) {
     jac_th[, col_th:(col_th + n_alp - 1L)] <-
       -delta_star * M[v_slot, m.alpha.idx, drop = FALSE]
    }
  }

  # categorical + conditional.x: build jac_pi (slope rows).
  #   pi_obs[r, k] = Delta[r] * (LAMBDA(I-B)^{-1} GAMMA)[r, k]
  # vec(PI) row layout is column-major: index (r, k) -> (k-1)*nvar + r.
  # Free params (in this column order): lambda | beta | gamma | delta.
  jac_pi <- NULL
  n_pi_rows <- 0L
  if (categorical && conditional.x && nexo_g > 0L) {
    n_pi_rows <- nvar * nexo_g

    gamma_IB <- if (beta.flag) A %*% MLIST$gamma else MLIST$gamma
    delta_var <- if (delta.flag) delta else rep(1, nvar)

    n_free_pi <- n_lam + n_bet + n_gam + n_del
    jac_pi <- matrix(0, n_pi_rows, n_free_pi)
    col_pi <- 1L

    # lambda[i,j]: dpi[r,k] = Delta[r] * I(r==i) * gamma_IB[j, k]
    if (n_lam > 0L) {
      rc <- vec2rc(m.lambda.idx, nvar); i_v <- rc[, 1L]; j_v <- rc[, 2L]
      for (k in seq_len(nexo_g)) {
        jac_pi[cbind((k - 1L) * nvar + i_v,
                     col_pi + seq_len(n_lam) - 1L)] <-
          delta_var[i_v] * gamma_IB[j_v, k]
      }
      col_pi <- col_pi + n_lam
    }

    # beta[i,j]: dpi[r,k] = Delta[r] * M[r,i] * gamma_IB[j, k]
    if (n_bet > 0L) {
      rc <- vec2rc(m.beta.idx, nfac); i_v <- rc[, 1L]; j_v <- rc[, 2L]
      for (k in seq_len(nexo_g)) {
        rows_k <- ((k - 1L) * nvar + 1L):(k * nvar)
        jac_pi[rows_k, col_pi:(col_pi + n_bet - 1L)] <-
          delta_var * M[, i_v, drop = FALSE] *
            rep(gamma_IB[j_v, k], each = nvar)
      }
      col_pi <- col_pi + n_bet
    }

    # gamma[i,k']: dpi[r,k] = Delta[r] * M[r,i] * I(k==k')
    if (n_gam > 0L) {
      rc <- vec2rc(m.gamma.idx, nfac); i_v <- rc[, 1L]; k_v <- rc[, 2L]
      for (c2 in seq_len(n_gam)) {
        rows_k <- ((k_v[c2] - 1L) * nvar + 1L):(k_v[c2] * nvar)
        jac_pi[rows_k, col_pi + c2 - 1L] <- delta_var * M[, i_v[c2]]
      }
      col_pi <- col_pi + n_gam
    }

    # delta[r] (delta-parameterization free): dpi[r,k] = (M GAMMA)[r, k] (unscaled).
    if (n_del > 0L) {
      M_GAMMA <- M %*% MLIST$gamma  # nvar x nexo
      for (c2 in seq_len(n_del)) {
        kd <- m.delta.idx[c2]
        for (k in seq_len(nexo_g)) {
          jac_pi[(k - 1L) * nvar + kd, col_pi + c2 - 1L] <- M_GAMMA[kd, k]
        }
      }
    }
  }

  # sigma
  out <- matrix(0, nrow = nrow(jac_sigma), ncol = nx.free)
  el.idx_sigma <- if (!wmat.flag) {
    c(x.lambda.idx, x.beta.idx, x.psi.idx, x.theta.idx, x.delta.idx)
  } else {
    c(x.lambda.idx, x.wmat.idx, x.beta.idx, x.psi.idx, x.theta.idx, x.delta.idx)
  }
  out[, el.idx_sigma] <- jac_sigma

  # meanstructure / conditional.x
  if (!categorical && conditional.x && !is.null(jac_beta)) {
    el.idx_beta <- c(x.nu.idx, x.lambda.idx, x.beta.idx, x.alpha.idx,
                     x.gamma.idx)
    outm <- matrix(0, nrow = nrow(jac_beta), ncol = nx.free)
    outm[, el.idx_beta] <- jac_beta
    out <- rbind(outm, out)
  } else if (!categorical && meanstructure) {
    # right order
    el.idx <- c(x.nu.idx, x.lambda.idx, x.beta.idx, x.alpha.idx)
    outm <- matrix(0, nrow = nrow(jac_mean), ncol = nx.free)
    outm[, el.idx] <- jac_mean
    out <- rbind(outm, out)
  } else if (categorical && !is.null(jac_th)) {
    el.idx_th <- c(
      if (n_th  > 0L) x.th.idx     else integer(0L),
      if (n_del > 0L) x.delta.idx  else integer(0L),
      if (n_nu  > 0L) x.nu.idx     else integer(0L),
      if (n_lam > 0L) x.lambda.idx else integer(0L),
      if (n_bet > 0L) x.beta.idx   else integer(0L),
      if (n_alp > 0L) x.alpha.idx  else integer(0L)
    )
    out_th <- matrix(0, nth_full, nx.free)
    out_th[, el.idx_th] <- jac_th

    # jac_diag_full only needed under theta parameterization (cached for
    # tau- and pi-side chain corrections).
    jac_diag_full <- NULL
    is_ord_v <- rep(TRUE, nvar)
    if (length(num.idx) > 0L) {
      is_ord_v[num.idx] <- FALSE
    }
    if (parameterization == "theta" && !is.null(jac_diag_theta)) {
      jac_diag_full <- matrix(0, nvar, nx.free)
      jac_diag_full[, el.idx_sigma] <- jac_diag_theta
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
      th_obs <- delta_star * TH0 * is_ord_v[v_slot]
      out_th <- out_th +
                (-0.5 * th_obs) * jac_diag_full[v_slot, , drop = FALSE]
    }

    out <- rbind(out_th, out)

    # categorical + conditional.x: insert pi rows between th and sigma.
    if (!is.null(jac_pi)) {
      el.idx_pi <- c(
        if (n_lam > 0L) x.lambda.idx else integer(0L),
        if (n_bet > 0L) x.beta.idx   else integer(0L),
        if (n_gam > 0L) x.gamma.idx  else integer(0L),
        if (n_del > 0L) x.delta.idx  else integer(0L)
      )
      out_pi <- matrix(0, n_pi_rows, nx.free)
      out_pi[, el.idx_pi] <- jac_pi

      # theta parameterization: chain rule for pi-side moments.
      #   pi_obs[r, k] = Delta_r * unscaled_pi[r, k], with Delta_r = 1/s_r
      #   d pi_obs[r, k] / dphi += -0.5 * pi_obs[r, k] / Sigma*_rr * dSigma*_rr
      # only for ordinal r. Numeric variables have Delta_r = 1 (constant).
      if (!is.null(jac_diag_full)) {
        # construct vec(PI_obs) following the column-major row layout
        unscaled_pi <- M %*% MLIST$gamma
        pi_obs_mat  <- if (delta.flag) {
          unscaled_pi * delta
        } else unscaled_pi
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
  if (wmat.flag && length(clv.idx) > 0L) {
    n_lvc <- length(clv.idx)

    # current free parameter vector
    x0_chain <- numeric(nx.free)
    for (mm in seq_along(MLIST)) {
      if (length(m.free.idx[[mm]])) {
        x0_chain[x.free.idx[[mm]]] <- MLIST[[mm]][m.free.idx[[mm]]]
      }
    }

    ML.tmpl    <- MLIST
    m.free.cap <- m.free.idx
    x.free.cap <- x.free.idx
    clv.cap    <- clv.idx
    meanstr.cap <- meanstructure
    compute_lvc <- function(xx) {
      ML <- ML.tmpl
      for (mm in seq_along(ML)) {
        if (length(m.free.cap[[mm]])) {
          ML[[mm]][m.free.cap[[mm]]] <- xx[x.free.cap[[mm]]]
        }
      }
      ML <- lav_lisrel_comp_set_intresvar(ML)
      o <- ML$psi[cbind(clv.cap, clv.cap)]
      if (meanstr.cap && !is.null(ML$alpha)) {
        o <- c(o, ML$alpha[clv.cap, 1L])
      }
      o
    }
    J_lvc <- lav_func_jacobian_complex(func = compute_lvc, x = x0_chain)
    J_psi <- J_lvc[seq_len(n_lvc), , drop = FALSE]

    # dSigma/dpsi*[k,k] for k in clv.idx (closed form):
    #   H_c[r,k]*H_c[s,k] - sum_{k' in clv} Lambdafull[r,k']*Lambdafull[s,k']*A[k',k]^2
    H_c.lvc <- H_c[, clv.idx, drop = FALSE]
    M_psi   <- H_c.lvc[r_s, , drop = FALSE] *
                 H_c.lvc[c_s, , drop = FALSE]
    A_lvc_lvc <- A[clv.idx, clv.idx, drop = FALSE]
    M_psi   <- M_psi - XrsXss %*% (A_lvc_lvc^2)

    sigma_chain <- M_psi %*% J_psi

    if (!categorical && meanstructure) {
      sigma.row.offset <- nvar
    } else if (categorical && !is.null(jac_th)) {
      sigma.row.offset <- nth_full
    } else {
      sigma.row.offset <- 0L
    }
    sig.rows <- sigma.row.offset + seq_len(nrow(jac_sigma))
    out[sig.rows, ] <- out[sig.rows, , drop = FALSE] + sigma_chain

    # dMu/dalpha*[k] for k in clv.idx: column M[, k]
    if (meanstructure && !categorical) {
      J_alpha  <- J_lvc[n_lvc + seq_len(n_lvc), , drop = FALSE]
      mu_chain <- M[, clv.idx, drop = FALSE] %*% J_alpha
      out[seq_len(nvar), ] <- out[seq_len(nvar), , drop = FALSE] + mu_chain
    }
  }

  # group weight: prepend a row with 1.0 at the gw column
  if (group.w.free) {
    mm.gw.idx <- which(mnames == "gw")
    if (length(mm.gw.idx) > 0L) {
      x.gw.idx <- x.free.idx[[mm.gw.idx]]
      out_gw <- matrix(0, nrow = 1L, ncol = nx.free)
      out_gw[1L, x.gw.idx] <- 1.0
      out <- rbind(out_gw, out)
    }
  }

  out
}
