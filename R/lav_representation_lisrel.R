# and matrix-representation specific functions:
# - computeSigmaHat
# - computeMuHat
# - derivative.F

# initital version: YR 2011-01-21: LISREL stuff
# updates:          YR 2011-12-01: group specific extraction
#                   YR 2012-05-17: thresholds
#                   YR 2021-10-04: rename representation.LISREL -> lav_lisrel

lav_lisrel <- function(lavpartable = NULL,
                       target = NULL,
                       extra = FALSE,
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
  group.w.free <- any(lavpartable$lhs == "group" & lavpartable$op == "%")

  # gamma?only if conditional.x
  if (any(lavpartable$op %in% c("~", "<~") & lavpartable$exo == 1L)) {
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
      ov.names <- vnames(lavpartable, "ov.nox", block = g)
    } else {
      ov.names <- vnames(lavpartable, "ov", block = g)
    }
    nvar <- length(ov.names)
    lv.names <- vnames(lavpartable, "lv", block = g)
    nfac <- length(lv.names)
    ov.th <- vnames(lavpartable, "th", block = g)
    nth <- length(ov.th)
    ov.names.x <- vnames(lavpartable, "ov.x", block = g)
    nexo <- length(ov.names.x)
    ov.names.nox <- vnames(lavpartable, "ov.nox", block = g)

    # in this representation, we need to create 'phantom/dummy' latent
    # variables for all `x' and `y' variables not in lv.names
    # (only y if conditional.x = TRUE)

    # regression dummys
    if (gamma) {
      tmp.names <-
        unique(lavpartable$lhs[(lavpartable$op == "~" |
          lavpartable$op == "<~") &
          lavpartable$block == g])
      # new in 0.6-12: fix for multilevel + conditional.x: splitted ov.x
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

    # 2. "~" regressions
    if (gamma) {
      # gamma
      idx <- which(target$rhs %in% ov.names.x &
        target$block == g & (target$op == "~" |
        target$op == "<~"))
      tmp.mat[idx] <- "gamma"
      tmp.row[idx] <- match(target$lhs[idx], lv.names)
      tmp.col[idx] <- match(target$rhs[idx], ov.names.x)

      # beta
      idx <- which(!target$rhs %in% ov.names.x &
        target$block == g & (target$op == "~" |
        target$op == "<~"))
      tmp.mat[idx] <- "beta"
      tmp.row[idx] <- match(target$lhs[idx], lv.names)
      tmp.col[idx] <- match(target$rhs[idx], lv.names)
    } else {
      idx <- which(target$block == g & (target$op == "~" |
        target$op == "<~"))
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

    if (extra) {
      # mRows
      mmRows <- list(
        tau = nth,
        delta = nvar,
        nu = nvar,
        lambda = nvar,
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
      mmNames <- c("lambda", "theta", "psi")
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
computeEETA.LISREL <- function(MLIST = NULL, mean.x = NULL,
                               sample.mean = NULL,
                               ov.y.dummy.ov.idx = NULL,
                               ov.x.dummy.ov.idx = NULL,
                               ov.y.dummy.lv.idx = NULL,
                               ov.x.dummy.lv.idx = NULL) {
  BETA <- MLIST$beta
  GAMMA <- MLIST$gamma

  # ALPHA? (reconstruct, but no 'fix')
  ALPHA <- .internal_get_ALPHA(
    MLIST = MLIST, sample.mean = sample.mean,
    ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
    ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
    ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
    ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
  )

  # BETA?
  if (!is.null(BETA)) {
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
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
computeEETAx.LISREL <- function(MLIST = NULL, eXo = NULL, N = nrow(eXo),
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
  ALPHA <- .internal_get_ALPHA(
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
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
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
computeVETA.LISREL <- function(MLIST = NULL) {
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
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
    VETA <- tcrossprod(IB.inv %*% PSI, IB.inv)
  }

  VETA
}

# 4) VETAx
# compute V(ETA|x_i): variances/covariances of latent variables
#     V(ETA) = (I-B)^-1 PSI (I-B)^-T  + remove dummies
computeVETAx.LISREL <- function(MLIST = NULL, lv.dummy.idx = NULL) {
  PSI <- MLIST$psi
  BETA <- MLIST$beta

  # beta?
  if (is.null(BETA)) {
    VETA <- PSI
  } else {
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
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
# this is similar to computeMuHat but:
# - we ALWAYS compute NU+ALPHA, even if meanstructure=FALSE
# - never used if GAMMA, since we then have categorical variables, and the
#   'part 1' structure contains the (thresholds +) intercepts, not
#   the means
computeEY.LISREL <- function(MLIST = NULL, mean.x = NULL, sample.mean = NULL,
                             ov.y.dummy.ov.idx = NULL,
                             ov.x.dummy.ov.idx = NULL,
                             ov.y.dummy.lv.idx = NULL,
                             ov.x.dummy.lv.idx = NULL, delta = TRUE) {
  LAMBDA <- MLIST$lambda

  # get NU, but do not 'fix'
  NU <- .internal_get_NU(
    MLIST = MLIST, sample.mean = sample.mean,
    ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
    ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
    ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
    ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
  )

  # compute E(ETA)
  EETA <- computeEETA.LISREL(
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

# 2) EYx
# compute E(Y|x_i): expected value of observed, conditional on x_i
# E(Y|x_i) = NU + LAMBDA %*% E(eta|x_i)

# - if no eXo (and GAMMA):
#     E(ETA|x_i) = (I-B)^-1 ALPHA
#     we return a matrix of size [nobs x nfac] replicating E(ETA)
# - if eXo and GAMMA:
#     E(ETA|x_i) = (I-B)^-1 ALPHA + (I-B)^-1 GAMMA x_i
#     we return  a matrix of size [nobs x nfac]
#
# - we ALWAYS compute NU+ALPHA, even if meanstructure=FALSE
# - never used if GAMMA, since we then have categorical variables, and the
#   'part 1' structure contains the (thresholds +) intercepts, not
#   the means
computeEYx.LISREL <- function(MLIST = NULL,
                              eXo = NULL,
                              N = nrow(eXo),
                              sample.mean = NULL,
                              ov.y.dummy.ov.idx = NULL,
                              ov.x.dummy.ov.idx = NULL,
                              ov.y.dummy.lv.idx = NULL,
                              ov.x.dummy.lv.idx = NULL,
                              delta = TRUE) {
  LAMBDA <- MLIST$lambda

  # get NU, but do not 'fix'
  NU <- .internal_get_NU(
    MLIST = MLIST,
    sample.mean = sample.mean,
    ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
    ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
    ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
    ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
  )

  # compute E(ETA|x_i)
  EETAx <- computeEETAx.LISREL(
    MLIST = MLIST,
    eXo = eXo,
    N = N,
    sample.mean = sample.mean,
    ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
    ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
    ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
    ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
  )

  # EYx
  EYx <- sweep(tcrossprod(EETAx, LAMBDA), 2L, STATS = NU, FUN = "+")

  # if delta, scale
  if (delta && !is.null(MLIST$delta)) {
    EYx <- sweep(EYx, 2L, STATS = MLIST$delta, FUN = "*")
  }

  EYx
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
computeEYetax.LISREL <- function(MLIST = NULL,
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
    ALPHA <- .internal_get_ALPHA(
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
  NU <- .internal_get_NU(
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
computeEYetax2.LISREL <- function(MLIST = NULL,
                                  ETA = NULL,
                                  sample.mean = NULL,
                                  ov.y.dummy.ov.idx = NULL,
                                  ov.x.dummy.ov.idx = NULL,
                                  ov.y.dummy.lv.idx = NULL,
                                  ov.x.dummy.lv.idx = NULL,
                                  delta = TRUE) {
  LAMBDA <- MLIST$lambda
  BETA <- MLIST$beta


  # only if we have dummy ov.y, we need to compute the 'yhat' values
  # beforehand, and impute them in ETA[,ov.y]
  if (length(ov.y.dummy.lv.idx) > 0L) {
    # ALPHA? (reconstruct, but no 'fix')
    ALPHA <- .internal_get_ALPHA(
      MLIST = MLIST, sample.mean = sample.mean,
      ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
      ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
      ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
      ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
    )

    # keep all, but ov.y values
    OV.NOY <- ETA[, -ov.y.dummy.lv.idx, drop = FALSE]
    # ov.y rows, non-ov.y cols
    BETAY <- BETA[ov.y.dummy.lv.idx, -ov.y.dummy.lv.idx, drop = FALSE]
    # ov.y intercepts
    ALPHAY <- ALPHA[ov.y.dummy.lv.idx, , drop = FALSE]

    # impute ov.y values in ETA
    ETA[, ov.y.dummy.lv.idx] <-
      sweep(tcrossprod(OV.NOY, BETAY), 2L, STATS = ALPHAY, FUN = "+")
  }

  # get NU, but do not 'fix'
  NU <- .internal_get_NU(
    MLIST = MLIST,
    sample.mean = sample.mean,
    ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
    ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
    ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
    ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
  )

  # EYetax
  EYetax <- sweep(tcrossprod(ETA, LAMBDA), 2L, STATS = NU, FUN = "+")

  # if delta, scale
  if (delta && !is.null(MLIST$delta)) {
    EYetax <- sweep(EYetax, 2L, STATS = MLIST$delta, FUN = "*")
  }

  EYetax
}

# unconditional version
computeEYetax3.LISREL <- function(MLIST = NULL,
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
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST2)
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv
  }

  # compute model-implied means
  EY <- computeEY.LISREL(
    MLIST = MLIST, mean.x = mean.x,
    sample.mean = sample.mean,
    ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
    ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
    ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
    ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
  )

  EETA <- computeEETA.LISREL(
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
computeVY.LISREL <- function(MLIST = NULL) {
  LAMBDA <- MLIST$lambda
  THETA <- MLIST$theta

  VETA <- computeVETA.LISREL(MLIST = MLIST)
  VY <- tcrossprod(LAMBDA %*% VETA, LAMBDA) + THETA
  VY
}

# 5) VYx
# compute V(Y*|x_i) == model-implied covariance matrix
# this equals V(Y*) if no (explicit) eXo no GAMMA
computeVYx.LISREL <- computeSigmaHat.LISREL <- function(MLIST = NULL,
                                                        delta = TRUE) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  PSI <- MLIST$psi
  THETA <- MLIST$theta
  BETA <- MLIST$beta

  # beta?
  if (is.null(BETA)) {
    LAMBDA..IB.inv <- LAMBDA
  } else {
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv
  }

  # compute V(Y*|x_i)
  VYx <- tcrossprod(LAMBDA..IB.inv %*% PSI, LAMBDA..IB.inv) + THETA

  # if delta, scale
  if (delta && !is.null(MLIST$delta)) {
    DELTA <- diag(MLIST$delta[, 1L], nrow = nvar, ncol = nvar)
    VYx <- DELTA %*% VYx %*% DELTA
  }

  VYx
}

# 6) VYetax
# V(Y | eta_i, x_i) = THETA
computeVYetax.LISREL <- function(MLIST = NULL, delta = TRUE) {
  VYetax <- MLIST$theta
  nvar <- nrow(MLIST$theta)

  # if delta, scale
  if (delta && !is.null(MLIST$delta)) {
    DELTA <- diag(MLIST$delta[, 1L], nrow = nvar, ncol = nvar)
    VYetax <- DELTA %*% VYetax %*% DELTA
  }

  VYetax
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
computeMuHat.LISREL <- function(MLIST = NULL) {
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
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv
  }

  # compute Mu Hat
  Mu.hat <- NU + LAMBDA..IB.inv %*% ALPHA

  Mu.hat
}

# compute TH for a single block/group
computeTH.LISREL <- function(MLIST = NULL, th.idx = NULL, delta = TRUE) {
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
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
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
computePI.LISREL <- function(MLIST = NULL, delta = TRUE) {
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
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
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

computeLAMBDA.LISREL <- function(MLIST = NULL,
                                 ov.y.dummy.ov.idx = NULL,
                                 ov.x.dummy.ov.idx = NULL,
                                 ov.y.dummy.lv.idx = NULL,
                                 ov.x.dummy.lv.idx = NULL,
                                 remove.dummy.lv = FALSE) {
  lv.dummy.idx <- c(ov.y.dummy.lv.idx, ov.x.dummy.lv.idx)

  # fix LAMBDA
  LAMBDA <- MLIST$lambda
  if (length(ov.y.dummy.ov.idx) > 0L) {
    LAMBDA[ov.y.dummy.ov.idx, ] <- MLIST$beta[ov.y.dummy.lv.idx, ]
  }

  # remove dummy lv?
  if (remove.dummy.lv && length(lv.dummy.idx) > 0L) {
    LAMBDA <- LAMBDA[, -lv.dummy.idx, drop = FALSE]
  }

  LAMBDA
}

computeTHETA.LISREL <- function(MLIST = NULL,
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

computeNU.LISREL <- function(MLIST = NULL,
                             sample.mean = sample.mean,
                             ov.y.dummy.ov.idx = NULL,
                             ov.x.dummy.ov.idx = NULL,
                             ov.y.dummy.lv.idx = NULL,
                             ov.x.dummy.lv.idx = NULL) {
  # get NU, but do not 'fix'
  NU <- .internal_get_NU(
    MLIST = MLIST, sample.mean = sample.mean,
    ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
    ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
    ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
    ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
  )

  # ALPHA? (reconstruct, but no 'fix')
  ALPHA <- .internal_get_ALPHA(
    MLIST = MLIST, sample.mean = sample.mean,
    ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
    ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
    ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
    ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
  )

  ov.dummy.idx <- c(ov.y.dummy.ov.idx, ov.x.dummy.ov.idx)
  lv.dummy.idx <- c(ov.y.dummy.lv.idx, ov.x.dummy.lv.idx)

  # fix NU
  if (length(ov.dummy.idx) > 0L) {
    NU[ov.dummy.idx, 1] <- ALPHA[lv.dummy.idx, 1]
  }

  NU
}

# compute IB.inv
.internal_get_IB.inv <- function(MLIST = NULL) {
  BETA <- MLIST$beta
  nr <- nrow(MLIST$psi)

  if (!is.null(BETA)) {
    tmp <- -BETA
    tmp[lav_matrix_diag_idx(nr)] <- 1
    IB.inv <- solve(tmp)
  } else {
    IB.inv <- diag(nr)
  }

  IB.inv
}

# only if ALPHA=NULL but we need it anyway
# we 'reconstruct' ALPHA here (including dummy entries), no fixing
#
# without any dummy variables, this is just the zero vector
# but if we have dummy variables, we need to fill in their values
#
#
.internal_get_ALPHA <- function(MLIST = NULL, sample.mean = NULL,
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
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv
    LAMBDA..IB.inv.dummy <- LAMBDA..IB.inv[ov.dummy.idx, lv.dummy.idx]
    ALPHA[lv.dummy.idx] <-
      solve(LAMBDA..IB.inv.dummy) %*% sample.mean[ov.dummy.idx]
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
.internal_get_NU <- function(MLIST = NULL, sample.mean = NULL,
                             ov.y.dummy.ov.idx = NULL,
                             ov.x.dummy.ov.idx = NULL,
                             ov.y.dummy.lv.idx = NULL,
                             ov.x.dummy.lv.idx = NULL) {
  if (!is.null(MLIST$nu)) {
    return(MLIST$nu)
  }

  # if nexo > 0, substract lambda %*% EETA
  if (length(ov.x.dummy.ov.idx) > 0L) {
    EETA <- computeEETA.LISREL(MLIST,
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

.internal_get_KAPPA <- function(MLIST = NULL,
                                ov.y.dummy.ov.idx = NULL,
                                ov.x.dummy.ov.idx = NULL,
                                ov.y.dummy.lv.idx = NULL,
                                ov.x.dummy.lv.idx = NULL,
                                nexo = NULL) {
  nvar <- nrow(MLIST$lambda)
  if (!is.null(MLIST$gamma)) {
    this.nexo <- ncol(MLIST$gamma)
  } else if (!is.null(nexo)) {
    this.nexo <- nexo
  } else {
    lav_msg_stop(gettext("nexo not known"))
  }

  # create KAPPA
  KAPPA <- matrix(0, nvar, this.nexo)
  if (!is.null(MLIST$gamma)) {
    KAPPA[ov.y.dummy.ov.idx, ] <-
      MLIST$gamma[ov.y.dummy.lv.idx, , drop = FALSE]
  } else if (length(ov.x.dummy.ov.idx) > 0L) {
    KAPPA[ov.y.dummy.ov.idx, ] <-
      MLIST$beta[ov.y.dummy.lv.idx,
        ov.x.dummy.lv.idx,
        drop = FALSE
      ]
  }

  KAPPA
}


# old version of computeEYetax (using 'fixing')
computeYHATetax.LISREL <- function(MLIST = NULL, eXo = NULL, ETA = NULL,
                                   sample.mean = NULL,
                                   ov.y.dummy.ov.idx = NULL,
                                   ov.x.dummy.ov.idx = NULL,
                                   ov.y.dummy.lv.idx = NULL,
                                   ov.x.dummy.lv.idx = NULL,
                                   Nobs = 1L) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)
  lv.dummy.idx <- c(ov.y.dummy.lv.idx, ov.x.dummy.lv.idx)
  ov.dummy.idx <- c(ov.y.dummy.ov.idx, ov.x.dummy.ov.idx)

  # exogenous variables?
  if (is.null(eXo)) {
    nexo <- 0L
  } else {
    nexo <- ncol(eXo)
    # check ETA rows
    if (!(nrow(ETA) == 1L || nrow(ETA) == nrow(eXo))) {
      lav_msg_stop(gettext("!(nrow(ETA) == 1L || nrow(ETA) == nrow(eXo))"))
    }
  }

  # get NU
  NU <- .internal_get_NU(
    MLIST = MLIST, sample.mean = sample.mean,
    ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
    ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
    ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
    ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
  )

  # ALPHA? (reconstruct, but no 'fix')
  ALPHA <- .internal_get_ALPHA(
    MLIST = MLIST, sample.mean = sample.mean,
    ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
    ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
    ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
    ov.x.dummy.lv.idx = ov.x.dummy.lv.idx
  )

  # fix NU
  if (length(lv.dummy.idx) > 0L) {
    NU[ov.dummy.idx, 1L] <- ALPHA[lv.dummy.idx, 1L]
  }

  # fix LAMBDA (remove dummies) ## FIXME -- needed?
  LAMBDA <- MLIST$lambda
  if (length(lv.dummy.idx) > 0L) {
    LAMBDA <- LAMBDA[, -lv.dummy.idx, drop = FALSE]
    nfac <- ncol(LAMBDA)
    LAMBDA[ov.y.dummy.ov.idx, ] <-
      MLIST$beta[ov.y.dummy.lv.idx, seq_len(nfac), drop = FALSE]
  }

  # compute YHAT
  YHAT <- sweep(ETA %*% t(LAMBDA), MARGIN = 2, NU, "+")

  # Kappa + eXo?
  # note: Kappa elements are either in Gamma or in Beta
  if (nexo > 0L) {
    # create KAPPA
    KAPPA <- .internal_get_KAPPA(
      MLIST = MLIST,
      ov.y.dummy.ov.idx = ov.y.dummy.ov.idx,
      ov.x.dummy.ov.idx = ov.x.dummy.ov.idx,
      ov.y.dummy.lv.idx = ov.y.dummy.lv.idx,
      ov.x.dummy.lv.idx = ov.x.dummy.lv.idx,
      nexo = nexo
    )

    # expand YHAT if ETA only has 1 row
    if (nrow(YHAT) == 1L) {
      YHAT <- sweep(eXo %*% t(KAPPA), MARGIN = 2, STATS = YHAT, FUN = "+")
    } else {
      # add fixed part
      YHAT <- YHAT + (eXo %*% t(KAPPA))
    }

    # put back eXo
    if (length(ov.x.dummy.ov.idx) > 0L) {
      YHAT[, ov.x.dummy.ov.idx] <- eXo
    }
  } else {
    # duplicate?
    if (is.numeric(Nobs) && Nobs > 1L && nrow(YHAT) == 1L) {
      YHAT <- matrix(YHAT, Nobs, nvar, byrow = TRUE)
      # YHAT <- YHAT[ rep(1L, Nobs), ]
    }
  }

  # delta?
  # FIXME: not used here?
  # if(!is.null(DELTA)) {
  #    YHAT <- sweep(YHAT, MARGIN=2, DELTA, "*")
  # }

  YHAT
}


# deal with 'dummy' OV.X latent variables
# create additional matrices (eg GAMMA), and resize
# remove all ov.x related entries
MLIST2MLISTX <- function(MLIST = NULL,
                         ov.x.dummy.ov.idx = NULL,
                         ov.x.dummy.lv.idx = NULL) {
  lv.idx <- ov.x.dummy.lv.idx
  ov.idx <- ov.x.dummy.ov.idx
  if (length(lv.idx) == 0L) {
    return(MLIST)
  }
  if (!is.null(MLIST$gamma)) {
    nexo <- ncol(MLIST$gamma)
  } else {
    nexo <- length(ov.x.dummy.ov.idx)
  }
  nvar <- nrow(MLIST$lambda)
  nfac <- ncol(MLIST$lambda) - length(lv.idx)

  # copy
  MLISTX <- MLIST

  # fix LAMBDA:
  # - remove all ov.x related columns/rows
  MLISTX$lambda <- MLIST$lambda[-ov.idx, -lv.idx, drop = FALSE]

  # fix THETA:
  # - remove ov.x related columns/rows
  MLISTX$theta <- MLIST$theta[-ov.idx, -ov.idx, drop = FALSE]

  # fix PSI:
  # - remove ov.x related columns/rows
  MLISTX$psi <- MLIST$psi[-lv.idx, -lv.idx, drop = FALSE]

  # create GAMMA
  if (length(ov.x.dummy.lv.idx) > 0L) {
    MLISTX$gamma <- MLIST$beta[-lv.idx, lv.idx, drop = FALSE]
  }

  # fix BETA (remove if empty)
  if (!is.null(MLIST$beta)) {
    MLISTX$beta <- MLIST$beta[-lv.idx, -lv.idx, drop = FALSE]
    if (ncol(MLISTX$beta) == 0L) MLISTX$beta <- NULL
  }

  # fix NU
  if (!is.null(MLIST$nu)) {
    MLISTX$nu <- MLIST$nu[-ov.idx, 1L, drop = FALSE]
  }

  # fix ALPHA
  if (!is.null(MLIST$alpha)) {
    MLISTX$alpha <- MLIST$alpha[-lv.idx, 1L, drop = FALSE]
  }

  MLISTX
}


# create MLIST from MLISTX
MLISTX2MLIST <- function(MLISTX = NULL,
                         ov.x.dummy.ov.idx = NULL,
                         ov.x.dummy.lv.idx = NULL,
                         mean.x = NULL,
                         cov.x = NULL) {
  lv.idx <- ov.x.dummy.lv.idx
  ndum <- length(lv.idx)
  ov.idx <- ov.x.dummy.ov.idx
  if (length(lv.idx) == 0L) {
    return(MLISTX)
  }
  stopifnot(!is.null(cov.x), !is.null(mean.x))
  nvar <- nrow(MLISTX$lambda)
  nfac <- ncol(MLISTX$lambda)

  # copy
  MLIST <- MLISTX

  # resize matrices
  MLIST$lambda <- rbind(
    cbind(MLISTX$lambda, matrix(0, nvar, ndum)),
    matrix(0, ndum, nfac + ndum)
  )
  MLIST$psi <- rbind(
    cbind(MLISTX$psi, matrix(0, nfac, ndum)),
    matrix(0, ndum, nfac + ndum)
  )
  MLIST$theta <- rbind(
    cbind(MLISTX$theta, matrix(0, nvar, ndum)),
    matrix(0, ndum, nvar + ndum)
  )
  if (!is.null(MLISTX$beta)) {
    MLIST$beta <- rbind(
      cbind(MLISTX$beta, matrix(0, nfac, ndum)),
      matrix(0, ndum, nfac + ndum)
    )
  }
  if (!is.null(MLISTX$alpha)) {
    MLIST$alpha <- rbind(MLISTX$alpha, matrix(0, ndum, 1))
  }
  if (!is.null(MLISTX$nu)) {
    MLIST$nu <- rbind(MLISTX$nu, matrix(0, ndum, 1))
  }

  # fix LAMBDA:
  # - add columns for all dummy latent variables
  MLIST$lambda[cbind(ov.idx, lv.idx)] <- 1

  # fix PSI
  # - move cov.x elements to PSI
  MLIST$psi[lv.idx, lv.idx] <- cov.x

  # move (ov.x.dummy elements of) GAMMA to BETA
  MLIST$beta[seq_len(nfac), ov.x.dummy.lv.idx] <- MLISTX$gamma
  MLIST$gamma <- NULL

  # fix ALPHA
  if (!is.null(MLIST$alpha)) {
    MLIST$alpha[lv.idx] <- mean.x
  }

  MLIST
}

# if DELTA parameterization, compute residual elements (in theta, or psi)
# of observed categorical variables, as a function of other model parameters
setResidualElements.LISREL <- function(MLIST = NULL,
                                       num.idx = NULL,
                                       ov.y.dummy.ov.idx = NULL,
                                       ov.y.dummy.lv.idx = NULL) {
  # remove num.idx from ov.y.dummy.*
  if (length(num.idx) > 0L && length(ov.y.dummy.ov.idx) > 0L) {
    n.idx <- which(ov.y.dummy.ov.idx %in% num.idx)
    if (length(n.idx) > 0L) {
      ov.y.dummy.ov.idx <- ov.y.dummy.ov.idx[-n.idx]
      ov.y.dummy.lv.idx <- ov.y.dummy.lv.idx[-n.idx]
    }
  }

  # force non-numeric theta elements to be zero
  if (length(num.idx) > 0L) {
    diag(MLIST$theta)[-num.idx] <- 0.0
  } else {
    diag(MLIST$theta) <- 0.0
  }
  if (length(ov.y.dummy.ov.idx) > 0L) {
    MLIST$psi[cbind(ov.y.dummy.lv.idx, ov.y.dummy.lv.idx)] <- 0.0
  }

  # special case: PSI=0, and lambda=I (eg ex3.12)
  if (ncol(MLIST$psi) > 0L &&
    sum(diag(MLIST$psi)) == 0.0 && all(diag(MLIST$lambda) == 1)) {
    ### FIXME: more elegant/general solution??
    diag(MLIST$psi) <- 1
    Sigma.hat <- computeSigmaHat.LISREL(MLIST = MLIST, delta = FALSE)
    diag.Sigma <- diag(Sigma.hat) - 1.0
  } else if (ncol(MLIST$psi) == 0L) {
    diag.Sigma <- rep(0, ncol(MLIST$theta))
  } else {
    Sigma.hat <- computeSigmaHat.LISREL(MLIST = MLIST, delta = FALSE)
    diag.Sigma <- diag(Sigma.hat)
  }

  if (is.null(MLIST$delta)) {
    delta <- rep(1, length(diag.Sigma))
  } else {
    delta <- MLIST$delta
  }
  # theta = DELTA^(-2) - diag( LAMBDA (I-B)^-1 PSI (I-B)^-T t(LAMBDA) )
  RESIDUAL <- as.vector(1 / (delta * delta) - diag.Sigma)
  if (length(num.idx) > 0L) {
    diag(MLIST$theta)[-num.idx] <- RESIDUAL[-num.idx]
  } else {
    diag(MLIST$theta) <- RESIDUAL
  }

  # move ov.y.dummy 'RESIDUAL' elements from THETA to PSI
  if (length(ov.y.dummy.ov.idx) > 0L) {
    MLIST$psi[cbind(ov.y.dummy.lv.idx, ov.y.dummy.lv.idx)] <-
      MLIST$theta[cbind(ov.y.dummy.ov.idx, ov.y.dummy.ov.idx)]
    MLIST$theta[cbind(ov.y.dummy.ov.idx, ov.y.dummy.ov.idx)] <- 0.0
  }

  MLIST
}

# if THETA parameterization, compute delta elements
# of observed categorical variables, as a function of other model parameters
setDeltaElements.LISREL <- function(MLIST = NULL, num.idx = NULL) {
  Sigma.hat <- computeSigmaHat.LISREL(MLIST = MLIST, delta = FALSE)
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
computeCOV.LISREL <- function(MLIST = NULL, delta = TRUE) {
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
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
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
derivative.F.LISREL <- function(MLIST = NULL, Omega = NULL, Omega.mu = NULL) {
  LAMBDA <- MLIST$lambda
  PSI <- MLIST$psi
  BETA <- MLIST$beta
  ALPHA <- MLIST$alpha

  # beta?
  if (is.null(BETA)) {
    LAMBDA..IB.inv <- LAMBDA
  } else {
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
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
  } else {
    BETA.deriv <- NULL
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
  } else {
    NU.deriv <- NULL
    ALPHA.deriv <- NULL
  }

  if (group.w.free) {
    GROUP.W.deriv <- 0.0
  } else {
    GROUP.W.deriv <- NULL
  }

  list(
    lambda = LAMBDA.deriv,
    beta = BETA.deriv,
    theta = THETA.deriv,
    psi = PSI.deriv,
    nu = NU.deriv,
    alpha = ALPHA.deriv,
    gw = GROUP.W.deriv
  )
}

# dSigma/dx -- per model matrix
# note:
# we avoid using the duplication and elimination matrices
# for now (perhaps until we'll use the Matrix package)
derivative.sigma.LISREL_OLD <- function(m = "lambda",
                                        # all model matrix elements, or only a few?
                                        # NOTE: for symmetric matrices,
                                        # we assume that the have full size
                                        # (nvar*nvar) (but already correct for
                                        # symmetry)
                                        idx = seq_len(length(MLIST[[m]])),
                                        MLIST = NULL,
                                        delta = TRUE) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)
  PSI <- MLIST$psi

  # only lower.tri part of sigma (not same order as elimination matrix?)
  v.idx <- lav_matrix_vech_idx(nvar)
  pstar <- nvar * (nvar + 1) / 2

  # shortcut for gamma, nu, alpha and tau: empty matrix
  if (m == "nu" || m == "alpha" || m == "tau" || m == "gamma" || m == "gw" ||
    m == "cov.x" || m == "mean.x") {
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
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
  }

  # pre
  if (m == "lambda" || m == "beta") {
    IK <- diag(nvar * nvar) + lav_matrix_commutation(nvar, nvar)
  }
  if (m == "lambda" || m == "beta") {
    IB.inv..PSI..tIB.inv..tLAMBDA <-
      IB.inv %*% PSI %*% t(IB.inv) %*% t(LAMBDA)
  }
  if (m == "beta" || m == "psi") {
    LAMBDA..IB.inv <- LAMBDA %*% IB.inv
  }

  # here we go:
  if (m == "lambda") {
    DX <- IK %*% t(IB.inv..PSI..tIB.inv..tLAMBDA %x% diag(nvar))
    if (delta.flag) {
      DX <- DX * as.vector(DELTA %x% DELTA)
    }
  } else if (m == "beta") {
    DX <- IK %*% (t(IB.inv..PSI..tIB.inv..tLAMBDA) %x% LAMBDA..IB.inv)
    # this is not really needed (because we select idx=m.el.idx)
    # but just in case we need all elements of beta...
    DX[, lav_matrix_diag_idx(nfac)] <- 0.0
    if (delta.flag) {
      DX <- DX * as.vector(DELTA %x% DELTA)
    }
  } else if (m == "psi") {
    DX <- (LAMBDA..IB.inv %x% LAMBDA..IB.inv)
    # symmetry correction, but keeping all duplicated elements
    # since we depend on idx=m.el.idx
    # otherwise, we could simply postmultiply with the duplicationMatrix

    # we sum up lower.tri + upper.tri (but not the diagonal elements!)
    # imatrix <- matrix(1:nfac^2,nfac,nfac)
    # lower.idx <- imatrix[lower.tri(imatrix, diag=FALSE)]
    # upper.idx <- imatrix[upper.tri(imatrix, diag=FALSE)]
    lower.idx <- lav_matrix_vech_idx(nfac, diagonal = FALSE)
    upper.idx <- lav_matrix_vechru_idx(nfac, diagonal = FALSE)
    # NOTE YR: upper.idx (see 3 lines up) is wrong in MH patch!
    # fixed again 13/06/2012 after bug report of Mijke Rhemtulla.

    offdiagSum <- DX[, lower.idx] + DX[, upper.idx]
    DX[, c(lower.idx, upper.idx)] <- cbind(offdiagSum, offdiagSum)
    if (delta.flag) {
      DX <- DX * as.vector(DELTA %x% DELTA)
    }
  } else if (m == "theta") {
    DX <- diag(nvar * nvar) # very sparse...
    # symmetry correction not needed, since all off-diagonal elements
    # are zero?
    if (delta.flag) {
      DX <- DX * as.vector(DELTA %x% DELTA)
    }
  } else if (m == "delta") {
    Omega <- computeSigmaHat.LISREL(MLIST, delta = FALSE)
    DD <- diag(DELTA[, 1], nvar, nvar)
    DD.Omega <- (DD %*% Omega)
    A <- DD.Omega %x% diag(nvar)
    B <- diag(nvar) %x% DD.Omega
    DX <- A[, lav_matrix_diag_idx(nvar), drop = FALSE] +
      B[, lav_matrix_diag_idx(nvar), drop = FALSE]
  } else {
    lav_msg_stop(gettext("wrong model matrix names:"), m)
  }

  DX <- DX[v.idx, idx, drop = FALSE]
  DX
}

# dSigma/dx -- per model matrix
derivative.sigma.LISREL <- function(m = "lambda",
                                    # all model matrix elements, or only a few?
                                    # NOTE: for symmetric matrices,
                                    # we assume that the have full size
                                    # (nvar*nvar) (but already correct for
                                    # symmetry)
                                    idx = seq_len(length(MLIST[[m]])),
                                    MLIST = NULL,
                                    vech = TRUE,
                                    delta = TRUE) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)
  PSI <- MLIST$psi

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
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
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
    KOL.idx <- matrix(1:(nfac * nfac), nfac, nfac, byrow = TRUE)[idx]
    DX <- (L1 %x% LAMBDA..IB.inv)[, idx, drop = FALSE] +
      (LAMBDA..IB.inv %x% L1)[, KOL.idx, drop = FALSE]
    # this is not really needed (because we select idx=m.el.idx)
    # but just in case we need all elements of beta...
    DX[, which(idx %in% lav_matrix_diag_idx(nfac))] <- 0.0
  } else if (m == "psi") {
    DX <- (LAMBDA..IB.inv %x% LAMBDA..IB.inv)
    # symmetry correction, but keeping all duplicated elements
    # since we depend on idx=m.el.idx
    lower.idx <- lav_matrix_vech_idx(nfac, diagonal = FALSE)
    upper.idx <- lav_matrix_vechru_idx(nfac, diagonal = FALSE)
    offdiagSum <- DX[, lower.idx] + DX[, upper.idx]
    DX[, c(lower.idx, upper.idx)] <- cbind(offdiagSum, offdiagSum)
    DX <- DX[, idx, drop = FALSE]
  } else if (m == "theta") {
    # DX <- diag(nvar*nvar) # very sparse...
    DX <- matrix(0, nvar * nvar, length(idx))
    DX[cbind(idx, seq_along(idx))] <- 1
    # symmetry correction not needed, since all off-diagonal elements
    # are zero?
  } else if (m == "delta") {
    Omega <- computeSigmaHat.LISREL(MLIST, delta = FALSE)
    DD <- diag(DELTA[, 1], nvar, nvar)
    DD.Omega <- (DD %*% Omega)
    A <- DD.Omega %x% diag(nvar)
    B <- diag(nvar) %x% DD.Omega
    DX <- A[, lav_matrix_diag_idx(nvar), drop = FALSE] +
      B[, lav_matrix_diag_idx(nvar), drop = FALSE]
    DX <- DX[, idx, drop = FALSE]
  } else {
    lav_msg_stop(gettext("wrong model matrix names:"), m)
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
derivative.mu.LISREL <- function(m = "alpha",
                                 # all model matrix elements, or only a few?
                                 idx = seq_len(length(MLIST[[m]])),
                                 MLIST = NULL) {
  LAMBDA <- MLIST$lambda
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)

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
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
  }

  if (m == "nu") {
    DX <- diag(nvar)
  } else if (m == "lambda") {
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
derivative.th.LISREL <- function(m = "tau",
                                 # all model matrix elements, or only a few?
                                 idx = seq_len(length(MLIST[[m]])),
                                 th.idx = NULL,
                                 MLIST = NULL,
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
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
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
derivative.pi.LISREL <- function(m = "lambda",
                                 # all model matrix elements, or only a few?
                                 idx = seq_len(length(MLIST[[m]])),
                                 MLIST = NULL) {
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
    IB.inv <- .internal_get_IB.inv(MLIST = MLIST)
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
derivative.gw.LISREL <- function(m = "gw",
                                 # all model matrix elements, or only a few?
                                 idx = seq_len(length(MLIST[[m]])),
                                 MLIST = NULL) {
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
derivative.lambda.LISREL <- function(m = "lambda",
                                     # all model matrix elements, or only a few?
                                     idx = seq_len(length(MLIST[[m]])),
                                     MLIST = NULL) {
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
derivative.psi.LISREL <- function(m = "psi",
                                  # all model matrix elements, or only a few?
                                  idx = seq_len(length(MLIST[[m]])),
                                  MLIST = NULL) {
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
derivative.theta.LISREL <- function(m = "theta",
                                    # all model matrix elements, or only a few?
                                    idx = seq_len(length(MLIST[[m]])),
                                    MLIST = NULL) {
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
derivative.beta.LISREL <- function(m = "beta",
                                   # all model matrix elements, or only a few?
                                   idx = seq_len(length(MLIST[[m]])),
                                   MLIST = NULL) {
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
derivative.gamma.LISREL <- function(m = "gamma",
                                    # all model matrix elements, or only a few?
                                    idx = seq_len(length(MLIST[[m]])),
                                    MLIST = NULL) {
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
derivative.nu.LISREL <- function(m = "nu",
                                 # all model matrix elements, or only a few?
                                 idx = seq_len(length(MLIST[[m]])),
                                 MLIST = NULL) {
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
derivative.tau.LISREL <- function(m = "tau",
                                  # all model matrix elements, or only a few?
                                  idx = seq_len(length(MLIST[[m]])),
                                  MLIST = NULL) {
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
derivative.alpha.LISREL <- function(m = "alpha",
                                    # all model matrix elements, or only a few?
                                    idx = seq_len(length(MLIST[[m]])),
                                    MLIST = NULL) {
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
# computeSigmaHat.LISREL <- lavaan:::computeSigmaHat.LISREL
# setDeltaElements.LISREL <- lavaan:::setDeltaElements.LISREL
TESTING_derivatives.LISREL <- function(MLIST = NULL,
                                       nvar = NULL, nfac = NULL, nexo = NULL,
                                       th.idx = NULL, num.idx = NULL,
                                       meanstructure = TRUE,
                                       th = TRUE, delta = TRUE, pi = TRUE,
                                       gw = FALSE, theta = FALSE,
                                       debug = FALSE) {
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
      mlist <- setDeltaElements.LISREL(MLIST = mlist, num.idx = num.idx)
    }
    lav_matrix_vech(computeSigmaHat.LISREL(mlist))
  }

  compute.mu <- function(x, mm = "lambda", MLIST = NULL) {
    mlist <- MLIST
    if (mm %in% c("psi", "theta")) {
      mlist[[mm]] <- lav_matrix_vech_reverse(x)
    } else {
      mlist[[mm]][, ] <- x
    }
    if (theta) {
      mlist <- setDeltaElements.LISREL(MLIST = mlist, num.idx = num.idx)
    }
    computeMuHat.LISREL(mlist)
  }

  compute.th2 <- function(x, mm = "tau", MLIST = NULL, th.idx) {
    mlist <- MLIST
    if (mm %in% c("psi", "theta")) {
      mlist[[mm]] <- lav_matrix_vech_reverse(x)
    } else {
      mlist[[mm]][, ] <- x
    }
    if (theta) {
      mlist <- setDeltaElements.LISREL(MLIST = mlist, num.idx = num.idx)
    }
    computeTH.LISREL(mlist, th.idx = th.idx)
  }

  compute.pi <- function(x, mm = "lambda", MLIST = NULL) {
    mlist <- MLIST
    if (mm %in% c("psi", "theta")) {
      mlist[[mm]] <- lav_matrix_vech_reverse(x)
    } else {
      mlist[[mm]][, ] <- x
    }
    if (theta) {
      mlist <- setDeltaElements.LISREL(MLIST = mlist, num.idx = num.idx)
    }
    computePI.LISREL(mlist)
  }

  compute.gw <- function(x, mm = "gw", MLIST = NULL) {
    mlist <- MLIST
    if (mm %in% c("psi", "theta")) {
      mlist[[mm]] <- lav_matrix_vech_reverse(x)
    } else {
      mlist[[mm]][, ] <- x
    }
    if (theta) {
      mlist <- setDeltaElements.LISREL(MLIST = mlist, num.idx = num.idx)
    }
    mlist$gw[1, 1]
  }

  # if theta, set MLIST$delta
  if (theta) {
    MLIST <- setDeltaElements.LISREL(MLIST = MLIST, num.idx = num.idx)
  }

  for (mm in names(MLIST)) {
    if (mm %in% c("psi", "theta")) {
      x <- lav_matrix_vech(MLIST[[mm]])
    } else {
      x <- lav_matrix_vec(MLIST[[mm]])
    }
    if (mm == "delta" && theta) next
    if (debug) {
      cat("### mm = ", mm, "\n")
    }

    # 1. sigma
    DX1 <- lav_func_jacobian_complex(func = compute.sigma, x = x, mm = mm, MLIST = MLIST)
    DX2 <- derivative.sigma.LISREL(
      m = mm, idx = seq_len(length(MLIST[[mm]])),
      MLIST = MLIST, delta = !theta
    )
    if (mm %in% c("psi", "theta")) {
      # remove duplicated columns of symmetric matrices
      idx <- lav_matrix_vechru_idx(sqrt(ncol(DX2)), diagonal = FALSE)
      if (length(idx) > 0L) DX2 <- DX2[, -idx]
    }
    if (theta) {
      sigma.hat <- computeSigmaHat.LISREL(MLIST = MLIST, delta = FALSE)
      R <- lav_deriv_cov2cor(sigma.hat, num.idx = num.idx)

      DX3 <- DX2
      DX2 <- R %*% DX2
    }
    if (debug) {
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
    DX2 <- derivative.mu.LISREL(
      m = mm, idx = seq_len(length(MLIST[[mm]])),
      MLIST = MLIST
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
    if (debug) {
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
      DX2 <- derivative.th.LISREL(
        m = mm, idx = seq_len(length(MLIST[[mm]])),
        MLIST = MLIST, th.idx = th.idx,
        delta = TRUE
      )
      if (theta) {
        # 1. compute dDelta.dx
        dxSigma <-
          derivative.sigma.LISREL(
            m = mm, idx = seq_len(length(MLIST[[mm]])),
            MLIST = MLIST, delta = !theta
          )
        var.idx <- which(!lav_matrix_vech_idx(nvar) %in%
          lav_matrix_vech_idx(nvar, diagonal = FALSE))
        sigma.hat <- computeSigmaHat.LISREL(MLIST = MLIST, delta = FALSE)
        dsigma <- diag(sigma.hat)
        # dy/ddsigma = -0.5/(ddsigma*sqrt(ddsigma))
        dDelta.dx <- dxSigma[var.idx, ] * -0.5 / (dsigma * sqrt(dsigma))

        # 2. compute dth.dDelta
        dth.dDelta <-
          derivative.th.LISREL(
            m = "delta",
            idx = seq_len(length(MLIST[["delta"]])),
            MLIST = MLIST, th.idx = th.idx
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
      if (debug) {
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
      DX2 <- derivative.pi.LISREL(
        m = mm, idx = seq_len(length(MLIST[[mm]])),
        MLIST = MLIST
      )
      if (mm %in% c("psi", "theta")) {
        # remove duplicated columns of symmetric matrices
        idx <- lav_matrix_vechru_idx(sqrt(ncol(DX2)), diagonal = FALSE)
        if (length(idx) > 0L) DX2 <- DX2[, -idx]
      }
      if (theta) {
        # 1. compute dDelta.dx
        dxSigma <-
          derivative.sigma.LISREL(
            m = mm, idx = seq_len(length(MLIST[[mm]])),
            MLIST = MLIST, delta = !theta
          )
        if (mm %in% c("psi", "theta")) {
          # remove duplicated columns of symmetric matrices
          idx <- lav_matrix_vechru_idx(sqrt(ncol(dxSigma)), diagonal = FALSE)
          if (length(idx) > 0L) dxSigma <- dxSigma[, -idx]
        }
        var.idx <- which(!lav_matrix_vech_idx(nvar) %in%
          lav_matrix_vech_idx(nvar, diagonal = FALSE))
        sigma.hat <- computeSigmaHat.LISREL(MLIST = MLIST, delta = FALSE)
        dsigma <- diag(sigma.hat)
        # dy/ddsigma = -0.5/(ddsigma*sqrt(ddsigma))
        dDelta.dx <- dxSigma[var.idx, ] * -0.5 / (dsigma * sqrt(dsigma))

        # 2. compute dpi.dDelta
        dpi.dDelta <-
          derivative.pi.LISREL(
            m = "delta",
            idx = seq_len(length(MLIST[["delta"]])),
            MLIST = MLIST
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
      if (debug) {
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
      DX2 <- derivative.gw.LISREL(
        m = mm, idx = seq_len(length(MLIST[[mm]])),
        MLIST = MLIST
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
      if (debug) {
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
