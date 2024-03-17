# compute model implied statistics
# per block

# YR 7 May 2022: add cov.x and mean.x if conditional.x (so that we do
#                no longer depend on SampleStats)

lav_model_implied <- function(lavmodel = NULL, GLIST = NULL, delta = TRUE) {
  stopifnot(inherits(lavmodel, "lavModel"))

  # state or final?
  if (is.null(GLIST)) GLIST <- lavmodel@GLIST

  # model-implied variance/covariance matrix ('sigma hat')
  Sigma.hat <- computeSigmaHat(
    lavmodel = lavmodel, GLIST = GLIST,
    delta = delta
  )

  # model-implied mean structure ('mu hat')
  if (lavmodel@meanstructure) {
    Mu.hat <- computeMuHat(lavmodel = lavmodel, GLIST = GLIST)
  } else {
    Mu.hat <- vector("list", length = lavmodel@nblocks)
  }

  # if conditional.x, slopes, cov.x, mean.x
  if (lavmodel@conditional.x) {
    SLOPES <- computePI(lavmodel = lavmodel, GLIST = GLIST)

    # per block, because for some blocks, cov.x may not exist
    COV.X <- vector("list", lavmodel@nblocks)
    MEAN.X <- vector("list", lavmodel@nblocks)
    for (b in seq_len(lavmodel@nblocks)) {
      mm.in.block <- (seq_len(lavmodel@nmat[b]) +
        cumsum(c(0, lavmodel@nmat))[b])
      MLIST <- lavmodel@GLIST[mm.in.block]
      cov.x.idx <- which(names(MLIST) == "cov.x")
      if (length(cov.x.idx) > 0L) {
        COV.X[[b]] <- MLIST[[cov.x.idx]]
      } else {
        COV.X[[b]] <- matrix(0, 0L, 0L)
      }
      mean.x.idx <- which(names(MLIST) == "mean.x")
      if (length(mean.x.idx) > 0L) {
        MEAN.X[[b]] <- MLIST[[mean.x.idx]]
      } else {
        MEAN.X[[b]] <- matrix(0, 0L, 1L)
      }
    }
  } else {
    SLOPES <- vector("list", length = lavmodel@nblocks)
  }

  # if categorical, model-implied thresholds
  if (lavmodel@categorical) {
    TH <- computeTH(lavmodel = lavmodel, GLIST = GLIST)
  } else {
    TH <- vector("list", length = lavmodel@nblocks)
  }

  if (lavmodel@group.w.free) {
    w.idx <- which(names(lavmodel@GLIST) == "gw")
    GW <- unname(GLIST[w.idx])
    GW <- lapply(GW, as.numeric)
  } else {
    GW <- vector("list", length = lavmodel@nblocks)
  }

  if (lavmodel@conditional.x) {
    implied <- list(
      res.cov = Sigma.hat, res.int = Mu.hat,
      res.slopes = SLOPES, cov.x = COV.X, mean.x = MEAN.X,
      res.th = TH, group.w = GW
    )
  } else {
    implied <- list(cov = Sigma.hat, mean = Mu.hat, th = TH, group.w = GW)
  }

  implied
}

# convert 'conditional.x = TRUE' to 'conditional.x = FALSE'
lav_model_implied_cond2uncond <- function(lavimplied) {
  # check for res.cov
  if (is.null(lavimplied$res.cov[[1]])) {
    # already unconditional
    return(lavimplied)
  } else {
    nblocks <- length(lavimplied$res.cov)
  }

  COV <- vector("list", length = nblocks)
  MEAN <- vector("list", length = nblocks)

  # reconstruct COV/MEAN per block
  for (b in seq_len(nblocks)) {
    res.Sigma <- lavimplied$res.cov[[b]]
    res.slopes <- lavimplied$res.slopes[[b]]
    res.int <- lavimplied$res.int[[b]]
    S.xx <- lavimplied$cov.x[[b]]
    M.x <- lavimplied$mean.x[[b]]

    S.yx <- res.slopes %*% S.xx
    S.xy <- t(S.yx)
    S.yy <- res.Sigma + tcrossprod(S.yx, res.slopes)
    COV[[b]] <- rbind(cbind(S.yy, S.yx), cbind(S.xy, S.xx))

    Mu.y <- as.vector(res.int + res.slopes %*% M.x)
    Mu.x <- as.vector(M.x)
    MEAN[[b]] <- matrix(c(Mu.y, Mu.x), ncol = 1L)
  }

  # we ignore res.th for now, as we do not support categorical data
  # in the two-level setting anyway
  implied <- list(
    cov = COV, mean = MEAN, th = lavimplied$res.th,
    group.w = lavimplied$group.w
  )

  implied
}
