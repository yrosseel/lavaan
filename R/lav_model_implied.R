# compute model implied statistics
# per block

# YR 7 May 2022: add cov.x and mean.x if conditional.x (so that we do
#                no longer depend on SampleStats)

lav_model_implied <- function(lavmodel = NULL, GLIST = NULL, delta = TRUE) { # nolint
  stopifnot(inherits(lavmodel, "lavModel"))

  # state or final?
  glist <- GLIST
  if (is.null(glist)) glist <- lavmodel@GLIST

  implied_fast <- lav_model_implied_fast(
    lavmodel = lavmodel, glist = glist,
    need_sigma = TRUE,
    need_mu = lavmodel@meanstructure,
    need_pi = lavmodel@conditional.x,
    need_th = lavmodel@categorical,
    delta = delta
  )
  # model-implied variance/covariance matrix ('sigma hat')
  sigma_hat <- implied_fast$sigma

  # model-implied mean structure ('mu hat')
  if (lavmodel@meanstructure) {
    mu_hat <- implied_fast$mu
  } else {
    mu_hat <- vector("list", length = lavmodel@nblocks)
  }

  # if conditional.x, slopes, cov.x, mean.x
  if (lavmodel@conditional.x) {
    slopes <- implied_fast$pi

    # per block, because for some blocks, cov.x may not exist
    cov_x <- vector("list", lavmodel@nblocks)
    mean_x <- vector("list", lavmodel@nblocks)
    for (b in seq_len(lavmodel@nblocks)) {
      mm_in_block <- (seq_len(lavmodel@nmat[b]) +
        cumsum(c(0, lavmodel@nmat))[b])
      mlist <- lavmodel@GLIST[mm_in_block]
      cov_x_idx <- which(names(mlist) == "cov.x")
      if (length(cov_x_idx) > 0L) {
        cov_x[[b]] <- mlist[[cov_x_idx]]
      } else {
        cov_x[[b]] <- matrix(0, 0L, 0L)
      }
      mean_x_idx <- which(names(mlist) == "mean.x")
      if (length(mean_x_idx) > 0L) {
        mean_x[[b]] <- mlist[[mean_x_idx]]
      } else {
        mean_x[[b]] <- matrix(0, 0L, 1L)
      }
    }
  } else {
    slopes <- vector("list", length = lavmodel@nblocks)
  }

  # if categorical, model-implied thresholds
  if (lavmodel@categorical) {
    th <- implied_fast$th
  } else {
    th <- vector("list", length = lavmodel@nblocks)
  }

  if (lavmodel@group.w.free) {
    w_idx <- which(names(lavmodel@GLIST) == "gw")
    gw <- unname(glist[w_idx])
    gw <- lapply(gw, as.numeric)
  } else {
    gw <- vector("list", length = lavmodel@nblocks)
  }

  if (lavmodel@conditional.x) {
    implied <- list(
      res.cov = sigma_hat, res.int = mu_hat,
      res.slopes = slopes, cov.x = cov_x, mean.x = mean_x,
      res.th = th, group.w = gw
    )
  } else {
    implied <- list(cov = sigma_hat, mean = mu_hat, th = th, group.w = gw)
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

  cov_1 <- vector("list", length = nblocks)
  mean_1 <- vector("list", length = nblocks)

  # reconstruct COV/MEAN per block
  for (b in seq_len(nblocks)) {
    res_sigma <- lavimplied$res.cov[[b]]
    res_slopes <- lavimplied$res.slopes[[b]]
    res_int <- lavimplied$res.int[[b]]
    s_xx <- lavimplied$cov.x[[b]]
    m_x <- lavimplied$mean.x[[b]]

    s_yx <- res_slopes %*% s_xx
    s_xy <- t(s_yx)
    s_yy <- res_sigma + tcrossprod(s_yx, res_slopes)
    cov_1[[b]] <- rbind(cbind(s_yy, s_yx), cbind(s_xy, s_xx))

    mu_y <- as.vector(res_int + res_slopes %*% m_x)
    mu_x <- as.vector(m_x)
    mean_1[[b]] <- matrix(c(mu_y, mu_x), ncol = 1L)
  }

  # we ignore res.th for now, as we do not support categorical data
  # in the two-level setting anyway
  implied <- list(
    cov = cov_1, mean = mean_1, th = lavimplied$res.th,
    group.w = lavimplied$group.w
  )

  implied
}
