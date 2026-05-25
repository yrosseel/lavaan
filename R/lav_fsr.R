# compute the jacobian: dtheta_2/dtheta_1:
#
# theta_2: - in the rows
#          - the croon corrections, expressed as
#              1) scaled offsets (scoffset), and
#              2) scaling factors
# theta_1: - in the columns
#          - the free parameters of the measurement model
#
lav_fsr_delta21 <- function(object, fsm = NULL) {
  lavmodel <- object@Model
  nmat <- lavmodel@nmat

  n_col <- lavmodel@nx.free
  m_el_idx <- x_el_idx <- vector("list", length = length(lavmodel@GLIST))
  for (mm in seq_along(lavmodel@GLIST)) {
    m_el_idx[[mm]] <- lavmodel@m.free.idx[[mm]]
    x_el_idx[[mm]] <- lavmodel@x.free.idx[[mm]]
    # handle symmetric matrices
    if (lavmodel@isSymmetric[mm]) {
      # since we use 'x.free.idx', only symmetric elements
      # are duplicated (not the equal ones, only in x.free.free)
      dix <- duplicated(x_el_idx[[mm]])
      if (any(dix)) {
        m_el_idx[[mm]] <- m_el_idx[[mm]][!dix]
        x_el_idx[[mm]] <- x_el_idx[[mm]][!dix]
      }
    }
  }

  # Delta per group (or block?)
  delta_1 <- vector("list", length = lavmodel@ngroups)

  for (g in 1:lavmodel@ngroups) {
    fsm_1 <- fsm[[g]]

    # which mm belong to group g?
    mm_in_group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    mlist <- lavmodel@GLIST[mm_in_group]

    nrow_scoffset <- ncol(mlist$lambda)
    nrow_scale <- ncol(mlist$lambda)
    n_row <- nrow_scoffset + nrow_scale
    delta_group <- matrix(0, nrow = n_row, ncol = n_col)

    # prepare some computations
    al_inv <- solve(fsm_1 %*% mlist$lambda)
    ata <- fsm_1 %*% mlist$theta %*% t(fsm_1)

    for (mm in mm_in_group) {
      mname <- names(lavmodel@GLIST)[mm]

      # skip empty ones
      if (!length(m_el_idx[[mm]])) next

      if (mname == "lambda") {
        d_l <- (-1 * (ata %*% al_inv + al_inv %*% ata) %*%
          (al_inv %x% al_inv) %*% fsm_1)

        delta_scoffset <- d_l
        delta_scale <- fsm_1 ## only ok for 1 row!!!
        delta <- rbind(delta_scoffset, delta_scale)

        delta_group[, x_el_idx[[mm]]] <- delta[, m_el_idx[[mm]]]
      } else if (mname == "theta") {
        d_t <- lav_mat_vec((t(al_inv) %*% fsm_1) %x%
          (t(fsm_1) %*% al_inv))
        delta_scoffset <- d_t
        delta_scale <- matrix(0,
          nrow = nrow_scale,
          ncol = length(mlist$theta)
        )
        delta <- rbind(delta_scoffset, delta_scale)

        delta_group[, x_el_idx[[mm]]] <- delta[, m_el_idx[[mm]]]
      } else if (mname %in% c("psi", "nu", "alpha")) {
        # zero
        next
      } else {
        lav_msg_stop(gettextf(
          "model matrix %s is not lambda/theta/psi", mname))
      }
    } # mm

    delta_1[[g]] <- delta_group
  } # g

  delta_1
}

lav_fsr_pa2si <- function(pt_1 = NULL, lvinfo) {
  # pt_orig <- pt_1

  # remove se column (if any)
  if (!is.null(pt_1$se)) {
    pt_1$se <- NULL
  }

  # ngroups
  ngroups <- lav_pt_ngroups(pt_1)

  lhs <- rhs <- op <- character(0)
  group <- block <- free <- exo <- integer(0)
  ustart <- est <- start <- numeric(0)

  for (g in seq_len(ngroups)) {
    n_mm <- length(lvinfo[[g]])
    for (mm in seq_len(n_mm)) {
      lvinfo_1 <- lvinfo[[g]][[mm]]
      lv_names <- lvinfo_1$lv.names

      nfac <- length(lv_names)
      if (nfac > 1L) {
        lav_msg_stop(gettext("more than 1 factor in measurement block"))
      }

      lv <- lv_names
      ind <- paste(lv, ".si", sep = "")
      scoffset <- lvinfo_1$scoffset[1, 1]
      scale <- lvinfo_1$scale[1, 1]

      lhs <- c(lhs, lv, ind, ind, ind)
      op <- c(op, "=~", "~~", "~*~", "~1")
      rhs <- c(rhs, ind, ind, ind, "")
      block <- c(block, rep(g, 4L))
      free <- c(free, 0L, 1L, 1L, 0L)
      ustart <- c(ustart, 1, scoffset, scale, 0)
      exo <- c(exo, rep(0L, 4L))
      group <- c(group, rep(g, 4L))
      start <- c(start, 1, scoffset, scale, 0)
      est <- c(est, 1, scoffset, scale, 0)
    }
  }

  # free counter
  idx_free <- which(free > 0)
  free[idx_free] <- max(pt_1$free) + seq_along(idx_free)

  list_1 <- list(
    id = max(pt_1$id) + seq_along(lhs),
    lhs = lhs,
    op = op,
    rhs = rhs,
    user = rep(10L, length(lhs)),
    block = block,
    group = group,
    level = rep(1L, length(lhs)),
    free = free,
    ustart = ustart,
    exo = exo,
    start = start,
    est = est
  )

  pt_si <- lav_pt_merge(pt_1, list_1)

  pt_si
}
