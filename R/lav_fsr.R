# compute the jacobian: dtheta_2/dtheta_1:
#
# theta_2: - in the rows
#          - the croon corrections, expressed as
#              1) scaled offsets (scoffset), and
#              2) scaling factors
# theta_1: - in the columns
#          - the free parameters of the measurement model
#
lav_fsr_delta21 <- function(object, FSM = NULL) {
  lavmodel <- object@Model
  nmat <- lavmodel@nmat

  NCOL <- lavmodel@nx.free
  m.el.idx <- x.el.idx <- vector("list", length = length(lavmodel@GLIST))
  for (mm in seq_len(length(lavmodel@GLIST))) {
    m.el.idx[[mm]] <- lavmodel@m.free.idx[[mm]]
    x.el.idx[[mm]] <- lavmodel@x.free.idx[[mm]]
    # handle symmetric matrices
    if (lavmodel@isSymmetric[mm]) {
      # since we use 'x.free.idx', only symmetric elements
      # are duplicated (not the equal ones, only in x.free.free)
      dix <- duplicated(x.el.idx[[mm]])
      if (any(dix)) {
        m.el.idx[[mm]] <- m.el.idx[[mm]][!dix]
        x.el.idx[[mm]] <- x.el.idx[[mm]][!dix]
      }
    }
  }

  # Delta per group (or block?)
  Delta <- vector("list", length = lavmodel@ngroups)

  for (g in 1:lavmodel@ngroups) {
    fsm <- FSM[[g]]

    # which mm belong to group g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- lavmodel@GLIST[mm.in.group]

    nrow.scoffset <- ncol(MLIST$lambda)
    nrow.scale <- ncol(MLIST$lambda)
    NROW <- nrow.scoffset + nrow.scale
    Delta.group <- matrix(0, nrow = NROW, ncol = NCOL)

    # prepare some computations
    AL.inv <- solve(fsm %*% MLIST$lambda)
    ATA <- fsm %*% MLIST$theta %*% t(fsm)

    for (mm in mm.in.group) {
      mname <- names(lavmodel@GLIST)[mm]

      # skip empty ones
      if (!length(m.el.idx[[mm]])) next

      if (mname == "lambda") {
        dL <- (-1 * (ATA %*% AL.inv + AL.inv %*% ATA) %*%
          (AL.inv %x% AL.inv) %*% fsm)

        delta.scoffset <- dL
        delta.scale <- fsm ## only ok for 1 row!!!
        delta <- rbind(delta.scoffset, delta.scale)

        Delta.group[, x.el.idx[[mm]]] <- delta[, m.el.idx[[mm]]]
      } else if (mname == "theta") {
        dT <- lav_matrix_vec((t(AL.inv) %*% fsm) %x%
          (t(fsm) %*% AL.inv))
        delta.scoffset <- dT
        delta.scale <- matrix(0,
          nrow = nrow.scale,
          ncol = length(MLIST$theta)
        )
        delta <- rbind(delta.scoffset, delta.scale)

        Delta.group[, x.el.idx[[mm]]] <- delta[, m.el.idx[[mm]]]
      } else if (mname %in% c("psi", "nu", "alpha")) {
        # zero
        next
      } else {
        lav_msg_stop(gettextf(
          "model matrix %s is not lambda/theta/psi", mname))
      }
    } # mm

    Delta[[g]] <- Delta.group
  } # g

  Delta
}

lav_fsr_pa2si <- function(PT = NULL, LVINFO) {
  PT.orig <- PT

  # remove se column (if any)
  if (!is.null(PT$se)) {
    PT$se <- NULL
  }

  # ngroups
  ngroups <- lav_partable_ngroups(PT)

  lhs <- rhs <- op <- character(0)
  group <- block <- level <- free <- exo <- integer(0)
  ustart <- est <- start <- numeric(0)

  for (g in seq_len(ngroups)) {
    nMM <- length(LVINFO[[g]])
    for (mm in seq_len(nMM)) {
      lvinfo <- LVINFO[[g]][[mm]]
      lv.names <- lvinfo$lv.names

      nfac <- length(lv.names)
      if (nfac > 1L) {
        lav_msg_stop(gettext("more than 1 factor in measurement block"))
      }

      LV <- lv.names
      ind <- paste(LV, ".si", sep = "")
      scoffset <- lvinfo$scoffset[1, 1]
      scale <- lvinfo$scale[1, 1]

      lhs <- c(lhs, LV, ind, ind, ind)
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

  # ree counter
  idx.free <- which(free > 0)
  free[idx.free] <- max(PT$free) + 1:length(idx.free)

  LIST <- list(
    id = max(PT$id) + 1:length(lhs),
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

  PT.si <- lav_partable_merge(PT, LIST)

  PT.si
}
