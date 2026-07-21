# lav_model utility functions

# initial version: YR 25/03/2009: `methods' for the Model class
# - YR 14 Jan 2014: rename object -> lavmodel, all functions as lav_model_*
# - YR 20 Nov 2021: add lav_model_dmmdpar

lav_model_get_parameters <- function(lavmodel = NULL, glist = NULL,   # nolint start
                                     type = "free", extra = TRUE, ...) {   # nolint end
  # type == "free": only non-redundant free parameters (x)
  # type == "user": all parameters listed in User model
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, NULL)
  # state or final?
  if (is.null(glist)) glist <- lavmodel@GLIST
                  # if this changes, tag @TDJorgensen in commit message

  if (type == "free") {
    n <- lavmodel@nx.free
    # } else if(type == "unco") {
    #    N <- lavmodel@nx.unco
  } else if (type == "user") {
    n <- lavmodel@nx.user
  }
  x <- numeric(n)

  for (mm in seq_along(lavmodel@GLIST)) {
    if (type == "free") {
      m_idx <- lavmodel@m.free.idx[[mm]]
      x_idx <- lavmodel@x.free.idx[[mm]]
      # } else if(type == "unco") {
      #    m.idx <- lavmodel@m.unco.idx[[mm]]
      #    x.idx <- lavmodel@x.unco.idx[[mm]]
    } else if (type == "user") {
      m_idx <- lavmodel@m.user.idx[[mm]]
      x_idx <- lavmodel@x.user.idx[[mm]]
    }
    x[x_idx] <- glist[[mm]][m_idx]
  }

  if (type == "user" && extra && sum(
    lavmodel@x.def.idx,
    lavmodel@x.ceq.idx,
    lavmodel@x.cin.idx
  ) > 0L) {
    # we need 'free' x
    x_free <- lav_model_get_parameters(
      lavmodel = lavmodel, glist = glist,
      type = "free"
    )
    if (length(lavmodel@x.def.idx) > 0L) {
      x[lavmodel@x.def.idx] <- lavmodel@def.function(x_free)
    }
    if (length(lavmodel@x.ceq.idx) > 0L) {
      x[lavmodel@x.ceq.idx] <- lavmodel@ceq.function(x_free)
    }
    if (length(lavmodel@x.cin.idx) > 0L) {
      tmp <- lavmodel@cin.function(x_free)
      # remove lower/upper bound values (if any)
      bound_idx <- attr(tmp, "bound.idx")
      if (length(bound_idx) > 0L) {
        tmp <- tmp[-bound_idx]
      }
      x[lavmodel@x.cin.idx] <- tmp
    }
  }

  x
}

# warning: this will make a copy of lavmodel
# warning: if categorical/correlation: 'delta' parameterization does
#          not work properly if we have 'mediators' (where x is not fixed)
#          that are observed (residuals are in PSI, and are not additive)
#          Note: fixed in 0.6-20 for recursive models
lav_model_set_parameters <- function(lavmodel = NULL, x = NULL) {
  tmp <- lavmodel@GLIST
  for (mm in seq_along(lavmodel@GLIST)) {
    m_free_idx <- lavmodel@m.free.idx[[mm]]
    x_free_idx <- lavmodel@x.free.idx[[mm]]
    tmp[[mm]][m_free_idx] <- x[x_free_idx]
  }

  correlation <- lavmodel@correlation

  # categorical? set categorical theta elements (if any)
  if (lavmodel@categorical || correlation) {
    nmat <- lavmodel@nmat
    if (lavmodel@representation == "LISREL") {
      for (g in 1:lavmodel@nblocks) {
        # which mm belong to group g?
        mm_in_group <- 1:nmat[g] + cumsum(c(0L, nmat))[g]

        if (lavmodel@estimator %in% c("MML", "FML")) {
          #  ttt <- diag(tmp[mm.in.group]$theta)
          #  diag(tmp[mm.in.group]$theta) <- as.numeric(NA)
          #  if(length(lavmodel@num.idx[[g]]) > 0L) {
          #      diag(tmp[mm.in.group]$theta)[ lavmodel@num.idx[[g]] ] <-
          #          ttt[ lavmodel@num.idx[[g]] ]
          #  }
        } else {
          if (lavmodel@parameterization == "delta") {
            tmp[mm_in_group] <-
              lav_lisrel_residual_variances(
                mlist = tmp[mm_in_group],
                num_idx = lavmodel@num.idx[[g]],
                ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[g]],
                ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[g]],
                marginal = correlation && lavmodel@conditional.x
              )
          } else if (lavmodel@parameterization == "theta") {
            tmp[mm_in_group] <-
              lav_lisrel_delta(
                mlist = tmp[mm_in_group],
                num_idx = lavmodel@num.idx[[g]]
              )
          }
        }
      }
    } else {
      cat("FIXME: deal with theta elements in the categorical case (RAM)")
    }
  }

  if (lavmodel@composites) {
    # for package stdmod only! (vignette stdmod_lavaan uses old fit object)
    # if (.hasSlot(lavmodel, "composites") && lavmodel@composites) {
    nmat <- lavmodel@nmat
    if (lavmodel@representation == "LISREL") {
      for (g in 1:lavmodel@nblocks) {
        # which mm belong to group g?
        mm_in_group <- 1:nmat[g] + cumsum(c(0L, nmat))[g]

        tmp[mm_in_group] <-
          lav_lisrel_comp_set_intresvar(mlist = tmp[mm_in_group])
      }
    } else {
      cat("FIXME: deal with Composites if representation = RAM")
    }
  }

  lavmodel@GLIST <- tmp

  lavmodel
}

# create a standalone GLIST, filled with (new) x values
# (avoiding a copy of lavmodel)
lav_model_x2glist <- function(lavmodel = NULL, x = NULL,
                              type = "free", set_delta = TRUE,
                              m_el_idx = NULL, x_el_idx = NULL) {
  correlation <- lavmodel@correlation

  glist <- lavmodel@GLIST
  for (mm in seq_along(glist)) {
    # skip empty matrix
    if (nrow(glist[[mm]]) == 0L) {
      next
    }
    if (type == "free") {
      m_el_idx_1 <- lavmodel@m.free.idx[[mm]]
      x_el_idx_1 <- lavmodel@x.free.idx[[mm]]
    } else if (type == "unco") {
      m_el_idx_1 <- lavmodel@m.free.idx[[mm]]
      x_el_idx_1 <- lavmodel@x.unco.idx[[mm]]
    } else if (type == "full") {
      if (lavmodel@isSymmetric[mm]) {
        n <- ncol(glist[[mm]])
        m_el_idx_1 <- lav_mat_vech_idx(n)
      } else {
        m_el_idx_1 <- seq_along(glist[[mm]])
      }
      x_el_idx_1 <- seq_along(m_el_idx)
      if (mm > 1) x_el_idx_1 <- x_el_idx_1 + sum(lavmodel@mmSize[1:(mm - 1)])
    } else if (type == "custom") {
      # nothing to do, m.el.idx and x.el.idx should be given
      m_el_idx_1 <- m_el_idx[[mm]]
      x_el_idx_1 <- x_el_idx[[mm]]
    }

    # assign
    glist[[mm]][m_el_idx_1] <- x[x_el_idx_1]

    # make symmetric (if full)
    if (type == "full" && lavmodel@isSymmetric[mm]) {
      t_1 <- t(glist[[mm]])
      glist[[mm]][upper.tri(glist[[mm]])] <- t_1[upper.tri(t_1)]
    }
  }

  #    # theta parameterization: delta must be reset!
  #    if((lavmodel@categorical || correlation) && setDelta &&
  #       lavmodel@parameterization == "theta") {
  #        nmat <- lavmodel@nmat
  #        for(g in 1:lavmodel@nblocks) {
  #            # which mm belong to group g?
  #            mm.in.group <- 1:nmat[g] + cumsum(c(0L,nmat))[g]
  #            GLIST[mm.in.group] <-
  #                lav_lisrel_delta(MLIST = GLIST[mm.in.group],
  #                    num.idx = lavmodel@num.idx[[g]])
  #        }
  #    }

  # in 0.6-13: we always set theta/delta
  if ((lavmodel@categorical || correlation) && set_delta) {
    nmat <- lavmodel@nmat
    if (lavmodel@representation == "LISREL") {
      for (g in 1:lavmodel@nblocks) {
        # which mm belong to group g?
        mm_in_group <- 1:nmat[g] + cumsum(c(0L, nmat))[g]

        if (lavmodel@parameterization == "delta") {
          glist[mm_in_group] <-
            lav_lisrel_residual_variances(
              mlist = glist[mm_in_group],
              num_idx = lavmodel@num.idx[[g]],
              ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[g]],
              ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[g]],
              marginal = correlation && lavmodel@conditional.x
            )
        } else if (lavmodel@parameterization == "theta") {
          glist[mm_in_group] <-
            lav_lisrel_delta(
              mlist = glist[mm_in_group],
              num_idx = lavmodel@num.idx[[g]]
            )
        }
      } # blocks
    } else {
      cat("FIXME: deal with theta elements in the categorical case (RAM)")
    }
  }

  if (lavmodel@composites) {
    nmat <- lavmodel@nmat
    if (lavmodel@representation == "LISREL") {
      for (g in 1:lavmodel@nblocks) {
        # which mm belong to group g?
        mm_in_group <- 1:nmat[g] + cumsum(c(0L, nmat))[g]

        glist[mm_in_group] <-
          lav_lisrel_comp_set_intresvar(mlist = glist[mm_in_group])
      }
    } else {
      cat("FIXME: deal with Composites when representation = RAM")
    }
  }

  glist
}

# backwards compatibility
# getModelParameters <- lav_model_get_parameters
# setModelParameters <- lav_model_set_parameters
# x2GLIST            <- lav_model_x2glist
