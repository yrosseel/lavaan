lav_standardize_lv_x <- function(x, lavobject, partable = NULL, cov.std = TRUE,
                                 lv.var = NULL,
                                 rotation = FALSE) {
  # set new values for x
  lavmodel <- lav_model_set_parameters(lavmodel = lavobject@Model, x = x)

  if (rotation) {
    x.unrotated <- x
    lavmodel@GLIST <- lavTech(lavobject, "est.unrotated") # unrotated!
    est.rot <- lav_model_efa_rotate_x(
      x = x.unrotated,
      lavmodel = lavmodel, # unrotated!
      lavoptions = lavobject@Options,
      init.rot = lavmodel@H,
      type = "user",
      extra = TRUE
    )
    GLIST <- attr(est.rot, "extra")$GLIST
    attributes(est.rot) <- NULL
    est <- est.rot
  } else {
    GLIST <- lavmodel@GLIST # if this changes, tag @TDJorgensen in commit message
    est <- lav_model_get_parameters(lavmodel, type = "user")
  }

  x.stand.user <- lav_standardize_lv(
    lavobject = lavobject,
    partable = partable, est = est,
    GLIST = GLIST, cov.std = cov.std,
    lv.var = lv.var
  )

  x.stand.user
}

lav_standardize_all_x <- function(x, lavobject, partable = NULL, cov.std = TRUE,
                                  rotation = FALSE) {
  lavmodel <- lav_model_set_parameters(lavmodel = lavobject@Model, x = x)

  if (rotation) {
    x.unrotated <- x
    lavmodel@GLIST <- lavTech(lavobject, "est.unrotated") # unrotated!
    est.rot <- lav_model_efa_rotate_x(
      x = x.unrotated,
      lavmodel = lavmodel, # unrotated!
      lavoptions = lavobject@Options,
      init.rot = lavmodel@H,
      type = "user",
      extra = TRUE
    )
    GLIST <- attr(est.rot, "extra")$GLIST
    attributes(est.rot) <- NULL
    est <- est.rot
  } else {
    GLIST <- lavmodel@GLIST # if this changes, tag @TDJorgensen in commit message
    est <- lav_model_get_parameters(lavmodel, type = "user")
  }

  x.stand.user <- lav_standardize_all(
    lavobject = lavobject,
    partable = partable, est = est,
    est.std = NULL, GLIST = GLIST,
    cov.std = cov.std
  )
  x.stand.user
}

lav_standardize_all_nox_x <- function(x, lavobject, partable = NULL,
                                      cov.std = TRUE, rotation = FALSE) {
  lavmodel <- lav_model_set_parameters(lavmodel = lavobject@Model, x = x)

  if (rotation) {
    x.unrotated <- x
    lavmodel@GLIST <- lavTech(lavobject, "est.unrotated") # unrotated!
    est.rot <- lav_model_efa_rotate_x(
      x = x.unrotated,
      lavmodel = lavmodel, # unrotated!
      lavoptions = lavobject@Options,
      init.rot = lavmodel@H,
      type = "user",
      extra = TRUE
    )
    GLIST <- attr(est.rot, "extra")$GLIST
    attributes(est.rot) <- NULL
    est <- est.rot
  } else {
    GLIST <- lavmodel@GLIST # if this changes, tag @TDJorgensen in commit message
    est <- lav_model_get_parameters(lavmodel, type = "user")
  }

  x.stand.user <- lav_standardize_all_nox(
    lavobject = lavobject,
    partable = partable, est = est,
    est.std = NULL, GLIST = GLIST,
    cov.std = cov.std
  )
  x.stand.user
}

lav_unstandardize_ov_x <- function(x, lavobject) {
  partable <- lavobject@ParTable
  partable$ustart <- x
  lav_unstandardize_ov(
    partable = partable,
    ov.var = lavobject@SampleStats@var,
    cov.std = TRUE
  )
}


lav_standardize_lv <- function(lavobject = NULL,
                               partable = NULL, est = NULL, GLIST = NULL,
                               cov.std = TRUE, lv.var = NULL,
                               lavmodel = NULL, lavpartable = NULL) {
  if (is.null(lavobject)) {
    stopifnot(!is.null(lavmodel))
    stopifnot(!is.null(lavpartable))
    if (is.null(est)) {
      if (!is.null(lavpartable$est)) {
        est <- lavpartable$est # if this changes, tag @TDJorgensen in commit message
      } else {
        lav_msg_stop(gettext("could not find `est' in lavpartable"))
      }
    }
  } else {
    lavmodel <- lavobject@Model
    lavpartable <- lavobject@ParTable
    if (is.null(est)) {
      est <- lav_object_inspect_est(lavobject)
    }
  }

  if (is.null(partable)) {
    partable <- lavpartable
  }
  if (is.null(GLIST)) {
    GLIST <- lavmodel@GLIST
  }

  out <- est
  N <- length(est)
  stopifnot(N == length(partable$lhs))

  nmat <- lavmodel@nmat

  # compute ETA
  if (is.null(lv.var)) {
    LV.ETA <- computeVETA(
      lavmodel = lavmodel,
      GLIST = GLIST
    )
  }

  for (g in 1:lavmodel@nblocks) {
    ov.names <- vnames(lavpartable, "ov", block = g) # not user,
    # which may be incomplete
    lv.names <- vnames(lavpartable, "lv", block = g)

    # shortcut: no latents in this block, nothing to do
    if (length(lv.names) == 0L) {
      next
    }

    # which mm belong to block g?
    mm.in.group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    MLIST <- GLIST[mm.in.group]

    if (is.null(lv.var)) {
      ETA2 <- diag(LV.ETA[[g]])
    } else {
      ETA2 <- lv.var[[g]]
    }
    # change negative values to NA
    ETA2[ETA2 < 0] <- as.numeric(NA)
    ETA <- sqrt(ETA2)

    # 1a. "=~" regular indicators
    idx <- which(partable$op == "=~" & !(partable$rhs %in% lv.names) &
      partable$block == g)
    out[idx] <- out[idx] * ETA[match(partable$lhs[idx], lv.names)]

    # 1b. "=~" regular higher-order lv indicators
    idx <- which(partable$op == "=~" & !(partable$rhs %in% ov.names) &
      partable$block == g)
    out[idx] <- (out[idx] * ETA[match(partable$lhs[idx], lv.names)]
      / ETA[match(partable$rhs[idx], lv.names)])

    # 1c. "=~" indicators that are both in ov and lv
    # idx <- which(partable$op == "=~" & partable$rhs %in% ov.names
    #                             & partable$rhs %in% lv.names &
    #             partable$block == g)

    # 2. "~" regressions (and "<~")
    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$lhs %in% lv.names &
      partable$block == g)
    out[idx] <- out[idx] / ETA[match(partable$lhs[idx], lv.names)]

    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$rhs %in% lv.names &
      partable$block == g)
    out[idx] <- out[idx] * ETA[match(partable$rhs[idx], lv.names)]

    # 3a. "~~" ov
    # idx <- which(partable$op == "~~" & !(partable$lhs %in% lv.names) &
    #             partable$block == g)

    # 3b. "~~" lv
    # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances
    #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
    #            were i.var and j.var where diagonal elements of ETA
    #
    #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
    #            elements are the 'PSI' diagonal elements!!

    # variances
    rv.idx <- which(partable$op == "~~" & partable$rhs %in% lv.names &
      partable$lhs == partable$rhs &
      partable$block == g)
    out[rv.idx] <- (out[rv.idx] / ETA[match(partable$lhs[rv.idx], lv.names)]
      / ETA[match(partable$rhs[rv.idx], lv.names)])

    # covariances lv
    # three types:
    # - only lhs is LV (and fixed.x = FALSE)
    # - only rhs is LV (and fixed.x = FALSE)
    # - both lhs and rhs are LV (regular case)
    if (cov.std) {
      if (!is.complex(est[rv.idx])) {
        RV <- sqrt(abs(est[rv.idx])) # abs in case of heywood cases
      } else {
        RV <- sqrt(est[rv.idx])
      }
      rv.names <- partable$lhs[rv.idx]
    }

    # left
    idx.lhs <- which(partable$op == "~~" &
      partable$lhs %in% lv.names &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx.lhs) > 0L) {
      if (cov.std == FALSE) {
        out[idx.lhs] <-
          (out[idx.lhs] / ETA[match(partable$lhs[idx.lhs], lv.names)])
      } else {
        out[idx.lhs] <-
          (out[idx.lhs] / RV[match(partable$lhs[idx.lhs], rv.names)])
      }
    }

    # right
    idx.rhs <- which(partable$op == "~~" &
      partable$rhs %in% lv.names &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx.rhs) > 0L) {
      if (cov.std == FALSE) {
        out[idx.rhs] <-
          (out[idx.rhs] / ETA[match(partable$rhs[idx.rhs], lv.names)])
      } else {
        out[idx.rhs] <-
          (out[idx.rhs] / RV[match(partable$rhs[idx.rhs], rv.names)])
      }
    }


    # 4a. "~1" ov
    # idx <- which(partable$op == "~1" & !(partable$lhs %in% lv.names) &
    #             partable$block == g)

    # 4b. "~1" lv
    idx <- which(partable$op == "~1" & partable$lhs %in% lv.names &
      partable$block == g)
    out[idx] <- out[idx] / ETA[match(partable$lhs[idx], lv.names)]
  }

  # 5a ":="
  idx <- which(partable$op == ":=")
  if (length(idx) > 0L) {
    x <- out[partable$free & !duplicated(partable$free)]
    out[idx] <- lavmodel@def.function(x)
  }

  # 5b "=="
  idx <- which(partable$op == "==")
  if (length(idx) > 0L) {
    x <- out[partable$free & !duplicated(partable$free)]
    out[idx] <- lavmodel@ceq.function(x)
  }

  # 5c. "<" or ">"
  idx <- which((partable$op == "<" | partable$op == ">"))
  if (length(idx) > 0L) {
    x <- out[partable$free & !duplicated(partable$free)]
    out[idx] <- lavmodel@cin.function(x)
  }

  out
}

lav_standardize_all <- function(lavobject = NULL,
                                partable = NULL, est = NULL, est.std = NULL,
                                GLIST = NULL, cov.std = TRUE, ov.var = NULL,
                                lv.var = NULL,
                                lavmodel = NULL, lavpartable = NULL,
                                cov.x = NULL) {
  if (is.null(lavobject)) {
    stopifnot(!is.null(lavmodel))
    stopifnot(!is.null(lavpartable))
    if (is.null(est)) {
      if (!is.null(lavpartable$est)) {
        est <- lavpartable$est # if this changes, tag @TDJorgensen in commit message
      } else {
        lav_msg_stop(gettext("could not find `est' in lavpartable"))
      }
    }
  } else {
    lavmodel <- lavobject@Model
    lavpartable <- lavobject@ParTable
    if (is.null(est)) {
      est <- lav_object_inspect_est(lavobject)
    }
    if (lavmodel@conditional.x) {
      if (is.null(cov.x)) {
        # try SampleStats slot
        # if("SampleStats" %in% slotNames(lavobject)) {
        #    cov.x <- lavobject@SampleStats@cov.x
        if (!is.null(lavobject@implied$cov.x[[1]])) {
          cov.x <- lavobject@implied$cov.x # if this changes, tag @TDJorgensen in commit message
        } else {
          # perhaps lavaanList object
          # extract it from GLIST per block
          cov.x <- vector("list", length = lavmodel@nblocks)
          for (b in seq_len(lavmodel@nblocks)) {
            # which mm belong to block b?
            mm.in.block <- (seq_len(lavmodel@nmat[b]) +
              cumsum(c(0, lavmodel@nmat))[b])
            MLIST <- lavmodel@GLIST[mm.in.block]
            cov.x[[b]] <- MLIST[["cov.x"]]
          }
        }
      }
    }
  }

  if (is.null(partable)) {
    partable <- lavpartable
  }
  if (is.null(GLIST)) {
    GLIST <- lavmodel@GLIST
  }
  if (is.null(est.std)) {
    est.std <- lav_standardize_lv(
      lavobject = lavobject,
      partable = partable, est = est, GLIST = GLIST,
      cov.std = cov.std, lv.var = lv.var, lavmodel = lavmodel,
      lavpartable = lavpartable
    )
  }

  out <- est.std
  N <- length(est.std)
  stopifnot(N == length(partable$lhs))

  VY <- computeVY(
    lavmodel = lavmodel, GLIST = GLIST,
    diagonal.only = TRUE
  )


  for (g in 1:lavmodel@nblocks) {
    ov.names <- vnames(lavpartable, "ov", block = g) # not user
    lv.names <- vnames(lavpartable, "lv", block = g)

    if (is.null(ov.var)) {
      OV2 <- VY[[g]]
      # replace zero values by NA (but keep negative values)
      zero.idx <- which(abs(OV2) < .Machine$double.eps)
      if (length(zero.idx) > 0L) {
        OV2[zero.idx] <- as.numeric(NA)
      }

      # replace negative values by NA (for sqrt)
      tmp.OV2 <- OV2
      neg.idx <- which(tmp.OV2 < 0)
      if (length(neg.idx) > 0L) {
        tmp.OV2[neg.idx] <- as.numeric(NA)
      }
      OV <- sqrt(tmp.OV2)
    } else {
      OV2 <- ov.var[[g]]
      OV <- sqrt(OV2)
    }

    if (lavmodel@conditional.x) {
      # extend OV with ov.names.x
      ov.names.x <- vnames(lavpartable, "ov.x", block = g)
      ov.names.nox <- vnames(lavpartable, "ov.nox", block = g)
      ov.names <- c(ov.names.nox, ov.names.x)
      OV2 <- c(OV2, diag(cov.x[[g]]))
      OV <- c(OV, sqrt(diag(cov.x[[g]])))
    }

    # 1a. "=~" regular indicators
    idx <- which(partable$op == "=~" & !(partable$rhs %in% lv.names) &
      partable$block == g)
    out[idx] <- out[idx] / OV[match(partable$rhs[idx], ov.names)]

    # 1b. "=~" regular higher-order lv indicators

    # 1c. "=~" indicators that are both in ov and lv
    # idx <- which(partable$op == "=~" & partable$rhs %in% ov.names
    #                             & partable$rhs %in% lv.names &
    #             partable$block == g)

    # 2. "~" regressions (and "<~")
    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$lhs %in% ov.names &
      partable$block == g)
    out[idx] <- out[idx] / OV[match(partable$lhs[idx], ov.names)]

    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$rhs %in% ov.names &
      partable$block == g)
    out[idx] <- out[idx] * OV[match(partable$rhs[idx], ov.names)]

    # 3a. "~~" ov
    # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances
    #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
    #            were i.var and j.var where diagonal elements of OV
    #
    #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
    #            elements are the 'THETA' diagonal elements!!

    # variances
    rv.idx <- which(partable$op == "~~" & !(partable$lhs %in% lv.names) &
      partable$lhs == partable$rhs &
      partable$block == g)
    # out[rv.idx] <- ( out[rv.idx] / OV[ match(partable$lhs[rv.idx], ov.names) ]
    #                             / OV[ match(partable$rhs[rv.idx], ov.names) ] )
    out[rv.idx] <- (out[rv.idx] /
      OV2[match(partable$lhs[rv.idx], ov.names)])

    # covariances ov
    # three types:
    # - only lhs is OV (and fixed.x = FALSE)
    # - only rhs is OV (and fixed.x = FALSE)
    # - both lhs and rhs are OV (regular case)
    if (cov.std) {
      if (!is.complex(est[rv.idx])) {
        RV <- sqrt(abs(est[rv.idx]))
      } else {
        RV <- sqrt(est[rv.idx])
      }
      rv.names <- partable$lhs[rv.idx]
    }

    # left
    idx.lhs <- which(partable$op == "~~" &
      !(partable$lhs %in% lv.names) &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx.lhs) > 0L) {
      if (cov.std == FALSE) {
        out[idx.lhs] <-
          (out[idx.lhs] / OV[match(partable$lhs[idx.lhs], ov.names)])
      } else {
        out[idx.lhs] <-
          (out[idx.lhs] / RV[match(partable$lhs[idx.lhs], rv.names)])
      }
    }

    # right
    idx.rhs <- which(partable$op == "~~" &
      !(partable$rhs %in% lv.names) &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx.rhs) > 0L) {
      if (cov.std == FALSE) {
        out[idx.rhs] <-
          (out[idx.rhs] / OV[match(partable$rhs[idx.rhs], ov.names)])
      } else {
        out[idx.rhs] <-
          (out[idx.rhs] / RV[match(partable$rhs[idx.rhs], rv.names)])
      }
    }

    # 3b. "~~" lv
    # idx <- which(partable$op == "~~" & partable$rhs %in% lv.names &
    #             partable$block == g)

    # 4a. "~1" ov
    idx <- which(partable$op == "~1" & !(partable$lhs %in% lv.names) &
      partable$block == g)
    out[idx] <- out[idx] / OV[match(partable$lhs[idx], ov.names)]

    # 4b. "~1" lv
    # idx <- which(partable$op == "~1" & partable$lhs %in% lv.names &
    #             partable$block == g)

    # 4c. "|" thresholds
    idx <- which(partable$op == "|" & !(partable$lhs %in% lv.names) &
      partable$block == g)
    out[idx] <- out[idx] / OV[match(partable$lhs[idx], ov.names)]

    # 4d. "~*~" scales
    idx <- which(partable$op == "~*~" & !(partable$lhs %in% lv.names) &
      partable$block == g)
    out[idx] <- 1.0
  }

  # 5a ":="
  idx <- which(partable$op == ":=")
  if (length(idx) > 0L) {
    x <- out[partable$free & !duplicated(partable$free)]
    out[idx] <- lavmodel@def.function(x)
  }

  # 5b "=="
  idx <- which(partable$op == "==")
  if (length(idx) > 0L) {
    x <- out[partable$free & !duplicated(partable$free)]
    out[idx] <- lavmodel@ceq.function(x)
  }

  # 5c. "<" or ">"
  idx <- which((partable$op == "<" | partable$op == ">"))
  if (length(idx) > 0L) {
    x <- out[partable$free & !duplicated(partable$free)]
    out[idx] <- lavmodel@cin.function(x)
  }

  out
}


lav_standardize_all_nox <- function(lavobject = NULL,
                                    partable = NULL, est = NULL, est.std = NULL,
                                    GLIST = NULL, cov.std = TRUE, ov.var = NULL,
                                    lv.var = NULL,
                                    lavmodel = NULL, lavpartable = NULL,
                                    cov.x = NULL) {
  if (is.null(lavobject)) {
    stopifnot(!is.null(lavmodel))
    stopifnot(!is.null(lavpartable))
    if (is.null(est)) {
      if (!is.null(lavpartable$est)) {
        est <- lavpartable$est # if this changes, tag @TDJorgensen in commit message
      } else {
        lav_msg_stop(gettext("could not find `est' in lavpartable"))
      }
    }
  } else {
    lavmodel <- lavobject@Model
    lavpartable <- lavobject@ParTable
    if (is.null(est)) {
      est <- lav_object_inspect_est(lavobject)
    }
    if (lavmodel@conditional.x) {
      if (is.null(cov.x)) {
        # try SampleStats slot
        # if("SampleStats" %in% slotNames(lavobject)) {
        #    cov.x <- lavobject@SampleStats@cov.x
        if (!is.null(lavobject@implied$cov.x[[1]])) {
          cov.x <- lavobject@implied$cov.x # if this changes, tag @TDJorgensen in commit message
        } else {
          # perhaps lavaanList object
          # extract it from GLIST per block
          cov.x <- vector("list", length = lavmodel@nblocks)
          for (b in seq_len(lavmodel@nblocks)) {
            # which mm belong to block b?
            mm.in.block <- (seq_len(lavmodel@nmat[b]) +
              cumsum(c(0, lavmodel@nmat))[b])
            MLIST <- lavmodel@GLIST[mm.in.block]
            cov.x[[b]] <- MLIST[["cov.x"]]
          }
        }
      }
    }
  }

  if (is.null(partable)) {
    partable <- lavpartable
  }
  if (is.null(GLIST)) {
    GLIST <- lavmodel@GLIST
  }
  if (is.null(est.std)) {
    est.std <- lav_standardize_lv(
      lavobject = lavobject,
      partable = partable, est = est, GLIST = GLIST,
      cov.std = cov.std, lv.var = lv.var, lavmodel = lavmodel,
      lavpartable = lavpartable
    )
  }


  out <- est.std
  N <- length(est.std)
  stopifnot(N == length(partable$lhs))

  VY <- computeVY(
    lavmodel = lavmodel, GLIST = GLIST,
    diagonal.only = TRUE
  )


  for (g in 1:lavmodel@nblocks) {
    ov.names <- vnames(lavpartable, "ov", block = g)
    ov.names.x <- vnames(lavpartable, "ov.x", block = g)
    ov.names.nox <- vnames(lavpartable, "ov.nox", block = g)
    lv.names <- vnames(lavpartable, "lv", block = g)

    if (is.null(ov.var)) {
      OV2 <- VY[[g]]
      # replace zero values by NA (but keep negative values)
      zero.idx <- which(abs(OV2) < .Machine$double.eps)
      if (length(zero.idx) > 0L) {
        OV2[zero.idx] <- as.numeric(NA)
      }

      # replace negative values by NA (for sqrt)
      tmp.OV2 <- OV2
      neg.idx <- which(tmp.OV2 < 0)
      if (length(neg.idx) > 0L) {
        tmp.OV2[neg.idx] <- as.numeric(NA)
      }
      OV <- sqrt(tmp.OV2)
    } else {
      OV2 <- ov.var[[g]]
      OV <- sqrt(OV2)
    }


    if (lavmodel@conditional.x) {
      # extend OV with ov.names.x
      ov.names.x <- vnames(lavpartable, "ov.x", block = g)
      ov.names <- c(ov.names.nox, ov.names.x)
      OV2 <- c(OV2, diag(cov.x[[g]]))
      OV <- c(OV, sqrt(diag(cov.x[[g]])))
    }

    # 1a. "=~" regular indicators
    idx <- which(partable$op == "=~" & !(partable$rhs %in% lv.names) &
      partable$block == g)
    out[idx] <- out[idx] / OV[match(partable$rhs[idx], ov.names)]

    # 1b. "=~" regular higher-order lv indicators

    # 1c. "=~" indicators that are both in ov and lv
    # idx <- which(partable$op == "=~" & partable$rhs %in% ov.names
    #                             & partable$rhs %in% lv.names &
    #             partable$block == g)

    # 2. "~" regressions (and "<~")
    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$lhs %in% ov.names &
      partable$block == g)
    out[idx] <- out[idx] / OV[match(partable$lhs[idx], ov.names)]

    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$rhs %in% ov.names.nox &
      partable$block == g)
    out[idx] <- out[idx] * OV[match(partable$rhs[idx], ov.names.nox)]

    # 3a. "~~" ov
    # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances
    #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
    #            were i.var and j.var where diagonal elements of OV
    #
    #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
    #            elements are the 'THETA' diagonal elements!!

    # variances
    rv.idx <- which(partable$op == "~~" & !(partable$lhs %in% lv.names) &
      !(partable$lhs %in% ov.names.x) &
      partable$lhs == partable$rhs &
      partable$block == g)
    # out[rv.idx] <- ( out[rv.idx] / OV[ match(partable$lhs[rv.idx], ov.names) ]
    #                            / OV[ match(partable$rhs[rv.idx], ov.names) ] )
    out[rv.idx] <- (out[rv.idx] /
      OV2[match(partable$lhs[rv.idx], ov.names)])

    # covariances ov
    # three types:
    # - only lhs is OV (and fixed.x = FALSE)
    # - only rhs is OV (and fixed.x = FALSE)
    # - both lhs and rhs are OV (regular case)
    if (cov.std) {
      if (!is.complex(est[rv.idx])) {
        RV <- sqrt(abs(est[rv.idx]))
      } else {
        RV <- sqrt(est[rv.idx])
      }
      rv.names <- partable$lhs[rv.idx]
    }

    # left
    idx.lhs <- which(partable$op == "~~" &
      !(partable$lhs %in% lv.names) &
      !(partable$lhs %in% ov.names.x) &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx.lhs) > 0L) {
      if (cov.std == FALSE) {
        out[idx.lhs] <-
          (out[idx.lhs] / OV[match(partable$lhs[idx.lhs], ov.names)])
      } else {
        out[idx.lhs] <-
          (out[idx.lhs] / RV[match(partable$lhs[idx.lhs], rv.names)])
      }
    }

    # right
    idx.rhs <- which(partable$op == "~~" &
      !(partable$rhs %in% lv.names) &
      !(partable$rhs %in% ov.names.x) &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx.rhs) > 0L) {
      if (cov.std == FALSE) {
        out[idx.rhs] <-
          (out[idx.rhs] / OV[match(partable$rhs[idx.rhs], ov.names)])
      } else {
        out[idx.rhs] <-
          (out[idx.rhs] / RV[match(partable$rhs[idx.rhs], rv.names)])
      }
    }

    # 3b. "~~" lv
    # idx <- which(partable$op == "~~" & partable$rhs %in% lv.names &
    #             partable$block == g)

    # 4a. "~1" ov
    idx <- which(partable$op == "~1" & !(partable$lhs %in% lv.names) &
      !(partable$lhs %in% ov.names.x) &
      partable$block == g)
    out[idx] <- out[idx] / OV[match(partable$lhs[idx], ov.names)]

    # 4b. "~1" lv
    # idx <- which(partable$op == "~1" & partable$lhs %in% lv.names &
    #             partable$block == g)

    # 4c. "|" thresholds
    idx <- which(partable$op == "|" & !(partable$lhs %in% lv.names) &
      partable$block == g)
    out[idx] <- out[idx] / OV[match(partable$lhs[idx], ov.names)]

    # 4d. "~*~" scales
    idx <- which(partable$op == "~*~" & !(partable$lhs %in% lv.names) &
      partable$block == g)
    out[idx] <- 1.0
  }

  # 5a ":="
  idx <- which(partable$op == ":=")
  if (length(idx) > 0L) {
    x <- out[partable$free & !duplicated(partable$free)]
    out[idx] <- lavmodel@def.function(x)
  }

  # 5b "=="
  idx <- which(partable$op == "==")
  if (length(idx) > 0L) {
    x <- out[partable$free & !duplicated(partable$free)]
    out[idx] <- lavmodel@ceq.function(x)
  }

  # 5c. "<" or ">"
  idx <- which((partable$op == "<" | partable$op == ">"))
  if (length(idx) > 0L) {
    x <- out[partable$free & !duplicated(partable$free)]
    out[idx] <- lavmodel@cin.function(x)
  }

  out
}

lav_unstandardize_ov <- function(partable, ov.var = NULL, cov.std = TRUE) {
  # check if ustart is missing; if so, look for est
  if (is.null(partable$ustart)) {
    partable$ustart <- partable$est
  }

  # check if block is missing
  if (is.null(partable$block)) {
    partable$block <- rep(1L, length(partable$ustart))
  }

  stopifnot(!any(is.na(partable$ustart)))
  est <- out <- partable$ustart
  N <- length(est)

  # nblocks
  nblocks <- lav_partable_nblocks(partable)

  # if ov.var is NOT a list, make a list
  if (!is.list(ov.var)) {
    tmp <- ov.var
    ov.var <- vector("list", length = nblocks)
    ov.var[1:nblocks] <- list(tmp)
  }

  for (g in 1:nblocks) {
    ov.names <- vnames(partable, "ov", block = g) # not user
    lv.names <- vnames(partable, "lv", block = g)

    OV <- sqrt(ov.var[[g]])

    # 1a. "=~" regular indicators
    idx <- which(partable$op == "=~" & !(partable$rhs %in% lv.names) &
      partable$block == g)
    out[idx] <- out[idx] * OV[match(partable$rhs[idx], ov.names)]

    # 1b. "=~" regular higher-order lv indicators

    # 1c. "=~" indicators that are both in ov and lv
    # idx <- which(partable$op == "=~" & partable$rhs %in% ov.names
    #                             & partable$rhs %in% lv.names &
    #             partable$block == g)

    # 2. "~" regressions (and "<~")
    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$lhs %in% ov.names &
      partable$block == g)
    out[idx] <- out[idx] * OV[match(partable$lhs[idx], ov.names)]

    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$rhs %in% ov.names &
      partable$block == g)
    out[idx] <- out[idx] / OV[match(partable$rhs[idx], ov.names)]

    # 3a. "~~" ov
    # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances
    #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
    #            were i.var and j.var where diagonal elements of OV
    #
    #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
    #            elements are the 'THETA' diagonal elements!!

    # variances
    rv.idx <- which(partable$op == "~~" & !(partable$lhs %in% lv.names) &
      partable$lhs == partable$rhs &
      partable$block == g)
    out[rv.idx] <- (out[rv.idx] * OV[match(partable$lhs[rv.idx], ov.names)]
      * OV[match(partable$rhs[rv.idx], ov.names)])

    # covariances
    idx <- which(partable$op == "~~" & !(partable$lhs %in% lv.names) &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx) > 0L) {
      if (cov.std == FALSE) {
        out[idx] <- (out[idx] * OV[match(partable$lhs[idx], ov.names)]
          * OV[match(partable$rhs[idx], ov.names)])
      } else {
        if (!is.complex(out[rv.idx])) {
          RV <- sqrt(abs(out[rv.idx]))
        } else {
          RV <- sqrt(out[rv.idx])
        }
        rv.names <- partable$lhs[rv.idx]
        out[idx] <- (out[idx] * RV[match(partable$lhs[idx], rv.names)]
          * RV[match(partable$rhs[idx], rv.names)])
      }
    }

    # 3b. "~~" lv
    # idx <- which(partable$op == "~~" & partable$rhs %in% lv.names &
    #             partable$block == g)

    # 4a. "~1" ov
    idx <- which(partable$op == "~1" & !(partable$lhs %in% lv.names) &
      partable$block == g)
    out[idx] <- out[idx] * OV[match(partable$lhs[idx], ov.names)]

    # 4b. "~1" lv
    # idx <- which(partable$op == "~1" & partable$lhs %in% lv.names &
    #             partable$block == g)
  }

  # 5a ":="
  # 5b "=="
  # 5c. "<" or ">"

  out
}
