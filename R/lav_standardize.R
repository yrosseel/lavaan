lav_standardize_lv_x <- function(x, lavobject, partable = NULL, cov_std = TRUE,
                                 lv_var = NULL,
                                 rotation = FALSE) {
  # set new values for x
  lavmodel <- lav_model_set_parameters(lavmodel = lavobject@Model, x = x)

  if (rotation) {
    x_unrotated <- x
    lavmodel@GLIST <- lavTech(lavobject, "est.unrotated") # unrotated!
    est_rot <- lav_model_efa_rotate_x(
      x = x_unrotated,
      lavmodel = lavmodel, # unrotated!
      lavoptions = lavobject@Options,
      init_rot = lavmodel@H,
      type = "user",
      extra = TRUE
    )
    glist <- attr(est_rot, "extra")$glist
    attributes(est_rot) <- NULL
    est <- est_rot
  } else {
    glist <- lavmodel@GLIST
                    # if this changes, tag @TDJorgensen in commit message
    est <- lav_model_get_parameters(lavmodel, type = "user")
  }

  x_stand_user <- lav_standardize_lv(
    lavobject = lavobject,
    partable = partable, est = est,
    glist = glist, cov_std = cov_std,
    lv_var = lv_var
  )

  x_stand_user
}

lav_standardize_all_x <- function(x, lavobject, partable = NULL, cov_std = TRUE,
                                  rotation = FALSE) {
  lavmodel <- lav_model_set_parameters(lavmodel = lavobject@Model, x = x)

  if (rotation) {
    x_unrotated <- x
    lavmodel@GLIST <- lavTech(lavobject, "est.unrotated") # unrotated!
    est_rot <- lav_model_efa_rotate_x(
      x = x_unrotated,
      lavmodel = lavmodel, # unrotated!
      lavoptions = lavobject@Options,
      init_rot = lavmodel@H,
      type = "user",
      extra = TRUE
    )
    glist <- attr(est_rot, "extra")$glist
    attributes(est_rot) <- NULL
    est <- est_rot
  } else {
    glist <- lavmodel@GLIST
           # if this changes, tag @TDJorgensen in commit message
    est <- lav_model_get_parameters(lavmodel, type = "user")
  }

  x_stand_user <- lav_standardize_all(
    lavobject = lavobject,
    partable = partable, est = est,
    est_std = NULL, glist = glist,
    cov_std = cov_std
  )
  x_stand_user
}

lav_standardize_all_nox_x <- function(x, lavobject, partable = NULL,
                                      cov_std = TRUE, rotation = FALSE) {
  lavmodel <- lav_model_set_parameters(lavmodel = lavobject@Model, x = x)

  if (rotation) {
    x_unrotated <- x
    lavmodel@GLIST <- lavTech(lavobject, "est.unrotated") # unrotated!
    est_rot <- lav_model_efa_rotate_x(
      x = x_unrotated,
      lavmodel = lavmodel, # unrotated!
      lavoptions = lavobject@Options,
      init_rot = lavmodel@H,
      type = "user",
      extra = TRUE
    )
    glist <- attr(est_rot, "extra")$glist
    attributes(est_rot) <- NULL
    est <- est_rot
  } else {
    glist <- lavmodel@GLIST
                # if this changes, tag @TDJorgensen in commit message
    est <- lav_model_get_parameters(lavmodel, type = "user")
  }

  x_stand_user <- lav_standardize_all_nox(
    lavobject = lavobject,
    partable = partable, est = est,
    est_std = NULL, glist = glist,
    cov_std = cov_std
  )
  x_stand_user
}


lav_standardize_lv <- function(lavobject = NULL,
                               partable = NULL, est = NULL, glist = NULL,
                               cov_std = TRUE, lv_var = NULL,
                               lavmodel = NULL, lavpartable = NULL) {
  if (is.null(lavobject)) {
    stopifnot(!is.null(lavmodel))
    stopifnot(!is.null(lavpartable))
    if (is.null(est)) {
      if (!is.null(lavpartable$est)) {
        est <- lavpartable$est
                # if this changes, tag @TDJorgensen in commit message
      } else {
        lav_msg_stop(gettext("could not find `est' in lavpartable"))
      }
    }
  } else {
    lavmodel <- lavobject@Model
    lavpartable <- lavobject@ParTable
    if (is.null(est)) {
      est <- lav_inspect_est(lavobject)
    }
  }

  if (is.null(partable)) {
    partable <- lavpartable
  }
  if (is.null(glist)) {
    glist <- lavmodel@GLIST
  }

  out <- est
  n <- length(est)
  stopifnot(n == length(partable$lhs))

  nmat <- lavmodel@nmat

  # compute ETA
  if (is.null(lv_var)) {
    lv_eta <- lav_model_veta(
      lavmodel = lavmodel,
      glist = glist
    )

    lv_eeta <- lav_model_eeta(
      lavmodel = lavmodel,
      glist = glist,
      lavsamplestats = lavobject@SampleStats
    )
  }

  for (g in 1:lavmodel@nblocks) {
    ov_names <- lav_pt_vnames(lavpartable, "ov", block = g) # not user,
    # which may be incomplete
    lv_names <- lav_pt_vnames(lavpartable, "lv", block = g)

    # shortcut: no latents in this block, nothing to do
    if (length(lv_names) == 0L) {
      next
    }

    # which mm belong to block g?
    mm_in_group <- 1:nmat[g] + cumsum(c(0, nmat))[g]
    # mlist <- glist[mm_in_group]

    if (is.null(lv_var)) {
      eta2 <- diag(lv_eta[[g]])
      eeta <- lv_eeta[[g]]
    } else {
      eta2 <- lv_var[[g]]
      eeta <- numeric(length(eta2))
    }
    # change negative values to NA
    eta2[eta2 < 0] <- as.numeric(NA)
    eta <- sqrt(eta2)

    # Interaction/quadratic term correction (FV)
    # (based on Kelava & Brandt, 2022; Brandt et al., 2015)
    # For interaction terms A:B, the standardized coefficient should use
    # SD(A)*SD(B) instead of SD(A:B) (Eq. 13 in Brandt et al., 2015).
    lv_int_names <- lav_pt_vnames(lavpartable, "lv.interaction",
      block = g
    )

    if (length(lv_int_names) > 0L) {

      for (int_name in lv_int_names) {
        components <- strsplit(int_name, ":", fixed = TRUE)[[1L]]
        a <- components[1]
        b <- components[2]

        idx_int <- match(int_name, lv_names)
        idx_a   <- match(a, lv_names)
        idx_b   <- match(b, lv_names)

        # Need a, b and interaction term before moving on
        if (is.na(idx_a) || is.na(idx_b) || is.na(idx_int)) next

        exp_a <- eeta[idx_a]
        exp_b <- eeta[idx_b]

        # Do we need to shift simple main effects?
        if (exp_a != 0 || exp_b != 0) {
          # KS, 08/05/26
          #
          # This covers the general case of y ~ x + z + w + x:z + x:x
          # This should naturally extend to other conditions with more
          # interaction terms. e.g., y ~ x + z + w + x:z + x:x + w:z
          #
          # For y = b0 + b1*x + b2*z + b3*xz + b4x^2
          # b1 and b2 are the slope at the x=0 and z=0. When mean-centering
          # we have to re-define b1 and b2 to the expected slopes
          # (i.e., the partial derivatives) at x=mean(x), z=mean(z).
          #
          # Computing the partial derivatives we get
          #
          #   dy/dx = b1 + b3*z + 2*b4*x
          #   dy/dz = b2 + b3*x
          #
          # With x=mean(x) and z=mean(x), we get
          #
          #   dy/dx = b1 + b3 * mean(z) + 2 * b4 * mean(x)
          #   dy/dz = b2 + b3 * mean(x)
          #
          # Since we iteratively update out[<idx>] for b1 and b2, we carry
          # out all the necessary conditions, by the end of the loop.
          #
          # Iteration: 1, interaction term: x:z
          #   b1 <- b1 + b3 * mean(z)
          #   b2 <- b2 + b3 * mean(x)
          #
          # Iteration: 2, quadratic term: x:x
          #   b1 <- b1 + b4 * mean(x)
          #   b1 <- b1 + b4 * mean(x)
          #
          #   The above is thus equivalent to b1 + 2 * b4 * mean(x)
          #
          # Putting it all together we get
          #
          #   b1 <- b1 + b3 * mean(z) + 2 * b4 * mean(x)
          #   b2 <- b2 + b3 * mean(x)
          #

          # Find dependent variables
          lv_int_dep <- unique(partable$lhs[
            partable$op == "~" & partable$rhs == int_name
          ])

          for (dep in lv_int_dep) {
            idx_beta_a <- which(
              partable$lhs == dep & partable$op == "~" &
              partable$rhs == a & partable$group == g
            )

            idx_beta_b <- which(
              partable$lhs == dep & partable$op == "~" &
              partable$rhs == b & partable$group == g
            )

            beta_ab <- partable$est[
              partable$lhs == dep & partable$op == "~" &
              partable$rhs == int_name & partable$group == g
            ]

            out[idx_beta_a] <- out[idx_beta_a] + beta_ab * exp_b
            out[idx_beta_b] <- out[idx_beta_b] + beta_ab * exp_a
          }
        }

        # Replace ETA for interaction terms with SD(A)*SD(B)
        eta[idx_int] <- eta[idx_a] * eta[idx_b]
      }
    }

    # End interaction/quadratic term correction for now (08/05/26;KS)
    # The next step is to correctly scale the variances of the interaction terms

    # 1a. "=~" regular indicators
    idx <- which(partable$op == "=~" & !(partable$rhs %in% lv_names) &
      partable$block == g)
    out[idx] <- out[idx] * eta[match(partable$lhs[idx], lv_names)]

    # 1b. "=~" regular higher-order lv indicators
    idx <- which(partable$op == "=~" & !(partable$rhs %in% ov_names) &
      partable$block == g)
    out[idx] <- (out[idx] * eta[match(partable$lhs[idx], lv_names)]
      / eta[match(partable$rhs[idx], lv_names)])

    # 1c. "=~" indicators that are both in ov and lv
    # idx <- which(partable$op == "=~" & partable$rhs %in% ov.names
    #                             & partable$rhs %in% lv.names &
    #             partable$block == g)

    # 2. "~" regressions (and "<~")
    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$lhs %in% lv_names &
      partable$block == g)
    out[idx] <- out[idx] / eta[match(partable$lhs[idx], lv_names)]

    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$rhs %in% lv_names &
      partable$block == g)
    out[idx] <- out[idx] * eta[match(partable$rhs[idx], lv_names)]

    # 3a. "~~" ov
    # idx <- which(partable$op == "~~" & !(partable$lhs %in% lv.names) &
    #             partable$block == g)

    # 3b. "~~" lv
    # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances
    #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
    #            where i.var and j.var were diagonal elements of ETA
    #
    #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
    #            elements are the 'PSI' diagonal elements!!

    # variances
    rv_idx <- which(partable$op == "~~" & partable$rhs %in% lv_names &
      partable$lhs == partable$rhs &
      partable$block == g)
    out[rv_idx] <- (out[rv_idx] / eta[match(partable$lhs[rv_idx], lv_names)]
      / eta[match(partable$rhs[rv_idx], lv_names)])

    # covariances lv
    # three types:
    # - only lhs is LV (and fixed.x = FALSE)
    # - only rhs is LV (and fixed.x = FALSE)
    # - both lhs and rhs are LV (regular case)
    if (cov_std) {
      if (!is.complex(est[rv_idx])) {
        rv <- sqrt(abs(est[rv_idx])) # abs in case of heywood cases
      } else {
        rv <- sqrt(est[rv_idx])
      }
      rv_names <- partable$lhs[rv_idx]
    }

    # left
    idx_lhs <- which(partable$op == "~~" &
      partable$lhs %in% lv_names &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx_lhs) > 0L) {
      if (cov_std == FALSE) {
        out[idx_lhs] <-
          (out[idx_lhs] / eta[match(partable$lhs[idx_lhs], lv_names)])
      } else {
        out[idx_lhs] <-
          (out[idx_lhs] / rv[match(partable$lhs[idx_lhs], rv_names)])
      }
    }

    # right
    idx_rhs <- which(partable$op == "~~" &
      partable$rhs %in% lv_names &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx_rhs) > 0L) {
      if (cov_std == FALSE) {
        out[idx_rhs] <-
          (out[idx_rhs] / eta[match(partable$rhs[idx_rhs], lv_names)])
      } else {
        out[idx_rhs] <-
          (out[idx_rhs] / rv[match(partable$rhs[idx_rhs], rv_names)])
      }
    }


    # 4a. "~1" ov
    # idx <- which(partable$op == "~1" & !(partable$lhs %in% lv.names) &
    #             partable$block == g)

    # 4b. "~1" lv
    idx <- which(partable$op == "~1" & partable$lhs %in% lv_names &
      partable$block == g)
    out[idx] <- out[idx] / eta[match(partable$lhs[idx], lv_names)]
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
                                partable = NULL, est = NULL, est_std = NULL,
                                glist = NULL, cov_std = TRUE, ov_var = NULL,
                                lv_var = NULL,
                                lavmodel = NULL, lavpartable = NULL,
                                cov_x = NULL) {
  if (is.null(lavobject)) {
    stopifnot(!is.null(lavmodel))
    stopifnot(!is.null(lavpartable))
    if (is.null(est)) {
      if (!is.null(lavpartable$est)) {
        est <- lavpartable$est
             # if this changes, tag @TDJorgensen in commit message
      } else {
        lav_msg_stop(gettext("could not find `est' in lavpartable"))
      }
    }
  } else {
    lavmodel <- lavobject@Model
    lavpartable <- lavobject@ParTable
    if (is.null(est)) {
      est <- lav_inspect_est(lavobject)
    }
    if (lavmodel@conditional.x) {
      if (is.null(cov_x)) {
        # try SampleStats slot
        # if("SampleStats" %in% slotNames(lavobject)) {
        #    cov.x <- lavobject@SampleStats@cov.x
        if (!is.null(lavobject@implied$cov.x[[1]])) {
          cov_x <- lavobject@implied$cov.x
             # if this changes, tag @TDJorgensen in commit message
        } else {
          # perhaps lavaanList object
          # extract it from GLIST per block
          cov_x <- vector("list", length = lavmodel@nblocks)
          for (b in seq_len(lavmodel@nblocks)) {
            # which mm belong to block b?
            mm_in_block <- (seq_len(lavmodel@nmat[b]) +
              cumsum(c(0, lavmodel@nmat))[b])
            mlist <- lavmodel@GLIST[mm_in_block]
            cov_x[[b]] <- mlist[["cov.x"]]
          }
        }
      }
    }
  }

  if (is.null(partable)) {
    partable <- lavpartable
  }
  if (is.null(glist)) {
    glist <- lavmodel@GLIST
  }
  if (is.null(est_std)) {
    est_std <- lav_standardize_lv(
      lavobject = lavobject,
      partable = partable, est = est, glist = glist,
      cov_std = cov_std, lv_var = lv_var, lavmodel = lavmodel,
      lavpartable = lavpartable
    )
  }

  out <- est_std
  n <- length(est_std)
  stopifnot(n == length(partable$lhs))

  vy <- lav_model_vy(
    lavmodel = lavmodel, glist = glist,
    diagonal_only = TRUE
  )


  for (g in 1:lavmodel@nblocks) {
    ov_names <- lav_pt_vnames(lavpartable, "ov", block = g) # not user
    lv_names <- lav_pt_vnames(lavpartable, "lv", block = g)

    if (is.null(ov_var)) {
      ov2 <- vy[[g]]
      # replace zero values by NA (but keep negative values)
      zero_idx <- which(abs(ov2) < .Machine$double.eps)
      if (length(zero_idx) > 0L) {
        ov2[zero_idx] <- as.numeric(NA)
      }

      # replace negative values by NA (for sqrt)
      tmp_ov2 <- ov2
      neg_idx <- which(tmp_ov2 < 0)
      if (length(neg_idx) > 0L) {
        tmp_ov2[neg_idx] <- as.numeric(NA)
      }
      ov <- sqrt(tmp_ov2)
    } else {
      ov2 <- ov_var[[g]]
      ov <- sqrt(ov2)
    }

    if (lavmodel@conditional.x) {
      # extend OV with ov.names.x
      ov_names_x <- lav_pt_vnames(lavpartable, "ov.x", block = g)
      ov_names_nox <- lav_pt_vnames(lavpartable, "ov.nox", block = g)
      ov_names <- c(ov_names_nox, ov_names_x)
      ov2 <- c(ov2, diag(cov_x[[g]]))
      ov <- c(ov, sqrt(diag(cov_x[[g]])))
    }

    # 1a. "=~" regular indicators
    idx <- which(partable$op == "=~" & !(partable$rhs %in% lv_names) &
      partable$block == g)
    out[idx] <- out[idx] / ov[match(partable$rhs[idx], ov_names)]

    # 1b. "=~" regular higher-order lv indicators

    # 1c. "=~" indicators that are both in ov and lv
    # idx <- which(partable$op == "=~" & partable$rhs %in% ov.names
    #                             & partable$rhs %in% lv.names &
    #             partable$block == g)

    # 2. "~" regressions (and "<~")
    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$lhs %in% ov_names &
      partable$block == g)
    out[idx] <- out[idx] / ov[match(partable$lhs[idx], ov_names)]

    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$rhs %in% ov_names &
      partable$block == g)
    out[idx] <- out[idx] * ov[match(partable$rhs[idx], ov_names)]

    # 3a. "~~" ov
    # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances
    #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
    #            where i.var and j.var were diagonal elements of OV
    #
    #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
    #            elements are the 'THETA' diagonal elements!!

    # variances
    rv_idx <- which(partable$op == "~~" & !(partable$lhs %in% lv_names) &
      partable$lhs == partable$rhs &
      partable$block == g)
    # out[rv.idx] <- ( out[rv.idx] / OV[ match(partable$lhs[rv.idx], ov.names) ]
    #                          / OV[ match(partable$rhs[rv.idx], ov.names) ] )
    out[rv_idx] <- (out[rv_idx] /
      ov2[match(partable$lhs[rv_idx], ov_names)])

    # covariances ov
    # three types:
    # - only lhs is OV (and fixed.x = FALSE)
    # - only rhs is OV (and fixed.x = FALSE)
    # - both lhs and rhs are OV (regular case)
    if (cov_std) {
      if (!is.complex(est[rv_idx])) {
        rv <- sqrt(abs(est[rv_idx]))
      } else {
        rv <- sqrt(est[rv_idx])
      }
      rv_names <- partable$lhs[rv_idx]
    }

    # left
    idx_lhs <- which(partable$op == "~~" &
      !(partable$lhs %in% lv_names) &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx_lhs) > 0L) {
      if (cov_std == FALSE) {
        out[idx_lhs] <-
          (out[idx_lhs] / ov[match(partable$lhs[idx_lhs], ov_names)])
      } else {
        out[idx_lhs] <-
          (out[idx_lhs] / rv[match(partable$lhs[idx_lhs], rv_names)])
      }
    }

    # right
    idx_rhs <- which(partable$op == "~~" &
      !(partable$rhs %in% lv_names) &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx_rhs) > 0L) {
      if (cov_std == FALSE) {
        out[idx_rhs] <-
          (out[idx_rhs] / ov[match(partable$rhs[idx_rhs], ov_names)])
      } else {
        out[idx_rhs] <-
          (out[idx_rhs] / rv[match(partable$rhs[idx_rhs], rv_names)])
      }
    }

    # 3b. "~~" lv
    # idx <- which(partable$op == "~~" & partable$rhs %in% lv.names &
    #             partable$block == g)

    # 4a. "~1" ov
    idx <- which(partable$op == "~1" & !(partable$lhs %in% lv_names) &
      partable$block == g)
    out[idx] <- out[idx] / ov[match(partable$lhs[idx], ov_names)]

    # 4b. "~1" lv
    # idx <- which(partable$op == "~1" & partable$lhs %in% lv.names &
    #             partable$block == g)

    # 4c. "|" thresholds
    idx <- which(partable$op == "|" & !(partable$lhs %in% lv_names) &
      partable$block == g)
    out[idx] <- out[idx] / ov[match(partable$lhs[idx], ov_names)]

    # 4d. "~*~" scales
    idx <- which(partable$op == "~*~" & !(partable$lhs %in% lv_names) &
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
                                    partable = NULL, est = NULL, est_std = NULL,
                                    glist = NULL, cov_std = TRUE, ov_var = NULL,
                                    lv_var = NULL,
                                    lavmodel = NULL, lavpartable = NULL,
                                    cov_x = NULL) {
  if (is.null(lavobject)) {
    stopifnot(!is.null(lavmodel))
    stopifnot(!is.null(lavpartable))
    if (is.null(est)) {
      if (!is.null(lavpartable$est)) {
        est <- lavpartable$est
                   # if this changes, tag @TDJorgensen in commit message
      } else {
        lav_msg_stop(gettext("could not find `est' in lavpartable"))
      }
    }
  } else {
    lavmodel <- lavobject@Model
    lavpartable <- lavobject@ParTable
    if (is.null(est)) {
      est <- lav_inspect_est(lavobject)
    }
    if (lavmodel@conditional.x) {
      if (is.null(cov_x)) {
        # try SampleStats slot
        # if("SampleStats" %in% slotNames(lavobject)) {
        #    cov.x <- lavobject@SampleStats@cov.x
        if (!is.null(lavobject@implied$cov.x[[1]])) {
          cov_x <- lavobject@implied$cov.x
                      # if this changes, tag @TDJorgensen in commit message
        } else {
          # perhaps lavaanList object
          # extract it from GLIST per block
          cov_x <- vector("list", length = lavmodel@nblocks)
          for (b in seq_len(lavmodel@nblocks)) {
            # which mm belong to block b?
            mm_in_block <- (seq_len(lavmodel@nmat[b]) +
              cumsum(c(0, lavmodel@nmat))[b])
            mlist <- lavmodel@GLIST[mm_in_block]
            cov_x[[b]] <- mlist[["cov.x"]]
          }
        }
      }
    }
  }

  if (is.null(partable)) {
    partable <- lavpartable
  }
  if (is.null(glist)) {
    glist <- lavmodel@GLIST
  }
  if (is.null(est_std)) {
    est_std <- lav_standardize_lv(
      lavobject = lavobject,
      partable = partable, est = est, glist = glist,
      cov_std = cov_std, lv_var = lv_var, lavmodel = lavmodel,
      lavpartable = lavpartable
    )
  }

  out <- est_std
  n <- length(est_std)
  stopifnot(n == length(partable$lhs))

  vy <- lav_model_vy(
    lavmodel = lavmodel, glist = glist,
    diagonal_only = TRUE
  )


  for (g in 1:lavmodel@nblocks) {
    ov_names <- lav_pt_vnames(lavpartable, "ov", block = g)
    ov_names_x <- lav_pt_vnames(lavpartable, "ov.x", block = g)
    ov_names_nox <- lav_pt_vnames(lavpartable, "ov.nox", block = g)
    lv_names <- lav_pt_vnames(lavpartable, "lv", block = g)

    if (is.null(ov_var)) {
      ov2 <- vy[[g]]
      # replace zero values by NA (but keep negative values)
      zero_idx <- which(abs(ov2) < .Machine$double.eps)
      if (length(zero_idx) > 0L) {
        ov2[zero_idx] <- as.numeric(NA)
      }

      # replace negative values by NA (for sqrt)
      tmp_ov2 <- ov2
      neg_idx <- which(tmp_ov2 < 0)
      if (length(neg_idx) > 0L) {
        tmp_ov2[neg_idx] <- as.numeric(NA)
      }
      ov <- sqrt(tmp_ov2)
    } else {
      ov2 <- ov_var[[g]]
      ov <- sqrt(ov2)
    }


    if (lavmodel@conditional.x) {
      # extend OV with ov.names.x
      ov_names_x <- lav_pt_vnames(lavpartable, "ov.x", block = g)
      ov_names <- c(ov_names_nox, ov_names_x)
      ov2 <- c(ov2, diag(cov_x[[g]]))
      ov <- c(ov, sqrt(diag(cov_x[[g]])))
    }

    # 1a. "=~" regular indicators
    idx <- which(partable$op == "=~" & !(partable$rhs %in% lv_names) &
      partable$block == g)
    out[idx] <- out[idx] / ov[match(partable$rhs[idx], ov_names)]

    # 1b. "=~" regular higher-order lv indicators

    # 1c. "=~" indicators that are both in ov and lv
    # idx <- which(partable$op == "=~" & partable$rhs %in% ov.names
    #                             & partable$rhs %in% lv.names &
    #             partable$block == g)

    # 2. "~" regressions (and "<~")
    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$lhs %in% ov_names &
      partable$block == g)
    out[idx] <- out[idx] / ov[match(partable$lhs[idx], ov_names)]

    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$rhs %in% ov_names_nox &
      partable$block == g)
    out[idx] <- out[idx] * ov[match(partable$rhs[idx], ov_names_nox)]

    # 3a. "~~" ov
    # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances
    #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
    #            where i.var and j.var were diagonal elements of OV
    #
    #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
    #            elements are the 'THETA' diagonal elements!!

    # variances
    rv_idx <- which(partable$op == "~~" & !(partable$lhs %in% lv_names) &
      !(partable$lhs %in% ov_names_x) &
      partable$lhs == partable$rhs &
      partable$block == g)
    # out[rv.idx] <- ( out[rv.idx] / OV[ match(partable$lhs[rv.idx], ov.names) ]
    #                            / OV[ match(partable$rhs[rv.idx], ov.names) ] )
    out[rv_idx] <- (out[rv_idx] /
      ov2[match(partable$lhs[rv_idx], ov_names)])

    # covariances ov
    # three types:
    # - only lhs is OV (and fixed.x = FALSE)
    # - only rhs is OV (and fixed.x = FALSE)
    # - both lhs and rhs are OV (regular case)
    if (cov_std) {
      if (!is.complex(est[rv_idx])) {
        rv <- sqrt(abs(est[rv_idx]))
      } else {
        rv <- sqrt(est[rv_idx])
      }
      rv_names <- partable$lhs[rv_idx]
    }

    # left
    idx_lhs <- which(partable$op == "~~" &
      !(partable$lhs %in% lv_names) &
      !(partable$lhs %in% ov_names_x) &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx_lhs) > 0L) {
      if (cov_std == FALSE) {
        out[idx_lhs] <-
          (out[idx_lhs] / ov[match(partable$lhs[idx_lhs], ov_names)])
      } else {
        out[idx_lhs] <-
          (out[idx_lhs] / rv[match(partable$lhs[idx_lhs], rv_names)])
      }
    }

    # right
    idx_rhs <- which(partable$op == "~~" &
      !(partable$rhs %in% lv_names) &
      !(partable$rhs %in% ov_names_x) &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx_rhs) > 0L) {
      if (cov_std == FALSE) {
        out[idx_rhs] <-
          (out[idx_rhs] / ov[match(partable$rhs[idx_rhs], ov_names)])
      } else {
        out[idx_rhs] <-
          (out[idx_rhs] / rv[match(partable$rhs[idx_rhs], rv_names)])
      }
    }

    # 3b. "~~" lv
    # idx <- which(partable$op == "~~" & partable$rhs %in% lv.names &
    #             partable$block == g)

    # 4a. "~1" ov
    idx <- which(partable$op == "~1" & !(partable$lhs %in% lv_names) &
      !(partable$lhs %in% ov_names_x) &
      partable$block == g)
    out[idx] <- out[idx] / ov[match(partable$lhs[idx], ov_names)]

    # 4b. "~1" lv
    # idx <- which(partable$op == "~1" & partable$lhs %in% lv.names &
    #             partable$block == g)

    # 4c. "|" thresholds
    idx <- which(partable$op == "|" & !(partable$lhs %in% lv_names) &
      partable$block == g)
    out[idx] <- out[idx] / ov[match(partable$lhs[idx], ov_names)]

    # 4d. "~*~" scales
    idx <- which(partable$op == "~*~" & !(partable$lhs %in% lv_names) &
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

lav_unstandardize_ov <- function(partable, ov_var = NULL, cov_std = TRUE) {
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
  # n <- length(est)

  # nblocks
  nblocks <- lav_pt_nblocks(partable)

  # if ov.var is NOT a list, make a list
  if (!is.list(ov_var)) {
    tmp <- ov_var
    ov_var <- vector("list", length = nblocks)
    ov_var[1:nblocks] <- list(tmp)
  }

  for (g in 1:nblocks) {
    ov_names <- lav_pt_vnames(partable, "ov", block = g) # not user
    lv_names <- lav_pt_vnames(partable, "lv", block = g)

    ov <- sqrt(ov_var[[g]])

    # 1a. "=~" regular indicators
    idx <- which(partable$op == "=~" & !(partable$rhs %in% lv_names) &
      partable$block == g)
    out[idx] <- out[idx] * ov[match(partable$rhs[idx], ov_names)]

    # 1b. "=~" regular higher-order lv indicators

    # 1c. "=~" indicators that are both in ov and lv
    # idx <- which(partable$op == "=~" & partable$rhs %in% ov.names
    #                             & partable$rhs %in% lv.names &
    #             partable$block == g)

    # 2. "~" regressions (and "<~")
    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$lhs %in% ov_names &
      partable$block == g)
    out[idx] <- out[idx] * ov[match(partable$lhs[idx], ov_names)]

    idx <- which((partable$op == "~" | partable$op == "<~") &
      partable$rhs %in% ov_names &
      partable$block == g)
    out[idx] <- out[idx] / ov[match(partable$rhs[idx], ov_names)]

    # 3a. "~~" ov
    # ATTENTION: in Mplus 4.1, the off-diagonal residual covariances
    #            were computed by the formula cov(i,j) / sqrt(i.var*j.var)
    #            where i.var and j.var were diagonal elements of OV
    #
    #            in Mplus 6.1 (but also AMOS and EQS), the i.var and j.var
    #            elements are the 'THETA' diagonal elements!!

    # variances
    rv_idx <- which(partable$op == "~~" & !(partable$lhs %in% lv_names) &
      partable$lhs == partable$rhs &
      partable$block == g)
    out[rv_idx] <- (out[rv_idx] * ov[match(partable$lhs[rv_idx], ov_names)]
      * ov[match(partable$rhs[rv_idx], ov_names)])

    # covariances
    idx <- which(partable$op == "~~" & !(partable$lhs %in% lv_names) &
      partable$lhs != partable$rhs &
      partable$block == g)
    if (length(idx) > 0L) {
      if (cov_std == FALSE) {
        out[idx] <- (out[idx] * ov[match(partable$lhs[idx], ov_names)]
          * ov[match(partable$rhs[idx], ov_names)])
      } else {
        if (!is.complex(out[rv_idx])) {
          rv <- sqrt(abs(out[rv_idx]))
        } else {
          rv <- sqrt(out[rv_idx])
        }
        rv_names <- partable$lhs[rv_idx]
        out[idx] <- (out[idx] * rv[match(partable$lhs[idx], rv_names)]
          * rv[match(partable$rhs[idx], rv_names)])
      }
    }

    # 3b. "~~" lv
    # idx <- which(partable$op == "~~" & partable$rhs %in% lv.names &
    #             partable$block == g)

    # 4a. "~1" ov
    idx <- which(partable$op == "~1" & !(partable$lhs %in% lv_names) &
      partable$block == g)
    out[idx] <- out[idx] * ov[match(partable$lhs[idx], ov_names)]

    # 4b. "~1" lv
    # idx <- which(partable$op == "~1" & partable$lhs %in% lv.names &
    #             partable$block == g)
  }

  # 5a ":="
  # 5b "=="
  # 5c. "<" or ">"

  out
}
