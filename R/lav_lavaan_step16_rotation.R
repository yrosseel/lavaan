lav_lavaan_step16_rotation <- function(lavoptions = NULL,
                                       lavmodel = NULL,
                                       lavpartable = NULL,
                                       lavh1 = NULL,
                                       lavdata = NULL,
                                       x = NULL,
                                       lavvcov = NULL,
                                       vcov = NULL,
                                       lavcache = NULL,
                                       lavimplied = NULL,
                                       lavsamplestats = NULL) {
  # # # # # # # # # # #
  # #  16. rotation # #
  # # # # # # # # # # #

  # if lavmodel@nefa > 0L and lavoptions$rotation not "none"
  #   store unrotated  solution in partable (column est.unrotated)
  #   rotate lavmodel via lav_model_efa_rotate and overwrite column est
  #     in partable
  #   if lavoptions$se not in none, bootstrap, external, twostep
  #     if lavoptions$rotation.se == "delta"
  #       re-compute vcov with delta rule (*)
  #       re-compute SE and store them in lavpartable (*)
  #     else if lavoptions$rotation.se == "bordered"
  #       create 'new' partable where the user = 7/77 parameters are free (*)
  #
  # (*) code too complicated to summarize here

  if (lavmodel@nefa > 0L &&
    (lavoptions$rotation != "none")) {
    # store unrotated solution in partable
    lavpartable$est.unrotated <- lavpartable$est
    lavpartable$se.unrotated <- lavpartable$se

    # rotate, and create new lavmodel
    if (lav_verbose()) {
      cat(
        "rotating EFA factors using rotation method =",
        toupper(lavoptions$rotation), "..."
      )
    }
    x_unrotated <- as.numeric(x)
    lavmodel_unrot <- lavmodel
    efa_out <- lav_model_efa_rotate(
      lavmodel = lavmodel,
      lavoptions = lavoptions
    )

    # adapt partable:
    # - change 'free' column to reflect that user = 7/77 parameters are free
    # - save unrotated free column in free.unrotated
    lavpartable$free.unrotated <- lavpartable$free
    user7_idx <- which((lavpartable$user == 7L | lavpartable$user == 77L) &
      lavpartable$free == 0L)
    lavpartable$free[user7_idx] <- 1L
    lavpartable$free[lavpartable$free > 0L] <-
      seq_len(sum(lavpartable$free > 0L))
    # avoid cin.simple entries for these user=7 parameters
    if (!is.null(lavpartable$lower)) {
      lavpartable$lower[user7_idx] <- -Inf
    }
    if (!is.null(lavpartable$upper)) {
      lavpartable$upper[user7_idx] <- +Inf
    }

    # create 'rotated' lavmodel, reflecting the 'new' free parameters
    lavmodel <- lav_model(
      lavpartable = lavpartable,
      lavoptions = lavoptions,
      th.idx = lavmodel@th.idx
    )

    # add rotated information
    lavmodel@H <- efa_out$H
    lavmodel@lv.order <- efa_out$lv.order
    lavmodel@GLIST <- efa_out$GLIST

    # add con.jac information (if any)
    lavmodel@con.lambda <- lavmodel_unrot@con.lambda
    if (nrow(lavmodel_unrot@con.jac) > 0L) {
      con_jac <- rbind(lavmodel@ceq.JAC, lavmodel@cin.JAC)
      attr(con_jac, "inactive.idx") <-
        attr(lavmodel_unrot@con.jac, "inactive.idx")
      attr(con_jac, "cin.idx") <- attr(lavmodel_unrot@con.jac, "cin.idx")
      # use nrow(lavmodel@ceq.JAC), not lavmodel.unrot (which may have had
      # zero rows removed by lav_constraints_parse for temporarily-fixed
      # user=7 EFA identification parameters)
      attr(con_jac, "ceq.idx") <- seq_len(nrow(lavmodel@ceq.JAC))
      lavmodel@con.jac <- con_jac
    }

    # overwrite parameters in @ParTable$est
    lavpartable$est <- lav_model_get_parameters(
      lavmodel = lavmodel,
      type = "user", extra = TRUE
    )
    if (lav_verbose()) {
      cat(" done.\n")
    }

    # VCOV rotated parameters
    if (!lavoptions$se %in% c("none", "bootstrap", "external", "two.step")) {
      if (lav_verbose()) {
        cat(
          "computing VCOV for se =", lavoptions$se,
          "and rotation.se =", lavoptions$rotation.se, "..."
        )
      }

      # use delta rule to recompute vcov
      if (lavoptions$rotation.se == "delta") {
        # Jacobian
        jac <- numDeriv::jacobian( # nolint
          func = lav_model_efa_rotate_x,
          x = x_unrotated, lavmodel = lavmodel_unrot,
          init.rot = lavmodel@H, lavoptions = lavoptions,
          type = "user", extra = FALSE,
          method.args = list(eps = 0.0050),
          method = "simple"
        ) # important!

        # force VCOV to be pd, before we transform (not very elegant)
        vcov_in <- lav_matrix_symmetric_force_pd(lavvcov$vcov, # nolint
          tol = 1e-10
        )
        # apply Delta rule
        vcov_user <- jac %*% vcov_in %*% t(jac) # nolint

        # re-compute SE and store them in lavpartable
        tmp <- diag(vcov_user)
        min_idx <- which(tmp < 0)
        if (length(min_idx) > 0L) {
          tmp[min_idx] <- as.numeric(NA)
        }
        tmp <- sqrt(tmp)
        # catch near-zero SEs  (was ^(1/2) < 0.6)
        zero_idx <- which(tmp < .Machine$double.eps^(1 / 3))
        if (length(zero_idx) > 0L) {
          tmp[zero_idx] <- 0.0
        }
        lavpartable$se <- tmp

        # store rotated VCOV
        # lavvcov$vcov.unrotated <- lavvcov$vcov
        if (lavmodel@ceq.simple.only) {
          free_idx <- which(lavpartable$free > 0L &&
            !duplicated(lavpartable$free))
        } else {
          free_idx <- which(lavpartable$free > 0L)
        }
        lavvcov$vcov <- vcov_user[free_idx, free_idx, drop = FALSE]

        # rotation.se = "bordered" is the default
      } else if (lavoptions$rotation.se == "bordered") {
        # create 'border' for augmented information matrix
        x_rot <- lav_model_get_parameters(lavmodel)
        jac <- numDeriv::jacobian( # nolint
          func = lav_model_efa_rotate_border_x,
          x = x_rot, lavmodel = lavmodel,
          lavoptions = lavoptions,
          lavpartable = lavpartable,
          # method.args = list(eps = 0.0005),
          # method = "simple")
          method = "Richardson"
        )
        # store JAC
        lavmodel@ceq.efa.JAC <- jac

        # no other constraints
        if (nrow(lavmodel@con.jac) == 0L) {
          lavmodel@con.jac <- jac
          attr(lavmodel@con.jac, "inactive.idx") <- integer(0L)
          attr(lavmodel@con.jac, "ceq.idx") <- seq_len(nrow(jac))
          attr(lavmodel@con.jac, "cin.idx") <- integer(0L)
          lavmodel@con.lambda <- rep(0, nrow(jac))

          # other constraints
        } else {
          inactive_idx <- attr(lavmodel@con.jac, "inactive.idx")
          ceq_idx <- attr(lavmodel@con.jac, "ceq.idx")
          cin_idx <- attr(lavmodel@con.jac, "cin.idx")
          lambda <- lavmodel@con.lambda
          nbord <- nrow(jac)

          # reconstruct con.jac
          con_jac_1 <- rbind(jac, lavmodel@ceq.JAC, lavmodel@cin.JAC)
          attr(con_jac_1, "cin.idx") <- cin_idx + nbord
          attr(con_jac_1, "ceq.idx") <- c(1:nbord, ceq_idx + nbord)
          attr(con_jac_1, "inactive.idx") <- inactive_idx + nbord

          lavmodel@con.jac <- con_jac_1
          lavmodel@con.lambda <- c(rep(0, nbord), lambda)
        }

        # compute VCOV, taking 'rotation constraints' into account
        vcov <- lav_model_vcov(
          lavmodel = lavmodel,
          lavsamplestats = lavsamplestats,
          lavoptions = lavoptions,
          lavdata = lavdata,
          lavpartable = lavpartable,
          lavcache = lavcache,
          lavimplied = lavimplied,
          lavh1 = lavh1
        )

        # compute SE and store them in lavpartable
        tmp <- lav_model_vcov_se(
          lavmodel = lavmodel,
          lavpartable = lavpartable, VCOV = vcov
        )
        lavpartable$se <- tmp

        # store rotated VCOV in lavvcov
        tmp_attr <- attributes(vcov)
        vcov1 <- vcov
        attributes(vcov1) <- tmp_attr["dim"]
        # lavvcov$vcov.unrotated <- lavvcov$vcov
        lavvcov$vcov <- vcov1
      } # bordered

      if (lav_verbose()) {
        cat(" done.\n")
      }
    } # vcov
  } # efa

  list(
    lavpartable = lavpartable,
    lavmodel    = lavmodel,
    lavvcov     = lavvcov
  )
}
