lav_lavaan_step16_rotation <- function(lavoptions = NULL,
                                       lavmodel = NULL,
                                       lavpartable = NULL,
                                       lavh1 = NULL,
                                       lavdata = NULL,
                                       x = NULL,
                                       lavvcov = NULL,
                                       VCOV = NULL, # nolint
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

  if ((.hasSlot(lavmodel, "nefa")) && (lavmodel@nefa > 0L) &&
    (lavoptions$rotation != "none")) {
    # store unrotated solution in partable
    lavpartable$est.unrotated <- lavpartable$est
    lavpartable$se.unrotated <- lavpartable$se

    # rotate, and create new lavmodel
    if (lavoptions$verbose) {
      cat(
        "rotating EFA factors using rotation method =",
        toupper(lavoptions$rotation), "..."
      )
    }
    x.unrotated <- as.numeric(x)
    lavmodel.unrot <- lavmodel
    efa.out <- lav_model_efa_rotate(
      lavmodel = lavmodel,
      lavoptions = lavoptions
    )

    # adapt partable:
    # - change 'free' column to reflect that user = 7/77 parameters are free
    # - save unrotated free column in free.unrotated
    lavpartable$free.unrotated <- lavpartable$free
    user7.idx <- which((lavpartable$user == 7L | lavpartable$user == 77L) &
      lavpartable$free == 0L)
    lavpartable$free[user7.idx] <- 1L
    lavpartable$free[lavpartable$free > 0L] <-
      seq_len(sum(lavpartable$free > 0L))

    # create 'rotated' lavmodel, reflecting the 'new' free parameters
    lavmodel <- lav_model(
      lavpartable = lavpartable,
      lavoptions = lavoptions,
      th.idx = lavmodel@th.idx
    )

    # add rotated information
    lavmodel@H <- efa.out$H
    lavmodel@lv.order <- efa.out$lv.order
    lavmodel@GLIST <- efa.out$GLIST

    # add con.jac information (if any)
    lavmodel@con.lambda <- lavmodel.unrot@con.lambda
    if (nrow(lavmodel.unrot@con.jac) > 0L) {
      con.jac <- rbind(lavmodel@ceq.JAC, lavmodel@cin.JAC)
      attr(con.jac, "inactive.idx") <-
        attr(lavmodel.unrot@con.jac, "inactive.idx")
      attr(con.jac, "cin.idx") <- attr(lavmodel.unrot@con.jac, "cin.idx")
      attr(con.jac, "ceq.idx") <- attr(lavmodel.unrot@con.jac, "ceq.idx")
      lavmodel@con.jac <- con.jac
    }

    # overwrite parameters in @ParTable$est
    lavpartable$est <- lav_model_get_parameters(
      lavmodel = lavmodel,
      type = "user", extra = TRUE
    )
    if (lavoptions$verbose) {
      cat(" done.\n")
    }

    # VCOV rotated parameters
    if (!lavoptions$se %in% c("none", "bootstrap", "external", "two.step")) {
      if (lavoptions$verbose) {
        cat(
          "computing VCOV for se =", lavoptions$se,
          "and rotation.se =", lavoptions$rotation.se, "..."
        )
      }

      # use delta rule to recompute vcov
      if (lavoptions$rotation.se == "delta") {
        # Jacobian
        JAC <- numDeriv::jacobian( # nolint
          func = lav_model_efa_rotate_x,
          x = x.unrotated, lavmodel = lavmodel.unrot,
          init.rot = lavmodel@H, lavoptions = lavoptions,
          type = "user", extra = FALSE,
          method.args = list(eps = 0.0050),
          method = "simple"
        ) # important!

        # force VCOV to be pd, before we transform (not very elegant)
        VCOV.in <- lav_matrix_symmetric_force_pd(lavvcov$vcov, # nolint
          tol = 1e-10
        )
        # apply Delta rule
        VCOV.user <- JAC %*% VCOV.in %*% t(JAC) # nolint

        # re-compute SE and store them in lavpartable
        tmp <- diag(VCOV.user)
        min.idx <- which(tmp < 0)
        if (length(min.idx) > 0L) {
          tmp[min.idx] <- as.numeric(NA)
        }
        tmp <- sqrt(tmp)
        # catch near-zero SEs  (was ^(1/2) < 0.6)
        zero.idx <- which(tmp < .Machine$double.eps^(1 / 3))
        if (length(zero.idx) > 0L) {
          tmp[zero.idx] <- 0.0
        }
        lavpartable$se <- tmp

        # store rotated VCOV
        # lavvcov$vcov.unrotated <- lavvcov$vcov
        if (.hasSlot(lavmodel, "ceq.simple.only") &&
          lavmodel@ceq.simple.only) {
          free.idx <- which(lavpartable$free > 0L &&
            !duplicated(lavpartable$free))
        } else {
          free.idx <- which(lavpartable$free > 0L)
        }
        lavvcov$vcov <- VCOV.user[free.idx, free.idx, drop = FALSE]

        # rotation.se = "bordered" is the default
      } else if (lavoptions$rotation.se == "bordered") {
        # create 'border' for augmented information matrix
        x.rot <- lav_model_get_parameters(lavmodel)
        JAC <- numDeriv::jacobian( # nolint
          func = lav_model_efa_rotate_border_x,
          x = x.rot, lavmodel = lavmodel,
          lavoptions = lavoptions,
          lavpartable = lavpartable,
          # method.args = list(eps = 0.0005),
          # method = "simple")
          method = "Richardson"
        )
        # store JAC
        lavmodel@ceq.efa.JAC <- JAC

        # no other constraints
        if (length(lavmodel@ceq.linear.idx) == 0L &&
          length(lavmodel@ceq.nonlinear.idx) == 0L &&
          length(lavmodel@cin.linear.idx) == 0L &&
          length(lavmodel@cin.nonlinear.idx) == 0L) {
          lavmodel@con.jac <- JAC
          attr(lavmodel@con.jac, "inactive.idx") <- integer(0L)
          attr(lavmodel@con.jac, "ceq.idx") <- seq_len(nrow(JAC))
          attr(lavmodel@con.jac, "cin.idx") <- integer(0L)
          lavmodel@con.lambda <- rep(0, nrow(JAC))

          # other constraints
        } else {
          inactive.idx <- attr(lavmodel@con.jac, "inactive.idx")
          ceq.idx <- attr(lavmodel@con.jac, "ceq.idx")
          cin.idx <- attr(lavmodel@con.jac, "cin.idx")
          lambda <- lavmodel@con.lambda
          nbord <- nrow(JAC)

          # reconstruct con.jac
          CON.JAC <- rbind(JAC, lavmodel@ceq.JAC, lavmodel@cin.JAC) # nolint
          attr(CON.JAC, "cin.idx") <- cin.idx + nbord # nolint
          attr(CON.JAC, "ceq.idx") <- c(1:nbord, ceq.idx + nbord) # nolint
          attr(CON.JAC, "inactive.idx") <- inactive.idx + nbord # nolint

          lavmodel@con.jac <- CON.JAC
          lavmodel@con.lambda <- c(rep(0, nbord), lambda)
        }

        # compute VCOV, taking 'rotation constraints' into account
        VCOV <- lav_model_vcov( # nolint
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
          lavpartable = lavpartable, VCOV = VCOV
        )
        lavpartable$se <- tmp

        # store rotated VCOV in lavvcov
        tmp.attr <- attributes(VCOV)
        VCOV1 <- VCOV # nolint
        attributes(VCOV1) <- tmp.attr["dim"] # nolint
        # lavvcov$vcov.unrotated <- lavvcov$vcov
        lavvcov$vcov <- VCOV1
      } # bordered

      if (lavoptions$verbose) {
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
