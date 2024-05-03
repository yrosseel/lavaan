# compute two-step standard errors for SAM models

lav_sam_step2_se <- function(FIT = NULL, JOINT = NULL,
                             STEP1 = NULL, STEP2 = NULL,
                             local.options = list()) {
  # current approach for se = "twostep":
  # - create 'global' model, only to get the 'joint' information matrix
  # - partition information matrix (step 1, step 2)
  # - apply two-step correction for second step
  # - 'insert' these corrected SEs (and vcov) in JOINT

  out <- list()
  Sigma.11 <- STEP1$Sigma.11
  step1.free.idx <- STEP1$step1.free.idx
  step2.free.idx <- STEP2$step2.free.idx
  lavoptions <- FIT@Options
  nlevels <- FIT@pta$nlevels
  FIT.PA <- STEP2$FIT.PA
  extra.id <- STEP2$extra.id

  # catch empty step2.free.idx
  if (length(step2.free.idx) == 0L) {
    # no (free) structural parameters at all!
    out <- list(
      V1 = matrix(0, 0, 0), V2 = matrix(0, 0, 0),
      VCOV = matrix(0, 0, 0)
    )
    return(out)
  }

  if (!lavoptions$se %in%
    c("none", "standard", "naive", "twostep", "twostep2")) {
    lav_msg_warn(gettext(
      "unknown se= argument: \"%s\". Switching to twostep.",
      lavoptions$se
    ))
  }

  if (lavoptions$se == "none") {
    return(out)
  }

  if (lavoptions$verbose) {
    cat("Computing ", lavoptions$se, " standard errors ... ", sep = "")
  }

  INFO <- lavInspect(JOINT, "information")
  I.12 <- INFO[step1.free.idx, step2.free.idx]
  I.22 <- INFO[step2.free.idx, step2.free.idx]
  I.21 <- INFO[step2.free.idx, step1.free.idx]

  # V2
  if (nlevels > 1L) {
    # FIXME: not ok for multigroup multilevel
    N <- FIT@Data@Lp[[1]]$nclusters[[2]] # first group only
  } else {
    N <- nobs(FIT)
  }

  # do we have 'extra' free parameter in FIT.PA that are not free in JOINT?
  step2.rm.idx <- integer(0L)
  if (length(extra.id) > 0L) {
    id.idx <- which(FIT.PA@ParTable$id %in% extra.id &
      FIT.PA@ParTable$free > 0L)
    step2.rm.idx <- FIT.PA@ParTable$free[id.idx]
  }

  # invert augmented information, for I.22 block only
  # new in 0.6-16 (otherwise, eq constraints in struc part are ignored)
  if (lavoptions$se != "naive") {
    I.22.inv <-
      lav_model_information_augment_invert(
        lavmodel = FIT.PA@Model,
        information = I.22,
        inverted = TRUE,
        rm.idx = step2.rm.idx
      )
    if (inherits(I.22.inv, "try-error")) {
      # hm, not good
      if (lavoptions$se != "naive") {
        lav_msg_warn(gettext(
          "problem inverting information matrix (I.22); -> switching
            to naive standard errors!"
        ))
        lavoptions$se <- "naive"
      }
    }
  } # se is not "naive", but based  on I.22

  # method below has the advantage that we can use a 'robust' vcov
  # for the joint model;
  # but does not work if we have equality constraints in the MM!
  # -> D will be singular
  # A <- JOINT@vcov$vcov[ step2.free.idx,  step2.free.idx]
  # B <- JOINT@vcov$vcov[ step2.free.idx, -step2.free.idx]
  # C <- JOINT@vcov$vcov[-step2.free.idx,  step2.free.idx]
  # D <- JOINT@vcov$vcov[-step2.free.idx, -step2.free.idx]
  # I.22.inv <- A - B %*% solve(D) %*% C

  # se = "standard"
  if (lavoptions$se == "standard") {
    VCOV <- 1 / N * I.22.inv
    out$VCOV <- VCOV

  # se = "naive"
  } else if (lavoptions$se == "naive") {
    if (is.null(FIT.PA@vcov$vcov)) {
      FIT.PA@Options$se <- "standard"
      VCOV.naive <- lavTech(FIT.PA, "vcov")
    } else {
      VCOV.naive <- FIT.PA@vcov$vcov
    }
    if (length(step2.rm.idx) > 0L) {
      VCOV.naive <- VCOV.naive[-step2.rm.idx, -step2.rm.idx]
    }
    out$VCOV <- VCOV.naive

  # se = "twostep" or "twostep2"
  } else if (lavoptions$se %in% c("twostep", "twostep2")) {
    V2 <- 1 / N * I.22.inv # not the same as FIT.PA@vcov$vcov!!

	if (lavoptions$se == "twostep" ) {
      V1 <- I.22.inv %*% I.21 %*% Sigma.11 %*% I.12 %*% I.22.inv
	} else {
	  stop("not ready yet")
	  # R.21 <- crossprod(J2, J1)/nrow(J2)
	  # V1 <- I.22.inv %*% (I.21 %*% Sigma.11 %*% t(I.21) -
	  #                     I.21 %*% Sigma.11 %*% t(R.21) -
      #                     R.21 %*% Sigma.11 %*% t(I.21)) %*% I.22.inv
	}

    # V for second step
    if (!is.null(local.options$alpha.correction) &&
      local.options$alpha.correction > 0) {
      alpha.N1 <- local.options$alpha.correction / (N - 1)
      if (alpha.N1 > 1.0) {
        alpha.N1 <- 1.0
      } else if (alpha.N1 < 0.0) {
        alpha.N1 <- 0.0
      }
      if (is.null(FIT.PA@vcov$vcov)) {
        FIT.PA@Options$se <- "standard"
        VCOV.naive <- lavTech(FIT.PA, "vcov")
      } else {
        VCOV.naive <- FIT.PA@vcov$vcov
      }
      if (length(step2.rm.idx) > 0L) {
        VCOV.naive <- VCOV.naive[-step2.rm.idx, -step2.rm.idx]
      }
      VCOV.corrected <- V2 + V1
      VCOV <- alpha.N1 * VCOV.naive + (1 - alpha.N1) * VCOV.corrected
    } else {
      VCOV <- V2 + V1
    }

    # store in out
    out$V2 <- V2
    out$V1 <- V1
    out$VCOV <- VCOV
  }

  if (lavoptions$verbose) {
    cat("done.\n")
  }

  out
}
