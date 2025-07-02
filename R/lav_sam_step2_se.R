# compute two-step standard errors for SAM models
#
# several possibilities:
# 1) se = "twostep": classic (but global) two-step corrected SEs
#     - create 'global' model, only to get the 'joint' information matrix
#     - partition information matrix (step 1, step 2)
#     - apply two-step correction for second step
#     - 'insert' these corrected SEs (and vcov) in JOINT
# 2) se = "standard": using I.22.inv (but without correction term)
# 3) se = "naive": grab (naive) VCOV from FIT.PA
# 4) se = "local": grab (robust) VCOV from FIT.PA

lav_sam_step2_se <- function(FIT = NULL, JOINT = NULL,
                             STEP1 = NULL, STEP2 = NULL,
                             local.options = list()) {
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
    c("none", "standard", "naive", "twostep", "twostep.robust",
      "local", "local.nt")) {
    lav_msg_warn(gettextf(
      "unknown se= argument: %s. Switching to twostep.",
      lavoptions$se
    ))
  }

  if (lavoptions$se == "none") {
    return(out)
  }

  if (lav_verbose()) {
    cat("Computing ", lavoptions$se, " standard errors ... ", sep = "")
  }

  if (lavoptions$se %in% c("naive", "twostep", "twostep.robust")) {
    INFO <- lavInspect(JOINT, "information")
    I.12 <- INFO[step1.free.idx, step2.free.idx]
    I.22 <- INFO[step2.free.idx, step2.free.idx]
    I.21 <- INFO[step2.free.idx, step1.free.idx]
  }

  # V2
  if (nlevels > 1L) {
    # FIXME: not ok for multigroup multilevel
    N <- FIT@Data@Lp[[1]]$nclusters[[2]] # first group only
  } else {
    N <- nobs(FIT)
  }

  # total number of free parameters STRUC
  if (FIT.PA@Model@ceq.simple.only) {
    npar <- FIT.PA@Model@nx.unco
    PTS.free <- FIT.PA@ParTable$free
    PTS.free[PTS.free > 0] <- seq_len(npar)
  } else {
    npar <- FIT.PA@Model@nx.free
    PTS.free <- FIT.PA@ParTable$free
  }

  # do we have 'extra' free parameter in FIT.PA that are not free in JOINT?
  step2.rm.idx <- integer(0L)
  if (length(extra.id) > 0L) {
    id.idx <- which(FIT.PA@ParTable$id %in% extra.id &
      FIT.PA@ParTable$free > 0L)
    step2.rm.idx <- PTS.free[id.idx]
  }

  # invert augmented information, for I.22 block only
  # new in 0.6-16 (otherwise, eq constraints in struc part are ignored)
  if (lavoptions$se %in% c("standard", "twostep", "twostep.robust")) {
    I.22.inv <-
      lav_model_information_augment_invert(
        lavmodel = FIT.PA@Model,
        information = I.22,
        inverted = TRUE,
        use.ginv = FALSE, # if interaction, SEs end up smaller than naive...
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
  } # se needs I.22.inv

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

  # se = "naive" or "local": grab VCOV directly from FIT.PA
  } else if (lavoptions$se %in% c("naive", "local", "local.nt")) {
    if (is.null(FIT.PA@vcov$vcov)) {
      FIT.PA@Options$se <- "standard"
      VCOV <- lavTech(FIT.PA, "vcov")
    } else {
      VCOV <- FIT.PA@vcov$vcov
    }
    if (length(step2.rm.idx) > 0L) {
      VCOV <- VCOV[-step2.rm.idx, -step2.rm.idx]
    }
    # order rows/cols of VCOV, so that they correspond with the (step 2)
    # parameters of the JOINT model
    idx <- sort.int(STEP2$pt.idx, index.return = TRUE)$ix
    VCOV <- VCOV[idx, idx]

    out$VCOV <- VCOV

  # se = "twostep" or "twostep.robust"
  } else if (lavoptions$se  == "twostep" || lavoptions$se == "twostep.robust") {

    if (lavoptions$se  == "twostep") {
      V2 <- 1 / N * I.22.inv # not the same as FIT.PA@vcov$vcov!!
      V1 <- I.22.inv %*% I.21 %*% Sigma.11 %*% I.12 %*% I.22.inv
    } else if(lavoptions$se == "twostep.robust") {
      # following Yuan & Chan 2002, eqs 4, 10, 11, 12, 13 and 14
      # but for V11, V12, V21, V22: we use index '1' for step1, and '2'
      # for step 2!!

      A <- -1 * INFO[step2.free.idx, step2.free.idx, drop = FALSE]
      B <- -1 * INFO[step2.free.idx, step1.free.idx, drop = FALSE]

      # get P (for a single group!! for now)
      P <- lav_sam_step1_local_jac(STEP1 = STEP1, FIT = FIT, P.only = TRUE)

      # get V22
      if (is.null(JOINT@SampleStats@NACOV[[1]])) {
        JOINT@SampleStats@NACOV <- lavTech(JOINT, "gamma")
      }
      tmp <- lav_model_nvcov_robust_sem(
        lavmodel = JOINT@Model, lavsamplestats = JOINT@SampleStats,
        lavcache = JOINT@cache, lavdata = JOINT@Data,
        lavimplied = JOINT@implied, lavh1 = JOINT@h1,
        lavoptions = JOINT@Options, use.ginv = FALSE,
        attr.Delta = TRUE, attr.tDVGVD = TRUE, attr.E.inv = TRUE,
        attr.WLS.V = TRUE)
      NVarCov <- tmp[,] # remove attributes
      Delta <- attr(tmp, "Delta")
      E.inv <- attr(tmp, "E.inv")
      WLS.V <- attr(tmp, "WLS.V")
      tDVGVD <- attr(tmp, "tDVGVD")

      V22 <- tDVGVD[ step2.free.idx, step2.free.idx, drop = FALSE] # ok

      # FIXME: for a single group only:
      V11 <- P %*% lavTech(JOINT, "gamma")[[1]] %*% t(P)
      V21 <- t(Delta[[1]][, step2.free.idx]) %*% WLS.V[[1]] %*% lavTech(JOINT, "gamma")[[1]] %*% t(P)
      V12 <- t(V21)

      #V11 <- NVarCov[step1.free.idx, step1.free.idx, drop = FALSE]
      #V21 <- tmp2[   step2.free.idx, step1.free.idx, drop = FALSE]
      #PI <- V22 + B %*% V12 + V21 %*% t(B) + B %*% V11 %*% t(B)
      #A.inv <- solve(A)
      A.inv <-  -1 * I.22.inv
      #VCOV <- A.inv %*% PI %*% t(A.inv)
      V2 <- 1/N * (A.inv %*% V22 %*% A.inv)
      V1 <- 1/N * (A.inv %*% (B %*% V12 + V21 %*% t(B) + B %*% V11 %*% t(B)) %*% A.inv)
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
      # no alpha correction
      VCOV <- V2 + V1
    }

    # store in out
    out$V2 <- V2
    out$V1 <- V1
    out$VCOV <- VCOV
  } # twostep

  # store se
  out$se <- lavoptions$se # in case it changed

  if (lav_verbose()) {
    cat("done.\n")
  }

  out
}
