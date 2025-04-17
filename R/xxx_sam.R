# SAM: a Structural After Measurement approach
#
# Yves Rosseel & Wen-Wei Loh, Feb-May 2019

# local vs global sam
# local sam = alternative for FSR+Croon
# - but no need to compute factor scores or corrections
# global sam = (old) twostep
# - but we can also take a 'local' perspective

# restrictions:
#
# local and global:
#  - all (measured) latent variables must have indicators that are observed
#    update: higher-order measurement models are supported in local SAM (0.6-20)
# local:
#  - only if LAMBDA is of full column rank (eg no SRM, no bi-factor, no MTMM)
#  - if multiple groups: each group has the same set of latent variables!
#  - global approach is used to compute corrected two-step standard errors
#    update: se = "local" is truely local, and uses 'Gamma' as an additional
#            ingredient; Gamma reflects the sampling variability of the
#            sample statistics (VETA and EETA)

# YR 12 May 2019 - first version
# YR 22 May 2019 - merge sam/twostep (call it 'local' vs 'global' sam)

# YR 27 June 2021 - prepare for `public' release
#                 - add Fuller (1987) 'lambda' correction if (MSM - MTM) is not
#                   positive definite
#                 - se = "none" now works
#                 - store 'local' information in @internal slot (for printing)

# YR 16 Oct 2021 - if an indicator is also a predictor/outcome in the
#                  structural part, treat it as an observed predictor
#                  without measurement error in the second step
#                  (ie, set THETA element to zero)

# YR 03 Dec 2022 - allow for sam.method = "fsr" and se = "naive"
#                - add alpha.correction= argument (for small sample correction)

# YR 21 May 2023 - allow for latent quadratic/interaction terms in the
#                  structural part (assuming the errors are normal, for now)

# YR 25 May 2023 - restructure code into multiple files
#                - rename veta.force.pd -> lambda.correction
#                - move alpha.correction= argument to local.options

# YR 09 Nov 2024 - add se = "bootstrap"
# YR 14 Nov 2024 - add se = "local"
# YR 01 Mar 2025 - allow for higher-order measurement models in local SAM
# YR 25 Mar 2025 - add sam.method = "cFSR" (using correlation-preserving FS)
# YR 15 Apr 2025 - add se = "twostep.robust"

# twostep = wrapper for global sam
twostep <- function(model = NULL, data = NULL, cmd = "sem",
                    mm.list = NULL, mm.args = list(), struc.args = list(),
                    ..., # global options
                    output = "lavaan") {
  sam(
    model = model, data = data, cmd = cmd, mm.list = mm.list,
    mm.args = mm.args, struc.args = struc.args,
    sam.method = "global", # or global
    ..., # global options
    output = output
  )
}

# fsr = wrapper for local sam
# TODO


sam <- function(model = NULL,
                data = NULL,
                cmd = "sem",
                se = "twostep",
                mm.list = NULL,
                mm.args = list(bounds = "wide.zerovar"),
                struc.args = list(estimator = "ML"),
                sam.method = "local", # or "global", or "fsr", or "cfsr"
                ..., # common options
                local.options = list(
                  M.method = "ML", # mapping matrix
                  lambda.correction = TRUE,
                  alpha.correction = 0L, # 0 -> (N-1)
                  twolevel.method = "h1"
                ),
                # h1, anova, mean
                global.options = list(), # not used for now
                bootstrap.args = list(R = 1000L, type = "ordinary",
                                      show.progress = FALSE),
                output = "lavaan") {

  # check model= argument
  has.sam.object.flag <- FALSE
  if (inherits(model, "lavaan") && !is.null(model@internal$sam.method)) {
    has.sam.object.flag <- TRUE
  }

  # check sam.method
  sam.method <- tolower(sam.method)
  if (!sam.method %in% c("local", "global", "fsr", "cfsr")) {
    lav_msg_stop(gettextf("unknown option for sam.method: [%s]",
                          "available options are local, global, fsr and cfs."))
  }

  # ------------- handling of warn/debug/verbose switches ----------
  dotdotdot <- list(...)
  if( length(dotdotdot) > 0L) {
    if (!is.null(dotdotdot$debug)) {
     current.debug <- lav_debug()
      if (lav_debug(dotdotdot$debug))
        on.exit(lav_debug(current.debug), TRUE)
      dotdotdot$debug <- NULL
      if (lav_debug()) {
        dotdotdot$warn <- TRUE       # force warnings if debug
        dotdotdot$verbose <- TRUE    # force verbose if debug
      }
    }
    if (!is.null(dotdotdot$warn)) {
      current.warn <- lav_warn()
      if (lav_warn(dotdotdot$warn))
        on.exit(lav_warn(current.warn), TRUE)
      dotdotdot$warn <- NULL
    }
    if (!is.null(dotdotdot$verbose)) {
      current.verbose <- lav_verbose()
      if (lav_verbose(dotdotdot$verbose))
        on.exit(lav_verbose(current.verbose), TRUE)
      dotdotdot$verbose <- NULL
    }
    # check for conditional.x= argument
    # if (!sam.method == "global" && !is.null(dotdotdot$conditional.x) &&
    #     dotdotdot$conditional.x) {
    #   lav_msg_warn(gettext(
    #     "local sam() does not support conditional.x = TRUE (yet) -> switching to
    #       conditional.x = FALSE"))
    #   dotdotdot$conditional.x <- FALSE
    # }
    # check for orthogonal= argument
    if (!is.null(dotdotdot$orthogonal) &&
        dotdotdot$orthogonal &&
        sam.method != "global") {
      lav_msg_warn(gettext(
        "local sam does not support orthogonal = TRUE -> switching to
         global sam"))
     sam.method <- "global"
    }
  } # length(dotdotdot) > 0L

  # check output= argument
  output <- tolower(output)
  if (output == "list" || output == "lavaan") {
    # nothing to do
  } else {
    lav_msg_stop(gettext("output should be \"list\" or \"lavaan.\""))
  }

  #

  # check se= argument
  if (!missing(se)) {
    se <- tolower(se)
    # aliases
    if (se %in% c("two-step", "two_step", "two.step")) {
      se <- "twostep"
    } else if(se %in% c("two-step-robust", "two-step.robust",
                        "two_step_robust", "two.step.robust",
                        "twostep.robust")) {
      se <- "twostep.robust"
    }
    else if (se %in% c("ij", "local")) {
      se <- "local"
    }
    # check if valid
    if (!se %in% c("standard", "naive", "twostep", "local", "local.nt",
                   "twostep.robust", "bootstrap", "none")) {
      lav_msg_stop(gettext(
        "se= argument must be twostep, twostep.robust, bootstrap, or local"))
    }
    # check for local
    if (se %in% c("local", "local.nt")) {
      if (!sam.method %in% c("local", "fsr", "cfsr")) { # for now
        lav_msg_stop(gettext("local se only available if sam.method is local, fsr, or cfsr"))
      }
    }
  }
  # default is twostep


  ###############################################
  # STEP 0: process full model, without fitting #
  ###############################################
  if (has.sam.object.flag) {
    FIT <- model
    # restore options
    FIT@Options <- FIT@internal$sam.lavoptions
    # extract other argments from FIT@internal, unless specified as arguments
    if (missing(mm.list)){
      mm.list <- FIT@internal$sam.mm.list
    }
    if (missing(mm.args)){
      mm.args <- FIT@internal$sam.mm.args
    }
    if (missing(struc.args)) {
      struc.args <- FIT@internal$sam.struc.args
    }
    if (missing(sam.method)) {
      sam.method <- FIT@internal$sam.method
    }
    if (missing(local.options)) {
      local.options <- FIT@internal$sam.local.options
    }
    if (missing(global.options)) {
      global.options <- FIT@internal$sam.global.options
    }
    if (missing(se)) {
      se <- FIT@Options$se
    } else {
      FIT@Options$se <- se
    }
    if (missing(cmd)) {
      cmd <- FIT@internal$sam.cmd
    }
    # remove @internal slot
    FIT@internal <- list()
  } else {
    FIT <- lav_sam_step0(
      cmd = cmd, model = model, data = data, se = se,
      sam.method = sam.method, dotdotdot = dotdotdot
    )

    # check for data.type == "none"
    if (FIT@Data@data.type == "none") {
      # we are done; perhaps we only wished to create a FIT object?
      return(FIT)
      #lav_msg_stop(gettext("no data or sample statistics are provided."))
    }

    # check if we have categorical data
    if (FIT@Model@categorical) {
      # switch to M.method = "ULS"
      local.options[["M.method"]] <- "ULS"
      # if sam.method = "global", force estimator to DWLS in struc par
  	  if (sam.method == "global" &&
  	      !is.null(struc.args[["estimator"]]) &&
    	    struc.args[["estimator"]] == "ML") {
          struc.args[["estimator"]] <- "DWLS"
      }
    }
  }

  lavoptions <- lavInspect(FIT, "options")
  if (lav_verbose()) {
    cat("This is sam using sam.method = ", sam.method, ".\n", sep = "")
  }

  ##############################################
  # STEP 1: fit each measurement model (block) #
  ##############################################
  if (lav_verbose()) {
    cat("Fitting the measurement part:\n")
  }
  STEP1 <- lav_sam_step1(
    cmd = cmd, mm.list = mm.list, mm.args = mm.args,
    FIT = FIT, sam.method = sam.method
  )

  ##################################################
  # STEP 1b: compute Var(eta) and E(eta) per block #
  #          only needed for local approach!       #
  ##################################################
  if (sam.method %in% c("local", "fsr", "cfsr")) {
    # default local.options
    local.opt <- list(
      M.method = "ML",
      lambda.correction = TRUE,
      alpha.correction = 0L,
      twolevel.method = "h1"
    )
    local.options <- modifyList(local.opt, local.options,
      keep.null = FALSE
    )

    # collect COV/YBAR sample statistics per block from FIT
    out <- lav_sam_get_cov_ybar(FIT = FIT, local.options = local.options)
    STEP1$COV  <- out$COV
    STEP1$YBAR <- out$YBAR

    # compute EETA/VETA
    STEP1 <- lav_sam_step1_local(
      STEP1 = STEP1, FIT = FIT,
      sam.method = sam.method,
      local.options = local.options
    )
  }

  ##################################################
  # STEP 1c: jacobian of vech(VETA) = f(vech(S))   #
  #          only needed for local approach!       #
  #          only if se = "local"                  #
  ##################################################
  if (se %in% c("local", "local.nt")) {
    JAC <- lav_sam_step1_local_jac(STEP1 = STEP1, FIT = FIT)
    if (se == "local") {
      Gamma <- FIT@SampleStats@NACOV
    } else if (se == "local.nt") {
      Gamma <- lav_object_gamma(lavobject = FIT, ADF = FALSE)
    }

    Gamma.eta <- vector("list", length = FIT@Data@ngroups)
    for (g in seq_len(FIT@Data@ngroups)) {
      Gamma.eta[[g]] <- JAC[[g]] %*% Gamma[[g]] %*% t(JAC[[g]])
    }
    STEP1$JAC <- JAC
    STEP1$Gamma.eta <- Gamma.eta
  }

  ####################################
  # STEP 2: estimate structural part #
  ####################################

  STEP2 <- lav_sam_step2(
    STEP1 = STEP1, FIT = FIT,
    sam.method = sam.method, struc.args = struc.args
  )

  # make sure step1.free.idx and step2.free.idx are disjoint
  both.idx <- which(STEP1$step1.free.idx %in% STEP2$step2.free.idx)
  if (length(both.idx) > 0L) {
    STEP1$step1.free.idx <- STEP1$step1.free.idx[-both.idx]
    # STEP1$Sigma.11[both.idx,] <- 0
    # STEP1$Sigma.11[,both.idx] <- 0
    STEP1$Sigma.11 <- STEP1$Sigma.11[-both.idx, -both.idx]
  }

  if (output == "list" && lavoptions$se == "none") {
    return(c(STEP1, STEP2))
  }

  ################################################################
  # Step 3: assemble results in a 'dummy' JOINT model for output #
  ################################################################
  if (lav_verbose()) {
    cat("Creating JOINT lavaan object ... ")
  }
  JOINT <- lav_sam_step3_joint(
    FIT = FIT, PT = STEP2$PT,
    sam.method = sam.method
  )
  # fill information from FIT.PA
  JOINT@Options$optim.method <- STEP2$FIT.PA@Options$optim.method
  JOINT@Model@estimator <- FIT@Options$estimator # could be DWLS!
  if (sam.method %in% c("local", "fsr", "cfsr")) {
    JOINT@optim <- STEP2$FIT.PA@optim
    JOINT@test <- STEP2$FIT.PA@test
  }
  # fill in vcov/se information from step 1
  if (!lavoptions$se %in% c("none", "bootstrap")) {
    JOINT@Options$se <- lavoptions$se # naive/twostep/none/local/bootstrap
    if (JOINT@Model@ceq.simple.only) {
      VCOV.ALL <- matrix(
        0, JOINT@Model@nx.unco,
        JOINT@Model@nx.unco
      )
    } else {
      VCOV.ALL <- matrix(
        0, JOINT@Model@nx.free,
        JOINT@Model@nx.free
      )
    }
    VCOV.ALL[STEP1$step1.free.idx, STEP1$step1.free.idx] <- STEP1$Sigma.11
    JOINT@vcov <- list(
      se = lavoptions$se,
      information = lavoptions$information[1],
      vcov = VCOV.ALL
    )
    # no need to fill @ParTable$se, as step1 SE values should already
    # be in place
  }

  if (lav_verbose()) {
    cat("done.\n")
  }

  ##############################################
  # Step 4: compute standard errors for step 2 #
  ##############################################

  if (lavoptions$se == "bootstrap") {
    if (lav_verbose()) {
      cat("computing VCOV for      se = bootstrap ...")
    }
    # construct temporary sam object, so that lav_bootstrap_internal() can
    # use it
    SAM <- lav_sam_table(
      JOINT = JOINT, STEP1 = STEP1,
      FIT.PA = STEP2$FIT.PA,
      cmd = cmd, lavoptions = FIT@Options,
      mm.args = mm.args,
      struc.args = struc.args,
      sam.method = sam.method,
      local.options = local.options,
      global.options = global.options
    )
    sam_object <- JOINT
    sam_object@internal <- SAM
    default.args <- list(R = 1000L, type = "ordinary",
                         show.progress = FALSE,
                         check.post = TRUE, keep.idx = FALSE)
    this.args <- modifyList(default.args, bootstrap.args)
    COEF <- lav_bootstrap_internal(object = sam_object,
      R = this.args$R, show.progress = this.args$show.progress,
      type = this.args$type, FUN = "coef",
      check.post = this.args$check.post, keep.idx = this.args$keep.idx)
    COEF.orig <- COEF
    error.idx <- attr(COEF, "error.idx")
    nfailed <- length(error.idx) # zero if NULL
    if (nfailed > 0L) {
      lav_msg_warn(gettextf(
        "%s bootstrap runs failed or did not converge.", nfailed))
    }
    notok <- length(attr(COEF, "nonadmissible")) # zero if NULL
    if (notok > 0L) {
      lav_msg_warn(gettextf(
        "%s bootstrap runs resulted in nonadmissible solutions.", notok))
    }
    if (length(error.idx) > 0L) {
      # new in 0.6-13: we must still remove them!
      COEF <- COEF[-error.idx, , drop = FALSE]
      # this also drops the attributes
    }
    nboot <- nrow(COEF)
    VarCov <- cov(COEF) * (nboot - 1) / nboot
    JOINT@boot$coef <- COEF.orig
    JOINT@Options$bootstrap <- this.args$R
    VCOV <- list(se = "bootstrap", VCOV = VarCov)

    if (lav_verbose()) {
      cat("done.\n")
    }

  # analytic twostep/standard/naive/local
  } else {
    VCOV <- lav_sam_step2_se(
      FIT = FIT, JOINT = JOINT, STEP1 = STEP1,
      STEP2 = STEP2, local.options = local.options
    )
  }

  # fill in twostep standard errors
  if (lavoptions$se != "none") {
    PT <- JOINT@ParTable
    JOINT@Options$se <- VCOV$se
    JOINT@vcov$se <- VCOV$se
    if (lavoptions$se == "bootstrap") {
      JOINT@vcov$vcov <- VCOV$VCOV
    } else {
      JOINT@vcov$vcov[STEP2$step2.free.idx, STEP2$step2.free.idx] <- VCOV$VCOV
    }
    PT$se <- lav_model_vcov_se(
      lavmodel = JOINT@Model,
      lavpartable = PT,
      VCOV = JOINT@vcov$vcov
    )
    JOINT@ParTable <- PT
  }



  ##################
  # Step 5: Output #
  ##################

  # assemble pieces to assemble final lavaan object
  if (output == "lavaan") {
    if (lav_verbose()) {
      cat("Assembling results for output ... ")
    }
    SAM <- lav_sam_table(
      JOINT = JOINT, STEP1 = STEP1,
      FIT.PA = STEP2$FIT.PA,
      cmd = cmd, lavoptions = FIT@Options,
      mm.args = mm.args,
      struc.args = struc.args,
      sam.method = sam.method,
      local.options = local.options,
      global.options = global.options
    )
    res <- JOINT
    res@internal <- SAM
    if (lav_verbose()) {
      cat("done.\n")
    }
  } else {
    res <- c(STEP1, STEP2, VCOV)
  }

  if (lav_verbose()) {
    cat("End of sam.\n")
  }

  res
}
