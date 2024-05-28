# SAM: a Structural After Measurement approach
#
# Yves Rosseel & Wen-Wei Loh, Feb-May 2019

# local vs global sam
# local sam = alternative for FSR+Croon
# - but no need to compute factor scores or corrections
# gloal sam = (old) twostep
# - but we can also take a 'local' perspective

# restrictions:
#
# local and global:
#  - all (measured) latent variables must have indicators that are observed
# local:
#  - only if LAMBDA is of full column rank (eg no SRM, no bi-factor, no MTMM)
#  - if multiple groups: each group has the same set of latent variables!
#  - global approach is used to compute corrected two-step standard errors

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
                sam.method = "local", # or "global", or "fsr"
                ..., # common options
                local.options = list(
                  M.method = "ML", # mapping matrix
                  lambda.correction = TRUE,
                  alpha.correction = 0L, # 0 -> (N-1)
                  twolevel.method = "h1"
                ),
                # h1, anova, mean
                global.options = list(), # not used for now
                output = "lavaan") {
  # output
  output <- tolower(output)
  if (output == "list" || output == "lavaan") {
    # nothing to do
  } else {
    lav_msg_stop(gettext("output should be \"list\" or \"lavaan.\""))
  }

  # check se= argument
  if (!se %in% c("standard", "naive", "twostep", "twostep2",
                 "bootstrap", "none")) {
    lav_msg_stop(gettext(
      "se= argument must be twostep, bootstrap, naive, standard or none."))
  }

  # handle dot dot dot
  dotdotdot <- list(...)

  if (!is.null(dotdotdot$conditional.x)) {
    lav_msg_warn(gettext(
      "sam() does not support conditional.x = TRUE (yet) -> switching to
      conditional.x = FALSE"))
    dotdotdot$conditional.x <- FALSE
  }

  ###############################################
  # STEP 0: process full model, without fitting #
  ###############################################
  FIT <- lav_sam_step0(
    cmd = cmd, model = model, data = data, se = se,
    sam.method = sam.method, dotdotdot = dotdotdot
  )

  # check for data.type == "none"
  if (FIT@Data@data.type == "none") {
    lav_msg_stop(gettext("no data or sample statistics are provided."))
  }

  lavoptions <- lavInspect(FIT, "options")
  if (lavoptions$verbose) {
    cat("This is sam using sam.method = ", sam.method, ".\n", sep = "")
  }

  ##############################################
  # STEP 1: fit each measurement model (block) #
  ##############################################
  if (lavoptions$verbose) {
    cat("Fitting the measurement part:\n")
  }
  STEP1 <- lav_sam_step1(
    cmd = cmd, mm.list = mm.list, mm.args = mm.args,
    FIT = FIT, data = data, sam.method = sam.method
  )

  ##################################################
  # STEP 1b: compute Var(eta) and E(eta) per block #
  #          only needed for local approach!       #
  ##################################################
  if (sam.method %in% c("local", "fsr")) {
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

    STEP1 <- lav_sam_step1_local(
      STEP1 = STEP1, FIT = FIT,
      sam.method = sam.method,
      local.options = local.options
    )
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
  if (lavoptions$verbose) {
    cat("Assembling results for output ... ")
  }
  JOINT <- lav_sam_step3_joint(
    FIT = FIT, PT = STEP2$PT,
    sam.method = sam.method
  )
  # fill information from FIT.PA
  JOINT@Options$optim.method <- STEP2$FIT.PA@Options$optim.method
  JOINT@Model@estimator <- FIT@Options$estimator # could be DWLS!
  if (sam.method %in% c("local", "fsr")) {
    JOINT@optim <- STEP2$FIT.PA@optim
    JOINT@test <- STEP2$FIT.PA@test
  }
  # fill in vcov/se information from step 1
  if (lavoptions$se != "none") {
    JOINT@Options$se <- lavoptions$se # naive/twostep/none
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

  if (lavoptions$verbose) {
    cat("done.\n")
  }


  ##############################################
  # Step 4: compute standard errors for step 2 #
  ##############################################

  if (lavoptions$se == "bootstrap") {
    #
	# TODO
	#
    #VCOV <- lav_sam_boot_se(FIT = FIT, JOINT = JOINT, STEP1 = STEP1,
    #  STEP2 = STEP2, local.options = local.options
    #)
	#
  } else {
    VCOV <- lav_sam_step2_se(
      FIT = FIT, JOINT = JOINT, STEP1 = STEP1,
      STEP2 = STEP2, local.options = local.options
    )
    # fill in twostep standard errors
    if (lavoptions$se != "none") {
      PT <- JOINT@ParTable
      JOINT@Options$se <- lavoptions$se
      JOINT@vcov$se <- lavoptions$se
      JOINT@vcov$vcov[STEP2$step2.free.idx, STEP2$step2.free.idx] <- VCOV$VCOV
      PT$se <- lav_model_vcov_se(
        lavmodel = JOINT@Model,
        lavpartable = PT,
        VCOV = JOINT@vcov$vcov
      )
      JOINT@ParTable <- PT
    }
  }



  ##################
  # Step 5: Output #
  ##################

  # assemble pieces to assemble final lavaan object
  if (output == "lavaan") {
    if (lavoptions$verbose) {
      cat("Assembling results for output ... ")
    }
    SAM <- lav_sam_table(
      JOINT = JOINT, STEP1 = STEP1,
      FIT.PA = STEP2$FIT.PA,
      mm.args = mm.args,
      struc.args = struc.args,
      sam.method = sam.method,
      local.options = local.options,
      global.options = global.options
    )
    res <- JOINT
    res@internal <- SAM
    if (lavoptions$verbose) {
      cat("done.\n")
    }
  } else {
    res <- c(STEP1, STEP2, VCOV)
  }

  if (lavoptions$verbose) {
    cat("End of sam.\n")
  }

  res
}
