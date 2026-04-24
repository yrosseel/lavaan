
# fsr = wrapper for local sam
# TODO


sam <- function(model = NULL,                      # nolint start
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
                bootstrap = list(R = 1000L, type = "ordinary",
                                      show.progress = FALSE),
                output = "lavaan",
                bootstrap.args = bootstrap) {        # nolint end
  # "bootstrap" is the new way to specify arguments, replacing bootstrap.args
  if (!missing(bootstrap.args)) {
    lav_msg_warn(gettext(
      "'bootstrap.args' is deprecated; please use 'bootstrap' instead."))
    bootstrap <- bootstrap.args
  }

  # check model= argument
  has_sam_object_flag <- FALSE
  if (inherits(model, "lavaan") && !is.null(model@internal$sam.method)) {
    has_sam_object_flag <- TRUE
  }

  # check sam.method
  sam_method <- tolower(sam.method)
  if (!sam_method %in% c("local", "global", "fsr", "cfsr")) {
    lav_msg_stop(gettextf("unknown option for sam_method: [%s]",
                          "available options are local, global, fsr and cfs."))
  }

  # ------------- handling of warn/debug/verbose switches ----------
  dotdotdot <- list(...)
  if (length(dotdotdot) > 0L) {
    if (!is.null(dotdotdot$debug)) {
     current_debug <- lav_debug()
      if (lav_debug(dotdotdot$debug))
        on.exit(lav_debug(current_debug), TRUE)
      dotdotdot$debug <- NULL
      if (lav_debug()) {
        dotdotdot$warn <- TRUE       # force warnings if debug
        dotdotdot$verbose <- TRUE    # force verbose if debug
      }
    }
    if (!is.null(dotdotdot$warn)) {
      current_warn <- lav_warn()
      if (lav_warn(dotdotdot$warn))
        on.exit(lav_warn(current_warn), TRUE)
      dotdotdot$warn <- NULL
    }
    if (!is.null(dotdotdot$verbose)) {
      current_verbose <- lav_verbose()
      if (lav_verbose(dotdotdot$verbose))
        on.exit(lav_verbose(current_verbose), TRUE)
      dotdotdot$verbose <- NULL
    }
    # check for conditional.x= argument
    # if (!sam_method == "global" && !is.null(dotdotdot$conditional.x) &&
    #     dotdotdot$conditional.x) {
    #   lav_msg_warn(gettext(
    #     "local sam() does not support conditional.x = TRUE (yet) ->
    #             switching to conditional.x = FALSE"))
    #   dotdotdot$conditional.x <- FALSE
    # }
    # check for orthogonal= argument
    if (!is.null(dotdotdot$orthogonal) &&
        dotdotdot$orthogonal &&
        sam_method != "global") {
      lav_msg_warn(gettext(
        "local sam does not support orthogonal = TRUE -> switching to
         global sam"))
     sam_method <- "global"
    }
  } # length(dotdotdot) > 0L

  # check output= argument
  output <- tolower(output)
  if (output %in% c("list", "list.step1.only", "lavaan")) {
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
    } else if (se %in% c("two-step-robust", "two-step.robust",
                        "two_step_robust", "two.step.robust",
                        "twostep.robust")) {
      se <- "twostep.robust"
    } else if (se %in% c("ij", "local")) {
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
      if (!sam_method %in% c("local", "fsr", "cfsr")) { # for now
        lav_msg_stop(gettext(
          "local se only available if sam.method is local, fsr, or cfsr"))
      }
    }
  }
  # default is twostep

  # check for gamma.unbiased
  if (is.null(dotdotdot$gamma.unbiased)) {
     # put in TRUE in dotdotdot # lavaan default is still FALSE
     dotdotdot$gamma.unbiased <- TRUE
  }


  ###############################################
  # STEP 0: process full model, without fitting #
  ###############################################
  if (has_sam_object_flag) {
    fit <- model
    # restore options
    fit@Options <- fit@internal$sam.lavoptions
    # extract other argments from FIT@internal, unless specified as arguments
    if (missing(mm.list)) {
      mm_list <- fit@internal$sam.mm.list
    } else {
      mm_list <- mm.list
    }
    if (missing(mm.args)) {
      mm_args <- fit@internal$sam.mm.args
    } else {
      mm_args <- mm.args
    }
    if (missing(struc.args)) {
      struc_args <- fit@internal$sam.struc.args
    } else {
      struc_args <- struc.args
    }
    if (missing(sam.method)) {
      sam_method <- fit@internal$sam.method
    } else {
      sam_method <- sam.method
    }
    if (missing(local.options)) {
      local_options <- fit@internal$sam.local.options
    } else {
      local_options <- local.options
    }
    if (missing(global.options)) {
      global_options <- fit@internal$sam.global.options
    } else {
      global_options <- global.options
    }
    if (missing(se)) {
      se <- fit@Options$se
    } else {
      fit@Options$se <- se
    }
    if (missing(cmd)) {
      cmd <- fit@internal$sam.cmd
    }
    # remove @internal slot
    fit@internal <- list()
  } else {
    mm_list <- mm.list
    mm_args <- mm.args
    struc_args <- struc.args
    local_options <- local.options
    global_options <- global.options
    sam_method <- sam.method
    fit <- lav_sam_step0(
      cmd = cmd, model = model, data = data, se = se,
      sam.method = sam_method, dotdotdot = dotdotdot
    )

    # check for data.type == "none"
    if (fit@Data@data.type == "none") {
      # we are done; perhaps we only wished to create a FIT object?
      return(fit)
      #lav_msg_stop(gettext("no data or sample statistics are provided."))
    }

    # check if we have categorical data
    if (fit@Model@categorical) {
      # switch to M.method = "ULS"
      local_options[["M.method"]] <- "ULS"
      # if sam_method = "global", force estimator to DWLS in struc par
      if (sam_method == "global" &&
          !is.null(struc_args[["estimator"]]) &&
          struc_args[["estimator"]] == "ML") {
          struc_args[["estimator"]] <- "DWLS"
      }
    }

    # check if we have latent interactions
    lv_interaction_flag <- FALSE
    if (length(unlist(fit@pta$vnames$lv.interaction)) > 0L) {
      lv_interaction_flag <- TRUE
      if (!se %in% c("none", "bootstrap")) {
        se <- "local"
        fit@Options$se <- se
      }
    }

    # check for multiple groups/blocks
    if (fit@Model@nblocks > 1L && se == "local") {
      lav_msg_stop(gettext("se = \"local\" not available (yet) if multiple
                            blocks (groups, levels) are involved."))
    }
  }

  lavoptions <- lavInspect(fit, "options")
  if (lav_verbose()) {
    cat("This is sam using sam.method = ", sam_method, ".\n", sep = "")
  }

  ##############################################
  # STEP 1: fit each measurement model (block) #
  ##############################################
  if (lav_verbose()) {
    cat("Fitting the measurement part:\n")
  }
  step1 <- lav_sam_step1(
    cmd = cmd, mm.list = mm_list, mm.args = mm_args,
    FIT = fit, sam.method = sam_method
  )

  ##################################################
  # STEP 1b: compute Var(eta) and E(eta) per block #
  #          only needed for local approach!       #
  ##################################################
  if (sam_method %in% c("local", "fsr", "cfsr")) {
    # default local_options
    local_opt <- list(
      M.method = "ML",
      lambda.correction = TRUE,
      alpha.correction = 0L,
      twolevel.method = "h1"
    )
    local_options <- modifyList(local_opt, local_options,
      keep.null = FALSE
    )

    # collect COV/YBAR sample statistics per block from FIT
    out <- lav_sam_get_cov_ybar(FIT = fit, local.options = local_options)
    step1$COV  <- out$COV
    step1$YBAR <- out$YBAR

    # compute EETA/VETA
    step1 <- lav_sam_step1_local(
      STEP1 = step1, FIT = fit,
      sam.method = sam_method,
      local.options = local_options,
      return.cov.iveta2 = (se %in% c("local", "local.nt"))
    )
  }

  ##################################################
  # STEP 1c: jacobian of vech(VETA) = f(vech(S))   #
  #          only needed for local approach!       #
  #          only if se = "local"                  #
  ##################################################
  if (se %in% c("local", "local.nt")) {
    gamma_eta <- vector("list", length = fit@Data@ngroups)
    if (lv_interaction_flag) {
      for (g in seq_len(fit@Data@ngroups)) { # group or block
        # initial Gamma.eta
        gamma_eta_init <- step1$COV.IVETA2[[g]]
        # compute 'additional variability' due to step1
        gamma_eta_add <- lav_sam_gamma_add(STEP1 = step1, FIT = fit, group = g)
        gamma_eta[[g]] <- gamma_eta_init + gamma_eta_add
      }
    } else {
      jac <- lav_sam_step1_local_jac(STEP1 = step1, FIT = fit)
      if (se == "local") {
        gamma <- fit@SampleStats@NACOV
      } else if (se == "local.nt") {
        gamma <- lav_object_gamma(lavobject = fit, ADF = FALSE)
      }

      for (g in seq_len(fit@Data@ngroups)) {
        gamma_eta[[g]] <- jac[[g]] %*% gamma[[g]] %*% t(jac[[g]])
      }
      step1$JAC <- jac
    } # no lv-interaction
    step1$Gamma.eta <- gamma_eta
  }

  if (output == "list.step1.only") {
    # stop here, return interim results
    return(step1)
  }

  ####################################
  # STEP 2: estimate structural part #
  ####################################

  step2 <- lav_sam_step2(
    STEP1 = step1, FIT = fit,
    sam.method = sam_method, struc.args = struc_args
  )

  # make sure step1.free.idx and step2.free.idx are disjoint
  both_idx <- which(step1$step1.free.idx %in% step2$step2.free.idx)
  if (length(both_idx) > 0L) {
    step1$step1.free.idx <- step1$step1.free.idx[-both_idx]
    # STEP1$Sigma.11[both.idx,] <- 0
    # STEP1$Sigma.11[,both.idx] <- 0
    if (!is.null(step1$Sigma.11)) { # maybe se = "none"...
      step1$Sigma.11 <- step1$Sigma.11[-both_idx, -both_idx]
    }
  }

  if (output == "list" && lavoptions$se == "none") {
    return(c(step1, step2))
  }

  ################################################################
  # Step 3: assemble results in a 'dummy' JOINT model for output #
  ################################################################
  if (lav_verbose()) {
    cat("Creating JOINT lavaan object ... ")
  }
  joint <- lav_sam_step3_joint(
    FIT = fit, PT = step2$PT,
    sam.method = sam_method
  )
  # fill information from FIT.PA
  joint@Options$optim.method <- step2$FIT.PA@Options$optim.method
  joint@Model@estimator <- fit@Options$estimator # could be DWLS!
  if (sam_method %in% c("local", "fsr", "cfsr")) {
    joint@optim <- step2$FIT.PA@optim
    joint@test <- step2$FIT.PA@test
  }
  # fill in vcov/se information from step 1
  if (!lavoptions$se %in% c("none", "bootstrap")) {
    joint@Options$se <- lavoptions$se # naive/twostep/none/local/bootstrap
    if (joint@Model@ceq.simple.only) {
      vcov_all <- matrix(
        0, joint@Model@nx.unco,
        joint@Model@nx.unco
      )
    } else {
      vcov_all <- matrix(
        0, joint@Model@nx.free,
        joint@Model@nx.free
      )
    }
    vcov_all[step1$step1.free.idx, step1$step1.free.idx] <- step1$Sigma.11
    joint@vcov <- list(
      se = lavoptions$se,
      information = lavoptions$information[1],
      vcov = vcov_all
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
    sam_1 <- lav_sam_table(
      JOINT = joint, STEP1 = step1,
      FIT.PA = step2$FIT.PA,
      cmd = cmd, lavoptions = fit@Options,
      mm.args = mm_args,
      struc.args = struc_args,
      sam.method = sam_method,
      local.options = local_options,
      global.options = global_options
    )
    sam_object <- joint
    sam_object@internal <- sam_1
    default_args <- list(R = 1000L, type = "ordinary",
                         show.progress = FALSE,
                         check.post = TRUE, keep.idx = FALSE)
    this_args <- modifyList(default_args, bootstrap)
    coef_1 <- lav_bootstrap_internal(object = sam_object,
      r = this_args$R, show_progress = this_args$show.progress,
      type = this_args$type, fun = "coef",
      check_post = this_args$check.post, keep_idx = this_args$keep.idx)
    coef_orig <- coef_1
    error_idx <- attr(coef_1, "error.idx")
    nfailed <- length(error_idx) # zero if NULL
    if (nfailed > 0L) {
      lav_msg_warn(gettextf(
        "%s bootstrap runs failed or did not converge.", nfailed))
    }
    notok <- length(attr(coef_1, "nonadmissible")) # zero if NULL
    if (notok > 0L) {
      lav_msg_warn(gettextf(
        "%s bootstrap runs resulted in nonadmissible solutions.", notok))
    }
    if (length(error_idx) > 0L) {
      # new in 0.6-13: we must still remove them!
      coef_1 <- coef_1[-error_idx, , drop = FALSE]
      # this also drops the attributes
    }
    nboot <- nrow(coef_1)
    var_cov <- cov(coef_1) * (nboot - 1) / nboot
    joint@boot$coef <- coef_orig
    joint@Options$bootstrap <- this_args$R
    vcov_1 <- list(se = "bootstrap", VCOV = var_cov)

    if (lav_verbose()) {
      cat("done.\n")
    }

  # analytic twostep/standard/naive/local
  } else {
    vcov_1 <- lav_sam_step2_se(
      FIT = fit, JOINT = joint, STEP1 = step1,
      STEP2 = step2, local.options = local_options
    )
  }

  # fill in twostep standard errors
  if (lavoptions$se != "none") {
    pt_1 <- joint@ParTable
    joint@Options$se <- vcov_1$se
    joint@vcov$se <- vcov_1$se
    if (lavoptions$se == "bootstrap") {
      joint@vcov$vcov <- vcov_1$VCOV
    } else {
      joint@vcov$vcov[step2$step2.free.idx, step2$step2.free.idx] <- vcov_1$VCOV
    }
    pt_1$se <- lav_model_vcov_se(
      lavmodel = joint@Model,
      lavpartable = pt_1,
      VCOV = joint@vcov$vcov
    )
    joint@ParTable <- pt_1
  }



  ##################
  # Step 5: Output #
  ##################

  # assemble pieces to assemble final lavaan object
  if (output == "lavaan") {
    if (lav_verbose()) {
      cat("Assembling results for output ... ")
    }
    sam_1 <- lav_sam_table(
      JOINT = joint, STEP1 = step1,
      FIT.PA = step2$FIT.PA,
      cmd = cmd, lavoptions = fit@Options,
      mm.args = mm_args,
      struc.args = struc_args,
      sam.method = sam_method,
      local.options = local_options,
      global.options = global_options
    )
    res <- joint
    res@internal <- sam_1
    if (lav_verbose()) {
      cat("done.\n")
    }
  } else {
    res <- c(step1, step2, vcov_1)
  }

  if (lav_verbose()) {
    cat("End of sam.\n")
  }

  res
}
