# sam: fit a SEM in two (or more) steps: Structural After Measurement
#
# sam.method = "local"  : compute (corrected) VETA/EETA per block, and use
#                         these summary statistics to fit the structural part
# sam.method = "global" : fix the measurement parameters, and estimate the
#                         structural parameters in the full (joint) model
# sam.method = "fsr"    : naive (uncorrected) factor score regression
# sam.method = "cfsr"   : Croon-corrected factor score regression

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
    lav_msg_stop(gettextf(
      "unknown option for sam.method: [%s]; available options are
       local, global, fsr and cfsr.", sam.method))
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
  if (!output %in% c("list", "list.step1.only", "lavaan")) {
    lav_msg_stop(gettext(
      "output should be \"lavaan\", \"list\" or \"list.step1.only\"."))
  }

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
    # extract other arguments from FIT@internal, unless specified as arguments
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
    } # else: sam_method already holds the validated sam.method
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
    fit <- lav_sam_step0(
      cmd = cmd, model = model, data = data, se = se,
      sam_method = sam_method, dotdotdot = dotdotdot
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

    # check for multiple blocks: multigroup (single level) is supported for
    # the local approach as long as the measurement model has no across-group
    # constraints (checked later in lav_sam_step1_local_jac_mg()); multilevel
    # local SEs are not available yet
    if (fit@Data@nlevels > 1L && se %in% c("local", "local.nt")) {
      lav_msg_stop(gettext("se = \"local\" not available (yet) if multiple
                            levels are involved."))
    }
  }

  lavoptions <- lavInspect(fit, "options")
  if (lav_verbose()) {
    cat("This is sam using sam.method = ", sam_method, ".\n", sep = "")
  }

  ##############################################
  # STEP 1: fit each measurement model (block) #
  ##############################################
  step1 <- lav_sam_step1(
    cmd = cmd, mm_list = mm_list, mm_args = mm_args,
    fit = fit, sam_method = sam_method
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
    out <- lav_sam_get_cov_ybar(fit = fit, local_options = local_options)
    step1$COV  <- out$COV
    step1$YBAR <- out$YBAR

    # compute EETA/VETA
    step1 <- lav_sam_step1_local(
      step1 = step1, fit = fit,
      sam_method = sam_method,
      local_options = local_options,
      return_cov_iveta2 = (se %in% c("local", "local.nt"))
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
        gamma_eta_add <- lav_sam_gamma_add(step1 = step1, fit = fit, group = g)
        gamma_eta[[g]] <- gamma_eta_init + gamma_eta_add
      }
    } else {
      jac <- lav_sam_step1_local_jac(step1 = step1, fit = fit)
      if (se == "local") {
        gamma <- fit@SampleStats@NACOV
      } else if (se == "local.nt") {
        gamma <- lav_object_gamma(lavobject = fit, adf = FALSE)
      }

      for (g in seq_len(fit@Data@ngroups)) {
        gamma_eta[[g]] <- jac[[g]] %*% gamma[[g]] %*% t(jac[[g]])
      }
      step1$JAC <- jac

      # Case B (across-group constraints in the measurement model): the
      # per-group Gamma.eta blocks above are incomplete, because vech(VETA_g)
      # also depends on the other groups' sample statistics. Build the full
      # (cross-group) covariance of the stacked structural statistics so that
      # step 2 can compute a proper cross-group sandwich. Cov(vech(S_g)) =
      # gamma[[g]] / n_g, so the actual covariance of the stacked structural
      # statistics is jac.full %*% bdiag(gamma_g / n_g) %*% t(jac.full).
      if (isTRUE(attr(jac, "caseB"))) {
        jac_full <- attr(jac, "jac.full")
        nobs_g <- unlist(fit@SampleStats@nobs)
        gscaled <- lapply(seq_len(fit@Data@ngroups),
                          function(g) gamma[[g]] / nobs_g[g])
        gamma_obs_full <- lav_mat_bdiag(gscaled)
        step1$Gamma.eta.full <- jac_full %*% gamma_obs_full %*% t(jac_full)
        step1$caseB <- TRUE
      }
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
    step1 = step1, fit = fit,
    sam_method = sam_method, struc_args = struc_args
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
    fit = fit, pt_1 = step2$PT,
    sam_method = sam_method
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

  # assemble the extra information for the @internal slot
  # (needed for output = "lavaan", and for se = "bootstrap", where
  #  lav_bootstrap_internal() needs a full sam object)
  if (output == "lavaan" || lavoptions$se == "bootstrap") {
    sam_1 <- lav_sam_table(
      joint = joint, step1 = step1,
      fit_pa = step2$FIT.PA,
      cmd = cmd, lavoptions = fit@Options,
      mm_args = mm_args,
      struc_args = struc_args,
      sam_method = sam_method,
      local_options = local_options,
      global_options = global_options
    )
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
    sam_object <- joint
    sam_object@internal <- sam_1
    boot_out <- lav_sam_step2_se_bootstrap(
      sam_object = sam_object, bootstrap = bootstrap
    )
    joint@boot$coef <- boot_out$boot.coef
    joint@Options$bootstrap <- boot_out$R
    vcov_1 <- list(se = "bootstrap", VCOV = boot_out$VCOV)

    if (lav_verbose()) {
      cat("done.\n")
    }

  # analytic twostep/standard/naive/local
  } else {
    vcov_1 <- lav_sam_step2_se(
      fit = fit, joint = joint, step1 = step1,
      step2 = step2, local_options = local_options
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

    # Monte Carlo samples for defined parameter CIs (Preacher & Selig 2012)
    mc_coef <- NULL
    if (sum(pt_1$op == ":=") > 0L && lavoptions$se != "bootstrap" &&
      lav_model_vcov_se_mc_active(lavoptions)) {
      mc_coef <- lav_model_vcov_mc(lavmodel = joint@Model,
                                   VCOV = joint@vcov$vcov,
                                   lavoptions = lavoptions)
      joint@internal$monte.carlo <- list(coef = mc_coef)
    }

    pt_1$se <- lav_model_vcov_se(
      lavmodel = joint@Model,
      lavpartable = pt_1,
      VCOV = joint@vcov$vcov,
      MC = mc_coef,
      lavoptions = lavoptions
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
    res <- joint
    mc_keep <- joint@internal$monte.carlo
    res@internal <- sam_1
    if (!is.null(mc_keep)) {
      res@internal$monte.carlo <- mc_keep
    }
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
