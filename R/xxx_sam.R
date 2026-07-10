# sam: fit a SEM in two (or more) steps: Structural After Measurement
#
# sam_method = "local"  : compute (corrected) VETA/EETA per block, and use
#                         these summary statistics to fit the structural part
# sam_method = "global" : fix the measurement parameters, and estimate the
#                         structural parameters in the full (joint) model
# sam_method = "fsr"    : naive (uncorrected) factor score regression
# sam_method = "cfsr"   : Croon-corrected factor score regression

sam <- function(model = NULL,
                data = NULL,
                aux = NULL,
                cmd = "sem",
                se = "twostep",
                mm_list = NULL,
                mm_args = list(bounds = "wide.zerovar"),
                struc_args = list(estimator = "ML"),
                sam_method = "local", # or "global", or "fsr", or "cfsr"
                ..., # common options
                local_options = list(
                  M.method = "ML", # mapping matrix
                  lambda.correction = TRUE,
                  alpha.correction = 0L, # 0 -> (N-1)
                  twolevel.method = "h1"
                ),
                # h1, anova, mean
                global_options = list(), # not used for now
                bootstrap = list(R = 1000L, type = "ordinary",
                                      show.progress = FALSE),
                output = "lavaan",
                bootstrap_args = bootstrap) {
  # capture the user's call; it is stored in the @call slot of the returned
  # object (if output = "lavaan"), so that update() and getCall() work
  mc <- match.call()
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, TRUE)

  # auxiliary variables: forward via dotdotdot, so they reach the
  # underlying measurement-block (and structural) lavaan() calls
  if (!is.null(aux)) {
    dotdotdot$aux <- aux
  }

  # "bootstrap" is the new way to specify arguments, replacing bootstrap_args
  if (!missing(bootstrap_args)) {
    lav_msg_warn(gettext(
      "'bootstrap_args' is deprecated; please use 'bootstrap' instead."))
    bootstrap <- bootstrap_args
  }
  # plugin compatibility with sem()/cfa(): bootstrap = <number of draws>
  if (!is.list(bootstrap) && is.numeric(bootstrap) &&
      length(bootstrap) == 1L && is.finite(bootstrap)) {
    bootstrap <- list(R = as.integer(bootstrap))
  }
  if (!is.list(bootstrap)) {
    lav_msg_stop(gettext(
      "the bootstrap= argument must be a list (eg list(R = 1000L)) or a
       single number (the number of bootstrap draws)."))
  }

  # check model= argument
  has_sam_object_flag <- FALSE
  if (inherits(model, "lavaan") && !is.null(model@internal$sam.method)) {
    has_sam_object_flag <- TRUE
  }

  # check sam_method
  sam_method <- tolower(sam_method)
  if (!sam_method %in% c("local", "global", "fsr", "cfsr")) {
    lav_msg_stop(gettextf(
      "unknown option for sam_method: [%s]; available options are
       local, global, fsr and cfsr.", sam_method))
  }

  # ------------- handling of warn/debug/verbose switches ----------
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
    } else if (se %in% c("robust", "robust.sem", "robust.huber.white",
                         "robust.cluster", "robust.cluster.sem")) {
      # plugin compatibility with sem()/cfa(): map the robust se= variants
      # to the SAM analogue (the robust two-step correction); clustering is
      # picked up automatically when cluster= is set
      lav_msg_note(gettextf(
        "se = \"%s\" is not a sam() option; using the SAM analogue
         se = \"twostep.robust\" instead.", se))
      se <- "twostep.robust"
    }
    # check if valid
    if (!se %in% c("standard", "naive", "twostep", "local", "local.nt",
                   "twostep.robust", "bootstrap", "none")) {
      lav_msg_stop(gettextf(
        "invalid se= argument (%1$s) for sam(); valid options are %2$s.",
        se,
        lav_msg_view(c("twostep", "twostep.robust", "local", "local.nt",
                       "naive", "standard", "bootstrap", "none"),
                     log_sep = "or")))
    }
    # check for local
    if (se %in% c("local", "local.nt")) {
      if (!sam_method %in% c("local", "fsr", "cfsr")) { # for now
        lav_msg_stop(gettext(
          "local se only available if sam_method is local, fsr, or cfsr"))
      }
    }
  }
  # default is twostep

  # check for gamma.unbiased
  if (is.null(dotdotdot$gamma.unbiased)) {
     # put in TRUE in dotdotdot # lavaan default is still FALSE
     # but not for clustered data or conditional.x = TRUE: the unbiased Gamma
     # is not available there; the (biased) cluster-robust/conditional Gamma
     # is used instead
     dotdotdot$gamma.unbiased <- is.null(dotdotdot$cluster) &&
       !isTRUE(dotdotdot$conditional.x)
  }


  ###############################################
  # STEP 0: process full model, without fitting #
  ###############################################
  if (has_sam_object_flag) {
    fit <- model
    # restore options
    fit@Options <- fit@internal$sam.lavoptions
    # extract other arguments from FIT@internal, unless specified as arguments
    if (missing(mm_list)) {
      mm_list <- fit@internal$sam.mm.list
    } else {
      mm_list <- mm_list
    }
    if (missing(mm_args)) {
      mm_args <- fit@internal$sam.mm.args
    } else {
      mm_args <- mm_args
    }
    if (missing(struc_args)) {
      struc_args <- fit@internal$sam.struc.args
    } else {
      struc_args <- struc_args
    }
    if (missing(sam_method)) {
      sam_method <- fit@internal$sam.method
    } # else: sam_method already holds the validated sam_method
    if (missing(local_options)) {
      local_options <- fit@internal$sam.local.options
    } else {
      local_options <- local_options
    }
    if (missing(global_options)) {
      global_options <- fit@internal$sam.global.options
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
    # local SEs are supported for the single-group continuous complete-data
    # setting (see lav_sam_gamma_eta_2l()), but not the normal-theory variant
    if (fit@Data@nlevels > 1L && se == "local.nt") {
      lav_msg_stop(gettext("se = \"local.nt\" not available if multiple
                            levels are involved; use se = \"local\"."))
    }

    # single-level clustered data: the normal-theory Gamma used by
    # se = "local.nt" does not account for clustering
    if (se == "local.nt" && fit@Data@nlevels == 1L &&
        length(fit@Data@cluster) > 0L) {
      lav_msg_warn(gettext("se = \"local.nt\" uses the normal-theory Gamma,
        which does not account for clustering; consider se = \"local\"
        instead."))
    }
  }

  lavoptions <- lavInspect(fit, "options")
  if (lav_verbose()) {
    cat("This is sam using sam_method = ", sam_method, ".\n", sep = "")
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
    if (fit@Model@conditional.x) {
      step1$SLOPES <- out$SLOPES # res.slopes of y on x
    }

    # compute EETA/VETA
    step1 <- lav_sam_step1_local(
      step1 = step1, fit = fit,
      sam_method = sam_method,
      local_options = local_options,
      return_cov_iveta2 = (se %in% c("local", "local.nt", "twostep.robust"))
    )
  }

  ##################################################
  # STEP 1c: jacobian of vech(VETA) = f(vech(S))   #
  #          Gamma.eta = NACOV of vech(VETA), the   #
  #          ingredient for (1) the local SEs and   #
  #          (2) the corrected two-step STRUCTURAL  #
  #          test (Satorra-Bentler), which is the   #
  #          DEFAULT test for sam.method = local    #
  #          REGARDLESS of the requested SE. So we  #
  #          compute Gamma.eta for every analytic   #
  #          SE of the local family (incl. the      #
  #          default se = "twostep" and "naive").   #
  #          For SEs that do not strictly need it   #
  #          (twostep/twostep.robust/naive) a        #
  #          failure is NON-fatal: we just skip the #
  #          corrected test. local/local.nt require #
  #          it (re-raise on failure).              #
  ##################################################
  if (se %in% c("local", "local.nt", "twostep", "twostep.robust", "naive")) {
    gamma_eta_required <- se %in% c("local", "local.nt")
    ge_try <- tryCatch({
      gamma_eta <- vector("list", length = fit@Data@ngroups)
      if (lv_interaction_flag) {
        for (g in seq_len(fit@Data@ngroups)) { # group or block
          # initial Gamma.eta
          gamma_eta_init <- step1$COV.IVETA2[[g]]
          # compute 'additional variability' due to step1
          gamma_eta_add <- lav_sam_gamma_add(step1 = step1, fit = fit, group = g)
          gamma_eta[[g]] <- gamma_eta_init + gamma_eta_add
        }
      } else if (fit@Data@nlevels > 1L) {
        # two-level: Gamma.eta per level (one entry per pseudo-group of the
        # structural refit) from the cluster-sandwich Gamma of the h1
        # estimates; the cross-level covariance of the stacked statistics
        # goes through the Case B (full-covariance) sandwich in step 2
        tmp_2l <- lav_sam_gamma_eta_2l(step1 = step1, fit = fit)
        gamma_eta <- tmp_2l$Gamma.eta
        step1$JAC <- tmp_2l$JAC
        step1$Gamma.eta.full <- tmp_2l$Gamma.eta.full
        step1$caseB <- TRUE
      } else {
        jac <- lav_sam_step1_local_jac(step1 = step1, fit = fit)

        # build Gamma.eta = JAC %*% Gamma %*% t(JAC) in a memory-lean form that
        # never materializes the p* x p* observed Gamma (see
        # lav_sam_gamma_eta_g()). Falls back to the explicit observed Gamma for
        # settings where that form does not apply (categorical, conditional.x,
        # missing data, ...) and for se = "local.nt" (normal-theory Gamma). The
        # ADF/influence path applies to all SEs except the normal-theory one.
        use_influence <- (se != "local.nt")
        if (use_influence) {
          for (g in seq_len(fit@Data@ngroups)) {
            ge_g <- lav_sam_gamma_eta_g(fit = fit, jac_g = jac[[g]], g = g)
            if (is.null(ge_g)) {
              use_influence <- FALSE
              break
            }
            gamma_eta[[g]] <- ge_g
          }
        }

        gamma <- NULL
        if (!use_influence) {
          if (se == "local.nt") {
            gamma <- lav_object_gamma(lavobject = fit, adf = FALSE)
          } else {
            # stored NACOV, else recompute with the fit-time settings. The
            # (unbiased) ADF Gamma is unavailable for fixed.x / conditional.x
            # (lav_samp_gamma() stops), which the default se = "twostep" can
            # hit (its default is fixed.x = TRUE); fall back to the *biased*
            # ADF Gamma there, exactly as the global yuan.chan test does.
            gamma <- tryCatch(lav_gamma_used(fit), error = function(e) NULL)
            if (is.null(gamma) || is.null(gamma[[1]])) {
              opts <- fit@Options
              opts$gamma.unbiased <- FALSE
              gamma <- lav_gamma_used(
                lavdata = fit@Data, lavoptions = opts,
                lavsamplestats = fit@SampleStats, lavh1 = fit@h1,
                lavimplied = fit@implied
              )
            }
          }
          for (g in seq_len(fit@Data@ngroups)) {
            gamma_eta[[g]] <- jac[[g]] %*% gamma[[g]] %*% t(jac[[g]])
          }
        }
        step1$JAC <- jac

        # Case B (across-group constraints in the measurement model): the
        # per-group Gamma.eta blocks above are incomplete, because vech(VETA_g)
        # also depends on the other groups' sample statistics. Build the full
        # (cross-group) covariance of the stacked structural statistics so that
        # step 2 can compute a proper cross-group sandwich. Cov(vech(S_g)) =
        # gamma[[g]] / n_g, so the actual covariance of the stacked structural
        # statistics is jac.full %*% bdiag(gamma_g / n_g) %*% t(jac.full). This
        # rare path needs the explicit per-group observed Gamma.
        if (isTRUE(attr(jac, "caseB"))) {
          jac_full <- attr(jac, "jac.full")
          nobs_g <- unlist(fit@SampleStats@nobs)
          if (is.null(gamma)) {
            gamma <- fit@SampleStats@NACOV
            if (is.null(gamma) || is.null(gamma[[1]])) {
              gamma <- lavTech(fit, "gamma")
            }
          }
          gscaled <- lapply(seq_len(fit@Data@ngroups),
                            function(g) gamma[[g]] / nobs_g[g])
          gamma_obs_full <- lav_mat_bdiag(gscaled)
          step1$Gamma.eta.full <- jac_full %*% gamma_obs_full %*% t(jac_full)
          step1$caseB <- TRUE
        }
      } # no lv-interaction
      step1$Gamma.eta <- gamma_eta
      TRUE
    }, error = function(e) e)
    if (inherits(ge_try, "error")) {
      if (gamma_eta_required) {
        # the local / local.nt SEs cannot be computed without Gamma.eta
        lav_msg_stop(gettextf(
          "Gamma.eta (needed for se = %s) could not be computed: %s",
          dQuote(se, q = FALSE), conditionMessage(ge_try)))
      }
      # twostep / twostep.robust / naive: the SEs do not need Gamma.eta, so a
      # failure is non-fatal -- we simply do not report the corrected structural
      # test (the naive structural test is reported instead).
      step1$Gamma.eta <- NULL
    }
  }

  if (output == "list.step1.only") {
    # stop here, return interim results
    return(step1)
  }

  # If the structural model uses only a SUBSET of the variables in VETA, restrict
  # VETA / EETA / Gamma.eta to that sub-block. A bi-factor's specific factors, for
  # instance, are pure measurement factors with no structural role: they are part
  # of VETA (and hence of Gamma.eta, the NACOV of vech(VETA)) but absent from the
  # structural model, so the structural fit's Delta/WLS.V would not conform with
  # the full Gamma.eta. For ordinary models every latent appears in the structural
  # model, so 'keep' is everything and this is a no-op.
  if (sam_method %in% c("local", "fsr", "cfsr") &&
      !is.null(step1$VETA) && !is.null(colnames(step1$VETA[[1]])) &&
      !isTRUE(step1$caseB)) {
    keep <- tryCatch({
      struc_pt <- lav_pt_subset_sm(step1$PT, add_exo_cov = TRUE,
                    fixed_x = FALSE,
                    conditional_x = fit@Options$conditional.x,
                    free_fixed_var = TRUE,
                    meanstructure = fit@Options$meanstructure)
      struc_vars <- unique(c(struc_pt$lhs, struc_pt$rhs))
      which(colnames(step1$VETA[[1]]) %in% struc_vars)
    }, error = function(e) integer(0L))
    p_full <- ncol(step1$VETA[[1]])
    if (length(keep) > 0L && length(keep) < p_full) {
      meanstr <- length(step1$EETA[[1]]) > 0L
      # vech positions of the kept sub-block within the full vech(VETA)
      pos <- matrix(0L, p_full, p_full)
      pos[lower.tri(pos, diag = TRUE)] <- seq_len(p_full * (p_full + 1L) / 2L)
      pos[upper.tri(pos)] <- t(pos)[upper.tri(pos)]
      sub_pos <- pos[keep, keep, drop = FALSE]
      vech_keep <- sub_pos[lower.tri(sub_pos, diag = TRUE)]
      # Gamma.eta row order per group is c(EETA (means), vech(VETA));
      # under conditional.x it is c(vecr(cbind(EETA, RS_eta)), vech(VETA)),
      # i.e. (q+1) interleaved entries per variable, then vech(VETA)
      res_slopes_attr <- attr(step1$VETA, "res.slopes")
      if (!is.null(res_slopes_attr)) {
        q_x <- ncol(res_slopes_attr[[1]])
        ge_mean <- unlist(lapply(keep, function(j) {
          (j - 1L) * (q_x + 1L) + seq_len(q_x + 1L)
        }))
        ge_keep <- c(ge_mean, p_full * (q_x + 1L) + vech_keep)
      } else {
        ge_keep <- if (meanstr) c(keep, p_full + vech_keep) else vech_keep
      }
      for (g in seq_along(step1$VETA)) {
        step1$VETA[[g]] <- step1$VETA[[g]][keep, keep, drop = FALSE]
        if (meanstr) {
          step1$EETA[[g]] <- step1$EETA[[g]][keep]
        }
        # conditional.x: keep only the slopes of the retained variables
        # (the exogenous covariates are unaffected)
        if (!is.null(res_slopes_attr)) {
          res_slopes_attr[[g]] <- res_slopes_attr[[g]][keep, , drop = FALSE]
        }
        if (!is.null(step1$Gamma.eta) && !is.null(step1$Gamma.eta[[g]])) {
          step1$Gamma.eta[[g]] <-
            step1$Gamma.eta[[g]][ge_keep, ge_keep, drop = FALSE]
        }
      }
      if (!is.null(res_slopes_attr)) {
        attr(step1$VETA, "res.slopes") <- res_slopes_attr
      }
    }
  }

  # conditional.x: permute Gamma.eta to the structural model's internal order.
  # VETA/EETA are matched by name in the step-2 lavaan() call, but the NACOV
  # is consumed as-is, and the structural model's internal variable order
  # (y-variables first) generally differs from the measurement (VETA) order.
  # The permutation mirrors the lav_pt_subset_sm() call of lav_sam_step2().
  if (fit@Model@conditional.x &&
      sam_method %in% c("local", "fsr", "cfsr") &&
      !is.null(step1$Gamma.eta) && !is.null(step1$Gamma.eta[[1]])) {
    perm_ok <- tryCatch({
      struc_pt <- lav_pt_subset_sm(step1$PT, add_exo_cov = TRUE,
                    fixed_x = fit@Options$fixed.x,
                    conditional_x = TRUE,
                    free_fixed_var = TRUE,
                    meanstructure = fit@Options$meanstructure)
      rs <- attr(step1$VETA, "res.slopes")
      perm <- lav_sam_condx_perm_idx(
        eta_names  = colnames(step1$VETA[[1]]),
        x_names    = colnames(rs[[1]]),
        target_eta = unique(unlist(lav_pt_vnames(struc_pt, type = "ov.nox"))),
        target_x   = unique(unlist(lav_pt_vnames(struc_pt, type = "ov.x"))))
      for (g in seq_along(step1$Gamma.eta)) {
        step1$Gamma.eta[[g]] <- step1$Gamma.eta[[g]][perm, perm, drop = FALSE]
      }
      TRUE
    }, error = function(e) FALSE)
    if (!isTRUE(perm_ok)) {
      if (se %in% c("local", "local.nt")) {
        lav_msg_stop(gettextf(
          "Gamma.eta (needed for se = %s) could not be aligned with the
          structural model under conditional.x = TRUE.",
          dQuote(se, q = FALSE)))
      }
      step1$Gamma.eta <- NULL # no corrected structural test
    }
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
    # keep @Options$test consistent with the (structural) tests that ended up
    # in @test, just like the global branch below does (the joint itself was
    # fitted with test = "none")
    joint@Options$test <- unname(vapply(joint@test, `[[`,
                                        character(1L), "test"))
  } else {
    # global SAM: if requested (the default for sam.method = "global"), augment
    # the standard full-model test with the Yuan & Chan (2002) rescaled GLOBAL
    # chi-square, which corrects the full-model fit statistic for the separate
    # estimation of the measurement parameters in step 1 and for non-normality.
    # @test then holds both "standard" and "yuan.chan", exactly as sem() holds
    # "standard" and "satorra.bentler". See lav_sam_global_test().
    if (any(fit@Options$test == "yuan.chan")) {
      yc_test <- lav_sam_global_test(
        joint = joint, step1 = step1, step2 = step2, fit = fit
      )
      joint@test <- yc_test$test
      if (!is.null(yc_test$baseline.test)) {
        joint@baseline$test <- yc_test$baseline.test
      }
    }
    # the joint was fitted internally with "yuan.chan" stripped from its test
    # option (lavaan's core cannot compute it). Report exactly the test(s) that
    # ended up in @test, so @Options$test stays consistent with @test (as sem()
    # does for "satorra.bentler"). If the Yuan-Chan test could not be computed
    # (lav_sam_global_test() warned and fell back), this honestly reports just
    # "standard" instead of misleadingly claiming "yuan.chan".
    joint@Options$test <- unname(vapply(joint@test, `[[`, character(1L), "test"))
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
    joint@Options$bootstrap <- list(R = boot_out$R)
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
      # sanity check: the step-2 vcov must match the number of structural
      # parameters; a mismatch means the structural vcov could not be
      # computed properly (typically an identification problem that
      # slipped through, eg an empirically under-identified structural
      # part) -- give a clear error instead of a cryptic assignment
      # failure
      n2 <- length(step2$step2.free.idx)
      if (is.null(vcov_1$VCOV) || is.null(dim(vcov_1$VCOV)) ||
          !all(dim(vcov_1$VCOV) == c(n2, n2))) {
        lav_msg_stop(gettext(
          "the variance matrix of the structural (step 2) parameters could
           not be computed; the structural part may not be identified given
           the (fixed) measurement part. Consider simplifying the
           structural part, or using sem() instead."))
      }
      joint@vcov$vcov[step2$step2.free.idx, step2$step2.free.idx] <- vcov_1$VCOV
    }

    # Monte Carlo samples for defined parameter CIs (Preacher & Selig 2012)
    mc_coef <- NULL
    if (sum(pt_1$op == ":=") > 0L && lavoptions$se != "bootstrap" &&
      lav_model_vcov_se_mc_active(lavoptions)) {
      mc_coef <- lav_model_vcov_mc(lavmodel = joint@Model,
                                   vcov = joint@vcov$vcov,
                                   lavoptions = lavoptions)
      joint@internal$monte.carlo <- list(coef = mc_coef)
    }

    pt_1$se <- lav_model_vcov_se(
      lavmodel = joint@Model,
      lavpartable = pt_1,
      vcov = joint@vcov$vcov,
      mc = mc_coef,
      lavoptions = lavoptions
    )

    # surface SEs for latent (residual) variances that are 'fixed' in the
    # joint model but data-derived (eg std.lv = TRUE: the residual variance of
    # an endogenous factor = total - explained). These are free = 0 in the
    # joint model, so lav_model_vcov_se() set their se to 0; recover the proper
    # (delta-method) se from the joint vcov. See lav_sam_step2_se_lv_var().
    lv_var_se <- lav_sam_step2_se_lv_var(joint = joint)
    nz_idx <- which(lv_var_se > 0)
    if (length(nz_idx) > 0L) {
      pt_1$se[nz_idx] <- lv_var_se[nz_idx]
    }

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
    # replace the internal 'joint' lavaan call by the original sam() call,
    # so that update() and getCall() work (issue #514)
    res@call <- mc
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
