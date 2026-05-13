# here, we generate new models based on the original model in lavobject
# 1. the independence model
# 2. the unrestricted model
# 3. model + extra parameters (for modindices/lavTestScore)
# 4. catML fit based on DWLS fit (for robust RMSEA/CFI)


# 1. fit an 'independence/baseline' model
#    note that for ML (and ULS and DWLS), the 'estimates' of the
#    independence model are simply the observed variances
#    but for GLS and WLS, this is not the case!!
#
#  - YR 01 Feb 2026: allow for baseline.type = "nested"
lav_object_independence <- lav_object_baseline <- function(object = NULL,
                                                         # or
                                                         lavsamplestats = NULL,
                                                         lavdata = NULL,
                                                         lavcache = NULL,
                                                         lavoptions = NULL,
                                                         lavpartable = NULL,
                                                         lavh1 = NULL,
                                                         # local options
                                                         se = FALSE) {
  # object or slots?
  if (!is.null(object)) {
    stopifnot(inherits(object, "lavaan"))
    object <- lav_object_check_version(object)

    # extract needed slots
    lavsamplestats <- object@SampleStats
    lavdata <- object@Data
    lavcache <- object@Cache
    lavoptions <- object@Options
    lavpartable <- object@ParTable
    lavpta <- object@pta
    lavh1 <- object@h1
    if (is.null(lavoptions$estimator.args)) {
      lavoptions$estimator.args <- list()
    }
  } else {
    lavpta <- lav_partable_attributes(lavpartable)
    lavpartable <- lav_partable_set_cache(lavpartable, lavpta)
  }

  # if two-level, force conditional.x = FALSE (for now)
  if (lavdata@nlevels > 1L && lavoptions$conditional.x) {
    lavoptions$conditional.x <- FALSE
  }

  # construct parameter table for independence model
  if (!is.null(lavoptions$baseline.type) &&
    lavoptions$baseline.type == "nested") {
    lavpartable <- lav_partable_baseline(
      lavobject = NULL,
      lavpartable = lavpartable, lavh1 = lavh1
    )
  } else {
    lavpartable <- lav_partable_indep_or_unrestricted(
      lavobject = NULL,
      lavdata = lavdata, lavpta = lavpta, lavoptions = lavoptions,
      lavsamplestats = lavsamplestats, lavh1 = lavh1, independent = TRUE
    )
  }

  # new in 0.6-6: add lower bounds for ov.var
  if (!is.null(lavoptions$optim.bounds)) {
    lavoptions$bounds <- "doe.maar"
    lavoptions$effect.coding <- "" # to avoid warning
    lavoptions$optim.bounds <- list(lower = "ov.var")
    lavpartable <- lav_partable_add_bounds(
      partable = lavpartable,
      lavh1 = lavh1, lavdata = lavdata,
      lavsamplestats = lavsamplestats, lavoptions = lavoptions
    )
  }

  # new in 0.6-8: if DLS, change to sample-based
  if (lavoptions$estimator == "DLS") {
    if (lavoptions$estimator.args$dls.GammaNT == "sample") {
      # nothing to do
    } else {
      lavoptions$estimator.args$dls.GammaNT <- "sample"
      dls_a <- lavoptions$estimator.args$dls.a
      for (g in 1:lavsamplestats@ngroups) {
        gamma_nt <- lav_samplestats_gamma_nt(
          m_cov          = lavsamplestats@cov[[g]],
          m_mean         = lavsamplestats@mean[[g]],
          x_idx          = lavsamplestats@x.idx[[g]],
          fixed_x        = lavoptions$fixed.x,
          conditional_x  = lavoptions$conditional.x,
          meanstructure  = lavoptions$meanstructure,
          slopestructure = lavoptions$conditional.x
        )
        w_dls <- (1 - dls_a) * lavsamplestats@NACOV[[g]] + dls_a * gamma_nt
        # overwrite
        lavsamplestats@WLS.V[[g]] <- lav_matrix_symmetric_inverse(w_dls)
      }
    }
  }

  # se
  if (se) {
    if (lavoptions$se == "none") {
      lavoptions$se <- "standard"
    }
  } else {
    if (identical(lavoptions$test, "standard")) {
      lavoptions$se <- "none"
    }

    # 0.6-20 -- except if se = "bootstrap" -> "none"
    if (lavoptions$se == "bootstrap") {
      lavoptions$se <- "none"
    }

    ## FIXME: if test = scaled, we need it anyway?
    # if(lavoptions$missing %in% c("two.stage", "two.stage.robust")) {
    # don't touch it
    # } else {
    #    lavoptions$se <- "none"
    # }
  }

  # change options
  lavoptions$h1 <- FALSE # already provided by lavh1
  lavoptions$baseline <- FALSE # of course
  lavoptions$loglik <- TRUE # eg for multilevel
  lavoptions$implied <- TRUE # needed for loglik (multilevel)
  lavoptions$check.start <- FALSE
  lavoptions$check.gradient <- FALSE
  lavoptions$check.post <- FALSE
  lavoptions$check.vcov <- FALSE
  lavoptions$optim.bounds <- list() # we already have the bounds
  lavoptions$start <- "default" # don't re-use user-specified starting values
  lavoptions$rstarts <- 0L # no random starts

  # ALWAYS do.fit and set optim.method = "nlminb" (if npar > 0)
  npar <- lav_partable_npar(lavpartable)
  if (npar > 0L) {
    lavoptions$do.fit <- TRUE
    if (lavoptions$optim.method != "noniter") {
      lavoptions$optim.method <- "nlminb"
    }
  } else {
    # perhaps a correlation structure?
    lavoptions$optim.method <- "none"
    lavoptions$optim.force.converged <- TRUE
  }

  # needed?
  if (any(lavpartable$op == "~1")) lavoptions$meanstructure <- TRUE

  # FIXME: it is crucial that the order of the ov's, as returned by
  # lav_object_vnames() remains the same
  # so lav_object_vnames(object) should equal lav_object_vnames(lavpartable)
  # otherwise, we will use the wrong sample statistics!!!
  #
  # this seems ok now, because we first generate the covariances in
  # lavpartable, and they should be in the right order (unlike the
  # intercepts)

  fit <- lavaan(lavpartable,
    slotOptions     = lavoptions,
    slotSampleStats = lavsamplestats,
    slotData        = lavdata,
    slotCache       = lavcache,
    sloth1          = lavh1
  )

  fit
}


# 3. extended model (used for modification indices)
lav_object_extended <- function(object, add = NULL,
                                remove_duplicated = TRUE,
                                all_free = FALSE,
                                do_fit = FALSE) {
  object <- lav_object_check_version(object)
  # partable original model
  partable <- object@ParTable[c(
    "lhs", "op", "rhs", "free", "exo", "label",
    "plabel"
  )]

  # new in 0.6-3: check for non-parameters
  nonpar_idx <- which(partable$op %in% c("==", ":=", "<", ">"))

  # always add block/group/level
  if (!is.null(object@ParTable$group)) {
    partable$group <- object@ParTable$group
  } else {
    partable$group <- rep(1L, length(partable$lhs))
    if (length(nonpar_idx) > 0L) {
      partable$group[nonpar_idx] <- 0L
    }
  }
  if (!is.null(object@ParTable$level)) {
    partable$level <- object@ParTable$level
  } else {
    partable$level <- rep(1L, length(partable$lhs))
    if (length(nonpar_idx) > 0L) {
      partable$level[nonpar_idx] <- 0L
    }
  }
  if (!is.null(object@ParTable$block)) {
    partable$block <- object@ParTable$block
  } else {
    partable$block <- rep(1L, length(partable$lhs))
    if (length(nonpar_idx) > 0L) {
      partable$block[nonpar_idx] <- 0L
    }
  }

  # TDJ: Added to prevent error when lav_partable_merge() is called below.
  #      Problematic if object@ParTable is missing one of the requested slots,
  #      which returns a NULL slot with a missing <NA> name.  For example:
  #        example(cfa)
  #        lav_partable_independence(lavdata = fit@Data, lavpta = fit@pta,
  #                                  lavoptions = lavInspect(fit, "options"))
  #     Has no "label" or "plabel" elements.
  empties <- which(sapply(partable, is.null))
  if (length(empties)) {
    partable[empties] <- NULL
  }

  if (all_free) {
    partable$user <- rep(1L, length(partable$lhs))
    non_free_idx <- which(partable$free == 0L & partable$op != "==" &
      partable$op != ":=" & partable$op != "<" &
      partable$op != ">")
    partable$free[non_free_idx] <- 1L
    partable$user[non_free_idx] <- 10L
  }

  # replace 'start' column, since lav_model will fill these in in GLIST
  partable$start <- lavParameterEstimates(object,
    remove.system.eq = FALSE,
    remove.def = FALSE,
    remove.eq = FALSE,
    remove.ineq = FALSE,
    remove.nonfree = FALSE,
    remove.unused = FALSE
  )$est

  # add new parameters, extend model
  if (is.list(add)) {
    stopifnot(
      !is.null(add$lhs),
      !is.null(add$op),
      !is.null(add$rhs)
    )
    add_1 <- add
  } else if (is.character(add)) {
    ngroups <- lav_partable_ngroups(partable)
    add_orig <- lav_model_partable(add, ngroups = ngroups)
    add_1 <- add_orig[, c("lhs", "op", "rhs", "user", "label")] # minimum

    # always add block/group/level
    if (!is.null(add_orig$group)) {
      add_1$group <- add_orig$group
    } else {
      add_1$group <- rep(1L, length(add_1$lhs))
    }
    if (!is.null(add_orig$level)) {
      add_1$level <- add_orig$level
    } else {
      add_1$level <- rep(1L, length(add_1$lhs))
    }
    if (!is.null(add_orig$block)) {
      add_1$block <- add_orig$block
    } else {
      add_1$block <- rep(1L, length(add_1$lhs))
    }

    remove_idx <- which(add_1$user == 0)
    if (length(remove_idx) > 0L) {
      add_1 <- add_1[-remove_idx, ]
    }
    add_1$start <- rep(0, nrow(add_1))
    add_1$free <- rep(1, nrow(add_1))
    add_1$user <- rep(10, nrow(add_1))
  }

  # merge
  list_1 <- lav_partable_merge(partable, add_1,
    remove.duplicated = remove_duplicated,
    warn = FALSE
  )

  # remove nonpar?
  # if(remove.nonpar) {
  #    nonpar.idx <- which(LIST$op %in% c("==", ":=", "<", ">"))
  #    if(length(nonpar.idx) > 0L) {
  #        LIST <- LIST[-nonpar.idx,]
  #    }
  # }

  # redo 'free'
  free_idx <- which(list_1$free > 0)
  list_1$free[free_idx] <- seq_along(free_idx)

  # adapt options
  lavoptions <- object@Options

  # do.fit?
  lavoptions$do.fit <- do_fit

  # needed?
  if (any(list_1$op == "~1")) lavoptions$meanstructure <- TRUE

  # switch off 'consider switching to parameterization = theta' warning
  # (modindices!)
  lavoptions$check.delta.cat.mediator <- FALSE

  fit <- lavaan(list_1,
    slotOptions     = lavoptions,
    slotSampleStats = object@SampleStats,
    slotData        = object@Data,
    slotCache       = object@Cache,
    sloth1          = object@h1
  )

  fit
}

# 4. catml model
lav_object_catml <- function(lavobject = NULL) {
  stopifnot(inherits(lavobject, "lavaan"))
  stopifnot(lavobject@Model@categorical)

  # extract slots
  lavdata <- lavobject@Data
  lavsamplestats <- lavobject@SampleStats
  lavoptions <- lavobject@Options
  lavpta <- lavobject@pta

  # if only categorical variables: remove thresholds and intercepts
  refit <- FALSE
  partable_catml <- parTable(lavobject)
  rm_idx <- which(partable_catml$op %in% c("|", "~1"))
  # remove also any constraints that refer to the labels that belong
  # to these intercepts/thresholds (new in 0.6-20)
  th_int_labels <- unique(
    partable_catml$label[rm_idx],
    partable_catml$plabel[rm_idx]
  )
  if (any(nchar(th_int_labels) == 0L)) {
    th_int_labels <- th_int_labels[-which(nchar(th_int_labels) == 0L)]
  }
  con_idx <- which(partable_catml$op %in% c("==", "<", ">", ":=") &
    (partable_catml$lhs %in% th_int_labels |
      partable_catml$rhs %in% th_int_labels))
  partable_catml <- partable_catml[-c(rm_idx, con_idx), ]
  # never rotate (new in 0.6-19), as we only need fit measures
  if (!is.null(partable_catml$efa)) {
    partable_catml$efa <- NULL
    partable_catml$free <- partable_catml$free.unrotated
  }

  if (!all(lavdata@ov$type == "ordered")) {
    refit <- TRUE
    partable_catml$start <- partable_catml$est
    partable_catml$se <- NULL
    partable_catml$ustart <- partable_catml$est
    for (b in seq_len(lavpta$nblocks)) {
      ov_names_num <- lavpta$vnames$ov.num[[b]]
      ov_var_idx <- which(partable_catml$op == "~~" &
        partable_catml$lhs %in% ov_names_num &
        partable_catml$lhs == partable_catml$rhs)
      partable_catml$free[ov_var_idx] <- 0L
    }
  }
  partable_catml <- lav_partable_complete(partable_catml)

  # adapt lavsamplestats
  for (g in seq_len(lavdata@ngroups)) {
    lavsamplestats@WLS.V[[g]] <- NULL
    lavsamplestats@WLS.VD[[g]] <- NULL
    cor_1 <- lavsamplestats@cov[[g]]
    # check if COV is pd or not
    ev <- eigen(cor_1, symmetric = TRUE, only.values = TRUE)$values
    if (any(ev < .Machine$double.eps^(1 / 2))) {
      # not PD!
      cov_1 <- cov2cor(lav_matrix_symmetric_force_pd(cor_1, tol = 1e-04))
      lavsamplestats@cov[[g]] <- cov_1
      lavsamplestats@var[[g]] <- diag(cov_1)
      refit <- TRUE
    } else {
      cov_1 <- cor_1
    }

    current_warn <- lav_warn()
    if (lav_warn(FALSE)) {
      on.exit(lav_warn(current_warn), TRUE)
    }
    out <- lav_samplestats_icov(
      cov_1 = cov_1, ridge = 1e-05,
      x_idx = lavsamplestats@x.idx[[g]],
      ngroups = lavdata@ngroups, g = g
    )
    lav_warn(current_warn)
    lavsamplestats@icov[[g]] <- out$icov
    lavsamplestats@cov.log.det[[g]] <- out$cov.log.det

    # NACOV <- lavsamplestats@NACOV[[g]]
    # nvar   <- nrow(COV)
    # ntotal <- nrow(NACOV)
    # pstar <- nvar*(nvar-1)/2
    # nocor <- ntotal - pstar
    # if(length(nocor) > 0L) {
    #    lavsamplestats@NACOV[[g]] <- NACOV[-seq_len(nocor),
    #                                       -seq_len(nocor)]
    # }
  }

  # adapt lavoptions
  lavoptions$estimator <- "catML"
  lavoptions$.categorical <- FALSE
  lavoptions$categorical <- FALSE
  lavoptions$correlation <- TRUE
  lavoptions$meanstructure <- FALSE
  lavoptions$conditional.x <- FALSE # fixme
  lavoptions$information <- c("expected", "expected")
  lavoptions$h1.information <- c("structured", "structured") # unlike DWLS
  lavoptions$se <- "none"
  lavoptions$test <- "standard" # always for now
  lavoptions$rotation <- "none" # new in 0.6-19
  lavoptions$baseline <- TRUE # also for RMSEA?
  if (!refit) {
    lavoptions$optim.method <- "none"
    lavoptions$optim.force.converged <- TRUE
  } else {
    lavoptions$optim.gradient <- "numerical"
  }

  # dummy fit
  fit <- lavaan(
    slotParTable = partable_catml,
    slotSampleStats = lavsamplestats,
    slotData = lavdata,
    slotOptions = lavoptions
  )

  fit
}
