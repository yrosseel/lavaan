# here, we generate new models based on the original model in lavobject
# 1. the independence model
# 2. the unrestricted model
# 3. model + extra parameters (for modindices/lavTestScore)
# 4. catML fit based on DWLS fit (for robust RMSEA/CFI)


# 1. fit an 'independence' model
#    note that for ML (and ULS and DWLS), the 'estimates' of the
#    independence model are simply the observed variances
#    but for GLS and WLS, this is not the case!!
lav_object_independence <- function(object = NULL,
                                    # or
                                    lavsamplestats = NULL,
                                    lavdata = NULL,
                                    lavcache = NULL,
                                    lavoptions = NULL,
                                    lavpartable = NULL,
                                    lavh1 = NULL,
                                    # local options
                                    se = FALSE,
                                    verbose = FALSE,
                                    warn = FALSE) {
  # object or slots?
  if (!is.null(object)) {
    stopifnot(inherits(object, "lavaan"))

    # extract needed slots
    lavsamplestats <- object@SampleStats
    lavdata <- object@Data
    lavcache <- object@Cache
    lavoptions <- object@Options
    lavpta <- object@pta
    if (.hasSlot(object, "h1")) {
      lavh1 <- object@h1
    } else {
      lavh1 <- lav_h1_implied_logl(
        lavdata = object@Data,
        lavsamplestats = object@SampleStats,
        lavoptions = object@Options
      )
    }
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
  lavpartable <- lav_partable_indep_or_unrestricted(
    lavobject = NULL,
    lavdata = lavdata, lavpta = lavpta, lavoptions = lavoptions,
    lavsamplestats = lavsamplestats, lavh1 = lavh1, independent = TRUE
  )

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
      dls.a <- lavoptions$estimator.args$dls.a
      for (g in 1:lavsamplestats@ngroups) {
        GammaNT <- lav_samplestats_Gamma_NT(
          COV            = lavsamplestats@cov[[g]],
          MEAN           = lavsamplestats@mean[[g]],
          rescale        = FALSE,
          x.idx          = lavsamplestats@x.idx[[g]],
          fixed.x        = lavoptions$fixed.x,
          conditional.x  = lavoptions$conditional.x,
          meanstructure  = lavoptions$meanstructure,
          slopestructure = lavoptions$conditional.x
        )
        W.DLS <- (1 - dls.a) * lavsamplestats@NACOV[[g]] + dls.a * GammaNT
        # overwrite
        lavsamplestats@WLS.V[[g]] <- lav_matrix_symmetric_inverse(W.DLS)
      }
    }
  }

  # se
  if (se) {
    if (lavoptions$se == "none") {
      lavoptions$se <- "standard"
    }
  } else {
    # 0.6-18: slower, but safer to just keep it

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
  lavoptions$rstarts <- 0L # no random starts

  # ALWAYS do.fit and set optim.method = "nlminb" (if npar > 0)
  npar <- lav_partable_npar(lavpartable)
  if (npar > 0L) {
    lavoptions$do.fit <- TRUE
    lavoptions$optim.method <- "nlminb"
  } else {
    # perhaps a correlation structure?
    lavoptions$optim.method <- "none"
    lavoptions$optim.force.converged <- TRUE
  }

  # verbose?
  lavoptions$verbose <- verbose

  # warn?
  lavoptions$warn <- warn

  # needed?
  if (any(lavpartable$op == "~1")) lavoptions$meanstructure <- TRUE

  # FIXME: it is crucial that the order of the ov's, as returned by
  # lavNames() remains the same
  # so lavNames(object) should equal lavNames(lavpartable)
  # otherwise, we will use the wrong sample statistics!!!
  #
  # this seems ok now, because we first generate the covariances in
  # lavpartable, and they should be in the right order (unlike the
  # intercepts)

  FIT <- lavaan(lavpartable,
    slotOptions     = lavoptions,
    slotSampleStats = lavsamplestats,
    slotData        = lavdata,
    slotCache       = lavcache,
    sloth1          = lavh1
  )

  FIT
}


# 2. unrestricted model
lav_object_unrestricted <- function(object, se = FALSE, verbose = FALSE,
                                    warn = FALSE) {
  # construct parameter table for unrestricted model
  lavpartable <- lav_partable_unrestricted(object)

  # adapt options
  lavoptions <- object@Options

  # se
  if (se) {
    if (lavoptions$se == "none") {
      lavoptions$se <- "standard"
    }
  } else {
    ## FIXME: if test = scaled, we need it anyway?
    lavoptions$se <- "none"
  }

  # ALWAYS do.fit
  lavoptions$do.fit <- TRUE

  # verbose?
  if (verbose) {
    lavoptions$verbose <- TRUE
  } else {
    lavoptions$verbose <- FALSE
  }

  # warn?
  if (warn) {
    lavoptions$warn <- TRUE
  } else {
    lavoptions$warn <- FALSE
  }

  # needed?
  if (any(lavpartable$op == "~1")) lavoptions$meanstructure <- TRUE

  if (.hasSlot(object, "h1")) {
    lavh1 <- object@h1
  } else {
    lavh1 <- lav_h1_implied_logl(
      lavdata = object@Data,
      lavsamplestats = object@SampleStats,
      lavoptions = object@Options
    )
  }

  FIT <- lavaan(lavpartable,
    slotOptions     = lavoptions,
    slotSampleStats = object@SampleStats,
    slotData        = object@Data,
    slotCache       = object@Cache,
    sloth1          = lavh1
  )

  FIT
}


# 3. extended model
lav_object_extended <- function(object, add = NULL,
                                remove.duplicated = TRUE,
                                all.free = FALSE,
                                verbose = FALSE, warn = FALSE,
                                do.fit = FALSE) {
  # partable original model
  partable <- object@ParTable[c(
    "lhs", "op", "rhs", "free", "exo", "label",
    "plabel"
  )]

  # new in 0.6-3: check for non-parameters
  nonpar.idx <- which(partable$op %in% c("==", ":=", "<", ">"))

  # always add block/group/level
  if (!is.null(object@ParTable$group)) {
    partable$group <- object@ParTable$group
  } else {
    partable$group <- rep(1L, length(partable$lhs))
    if (length(nonpar.idx) > 0L) {
      partable$group[nonpar.idx] <- 0L
    }
  }
  if (!is.null(object@ParTable$level)) {
    partable$level <- object@ParTable$level
  } else {
    partable$level <- rep(1L, length(partable$lhs))
    if (length(nonpar.idx) > 0L) {
      partable$level[nonpar.idx] <- 0L
    }
  }
  if (!is.null(object@ParTable$block)) {
    partable$block <- object@ParTable$block
  } else {
    partable$block <- rep(1L, length(partable$lhs))
    if (length(nonpar.idx) > 0L) {
      partable$block[nonpar.idx] <- 0L
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

  if (all.free) {
    partable$user <- rep(1L, length(partable$lhs))
    non.free.idx <- which(partable$free == 0L & partable$op != "==" &
      partable$op != ":=" & partable$op != "<" &
      partable$op != ">")
    partable$free[non.free.idx] <- 1L
    partable$user[non.free.idx] <- 10L
  }

  # replace 'start' column, since lav_model will fill these in in GLIST
  partable$start <- parameterEstimates(object,
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
    ADD <- add
  } else if (is.character(add)) {
    ngroups <- lav_partable_ngroups(partable)
    ADD.orig <- lavaanify(add, ngroups = ngroups)
    ADD <- ADD.orig[, c("lhs", "op", "rhs", "user", "label")] # minimum

    # always add block/group/level
    if (!is.null(ADD.orig$group)) {
      ADD$group <- ADD.orig$group
    } else {
      ADD$group <- rep(1L, length(ADD$lhs))
    }
    if (!is.null(ADD.orig$level)) {
      ADD$level <- ADD.orig$level
    } else {
      ADD$level <- rep(1L, length(ADD$lhs))
    }
    if (!is.null(ADD.orig$block)) {
      ADD$block <- ADD.orig$block
    } else {
      ADD$block <- rep(1L, length(ADD$lhs))
    }

    remove.idx <- which(ADD$user == 0)
    if (length(remove.idx) > 0L) {
      ADD <- ADD[-remove.idx, ]
    }
    ADD$start <- rep(0, nrow(ADD))
    ADD$free <- rep(1, nrow(ADD))
    ADD$user <- rep(10, nrow(ADD))
  }

  # merge
  LIST <- lav_partable_merge(partable, ADD,
    remove.duplicated = remove.duplicated,
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
  free.idx <- which(LIST$free > 0)
  LIST$free[free.idx] <- 1:length(free.idx)

  # adapt options
  lavoptions <- object@Options

  # verbose?
  lavoptions$verbose <- verbose

  # warn?
  lavoptions$warn <- warn

  # do.fit?
  lavoptions$do.fit <- do.fit

  # needed?
  if (any(LIST$op == "~1")) lavoptions$meanstructure <- TRUE

  if (.hasSlot(object, "h1")) {
    lavh1 <- object@h1
  } else {
    # old object -- for example 'usemmodelfit' in package 'pompom'

    # add a few fields
    lavoptions$h1 <- FALSE
    lavoptions$implied <- FALSE
    lavoptions$baseline <- FALSE
    lavoptions$loglik <- FALSE
    lavoptions$estimator.args <- list()

    # add a few slots
    object@Data@weights <- vector("list", object@Data@ngroups)
    object@Model@estimator <- object@Options$estimator
    object@Model@estimator.args <- list()

    lavh1 <- lav_h1_implied_logl(
      lavdata = object@Data,
      lavsamplestats = object@SampleStats,
      lavoptions = object@Options
    )
  }

  FIT <- lavaan(LIST,
    slotOptions     = lavoptions,
    slotSampleStats = object@SampleStats,
    slotData        = object@Data,
    slotCache       = object@Cache,
    sloth1          = lavh1
  )

  FIT
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
  if (all(lavdata@ov$type == "ordered")) {
    partable.catml <- parTable(lavobject)
    rm.idx <- which(partable.catml$op %in% c("|", "~1"))
    partable.catml <- partable.catml[-rm.idx, ]
    partable.catml <- lav_partable_complete(partable.catml)
  } else {
    refit <- TRUE
    partable.catml <- parTable(lavobject)
    partable.catml$start <- partable.catml$est
    partable.catml$se <- NULL
    rm.idx <- which(partable.catml$op %in% c("|", "~1"))
    partable.catml <- partable.catml[-rm.idx, ]
    partable.catml$ustart <- partable.catml$est
    for (b in seq_len(lavpta$nblocks)) {
      ov.names.num <- lavpta$vnames$ov.num[[b]]
      ov.var.idx <- which(partable.catml$op == "~~" &
        partable.catml$lhs %in% ov.names.num &
        partable.catml$lhs == partable.catml$rhs)
      partable.catml$free[ov.var.idx] <- 0L
    }
    partable.catml <- lav_partable_complete(partable.catml)
  }

  # adapt lavsamplestats
  for (g in seq_len(lavdata@ngroups)) {
    lavsamplestats@WLS.V[[g]] <- NULL
    lavsamplestats@WLS.VD[[g]] <- NULL
    COR <- lavsamplestats@cov[[g]]
    # check if COV is pd or not
    ev <- eigen(COR, symmetric = TRUE, only.values = TRUE)$values
    if (any(ev < .Machine$double.eps^(1 / 2))) {
      # not PD!
      COV <- cov2cor(lav_matrix_symmetric_force_pd(COR, tol = 1e-04))
      lavsamplestats@cov[[g]] <- COV
      lavsamplestats@var[[g]] <- diag(COV)
      refit <- TRUE
    } else {
      COV <- COR
    }

    out <- lav_samplestats_icov(
      COV = COV, ridge = 1e-05,
      x.idx = lavsamplestats@x.idx[[g]],
      ngroups = lavdata@ngroups, g = g,
      warn = FALSE
    )
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
  lavoptions$verbose <- FALSE
  lavoptions$debug <- FALSE
  lavoptions$.categorical <- FALSE
  lavoptions$categorical <- FALSE
  lavoptions$correlation <- TRUE
  lavoptions$meanstructure <- FALSE
  lavoptions$conditional.x <- FALSE # fixme
  lavoptions$information <- c("expected", "expected")
  lavoptions$h1.information <- c("structured", "structured") # unlike DWLS
  lavoptions$se <- "none"
  lavoptions$test <- "standard" # always for now
  lavoptions$baseline <- TRUE
  if (!refit) {
    lavoptions$optim.method <- "none"
    lavoptions$optim.force.converged <- TRUE
  } else {
    lavoptions$optim.gradient <- "numerical"
  }

  # dummy fit
  FIT <- lavaan(
    slotParTable = partable.catml,
    slotSampleStats = lavsamplestats,
    slotData = lavdata,
    slotOptions = lavoptions
  )

  FIT
}
