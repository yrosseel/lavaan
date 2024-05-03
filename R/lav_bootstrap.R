# main function used by various bootstrap related functions
# this function draws the bootstrap samples, and estimates the
# free parameters for each bootstrap sample
#
# return COEF matrix of size R x npar (R = number of bootstrap samples)
#
# Ed. 9 mar 2012
#
# Notes: - faulty runs are simply ignored (with a warning)
#        - default R=1000
#
# Updates: - now we have a separate @Data slot, we only need to transform once
#            for the bollen.stine bootstrap (13 dec 2011)
#          - bug fix: we need to 'update' the fixed.x variances/covariances
#            for each bootstrap draw!!
#
# Question: if fixed.x=TRUE, should we not keep X fixed, and bootstrap Y
#           only, conditional on X?? How to implement the conditional part?

# YR 27 Aug: - add keep.idx argument
#            - always return 'full' set of bootstrap results, including
#              failed runs (as NAs)
#            - idx nonadmissible/error solutions as an attribute
#            - thanks to keep.idx, it is easy to replicate/investigate these
#              cases if needed


bootstrapLavaan <- function(object,
                            R = 1000L,
                            type = "ordinary",
                            verbose = FALSE,
                            FUN = "coef",
                            # return.boot = FALSE, # no use, as boot stores
                            #                      # sample indices differently
                            keep.idx = FALSE,
                            parallel = c("no", "multicore", "snow"),
                            ncpus = 1L,
                            cl = NULL,
                            iseed = NULL,
                            h0.rmsea = NULL,
                            ...) {
  # checks
  type. <- tolower(type) # overwritten if nonparametric
  stopifnot(
    inherits(object, "lavaan"),
    type. %in% c(
      "nonparametric", "ordinary",
      "bollen.stine", "parametric", "yuan"
    )
  )
  if (type. == "nonparametric") {
    type. <- "ordinary"
  }
  if (missing(parallel)) {
    parallel <- "no"
  }
  parallel <- match.arg(parallel)

  # check if options$se is not bootstrap, otherwise, we get an infinite loop
  if (object@Options$se == "bootstrap") {
    object@Options$se <- "standard"
  }

  # check if options$test is not bollen.stine
  if ("bollen.stine" %in% object@Options$test) {
    object@Options$test <- "standard"
  }

  # check for conditional.x = TRUE
  if (object@Model@conditional.x) {
    lav_msg_stop(gettext(
      "this function is not (yet) available if conditional.x = TRUE"))
  }

  lavoptions. <- list(
    parallel = parallel, ncpus = ncpus, cl = cl,
    iseed = iseed
  )

  out <- lav_bootstrap_internal(
    object = object,
    lavdata. = NULL,
    lavmodel. = NULL,
    lavsamplestats. = NULL,
    lavoptions. = lavoptions.,
    lavpartable. = NULL,
    R = R,
    type = type.,
    verbose = verbose,
    FUN = FUN,
    keep.idx = keep.idx,
    h0.rmsea = h0.rmsea,
    ...
  )

  # new in 0.6-12: always warn for failed and nonadmissible runs
  nfailed <- length(attr(out, "error.idx")) # zero if NULL
  if (nfailed > 0L && object@Options$warn) {
    lav_msg_warn(gettextf(
      "%s bootstrap runs failed or did not converge.", nfailed))
  }

  notok <- length(attr(out, "nonadmissible")) # zero if NULL
  if (notok > 0L && object@Options$warn) {
     lav_msg_warn(gettextf(
      "%s bootstrap runs resulted in nonadmissible solutions.", notok))
  }

  out
}

# we need an internal version to be called from VCOV and lav_model_test
# when there is no lavaan object yet!
lav_bootstrap_internal <- function(object = NULL,
                                   lavdata. = NULL,
                                   lavmodel. = NULL,
                                   lavsamplestats. = NULL,
                                   lavoptions. = NULL,
                                   lavpartable. = NULL,
                                   R = 1000L,
                                   type = "ordinary",
                                   verbose = FALSE,
                                   FUN = "coef",
                                   # warn            = -1L, # not used anymore!
                                   check.post = TRUE,
                                   keep.idx = FALSE,
                                   # return.boot     = FALSE,
                                   h0.rmsea = NULL,
                                   ...) {
  # warning: avoid use of 'options', 'sample' (both are used as functions
  # below...
  # options -> opt
  # sample -> samp

  mc <- match.call()

  # object slots
  FUN.orig <- FUN
  if (!is.null(object)) {
    lavdata <- object@Data
    lavmodel <- object@Model
    lavsamplestats <- object@SampleStats
    lavoptions <- object@Options
    if (!is.null(lavoptions.)) {
      lavoptions$parallel <- lavoptions.$parallel
      lavoptions$ncpus <- lavoptions.$ncpus
      lavoptions$cl <- lavoptions.$cl
      lavoptions$iseed <- lavoptions.$iseed
    }
    lavpartable <- object@ParTable
    FUN <- match.fun(FUN)
    t0 <- FUN(object, ...)
    t.star <- matrix(as.numeric(NA), R, length(t0))
    colnames(t.star) <- names(t0)
  } else {
    # internal version!
    lavdata <- lavdata.
    lavmodel <- lavmodel.
    lavsamplestats <- lavsamplestats.
    lavoptions <- lavoptions.
    lavpartable <- lavpartable.
    if (FUN == "coef") {
      t.star <- matrix(as.numeric(NA), R, lavmodel@nx.free)
      lavoptions$test <- "none"
    } else if (FUN == "test") {
      t.star <- matrix(as.numeric(NA), R, 1L)
      lavoptions$test <- "standard"
    } else if (FUN == "coeftest") {
      t.star <- matrix(as.numeric(NA), R, lavmodel@nx.free + 1L)
      lavoptions$test <- "standard"
    }
  }

  # always shut off some options:
  lavoptions$verbose <- FALSE
  lavoptions$check.start <- FALSE
  lavoptions$check.post <- FALSE
  lavoptions$optim.attempts <- 1L # can save a lot of time

  # if internal or FUN == "coef", we can shut off even more
  if (is.null(object) || (is.character(FUN.orig) && FUN.orig == "coef")) {
    lavoptions$baseline <- FALSE
    lavoptions$h1 <- FALSE
    lavoptions$loglik <- FALSE
    lavoptions$implied <- FALSE
    lavoptions$store.vcov <- FALSE
    lavoptions$se <- "none"
    if (FUN.orig == "coef") {
      lavoptions$test <- "none"
    }
  }

  # bollen.stine, yuan, or parametric: we need the Sigma.hat values
  if (type == "bollen.stine" || type == "parametric" || type == "yuan") {
    Sigma.hat <- computeSigmaHat(lavmodel = lavmodel)
    Mu.hat <- computeMuHat(lavmodel = lavmodel)
  }

  # can we use the original data, or do we need to transform it first?
  if (type == "bollen.stine" || type == "yuan") {
    # check if data is continuous
    if (lavmodel@categorical) {
      lav_msg_stop(gettext(
        "bollen.stine/yuan bootstrap not available for categorical/ordinal data"
        ))
    }
    # check if data is complete
    if (lavoptions$missing != "listwise") {
      lav_msg_stop(gettext(
        "bollen.stine/yuan bootstrap not available for missing data"))
    }
    dataX <- vector("list", length = lavdata@ngroups)
  } else {
    dataX <- lavdata@X
  }

  # if bollen.stine, transform data here
  if (type == "bollen.stine") {
    for (g in 1:lavsamplestats@ngroups) {
      sigma.sqrt <- lav_matrix_symmetric_sqrt(Sigma.hat[[g]])
      S.inv.sqrt <- lav_matrix_symmetric_sqrt(lavsamplestats@icov[[g]])

      # center (needed???)
      X <- scale(lavdata@X[[g]], center = TRUE, scale = FALSE)

      # transform
      X <- X %*% S.inv.sqrt %*% sigma.sqrt

      # add model-based mean
      if (lavmodel@meanstructure) {
        X <- scale(X, center = (-1 * Mu.hat[[g]]), scale = FALSE)
      }

      # transformed data
      dataX[[g]] <- X
    }

    # if yuan, transform data here
  } else if (type == "yuan") {
    # page numbers refer to Yuan et al, 2007
    # Define a function to find appropriate value of a
    # (p. 272); code supplied 16 jun 2016 by Cheng & Wu
    search.a <- function(F0, d, p) {
      if (F0 == 0) {
        a0 <- 0
        return(a0)
      }
      max.a <- 1 / (1 - min(d)) - 1e-3
      # starting value; Yuan p. 272
      a0 <- min(sqrt(2 * F0 / sum((d - 1)^2)), max.a)

      # See Yuan p. 280
      for (i in 1:50) {
        dia <- a0 * d + (1 - a0)
        g1 <- -sum(log(dia)) + sum(dia) - p
        dif <- g1 - F0
        if (abs(dif) < 1e-6) {
          return(a0)
        }
        g2 <- a0 * sum((d - 1)^2 / dia)
        a0 <- min(max(a0 - dif / g2, 0), max.a)
      }
      # if search fails to converge in 50 iterations
      lav_msg_warn(gettext("yuan bootstrap search for `a` did not converge.
                           h0.rmsea may be too large."))
      a0
    }

    # Now use g.a within each group
    for (g in 1:lavsamplestats@ngroups) {
      S <- lavsamplestats@cov[[g]]
      # test is in Fit slot
      ghat <- object@test[[1]]$stat.group[[g]]
      df <- object@test[[1]]$df
      Sigmahat <- Sigma.hat[[g]]
      nmv <- nrow(Sigmahat)
      n <- nrow(lavdata@X[[g]])

      # Calculate tauhat_1, middle p. 267.
      # Yuan et al note that tauhat_1 could be negative;
      # if so, we need to let S.a = Sigmahat. (see middle p 275)
      ifelse(length(h0.rmsea) == 0,
        tau.hat <- (ghat - df) / (n - 1), # middle p 267
        tau.hat <- df * (h0.rmsea * h0.rmsea)
      ) # middle p 273

      if (tau.hat >= 0) {
        # from Cheng and Wu
        EL <- t(chol(Sigmahat))
        ESE <- forwardsolve(EL, t(forwardsolve(EL, S)))
        d <- eigen(ESE, symmetric = TRUE, only.values = TRUE)$values
        if ("a" %in% names(list(...))) {
          a <- list(...)$a
        } else {
          # Find a to minimize g.a
          a <- search.a(tau.hat, d, nmv)
        }
        # Calculate S_a (p. 267)
        S.a <- a * S + (1 - a) * Sigmahat
      } else {
        S.a <- Sigmahat
      }

      # Transform the data (p. 263)
      S.a.sqrt <- lav_matrix_symmetric_sqrt(S.a)
      S.inv.sqrt <- lav_matrix_symmetric_sqrt(lavsamplestats@icov[[g]])

      X <- lavdata@X[[g]]
      X <- X %*% S.inv.sqrt %*% S.a.sqrt

      # transformed data
      dataX[[g]] <- X
    }
  }

  # run bootstraps
  fn <- function(b) {
    # create bootstrap sample, and generate new 'data' object
    if (type == "bollen.stine" || type == "ordinary" || type == "yuan") {
      # take a bootstrap sample for each group
      BOOT.idx <- vector("list", length = lavdata@ngroups)
      # Note: we generate the bootstrap indices separately for each
      #       group, in order to ensure the group sizes do not change!
      for (g in 1:lavdata@ngroups) {
        stopifnot(nrow(lavdata@X[[g]]) > 1L)
        boot.idx <- sample.int(nrow(lavdata@X[[g]]), replace = TRUE)
        BOOT.idx[[g]] <- boot.idx
        dataX[[g]] <- dataX[[g]][boot.idx, , drop = FALSE]
      }
      newData <- lav_data_update(
        lavdata = lavdata, newX = dataX,
        BOOT.idx = BOOT.idx,
        lavoptions = lavoptions
      )
    } else { # parametric!
      for (g in 1:lavdata@ngroups) {
        dataX[[g]] <- MASS::mvrnorm(
          n = lavdata@nobs[[g]],
          Sigma = Sigma.hat[[g]],
          mu = Mu.hat[[g]]
        )
      }
      newData <- lav_data_update(
        lavdata = lavdata, newX = dataX,
        lavoptions = lavoptions
      )
    }

    # verbose
    if (verbose) cat("  ... bootstrap draw number:", sprintf("%4d", b))
    bootSampleStats <- try(lav_samplestats_from_data(
      lavdata       = newData,
      lavoptions    = lavoptions
    ), silent = TRUE)
    if (inherits(bootSampleStats, "try-error")) {
      if (verbose) {
        cat("     FAILED: creating sample statistics\n")
        cat(bootSampleStats[1])
      }
      out <- as.numeric(NA)
      attr(out, "nonadmissible.flag") <- TRUE
      if (keep.idx) {
        attr(out, "BOOT.idx") <- BOOT.idx
      }
      return(out)
    }

    # do we need to update Model slot? only if we have fixed exogenous
    # covariates, as their variances/covariances are stored in GLIST
    if (lavmodel@fixed.x && length(vnames(lavpartable, "ov.x")) > 0L) {
      model.boot <- NULL
    } else {
      model.boot <- lavmodel
    }

    # override option

    # fit model on bootstrap sample
    fit.boot <- suppressWarnings(lavaan(
      slotOptions = lavoptions,
      slotParTable = lavpartable,
      slotModel = model.boot,
      slotSampleStats = bootSampleStats,
      slotData = lavdata
    ))
    if (!fit.boot@optim$converged) {
      if (verbose) cat("     FAILED: no convergence\n")
      out <- as.numeric(NA)
      attr(out, "nonadmissible.flag") <- TRUE
      if (keep.idx) {
        attr(out, "BOOT.idx") <- BOOT.idx
      }
      return(out)
    }

    # extract information we need
    if (is.null(object)) { # internal use only!
      if (FUN == "coef") {
        out <- fit.boot@optim$x
      } else if (FUN == "test") {
        out <- fit.boot@test[[1L]]$stat
      } else if (FUN == "coeftest") {
        out <- c(fit.boot@optim$x, fit.boot@test[[1L]]$stat)
      }
    } else { # general use
      out <- try(as.numeric(FUN(fit.boot, ...)), silent = TRUE)
    }
    if (inherits(out, "try-error")) {
      if (verbose) cat("     FAILED: applying FUN to fit.boot\n")
      out <- as.numeric(NA)
      attr(out, "nonadmissible.flag") <- TRUE
      if (keep.idx) {
        attr(out, "BOOT.idx") <- BOOT.idx
      }
      return(out)
    }

    # check if the solution is admissible
    admissible.flag <- suppressWarnings(lavInspect(fit.boot, "post.check"))
    attr(out, "nonadmissible.flag") <- !admissible.flag

    if (verbose) {
      cat(
        "   OK -- niter = ",
        sprintf("%3d", fit.boot@optim$iterations), " fx = ",
        sprintf("%11.9f", fit.boot@optim$fx),
        if (admissible.flag) " " else "n", "\n"
      )
    }

    if (keep.idx) {
      # add BOOT.idx (for all groups)
      attr(out, "BOOT.idx") <- BOOT.idx
    }

    out
  } # end-of-fn

  # get parallelization options
  parallel <- lavoptions$parallel[1]
  ncpus <- lavoptions$ncpus
  cl <- lavoptions[["cl"]] # often NULL
  iseed <- lavoptions[["iseed"]] # often NULL

  # the next 10 lines are borrowed from the boot package
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") {
      have_mc <- .Platform$OS.type != "windows"
    } else if (parallel == "snow") have_snow <- TRUE
    if (!have_mc && !have_snow) ncpus <- 1L
    loadNamespace("parallel") # before recording seed!
  }

  # iseed:
  # this follows a proposal of Shu Fai Cheung (see github issue #240)
  # - iseed is used for both serial and parallel
  # - if iseed is not set, iseed is generated + .Random.seed created/updated
  #     -> tmp.seed <- NA
  # - if iseed is set: don't touch .Random.seed (if it exists)
  #     -> tmp.seed <- .Random.seed (if it exists)
  #     -> tmp.seed <- NULL (if it does not exist)
  if (is.null(iseed)) {
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      runif(1)
    }
    # identical(temp.seed, NA): Will not change .Random.seed in GlobalEnv
    temp.seed <- NA
    iseed <- runif(1, 0, 999999999)
  } else {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      temp.seed <-
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      # is.null(temp.seed): Will remove .Random.seed in GlobalEnv
      #                     if serial.
      #                     If parallel, .Random.seed will not be touched.
      temp.seed <- NULL
    }
  }
  if (!(ncpus > 1L && (have_mc || have_snow))) { # Only for serial
    set.seed(iseed)
  }


  # this is adapted from the boot function in package boot
  RR <- R
  if (verbose) {
    cat("\n")
  }
  res <- if (ncpus > 1L && (have_mc || have_snow)) {
    if (have_mc) {
      RNGkind_old <- RNGkind() # store current kind
      RNGkind("L'Ecuyer-CMRG") # to allow for reproducible results
      set.seed(iseed)
      parallel::mclapply(seq_len(RR), fn, mc.cores = ncpus)
    } else if (have_snow) {
      list(...) # evaluate any promises
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
        # # No need for
        # if(RNGkind()[1L] == "L'Ecuyer-CMRG")
        # clusterSetRNGStream() always calls `RNGkind("L'Ecuyer-CMRG")`
        parallel::clusterSetRNGStream(cl, iseed = iseed)
        res <- parallel::parLapply(cl, seq_len(RR), fn)
        parallel::stopCluster(cl)
        res
      } else {
        parallel::parLapply(cl, seq_len(RR), fn)
      }
    }
  } else {
    lapply(seq_len(RR), fn)
  }

  # restore old RNGkind()
  if (ncpus > 1L && have_mc) {
    RNGkind(RNGkind_old[1], RNGkind_old[2], RNGkind_old[3])
  }

  # fill in container
  t.star[] <- do.call("rbind", res)

  # handle errors
  error.idx <- which(sapply(res, function(x) is.na(x[1L])))
  attr(t.star, "error.idx") <- error.idx # could be integer(0L)

  # handle nonadmissible solutions
  if (check.post) {
    notok <- which(sapply(res, attr, "nonadmissible.flag"))
    if (length(error.idx) > 0L) {
      notok <- notok[-which(notok %in% error.idx)]
    }
    attr(t.star, "nonadmissible") <- notok
  }

  # store iseed
  attr(t.star, "seed") <- iseed

  # handle temp.seed
  if (!is.null(temp.seed) && !identical(temp.seed, NA)) {
    assign(".Random.seed", temp.seed, envir = .GlobalEnv)
  } else if (is.null(temp.seed) && !(ncpus > 1L && (have_mc || have_snow))) {
    # serial
    rm(.Random.seed, pos = 1)
  } else if (is.null(temp.seed) && (ncpus > 1L && have_mc)) {
    # parallel/multicore only
    rm(.Random.seed, pos = 1) # because set used set.seed()
  }

  # store BOOT.idx per group
  if (keep.idx) {
    BOOT.idx <- vector("list", length = lavsamplestats@ngroups)
    for (g in 1:lavsamplestats@ngroups) {
      # note that failed runs (NULL) are removed (for now)
      BOOT.idx[[g]] <- do.call(
        "rbind",
        lapply(res, function(x) attr(x, "BOOT.idx")[[g]])
      )
    }
    attr(t.star, "boot.idx") <- BOOT.idx
  }

  #    # No use, as boot package stores the sample indices differently
  #    # See boot:::boot.array() versus lav_utils_bootstrap_indices()
  #    if(return.boot) {
  #        # mimic output boot function
  #
  #        if(is.null(object)) {
  #            stop("lavaan ERROR: return.boot = TRUE requires a full lavaan object")
  #        }
  #
  #        # we start with ordinary only for now
  #        stopifnot(type == "ordinary")
  #
  #        if(! type %in% c("ordinary", "parametric")) {
  #            stop("lavaan ERROR: only ordinary and parametric bootstrap are supported if return.boot = TRUE")
  #        } else {
  #            sim <- type
  #        }
  #
  #        statistic. <- function(data, idx) {
  #                         data.boot <- data[idx,]
  #                         fit.boot <- update(object, data = data.boot)
  #                         out <- try(FUN(fit.boot, ...), silent = TRUE)
  #                         if(inherits(out, "try-error")) {
  #                             out <- rep(as.numeric(NA), length(t0))
  #                         }
  #                         out
  #                     }
  #        attr(t.star, "seed") <- NULL
  #        attr(t.star, "nonadmissible") <- NULL
  #        out <- list(t0 = t0, t = t.star, R = RR,
  #                    data = lavInspect(object, "data"),
  #                    seed = iseed, statistic = statistic.,
  #                    sim = sim, call = mc)
  #
  #        #if(sim == "parametric") {
  #        #    ran.gen. <- function() {} # TODO
  #        #    out <- c(out, list(ran.gen = ran.gen, mle = mle))
  #        #} else if(sim == "ordinary") {
  #            stype <- "i"
  #            strata <- rep(1, nobs(object))
  #            weights <- 1/tabulate(strata)[strata]
  #            out <- c(out, list(stype = stype, strata = strata,
  #                                weights = weights))
  #        #}
  #
  #        class(out) <- "boot"
  #        return(out)
  #    }
  #
  t.star
}
