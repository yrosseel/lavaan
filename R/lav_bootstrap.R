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

# YR 27 Aug 2022: - add keep.idx argument
#                 - always return 'full' set of bootstrap results, including
#                   failed runs (as NAs)
#                 - idx nonadmissible/error solutions as an attribute
#                 - thanks to keep.idx, it is easy to replicate/investigate
#                   these cases if needed

# YR 10 Nov 2024: - detect sam object

lavBootstrap <- function(object,                                   # nolint start
                            R = 1000L,
                            type = "ordinary",
                            verbose = FALSE,
                            FUN = "coef",
                            # return.boot = FALSE, # no use, as boot stores
                            #                      # sample indices differently
                            keep.idx = FALSE,
                            parallel = c("no", "multicore", "snow"),
                            ncpus = max(1L, parallel::detectCores() - 2L),
                            cl = NULL,
                            iseed = NULL,
                            h0.rmsea = NULL,
                            ...) {                                 # nolint end

  # check object
  object <- lav_object_check_version(object)

  # checks
  type_1 <- tolower(type) # overwritten if nonparametric
  stopifnot(
    inherits(object, "lavaan"),
    type_1 %in% c(
      "nonparametric", "ordinary",
      "bollen.stine", "parametric", "yuan"
    )
  )
  if (!missing(verbose)) {
    current_verbose <- lav_verbose()
    if (lav_verbose(verbose))
      on.exit(lav_verbose(current_verbose), TRUE)
  }
  if (type_1 == "nonparametric") {
    type_1 <- "ordinary"
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

  lavoptions <- list(
    parallel = parallel, ncpus = ncpus, cl = cl,
    iseed = iseed
  )

  out <- lav_bootstrap_internal(
    object = object,
    lavdata = NULL,
    lavmodel = NULL,
    lavsamplestats = NULL,
    lavoptions = lavoptions,
    lavpartable = NULL,
    r = R,
    show_progress = verbose,
    type = type_1,
    fun = FUN,
    keep_idx = keep.idx,
    h0_rmsea = h0.rmsea,
    ...
  )

  # new in 0.6-12: always warn for failed and nonadmissible runs
  nfailed <- length(attr(out, "error.idx")) # zero if NULL
  if (nfailed > 0L) {
    lav_msg_warn(gettextf(
      "%s bootstrap runs failed or did not converge.", nfailed))
  }

  notok <- length(attr(out, "nonadmissible")) # zero if NULL
  if (notok > 0L) {
     lav_msg_warn(gettextf(
      "%s bootstrap runs resulted in nonadmissible solutions.", notok))
  }

  out
}
lavBootstrap <- lavBootstrap  # synonym                                  #nolint

# we need an internal version to be called from VCOV and lav_model_test
# when there is no lavaan object yet!
lav_bootstrap_internal <- function(object = NULL,
                                   lavdata = NULL,
                                   lavmodel = NULL,
                                   lavsamplestats = NULL,
                                   lavoptions = NULL,
                                   lavpartable = NULL,
                                   r = 1000L,
                                   show_progress = FALSE,
                                   type = "ordinary",
                                   fun = "coef",
                                   check_post = TRUE,
                                   keep_idx = FALSE,
                                   # return.boot     = FALSE,
                                   h0_rmsea = NULL,
                                   ...) {
  # warning: avoid use of 'options', 'sample' (both are used as functions
  # below...
  # options -> opt
  # sample -> samp

  # object slots
  fun_orig <- fun
  has_sam_object_flag <- FALSE
  if (!is.null(object)) {
    stopifnot(inherits(object, "lavaan"))
    # check for sam object
    if (!is.null(object@internal$sam.method)) {
      has_sam_object_flag <- TRUE
    }
    lavdata_1 <- object@Data
    lavmodel_1 <- object@Model
    lavsamplestats_1 <- object@SampleStats
    lavoptions_1 <- object@Options
    if (!is.null(lavoptions)) {
      lavoptions_1$parallel <- lavoptions$parallel
      lavoptions_1$ncpus <- lavoptions$ncpus
      lavoptions_1$cl <- lavoptions$cl
      lavoptions_1$iseed <- lavoptions$iseed
    }
    lavpartable_1 <- object@ParTable
    fun <- match.fun(fun)
    t0 <- fun(object, ...)
    t_star <- matrix(as.numeric(NA), r, length(t0))
    colnames(t_star) <- names(t0)
  } else {
    # internal version!
    lavdata_1 <- lavdata
    lavmodel_1 <- lavmodel
    lavsamplestats_1 <- lavsamplestats
    lavoptions_1 <- lavoptions
    lavpartable_1 <- lavpartable
    if (fun == "coef") {
      t_star <- matrix(as.numeric(NA), r, lavmodel_1@nx.free)
      lavoptions_1$test <- "none"
    } else if (fun == "test") {
      t_star <- matrix(as.numeric(NA), r, 1L)
      lavoptions_1$test <- "standard"
    } else if (fun == "coeftest") {
      t_star <- matrix(as.numeric(NA), r, lavmodel_1@nx.free + 1L)
      lavoptions_1$test <- "standard"
    }
  }

  # always shut off some options:
  current_verbose <- lav_verbose()
  if (missing(show_progress)) {
    show_progress <- current_verbose
  }
  if (lav_verbose(FALSE)) on.exit(lav_verbose(current_verbose))
  lavoptions_1$check.start <- FALSE
  lavoptions_1$check.post <- FALSE
  lavoptions_1$optim.attempts <- 1L # can save a lot of time

  # if internal or fun == "coef", we can shut off even more
  if (is.null(object) || (is.character(fun_orig) && fun_orig == "coef")) {
    lavoptions_1$baseline <- FALSE
    lavoptions_1$h1 <- FALSE
    lavoptions_1$loglik <- FALSE
    lavoptions_1$implied <- FALSE
    lavoptions_1$store.vcov <- FALSE
    lavoptions_1$se <- "none"
    if (fun_orig == "coef") {
      lavoptions_1$test <- "none"
    }
  }

  # bollen.stine, yuan, or parametric: we need the Sigma.hat values
  if (type == "bollen.stine" || type == "parametric" || type == "yuan") {
    sigma_hat <- lav_model_sigma(lavmodel = lavmodel_1)
    mu_hat <- lav_model_mu(lavmodel = lavmodel_1)
  }

  # can we use the original data, or do we need to transform it first?
  if (type == "bollen.stine" || type == "yuan") {
    # check if data is continuous
    if (lavmodel_1@categorical) {
      lav_msg_stop(gettext(
        "bollen.stine/yuan bootstrap not available for categorical/ordinal data"
        ))
    }
    # check if data is complete
    if (lavoptions_1$missing != "listwise") {
      lav_msg_stop(gettext(
        "bollen.stine/yuan bootstrap not available for missing data"))
    }
    data_x <- vector("list", length = lavdata_1@ngroups)
  } else {
    data_x <- lavdata_1@X
  }

  # if bollen.stine, transform data here
  if (type == "bollen.stine") {
    for (g in 1:lavsamplestats_1@ngroups) {
      sigma_sqrt <- lav_matrix_symmetric_sqrt(sigma_hat[[g]])
      s_inv_sqrt <- lav_matrix_symmetric_sqrt(lavsamplestats_1@icov[[g]])

      # center (needed???)
      x_1 <- scale(lavdata_1@X[[g]], center = TRUE, scale = FALSE)

      # transform
      x_1 <- x_1 %*% s_inv_sqrt %*% sigma_sqrt

      # add model-based mean
      if (lavmodel_1@meanstructure) {
        x_1 <- scale(x_1, center = (-1 * mu_hat[[g]]), scale = FALSE)
      }

      # transformed data
      data_x[[g]] <- x_1
    }

    # if yuan, transform data here
  } else if (type == "yuan") {
    # page numbers refer to Yuan et al, 2007
    # Define a function to find appropriate value of a
    # (p. 272); code supplied 16 jun 2016 by Cheng & Wu
    search_a <- function(f0, d, p) {
      if (f0 == 0) {
        a0 <- 0
        return(a0)
      }
      max_a <- 1 / (1 - min(d)) - 1e-3
      # starting value; Yuan p. 272
      a0 <- min(sqrt(2 * f0 / sum((d - 1)^2)), max_a)

      # See Yuan p. 280
      for (i in 1:50) {
        dia <- a0 * d + (1 - a0)
        g1 <- -sum(log(dia)) + sum(dia) - p
        dif <- g1 - f0
        if (abs(dif) < 1e-6) {
          return(a0)
        }
        g2 <- a0 * sum((d - 1)^2 / dia)
        a0 <- min(max(a0 - dif / g2, 0), max_a)
      }
      # if search fails to converge in 50 iterations
      lav_msg_warn(gettext("yuan bootstrap search for `a` did not converge.
                           h0.rmsea may be too large."))
      a0
    }

    # Now use g.a within each group
    for (g in 1:lavsamplestats_1@ngroups) {
      s <- lavsamplestats_1@cov[[g]]
      # test is in Fit slot
      ghat <- object@test[[1]]$stat.group[[g]]
      df <- object@test[[1]]$df
      sigmahat <- sigma_hat[[g]]
      nmv <- nrow(sigmahat)
      n <- nrow(lavdata_1@X[[g]])

      # Calculate tauhat_1, middle p. 267.
      # Yuan et al note that tauhat_1 could be negative;
      # if so, we need to let S.a = Sigmahat. (see middle p 275)
      ifelse(length(h0_rmsea) == 0,
        tau_hat <- (ghat - df) / (n - 1), # middle p 267
        tau_hat <- df * (h0_rmsea * h0_rmsea)
      ) # middle p 273

      if (tau_hat >= 0) {
        # from Cheng and Wu
        el_1 <- t(chol(sigmahat))
        ese <- forwardsolve(el_1, t(forwardsolve(el_1, s)))
        d <- eigen(ese, symmetric = TRUE, only.values = TRUE)$values
        if ("a" %in% names(list(...))) {
          a <- list(...)$a
        } else {
          # Find a to minimize g.a
          a <- search_a(tau_hat, d, nmv)
        }
        # Calculate S_a (p. 267)
        s_a <- a * s + (1 - a) * sigmahat
      } else {
        s_a <- sigmahat
      }

      # Transform the data (p. 263)
      s_a_sqrt <- lav_matrix_symmetric_sqrt(s_a)
      s_inv_sqrt <- lav_matrix_symmetric_sqrt(lavsamplestats_1@icov[[g]])

      x_1 <- lavdata_1@X[[g]]
      x_1 <- x_1 %*% s_inv_sqrt %*% s_a_sqrt

      # transformed data
      data_x[[g]] <- x_1
    }
  }

  # run bootstraps
  fn <- function(b) {
    # create bootstrap sample, and generate new 'data' object
    if (type == "bollen.stine" || type == "ordinary" || type == "yuan") {
      # take a bootstrap sample for each group
      boot_idx_1 <- vector("list", length = lavdata_1@ngroups)
      # Note: we generate the bootstrap indices separately for each
      #       group, in order to ensure the group sizes do not change!
      for (g in 1:lavdata_1@ngroups) {
        stopifnot(nrow(lavdata_1@X[[g]]) > 1L)
        boot_idx <- sample.int(nrow(lavdata_1@X[[g]]), replace = TRUE)
        boot_idx_1[[g]] <- boot_idx
        data_x[[g]] <- data_x[[g]][boot_idx, , drop = FALSE]
      }
      new_data <- lav_data_update(
        lavdata = lavdata_1, newX = data_x,
        BOOT.idx = boot_idx_1,
        lavoptions = lavoptions_1
      )
    } else { # parametric! (using sign-invariant method for reproducibility)
      for (g in 1:lavdata_1@ngroups) {
        data_x[[g]] <- lav_mvrnorm(
          n = lavdata_1@nobs[[g]],
          Sigma = sigma_hat[[g]],
          mu = mu_hat[[g]]
        )
      }
      new_data <- lav_data_update(
        lavdata = lavdata_1, newX = data_x,
        lavoptions = lavoptions_1
      )
    }

    # show progress?
    if (show_progress) {
      cat(if (interactive()) "\r" else "",
          "  ... bootstrap draw number:", sprintf("%4d", b))
    }
    boot_sample_stats <- try(lav_samplestats_from_data(
      lavdata       = new_data,
      lavoptions    = lavoptions_1
    ), silent = TRUE)
    if (inherits(boot_sample_stats, "try-error")) {
      if (show_progress) {
        cat("     FAILED: creating sample statistics\n")
        cat(boot_sample_stats[1])
      }
      out <- as.numeric(NA)
      attr(out, "nonadmissible.flag") <- TRUE
      if (keep_idx) {
        attr(out, "BOOT.idx") <- boot_idx_1
      }
      return(out)
    }
    if (has_sam_object_flag) {
      # also need h1
      booth1 <- lav_h1_implied_logl(lavdata = new_data,
        lavsamplestats = boot_sample_stats, lavpartable = lavpartable_1,
        lavoptions = lavoptions_1)
    }

    # do we need to update Model slot? only if we have fixed exogenous
    # covariates, as their variances/covariances are stored in GLIST
    if (lavmodel_1@fixed.x &&
        length(lav_partable_vnames(lavpartable_1, "ov.x")) > 0L) {
      model_boot <- NULL
    } else {
      model_boot <- lavmodel_1
    }

    # fit model on bootstrap sample
    if (has_sam_object_flag) {
      new_object <- object
      new_object@Data <- new_data
      new_object@SampleStats <- boot_sample_stats
      new_object@h1 <- booth1
      # what about lavoptions?
      fit_boot <- suppressWarnings(try(sam(new_object, se = "none"),
                                       silent = FALSE)) # show what is wrong
    } else {
      fit_boot <- suppressWarnings(try(lavaan(
        slotOptions = lavoptions_1,
        slotParTable = lavpartable_1,
        slotModel = model_boot,
        slotSampleStats = boot_sample_stats,
        slotData = new_data
      ), silent = FALSE))
    }
    if (inherits(fit_boot, "try-error")) {
      if (show_progress) {
        cat("     FAILED: with ERROR message\n")
      }
      out <- as.numeric(NA)
      attr(out, "nonadmissible.flag") <- TRUE
      if (keep_idx) {
        attr(out, "BOOT.idx") <- boot_idx_1
      }
      return(out)
    }
    if (!fit_boot@optim$converged) {
      if (show_progress) {
        cat("     FAILED: no convergence\n")
      }
      out <- as.numeric(NA)
      attr(out, "nonadmissible.flag") <- TRUE
      if (keep_idx) {
        attr(out, "BOOT.idx") <- boot_idx_1
      }
      return(out)
    }

    # extract information we need
    if (is.null(object)) { # internal use only!
      if (fun == "coef") {
        out <- fit_boot@optim$x
      } else if (fun == "test") {
        out <- fit_boot@test[[1L]]$stat
      } else if (fun == "coeftest") {
        out <- c(fit_boot@optim$x, fit_boot@test[[1L]]$stat)
      }
    } else { # general use
      out <- try(as.numeric(fun(fit_boot, ...)), silent = TRUE)
    }
    if (inherits(out, "try-error")) {
      if (show_progress) {
        cat("     FAILED: applying fun to fit.boot\n")
      }
      out <- as.numeric(NA)
      attr(out, "nonadmissible.flag") <- TRUE
      if (keep_idx) {
        attr(out, "BOOT.idx") <- boot_idx_1
      }
      return(out)
    }

    # check if the solution is admissible
    admissible_flag <- suppressWarnings(lavInspect(fit_boot, "post.check"))
    attr(out, "nonadmissible.flag") <- !admissible_flag

    if (show_progress) {
      cat(
        "   OK -- niter = ",
        sprintf("%3d", fit_boot@optim$iterations), " fx = ",
        sprintf("%11.9f", fit_boot@optim$fx),
        # proposed by issue #519 Shu Fai Cheung
        if (admissible_flag) " " else "n", "",
        if (interactive()) "" else "\n"
      )
    }

    if (keep_idx) {
      # add BOOT.idx (for all groups)
      attr(out, "BOOT.idx") <- boot_idx_1
    }

    out
  } # end-of-fn

  # get parallelization options
  parallel <- lavoptions_1$parallel[1]
  ncpus <- lavoptions_1$ncpus
  cl <- lavoptions_1[["cl"]] # often NULL
  iseed <- lavoptions_1[["iseed"]] # often NULL

  # the next 8 lines are borrowed from the boot package
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") {
      have_mc <- .Platform$OS.type != "windows"
    } else if (parallel == "snow") {
      have_snow <- TRUE
    }
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
    temp_seed <- NA
    iseed <- runif(1, 0, 999999999)
  } else {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      temp_seed <-
        get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    } else {
      # is.null(temp.seed): Will remove .Random.seed in GlobalEnv
      #                     if serial.
      #                     If parallel, .Random.seed will not be touched.
      temp_seed <- NULL
    }
  }
  if (!(ncpus > 1L && (have_mc || have_snow))) { # Only for serial
    set.seed(iseed)
  }


  # this is adapted from the boot function in package boot
  rr <- r
  if (show_progress) {
    cat("\n")
  }
  res <- if (ncpus > 1L && (have_mc || have_snow)) {
    if (have_mc) {
      rngkind_old <- RNGkind() # store current kind
      RNGkind("L'Ecuyer-CMRG") # to allow for reproducible results
      set.seed(iseed)
      parallel::mclapply(seq_len(rr), fn, mc.cores = ncpus)
    } else if (have_snow) {
      list(...) # evaluate any promises
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
        # # No need for
        # if(RNGkind()[1L] == "L'Ecuyer-CMRG")
        # clusterSetRNGStream() always calls `RNGkind("L'Ecuyer-CMRG")`
        parallel::clusterSetRNGStream(cl, iseed = iseed)
        res <- parallel::parLapply(cl, seq_len(rr), fn)
        parallel::stopCluster(cl)
        res
      } else {
        parallel::parLapply(cl, seq_len(rr), fn)
      }
    }
  } else {
    lapply(seq_len(rr), fn)
  }
  if (show_progress) {
    cat("\n")
  }
  # restore old RNGkind()
  if (ncpus > 1L && have_mc) {
    RNGkind(rngkind_old[1], rngkind_old[2], rngkind_old[3])
  }

  # fill in container
  t_star[] <- do.call("rbind", res)

  # handle errors
  error_idx <- which(sapply(res, function(x) is.na(x[1L])))
  attr(t_star, "error.idx") <- error_idx # could be integer(0L)

  # handle nonadmissible solutions
  if (check_post) {
    notok <- which(sapply(res, attr, "nonadmissible.flag"))
    if (length(error_idx) > 0L) {
      notok <- notok[-which(notok %in% error_idx)]
    }
    attr(t_star, "nonadmissible") <- notok
  }

  # store iseed
  attr(t_star, "seed") <- iseed

  # handle temp.seed
  if (!is.null(temp_seed) && !identical(temp_seed, NA)) {
    assign(".Random.seed", temp_seed, envir = .GlobalEnv)                   # nolint start
  } else if (is.null(temp_seed) && !(ncpus > 1L && (have_mc || have_snow))) {
    # serial
    rm(.Random.seed, pos = 1)
  } else if (is.null(temp_seed) && (ncpus > 1L && have_mc)) {
    # parallel/multicore only
    rm(.Random.seed, pos = 1) # because set used set.seed()                 # nolint end
  }

  # store BOOT.idx per group
  if (keep_idx) {
    boot_idx_1 <- vector("list", length = lavsamplestats_1@ngroups)
    for (g in 1:lavsamplestats_1@ngroups) {
      # note that failed runs (NULL) are removed (for now)
      boot_idx_1[[g]] <- do.call(
        "rbind",
        lapply(res, function(x) attr(x, "BOOT.idx")[[g]])
      )
    }
    attr(t_star, "boot.idx") <- boot_idx_1
  }

  #    # No use, as boot package stores the sample indices differently
  #    # See boot:::boot.array() versus lav_bootstrap_indices()
  #    if(return.boot) {
  #        # mimic output boot function
  #
  #        if(is.null(object)) {
  #            stop("lavaan ERROR: return.boot = TRUE requires a full
  #                  lavaan object")
  #        }
  #
  #        # we start with ordinary only for now
  #        stopifnot(type == "ordinary")
  #
  #        if(! type %in% c("ordinary", "parametric")) {
  #            stop("lavaan ERROR: only ordinary and parametric bootstrap are
  #                  supported if return.boot = TRUE")
  #        } else {
  #            sim <- type
  #        }
  #
  #        statistic. <- function(data, idx) {
  #                         data.boot <- data[idx,]
  #                         fit.boot <- update(object, data = data.boot)
  #                         out <- try(fun(fit.boot, ...), silent = TRUE)
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
  t_star
}

# create matrix with indices to reconstruct the bootstrap samples
# per group
# (originally needed for BCa confidence intervals)
#
# rows are the (R) bootstrap runs
# columns are the (N) observations
#
# simple version: no strata, no weights
#
lav_bootstrap_indices <- function(r = 0L,
                                  nobs = list(0L), # per group
                                  parallel = "no",
                                  ncpus = 1L,
                                  cl = NULL,
                                  iseed = NULL,
                                  merge_groups = FALSE,
                                  return_freq = FALSE) {
  # iseed must be set!
  stopifnot(!is.null(iseed))

  if (return_freq && !merge_groups) {
    lav_msg_stop(gettext("return.freq only available if merge.groups = TRUE"))
  }

  if (is.integer(nobs)) {
    nobs <- list(nobs)
  }

  # number of groups
  ngroups <- length(nobs)

  # mimic 'random' sampling from lav_bootstrap_internal:

  # the next 7 lines are borrowed from the boot package
  have_mc <- have_snow <- FALSE
  parallel <- parallel[1]
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") {
      have_mc <- .Platform$OS.type != "windows"
    } else {
      if (parallel == "snow") have_snow <- TRUE
    }
    if (!have_mc && !have_snow) ncpus <- 1L
    loadNamespace("parallel") # before recording seed!
  }
  temp_seed <- NULL
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    temp_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  if (!(ncpus > 1L && (have_mc || have_snow))) { # Only for serial
    set.seed(iseed)
  }

  # fn() returns indices per group
  fn <- function(b) {
    boot_idx_1 <- vector("list", length = ngroups)
    offset_1 <- cumsum(c(0, unlist(nobs)))
    for (g in 1:ngroups) {
      stopifnot(nobs[[g]] > 1L)
      boot_idx <- sample.int(nobs[[g]], replace = TRUE)
      if (merge_groups) {
        boot_idx_1[[g]] <- boot_idx + offset_1[g]
      } else {
        boot_idx_1[[g]] <- boot_idx
      }
    }
    boot_idx_1
  }

  rr <- r
  res <- if (ncpus > 1L && (have_mc || have_snow)) {
    if (have_mc) {
      rngkind_old <- RNGkind() # store current kind
      RNGkind("L'Ecuyer-CMRG") # to allow for reproducible results
      set.seed(iseed)
      parallel::mclapply(seq_len(rr), fn, mc.cores = ncpus)
    } else if (have_snow) {
      # list(...) # evaluate any promises
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
        parallel::clusterSetRNGStream(cl, iseed = iseed)
        res <- parallel::parLapply(cl, seq_len(rr), fn)
        parallel::stopCluster(cl)
        res
      } else {
        parallel::parLapply(cl, seq_len(rr), fn)
      }
    }
  } else {
    lapply(seq_len(rr), fn)
  }

  # restore old RNGkind()
  if (ncpus > 1L && have_mc) {
    RNGkind(rngkind_old[1], rngkind_old[2], rngkind_old[3])
  }

  # handle temp.seed
  if (!is.null(temp_seed) && !identical(temp_seed, NA)) {            # nolint start
    assign(".Random.seed", temp_seed, envir = .GlobalEnv)
  } else if (is.null(temp_seed) && !(ncpus > 1L && (have_mc || have_snow))) {
    # serial
    rm(.Random.seed, pos = 1)
  } else if (is.null(temp_seed) && (ncpus > 1L && have_mc)) {
    # parallel/multicore only
    rm(.Random.seed, pos = 1) # because set used set.seed()          # nolint end
  }


  # assemble IDX
  boot_idx_1 <- vector("list", length = ngroups)
  for (g in 1:ngroups) {
    # FIXME: handle failed runs
    boot_idx_1[[g]] <- do.call("rbind", lapply(res, "[[", g))
  }

  # merge groups
  if (merge_groups) {
    out <- do.call("cbind", boot_idx_1)
  } else {
    out <- boot_idx_1
  }

  # NOTE: the order of the indices is different from the boot package!
  # we fill in the matrix 'row-wise' (1 row = sample(N, replace = TRUE)),
  # while boot fills in the matrix 'column-wise'
  # this also explains why we get different results with return.boot = TRUE
  # despite using the same iseed

  # return frequencies instead?
  if (return_freq && merge_groups) {
    out <- t(apply(out, 1L, tabulate, ncol(out)))
  }

  out
}
