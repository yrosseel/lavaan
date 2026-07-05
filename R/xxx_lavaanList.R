# lavaanList: fit the *same* model, on different datasets
# YR - 29 Jun 2016
# YR - 27 Jan 2017: change lavoptions; add dotdotdot to each call
# TDJ - 23 Aug 2018: change wrappers to preserve arguments from match.call()
# YR - 15 Oct 2024: add iseed (as in lavBootstrap)

lavaanList <- function(model = NULL, # model                    # nolint
                       data_list = NULL, # list of datasets
                       data_function = NULL, # generating function
                       data_function_args = list(), # optional arguments
                       ndat = length(data_list), # how many datasets?
                       cmd = "lavaan",
                       ...,
                       store_slots = c("partable"), # default is partable
                       fun = NULL, # arbitrary fun
                       show_progress = FALSE,
                       store_failed = FALSE,
                       parallel = c("no", "multicore", "snow"),
                       ncpus = max(1L, parallel::detectCores() - 1L, na.rm = TRUE),
                       cl = NULL,
                       iseed = NULL) {
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, FALSE)

  # store_slots call
  mc <- match.call()

  # store current random seed (if any) and method
  rngkind_old <- RNGkind()
  #cat("RNGkind = ", paste(RNGkind_old, collapse = " "), "\n")
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    init_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else {
    init_seed <- NULL
  }

  # initial seed (for the precomputations)
  if (!is.null(iseed)) {
    set.seed(min(1, iseed - 1L))
  }

  # check store_slots
  store_slots <- tolower(store_slots)
  if (length(store_slots) == 1L && store_slots == "all") {
    store_slots <- c(
      "timing", "partable", "data", "samplestats",
      "cache", "loglik", "h1", "baseline", "external",
      "vcov", "test", "optim", "implied"
    )
  }

  # data_list or function?
  if (is.function(data_function)) {
    if (ndat == 0L) {
      lav_msg_stop(gettext(
        "please specify number of requested datasets (ndat)"))
    }
    # here we already use the random generator
    first_data <- do.call(data_function, args = data_function_args)
    # data_list  <- vector("list", length = ndat)
  } else {
    first_data <- data_list[[1]]
  }

  # check data
  if (is.matrix(first_data)) {
    # check if we have column names?
    names_1 <- colnames(first_data)
    if (is.null(names_1)) {
      lav_msg_stop(gettext("data is a matrix without column names"))
    }
  } else if (inherits(first_data, "data.frame")) {
    # check?
  } else {
    lav_msg_stop(gettext("(generated) data is not a data.frame (or a matrix)"))
  }

  # if 'model' is a lavaan object (perhaps from lavSimulate), no need to
  # call `cmd'
  if (inherits(model, "lavaan")) {
    fit <- model
  } else {
    # adapt for FIT
    # dotdotdotFIT <- dotdotdot
    # dotdotdotFIT$do.fit  <- TRUE    # to get starting values
    # dotdotdotFIT$se      <- "none"
    # dotdotdotFIT$test    <- "none"

    # initial model fit, using first dataset
    fit <- do.call(cmd,
      args = c(list(
        model = model,
        data = first_data
      ), dotdotdot)
    )
  }

  lavoptions <- fit@Options
  lavmodel <- fit@Model
  lavpartable <- fit@ParTable
  lavpta <- fit@pta

  # remove any options in lavoptions from dotdotdot (sem/lavaan only)
  # because we can use the lavoptions slot
  if (length(dotdotdot) > 0L && cmd %in% c("lavaan", "sem", "cfa", "growth")) {
    rm_idx <- which(names(dotdotdot) %in% names(lavoptions))
    if (length(rm_idx) > 0L) {
      dotdotdot <- dotdotdot[-rm_idx]
    }
  }

  # remove start/est/se columns from lavpartable
  lavpartable$start <- lavpartable$est <- lavpartable$se <- NULL

  # empty slots
  timing_list <- par_table_list <- local_data_list <- sample_stats_list <-
    cache_list <- vcov_list <- test_list <- optim_list <-
    h1list <- loglik_list <- baseline_list <-
    implied_list <- fun_list <- list()

  # prepare store_slots slots
  if ("timing" %in% store_slots) {
    timing_list <- vector("list", length = ndat)
  }
  if ("partable" %in% store_slots) {
    par_table_list <- vector("list", length = ndat)
  }
  if ("data" %in% store_slots) {
    local_data_list <- vector("list", length = ndat)
  }
  if ("samplestats" %in% store_slots) {
    sample_stats_list <- vector("list", length = ndat)
  }
  if ("cache" %in% store_slots) {
    cache_list <- vector("list", length = ndat)
  }
  if ("vcov" %in% store_slots) {
    vcov_list <- vector("list", length = ndat)
  }
  if ("test" %in% store_slots) {
    test_list <- vector("list", length = ndat)
  }
  if ("optim" %in% store_slots) {
    optim_list <- vector("list", length = ndat)
  }
  if ("implied" %in% store_slots) {
    implied_list <- vector("list", length = ndat)
  }
  if ("loglik" %in% store_slots) {
    loglik_list <- vector("list", length = ndat)
  }
  if ("h1" %in% store_slots) {
    h1list <- vector("list", length = ndat)
  }
  if ("baseline" %in% store_slots) {
    baseline_list <- vector("list", length = ndat)
  }

  if (!is.null(fun)) {
    fun_list <- vector("list", length = ndat)
  }

  # single run
  fn <- function(i) {
    if (show_progress) {
      cat("   ... data set number:", sprintf("%4d", i))
    }

    # get new dataset
    if (i == 1L) {
      data_1 <- first_data
    } else if (is.function(data_function)) {
      data_1 <- do.call(data_function, args = data_function_args)
    } else if (is.list(data_list)) {
      data_1 <- data_list[[i]]
    }

    # if categorical, check if we have enough response categories
    # for each ordered variables in DATA
    data_ok_flag <- TRUE
    if (fit@Model@categorical) {
      # expected nlev
      ord_idx <- unique(unlist(fit@pta$vidx$ov.ord))
      nlev_exp <- fit@Data@ov$nlev[ord_idx]
      # observed nlev
      nlev_obs <- sapply(
        data_1[, unique(unlist(fit@pta$vnames$ov.ord)),
          drop = FALSE
        ],
        function(x) length(unique(na.omit(x)))
      )
      wrong_idx <- which(nlev_exp - nlev_obs != 0)
      if (length(wrong_idx) > 0L) {
        data_ok_flag <- FALSE
      }
    }

    # adapt lavmodel for this new dataset
    # - starting values will be different
    # - ov.x variances/covariances
    # FIXME: can we not make the changes internally?
    # if(lavmodel@fixed.x &&
    #     length(lav_pt_vnames(lavpartable, "ov.x")) > 0L) {
    # for(g in 1:FIT@Data@ngroups) {
    #
    # }
    lavmodel <- NULL
    # }

    # fit model with this (new) dataset
    if (data_ok_flag) {
      if (cmd %in% c("lavaan", "sem", "cfa", "growth")) {
        # lavoptions$start <- FIT # FIXME: needed?
        lavobject <- try(
          do.call("lavaan",
            args = c(
              list(
                slot_options = lavoptions,
                slot_par_table = lavpartable,
                slot_model = lavmodel,
                # start        = FIT,
                data = data_1
              ),
              dotdotdot
            )
          ),
          silent = TRUE
        )
      } else if (cmd == "fsr") {
        # extract fs.method and fsr.method from dotdotdot
        if (!is.null(dotdotdot$fs.method)) {
          fs_method <- dotdotdot$fs.method
        } else {
          fs_method <- formals(fsr)$fs_method # default
        }

        if (!is.null(dotdotdot$fsr.method)) {
          fsr_method <- dotdotdot$fsr.method
        } else {
          fsr_method <- formals(fsr)$fsr_method # default
        }

        lavoptions$start <- fit # FIXME: needed?
        lavobject <- try(
          do.call("fsr",
            args = c(
              list(
                slot_options = lavoptions,
                slot_par_table = lavpartable,
                slot_model = lavmodel,
                # start        = FIT,
                data = data_1,
                cmd = "lavaan",
                fs.method = fs_method,
                fsr.method = fsr_method
              ),
              dotdotdot
            )
          ),
          silent = TRUE
        )
      } else if (cmd == "sam") {
        lavobject <- try(
          do.call("sam",
            args = c(
              list(
                model = model,
                data = data_1
              ),
              dotdotdot
            )
          ),
          silent = TRUE
        )
      } else {
        lavobject <- try(do.call(cmd,
          args = c(
              list(
                model = model,
                data = data_1
              ),
              dotdotdot
            )
          ),
          silent = TRUE
        )
        #lav_msg_stop(gettext("unknown cmd:"), cmd)
      }
    } # data.ok.flag

    res <- list(
      ok = FALSE, timing = NULL, ParTable = NULL,
      Data = NULL, SampleStats = NULL, vcov = NULL,
      test = NULL, optim = NULL, implied = NULL,
      baseline = NULL, baseline.ok = FALSE, fun = NULL
    )

    if (data_ok_flag && inherits(lavobject, "lavaan") &&
      lavInspect(lavobject, "converged")) {
      res$ok <- TRUE

      if (show_progress) {
        cat(
          "   OK -- niter = ",
          sprintf("%3d", lavInspect(lavobject, "iterations")),
          "\n"
        )
      }

      # extract slots from fit
      if ("timing" %in% store_slots) {
        res$timing <- lavobject@timing
      }
      if ("partable" %in% store_slots) {
        res$ParTable <- lavobject@ParTable
        # if the (possibly asymmetric) confidence intervals are based on
        # the bootstrap or the Monte Carlo method, the internal ParTable
        # does not contain them; compute and store ci.lower/ci.upper here,
        # so summary() (eg of a lavSimulate() object) can evaluate the
        # coverage of these intervals. The intervals are computed at the
        # 95% level (the lavaan default).
        if (identical(lavobject@Options$se, "bootstrap") ||
            identical(lavobject@Options$se.def, "monte.carlo")) {
          pe_ci <- try(lavParameterEstimates(lavobject,
            se = TRUE, zstat = FALSE, pvalue = FALSE, ci = TRUE,
            level = 0.95, standardized = FALSE,
            remove_system_eq = FALSE, remove_eq = FALSE,
            remove_ineq = FALSE, remove_def = FALSE,
            remove_nonfree = FALSE, remove_step1 = FALSE,
            remove_unused = FALSE, output = "data.frame"),
            silent = TRUE)
          if (!inherits(pe_ci, "try-error") && !is.null(pe_ci$ci.lower) &&
              length(pe_ci$ci.lower) == length(res$ParTable$lhs)) {
            res$ParTable$ci.lower <- pe_ci$ci.lower
            res$ParTable$ci.upper <- pe_ci$ci.upper
          } else {
            # res$ParTable is a list (not a data.frame), so we must use a
            # full-length vector here, not a scalar
            res$ParTable$ci.lower <- rep(as.numeric(NA),
                                         length(res$ParTable$lhs))
            res$ParTable$ci.upper <- rep(as.numeric(NA),
                                         length(res$ParTable$lhs))
          }
        }
      }
      if ("data" %in% store_slots) {
        res$Data <- lavobject@Data
      }
      if ("samplestats" %in% store_slots) {
        res$SampleStats <- lavobject@SampleStats
      }
      if ("cache" %in% store_slots) {
        res$Cache <- lavobject@Cache
      }
      if ("vcov" %in% store_slots) {
        res$vcov <- lavobject@vcov
      }
      if ("test" %in% store_slots) {
        res$test <- lavobject@test
      }
      if ("optim" %in% store_slots) {
        res$optim <- lavobject@optim
      }
      if ("implied" %in% store_slots) {
        res$implied <- lavobject@implied
      }
      if ("loglik" %in% store_slots) {
        res$loglik <- lavobject@loglik
      }
      if ("h1" %in% store_slots) {
        res$h1 <- lavobject@h1
      }
      if ("baseline" %in% store_slots) {
        res$baseline <- lavobject@baseline
        if (length(lavobject@baseline) > 0L) {
          res$baseline.ok <- TRUE
        }
      }

      # custom fun
      if (!is.null(fun)) {
        res$fun <- fun(lavobject)
      }
    } else { # failed!
      if (show_progress) {
        if (data_ok_flag) {
          if (inherits(lavobject, "lavaan")) {
            cat("   FAILED: no convergence\n")
          } else {
            cat("   FAILED: could not construct lavobject\n")
            print(lavobject)
          }
        } else {
          cat("   FAILED: nlev too low for some vars\n")
        }
      }
      if ("partable" %in% store_slots) {
        res$ParTable <- lavpartable
        res$ParTable$est <- res$ParTable$start
        res$ParTable$est[res$ParTable$free > 0] <- as.numeric(NA)
        res$ParTable$se <- numeric(length(lavpartable$lhs))
        res$ParTable$se[res$ParTable$free > 0] <- as.numeric(NA)
        # keep the ci.lower/ci.upper columns in sync with the converged
        # runs (see above), so all stored ParTables have the same columns.
        # res$ParTable is a list, so use full-length vectors, not scalars.
        if (identical(lavoptions$se, "bootstrap") ||
            identical(lavoptions$se.def, "monte.carlo")) {
          res$ParTable$ci.lower <- rep(as.numeric(NA),
                                       length(lavpartable$lhs))
          res$ParTable$ci.upper <- rep(as.numeric(NA),
                                       length(lavpartable$lhs))
        }
      }
      if (store_failed) {
        tmpfile <- tempfile(pattern = "lavaanListData")
        datfile <- paste0(tmpfile, ".csv")
        write.csv(data_1, file = datfile, row.names = FALSE)
        if (data_ok_flag) {
          # or only if lavobject is of class lavaan?
          objfile <- paste0(tmpfile, ".RData")
          save(lavobject, file = objfile)
        }
      }
    }

    res
  }


  # the next 8 lines are borrowed from the boot package
  have_mc <- have_snow <- FALSE
  if (missing(parallel)) {
    parallel <- "no"
  }
  parallel <- match.arg(parallel)
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
  res <- if (ncpus > 1L && (have_mc || have_snow)) {
    if (have_mc) {
      #RNGkind_old <- RNGkind() # store current kind
      RNGkind("L'Ecuyer-CMRG") # to allow for reproducible results
      set.seed(iseed)
      parallel::mclapply(seq_len(ndat), fn, mc.cores = ncpus)
    } else if (have_snow) {
      list(...) # evaluate any promises
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
        # # No need for
        # if(RNGkind()[1L] == "L'Ecuyer-CMRG")
        # clusterSetRNGStream() always calls `RNGkind("L'Ecuyer-CMRG")`
        parallel::clusterSetRNGStream(cl, iseed = iseed)
        res <- parallel::parLapply(cl, seq_len(ndat), fn)
        parallel::stopCluster(cl)
        res
      } else {
        parallel::parLapply(cl, seq_len(ndat), fn)
      }
    }
  } else {
    lapply(seq_len(ndat), fn)
  }


  # restructure
  if ("baseline" %in% store_slots) {
    meta <- list(
      ndat = ndat, ok = sapply(res, "[[", "ok"),
      baseline.ok = sapply(res, "[[", "baseline.ok"),
      store.slots = store_slots
    )
  } else {
    meta <- list(
      ndat = ndat, ok = sapply(res, "[[", "ok"),
      store.slots = store_slots
    )
  }

  # extract store_slots slots
  if ("timing" %in% store_slots) {
    timing_list <- lapply(res, "[[", "timing")
  }
  if ("partable" %in% store_slots) {
    par_table_list <- lapply(res, "[[", "ParTable")
  }
  if ("data" %in% store_slots) {
    local_data_list <- lapply(res, "[[", "Data")
  }
  if ("samplestats" %in% store_slots) {
    sample_stats_list <- lapply(res, "[[", "SampleStats")
  }
  if ("cache" %in% store_slots) {
    cache_list <- lapply(res, "[[", "Cache")
  }
  if ("vcov" %in% store_slots) {
    vcov_list <- lapply(res, "[[", "vcov")
  }
  if ("test" %in% store_slots) {
    test_list <- lapply(res, "[[", "test")
  }
  if ("optim" %in% store_slots) {
    optim_list <- lapply(res, "[[", "optim")
  }
  if ("implied" %in% store_slots) {
    implied_list <- lapply(res, "[[", "implied")
  }
  if ("h1" %in% store_slots) {
    h1list <- lapply(res, "[[", "h1")
  }
  if ("loglik" %in% store_slots) {
    loglik_list <- lapply(res, "[[", "loglik")
  }
  if ("baseline" %in% store_slots) {
    baseline_list <- lapply(res, "[[", "baseline")
  }
  if (!is.null(fun)) {
    fun_list <- lapply(res, "[[", "fun")
  }

  # create lavaanList object
  lavaan_list <- new("lavaanList",
    version = packageDescription("lavaan", fields = "Version"),
    call = mc,
    Options = lavoptions,
    ParTable = lavpartable,
    pta = lavpta,
    Model = lavmodel,
    Data = fit@Data,

    # meta
    meta = meta,

    # per dataset
    timingList = timing_list,
    ParTableList = par_table_list,
    DataList = local_data_list,
    SampleStatsList = sample_stats_list,
    CacheList = cache_list,
    vcovList = vcov_list,
    testList = test_list,
    optimList = optim_list,
    impliedList = implied_list,
    h1List = h1list,
    loglikList = loglik_list,
    baselineList = baseline_list,
    funList = fun_list,
    external = list()
  )

  # restore random seed
  if (!is.null(init_seed)) {
    assign(".Random.seed", init_seed, envir = .GlobalEnv)    # nolint
  } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    # initially there was no .Random.seed, but we created one along the way
    # clean up
    remove(".Random.seed", envir = .GlobalEnv)
  }
  # restore random generator method (if changed)
  RNGkind(rngkind_old[1], rngkind_old[2], rngkind_old[3])

  lavaan_list
}

semList <- function(model = NULL,                            # nolint
                    data_list = NULL,
                    data_function = NULL,
                    data_function_args = list(),
                    ndat = length(data_list),
                    ...,
                    store_slots = c("partable"),
                    fun = NULL,
                    show_progress = FALSE,
                    store_failed = FALSE,
                    parallel = c("no", "multicore", "snow"),
                    ncpus = max(1L, parallel::detectCores() - 1L, na.rm = TRUE),
                    cl = NULL,
                    iseed = NULL) {
  mc <- match.call(expand.dots = TRUE)
  mc$cmd <- "sem"
  mc[[1L]] <- quote(lavaan::lavaanList)
  eval(mc, parent.frame())
}

cfaList <- function(model = NULL,                         # nolint
                    data_list = NULL,
                    data_function = NULL,
                    data_function_args = list(),
                    ndat = length(data_list),
                    ...,
                    store_slots = c("partable"),
                    fun = NULL,
                    show_progress = FALSE,
                    store_failed = FALSE,
                    parallel = c("no", "multicore", "snow"),
                    ncpus = max(1L, parallel::detectCores() - 1L, na.rm = TRUE),
                    cl = NULL,
                    iseed = NULL) {
  mc <- match.call(expand.dots = TRUE)
  mc$cmd <- "cfa"
  mc[[1L]] <- quote(lavaan::lavaanList)
  eval(mc, parent.frame())
}
