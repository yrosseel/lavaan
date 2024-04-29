# lavaanList: fit the *same* model, on different datasets
# YR - 29 Jun 2016
# YR - 27 Jan 2017: change lavoptions; add dotdotdot to each call
# TDJ - 23 Aug 2018: change wrappers to preserve arguments from match.call()

lavaanList <- function(model = NULL, # model
                       dataList = NULL, # list of datasets
                       dataFunction = NULL, # generating function
                       dataFunction.args = list(), # optional arguments
                       ndat = length(dataList), # how many datasets?
                       cmd = "lavaan",
                       ...,
                       store.slots = c("partable"), # default is partable
                       FUN = NULL, # arbitrary FUN
                       show.progress = FALSE,
                       store.failed = FALSE,
                       parallel = c("no", "multicore", "snow"),
                       ncpus = max(1L, parallel::detectCores() - 1L),
                       cl = NULL,
                       iseed = NULL) {
  # store.slots call
  mc <- match.call()

  # check store.slots
  store.slots <- tolower(store.slots)
  if (length(store.slots) == 1L && store.slots == "all") {
    store.slots <- c(
      "timing", "partable", "data", "samplestats",
      "cache", "loglik", "h1", "baseline", "external",
      "vcov", "test", "optim", "implied"
    )
  }

  # dataList or function?
  if (is.function(dataFunction)) {
    if (ndat == 0L) {
      lav_msg_stop(gettext("please specify number of requested datasets (ndat)"))
    }
    firstData <- do.call(dataFunction, args = dataFunction.args)
    # dataList  <- vector("list", length = ndat)
  } else {
    firstData <- dataList[[1]]
  }

  # check data
  if (is.matrix(firstData)) {
    # check if we have column names?
    NAMES <- colnames(firstData)
    if (is.null(NAMES)) {
      lav_msg_stop(gettext("data is a matrix without column names"))
    }
  } else if (inherits(firstData, "data.frame")) {
    # check?
  } else {
    lav_msg_stop(gettext("(generated) data is not a data.frame (or a matrix)"))
  }

  # parallel (see boot package)
  if (missing(parallel)) {
    # parallel <- getOption("boot.parallel", "no")
    parallel <- "no"
  }
  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") {
      have_mc <- .Platform$OS.type != "windows"
    } else if (parallel == "snow") {
      have_snow <- TRUE
    }
    if (!have_mc && !have_snow) {
      ncpus <- 1L
    }
    loadNamespace("parallel")
  }

  # dot dot dot
  dotdotdot <- list(...)

  # if 'model' is a lavaan object (perhaps from lavSimulate), no need to
  # call `cmd'
  if (inherits(model, "lavaan")) {
    FIT <- model
  } else {
    # adapt for FIT
    # dotdotdotFIT <- dotdotdot
    # dotdotdotFIT$do.fit  <- TRUE    # to get starting values
    # dotdotdotFIT$se      <- "none"
    # dotdotdotFIT$test    <- "none"
    # dotdotdotFIT$verbose <- FALSE
    # dotdotdotFIT$debug   <- FALSE

    # initial model fit, using first dataset
    FIT <- do.call(cmd,
      args = c(list(
        model = model,
        data = firstData
      ), dotdotdot)
    )
  }

  lavoptions <- FIT@Options
  lavmodel <- FIT@Model
  lavpartable <- FIT@ParTable
  lavpta <- FIT@pta

  # remove any options in lavoptions from dotdotdot
  if (length(dotdotdot) > 0L) {
    rm.idx <- which(names(dotdotdot) %in% names(lavoptions))
    if (length(rm.idx) > 0L) {
      dotdotdot <- dotdotdot[-rm.idx]
    }
  }

  # remove start/est/se columns from lavpartable
  lavpartable$start <- lavpartable$est <- lavpartable$se <- NULL

  # empty slots
  timingList <- ParTableList <- DataList <- SampleStatsList <-
    CacheList <- vcovList <- testList <- optimList <-
    h1List <- loglikList <- baselineList <-
    impliedList <- funList <- list()

  # prepare store.slotsd slots
  if ("timing" %in% store.slots) {
    timingList <- vector("list", length = ndat)
  }
  if ("partable" %in% store.slots) {
    ParTableList <- vector("list", length = ndat)
  }
  if ("data" %in% store.slots) {
    DataList <- vector("list", length = ndat)
  }
  if ("samplestats" %in% store.slots) {
    SampleStatsList <- vector("list", length = ndat)
  }
  if ("cache" %in% store.slots) {
    CacheList <- vector("list", length = ndat)
  }
  if ("vcov" %in% store.slots) {
    vcovList <- vector("list", length = ndat)
  }
  if ("test" %in% store.slots) {
    testList <- vector("list", length = ndat)
  }
  if ("optim" %in% store.slots) {
    optimList <- vector("list", length = ndat)
  }
  if ("implied" %in% store.slots) {
    impliedList <- vector("list", length = ndat)
  }
  if ("loglik" %in% store.slots) {
    loglikList <- vector("list", length = ndat)
  }
  if ("h1" %in% store.slots) {
    h1List <- vector("list", length = ndat)
  }
  if ("baseline" %in% store.slots) {
    baselineList <- vector("list", length = ndat)
  }

  if (!is.null(FUN)) {
    funList <- vector("list", length = ndat)
  }

  # single run
  fn <- function(i) {
    if (show.progress) {
      cat("   ... data set number:", sprintf("%4d", i))
    }

    # get new dataset
    if (i == 1L) {
      DATA <- firstData
    } else if (is.function(dataFunction)) {
      DATA <- do.call(dataFunction, args = dataFunction.args)
    } else if (is.list(dataList)) {
      DATA <- dataList[[i]]
    }

    # if categorical, check if we have enough response categories
    # for each ordered variables in DATA
    data.ok.flag <- TRUE
    if (FIT@Model@categorical) {
      # expected nlev
      ord.idx <- unique(unlist(FIT@pta$vidx$ov.ord))
      NLEV.exp <- FIT@Data@ov$nlev[ord.idx]
      # observed nlev
      NLEV.obs <- sapply(
        DATA[, unique(unlist(FIT@pta$vnames$ov.ord)),
          drop = FALSE
        ],
        function(x) length(unique(na.omit(x)))
      )
      wrong.idx <- which(NLEV.exp - NLEV.obs != 0)
      if (length(wrong.idx) > 0L) {
        data.ok.flag <- FALSE
      }
    }

    # adapt lavmodel for this new dataset
    # - starting values will be different
    # - ov.x variances/covariances
    # FIXME: can we not make the changes internally?
    # if(lavmodel@fixed.x && length(vnames(lavpartable, "ov.x")) > 0L) {
    # for(g in 1:FIT@Data@ngroups) {
    #
    # }
    lavmodel <- NULL
    # }

    # fit model with this (new) dataset
    if (data.ok.flag) {
      if (cmd %in% c("lavaan", "sem", "cfa", "growth")) {
        # lavoptions$start <- FIT # FIXME: needed?
        lavobject <- try(
          do.call("lavaan",
            args = c(
              list(
                slotOptions = lavoptions,
                slotParTable = lavpartable,
                slotModel = lavmodel,
                # start        = FIT,
                data = DATA
              ),
              dotdotdot
            )
          ),
          silent = TRUE
        )
      } else if (cmd == "fsr") {
        # extract fs.method and fsr.method from dotdotdot
        if (!is.null(dotdotdot$fs.method)) {
          fs.method <- dotdotdot$fs.method
        } else {
          fs.method <- formals(fsr)$fs.method # default
        }

        if (!is.null(dotdotdot$fsr.method)) {
          fsr.method <- dotdotdot$fsr.method
        } else {
          fsr.method <- formals(fsr)$fsr.method # default
        }

        lavoptions$start <- FIT # FIXME: needed?
        lavobject <- try(
          do.call("fsr",
            args = c(
              list(
                slotOptions = lavoptions,
                slotParTable = lavpartable,
                slotModel = lavmodel,
                # start        = FIT,
                data = DATA,
                cmd = "lavaan",
                fs.method = fs.method,
                fsr.method = fsr.method
              ),
              dotdotdot
            )
          ),
          silent = TRUE
        )
      } else {
        lav_msg_stop(gettext("unknown cmd:"), cmd)
      }
    } # data.ok.flag

    RES <- list(
      ok = FALSE, timing = NULL, ParTable = NULL,
      Data = NULL, SampleStats = NULL, vcov = NULL,
      test = NULL, optim = NULL, implied = NULL,
      baseline = NULL, baseline.ok = FALSE, fun = NULL
    )

    if (data.ok.flag && inherits(lavobject, "lavaan") &&
      lavInspect(lavobject, "converged")) {
      RES$ok <- TRUE

      if (show.progress) {
        cat(
          "   OK -- niter = ",
          sprintf("%3d", lavInspect(lavobject, "iterations")),
          "\n"
        )
      }

      # extract slots from fit
      if ("timing" %in% store.slots) {
        RES$timing <- lavobject@timing
      }
      if ("partable" %in% store.slots) {
        RES$ParTable <- lavobject@ParTable
      }
      if ("data" %in% store.slots) {
        RES$Data <- lavobject@Data
      }
      if ("samplestats" %in% store.slots) {
        RES$SampleStats <- lavobject@SampleStats
      }
      if ("cache" %in% store.slots) {
        RES$Cache <- lavobject@Cache
      }
      if ("vcov" %in% store.slots) {
        RES$vcov <- lavobject@vcov
      }
      if ("test" %in% store.slots) {
        RES$test <- lavobject@test
      }
      if ("optim" %in% store.slots) {
        RES$optim <- lavobject@optim
      }
      if ("implied" %in% store.slots) {
        RES$implied <- lavobject@implied
      }
      if ("loglik" %in% store.slots) {
        RES$loglik <- lavobject@loglik
      }
      if ("h1" %in% store.slots) {
        RES$h1 <- lavobject@h1
      }
      if ("baseline" %in% store.slots) {
        RES$baseline <- lavobject@baseline
        if (length(lavobject@baseline) > 0L) {
          RES$baseline.ok <- TRUE
        }
      }

      # custom FUN
      if (!is.null(FUN)) {
        RES$fun <- FUN(lavobject)
      }
    } else { # failed!
      if (show.progress) {
        if (data.ok.flag) {
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
      if ("partable" %in% store.slots) {
        RES$ParTable <- lavpartable
        RES$ParTable$est <- RES$ParTable$start
        RES$ParTable$est[RES$ParTable$free > 0] <- as.numeric(NA)
        RES$ParTable$se <- numeric(length(lavpartable$lhs))
        RES$ParTable$se[RES$ParTable$free > 0] <- as.numeric(NA)
      }
      if (store.failed) {
        tmpfile <- tempfile(pattern = "lavaanListData")
        datfile <- paste0(tmpfile, ".csv")
        write.csv(DATA, file = datfile, row.names = FALSE)
        if (data.ok.flag) {
          # or only if lavobject is of class lavaan?
          objfile <- paste0(tmpfile, ".RData")
          write(lavobject, file = objfile)
        }
      }
    }

    RES
  }

  # the next 20 lines are based on the boot package
  RES <- if (ncpus > 1L && (have_mc || have_snow)) {
    if (have_mc) {
      parallel::mclapply(seq_len(ndat), fn, mc.cores = ncpus)
    } else if (have_snow) {
      list(...) # evaluate any promises
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", ncpus))
        if (RNGkind()[1L] == "L'Ecuyer-CMRG") {
          parallel::clusterSetRNGStream(cl, iseed = iseed)
        }
        RES <- parallel::parLapply(cl, seq_len(ndat), fn)
        parallel::stopCluster(cl)
        RES
      } else {
        parallel::parLapply(cl, seq_len(ndat), fn)
      }
    }
  } else {
    lapply(seq_len(ndat), fn)
  }


  # restructure
  if ("baseline" %in% store.slots) {
    meta <- list(
      ndat = ndat, ok = sapply(RES, "[[", "ok"),
      baseline.ok = sapply(RES, "[[", "baseline.ok"),
      store.slots = store.slots
    )
  } else {
    meta <- list(
      ndat = ndat, ok = sapply(RES, "[[", "ok"),
      store.slots = store.slots
    )
  }

  # extract store.slots slots
  if ("timing" %in% store.slots) {
    timingList <- lapply(RES, "[[", "timing")
  }
  if ("partable" %in% store.slots) {
    ParTableList <- lapply(RES, "[[", "ParTable")
  }
  if ("data" %in% store.slots) {
    DataList <- lapply(RES, "[[", "Data")
  }
  if ("samplestats" %in% store.slots) {
    SampleStatsList <- lapply(RES, "[[", "SampleStats")
  }
  if ("cache" %in% store.slots) {
    CacheList <- lapply(RES, "[[", "Cache")
  }
  if ("vcov" %in% store.slots) {
    vcovList <- lapply(RES, "[[", "vcov")
  }
  if ("test" %in% store.slots) {
    testList <- lapply(RES, "[[", "test")
  }
  if ("optim" %in% store.slots) {
    optimList <- lapply(RES, "[[", "optim")
  }
  if ("implied" %in% store.slots) {
    impliedList <- lapply(RES, "[[", "implied")
  }
  if ("h1" %in% store.slots) {
    h1List <- lapply(RES, "[[", "h1")
  }
  if ("loglik" %in% store.slots) {
    loglikList <- lapply(RES, "[[", "loglik")
  }
  if ("baseline" %in% store.slots) {
    baselineList <- lapply(RES, "[[", "baseline")
  }
  if (!is.null(FUN)) {
    funList <- lapply(RES, "[[", "fun")
  }

  # create lavaanList object
  lavaanList <- new("lavaanList",
    version = packageDescription("lavaan", fields = "Version"),
    call = mc,
    Options = lavoptions,
    ParTable = lavpartable,
    pta = lavpta,
    Model = lavmodel,
    Data = FIT@Data,

    # meta
    meta = meta,

    # per dataset
    timingList = timingList,
    ParTableList = ParTableList,
    DataList = DataList,
    SampleStatsList = SampleStatsList,
    CacheList = CacheList,
    vcovList = vcovList,
    testList = testList,
    optimList = optimList,
    impliedList = impliedList,
    h1List = h1List,
    loglikList = loglikList,
    baselineList = baselineList,
    funList = funList,
    external = list()
  )

  lavaanList
}

semList <- function(model = NULL,
                    dataList = NULL,
                    dataFunction = NULL,
                    dataFunction.args = list(),
                    ndat = length(dataList),
                    ...,
                    store.slots = c("partable"),
                    FUN = NULL,
                    show.progress = FALSE,
                    store.failed = FALSE,
                    parallel = c("no", "multicore", "snow"),
                    ncpus = max(1L, parallel::detectCores() - 1L),
                    cl = NULL,
                    iseed = NULL) {
  mc <- match.call(expand.dots = TRUE)
  mc$cmd <- "sem"
  mc[[1L]] <- quote(lavaan::lavaanList)
  eval(mc, parent.frame())
}

cfaList <- function(model = NULL,
                    dataList = NULL,
                    dataFunction = NULL,
                    dataFunction.args = list(),
                    ndat = length(dataList),
                    ...,
                    store.slots = c("partable"),
                    FUN = NULL,
                    show.progress = FALSE,
                    store.failed = FALSE,
                    parallel = c("no", "multicore", "snow"),
                    ncpus = max(1L, parallel::detectCores() - 1L),
                    cl = NULL,
                    iseed = NULL) {
  mc <- match.call(expand.dots = TRUE)
  mc$cmd <- "cfa"
  mc[[1L]] <- quote(lavaan::lavaanList)
  eval(mc, parent.frame())
}
