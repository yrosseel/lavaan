# lavaanList: fit the *same* model, on different datasets
# YR - 29 June 2016

lavaanList <- function(model         = NULL,             # model
                       dataList      = NULL,             # list of datasets
                       dataFunction  = NULL,             # generating function
                       dataFunction.args = list(),       # optional arguments
                       ndat          = length(dataList), # how many datasets?
                       cmd           = "lavaan",
                       ...,
                       store.slots   = c("partable"),    # default is partable
                       FUN           = NULL,             # arbitrary FUN
                       show.progress = FALSE,
                       parallel      = c("no", "multicore", "snow"),
                       ncpus         = max(1L, parallel::detectCores() - 1L),
                       cl            = NULL) {

    # store.slots call
    mc  <- match.call()

    # check store.slots
    store.slots <- tolower(store.slots)
    if(length(store.slots) == 1L && store.slots == "all") {
        store.slots <- c("timing", "partable", "data", "samplestats",
                         "vcov", "test", "optim", "implied")
    }

    # dataList or function?
    if(is.function(dataFunction)) {
        if(ndat == 0L) {
            stop("lavaan ERROR: please specify number of requested datasets (ndat)")
        }
        firstData <- do.call(dataFunction, args = dataFunction.args)
        #dataList  <- vector("list", length = ndat)
    } else {
        firstData <- dataList[[1]]
    } 

    # check data
    if(is.matrix(firstData)) {
        # check if we have column names?
        NAMES <- colnames(firstData)
        if(is.null(NAMES)) {
            stop("lavaan ERROR: data is a matrix without column names")
        }
    } else if(inherits(firstData, "data.frame")) {
        # check?
    } else {
        stop("lavaan ERROR: (generated) data is not a data.frame (or a matrix)")
    }

    # parallel (see boot package)
    if (missing(parallel)) {
        #parallel <- getOption("boot.parallel", "no")
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

    # check dotdotdot
    dotdotdot$do.fit  <- TRUE    # to get starting values
    dotdotdot$se      <- "none"
    dotdotdot$test    <- "none"
    dotdotdot$verbose <- FALSE
    dotdotdot$debug   <- FALSE

    # initial model fit, using first dataset
    FIT <- do.call(cmd,
                   args = c(list(model  = model,
                                 data   = firstData), dotdotdot) )

    lavoptions  <- FIT@Options
    lavmodel    <- FIT@Model
    lavpartable <- FIT@ParTable

    # remove start/est/se columns from lavpartable
    lavpartable$start <- lavpartable$est
    lavpartable$est <- lavpartable$se <- NULL

    # empty slots
    timingList <- ParTableList <- DataList <- SampleStatsList <-
        vcovList <- testList <- optimList <- impliedList <- funList <- list()

    # prepare store.slotsd slots
    if("timing" %in% store.slots) {
        timingList <- vector("list", length = ndat)
    }
    if("partable" %in% store.slots) {
        ParTableList <- vector("list", length = ndat)
    }
    if("data" %in% store.slots) {
        DataList <- vector("list", length = ndat)
    }
    if("samplestats" %in% store.slots) {
        SampleStatsList <- vector("list", length = ndat)
    }
    if("vcov" %in% store.slots) {
        vcovList <- vector("list", length = ndat)
    }
    if("test" %in% store.slots) {
        testList <- vector("list", length = ndat)
    }
    if("optim" %in% store.slots) {
        optimList <- vector("list", length = ndat)
    }
    if("implied" %in% store.slots) {
        impliedList <- vector("list", length = ndat)
    }
    if(!is.null(FUN)) {
        funList <- vector("list", length = ndat)
    }

    # single run
    fn <- function(i) {

        if(show.progress) {
            cat("   ... data set number:", sprintf("%4d", i))
        }

        # get new dataset
        if(i == 1L) {
            DATA <- firstData
        } else if(is.function(dataFunction)) {
            DATA <- do.call(dataFunction, args = dataFunction.args)
        } else if(is.list(dataList)) {
            DATA <- dataList[[i]]
        }

        # adapt lavmodel for this new dataset
        # - starting values will be different
        # - ov.x variances/covariances
        # FIXME: can we not make the changes internally?
        #if(lavmodel@fixed.x && length(vnames(lavpartable, "ov.x")) > 0L) {
            #for(g in 1:FIT@Data@ngroups) {
            #    
            #}           
            lavmodel <- NULL
        #} 

        # fit model with this (new) dataset
        lavobject <- try(lavaan(slotOptions  = lavoptions,
                                slotParTable = lavpartable,
                                slotModel    = lavmodel,
                                start        = FIT,
                                data         = DATA), silent = TRUE)

        RES <- list(ok = FALSE, timing = NULL, ParTable = NULL,
                    Data = NULL, SampleStats = NULL, vcov = NULL,
                    test = NULL, optim = NULL, implied = NULL,
                    fun = NULL)

        if(inherits(lavobject, "lavaan") && 
           lavInspect(lavobject, "converged")) {
            RES$ok <- TRUE

            if(show.progress) {
                cat("   OK -- niter = ",
                        sprintf("%3d", lavInspect(lavobject, "iterations")), 
                        "\n")
            }

            # extract slots from fit
            if("timing" %in% store.slots) {
                RES$timing <- lavobject@timing
            }
            if("partable" %in% store.slots) {
                RES$ParTable <- lavobject@ParTable
            }
            if("data" %in% store.slots) {
                RES$Data <- lavobject@Data
            }
            if("samplestats" %in% store.slots) {
                RES$SampleStats <- lavobject@vcov
            }
            if("vcov" %in% store.slots) {
                RES$vcov <- lavobject@vcov
            }
            if("test" %in% store.slots) {
                RES$test <- lavobject@test
            }
            if("optim" %in% store.slots) {
                RES$optim <- lavobject@optim
            }
            if("implied" %in% store.slots) {
                RES$implied <- lavobject@implied
            }
 
            # custom FUN
            if(!is.null(FUN)) {
                RES$fun <- FUN(lavobject)
            }

            #if("coef" %in% output) {
            #    COEF[[i]] <- coef(lavobject)
            #}
        } else {
            if(show.progress) {
                cat("     FAILED: no convergence\n")
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
                if(RNGkind()[1L] == "L'Ecuyer-CMRG") {
                    parallel::clusterSetRNGStream(cl)
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
    meta <- list(ndat = ndat, ok = sapply(RES, "[[", "ok"),
                 store.slots = store.slots)

    # extract store.slots slots
    if("timing" %in% store.slots) {
        timingList <- lapply(RES, "[[", "timing")
    }
    if("partable" %in% store.slots) {
        ParTableList <- lapply(RES, "[[", "ParTable")
    }
    if("data" %in% store.slots) {
        DataList <- lapply(RES, "[[", "Data")
    }
    if("samplestats" %in% store.slots) {
        SampleStatsList <- lapply(RES, "[[", "SampleStats")
    }
    if("vcov" %in% store.slots) {
        vcovList <- lapply(RES, "[[", "vcov")
    }
    if("test" %in% store.slots) {
        testList <- lapply(RES, "[[", "test")
    }
    if("optim" %in% store.slots) {
        optimList <- lapply(RES, "[[", "optim")
    }
    if("implied" %in% store.slots) {
        impliedList <- lapply(RES, "[[", "implied")
    }
    if(!is.null(FUN)) {
        funList <- lapply(RES, "[[", "fun")
    }

    # create lavaanList object
    lavaanList <- new("lavaanList",
                      call         = mc,
                      Options      = lavoptions,
                      ParTable     = lavpartable,
                      Model        = lavmodel,
                      Data         = FIT@Data,

                      # meta
                      meta         = meta,

                      # per dataset
                      timingList      = timingList,
                      ParTableList    = ParTableList,
                      DataList        = DataList,
                      SampleStatsList = SampleStatsList,
                      vcovList        = vcovList,
                      testList        = testList,
                      optimList       = optimList,
                      impliedList     = impliedList,
                      funList         = funList,

                      external     = list()
                     )

    lavaanList
}

semList <- function(model         = NULL,
                    dataList      = NULL,
                    dataFunction  = NULL,
                    dataFunction.args = list(),
                    ndat          = length(dataList),
                    ...,
                    store.slots   = c("partable"),
                    FUN           = NULL,
                    show.progress = FALSE) {

    lavaanList(model = model, dataList = dataList, 
               dataFunction = dataFunction, 
               dataFunction.args = dataFunction.args, ndat = ndat, 
               cmd = "sem",
               ..., store.slots = store.slots, 
               FUN = FUN, show.progress = show.progress)
}

cfaList <- function(model         = NULL,             
                    dataList      = NULL,             
                    dataFunction  = NULL,             
                    dataFunction.args = list(),       
                    ndat          = length(dataList), 
                    ...,
                    store.slots   = c("partable"),    
                    FUN           = NULL,             
                    show.progress = FALSE) {

    lavaanList(model = model, dataList = dataList, 
               dataFunction = dataFunction, 
               dataFunction.args = dataFunction.args, ndat = ndat, 
               cmd = "cfa",
               ..., store.slots = store.slots, 
               FUN = FUN, show.progress = show.progress)
}

