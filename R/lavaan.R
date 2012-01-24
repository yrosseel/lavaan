# main user-visi ble cfa/sem/growth functions 
#
# initial version: YR 25/03/2009
# added lavaanOptions YR 02/08/2010
# major revision: YR 9/12/2010: - new workflow (since 0.4-5)
#                               - merge cfa/sem/growth functions

lavaan <- function(# user-specified model: can be syntax, parameter Table, ...
                   model           = NULL,
                   model.type      = "sem",
                
                   # model modifiers
                   meanstructure   = "default",
                   int.ov.free     = FALSE,
                   int.lv.free     = FALSE,
                   fixed.x         = "default", # or FALSE?
                   orthogonal      = FALSE,
                   std.lv          = FALSE,

                   auto.fix.first  = FALSE,
                   auto.fix.single = FALSE,
                   auto.var        = FALSE,
                   auto.cov.lv.x   = FALSE,
                   auto.cov.y      = FALSE,
                   
                   # full data
                   data            = NULL,
                   std.ov          = FALSE,
                   missing         = "default",

                   # summary data
                   sample.cov      = NULL,
                   sample.mean     = NULL,
                   sample.nobs     = NULL,

                   # multiple groups
                   group           = NULL,
                   group.equal     = '',
                   group.partial   = '',
             
                   # constraints
                   constraints     = '',

                   # estimation
                   estimator       = "default",
                   likelihood      = "default",
                   information     = "default",
                   se              = "default",
                   test            = "default",
                   bootstrap       = 1000L,
                   mimic           = "default",
                   representation  = "default",
                   do.fit          = TRUE,
                   control         = list(),

                   # starting values
                   start           = "default",

                   # full slots from previous fits
                   slotOptions     = NULL,
                   slotUser        = NULL,
                   slotSample      = NULL,
                   slotData        = NULL,
                   slotModel       = NULL,
  
                   # verbosity
                   verbose         = FALSE,
                   warn            = TRUE,
                   debug           = FALSE
                  )
{
    # 0. store call, start timer
    mc  <- match.call()
    start.time0 <- start.time <- proc.time()[3]
    timing <- list()

    # 1a. check data/sample.cov and get the number of groups
    if(!is.null(slotSample)) {
        stopifnot(class(slotSample) == "SampleStats")
        if(is.null(slotData)) {
            lavaanData <- list()
        } else {
            # need further checking?
            lavaanData <- slotData
        }
        ngroups <- slotSample@ngroups
        data.type = "sampleStats"
    } else if(!is.null(data)) {
        # good, we got a full data frame
        stopifnot(is.data.frame(data))
        if(!is.null(group)) {
            if(!(group %in% names(data))) {
                stop("lavaan ERROR: grouping variable `", group,
                     "' not found in names data:", names(data))
            }
            # note: we use the order as in the data; not as in levels(data)
            group.label <- unique(as.character(data[,group]))
            if(warn && any(is.na(group.label))) {
                cat("lavaan WARNING: group variable `", group, "` contains missing values\n",
                    sep="")
            }
            group.label <- group.label[!is.na(group.label)]
            ngroups     <- length(group.label)
        } else {
            group.label <- character(0)
            ngroups     <- 1L
        }
        data.type <- "full"
    } else if(!is.null(sample.cov)) {

        # we also need the number of observations (per group)
        if(is.null(sample.nobs))
            stop("lavaan ERROR: please specify number of observations")

        # if meanstructure=TRUE, we need sample.mean
        if(meanstructure == TRUE && is.null(sample.mean))
            stop("lavaan ERROR: please provide sample.mean if meanstructure=TRUE")
        # if group.equal contains "intercepts", we need sample.mean
        if("intercepts" %in% group.equal && is.null(sample.mean))
            stop("lavaan ERROR: please provide sample.mean if group.equal contains \"intercepts\"")


        # list?
        if(is.list(sample.cov)) {
            # multiple groups, multiple cov matrices
            if(!is.null(sample.mean)) {
                stopifnot(length(sample.mean) == length(sample.cov))
            }
            # multiple groups, multiple cov matrices
            ngroups     <- length(sample.cov)
            group.label <- names(sample.cov)
            if(is.null(group.label)) {
                group.label <- paste("Group ", 1:ngroups, sep="")
            }
        } else {
            if(!is.matrix(sample.cov))
                stop("lavaan ERROR: sample.cov must be a matrix or a list of matrices")
            ngroups <- 1L
            group.label <- character(0)
        }
        data.type <- "moment"
    } else {
        # both data and sample.cov are NULL
        # maybe we want an empty lavaan object to simulate?

        # number of groups
        if(!is.null(sample.nobs)) {
            ngroups <- length(unlist(sample.nobs))
        } else {
            ngroups <- 1L
        }
        if(ngroups > 1L) {
            group.label <- paste("Group ", 1:ngroups, sep="")
        } else {
            group.label <- character(0)
        }
        data.type <- "none"

        # no data? no fitting!
        do.fit <- FALSE
        se <- "none"
        test <- "none"
        start <- "simple"
    }


    # 1a. collect various options/flags and fill in `default' values
    #opt <- modifyList(formals(lavaan), as.list(mc)[-1])
    # force evaluation of `language` and/or `symbol` arguments
    #opt <- lapply(opt, function(x) if(typeof(x) %in% c("language", "symbol")) 
    #                                   eval(x, parent.frame()) else x)
    if(!is.null(slotOptions)) {
        lavaanOptions <- slotOptions
    } else {
        opt <- list(model = model, model.type = model.type,
            meanstructure = meanstructure, int.ov.free = int.ov.free,
            int.lv.free = int.lv.free, fixed.x = fixed.x, 
            orthogonal = orthogonal, std.lv = std.lv, 
            auto.fix.first = auto.fix.first, auto.fix.single = auto.fix.single,
            auto.var = auto.var, auto.cov.lv.x = auto.cov.lv.x, 
            auto.cov.y = auto.cov.y, std.ov = std.ov, missing = missing, 
            group = group, group.equal = group.equal, 
            group.partial = group.partial, 
            constraints = constraints,
            estimator = estimator, likelihood = likelihood,
            information = information, se = se, test = test, 
            bootstrap = bootstrap, mimic = mimic,
            representation = representation, do.fit = do.fit, verbose = verbose,
            warn = warn, debug = debug, data.type = data.type)
        lavaanOptions <- setLavaanOptions(opt)
    }
    timing$InitOptions <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
    

    # 2a. construct lavaan User list: description of the user-specified model
    if(!is.null(slotUser)) {
        lavaanUser <- slotUser
    } else if(is.character(model)) {
        lavaanUser <- 
            lavaanify(model.syntax    = model, 
                      meanstructure   = lavaanOptions$meanstructure, 
                      int.ov.free     = lavaanOptions$int.ov.free,
                      int.lv.free     = lavaanOptions$int.lv.free,
                      orthogonal      = lavaanOptions$orthogonal, 
                      fixed.x         = lavaanOptions$fixed.x,
                      std.lv          = lavaanOptions$std.lv,
                      constraints     = constraints,

                      auto.fix.first  = lavaanOptions$auto.fix.first,
                      auto.fix.single = lavaanOptions$auto.fix.single,
                      auto.var        = lavaanOptions$auto.var,
                      auto.cov.lv.x   = lavaanOptions$auto.cov.lv.x,
                      auto.cov.y      = lavaanOptions$auto.cov.y,

                      ngroups         = ngroups,
                      group.equal     = lavaanOptions$group.equal, 
                      group.partial   = lavaanOptions$group.partial,
                      debug           = lavaanOptions$debug,
                      warn            = lavaanOptions$warn,

                      as.data.frame.  = FALSE)
    } else if(is.list(model)) {
        # two possibilities: either model is already lavaanified
        # or it is a list of modeles (perhaps one for each group)
        if(!is.null(model$lhs) &&
           !is.null(model$op)  &&
           !is.null(model$rhs) &&
           !is.null(model$free)) {
            lavaanUser <- model
        } else if(is.character(model[[1]])) {
            # we lavaanify each model in term, and assume multiple groups
            # ... or not: we now allow multiple groups (with different
            # manifest variables within a single syntax... (dec 2011)
            stop("lavaan ERROR: model is a list, but not a parameterTable?")
        }
    } else {
        cat("model type: ", class(model), "\n")
        stop("lavaan ERROR: model is not of type character or list")
    }

    # 2b. change meanstructure flag?
    if(any(lavaanUser$op == "~1")) lavaanOptions$meanstructure <- TRUE
    timing$User <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 3. get sample statistics
    ov.names <- lapply(as.list(1:ngroups),
                       function(x) vnames(lavaanUser, type="ov", x))

    # 3a. handle full data
    if(data.type == "sampleStats") {
        lavaanSampleStats <- slotSample
    } else if(data.type == "full") {
        # ov.names User model
        ov.names <- lapply(as.list(1:ngroups),
                       function(x) vnames(lavaanUser, type="ov", x))   
        lavaanData <- getData(data        = data, 
                              ov.names    = ov.names,
                              std.ov      = std.ov,
                              group       = group,
                              group.label = group.label)

        lavaanMissing <- 
            getMissingPatterns(X       = lavaanData, 
                               missing = lavaanOptions$missing,
                               warn    = lavaanOptions$warn,
                               verbose = lavaanOptions$verbose)

        WLS.V <- list()
        if(lavaanOptions$estimator %in% c("GLS", "WLS")) {
            WLS.V <- getWLS.V(X             = lavaanData,
                              sample        = NULL,
                              estimator     = lavaanOptions$estimator,
                              mimic         = lavaanOptions$mimic,
                              meanstructure = lavaanOptions$meanstructure)
        }

        lavaanSampleStats <- getSampleStatsFromData(
                           X           = lavaanData,
                           M           = lavaanMissing,
                           rescale     = (lavaanOptions$estimator == "ML" &&
                                          lavaanOptions$likelihood == "normal"),
                           group.label = group.label,
                           WLS.V       = WLS.V)
                                                 
    } else if(data.type == "moment") {
        lavaanData <- list()
        lavaanSampleStats <- getSampleStatsFromMoments(
                           sample.cov  = sample.cov,
                           sample.mean = sample.mean,
                           sample.nobs = sample.nobs,
                           ov.names    = ov.names,
                           rescale     = (lavaanOptions$estimator == "ML" &&
                                          lavaanOptions$likelihood == "normal"),
                           group.label = group.label)

        if(lavaanOptions$estimator == "GLS") {
            WLS.V <- getWLS.V(X             = NULL,
                              sample        = lavaanSampleStats,
                              estimator     = lavaanOptions$estimator,
                              mimic         = lavaanOptions$mimic,
                              meanstructure = lavaanOptions$meanstructure)
            lavaanSampleStats@WLS.V <- WLS.V
        }
    } else {
        # no data
        lavaanData <- list()
        lavaanSampleStats <- new("SampleStats", ngroups=ngroups,
                                 ov.names=ov.names)
    } 
    if(debug) {
        print(str(lavaanData))
        print(str(lavaanSampleStats))
    }


    timing$Sample <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 4. compute some reasonable starting values 
    if(!is.null(slotModel)) {
        lavaanModel <- slotModel
        lavaanStart <- getModelParameters(lavaanModel, type="user")
        timing$Start <- (proc.time()[3] - start.time)
        start.time <- proc.time()[3]
        timing$Model <- (proc.time()[3] - start.time)
        start.time <- proc.time()[3]
    } else {
        lavaanStart <- 
            StartingValues(start.method = start,
                           user         = lavaanUser, 
                           sample       = lavaanSampleStats,
                           model.type   = lavaanOptions$model.type,
                           mimic        = lavaanOptions$mimic,
                           debug        = lavaanOptions$debug)
        timing$Start <- (proc.time()[3] - start.time)
        start.time <- proc.time()[3]

        # 5. construct internal model (S4) representation
        lavaanModel <- 
            Model(user           = lavaanUser, 
                  start          = lavaanStart, 
                  representation = lavaanOptions$representation,
                  debug          = lavaanOptions$debug)
    }

    # 6. estimate free parameters
    x <- NULL
    if(do.fit && lavaanModel@nx.free > 0L) {
        x <- estimateModel(lavaanModel,
                           sample  = lavaanSampleStats,
                           options = lavaanOptions,
                           control = control)
        lavaanModel <- setModelParameters(lavaanModel, x = x)
        if(!is.null(attr(x, "con.jac"))) 
            lavaanModel@con.jac <- attr(x, "con.jac")
        # check if model has converged or not
        if(!attr(x, "converged") && lavaanOptions$warn) {
           warning("lavaan WARNING: model has NOT converged!")
        }
    } else {
        x <- numeric(0L)
        attr(x, "iterations") <- 0L; attr(x, "converged") <- FALSE
        attr(x, "control") <- control
        attr(x, "fx") <- 
            computeObjective(lavaanModel, sample = lavaanSampleStats, 
                             estimator = lavaanOptions$estimator)
    }
    timing$Estimate <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]


    # 7. estimate vcov of free parameters (for standard errors)
    VCOV <- NULL
    if(lavaanOptions$se != "none" && lavaanModel@nx.free > 0L) {
        VCOV <- estimateVCOV(lavaanModel,
                             sample  = lavaanSampleStats,
                             options = lavaanOptions,
                             data    = lavaanData)
    }
    timing$VCOV <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 8. compute test statistic (chi-square and friends)
    TEST <- NULL
    if(lavaanOptions$test != "none") {
        TEST <- computeTestStatistic(lavaanModel,
                                     user    = lavaanUser,
                                     sample  = lavaanSampleStats,
                                     options = lavaanOptions,
                                     x       = x,
                                     VCOV    = VCOV,
                                     data    = lavaanData)
    }
    timing$TEST <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 9. collect information about model fit (S4)
    lavaanFit <- Fit(user  = lavaanUser, 
                     start = lavaanStart, 
                     model = lavaanModel,
                     x     = x, 
                     VCOV  = VCOV,
                     TEST  = TEST)
    timing$total <- (proc.time()[3] - start.time0)

    # 10. construct lavaan object
    lavaan <- new("lavaan",
                  call    = mc,                     # match.call
                  timing  = timing,                 # list
                  Options = lavaanOptions,          # list
                  User    = lavaanUser,             # list
                  Data    = lavaanData,             # list
                  Sample  = lavaanSampleStats,      # S4 class
                  Model   = lavaanModel,            # S4 class
                  Fit     = lavaanFit               # S4 class
                 )

    lavaan
}

# cfa + sem
cfa <- sem <- function(model = NULL,
    meanstructure = "default", fixed.x = "default",
    orthogonal = FALSE, std.lv = FALSE, data = NULL, std.ov = FALSE,
    missing = "default", sample.cov = NULL, sample.mean = NULL,
    sample.nobs = NULL, group = NULL, group.equal = "",
    group.partial = "", constraints = "",
    estimator = "default", likelihood = "default",
    information = "default", se = "default", test = "default",
    bootstrap = 1000L, mimic = "default", representation = "default",
    do.fit = TRUE, control = list(), start = "default", 
    slotOptions = NULL, slotUser = NULL, slotSample = NULL,
    slotData = NULL, slotModel = NULL,
    verbose = FALSE, warn = TRUE, debug = FALSE) {

    mc <- match.call()

    mc$model.type      = as.character( mc[[1L]] )
    if(length(mc$model.type) == 3L) mc$model.type <- mc$model.type[3L]
    mc$int.ov.free     = TRUE
    mc$int.lv.free     = FALSE
    mc$auto.fix.first  = !std.lv
    mc$auto.fix.single = TRUE
    mc$auto.var        = TRUE
    mc$auto.cov.lv.x   = TRUE
    mc$auto.cov.y      = TRUE
    mc[[1L]] <- as.name("lavaan")

    eval(mc, parent.frame())
}

# simple growth models
growth <- function(model = NULL,
    fixed.x = "default",
    orthogonal = FALSE, std.lv = FALSE, data = NULL, std.ov = FALSE,
    missing = "default", sample.cov = NULL, sample.mean = NULL,
    sample.nobs = NULL, group = NULL, group.equal = "",
    group.partial = "", constraints = "",
    estimator = "default", likelihood = "default",
    information = "default", se = "default", test = "default",
    bootstrap = 1000L, mimic = "default", representation = "default",
    do.fit = TRUE, control = list(), start = "default",
    slotOptions = NULL, slotUser = NULL, slotSample = NULL,
    slotData = NULL, slotModel = NULL,
    verbose = FALSE, warn = TRUE, debug = FALSE) {

    mc <- match.call()

    mc$model.type      = "growth"
    mc$int.ov.free     = FALSE
    mc$int.lv.free     = TRUE
    mc$auto.fix.first  = !std.lv
    mc$auto.fix.single = TRUE
    mc$auto.var        = TRUE
    mc$auto.cov.lv.x   = TRUE
    mc$auto.cov.y      = TRUE
    mc[[1L]] <- as.name("lavaan")

    eval(mc, parent.frame())
}
