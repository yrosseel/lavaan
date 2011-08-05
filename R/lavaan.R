# main user-visi ble cfa/sem/growth functions 
#
# initial version: YR 25/03/2009
# added lavaanOptions YR 02/08/2010
# major revision: YR 9/12/2010: - new workflow (since 0.4-5)
#                               - merge cfa/sem/growth functions

lavaan <- function(# user-specified model syntax
                   model.syntax    = '',
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
                   mimic           = "default",
                   representation  = "default",
                   do.fit          = TRUE,

                   # starting values
                   start           = "default",
  
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
    if(!is.null(data)) {
        # good, we got a full data frame
        stopifnot(is.data.frame(data))
        if(!is.null(group)) {
            if(!(group %in% names(data))) {
                stop("grouping variable `", group,
                     "' not found in names data:", names(data))
            }
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
    } else {
        # sample.cov is null?
        if(is.null(sample.cov))
            stop("both data and sample.cov are null")

        # we also need the number of observations (per group)
        if(is.null(sample.nobs))
            stop("please specify number of observations")

        # if meanstructure=TRUE, we need sample.mean
        if(meanstructure == TRUE && is.null(sample.mean))
            stop("please provide sample.mean if meanstructure=TRUE")

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
                stop("sample.cov must be a matrix or a list of matrices")
            ngroups <- 1L
            group.label <- character(0)
        }
        data.type <- "moment"
    }


    # 1a. collect various options/flags and fill in `default' values
    #opt <- modifyList(formals(lavaan), as.list(mc)[-1])
    # force evaluation of `language` and/or `symbol` arguments
    #opt <- lapply(opt, function(x) if(typeof(x) %in% c("language", "symbol")) 
    #                                   eval(x, parent.frame()) else x)
    opt <- list(model.syntax = model.syntax, model.type = model.type,
        meanstructure = meanstructure, int.ov.free = int.ov.free,
        int.lv.free = int.lv.free, fixed.x = fixed.x, orthogonal = orthogonal,
        std.lv = std.lv, auto.fix.first = auto.fix.first,
        auto.fix.single = auto.fix.single, auto.var = auto.var,
        auto.cov.lv.x = auto.cov.lv.x, auto.cov.y = auto.cov.y,
        std.ov = std.ov, missing = missing, group = group, 
        group.equal = group.equal, group.partial = group.partial, 
        constraints = constraints,
        estimator = estimator, likelihood = likelihood,
        information = information, se = se, test = test, mimic = mimic,
        representation = representation, do.fit = do.fit, verbose = verbose,
        warn = warn, debug = debug, data.type = data.type)
    lavaanOptions <- setLavaanOptions(opt)
    timing$InitOptions <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
    

    # 2a. construct lavaan User list: description of the user-specified model
    if(is.character(model.syntax)) {
        lavaanUser <- 
            lavaanify(model.syntax    = model.syntax, 
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
    } else if(is.list(model.syntax)) {
        # two possibilities: either model.syntax is already lavaanified
        # or it is a list of model.syntaxes (perhaps one for each group)
        if(!is.null(model.syntax$lhs) &&
           !is.null(model.syntax$op)  &&
           !is.null(model.syntax$rhs) &&
           !is.null(model.syntax$free)) {
            lavaanUser <- model.syntax
        } else if(is.character(model.syntax[[1]])) {
            # we lavaanify each model in term, and assume multiple groups
            stop("lavaan ERROR: list of model syntaxes: not implemented yet")
        }
    }

    # 2b. change meanstructure flag?
    if(any(lavaanUser$op == "~1")) lavaanOptions$meanstructure <- TRUE
    timing$User <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 3. construct lavaan Sample (S4) object: description of the data
    lavaanSample <- 
        Sample(data          = data, 
               group         = group,
               sample.cov    = sample.cov, 
               sample.mean   = sample.mean,
               sample.nobs   = sample.nobs, 
               std.ov        = std.ov,
               
               ov.names      = vnames(lavaanUser, "ov"),
               data.type     = data.type,
               ngroups       = ngroups,
               group.label   = group.label,
               estimator     = lavaanOptions$estimator,
               likelihood    = lavaanOptions$likelihood,
               mimic         = lavaanOptions$mimic,
               meanstructure = lavaanOptions$meanstructure,
               missing       = lavaanOptions$missing,
               warn          = lavaanOptions$warn,
               verbose       = lavaanOptions$verbose)
    timing$Sample <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 4. compute some reasonable starting values 
    lavaanStart <- 
        StartingValues(start.method = start,
                       user         = lavaanUser, 
                       sample       = lavaanSample, 
                       model.type   = lavaanOptions$model.type,
                       debug        = lavaanOptions$debug)
    timing$Start <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 5. construct internal model (S4) representation
    lavaanModel <- 
        Model(user           = lavaanUser, 
              start          = lavaanStart, 
              representation = lavaanOptions$representation,
              debug          = lavaanOptions$debug)
    timing$Model <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 6. estimate free parameters
    x <- NULL
    if(do.fit && lavaanModel@nx.free > 0L) {
        x <- estimateModel(lavaanModel,
                           sample  = lavaanSample, 
                           options = lavaanOptions)
        lavaanModel <- setModelParameters(lavaanModel, x = x)
        if(!is.null(attr(x, "con.jac"))) 
            lavaanModel@con.jac <- attr(x, "con.jac")
    } else {
        x <- numeric(0L)
        attr(x, "iterations") <- 0L; attr(x, "converged") <- FALSE
        attr(x, "fx") <- computeObjective(lavaanModel, sample = lavaanSample, 
                                          estimator = lavaanOptions$estimator)
    }
    timing$Estimate <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 7. estimate vcov of free parameters (for standard errors)
    VCOV <- NULL
    if(opt$se != "none" && lavaanModel@nx.free > 0L) {
        VCOV <- estimateVCOV(lavaanModel,
                             sample  = lavaanSample,
                             options = lavaanOptions)
    }
    timing$VCOV <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 8. compute test statistic (chi-square and friends)
    TEST <- NULL
    if(opt$test != "none") {
        TEST <- computeTestStatistic(lavaanModel,
                                     user    = lavaanUser,
                                     sample  = lavaanSample,
                                     options = lavaanOptions,
                                     x       = x,
                                     VCOV    = VCOV)
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
                  call    = mc,             # match.call
                  timing  = timing,         # list
                  Options = lavaanOptions,  # list
                  User    = lavaanUser,     # list
                  Sample  = lavaanSample,   # S4 class
                  Model   = lavaanModel,    # S4 class
                  Fit     = lavaanFit       # S4 class
                 )

    lavaan
}

# cfa + sem
cfa <- sem <- function(model.syntax = '', 
    meanstructure = "default", fixed.x = "default",
    orthogonal = FALSE, std.lv = FALSE, data = NULL, std.ov = FALSE,
    missing = "default", sample.cov = NULL, sample.mean = NULL,
    sample.nobs = NULL, group = NULL, group.equal = "",
    group.partial = "", constraints = "",
    estimator = "default", likelihood = "default",
    information = "default", se = "default", test = "default",
    mimic = "default", representation = "default",
    do.fit = TRUE, start = "default", 
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
growth <- function(model.syntax = '',
    fixed.x = "default",
    orthogonal = FALSE, std.lv = FALSE, data = NULL, std.ov = FALSE,
    missing = "default", sample.cov = NULL, sample.mean = NULL,
    sample.nobs = NULL, group = NULL, group.equal = "",
    group.partial = "", constraints = "",
    estimator = "default", likelihood = "default",
    information = "default", se = "default", test = "default",
    mimic = "default", representation = "default",
    do.fit = TRUE, start = "default",
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
