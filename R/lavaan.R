# main user-visi ble cfa/sem/growth functions 
#
# initial version: YR 25/03/2009
# added lavaanOptions YR 02/08/2010
# major revision: YR 9/12/2010: - new workflow (since 0.4-5)
#                               - merge cfa/sem/growth functions
# YR 25/02/2012: changed data slot (from list() to S4); data@X contains data

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
    # start timer
    start.time0 <- start.time <- proc.time()[3]; timing <- list()

    # 0a. store call
    mc  <- match.call()

    # 0b. get ov.names (per group) -- needed for getData()
    if(!is.null(slotUser)) {
        FLAT <- slotUser
    } else if(is.character(model)) {
        FLAT <- parseModelString(model)
    } else if(is.list(model)) {
        FLAT <- model
    }
    if(max(FLAT$group) < 2L) { # same model for all groups 
        ov.names <- vnames(FLAT, type="ov")
    } else { # different model per group
        ov.names <- lapply(1:max(FLAT$group),
                           function(x) vnames(FLAT, type="ov", group=x))
    }

    # 0c. get data.type
    if(!is.null(slotData)) {
        data.type = "slotData"
    } else if(!is.null(data)) {
        data.type = "full"
    } else if(!is.null(sample.cov)) {
        data.type = "moment"
    } else {
        data.type = "none"
        # no data? no fitting!
        do.fit <- FALSE; se <- "none"; test <- "none"; start <- "simple"
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

    # 1b. check data/sample.cov and get the number of groups
    if(data.type == "slotData") {
        lavaanData <- slotData
    } else if(data.type == "full") {
        stopifnot(is.data.frame(data)) ## FIXME!! we should also allow matrices
        lavaanData <- getData(data     = data,
                              group    = group,
                              ov.names = ov.names,
                              std.ov   = lavaanOptions$std.ov,
                              missing  = lavaanOptions$missing,
                              warn     = lavaanOptions$warn)
    } else if(data.type == "moment") {
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
            ngroups <- 1L; group.label <- character(0)
            if(!is.matrix(sample.cov))
                stop("lavaan ERROR: sample.cov must be a matrix or a list of matrices")
        }
        if(!is.list(ov.names)) {
            tmp <- ov.names; ov.names <- vector("list", length=ngroups)
            ov.names[1:ngroups] <- list(tmp)
        } else {
            if(length(ov.names) != ngroups)
                stop("lavaan ERROR: ov.names assumes ", length(ov.names),
                     " groups; data contains ", ngroups, " groups")
        }
        lavaanData <- new("lavaanData",
                          ngroups = ngroups, group.label = group.label,
                          nobs = as.list(sample.nobs),
                          norig = as.list(sample.nobs),
                          ov.names = ov.names)
    } else { # both data and sample.cov are NULL; simulating?
        if(is.null(sample.nobs)) sample.nobs <- 0L
        sample.nobs <- as.list(sample.nobs)
        ngroups <- length(unlist(sample.nobs))
        if(ngroups > 1L) 
            group.label <- paste("Group ", 1:ngroups, sep="")
        else
            group.label <- character(0)
        if(!is.list(ov.names)) {
            tmp <- ov.names; ov.names <- vector("list", length=ngroups)
            ov.names[1:ngroups] <- list(tmp)
        }
        lavaanData <- new("lavaanData",
                          ngroups = ngroups, group.label = group.label,
                          nobs = sample.nobs,
                          norig = sample.nobs,
                          ov.names = ov.names)
    }
    timing$InitData <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
    if(debug) print(str(lavaanData))


    # 2a. construct lavaan User list: description of the user-specified model
    if(!is.null(slotUser)) {
        lavaanUser <- slotUser
    } else if(is.character(model)) {
        lavaanUser <- 
            lavaanify(model           = FLAT,
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

                      ngroups         = lavaanData@ngroups,
                      group.equal     = lavaanOptions$group.equal, 
                      group.partial   = lavaanOptions$group.partial,
                      debug           = lavaanOptions$debug,
                      warn            = lavaanOptions$warn,

                      as.data.frame.  = FALSE)
    } else if(is.list(model)) {
        # two possibilities: either model is already lavaanified
        # or it is something else...
        if(!is.null(model$lhs) && !is.null(model$op)  &&
           !is.null(model$rhs) && !is.null(model$free)) {
            lavaanUser <- model
        } else if(is.character(model[[1]])) {
            stop("lavaan ERROR: model is a list, but not a parameterTable?")
        }
    } else {
        cat("model type: ", class(model), "\n")
        stop("lavaan ERROR: model is not of type character or list")
    }

    # 2b. change meanstructure flag?
    if(any(lavaanUser$op == "~1")) lavaanOptions$meanstructure <- TRUE

    # 2c. prepare constraints functions
    timing$User <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 3. get sample statistics
    # here we know the number of groups!
    ov.names <- lapply(as.list(1:lavaanData@ngroups),
                       function(x) vnames(lavaanUser, type="ov", x))

    if(!is.null(slotSample)) {
        lavaanSampleStats <- slotSample
    } else if(data.type == "full") {
        WLS.V <- list()
        if(lavaanOptions$estimator %in% c("GLS", "WLS")) {
            WLS.V <- getWLS.V(Data          = lavaanData,
                              sample        = NULL,
                              estimator     = lavaanOptions$estimator,
                              mimic         = lavaanOptions$mimic,
                              meanstructure = lavaanOptions$meanstructure)
        }

        lavaanSampleStats <- getSampleStatsFromData(
                           Data        = lavaanData,
                           rescale     = (lavaanOptions$estimator == "ML" &&
                                          lavaanOptions$likelihood == "normal"),
                           WLS.V       = WLS.V,
                           verbose     = lavaanOptions$verbose)
                                                 
    } else if(data.type == "moment") {
        lavaanSampleStats <- getSampleStatsFromMoments(
                           sample.cov  = sample.cov,
                           sample.mean = sample.mean,
                           sample.nobs = sample.nobs,
                           ov.names    = ov.names,
                           rescale     = (lavaanOptions$estimator == "ML" &&
                                          lavaanOptions$likelihood == "normal"))

        if(lavaanOptions$estimator == "GLS") {
            WLS.V <- getWLS.V(Data          = NULL,
                              sample        = lavaanSampleStats,
                              estimator     = lavaanOptions$estimator,
                              mimic         = lavaanOptions$mimic,
                              meanstructure = lavaanOptions$meanstructure)
            lavaanSampleStats@WLS.V <- WLS.V
        }
    } else {
        # no data
        lavaanSampleStats <- new("SampleStats", ngroups=lavaanData@ngroups,
                                 nobs=as.list(rep(0L, lavaanData@ngroups)),
                                 missing.flag=FALSE)
    } 
    timing$Sample <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
    if(debug) print(str(lavaanSampleStats))

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
        timing$Model <- (proc.time()[3] - start.time)
        start.time <- proc.time()[3]
    }

    # 6. estimate free parameters
    x <- NULL
    if(do.fit && lavaanModel@nx.free > 0L) {
        # catch simple linear regression models
        if(length(unique(lavaanUser$lhs[lavaanUser$op == "~"])) == 1L && 
           length(vnames(lavaanUser,   "lv")) == 0L &&
           length(lavaanData@X) > 0L &&
           lavaanData@ngroups == 1L &&
           lavaanOptions$fixed &&
           lavaanOptions$missing == "listwise") {
            # simple univariate regression
            ov.y.idx <- match(vnames(lavaanUser, "ov.y"), 
                              colnames(lavaanData@X[[1L]]))
            ov.x.idx <- match(vnames(lavaanUser, "ov.x"), 
                              colnames(lavaanData@X[[1L]]))
            YX <- lavaanData@X[[1L]]
            #print(head(YX))
            # constraints?
            if(sum(length(lavaanModel@x.ceq.idx) + 
                   length(lavaanModel@x.cin.idx)) == 0L) {
                out <- lm.fit(x=cbind(1,YX[,ov.x.idx]), 
                              y=YX[,ov.y.idx])
                x.beta <- out$coefficients
                y.rvar <- sum(out$residuals^2)/length(out$residuals) #ML?
                if(!lavaanOptions$meanstructure) {
                    x.beta <- x.beta[-1L]
                    x <- c(x.beta, y.rvar)
                } else {
                    x <- c(x.beta[-1L], y.rvar, x.beta[1L])
                }
                lavaanModel <- setModelParameters(lavaanModel, x = x)
                attr(x, "iterations") <- 1L; attr(x, "converged") <- TRUE
                attr(x, "control") <- control
                attr(x, "fx") <-
                computeObjective(lavaanModel, sample = lavaanSampleStats,
                                 estimator = lavaanOptions$estimator)
            } else if(checkLinearConstraints(lavaanModel) == TRUE) {
                require(quadprog)

                A.ceq <- A.cin <- matrix(0, lavaanModel@nx.free, 0)
                if(!is.null(body(lavaanModel@ceq.function)))
                    A.ceq <- t(jacobian(func=lavaanModel@ceq.function, 
                                        x=rep(0,lavaanModel@nx.free)))
                if(!is.null(body(lavaanModel@cin.function)))
                    A.cin <- t(jacobian(func=lavaanModel@cin.function, 
                                        x=rep(0,lavaanModel@nx.free)))
                A <- cbind(A.ceq, A.cin)
                # meanstructure? last row is intercept
                if(lavaanOptions$meanstructure) {
                    A <- rbind(A[nrow(A),,drop=FALSE], A[-nrow(A),,drop=FALSE])
                } else {
                    A <- rbind(rep(0,ncol(A)), A)
                }
                # remove residual variance (last row)
                A <- A[-nrow(A),,drop=FALSE]
                #A <- matrix( c( 0, -1, 1,  0,
                #                0, 0,  -1, 1), 4, 2)
                X <- cbind(1,YX[,ov.x.idx]); X.X <- crossprod(X)
                Y <- YX[,ov.y.idx]; X.Y <- crossprod(X, Y)
                out <- solve.QP(Dmat=X.X, dvec=X.Y, Amat=A, 
                                bvec=rep(0, NCOL(A)), 
                                meq=length(lavaanModel@x.ceq.idx))
                x.beta <- out$solution
                residuals <- Y - (X %*% x.beta)
                y.rvar <- sum(residuals^2)/length(residuals) #ML?
                if(!lavaanOptions$meanstructure) {
                    x.beta <- x.beta[-1L]
                } else {
                    x.beta <- c(x.beta[-1L], x.beta[1L])
                }
                x <- c(x.beta, y.rvar)
                lavaanModel <- setModelParameters(lavaanModel, x = x)
                attr(x, "iterations") <- 1L; attr(x, "converged") <- TRUE
                attr(x, "control") <- control
                attr(x, "fx") <-
                    computeObjective(lavaanModel, sample = lavaanSampleStats,
                                     estimator = lavaanOptions$estimator)
            } else {
                # regular estimation after all
                x <- estimateModel(lavaanModel,
                               sample  = lavaanSampleStats,
                               options = lavaanOptions,
                               control = control)
                lavaanModel <- setModelParameters(lavaanModel, x = x)
            }
        } else {
            #cat("REGULAR\n")
            # regular estimation
            x <- estimateModel(lavaanModel,
                               sample  = lavaanSampleStats,
                               options = lavaanOptions,
                               control = control)
            lavaanModel <- setModelParameters(lavaanModel, x = x)
        }
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
                             data    = lavaanData,
                             user    = lavaanUser,
                             control = control)
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
                                     data    = lavaanData,
                                     control = control)
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
