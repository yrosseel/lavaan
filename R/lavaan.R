# main user-visi ble cfa/sem/growth functions 
#
# initial version: YR 25/03/2009
# added lavaanOptions YR 02/08/2010
# major revision: YR 9/12/2010: - new workflow (since 0.4-5)
#                               - merge cfa/sem/growth functions
# YR 25/02/2012: changed data slot (from list() to S4); data@X contains data

lavaan <- function(# user-specified model: can be syntax, parameter Table, ...
                   model              = NULL,
                   data               = NULL,  # second argument, most used!
                   model.type         = "sem",
                
                   # model modifiers
                   meanstructure      = "default",
                   int.ov.free        = FALSE,
                   int.lv.free        = FALSE,
                   fixed.x            = "default", # or FALSE?
                   orthogonal         = FALSE,
                   std.lv             = FALSE,

                   auto.fix.first     = FALSE,
                   auto.fix.single    = FALSE,
                   auto.var           = FALSE,
                   auto.cov.lv.x      = FALSE,
                   auto.cov.y         = FALSE,
                   auto.th            = FALSE,
                   auto.delta         = FALSE,
                   
                   # full data
                   std.ov             = FALSE,
                   missing            = "default",
                   ordered            = NULL,

                   # summary data
                   sample.cov         = NULL,
                   sample.cov.rescale = "default",
                   sample.mean        = NULL,
                   sample.nobs        = NULL,

                   # multiple groups
                   group              = NULL,
                   group.label        = NULL,
                   group.equal        = '',
                   group.partial      = '',

                   # clusters
                   cluster            = NULL,
             
                   # constraints
                   constraints        = '',

                   # estimation
                   estimator          = "default",
                   likelihood         = "default",
                   information        = "default",
                   se                 = "default",
                   test               = "default",
                   bootstrap          = 1000L,
                   mimic              = "default",
                   representation     = "default",
                   do.fit             = TRUE,
                   control            = list(),
                   WLS.V              = NULL,
                   NACOV              = NULL,

                   # starting values
                   start              = "default",

                   # full slots from previous fits
                   slotOptions        = NULL,
                   slotParTable       = NULL,
                   slotSampleStats    = NULL,
                   slotData           = NULL,
                   slotModel          = NULL,
  
                   # verbosity
                   verbose            = FALSE,
                   warn               = TRUE,
                   debug              = FALSE
                  )
{

    # temporary block estimator = "PML"
    # if(estimator == "PML") stop("estimator PML is not available yet")

    # start timer
    start.time0 <- start.time <- proc.time()[3]; timing <- list()

    # 0a. store call
    mc  <- match.call()

    # 0b. get ov.names and ov.names.x (per group) -- needed for lavData()
    if(!is.null(slotParTable)) {
        FLAT <- slotParTable
    } else if(is.character(model)) {
        FLAT <- parseModelString(model)
    } else if(is.list(model)) {
        FLAT <- model
    }
    if(max(FLAT$group) < 2L) { # same model for all groups 
        ov.names   <- vnames(FLAT, type="ov")
    } else { # different model per group
        ov.names <- lapply(1:max(FLAT$group),
                           function(x) vnames(FLAT, type="ov", group=x))
    }

    # 0c categorical variables? -- needed for lavaanOptions
    categorical <- any(FLAT$op == "|")
    #if(!is.null(data) && length(ordered) > 0L) {
    #    categorical <- TRUE
    #}
    if(!is.null(data)) {
        if(length(ordered) > 0L) {
            # check 'ordered'
            not.in.ov <- which(!ordered %in% unlist(ov.names))
            if(length(not.in.ov) > 0L) {
                stop("ordered argument contains variable name(s) not in model: ",                     paste(ordered[not.in.ov], collapse=" "))
            }
            data[,ordered] <- lapply(data[,ordered,drop=FALSE], base::ordered)
            #
            # NOTE: we coerce these variables to 'ordered' here in
            # the 'data' data.frame; however, (at least in 2.15.1), this
            # creates (at least) 3 copies of the full data.frame...
            # 
            # we should try not to 'touch' the data.frame at all (giving us
            # a lot of housekeeping work in lavData) ... TODO!
            ov.types <- sapply(data[,unlist(ov.names)],
                               function(x) class(x)[1])
            categorical <- TRUE
        } else {
            ov.types <- sapply(data[,unlist(ov.names)], 
                               function(x) class(x)[1])
            if("ordered" %in% ov.types) 
                categorical <- TRUE
        }        
    }

    # if categorical, make a distinction between exo and the rest
    if(categorical) {
        if(max(FLAT$group) < 2L) { # same model for all groups 
            ov.names   <- vnames(FLAT, type="ov.nox")
            ov.names.x <- vnames(FLAT, type="ov.x")
        } else { # different model per group
            ov.names <- lapply(1:max(FLAT$group),
                               function(x) vnames(FLAT, type="ov.nox", group=x))
            ov.names.x <- lapply(1:max(FLAT$group),
                                 function(x) vnames(FLAT, type="ov.x", group=x))
        }
    } else {
        ov.names.x <- character(0)
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
            auto.cov.y = auto.cov.y, auto.th = auto.th, 
            auto.delta = auto.delta, missing = missing, 
            group = group, categorical = categorical,
            group.equal = group.equal, group.partial = group.partial, 
            constraints = constraints,
            estimator = estimator, likelihood = likelihood, 
            sample.cov.rescale = sample.cov.rescale,
            information = information, se = se, test = test, 
            bootstrap = bootstrap, mimic = mimic,
            representation = representation, do.fit = do.fit, verbose = verbose,
            warn = warn, debug = debug)
        lavaanOptions <- setLavaanOptions(opt)
    }
    timing$InitOptions <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # some additional checks for estimator="PML"
    if(lavaanOptions$estimator == "PML") {
        # 0. at least some variables must be ordinal
        if(!any(ov.types == "ordered")) {
            stop("lavaan ERROR: estimator=\"PML\" is only available if some variables are ordinal")
        }
        # 1. all variables must be ordinal (for now)
        #    (the mixed continuous/ordinal case will be added later)
        if(any(ov.types != "ordered")) {
            stop("lavaan ERROR: estimator=\"PML\" can not handle mixed continuous and ordinal data (yet)")
        }
        
        # 2. we can not handle exogenous covariates yet
        if(length(ov.names.x) > 0L) {
            stop("lavaan ERROR: estimator=\"PML\" can not handle exogenous covariates (yet)")
        }

        # 3. warn that this is still experimental
        message("Please note: the PML estimator is still under development.\n",
                "Future releases will improve speed, and allow for mixed ord/cont variables.\n",
                "Research on how to compute a proper goodness-of-fit test is ongoing.\n")
    }

    # 1b. check data/sample.cov and get the number of groups
    if(!is.null(slotData)) {
        lavaanData <- slotData
    } else {
        lavaanData <- lavData(data        = data,
                              group       = group,
                              group.label = group.label,
                              ov.names    = ov.names,
                              ordered     = ordered,
                              ov.names.x  = ov.names.x,
                              std.ov      = std.ov,
                              missing     = lavaanOptions$missing,
                              sample.cov  = sample.cov,
                              sample.mean = sample.mean,
                              sample.nobs = sample.nobs,
                              warn        = lavaanOptions$warn)

        # what have we learned from the data?
        if("ordered" %in% lavaanData@ov$type) {
            if(lavaanOptions$estimator == "ML")
                stop("lavaan ERROR: estimator ML for ordered data is not supported yet. Use WLSMV instead.")
            # Mplus style
            lavaanOptions$meanstructure <- TRUE
        }
    }
    if(lavaanData@data.type == "none") {
        do.fit <- FALSE; start <- "simple"
        lavaanOptions$se <- "none"; lavaanOptions$test <- "none"
    } else if(lavaanData@data.type == "moment") {
        # catch here some options that will not work with moments
        if(estimator %in% c("MLM", "MLMV", "MLR", "ULSM", "ULSMV") &&
           is.null(NACOV)) {
            stop("lavaan ERROR: estimator ", estimator, " requires full data or user-provided NACOV")
        }
        if(estimator %in% c("WLS", "WLSM", "WLSMV", "DWLS") &&
           is.null(WLS.V)) {
            stop("lavaan ERROR: estimator ", estimator, " requires full data or user-provided WLS.V") 
        }
    }
    timing$InitData <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
    if(debug) print(str(lavaanData))


    # 2a. construct paramter table: description of the user-specified model
    if(!is.null(slotParTable)) {
        lavaanParTable <- slotParTable
    } else if(is.character(model)) {

        # check FLAT before we proceed
        if(debug) print(as.data.frame(FLAT))

        # check 1: catch ~~ of fixed.x covariates if categorical
        if(lavaanOptions$categorical) {
            tmp <- vnames(FLAT, type="ov.x", ov.x.fatal=TRUE)
        }
        # check 2: catch ~1 of ordered variables (unless they are fixed)
        if(lavaanOptions$categorical) {
            int.idx <- which(FLAT$op == "~1" & FLAT$fixed == "")
            if(length(int.idx) > 0L) {
                INT <- FLAT$lhs[int.idx]
                ORD <- lavaanData@ov$name[ lavaanData@ov$type == "ordered" ]
                if(any(INT %in% ORD))        
                    stop("lavaan ERROR: model syntax contains free intercepts for ordinal dependent variable(s): [", paste(INT, collapse=" "), 
                         "];\n  Please remove them and try again.")
            }
        }

        lavaanParTable <- 
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
                      auto.th         = lavaanOptions$auto.th,
                      auto.delta      = lavaanOptions$auto.delta,

                      varTable        = lavaanData@ov,
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
            lavaanParTable <- model
        } else if(is.character(model[[1]])) {
            stop("lavaan ERROR: model is a list, but not a parameterTable?")
        }
    } else {
        cat("model type: ", class(model), "\n")
        stop("lavaan ERROR: model is not of type character or list")
    }

    # 2b. change meanstructure flag?
    if(any(lavaanParTable$op == "~1")) lavaanOptions$meanstructure <- TRUE

    # 2c. prepare constraints functions
    timing$ParTable <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 3. get sample statistics
    # here we know the number of groups!
    ov.names <- lapply(as.list(1:lavaanData@ngroups),
                       function(x) vnames(lavaanParTable, type="ov", x))

    if(!is.null(slotSampleStats)) {
        lavaanSampleStats <- slotSampleStats
    } else if(lavaanData@data.type == "full") {
        lavaanSampleStats <- lavSampleStatsFromData(
                       Data          = lavaanData,
                       rescale       = (lavaanOptions$estimator == "ML" &&
                                        lavaanOptions$likelihood == "normal"),
                       estimator     = lavaanOptions$estimator,
                       mimic         = lavaanOptions$mimic,
                       meanstructure = lavaanOptions$meanstructure,
                       missing.h1    = (lavaanOptions$missing != "listwise"),
                       WLS.V         = WLS.V,
                       NACOV         = NACOV,
                       verbose       = lavaanOptions$verbose)
                                                 
    } else if(lavaanData@data.type == "moment") {
        lavaanSampleStats <- lavSampleStatsFromMoments(
                           sample.cov    = sample.cov,
                           sample.mean   = sample.mean,
                           sample.nobs   = sample.nobs,
                           ov.names      = ov.names,
                           estimator     = lavaanOptions$estimator,
                           mimic         = lavaanOptions$mimic,
                           meanstructure = lavaanOptions$meanstructure,
                           WLS.V         = WLS.V,
                           NACOV         = NACOV,
                           rescale       = lavaanOptions$sample.cov.rescale)
                           
    } else {
        # no data
        lavaanSampleStats <- new("lavSampleStats", ngroups=lavaanData@ngroups,
                                 nobs=as.list(rep(0L, lavaanData@ngroups)),
                                 cov.x=vector("list", length=lavaanData@ngroups),
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
                           partable     = lavaanParTable, 
                           samplestats  = lavaanSampleStats,
                           model.type   = lavaanOptions$model.type,
                           mimic        = lavaanOptions$mimic,
                           debug        = lavaanOptions$debug)
        timing$Start <- (proc.time()[3] - start.time)
        start.time <- proc.time()[3]

        #lavaanParTable$start <- lavaanStart
        #print(as.data.frame(lavaanParTable))

        # 5. construct internal model (S4) representation
        lavaanModel <- 
            Model(partable       = lavaanParTable, 
                  start          = lavaanStart, 
                  representation = lavaanOptions$representation,
                  th.idx         = lavaanSampleStats@th.idx,
                  debug          = lavaanOptions$debug)
        timing$Model <- (proc.time()[3] - start.time)
        start.time <- proc.time()[3]
    }

    # prepare cache -- stuff needed for estimation, but also post-estimation
    lavaanCache <- vector("list", length=lavaanData@ngroups)
    if(estimator == "PML") {
        TH <- computeTH(lavaanModel)
        for(g in 1:lavaanData@ngroups) {
            nvar <- ncol(lavaanData@X[[g]])
            th.idx <- lavaanModel@th.idx[[g]]
            # pairwise tables, as a long vector
            PW <- pairwiseTables(data=lavaanData@X[[g]], no.x=nvar)$pairTables
            bifreq <- as.numeric(unlist(PW))
            ### FIXME!!! Check for zero cells!!
            #zero.idx <- which(bifreq == 0)
            #bifreq[zero.idx] <- 0.001 ####????
            LONG <- LongVecInd(no.x               = nvar,
                               all.thres          = TH[[g]],
                               index.var.of.thres = th.idx)
            lavaanCache[[g]] <- list(bifreq=bifreq, LONG=LONG)
        }
    }


    # 6. estimate free parameters
    x <- NULL
    if(do.fit && lavaanModel@nx.free > 0L) {
        # catch simple linear regression models
        if(length(unique(lavaanParTable$lhs[lavaanParTable$op == "~"])) == 1L && 
           length(vnames(lavaanParTable,   "lv")) == 0L &&
           #FALSE && # to debug
           sum(nchar(FLAT$fixed)) == 0 && # no fixed values in parTable
                                          # this includes intercepts!!
           !categorical &&
           !lavaanModel@eq.constraints &&
           length(lavaanData@X) > 0L &&
           lavaanData@ngroups == 1L &&
           lavaanOptions$fixed &&
           lavaanOptions$missing == "listwise") {
            # simple univariate regression
            ov.y.idx <- match(vnames(lavaanParTable, "ov.y"), 
                              lavaanData@ov.names[[1L]])
            ov.x.idx <- match(vnames(lavaanParTable, "ov.x"), 
                              lavaanData@ov.names[[1L]])
            YX <- lavaanData@X[[1L]]
            #print(head(YX))
            # constraints?
            if(sum(length(lavaanModel@x.ceq.idx) + 
                   length(lavaanModel@x.cin.idx)) == 0L) {

                ## forced zero intercept?
                #yvar <- vnames(lavaanParTable, "ov.y")
                #int.y.idx <- which(lavaanParTable$lhs == yvar &
                #                   lavaanParTable$op == "~1")
                #if(length(int.y.idx) > 0L && 
                #   !is.na(lavaanParTable$ustart[int.y.idx]) &&
                #   lavaanParTable$ustart[int.y.idx] == 0) {
                #    lm.intercept <- FALSE
                #} else {
                #    lm.intercept <- TRUE
                #}
                #if(lm.intercept) {
                    out <- lm.fit(x=cbind(1,YX[,ov.x.idx]), 
                                  y=YX[,ov.y.idx])
                    x.beta <- out$coefficients
                #} else {
                #     out <- lm.fit(x=YX[,ov.x.idx], y=YX[,ov.y.idx])
                #     x.beta <- c(0,out$coefficients)
                #}
                y.rvar <- sum(out$residuals^2)/length(out$residuals) #ML?
                if(!lavaanOptions$meanstructure) {
                    x <- numeric(1L + length(x.beta) - 1L)
                    x[lavaanParTable$unco[lavaanParTable$op == "~~" &
                                          lavaanParTable$unco]] <- y.rvar
                    x[lavaanParTable$unco[lavaanParTable$op == "~" &
                                          lavaanParTable$unco]] <- x.beta[-1L]
                } else {
                    x <- numeric(1L + length(x.beta))
                    x[lavaanParTable$unco[lavaanParTable$op == "~~" &
                                          lavaanParTable$unco]] <- y.rvar
                    x[lavaanParTable$unco[lavaanParTable$op == "~" &
                                          lavaanParTable$unco]] <- x.beta[-1L]
                    x[lavaanParTable$unco[lavaanParTable$op == "~1" &
                                          lavaanParTable$unco]] <- x.beta[1L]
                }
                lavaanModel <- setModelParameters(lavaanModel, x = x)
                attr(x, "iterations") <- 1L; attr(x, "converged") <- TRUE
                attr(x, "control") <- control
                attr(x, "fx") <-
                computeObjective(lavaanModel, samplestats = lavaanSampleStats,
                                 estimator = lavaanOptions$estimator)
            } else if(checkLinearConstraints(lavaanModel) == TRUE) {

                require(quadprog)

                A.ceq <- A.cin <- matrix(0, lavaanModel@nx.free, 0)
                if(!is.null(body(lavaanModel@ceq.function)))
                    A.ceq <- t(lavJacobianC(func=lavaanModel@ceq.function, 
                                            x=rep(0,lavaanModel@nx.free)))
                if(!is.null(body(lavaanModel@cin.function)))
                    A.cin <- t(lavJacobianC(func=lavaanModel@cin.function, 
                                            x=rep(0,lavaanModel@nx.free)))
                A <- cbind(A.ceq, A.cin)
                con.jac <- t(A)

                # meanstructure?
                rvar.idx <- lavaanParTable$unco[lavaanParTable$op == "~~" &
                                                lavaanParTable$unco]
                if(lavaanOptions$meanstructure) {
                    # where is the intercept?
                    int.idx <- lavaanParTable$unco[lavaanParTable$op == "~1" &
                                                   lavaanParTable$unco]
                    # first intercept, then coefficients, remove resvar
                    A <- rbind(A[int.idx,,drop=FALSE], 
                               A[-c(int.idx,rvar.idx),,drop=FALSE])
                } else {
                    # add intercept, then coefficients, remove resvar
                    A <- rbind(rep(0,ncol(A)), A[-c(rvar.idx),,drop=FALSE])
                }

                ## forced zero intercept?
                #yvar <- vnames(lavaanParTable, "ov.y")
                #int.y.idx <- which(lavaanParTable$lhs == yvar &
                #                   lavaanParTable$op == "~1")
                #if(length(int.y.idx) > 0L &&
                #   lavaanParTable$ustart[int.y.idx] == 0) {
                #    lm.intercept <- FALSE
                #} else {
                #    lm.intercept <- TRUE
                #}
                X <- cbind(1,YX[,ov.x.idx])
                Y <- YX[,ov.y.idx]; X.Y <- crossprod(X, Y)
                #if(lm.intercept) {
                    X <- cbind(1,YX[,ov.x.idx])
                #} else {
                #    X <- YX[,ov.x.idx]
                #}
                X.X <- crossprod(X)
                out <- solve.QP(Dmat=X.X, dvec=X.Y, Amat=A, 
                                bvec=rep(0, NCOL(A)), 
                                meq=length(lavaanModel@x.ceq.idx))
                x.beta <- out$solution
                residuals <- Y - (X %*% x.beta)
                y.rvar <- sum(residuals^2)/length(residuals) #ML?
                if(!lavaanOptions$meanstructure) {
                    x <- numeric(1L + length(x.beta) - 1L)
                    x[lavaanParTable$free[lavaanParTable$op == "~~" &
                                          lavaanParTable$free]] <- y.rvar
                    x[lavaanParTable$free[lavaanParTable$op == "~" &
                                          lavaanParTable$free]] <- x.beta[-1L]
                } else {
                    x <- numeric(1L + length(x.beta))
                    x[lavaanParTable$free[lavaanParTable$op == "~~" &
                                          lavaanParTable$free]] <- y.rvar
                    x[lavaanParTable$free[lavaanParTable$op == "~" &
                                          lavaanParTable$free]] <- x.beta[-1L]
                    x[lavaanParTable$free[lavaanParTable$op == "~1" &
                                          lavaanParTable$free]] <- x.beta[1L]
                }
                lavaanModel <- setModelParameters(lavaanModel, x = x)
                attr(x, "iterations") <- 1L; attr(x, "converged") <- TRUE
                attr(x, "control") <- control
                attr(x, "fx") <-
                    computeObjective(lavaanModel, 
                                     samplestats = lavaanSampleStats,
                                     estimator = lavaanOptions$estimator)
                # for VCOV
                attr(con.jac, "inactive.idx") <- integer(0) # FIXME!!
                attr(con.jac, "cin.idx") <- seq_len(ncol(A.cin)) + ncol(A.ceq)
                attr(con.jac, "ceq.idx") <- seq_len(ncol(A.ceq))
                attr(x, "con.jac") <- con.jac
                con.lambda <- rep(1, nrow(con.jac)) # FIXME!
                attr(x, "con.lambda") <- con.lambda
            } else {
                # regular estimation after all
                x <- estimateModel(lavaanModel,
                               samplestats  = lavaanSampleStats,
                               X            = lavaanData@X,
                               options      = lavaanOptions,
                               cache        = lavaanCache,
                               control      = control)
                lavaanModel <- setModelParameters(lavaanModel, x = x)
            }
        } else {
            #cat("REGULAR\n")
            # regular estimation
            x <- estimateModel(lavaanModel,
                               samplestats  = lavaanSampleStats,
                               X            = lavaanData@X,
                               options      = lavaanOptions,
                               cache        = lavaanCache,
                               control      = control)
            lavaanModel <- setModelParameters(lavaanModel, x = x)
        }
        if(!is.null(attr(x, "con.jac"))) 
            lavaanModel@con.jac <- attr(x, "con.jac")
        if(!is.null(attr(x, "con.lambda")))
            lavaanModel@con.lambda <- attr(x, "con.lambda")
        # check if model has converged or not
        if(!attr(x, "converged") && lavaanOptions$warn) {
           warning("lavaan WARNING: model has NOT converged!")
        }
    } else {
        x <- numeric(0L)
        attr(x, "iterations") <- 0L; attr(x, "converged") <- FALSE
        attr(x, "control") <- control
        attr(x, "fx") <- 
            computeObjective(lavaanModel, samplestats = lavaanSampleStats, 
                             X=lavaanData@X,
                             estimator = lavaanOptions$estimator)
    }
    timing$Estimate <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 7. estimate vcov of free parameters (for standard errors)
    VCOV <- NULL
    if(lavaanOptions$se != "none" && lavaanModel@nx.free > 0L &&
       attr(x, "converged")) {
        if(verbose) cat("Computing VCOV for se =", lavaanOptions$se, "...")
        VCOV <- estimateVCOV(lavaanModel,
                             samplestats  = lavaanSampleStats,
                             options      = lavaanOptions,
                             data         = lavaanData,
                             partable     = lavaanParTable,
                             cache        = lavaanCache,
                             control      = control)
        if(verbose) cat(" done.\n")
    }
    timing$VCOV <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 8. compute test statistic (chi-square and friends)
    TEST <- NULL
    if(lavaanOptions$test != "none" && attr(x, "converged")) {
        if(verbose) cat("Computing TEST for test =", lavaanOptions$test, "...")
        TEST <- computeTestStatistic(lavaanModel,
                                     partable    = lavaanParTable,
                                     samplestats = lavaanSampleStats,
                                     options  = lavaanOptions,
                                     x        = x,
                                     VCOV     = VCOV,
                                     data     = lavaanData,
                                     cache    = lavaanCache,
                                     control  = control)
        if(verbose) cat(" done.\n")
    } else {
        TEST <- list(list(test="none", stat=NA, 
                     stat.group=rep(NA, lavaanData@ngroups), df=NA, 
                     refdistr="unknown", pvalue=NA))
    }
    timing$TEST <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 9. collect information about model fit (S4)
    lavaanFit <- Fit(partable = lavaanParTable, 
                     start    = lavaanStart, 
                     model    = lavaanModel,
                     x        = x, 
                     VCOV     = VCOV,
                     TEST     = TEST)
    timing$total <- (proc.time()[3] - start.time0)

    # 9b. post-fitting checks
    if(attr(x, "converged")) { # only if estimation was successful
        # 1. check for Heywood cases, negative (residual) variances, ...
        var.idx <- which(lavaanParTable$op == "~~" &
                         lavaanParTable$lhs == lavaanParTable$rhs)
        if(length(var.idx) > 0L && any(lavaanFit@est[var.idx] < 0.0))
            warning("lavaan WARNING: some estimated variances are negative")
        
        # 2. is cov.lv (PSI) positive definite?
        if(length(vnames(lavaanParTable, type="lv")) > 0L) {
            ETA <- computeETA(lavaanModel, samplestats=lavaanSampleStats)
            for(g in 1:lavaanData@ngroups) {
                txt.group <- ifelse(lavaanData@ngroups > 1L,
                                    paste("in group", g, ".", sep=""), ".")
                eigvals <- eigen(ETA[[g]], symmetric=TRUE, 
                                 only.values=TRUE)$values
                if(any(eigvals < 0))
                    warning("lavaan WARNING: covariance matrix of latent variables is not positive definite;", txt.group, " use inspect(fit,\"cov.lv\") to investigate.")
            }
        }
    }

    # 10. construct lavaan object
    lavaan <- new("lavaan",
                  call         = mc,                     # match.call
                  timing       = timing,                 # list
                  Options      = lavaanOptions,          # list
                  ParTable     = lavaanParTable,         # list
                  Data         = lavaanData,             # S4 class
                  SampleStats  = lavaanSampleStats,      # S4 class
                  Model        = lavaanModel,            # S4 class
                  Cache        = lavaanCache,            # list
                  Fit          = lavaanFit               # S4 class
                 )

    lavaan
}

# cfa + sem
cfa <- sem <- function(model = NULL, data = NULL,
    meanstructure = "default", fixed.x = "default",
    orthogonal = FALSE, std.lv = FALSE, std.ov = FALSE,
    missing = "default", ordered = NULL, 
    sample.cov = NULL, sample.cov.rescale = "default", sample.mean = NULL,
    sample.nobs = NULL, group = NULL, group.label = NULL,
    group.equal = "", group.partial = "", cluster = NULL, constraints = "",
    estimator = "default", likelihood = "default", 
    information = "default", se = "default", test = "default",
    bootstrap = 1000L, mimic = "default", representation = "default",
    do.fit = TRUE, control = list(), WLS.V = NULL, NACOV = NULL,
    start = "default", verbose = FALSE, warn = TRUE, debug = FALSE) {

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
    mc$auto.th         = TRUE
    mc$auto.delta      = TRUE
    mc[[1L]] <- as.name("lavaan")

    eval(mc, parent.frame())
}

# simple growth models
growth <- function(model = NULL, data = NULL,
    fixed.x = "default",
    orthogonal = FALSE, std.lv = FALSE, std.ov = FALSE,
    missing = "default", ordered = NULL, 
    sample.cov = NULL, sample.cov.rescale = "default", sample.mean = NULL,
    sample.nobs = NULL, group = NULL, group.label = NULL,
    group.equal = "", group.partial = "", cluster = NULL, constraints = "",
    estimator = "default", likelihood = "default", 
    information = "default", se = "default", test = "default",
    bootstrap = 1000L, mimic = "default", representation = "default",
    do.fit = TRUE, control = list(), WLS.V = NULL, NACOV = NULL,
    start = "default", verbose = FALSE, warn = TRUE, debug = FALSE) {

    mc <- match.call()

    mc$model.type      = "growth"
    mc$int.ov.free     = FALSE
    mc$int.lv.free     = TRUE
    mc$auto.fix.first  = !std.lv
    mc$auto.fix.single = TRUE
    mc$auto.var        = TRUE
    mc$auto.cov.lv.x   = TRUE
    mc$auto.cov.y      = TRUE
    mc$auto.th         = TRUE
    mc$auto.delta      = TRUE
    mc[[1L]] <- as.name("lavaan")

    eval(mc, parent.frame())
}
