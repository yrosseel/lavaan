# main user-visi ble cfa/sem/growth functions 
#
# initial version: YR 25/03/2009
# added lavoptions YR 02/08/2010
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
                   parameterization   = "default",

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
                   ridge              = 1e-5,

                   # multiple groups
                   group              = NULL,
                   group.label        = NULL,
                   group.equal        = '',
                   group.partial      = '',
                   group.w.free       = FALSE,

                   # clusters
                   cluster            = NULL,
             
                   # constraints
                   constraints        = '',

                   # estimation
                   estimator          = "default",
                   likelihood         = "default",
                   link               = "default",
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

                   # zero values
                   zero.add           = "default",
                   zero.keep.margins  = "default",

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
    # start timer
    start.time0 <- start.time <- proc.time()[3]; timing <- list()

    # 0a. store call
    mc  <- match.call()

    # 0b. get ov.names and ov.names.x (per group) -- needed for lavData()
    if(!is.null(slotParTable)) {
        FLAT <- slotParTable
    } else if(is.character(model)) {
        FLAT <- lavParseModelString(model)
    } else if(is.list(model)) {
        FLAT <- model
    }
    if(max(FLAT$group) < 2L) { # same model for all groups 
        ov.names   <- vnames(FLAT, type="ov")
        ov.names.y <- vnames(FLAT, type="ov.nox")
        ov.names.x <- vnames(FLAT, type="ov.x")
    } else { # different model per group
        ov.names <- lapply(1:max(FLAT$group),
                           function(x) vnames(FLAT, type="ov", group=x))
        ov.names.y <- lapply(1:max(FLAT$group),
                           function(x) vnames(FLAT, type="ov.nox", group=x))
        ov.names.x <- lapply(1:max(FLAT$group),
                           function(x) vnames(FLAT, type="ov.x", group=x))
    }

    # 0c categorical variables? -- needed for lavoptions
    if(any(FLAT$op == "|")) {
        categorical <- TRUE
    } else if(!is.null(data) && length(ordered) > 0L) {
        categorical <- TRUE
    } else if(is.data.frame(data) && 
              lav_dataframe_check_ordered(frame=data, ov.names=ov.names.y)) {
        categorical <- TRUE
    } else {
        categorical <- FALSE
    }

    # 1a. collect various options/flags and fill in `default' values
    #opt <- modifyList(formals(lavaan), as.list(mc)[-1])
    # force evaluation of `language` and/or `symbol` arguments
    #opt <- lapply(opt, function(x) if(typeof(x) %in% c("language", "symbol")) 
    #                                   eval(x, parent.frame()) else x)
    if(!is.null(slotOptions)) {
        lavoptions <- slotOptions
    } else {
        opt <- list(model = model, model.type = model.type,
            meanstructure = meanstructure, int.ov.free = int.ov.free,
            int.lv.free = int.lv.free, fixed.x = fixed.x, 
            orthogonal = orthogonal, std.lv = std.lv, 
            parameterization = parameterization,
            auto.fix.first = auto.fix.first, auto.fix.single = auto.fix.single,
            auto.var = auto.var, auto.cov.lv.x = auto.cov.lv.x, 
            auto.cov.y = auto.cov.y, auto.th = auto.th, 
            auto.delta = auto.delta, missing = missing, 
            group = group, categorical = categorical,
            group.equal = group.equal, group.partial = group.partial, 
            group.w.free = group.w.free,
            constraints = constraints,
            estimator = estimator, likelihood = likelihood, link = link,
            sample.cov.rescale = sample.cov.rescale,
            information = information, se = se, test = test, 
            bootstrap = bootstrap, mimic = mimic,
            zero.add = zero.add, zero.keep.margins = zero.keep.margins,
            representation = representation, do.fit = do.fit, verbose = verbose,
            warn = warn, debug = debug)
        lavoptions <- lav_options_set(opt)
    }
    timing$InitOptions <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # some additional checks for estimator="PML"
    if(lavoptions$estimator == "PML") {
        ovy <- unique( unlist(ov.names.y) )
        ovx <- unique( unlist(ov.names.x) )
        if(!is.null(slotData)) {
            ov.types <- slotData@ov$type[ slotData@ov$name %in% ovy ]
        } else {
            ov.types <- lav_dataframe_check_vartype(data, ov.names=ov.names.y)
        }
        # ordered argument?
        if(length(ordered) > 0L) {
            ord.idx <- which(ovy %in% ordered)
            ov.types[ord.idx] <- "ordered"
        }
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
        if(length(ovx) > 0L) {
            stop("lavaan ERROR: estimator=\"PML\" can not handle exogenous covariates (yet)")
        }

        # 3. warn that this is still experimental
        message("Please note: the PML estimator is still under development.\n",
                "Future releases will improve speed, and allow for mixed ord/cont variables.\n",
                "Research on how to compute a proper goodness-of-fit test is ongoing.\n")
    }

    # 1b. check data/sample.cov and get the number of groups
    if(!is.null(slotData)) {
        lavdata <- slotData
    } else {
        if(categorical) {
            ov.names <- ov.names.y
        } 
        lavdata <- lavData(data        = data,
                           group       = group,
                           group.label = group.label,
                           ov.names    = ov.names,
                           ordered     = ordered,
                           ov.names.x  = ov.names.x,
                           std.ov      = std.ov,
                           missing     = lavoptions$missing,
                           sample.cov  = sample.cov,
                           sample.mean = sample.mean,
                           sample.nobs = sample.nobs,
                           warn        = lavoptions$warn)

        # what have we learned from the data?
        if("ordered" %in% lavdata@ov$type) {
            if(lavoptions$estimator == "ML")
                stop("lavaan ERROR: estimator ML for ordered data is not supported yet. Use WLSMV instead.")
            # Mplus style
            lavoptions$meanstructure <- TRUE
        }
    }
    if(lavdata@data.type == "none") {
        do.fit <- FALSE; start <- "simple"
        lavoptions$se <- "none"; lavoptions$test <- "none"
    } else if(lavdata@data.type == "moment") {
        # catch here some options that will not work with moments
        if(lavoptions$se == "bootstrap") {
            stop("lavaan ERROR: bootstrapping requires full data")
        }
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
    if(debug) print(str(lavdata))


    # 2a. construct paramter table: description of the user-specified model
    if(!is.null(slotParTable)) {
        lavpartable <- slotParTable
    } else if(is.character(model)) {

        # check FLAT before we proceed
        if(debug) print(as.data.frame(FLAT))

        # check 1: catch ~~ of fixed.x covariates if categorical
        if(lavoptions$categorical) {
            tmp <- vnames(FLAT, type="ov.x", ov.x.fatal=TRUE)
        }
        # check 2: catch ~1 of ordered variables (unless they are fixed)
        #if(lavoptions$categorical) {
        #    int.idx <- which(FLAT$op == "~1" & FLAT$fixed == "")
        #    if(length(int.idx) > 0L) {
        #        INT <- FLAT$lhs[int.idx]
        #        ORD <- lavdata@ov$name[ lavdata@ov$type == "ordered" ]
        #        if(any(INT %in% ORD))        
        #            stop("lavaan ERROR: model syntax contains free intercepts for ordinal dependent variable(s): [", paste(INT, collapse=" "), 
        #                 "];\n  Please remove them and try again.")
        #    }
        #}

        lavpartable <- 
            lavaanify(model            = FLAT,
                      meanstructure    = lavoptions$meanstructure, 
                      int.ov.free      = lavoptions$int.ov.free,
                      int.lv.free      = lavoptions$int.lv.free,
                      orthogonal       = lavoptions$orthogonal, 
                      fixed.x          = lavoptions$fixed.x,
                      std.lv           = lavoptions$std.lv,
                      parameterization = lavoptions$parameterization,
                      constraints      = constraints,

                      auto.fix.first   = lavoptions$auto.fix.first,
                      auto.fix.single  = lavoptions$auto.fix.single,
                      auto.var         = lavoptions$auto.var,
                      auto.cov.lv.x    = lavoptions$auto.cov.lv.x,
                      auto.cov.y       = lavoptions$auto.cov.y,
                      auto.th          = lavoptions$auto.th,
                      auto.delta       = lavoptions$auto.delta,

                      varTable         = lavdata@ov,
                      ngroups          = lavdata@ngroups,
                      group.equal      = lavoptions$group.equal, 
                      group.partial    = lavoptions$group.partial,
                      group.w.free     = lavoptions$group.w.free,
                      debug            = lavoptions$debug,
                      warn             = lavoptions$warn,

                      as.data.frame.   = FALSE)

    } else if(is.list(model)) {
        # two possibilities: either model is already lavaanified
        # or it is something else...
        if(!is.null(model$lhs) && !is.null(model$op)  &&
           !is.null(model$rhs) && !is.null(model$free)) {
            lavpartable <- as.list(model)
        } else if(is.character(model[[1]])) {
            stop("lavaan ERROR: model is a list, but not a parameterTable?")
        }
    } else {
        cat("model type: ", class(model), "\n")
        stop("lavaan ERROR: model is not of type character or list")
    }
    if(debug) print(as.data.frame(lavpartable))

    # 2b. change meanstructure flag?
    if(any(lavpartable$op == "~1")) lavoptions$meanstructure <- TRUE

    # 2c. get partable attributes
    lavpta <- lav_partable_attributes(lavpartable)
    timing$ParTable <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 3. get sample statistics
    if(!is.null(slotSampleStats)) {
        lavsamplestats <- slotSampleStats
    } else if(lavdata@data.type == "full") {
        lavsamplestats <- lav_samplestats_from_data(
                       Data          = lavdata,
                       missing       = lavoptions$missing,
                       rescale       = (lavoptions$estimator == "ML" &&
                                        lavoptions$likelihood == "normal"),
                       estimator     = lavoptions$estimator,
                       mimic         = lavoptions$mimic,
                       meanstructure = lavoptions$meanstructure,
                       group.w.free  = lavoptions$group.w.free,
                       missing.h1    = (lavoptions$missing != "listwise"),
                       WLS.V             = WLS.V,
                       NACOV             = NACOV,
                       ridge             = ridge,
                       zero.add          = lavoptions$zero.add,
                       zero.keep.margins = lavoptions$zero.keep.margins,
                       debug             = lavoptions$debug,
                       verbose           = lavoptions$verbose)
                                                 
    } else if(lavdata@data.type == "moment") {
        lavsamplestats <- lav_samplestats_from_moments(
                           sample.cov    = sample.cov,
                           sample.mean   = sample.mean,
                           sample.nobs   = sample.nobs,
                           ov.names      = lavpta$vnames$ov,
                           estimator     = lavoptions$estimator,
                           mimic         = lavoptions$mimic,
                           meanstructure = lavoptions$meanstructure,
                           group.w.free  = lavoptions$group.w.free,
                           WLS.V         = WLS.V,
                           NACOV         = NACOV,
                           ridge         = ridge,
                           rescale       = lavoptions$sample.cov.rescale)
    } else {
        # no data
        th.idx <- vector("list", length=lavdata@ngroups)
        for(g in 1:lavdata@ngroups) {
            th.idx[[g]] <- lav_partable_ov_idx(lavpartable, type="th")
        }
        lavsamplestats <- new("lavSampleStats", ngroups=lavdata@ngroups,
                                 nobs=as.list(rep(0L, lavdata@ngroups)),
                                 cov.x=vector("list",length=lavdata@ngroups),
                                 th.idx=th.idx,
                                 missing.flag=FALSE)
    }
    timing$Sample <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
    if(debug) print(str(lavsamplestats))

    # 4. compute some reasonable starting values 
    if(!is.null(slotModel)) {
        lavmodel <- slotModel
        lavaanStart <- lav_model_get_parameters(lavmodel, type="user")
        lavpartable$start <- lavaanStart
        timing$Start <- (proc.time()[3] - start.time)
        start.time <- proc.time()[3]
        timing$Model <- (proc.time()[3] - start.time)
        start.time <- proc.time()[3]
    } else {
        lavaanStart <- 
            lav_start(start.method = start,
                           partable     = lavpartable, 
                           samplestats  = lavsamplestats,
                           model.type   = lavoptions$model.type,
                           mimic        = lavoptions$mimic,
                           debug        = lavoptions$debug)
        timing$Start <- (proc.time()[3] - start.time)
        start.time <- proc.time()[3]

        lavpartable$start <- lavaanStart
        #print(as.data.frame(lavpartable))

        # 5. construct internal model (S4) representation
        lavmodel <- 
            lav_model(partable         = lavpartable,
                      representation   = lavoptions$representation,
                      th.idx           = lavsamplestats@th.idx,
                      parameterization = lavoptions$parameterization,
                      debug            = lavoptions$debug)
        timing$Model <- (proc.time()[3] - start.time)
        start.time <- proc.time()[3]
  
        # if no data, call lav_model_set_parameters once (for categorical case)
        if(lavdata@data.type == "none" && lavmodel@categorical) {
            lavmodel <- lav_model_set_parameters(lavmodel,
                                              x=lav_model_get_parameters(lavmodel),                                              estimator=lavoptions$estimator)
        }
    }

    # check for categorical
    #if(lavmodel@categorical && lavoptions$se == "bootstrap") {
    #    stop("lavaan ERROR: bootstrap not supported (yet) for categorical data")
    #}

    # prepare cache -- stuff needed for estimation, but also post-estimation
    lavcache <- vector("list", length=lavdata@ngroups)
    if(lavoptions$estimator == "PML") {
        TH <- computeTH(lavmodel)
        for(g in 1:lavdata@ngroups) {
            nvar <- ncol(lavdata@X[[g]])
            nobs <- nrow(lavdata@X[[g]])
            th.idx <- lavmodel@th.idx[[g]]
            # pairwise tables, as a long vector
            PW <- pairwiseTables(data=lavdata@X[[g]], no.x=nvar)$pairTables
            bifreq <- as.numeric(unlist(PW))
            ### FIXME!!! Check for zero cells!!
            #zero.idx <- which(bifreq == 0)
            #bifreq[zero.idx] <- 0.001 ####????
            LONG <- LongVecInd(no.x               = nvar,
                               all.thres          = TH[[g]],
                               index.var.of.thres = th.idx)
            lavcache[[g]] <- list(bifreq=bifreq, nobs=nobs, LONG=LONG)
        }
    }
    # copy response patterns to cache -- FIXME!! (data not included 
    # in Model only functions)
    if(lavdata@data.type == "full" && !is.null(lavdata@Rp[[1L]])) {
        for(g in 1:lavdata@ngroups) {
            lavcache[[g]]$pat <- lavdata@Rp[[g]]$pat
        }
    }

    # If estimator = MML, store Gauss-Hermite nodes/weights
    if(lavoptions$estimator == "MML") {
        if(!is.null(control$nGH)) {
            nGH <- control$nGH
        } else {
            nGH <- 21L
        }
        for(g in 1:lavdata@ngroups) {
            # count only the ones with non-normal indicators
            #nfac <- lavpta$nfac.nonnormal[[g]]
            nfac <- lavpta$nfac[[g]]
            lavcache[[g]]$GH <- 
                lav_gauss_hermite_xw_dnorm(n=nGH, revert=FALSE, ndim = nfac)
        }
    }

    # 6. estimate free parameters
    x <- NULL
    if(do.fit && lavoptions$estimator != "none" && 
       lavmodel@nx.free > 0L) {
        # catch simple linear regression models
        if(lavdata@data.type == "full" &&
           length(unique(lavpartable$lhs[lavpartable$op == "~"])) == 1L && 
           length(vnames(lavpartable,   "lv")) == 0L &&
           #FALSE && # to debug
           sum(nchar(FLAT$fixed)) == 0 && # no fixed values in parTable
                                          # this includes intercepts!!
           !categorical &&
           !lavmodel@eq.constraints &&
           length(lavdata@X) > 0L &&
           lavdata@ngroups == 1L &&
           lavoptions$fixed &&
           lavoptions$missing == "listwise") {
            # simple univariate regression
            ov.y.idx <- match(vnames(lavpartable, "ov.y"), 
                              lavdata@ov.names[[1L]])
            ov.x.idx <- match(vnames(lavpartable, "ov.x"), 
                              lavdata@ov.names[[1L]])
            YX <- lavdata@X[[1L]]
            #print(head(YX))
            # constraints?
            if(sum(length(lavmodel@x.ceq.idx) + 
                   length(lavmodel@x.cin.idx)) == 0L) {

                ## forced zero intercept?
                #yvar <- vnames(lavpartable, "ov.y")
                #int.y.idx <- which(lavpartable$lhs == yvar &
                #                   lavpartable$op == "~1")
                #if(length(int.y.idx) > 0L && 
                #   !is.na(lavpartable$ustart[int.y.idx]) &&
                #   lavpartable$ustart[int.y.idx] == 0) {
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
                if(!lavoptions$meanstructure) {
                    x <- numeric(1L + length(x.beta) - 1L)
                    x[lavpartable$unco[lavpartable$op == "~~" &
                                          lavpartable$unco]] <- y.rvar
                    x[lavpartable$unco[lavpartable$op == "~" &
                                          lavpartable$unco]] <- x.beta[-1L]
                } else {
                    x <- numeric(1L + length(x.beta))
                    x[lavpartable$unco[lavpartable$op == "~~" &
                                          lavpartable$unco]] <- y.rvar
                    x[lavpartable$unco[lavpartable$op == "~" &
                                          lavpartable$unco]] <- x.beta[-1L]
                    x[lavpartable$unco[lavpartable$op == "~1" &
                                          lavpartable$unco]] <- x.beta[1L]
                }
                lavmodel <- lav_model_set_parameters(lavmodel, x = x,
                                           estimator=lavoptions$estimator)
                attr(x, "iterations") <- 1L; attr(x, "converged") <- TRUE
                attr(x, "control") <- control
                FX <- try(lav_model_objective(lavmodel, 
                                           samplestats = lavsamplestats,
                                           estimator = lavoptions$estimator,
                                           link = lavoptions$link),
                          silent=TRUE)
                if(inherits(FX, "try-error")) {
                    # eg non-full rank design matrix
                    FX <- as.numeric(NA)
                    attr(FX, "fx.group") <- rep(as.numeric(NA), 
                                                lavdata@ngroups)
                    attr(x, "fx") <- as.numeric(NA)
                } 
                attr(x, "fx") <- FX
            } else if(lav_constraints_check_linear(lavmodel) == TRUE) {

                A.ceq <- A.cin <- matrix(0, lavmodel@nx.free, 0)
                if(!is.null(body(lavmodel@ceq.function)))
                    A.ceq <- t(lavJacobianC(func=lavmodel@ceq.function, 
                                            x=rep(0,lavmodel@nx.free)))
                if(!is.null(body(lavmodel@cin.function)))
                    A.cin <- t(lavJacobianC(func=lavmodel@cin.function, 
                                            x=rep(0,lavmodel@nx.free)))
                A <- cbind(A.ceq, A.cin)
                con.jac <- t(A)

                # meanstructure?
                rvar.idx <- lavpartable$unco[lavpartable$op == "~~" &
                                                lavpartable$unco]
                if(lavoptions$meanstructure) {
                    # where is the intercept?
                    int.idx <- lavpartable$unco[lavpartable$op == "~1" &
                                                   lavpartable$unco]
                    # first intercept, then coefficients, remove resvar
                    A <- rbind(A[int.idx,,drop=FALSE], 
                               A[-c(int.idx,rvar.idx),,drop=FALSE])
                } else {
                    # add intercept, then coefficients, remove resvar
                    A <- rbind(rep(0,ncol(A)), A[-c(rvar.idx),,drop=FALSE])
                }

                ## forced zero intercept?
                #yvar <- vnames(lavpartable, "ov.y")
                #int.y.idx <- which(lavpartable$lhs == yvar &
                #                   lavpartable$op == "~1")
                #if(length(int.y.idx) > 0L &&
                #   lavpartable$ustart[int.y.idx] == 0) {
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
                                bvec=rep(0, NCOL(A)), ### FIXME!!! always zero
                                meq=length(lavmodel@x.ceq.idx))
                x.beta <- out$solution
                residuals <- Y - (X %*% x.beta)
                y.rvar <- sum(residuals^2)/length(residuals) #ML?
                if(!lavoptions$meanstructure) {
                    x <- numeric(1L + length(x.beta) - 1L)
                    x[lavpartable$free[lavpartable$op == "~~" &
                                          lavpartable$free]] <- y.rvar
                    x[lavpartable$free[lavpartable$op == "~" &
                                          lavpartable$free]] <- x.beta[-1L]
                } else {
                    x <- numeric(1L + length(x.beta))
                    x[lavpartable$free[lavpartable$op == "~~" &
                                          lavpartable$free]] <- y.rvar
                    x[lavpartable$free[lavpartable$op == "~" &
                                          lavpartable$free]] <- x.beta[-1L]
                    x[lavpartable$free[lavpartable$op == "~1" &
                                          lavpartable$free]] <- x.beta[1L]
                }
                lavmodel <- lav_model_set_parameters(lavmodel, x = x,
                                       estimator=lavoptions$estimator)
                attr(x, "iterations") <- 1L; attr(x, "converged") <- TRUE
                attr(x, "control") <- control
                attr(x, "fx") <-
                    lav_model_objective(lavmodel, 
                                     samplestats = lavsamplestats,
                                     estimator = lavoptions$estimator,
                                     link = lavoptions$link)
                # for VCOV
                attr(con.jac, "inactive.idx") <- integer(0) # FIXME!!
                attr(con.jac, "cin.idx") <- seq_len(ncol(A.cin)) + ncol(A.ceq)
                attr(con.jac, "ceq.idx") <- seq_len(ncol(A.ceq))
                attr(x, "con.jac") <- con.jac
                con.lambda <- rep(1, nrow(con.jac)) # FIXME!
                attr(x, "con.lambda") <- con.lambda
            } else {
                # regular estimation after all
                x <- lav_model_estimate(lavmodel,
                               samplestats  = lavsamplestats,
                               X            = lavdata@X,
                               options      = lavoptions,
                               cache        = lavcache,
                               control      = control)
                lavmodel <- lav_model_set_parameters(lavmodel, x = x,
                     estimator=lavoptions$estimator)
            }
        } else {
            #cat("REGULAR\n")
            # regular estimation
            x <- lav_model_estimate(lavmodel,
                               samplestats  = lavsamplestats,
                               X            = lavdata@X,
                               options      = lavoptions,
                               cache        = lavcache,
                               control      = control)
            lavmodel <- lav_model_set_parameters(lavmodel, x = x,
                             estimator=lavoptions$estimator)
        }
        if(!is.null(attr(x, "con.jac"))) 
            lavmodel@con.jac <- attr(x, "con.jac")
        if(!is.null(attr(x, "con.lambda")))
            lavmodel@con.lambda <- attr(x, "con.lambda")
        # check if model has converged or not
        if(!attr(x, "converged") && lavoptions$warn) {
           warning("lavaan WARNING: model has NOT converged!")
        }
    } else {
        x <- numeric(0L)
        attr(x, "iterations") <- 0L; attr(x, "converged") <- FALSE
        attr(x, "control") <- control
        attr(x, "fx") <- 
            lav_model_objective(lavmodel, samplestats = lavsamplestats, 
                             X=lavdata@X, cache = lavcache,
                             estimator = lavoptions$estimator,
                             link = lavoptions$link)
    }
    timing$Estimate <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 7. estimate vcov of free parameters (for standard errors)
    VCOV <- NULL
    if(lavoptions$se != "none" && lavmodel@nx.free > 0L &&
       attr(x, "converged")) {
        if(verbose) cat("Computing VCOV for se =", lavoptions$se, "...")
        VCOV <- estimateVCOV(lavmodel,
                             samplestats  = lavsamplestats,
                             options      = lavoptions,
                             data         = lavdata,
                             partable     = lavpartable,
                             cache        = lavcache,
                             control      = control)
        if(verbose) cat(" done.\n")
    }
    timing$VCOV <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 8. compute test statistic (chi-square and friends)
    TEST <- NULL
    if(lavoptions$test != "none" && attr(x, "converged")) {
        if(verbose) cat("Computing TEST for test =", lavoptions$test, "...")
        TEST <- computeTestStatistic(lavmodel,
                                     partable    = lavpartable,
                                     samplestats = lavsamplestats,
                                     options  = lavoptions,
                                     x        = x,
                                     VCOV     = VCOV,
                                     data     = lavdata,
                                     cache    = lavcache,
                                     control  = control)
        if(verbose) cat(" done.\n")
    } else {
        TEST <- list(list(test="none", stat=NA, 
                     stat.group=rep(NA, lavdata@ngroups), df=NA, 
                     refdistr="unknown", pvalue=NA))
    }
    timing$TEST <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]

    # 9. collect information about model fit (S4)
    lavfit <- Fit(partable = lavpartable, 
                     model    = lavmodel,
                     x        = x, 
                     VCOV     = VCOV,
                     TEST     = TEST)
    timing$total <- (proc.time()[3] - start.time0)

    # 9b. post-fitting checks
    if(attr(x, "converged")) { # only if estimation was successful
        # 1. check for Heywood cases, negative (residual) variances, ...
        var.idx <- which(lavpartable$op == "~~" &
                         !lavpartable$lhs %in% unlist(lavpta$vnames$ov.ord) &
                         lavpartable$lhs == lavpartable$rhs)
        if(length(var.idx) > 0L && any(lavfit@est[var.idx] < 0.0))
            warning("lavaan WARNING: some estimated variances are negative")
        
        # 2. is cov.lv (PSI) positive definite?
        if(length(vnames(lavpartable, type="lv.regular")) > 0L) {
            ETA <- computeVETA(lavmodel, samplestats=lavsamplestats)
            for(g in 1:lavdata@ngroups) {
                txt.group <- ifelse(lavdata@ngroups > 1L,
                                    paste("in group", g, ".", sep=""), "")
                eigvals <- eigen(ETA[[g]], symmetric=TRUE, 
                                 only.values=TRUE)$values
                if(any(eigvals < -1 * .Machine$double.eps^(3/4)))
                    warning("lavaan WARNING: covariance matrix of latent variables is not positive definite;", txt.group, " use inspect(fit,\"cov.lv\") to investigate.")
            }
        }

        # 3. is THETA positive definite (but only for numeric variables)
        THETA <- computeTHETA(lavmodel)
        for(g in 1:lavdata@ngroups) { 
                num.idx <- lavmodel@num.idx[[g]]
                if(length(num.idx) > 0L) {
                    txt.group <- ifelse(lavdata@ngroups > 1L,
                                        paste("in group", g, ".", sep=""), "")
                    eigvals <- eigen(THETA[[g]][num.idx,num.idx,drop=FALSE], 
                                     symmetric=TRUE,
                                     only.values=TRUE)$values
                    if(any(eigvals < -1 * .Machine$double.eps^(3/4))) {
                        warning("lavaan WARNING: residual covariance matrix is not positive definite;", txt.group, " use inspect(fit,\"cov.ov\") to investigate.")
                    }
                }
        }
    }

    # 10. construct lavaan object
    lavaan <- new("lavaan",
                  call         = mc,                  # match.call
                  timing       = timing,              # list
                  Options      = lavoptions,          # list
                  ParTable     = lavpartable,         # list
                  pta          = lavpta,              # list
                  Data         = lavdata,             # S4 class
                  SampleStats  = lavsamplestats,      # S4 class
                  Model        = lavmodel,            # S4 class
                  Cache        = lavcache,            # list
                  Fit          = lavfit               # S4 class
                 )

    lavaan
}

# cfa + sem
cfa <- sem <- function(model = NULL, data = NULL,
    meanstructure = "default", fixed.x = "default",
    orthogonal = FALSE, std.lv = FALSE, 
    parameterization = "default", std.ov = FALSE,
    missing = "default", ordered = NULL, 
    sample.cov = NULL, sample.cov.rescale = "default", sample.mean = NULL,
    sample.nobs = NULL, ridge = 1e-5,
    group = NULL, group.label = NULL,
    group.equal = "", group.partial = "", group.w.free = FALSE,
    cluster = NULL, constraints = "",
    estimator = "default", likelihood = "default", link = "default",
    information = "default", se = "default", test = "default",
    bootstrap = 1000L, mimic = "default", representation = "default",
    do.fit = TRUE, control = list(), WLS.V = NULL, NACOV = NULL,
    zero.add = "default", zero.keep.margins = "default", start = "default",
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
    mc$auto.th         = TRUE
    mc$auto.delta      = TRUE
    mc[[1L]] <- quote(lavaan::lavaan)

    eval(mc, parent.frame())
}

# simple growth models
growth <- function(model = NULL, data = NULL,
    fixed.x = "default",
    orthogonal = FALSE, std.lv = FALSE, 
    parameterization = "default", std.ov = FALSE,
    missing = "default", ordered = NULL, 
    sample.cov = NULL, sample.cov.rescale = "default", sample.mean = NULL,
    sample.nobs = NULL, ridge = 1e-5,
    group = NULL, group.label = NULL,
    group.equal = "", group.partial = "", group.w.free = FALSE,
    cluster = NULL, constraints = "",
    estimator = "default", likelihood = "default", link = "default",
    information = "default", se = "default", test = "default",
    bootstrap = 1000L, mimic = "default", representation = "default",
    do.fit = TRUE, control = list(), WLS.V = NULL, NACOV = NULL,
    zero.add = "default", zero.keep.margins = "default", start = "default",
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
    mc$auto.th         = TRUE
    mc$auto.delta      = TRUE
    mc[[1L]] <- quote(lavaan::lavaan)

    eval(mc, parent.frame())
}
