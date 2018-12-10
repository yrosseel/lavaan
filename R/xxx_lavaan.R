# main user-visible cfa/sem/growth functions
#
# initial version: YR 25/03/2009
# added lavoptions YR 02/08/2010
# major revision: YR 9/12/2010: - new workflow (since 0.4-5)
#                               - merge cfa/sem/growth functions
# YR 25/02/2012: changed data slot (from list() to S4); data@X contains data

# YR 26 Jan 2017: use '...' to capture the never-ending list of options

lavaan <- function(# user-specified model: can be syntax, parameter Table, ...
                   model              = NULL,
                   # data (second argument, most used)
                   data               = NULL,

                   # variable information
                   ordered            = NULL,

                   # sampling weights
                   sampling.weights   = NULL,

                   # summary data
                   sample.cov         = NULL,
                   sample.mean        = NULL,
                   sample.th          = NULL,
                   sample.nobs        = NULL,

                   # multiple groups?
                   group              = NULL,

                   # multiple levels?
                   cluster            = NULL,

                   # constraints
                   constraints        = '',

                   # user-specified variance matrices
                   WLS.V              = NULL,
                   NACOV              = NULL,

                   # full slots from previous fits
                   slotOptions        = NULL,
                   slotParTable       = NULL,
                   slotSampleStats    = NULL,
                   slotData           = NULL,
                   slotModel          = NULL,
                   slotCache          = NULL,
                   sloth1             = NULL,

                   # options (dotdotdot)
                   ...
                  ) {
    # start timer
    start.time0 <- start.time <- proc.time()[3]; timing <- list()

    # 0a. store call
    mc  <- match.call(expand.dots = TRUE)

    # handle dotdotdot
    dotdotdot <- list(...)

    # backwards compatibility, control= argument (<0.5-23)
    if(!is.null(dotdotdot$control)) {
        # optim.method
        if(!is.null(dotdotdot$control$optim.method)) {
           dotdotdot$optim.method <- dotdotdot$control$optim.method
        }
        # cor.optim.method
        if(!is.null(dotdotdot$control$cor.optim.method)) {
            dotdotdot$optim.method.cor <- dotdotdot$control$cor.optim.method
        }
        # control$optim.force.converged
        if(!is.null(dotdotdot$control$optim.force.converged)) {
            dotdotdot$optim.force.converged <-
                dotdotdot$control$optim.force.converged
        }
        # gradient
        if(!is.null(dotdotdot$control$gradient)) {
            dotdotdot$optim.gradient <- dotdotdot$control$gradient
        }
        if(!is.null(dotdotdot$gradient)) {
            dotdotdot$optim.gradient <- dotdotdot$gradient
        }
        # init_nelder_mead
        if(!is.null(dotdotdot$control$init_nelder_mead)) {
            dotdotdot$optim.init_nelder_mead <-
                dotdotdot$control$init_nelder_mead
        }
    }


    ######################
    #### 1. ov.names  ####
    ######################
    # 1a. get ov.names and ov.names.x (per group) -- needed for lavData()
    if(!is.null(slotParTable)) {
        FLAT <- slotParTable
    } else if(is.character(model)) {
        FLAT <- lavParseModelString(model)
    } else if( inherits(model, "formula") ) {
        # two typical cases:
        # 1. regression type formula
        # 2. no quotes, eg f =~ x1 + x2 + x3
        tmp <- as.character(model)
        if(tmp[1] == "~" && length(tmp) == 2L) {
            # looks like an unquoted single factor model f =~ something
            warning("lavaan WARNING: model seems to be a formula; please enclose the model syntax between quotes")
            # create model and hope for the best
            model.bis <- paste("f =", paste(tmp, collapse= " "), sep = "")
            FLAT <- lavParseModelString(model.bis)
        } else if(tmp[1] == "~" && length(tmp) == 3L) {
            # looks like a (unquoted) regression formula
            warning("lavaan WARNING: model seems to be a formula; please enclose the model syntax between quotes")
            # create model and hope for the best
            model.bis <- paste(tmp[2], tmp[1], tmp[3])
            FLAT <- lavParseModelString(model.bis)
        } else {
            stop("lavaan ERROR: model seems to be a formula; please enclose the model syntax between quotes")
        }
    } else if(inherits(model, "lavaan")) {
        # hm, a lavaan model; let's try to extract the parameter table
        # and see what happens
        FLAT <- parTable(model)
    } else if(is.list(model)) {
        # two possibilities: either model is already lavaanified
        # or it is something else...

        # look for the bare minimum columns: lhs - op - rhs
        if(!is.null(model$lhs) && !is.null(model$op)  &&
           !is.null(model$rhs) && !is.null(model$free)) {

            # ok, we have something that looks like a parameter table
            # FIXME: we need to check for redundant arguments
            # (but if cfa/sem was used, we can not trust the call)
            # redundant <- c("meanstructure", "int.ov.free", "int.lv.free",
            #        "fixed.x", "orthogonal", "std.lv", "parameterization",
            #        "auto.fix.first", "auto.fix.single", "auto.var",
            #        "auto.cov.lv.x", "auto.cov.y", "auto.th", "auto.delta")
            FLAT <- model

            # fix semTools issue here? for auxiliary() which does not use
            # block column yet
            if(!is.null(FLAT$block)) {
                N <- length(FLAT$lhs)
                if(length(FLAT$block) != N) {
                    FLAT$block <- FLAT$group
                }
                if(any(is.na(FLAT$block))) {
                    FLAT$block <- FLAT$group
                }
            } else if(!is.null(FLAT$group)) {
                FLAT$block <- FLAT$group
            }

        } else {
            bare.minimum <- c("lhs", "op", "rhs", "free")
            missing.idx <- is.na(match(bare.minimum, names(model)))
            missing.txt <- paste(bare.minimum[missing.idx], collapse = ", ")
            stop("lavaan ERROR: model is a list, but not a parameterTable?",
                 "\n  lavaan  NOTE: ",
                 "missing column(s) in parameter table: [", missing.txt, "]")
        }
    } else if(is.null(model)) {
        stop("lavaan ERROR: model is NULL!")
    }

    # group blocks?
    if(any(FLAT$op == ":" & tolower(FLAT$lhs) == "group")) {
        # here, we only need to figure out:
        # - ngroups
        # - ov's per group
        # - FIXME: we need a more efficient way, avoiding lavaanify/vnames
        group.idx <- which(FLAT$op == ":" & tolower(FLAT$lhs) == "group")
        # replace by 'group' (in case we got 'Group'):
        FLAT$lhs[group.idx] <- "group"
        tmp.group.values <- unique(FLAT$rhs[group.idx])
        tmp.ngroups <- length(tmp.group.values)
        tmp.lav <- lavaanify(FLAT, ngroups = tmp.ngroups, warn = FALSE)
        ov.names <- ov.names.y <- ov.names.x <- vector("list",
                                                       length = tmp.ngroups)
        for(g in seq_len(tmp.ngroups)) {
            ov.names[[g]]   <- unique(unlist(lav_partable_vnames(tmp.lav,
                                type = "ov", group = tmp.group.values[g])))
            ov.names.y[[g]] <- unique(unlist(lav_partable_vnames(tmp.lav,
                                type = "ov.nox", group = tmp.group.values[g])))
            ov.names.x[[g]] <- unique(unlist(lav_partable_vnames(tmp.lav,
                                type = "ov.x", group = tmp.group.values[g])))
        }
    } else if(!is.null(FLAT$group)) {
        # user-provided full partable with group column!
        ngroups <- lav_partable_ngroups(FLAT)
        if(ngroups > 1L) {
            group.values <- lav_partable_group_values(FLAT)
            ov.names <- ov.names.y <- ov.names.x <- vector("list",
                                                       length = ngroups)
            for(g in seq_len(ngroups)) {
                ov.names[[g]]   <- lav_partable_vnames(FLAT, type = "ov",
                                                       group = group.values[g])
                ov.names.y[[g]] <- lav_partable_vnames(FLAT, type = "ov.nox",
                                                       group = group.values[g])
                ov.names.x[[g]] <- lav_partable_vnames(FLAT, type = "ov.x",
                                                       group = group.values[g])
            }
        } else {
            ov.names   <- lav_partable_vnames(FLAT, type = "ov")
            ov.names.y <- lav_partable_vnames(FLAT, type = "ov.nox")
            ov.names.x <- lav_partable_vnames(FLAT, type = "ov.x")
        }
    } else {
        ov.names   <- lav_partable_vnames(FLAT, type = "ov")
        ov.names.y <- lav_partable_vnames(FLAT, type = "ov.nox")
        ov.names.x <- lav_partable_vnames(FLAT, type = "ov.x")
    }

    # handle ov.names.l
    if(any(FLAT$op == ":" & tolower(FLAT$lhs) == "level")) {

        # check for cluster argument
        if(!is.null(data) && is.null(cluster)) {
            stop("lavaan ERROR: cluster argument is missing.")
        }

        # here, we only need to figure out:
        # - nlevels
        # - ov's per level
        # - FIXME: we need a more efficient way, avoiding lavaanify/vnames

        group.idx <- which(FLAT$op == ":" & FLAT$lhs == "group")
        tmp.group.values <- unique(FLAT$rhs[group.idx])
        tmp.ngroups <- max(c(length(tmp.group.values), 1))

        level.idx <- which(FLAT$op == ":" & tolower(FLAT$lhs) == "level")
        # replace by "level" (in case we got 'Level')
        FLAT$lhs[level.idx] <- "level"
        tmp.level.values <- unique(FLAT$rhs[level.idx])
        tmp.nlevels <- length(tmp.level.values)

        # we need at least 2 levels (for now)
        if(tmp.nlevels < 2L) {
            stop("lavaan ERROR: when data is clustered, you must specify a model\n", "  for each level in the model syntax (for now); see example(Demo.twolevel)")
        }

        tmp.lav <- lavaanify(FLAT, ngroups = tmp.ngroups, warn = FALSE)
        # check for empty levels
        if(max(tmp.lav$level) < 2L) {
            stop("lavaan ERROR: at least one level has no model syntax; you must specify a model for each level in the model syntax (for now); see example(Demo.twolevel)")
        }
        ov.names.l <- vector("list", length = tmp.ngroups) # per group

        for(g in seq_len(tmp.ngroups)) {
                ov.names.l[[g]] <- vector("list", length = tmp.nlevels)
            for(l in seq_len(tmp.nlevels)) {
                if(tmp.ngroups > 1L) {
                    ov.names.l[[g]][[l]] <-
                    unique(unlist(lav_partable_vnames(tmp.lav,
                                                  type = "ov",
                                                  group = tmp.group.values[g],
                                                  level = tmp.level.values[l])))
                } else {
                    ov.names.l[[g]][[l]] <-
                    unique(unlist(lav_partable_vnames(tmp.lav,
                                                  type = "ov",
                                                  level = tmp.level.values[l])))
                }
            } # levels
        } # groups
    } else {
        # perhaps model is already a parameter table
        nlevels <- lav_partable_nlevels(FLAT)
        if(nlevels > 1L) {

            # check for cluster argument (only if we have data)
            if(!is.null(data) && is.null(cluster)) {
                stop("lavaan ERROR: cluster argument is missing.")
            }

            ngroups <- lav_partable_ngroups(FLAT)
            ov.names.l <- vector("list", length = ngroups)
            for(g in 1:ngroups) {
                # note: lavNames() will return a list if any level:
                ov.names.l[[g]] <- lavNames(FLAT, "ov", group = g)
            }
        } else {
            # no level: in model syntax
            ov.names.l <- list()
        }
    }

    # sanity check ordered argument (just in case, add lhs variables names)
    ordered <- unique(c(ordered, lavNames(FLAT, "ov.ord")))



    #######################
    #### 2. lavoptions ####
    #######################
    if(!is.null(slotOptions)) {
        lavoptions <- slotOptions

        # but what if other 'options' are given anyway (eg 'start = ')?
        # give a warning!
        if(length(dotdotdot) > 0L) {
            dot.names <- names(dotdotdot)
            op.idx <- which(dot.names %in% names(slotOptions))
            warning("lavaan WARNING: the following argument(s) override(s) the options in slotOptions:\n\t\t", paste(dot.names[op.idx], collapse = " "))
            lavoptions[ dot.names[op.idx] ] <- dotdotdot[ op.idx ]
        }
    } else {

        # load default options
        opt <- lav_options_default()

        # catch unknown options
        ok.names <- names(opt)
        dot.names <- names(dotdotdot)
        wrong.idx <- which(!dot.names %in% ok.names)
        if(length(wrong.idx) > 0L) {
            idx <- wrong.idx[1L] # only show first one
            # stop or warning?? stop for now (there could be more)
            stop("lavaan ERROR: unknown argument `", dot.names[idx],"'")
        }

        # modifyList
        opt <- modifyList(opt, dotdotdot)

        # categorical mode?
        if(any(FLAT$op == "|")) {
            opt$categorical <- TRUE
        } else if(!is.null(data) && length(ordered) > 0L) {
            opt$categorical <- TRUE
        } else if(!is.null(sample.th)) {
            opt$categorical <- TRUE
        } else if(is.data.frame(data) &&
            lav_dataframe_check_ordered(frame = data, ov.names = ov.names.y)) {
            opt$categorical <- TRUE
        } else {
            opt$categorical <- FALSE
        }

        # clustered?
        if(length(cluster) > 0L) {
            opt$clustered <- TRUE
        } else {
            opt$clustered <- FALSE
        }

        # multilevel?
        if(length(ov.names.l) > 0L && length(ov.names.l[[1]]) > 1L) {
            opt$multilevel <- TRUE
        } else {
            opt$multilevel <- FALSE
        }

        # sampling weights? force MLR
        if(!is.null(sampling.weights)) {
            opt$estimator <- "MLR"
        }

        # constraints
        if(nchar(constraints) > 0L && opt$estimator %in% c("ML")) {
            opt$information <- "observed"
        }

        # meanstructure
        if(any(FLAT$op == "~1") || !is.null(sample.mean)) {
            opt$meanstructure <- TRUE
        }
        if(!is.null(group) && is.null(dotdotdot$meanstructure)) {
            opt$meanstructure <- TRUE
        }

        # conditional.x
        if( (is.list(ov.names.x) &&
             sum(sapply(ov.names.x, FUN = length)) == 0L) ||
            (is.character(ov.names.x) && length(ov.names.x) == 0L) ) {
            # if explicitly set to TRUE, give warning
            if(is.logical(dotdotdot$conditional.x) && dotdotdot$conditional.x) {
                warning("lavaan WARNING: no exogenous covariates; conditional.x will be set to FALSE")
            }
            opt$conditional.x <- FALSE
        }

        # fixed.x
        if( (is.list(ov.names.x) &&
             sum(sapply(ov.names.x, FUN = length)) == 0L) ||
            (is.character(ov.names.x) && length(ov.names.x) == 0L) ) {
            # if explicitly set to TRUE, give warning
            if(is.logical(dotdotdot$fixed.x) && dotdotdot$fixed.x) {
                # ok, we respect this: keep fixed.x = TRUE
            } else {
                opt$fixed.x <- FALSE
            }
        }

        # fill in remaining "default" values
        lavoptions <- lav_options_set(opt)
    }
    timing$Options <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]


    # fixed.x = FALSE? set ov.names.x = character(0L)
    # new in 0.6-1
    if(!lavoptions$fixed.x) {
        ov.names.x <- character(0L)
    }


    #####################
    #### 3. lavdata  ####
    #####################
    if(!is.null(slotData)) {
        lavdata <- slotData
    } else {
        # FIXME: ov.names should always contain both y and x!
        OV.NAMES <- if(lavoptions$conditional.x) { ov.names.y } else {ov.names}
        lavdata <- lavData(data             = data,
                           group            = group,
                           cluster          = cluster,
                           ov.names         = OV.NAMES,
                           ov.names.x       = ov.names.x,
                           ov.names.l       = ov.names.l,
                           ordered          = ordered,
                           sampling.weights = sampling.weights,
                           sample.cov       = sample.cov,
                           sample.mean      = sample.mean,
                           sample.th        = sample.th,
                           sample.nobs      = sample.nobs,
                           lavoptions       = lavoptions)
    }
    # what have we learned from the data?
    if(lavdata@data.type == "none") {
        lavoptions$do.fit <- FALSE
        lavoptions$start  <- "simple"
        lavoptions$se     <- "none"
        lavoptions$test   <- "none"
    } else if(lavdata@data.type == "moment") {

        # check user-specified options first
        if(!is.null(dotdotdot$estimator)) {
            if(dotdotdot$estimator %in%
                  c("MLM", "MLMV", "MLR", "MLR", "ULSM", "ULSMV", "ULSMVS") &&
               is.null(NACOV)) {
                stop("lavaan ERROR: estimator ", dotdotdot$estimator,
                 " requires full data or user-provided NACOV")
            } else if(dotdotdot$estimator %in%
                          c("WLS", "WLSM", "WLSMV", "WLSMVS", "DWLS") &&
                      is.null(WLS.V)) {
                stop("lavaan ERROR: estimator ", dotdotdot$estimator,
                 " requires full data or user-provided WLS.V and NACOV")
            }
        }

        # catch here some options that will not work with moments
        if(lavoptions$se == "bootstrap") {
            stop("lavaan ERROR: bootstrapping requires full data")
        }
        # more needed?
    }
    timing$Data <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
    if(lavoptions$debug) {
        print(str(lavdata))
    }

    # if lavdata@nlevels > 1L, adapt start option (for now)
    # until we figure out how to handle groups+blocks
    #if(lavdata@nlevels > 1L) {
    #   lavoptions$start <- "simple"
    #}








    ########################
    #### 4. lavpartable ####
    ########################
    if(!is.null(slotParTable)) {
        lavpartable <- slotParTable
    } else if(is.character(model) ||
              inherits(model, "formula")) {
        # check FLAT before we proceed
        if(lavoptions$debug) {
            print(as.data.frame(FLAT))
        }
        # catch ~~ of fixed.x covariates if fixed.x = TRUE
        if(lavoptions$fixed.x) {
            tmp <- try(vnames(FLAT, type = "ov.x", ov.x.fatal = TRUE),
                       silent = TRUE)
            if(inherits(tmp, "try-error")) {
                warning("lavaan WARNING: syntax contains parameters involving exogenous covariates; switching to fixed.x = FALSE")
                lavoptions$fixed.x <- FALSE
            }
        }
        if(lavoptions$conditional.x) {
            tmp <- vnames(FLAT, type = "ov.x", ov.x.fatal = TRUE)
        }

        lavpartable <-
            lavaanify(model            = FLAT,
                      constraints      = constraints,
                      varTable         = lavdata@ov,
                      ngroups          = lavdata@ngroups,

                      meanstructure    = lavoptions$meanstructure,
                      int.ov.free      = lavoptions$int.ov.free,
                      int.lv.free      = lavoptions$int.lv.free,
                      orthogonal       = lavoptions$orthogonal,
                      conditional.x    = lavoptions$conditional.x,
                      fixed.x          = lavoptions$fixed.x,
                      std.lv           = lavoptions$std.lv,
                      parameterization = lavoptions$parameterization,
                      auto.fix.first   = lavoptions$auto.fix.first,
                      auto.fix.single  = lavoptions$auto.fix.single,
                      auto.var         = lavoptions$auto.var,
                      auto.cov.lv.x    = lavoptions$auto.cov.lv.x,
                      auto.cov.y       = lavoptions$auto.cov.y,
                      auto.th          = lavoptions$auto.th,
                      auto.delta       = lavoptions$auto.delta,
                      group.equal      = lavoptions$group.equal,
                      group.partial    = lavoptions$group.partial,
                      group.w.free     = lavoptions$group.w.free,
                      debug            = lavoptions$debug,
                      warn             = lavoptions$warn,

                      as.data.frame.   = FALSE)

    } else if(inherits(model, "lavaan")) {
        lavpartable <- parTable(model)
    } else if(is.list(model)) {
        # we already checked this when creating FLAT
        # but we may need to complete it
        lavpartable <- as.list(FLAT) # in case model is a data.frame
        # complete table
        lavpartable <- lav_partable_complete(lavpartable)
    } else {
        stop("lavaan ERROR: model [type = ", class(model),
             "] is not of type character or list")
    }
    if(lavoptions$debug) {
        print(as.data.frame(lavpartable))
    }

    # at this point, we should check if the partable is complete
    # or not; this is especially relevant if the lavaan() function
    # was used, but the user has forgotten some variances/intercepts...
    junk <- lav_partable_check(lavpartable,
                               categorical = lavoptions$categorical,
                               warn = TRUE)

    # for EM only (for now), force fixed-to-zero (residual) variances
    # to be slightly larger than zero
    if(lavoptions$optim.method == "em") {
        zero.var.idx <- which(lavpartable$op == "~~" &
                              lavpartable$lhs == lavpartable$rhs &
                              lavpartable$free == 0L &
                              lavpartable$ustart == 0)
        if(length(zero.var.idx) > 0L) {
            lavpartable$ustart[zero.var.idx] <- lavoptions$em.zerovar.offset
        }
    }

    # 4b. get partable attributes
    lavpta <- lav_partable_attributes(lavpartable)
    timing$ParTable <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]




    ###########################
    #### 5. lavsamplestats ####
    ##########################
    if(!is.null(slotSampleStats)) {
        lavsamplestats <- slotSampleStats
    } else if(lavdata@data.type == "full") {
        lavsamplestats <- lav_samplestats_from_data(
                       lavdata       = lavdata,
                       missing       = lavoptions$missing,
                       rescale       =
                           (lavoptions$estimator %in% c("ML","REML","NTRLS") &&
                            lavoptions$likelihood == "normal"),
                       estimator     = lavoptions$estimator,
                       mimic         = lavoptions$mimic,
                       meanstructure = lavoptions$meanstructure,
                       conditional.x = lavoptions$conditional.x,
                       fixed.x       = lavoptions$fixed.x,
                       group.w.free  = lavoptions$group.w.free,
                       missing.h1    = (lavoptions$missing != "listwise"),
                       WLS.V             = WLS.V,
                       NACOV             = NACOV,
                       gamma.n.minus.one = lavoptions$gamma.n.minus.one,
                       se                = lavoptions$se,
                       information       = lavoptions$information,
                       ridge             = lavoptions$ridge,
                       optim.method      = lavoptions$optim.method.cor,
                       zero.add          = lavoptions$zero.add,
                       zero.keep.margins = lavoptions$zero.keep.margins,
                       zero.cell.warn    = lavoptions$zero.cell.warn,
                       debug             = lavoptions$debug,
                       verbose           = lavoptions$verbose)

    } else if(lavdata@data.type == "moment") {
        lavsamplestats <- lav_samplestats_from_moments(
                           sample.cov    = sample.cov,
                           sample.mean   = sample.mean,
                           sample.th     = sample.th,
                           sample.nobs   = sample.nobs,
                           ov.names      = ov.names,
                           ov.names.x    = ov.names.x,
                           estimator     = lavoptions$estimator,
                           mimic         = lavoptions$mimic,
                           group.w.free  = lavoptions$group.w.free,
                           WLS.V         = WLS.V,
                           NACOV         = NACOV,
                           ridge         = lavoptions$ridge,
                           rescale       = lavoptions$sample.cov.rescale)
    } else {
        # no data
        lavsamplestats <- new("lavSampleStats", ngroups=lavdata@ngroups,
                                 nobs=as.list(rep(0L, lavdata@ngroups)),
                                 #cov.x=vector("list",length=lavdata@ngroups),
                                 #mean.x=vector("list",length=lavdata@ngroups),
                                 th.idx=lavpta$th.idx,
                                 missing.flag=FALSE)
    }
    timing$SampleStats <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]
    if(lavoptions$debug) {
        print(str(lavsamplestats))
    }



    ###################
    #### 5b. lavh1 ####
    ###################
    if(!is.null(sloth1)) {
        lavh1 <- sloth1
    } else {
        lavh1 <- list()
        if(is.logical(lavoptions$h1) && lavoptions$h1) {
            if(length(lavsamplestats@ntotal) > 0L) { # lavsamplestats filled in

                # implied h1 statistics
                out <- lav_h1_implied_logl(lavdata        = lavdata,
                                           lavsamplestats = lavsamplestats,
                                           lavoptions     = lavoptions)
                h1.implied      <- out$implied
                h1.loglik       <- out$logl$loglik
                h1.loglik.group <- out$logl$loglik.group

                # collect in h1 list
                lavh1 <- list(implied      = h1.implied,
                              loglik       = h1.loglik,
                              loglik.group = h1.loglik.group)
            } else {
                # do nothing for now
            }
        } else {
            if(!is.logical(lavoptions$h1)) {
                stop("lavaan ERROR: argument `h1' must be logical (for now)")
            }
            # TODO: allow h1 to be either a model syntax, a parameter table,
            # or a fitted lavaan object
        }
    }
    timing$h1 <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]






    #####################
    #### 6. lavstart ####
    #####################
    if(!is.null(slotModel)) {
        lavmodel <- slotModel
        # FIXME
        #lavaanStart <- lav_model_get_parameters(lavmodel, type="user")
        #lavpartable$start <- lavaanStart
        timing$start <- (proc.time()[3] - start.time)
        start.time <- proc.time()[3]
        timing$Model <- (proc.time()[3] - start.time)
        start.time <- proc.time()[3]
    } else {
        # check if we have provide a full parameter table as model= input
        if(!is.null(lavpartable$est) && lavoptions$start == "default") {
            # check if all 'est' values look ok
            # this is not the case, eg, if partables have been merged eg, as
            # in semTools' auxiliary() function

            # check for zero free variances and NA values
            zero.idx <- which(lavpartable$free > 0L &
                              lavpartable$op == "~~" &
                              lavpartable$lhs == lavpartable$rhs &
                              lavpartable$est == 0)

            if(length(zero.idx) > 0L || any(is.na(lavpartable$est))) {
                lavpartable$start <- lav_start(start.method = lavoptions$start,
                                       lavpartable     = lavpartable,
                                       lavsamplestats  = lavsamplestats,
                                       model.type   = lavoptions$model.type,
                                       mimic        = lavoptions$mimic,
                                       debug        = lavoptions$debug)
            } else {
                lavpartable$start <- lavpartable$est
            }
        } else {
            START <- lav_start(start.method   = lavoptions$start,
                               lavpartable    = lavpartable,
                               lavsamplestats = lavsamplestats,
                               model.type     = lavoptions$model.type,
                               mimic          = lavoptions$mimic,
                               debug          = lavoptions$debug)

            # sanity check
            if(!is.null(lavoptions$check.start) && lavoptions$check.start) {
                START <- lav_start_check_cov(lavpartable = lavpartable,
                                             start       = START)
            }

            lavpartable$start <- START
        }


        timing$start <- (proc.time()[3] - start.time)
        start.time <- proc.time()[3]



    #####################
    #### 7. lavmodel ####
    #####################
        lavmodel <- lav_model(lavpartable      = lavpartable,
                              lavoptions       = lavoptions,
                              th.idx           = lavsamplestats@th.idx,
                              cov.x            = lavsamplestats@cov.x,
                              mean.x           = lavsamplestats@mean.x)
        timing$Model <- (proc.time()[3] - start.time)
        start.time <- proc.time()[3]

        # if no data, call lav_model_set_parameters once (for categorical case)
        if(lavdata@data.type == "none" && lavmodel@categorical) {
            lavmodel <- lav_model_set_parameters(lavmodel = lavmodel,
                            x = lav_model_get_parameters(lavmodel))
            # re-adjust parameter table
            lavpartable$start <- lav_model_get_parameters(lavmodel, type="user")

            # check/warn if theta/delta values make sense
            if(!all(lavpartable$start == lavpartable$ustart)) {
                if(lavmodel@parameterization == "delta") {
                    # did the user specify theta values?
                    user.var.idx <- which(lavpartable$op == "~~" &
                           lavpartable$lhs == lavpartable$rhs &
                           lavpartable$lhs %in% unlist(lavpta$vnames$ov.ord) &
                           lavpartable$user == 1L)
                    if(length(user.var.idx)) {
                        warning("lavaan WARNING: ",
              "variance (theta) values for categorical variables are ignored",
              "\n\t\t  if parameterization = \"delta\"!")
                    }
                } else if(lavmodel@parameterization == "theta") {
                    # did the user specify theta values?
                    user.delta.idx <- which(lavpartable$op == "~*~" &
                            lavpartable$lhs == lavpartable$rhs &
                            lavpartable$lhs %in% unlist(lavpta$vnames$ov.ord) &
                            lavpartable$user == 1L)
                    if(length(user.delta.idx)) {
                        warning("lavaan WARNING: ",
              "scaling (~*~) values for categorical variables are ignored",
              "\n\t\t  if parameterization = \"theta\"!")
                    }
                }
            }
        }

    } # slotModel

    #######################
    #### 7. seat belts #### # TODO: bounds,...
    #######################
    #if(lavoptions$safe.ov.var.ub) {
    #
    #}
    #if(lavoptions$safe.ov.var.lb) {
    #}







    #####################
    #### 8. lavcache ####
    #####################
    if(!is.null(slotCache)) {
        lavcache <- slotCache
    } else {
        # prepare cache -- stuff needed for estimation, but also post-estimation
        lavcache <- vector("list", length=lavdata@ngroups)

        # ov.types? (for PML check)
        ov.types <- lavdata@ov$type
        if(lavmodel@conditional.x && sum(lavmodel@nexo) > 0L) {
            # remove ov.x
            ov.x.idx <- unlist(lavpta$vidx$ov.x)
            ov.types <- ov.types[-ov.x.idx]
        }

        if(lavoptions$estimator == "PML" && all(ov.types == "ordered")) {
            TH <- computeTH(lavmodel)
            BI <- lav_tables_pairwise_freq_cell(lavdata)

            # handle option missing = "available.cases"
            if(lavoptions$missing == "available.cases" ||
               lavoptions$missing == "doubly.robust") {
                UNI <- lav_tables_univariate_freq_cell(lavdata)
            }

            # checks for missing = "double.robust"
            if (lavoptions$missing == "doubly.robust") {
                # check whether the probabilities pairwiseProbGivObs and
                # univariateProbGivObs are given by the user
                if(is.null(lavoptions$control$pairwiseProbGivObs)) {
                    stop("lavaan ERROR: could not find `pairwiseProbGivObs' in control() list")
                }
                if(is.null(lavoptions$control$univariateProbGivObs)) {
                    stop("lavaan ERROR: could not find `univariateProbGivObs' in control() list")
                }
            }

            for(g in 1:lavdata@ngroups) {

                if(is.null(BI$group) || max(BI$group) == 1L) {
                    bifreq <- BI$obs.freq
                    binobs  <- BI$nobs
                } else {
                    idx <- which(BI$group == g)
                    bifreq <- BI$obs.freq[idx]
                    binobs  <- BI$nobs[idx]
                }
                LONG  <- LongVecInd(no.x               = ncol(lavdata@X[[g]]),
                                    all.thres          = TH[[g]],
                                    index.var.of.thres = lavmodel@th.idx[[g]])
                lavcache[[g]] <- list(bifreq = bifreq,
                                      nobs   = binobs,
                                      LONG   = LONG)

                # available cases
                if(lavoptions$missing == "available.cases" ||
                   lavoptions$missing == "doubly.robust") {
                    if(is.null(UNI$group) || max(UNI$group) == 1L) {
                        unifreq <- UNI$obs.freq
                        uninobs <- UNI$nobs
                    } else {
                        idx <- which(UNI$group == g)
                        unifreq <- UNI$obs.freq[idx]
                        uninobs <- UNI$nobs[idx]
                    }
                    lavcache[[g]]$unifreq <- unifreq
                    lavcache[[g]]$uninobs <- uninobs

                    uniweights.casewise <- rowSums( is.na( lavdata@X[[g]] ) )
                    lavcache[[g]]$uniweights.casewise <- uniweights.casewise

                    #weights per response category per variable in the same
                    # order as unifreq; i.e. w_ia, i=1,...,p, (p variables),
                    # a=1,...,Ci, (Ci response categories for variable i),
                    # a running faster than i
                    tmp.uniweights <- apply(lavdata@X[[g]], 2,
                        function(x){
                            tapply(uniweights.casewise, as.factor(x), sum,
                                   na.rm=TRUE) } )
                    if( is.matrix(tmp.uniweights) ) {
                        lavcache[[g]]$uniweights <- c(tmp.uniweights)
                    }
                    if( is.list(tmp.uniweights) ) {
                        lavcache[[g]]$uniweights <- unlist(tmp.uniweights)
                    }
                } # "available.cases" or "double.robust"

                # doubly.robust only
                if (lavoptions$missing == "doubly.robust") {

                    # add the provided by the user probabilities
                    # pairwiseProbGivObs and univariateProbGivObs in Cache
                    lavcache[[g]]$pairwiseProbGivObs <-
                        lavoptions$control$pairwiseProbGivObs[[g]]
                    lavcache[[g]]$univariateProbGivObs <-
                        lavoptions$control$univariateProbGivObs[[g]]
                    # compute different indices vectors that will help to do
                    # calculations
                    ind.vec <- as.data.frame(LONG[1:5] )
                    ind.vec <-
                        ind.vec[ ((ind.vec$index.thres.var1.of.pair!=0) &
                                  (ind.vec$index.thres.var2.of.pair!=0)) , ]
                    idx.cat.y1 <- ind.vec$index.thres.var1.of.pair
                    idx.cat.y2 <- ind.vec$index.thres.var2.of.pair
                    idx.pairs  <- ind.vec$index.pairs.extended
                    lavcache[[g]]$idx.pairs <- idx.pairs

                    idx.cat.y1.split <- split(idx.cat.y1, idx.pairs)
                    idx.cat.y2.split <- split(idx.cat.y2, idx.pairs)
                    lavcache[[g]]$idx.cat.y1.split <- idx.cat.y1.split
                    lavcache[[g]]$idx.cat.y2.split <- idx.cat.y2.split

                    # generate the variables, categories indices vector which
                    # keep track to which variables and categories the
                    # elements of vector probY1Gy2 refer to
                    nlev <- lavdata@ov$nlev
                    nvar <- length(nlev)

                    idx.var.matrix <- matrix(1:nvar, nrow=nvar, ncol=nvar)
                    idx.diag <- diag( matrix(1:(nvar*nvar), nrow=nvar,
                                             ncol=nvar) )
                    idx.Y1Gy2.matrix <- rbind(t(idx.var.matrix)[-idx.diag],
                                                idx.var.matrix [-idx.diag])
                    no.pairs.Y1Gy2 <- ncol(idx.Y1Gy2.matrix)
                    idx.cat.Y1 <- unlist(lapply(1:no.pairs.Y1Gy2, function(x) {
                        rep(     1:nlev[ idx.Y1Gy2.matrix[1,x] ],
                            times= nlev[ idx.Y1Gy2.matrix[2,x] ]   )} ) )
                    idx.cat.Gy2 <- unlist(lapply(1:no.pairs.Y1Gy2, function(x) {
                        rep(     1:nlev[ idx.Y1Gy2.matrix[2,x] ],
                             each= nlev[ idx.Y1Gy2.matrix[1,x] ]   )} ) )
                    dim.pairs <- unlist(lapply(1:no.pairs.Y1Gy2, function(x) {
                        nlev[ idx.Y1Gy2.matrix[1,x] ] *
                        nlev[ idx.Y1Gy2.matrix[2,x] ] }) )
                    idx.Y1 <- unlist( mapply(rep, idx.Y1Gy2.matrix[1,],
                                             each=dim.pairs) )
                    idx.Gy2 <- unlist( mapply(rep, idx.Y1Gy2.matrix[2,],
                                              each=dim.pairs) )

                    lavcache[[g]]$idx.Y1      <- idx.Y1
                    lavcache[[g]]$idx.Gy2     <- idx.Gy2
                    lavcache[[g]]$idx.cat.Y1  <- idx.cat.Y1
                    lavcache[[g]]$idx.cat.Gy2 <- idx.cat.Gy2

                    # the vector below keeps track of the variable each column
                    # of the matrix univariateProbGivObs refers to
                    lavcache[[g]]$id.uniPrGivObs <-
                        sort( c( unique(lavmodel@th.idx[[g]]) ,
                                        lavmodel@th.idx[[g]] ) )
                } # doubly.robust



            } # g
        }
        # copy response patterns to cache -- FIXME!! (data not included
        # in Model only functions)
        if(lavdata@data.type == "full" && !is.null(lavdata@Rp[[1L]])) {
            for(g in 1:lavdata@ngroups) {
                lavcache[[g]]$pat <- lavdata@Rp[[g]]$pat
            }
        }
    }

    # If estimator = MML, store Gauss-Hermite nodes/weights
    if(lavoptions$estimator == "MML") {
        for(g in 1:lavdata@ngroups) {
            # count only the ones with non-normal indicators
            #nfac <- lavpta$nfac.nonnormal[[g]]
            nfac <- lavpta$nfac[[g]]
            lavcache[[g]]$GH <-
                lav_integration_gauss_hermite(n = lavoptions$integration.ngh,
                                              dnorm = TRUE,
                                              mean = 0, sd = 1,
                                              ndim = nfac)
            #lavcache[[g]]$DD <- lav_model_gradient_DD(lavmodel, group = g)
        }
    }

    timing$cache <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]








    ############################
    #### 10. est + lavoptim ####
    ############################
    x <- NULL
    if(lavoptions$do.fit && lavoptions$estimator != "none" &&
       lavmodel@nx.free > 0L) {

        if(lavoptions$optim.method == "em") {
            # multilevel only for now
            stopifnot(lavdata@nlevels > 1L)
            x <- lav_mvnorm_cluster_em_h0(lavsamplestats = lavsamplestats,
                               lavdata        = lavdata,
                               lavimplied     = NULL,
                               lavpartable    = lavpartable,
                               lavmodel       = lavmodel,
                               lavoptions     = lavoptions,
                               verbose        = lavoptions$verbose,
                               fx.tol         = lavoptions$em.fx.tol,
                               dx.tol         = lavoptions$em.dx.tol,
                               max.iter       = lavoptions$em.iter.max)
        } else {
            x <- lav_model_estimate(lavmodel        = lavmodel,
                                    lavpartable     = lavpartable,
                                    lavsamplestats  = lavsamplestats,
                                    lavdata         = lavdata,
                                    lavoptions      = lavoptions,
                                    lavcache        = lavcache)
        }
        lavmodel <- lav_model_set_parameters(lavmodel, x = x)
        # store parameters in @ParTable$est
        lavpartable$est <- lav_model_get_parameters(lavmodel = lavmodel,
                                                    type = "user", extra = TRUE)

        if(!is.null(attr(x, "con.jac")))
            lavmodel@con.jac <- attr(x, "con.jac")
        if(!is.null(attr(x, "con.lambda")))
            lavmodel@con.lambda <- attr(x, "con.lambda")
        # check if model has converged or not
        if(!attr(x, "converged") && lavoptions$warn) {
           warning("lavaan WARNING: the optimizer warns that a solution has NOT been found!")
        }
    } else {
        x <- numeric(0L)
        attr(x, "iterations") <- 0L; attr(x, "converged") <- FALSE
        attr(x, "control") <- lavoptions$control
        attr(x, "fx") <-
            lav_model_objective(lavmodel = lavmodel,
                lavsamplestats = lavsamplestats, lavdata = lavdata,
                lavcache = lavcache)

        lavpartable$est <- lavpartable$start
    }

    # should we fake/force convergence? (eg. to enforce the
    # computation of a test statistic)
    if(lavoptions$optim.force.converged) {
        attr(x, "converged") <- TRUE
    }

    # store optimization info in lavoptim
    lavoptim <- list()
    x2 <- x; attributes(x2) <- NULL
    lavoptim$x <- x2
    lavoptim$npar <- length(x)
    lavoptim$iterations <- attr(x, "iterations")
    lavoptim$converged  <- attr(x, "converged")
    lavoptim$parscale   <- attr(x, "parscale")
    fx.copy <- fx <- attr(x, "fx"); attributes(fx) <- NULL
    lavoptim$fx         <- fx
    lavoptim$fx.group   <- attr(fx.copy, "fx.group")
    if(!is.null(attr(fx.copy, "logl.group"))) {
        lavoptim$logl.group <- attr(fx.copy, "logl.group")
        lavoptim$logl       <- sum(lavoptim$logl.group)
    } else {
        lavoptim$logl.group <- as.numeric(NA)
        lavoptim$logl       <- as.numeric(NA)
    }
    lavoptim$control        <- attr(x, "control")

    timing$optim <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]








    ####################################
    #### 11. lavimplied + lavloglik ####
    ####################################
    lavimplied <- list()
    if(lavoptions$implied) {
         lavimplied <- lav_model_implied(lavmodel)
    }

    lavloglik <- list()
    if(lavoptions$loglik) {
         lavloglik <- lav_model_loglik(lavdata        = lavdata,
                                       lavsamplestats = lavsamplestats,
                                       lavimplied     = lavimplied,
                                       lavmodel       = lavmodel,
                                       lavoptions     = lavoptions)
    }

    timing$implied <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]












    ###############################
    #### 12. lavvcov + lavboot ####
    ###############################
    VCOV <- NULL
    if(lavoptions$se != "none" && lavoptions$se != "external" &&
       lavoptions$se != "twostep" &&
       lavmodel@nx.free > 0L && attr(x, "converged")) {
        if(lavoptions$verbose) {
            cat("Computing VCOV for se =", lavoptions$se, "...")
        }
        VCOV <- lav_model_vcov(lavmodel        = lavmodel,
                               lavsamplestats  = lavsamplestats,
                               lavoptions      = lavoptions,
                               lavdata         = lavdata,
                               lavpartable     = lavpartable,
                               lavcache        = lavcache,
                               lavimplied      = lavimplied,
                               lavh1           = lavh1)
        if(lavoptions$verbose) {
            cat(" done.\n")
        }
    }

    # extract bootstrap results (if any)
    if(!is.null(attr(VCOV, "BOOT.COEF"))) {
        lavboot <- list()
        lavboot$coef <- attr(VCOV, "BOOT.COEF")
    } else {
        lavboot <- list()
    }

    # store VCOV in vcov
    # strip all attributes but 'dim'
    tmp.attr <- attributes(VCOV)
    VCOV1 <- VCOV
    attributes(VCOV1) <- tmp.attr["dim"]
    lavvcov <- list(se = lavoptions$se, information = lavoptions$information,
                    vcov = VCOV1)

    # store se in partable
    if(lavoptions$se != "external" && lavoptions$se != "twostep") {
        lavpartable$se <- lav_model_vcov_se(lavmodel = lavmodel,
                                            lavpartable = lavpartable,
                                            VCOV = VCOV,
                                            BOOT = lavboot$coef)
    } else {
        if(is.null(lavpartable$se)) {
            lavpartable$se <- lav_model_vcov_se(lavmodel = lavmodel,
                                                lavpartable = lavpartable,
                                                VCOV = NULL, BOOT = NULL)
            warning("lavaan WARNING: se = \"external\" but parameter table does not contain a `se' column")
        }
    }

    timing$vcov <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]









    #####################
    #### 13. lavtest ####
    #####################
    TEST <- NULL
    if(lavoptions$test != "none" && attr(x, "converged")) {
        if(lavoptions$verbose) {
            cat("Computing TEST for test =", lavoptions$test, "...")
        }
        TEST <- lav_model_test(lavmodel            = lavmodel,
                               lavpartable         = lavpartable,
                               lavsamplestats      = lavsamplestats,
                               lavimplied          = lavimplied,
                               lavh1               = lavh1,
                               lavoptions          = lavoptions,
                               x                   = x,
                               VCOV                = VCOV,
                               lavdata             = lavdata,
                               lavcache            = lavcache,
                               lavloglik           = lavloglik)
        if(lavoptions$verbose) {
            cat(" done.\n")
        }
    } else {
        TEST <- list(list(test="none", stat=NA,
                     stat.group=rep(NA, lavdata@ngroups), df=NA,
                     refdistr="unknown", pvalue=NA))
    }

    # store test in lavtest
    lavtest <- TEST

    timing$test <- (proc.time()[3] - start.time)
    start.time <- proc.time()[3]










    ####################
    #### 14. lavfit ####
    ####################
    lavfit <- lav_model_fit(lavpartable = lavpartable,
                            lavmodel    = lavmodel,
                            x           = x,
                            VCOV        = VCOV,
                            TEST        = TEST)
    timing$total <- (proc.time()[3] - start.time0)




    ######################
    #### 15. baseline ####
    ######################
    lavbaseline <- list()
    if(is.logical(lavoptions$baseline) && lavoptions$baseline) {

    }


    ####################
    #### 16. lavaan ####
    ####################
    lavaan <- new("lavaan",
                  version      = as.character(packageVersion("lavaan")),
                  call         = mc,                  # match.call
                  timing       = timing,              # list
                  Options      = lavoptions,          # list
                  ParTable     = lavpartable,         # list
                  pta          = lavpta,              # list
                  Data         = lavdata,             # S4 class
                  SampleStats  = lavsamplestats,      # S4 class
                  Model        = lavmodel,            # S4 class
                  Cache        = lavcache,            # list
                  Fit          = lavfit,              # S4 class
                  boot         = lavboot,             # list
                  optim        = lavoptim,            # list
                  implied      = lavimplied,          # list
                  loglik       = lavloglik,           # list
                  vcov         = lavvcov,             # list
                  test         = lavtest,             # list
                  h1           = lavh1,               # list
                  baseline     = list(),              # list
                  external     = list()               # empty list
                 )



    # post-fitting check of parameters
    if(!is.null(lavoptions$check.post) && lavoptions$check.post &&
       lavTech(lavaan, "converged")) {
        lavInspect(lavaan, "post.check")
    }

    # new in 0.6-2
    # FIXME: not scale independent (should use solve(Hessian) %*% g)
    # but Hessian is not always available (or expensive to compute)
    hasExplicitConstraints <- FALSE
    if(is.character(constraints) && nchar(constraints) > 0L) {
        hasExplicitConstraints <- TRUE
    }
    hasNonLinearEqConstraints <- FALSE
    if(length(lavmodel@ceq.nonlinear.idx) > 0L) {
        hasNonLinearEqConstraints <- TRUE
    }
    hasIneqConstraints <- FALSE
    if(length(lavmodel@cin.linear.idx) > 0L ||
       length(lavmodel@cin.nonlinear.idx) > 0L) {
        hasIneqConstraints <- TRUE
    }
    if(!is.null(lavoptions$check.gradient) && lavoptions$check.gradient &&
       lavTech(lavaan, "converged") &&
       !hasExplicitConstraints      &&
       !hasNonLinearEqConstraints   &&
       !hasIneqConstraints) {
        grad <- lavInspect(lavaan, "optim.gradient")
        large.idx <- which(abs(grad) > 0.001)  # better 0.0001?
        if(length(large.idx) > 0L) {
            warning(
  "lavaan WARNING: not all elements of the gradient are (near) zero;\n",
"                  the optimizer may not have found a local solution;\n",
"                  use lavInspect(fit, \"optim.gradient\") to investigate")
        }
    }

    ########################
    #### post: baseline ####
    ########################
    #lavbaseline <- list()
    #if(is.logical(lavoptions$baseline) && lavoptions$baseline) {
    #    fit.indep <- try(lav_object_independence(lavaan), silent = TRUE)
    #    X2.null <- df.null <- as.numeric(NA)
    #    X2.null.scaled <- df.null.scaled <- as.numeric(NA)
    #    if(inherits(fit.indep, "try-error")) {
    #        warning("lavaan WARNING: estimation of the baseline model failed.")
    #    } else {
    #        X2.null <- fit.indep@test[[1]]$stat
    #        df.null <- fit.indep@test[[1]]$df
    #        if(length(fit.indep@test) > 1L) {
    #            X2.null.scaled <- fit.indep@test[[2]]$stat
    #            df.null.scaled <- fit.indep@test[[2]]$df
    #        }
    #    }
    #
    #    # store in list
    #    lavbaseline <- list(X2.null = X2.null,
    #                        df.null = df.null,
    #                        X2.null.scaled = X2.null.scaled,
    #                        df.null.scaled = df.null.scaled)
    #
    #    # add to lavaan object
    #    lavaan@baseline <- lavbaseline
    #}



    lavaan
}



# cfa + sem
cfa <- sem <- function(# user-specified model: can be syntax, parameter Table
                       model              = NULL,
                       # data (second argument, most used)
                       data               = NULL,

                       # variable information
                       ordered            = NULL,

                       # sampling weights
                       sampling.weights   = NULL,

                       # summary data
                       sample.cov         = NULL,
                       sample.mean        = NULL,
                       sample.th          = NULL,
                       sample.nobs        = NULL,

                       # multiple groups?
                       group              = NULL,

                       # multiple levels?
                       cluster            = NULL,

                       # constraints
                       constraints        = '',

                       # user-specified variance matrices
                       WLS.V              = NULL,
                       NACOV              = NULL,

                       # options (dotdotdot)
                       ...) {

    mc <- match.call(expand.dots = TRUE)

    # set model.type
    mc$model.type      = as.character( mc[[1L]] )
    if(length(mc$model.type) == 3L) {
        mc$model.type <- mc$model.type[3L]
    }

    dotdotdot <- list(...)
    if(!is.null(dotdotdot$std.lv)) {
        std.lv <- dotdotdot$std.lv
    } else {
        std.lv <- FALSE
    }

    # default options for sem/cfa call
    mc$int.ov.free     = TRUE
    mc$int.lv.free     = FALSE
    mc$auto.fix.first  = !std.lv
    mc$auto.fix.single = TRUE
    mc$auto.var        = TRUE
    mc$auto.cov.lv.x   = TRUE
    mc$auto.cov.y      = TRUE
    mc$auto.th         = TRUE
    mc$auto.delta      = TRUE

    # call mother function
    mc[[1L]] <- quote(lavaan::lavaan)
    eval(mc, parent.frame())
}

# simple growth models
growth <- function(# user-specified model: can be syntax, parameter Table
                   model              = NULL,
                   # data (second argument, most used)
                   data               = NULL,

                   # variable information
                   ordered            = NULL,

                   # sampling weights
                   sampling.weights   = NULL,

                   # summary data
                   sample.cov         = NULL,
                   sample.mean        = NULL,
                   sample.th          = NULL,
                   sample.nobs        = NULL,

                   # multiple groups?
                   group              = NULL,

                   # multiple levels?
                   cluster            = NULL,

                   # constraints
                   constraints        = '',

                   # user-specified variance matrices
                   WLS.V              = NULL,
                   NACOV              = NULL,

                   # options (dotdotdot)
                   ...) {

    mc <- match.call(expand.dots = TRUE)

    # set model.type
    mc$model.type      = as.character( mc[[1L]] )
    if(length(mc$model.type) == 3L) {
        mc$model.type <- mc$model.type[3L]
    }

    dotdotdot <- list(...)
    if(!is.null(dotdotdot$std.lv)) {
        std.lv <- dotdotdot$std.lv
    } else {
        std.lv <- FALSE
    }

    # default options for sem/cfa call
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

    # call mother function
    mc[[1L]] <- quote(lavaan::lavaan)
    eval(mc, parent.frame())
}
