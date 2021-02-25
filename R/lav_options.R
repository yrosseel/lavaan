# initial version YR 02/08/2010

# YR 28 Jan 2017: add lavOptions(), lav_options_default()

# public function
lavOptions <- function(x = NULL, default = NULL, mimic = "lavaan") {

    lavoptions <- lav_options_default(mimic = mimic)

    # selection only
    if(!is.null(x)) {
        if(is.character(x)) {
            # lower case only
            x <- tolower(x)

            # check if x is in names(lavoptions)
            not.ok <- which(!x %in% names(lavoptions))
            if(length(not.ok) > 0L) {
                # only warn if multiple options were requested
                if(length(x) > 1L) {
                    warning("lavaan WARNING: option `", x[not.ok],
                            "' not available")
                }
                x <- x[ -not.ok ]
            }

            # return requested option(s)
            if(length(x) == 0L) {
                return(default)
            } else {
                lavoptions[x]
            }
        } else {
            stop("lavaan ERROR: `x' must be a character string")
        }
    } else {
        lavoptions
    }
}

# set the default options (including unspecified values "default")
lav_options_default <- function(mimic = "lavaan") {

    opt <- list(model.type         = "sem",

                # global
                mimic              = "lavaan",

                # model modifiers
                meanstructure      = "default",
                int.ov.free        = FALSE,
                int.lv.free        = FALSE,
                conditional.x      = "default", # or FALSE?
                fixed.x            = "default", # or FALSE?
                orthogonal         = FALSE,
                orthogonal.x       = FALSE,
                orthogonal.y       = FALSE,
                std.lv             = FALSE,
                effect.coding      = FALSE,     # TRUE implies
                                                # c("loadings", "intercepts")
                parameterization   = "default",
                #ov.order           = "model",   # how to 'order' ov's in the
                #                                # partable

                auto.fix.first     = FALSE,
                auto.fix.single    = FALSE,
                auto.var           = FALSE,
                auto.cov.lv.x      = FALSE,
                auto.cov.y         = FALSE,
                auto.th            = FALSE,
                auto.delta         = FALSE,
                auto.efa           = FALSE,

                # seat belts
                safe.ov.var.ub     = FALSE,
                save.ov.var.lb     = FALSE,

                # rotation
                rotation           = "geomin",
                rotation.se        = "bordered", # "bordered" or "delta"
                rotation.args      = list(orthogonal     = FALSE,
                                          row.weights    = "default",
                                          std.ov         = TRUE,
                                          geomin.epsilon = 0.01,
                                          orthomax.gamma = 1,
                                          cf.gamma       = 0,
                                          oblimin.gamma  = 0,
                                          target         = matrix(0,0,0),
                                          target.mask    = matrix(0,0,0),
                                          rstarts        = 100L,
                                          algorithm      = "gpa",
                                          reflect        = TRUE,
                                          order.lv.by    = "index",
                                          gpa.tol        = 1e-05,
                                          tol            = 1e-08,
                                          warn           = FALSE,
                                          verbose        = FALSE,
                                          jac.init.rot   = TRUE,
                                          max.iter       = 10000L),

                # full data
                std.ov             = FALSE,
                missing            = "default",
                sampling.weights.normalization = "total",

                # summary data
                sample.cov.rescale = "default",
                ridge              = FALSE,
                ridge.x            = FALSE,
                ridge.constant     = "default",
                ridge.constant.x   = 1e-5,

                # multiple groups
                group.label        = NULL,
                group.equal        = '',
                group.partial      = '',
                group.w.free       = FALSE,

                # clusters
                level.label        = NULL,

                # estimation
                estimator              = "default",
                estimator.args         = list(),
                likelihood             = "default",
                link                   = "default",
                representation         = "default",
                do.fit                 = TRUE,
                bounds                 = "none", # new in 0.6-6

                # inference
                se                     = "default",
                test                   = "default",

                # information (se + test)
                information            = c("default",    "default"),
                h1.information         = c("structured", "structured"),
                observed.information   = c("hessian",    "default"),

                # information se only
                information.meat       = "default",
                h1.information.meat    = "default",

                # information for 'Omega' (yuan-benter test only)
                omega.information         = "default",
                omega.h1.information      = "default",
                omega.information.meat    = "default",
                omega.h1.information.meat = "default",

                bootstrap              = 1000L,
                gamma.n.minus.one      = FALSE,
                #gamma.unbiased         = FALSE,

                # optimization
                control                = list(),
                optim.method           = "nlminb",
                optim.attempts         = 4L,
                optim.force.converged  = FALSE,
                optim.gradient         = "analytic",
                optim.init_nelder_mead = FALSE,
                optim.var.transform    = "none",
                optim.parscale         = "none",
                optim.dx.tol           = 1e-03, # not too strict
                optim.bounds           = list(),
                em.iter.max            = 10000L,
                em.fx.tol              = 1e-08,
                em.dx.tol              = 1e-04,
                em.zerovar.offset      = 0.0001,
                optim.gn.iter.max      = 200L,
                optim.gn.tol.x         = 1e-07,

                # numerical integration
                integration.ngh        = 21L,

                # parallel
                parallel               = "no",
                ncpus                  = 1L,
                cl                     = NULL,
                iseed                  = NULL,

                # zero values
                zero.add               = "default",
                zero.keep.margins      = "default",
                zero.cell.warn         = FALSE, # since 0.6-1

                # starting values
                start                  = "default",

                # sanity checks
                check.start            = TRUE,
                check.post             = TRUE,
                check.gradient         = TRUE,
                check.vcov             = TRUE,
                check.lv.names         = TRUE,

                # more models/info
                h1                     = TRUE,
                baseline               = TRUE,
                baseline.conditional.x.free.slopes = TRUE,
                implied                = TRUE,
                loglik                 = TRUE,

                # storage of information
                store.vcov             = "default",

                # verbosity
                verbose                = FALSE,
                warn                   = TRUE,
                debug                  = FALSE)

    opt
}

# this function collects and checks the user-provided options/arguments,
# and fills in the "default" values, or changes them in an attempt to
# produce a consistent set of values...
#
# returns a list with the named options
lav_options_set <- function(opt = NULL) {

    if(opt$debug) { cat("lavaan DEBUG: lavaanOptions IN\n"); str(opt) }

    if(opt$debug) {
        opt$partrace <- TRUE
    } else {
        opt$partrace <- FALSE
    }

    # everything lowercase
    opt.old <- opt
    opt <- lapply(opt, function(x) { if(is.character(x)) tolower(x) else x})
    # except group,group.partial, which may contain capital letters
    opt$group.label <- opt.old$group.label
    opt$group.partial <- opt.old$group.partial


    # do.fit implies se="none and test="none" (unless not default)
    if(!opt$do.fit) {
        if(opt$se == "default") {
            opt$se <- "none"
        }
        if(length(opt$test) == 1L && opt$test == "default") {
            opt$test <- "none"
        }
    }

    # mimic
    if(opt$mimic == "default" || opt$mimic == "lavaan") {
        opt$mimic <- "lavaan"
    } else if(opt$mimic == "mplus") {
        opt$mimic <- "Mplus"
    } else if(opt$mimic == "eqs") {
        opt$mimic <- "EQS"
    } else if(opt$mimic == "lisrel") {
        cat("Warning: mimic=\"LISREL\" is not ready yet. Using EQS instead.\n")
        opt$mimic <- "EQS"
    } else {
        stop("lavaan ERROR: mimic must be \"lavaan\", \"Mplus\" or \"EQS\" \n")
    }

    # group.equal and group.partial
    if(opt$group.equal[1] == "none") {
        opt$group.equal <- character(0)
    } else if(is.null(opt$group.equal) || all(nchar(opt$group.equal) == 0L)) {
        if(opt$mimic == "Mplus") {
            if(opt$categorical) {
                opt$group.equal <- c("loadings", "thresholds")
            } else {
                if(is.logical(opt$meanstructure) && !opt$meanstructure) {
                    opt$group.equal <- "loadings"
                } else {
                    opt$group.equal <- c("loadings", "intercepts")
                }
            }
        } else {
            opt$group.equal <- character(0)
        }
    } else if(length(opt$group.equal) == 0) {
        # nothing to do
    } else if(all(opt$group.equal %in% c("loadings", "intercepts", "means",
                                         "composite.loadings",
                                         "regressions", "residuals",
                                         "residual.covariances", "thresholds",
                                         "lv.variances", "lv.covariances"))) {
        # nothing to do
    } else {
        wrong.idx <- which(!opt$group.equal %in%
            c("loadings", "intercepts", "means",
              "composite.loadings", "regressions", "residuals",
              "residual.covariances", "thresholds",
              "lv.variances", "lv.covariances"))
        stop("lavaan ERROR: unknown value for `group.equal' argument: ",
             sQuote(opt$group.equal[wrong.idx[1L]]), "\n")
    }
    if(is.null(opt$group.partial) || all(nchar(opt$group.partial) == 0L)) {
        opt$group.partial <- character(0)
    } else if(length(opt$group.partial) == 0) {
        # nothing to do
    } else {
        # strip white space
        opt$group.partial <- gsub("[[:space:]]+", "", opt$group.partial)
    }

    # if categorical, and group.equal contains "intercepts", also add
    # thresholds (and vice versa)
    if(opt$categorical && "intercepts" %in% opt$group.equal) {
        opt$group.equal <- unique(c(opt$group.equal, "thresholds"))
    }
    if(opt$categorical && "thresholds" %in% opt$group.equal) {
        opt$group.equal <- unique(c(opt$group.equal, "intercepts"))
    }

    # convenience (since 0.6-6)
    if(opt$se == "sandwich") {
        opt$se <- "robust.huber.white"
    }


    # representation
    if(opt$representation == "default") {
        opt$representation <- "LISREL"
    } else if(opt$representation == "lisrel") {
        opt$representation <- "LISREL"
    } else if(opt$representation == "eqs" ||
              opt$representation == "bentler-weeks") {
        opt$representation <- "EQS"
    } else {
        stop("lavaan ERROR: representation must be \"LISREL\" or \"EQS\" \n")
    }

    # clustered
    # brute-force override (for now)
    if(opt$clustered && !opt$multilevel) {
        opt$meanstructure <- TRUE
        #opt$missing <- "listwise"
        if(opt$missing == "ml") {
            optim.gradient = "numerical"
        }

        if(opt$estimator == "mlr") {
            opt$estimator <- "ml"
            opt$test <- "yuan.bentler.mplus"
            opt$se <- "robust.cluster"
        } else if(opt$estimator == "mlm") {
            opt$estimator <- "ml"
            opt$test <- "satorra.bentler"
            opt$se <- "robust.cluster.sem"
        }

        # test
        if(length(opt$test) == 1L && opt$test == "default") {
            opt$test <- "yuan.bentler.mplus"
        } else if(all(opt$test %in% c("none", "standard",
                                  "satorra.bentler",
                                  "yuan.bentler","yuan.bentler.mplus"))) {
            # nothing to do
        } else if(opt$se == "robust") {
            opt$test <- "yuan.bentler.mplus"
        } else {
            stop("lavaan ERROR: `test' argument must one of \"none\", \"yuan.bentler\", \"yuan.bentler.mplus\" or \"satorra.bentler\" in the clustered case")
        }

        # se
        if(opt$se == "default") {
            opt$se <- "robust.cluster"
        } else if(opt$se %in% c("none", "robust.cluster",
                                "robust.cluster.sem")) {
            # nothing to do
        } else if(opt$se == "robust") {
            opt$se <- "robust.cluster"
        } else {
            stop("lavaan ERROR: `se' argument must one of \"none\", \"robust.cluster\", or \"robust.cluster.sem\" in the clustered case")
        }

        # information
        if(opt$information[1] == "default") {
            if(opt$se == "robust.cluster") {
                opt$information[1] <- "observed"
            } else {
                opt$information[1] <- "expected"
            }
        }
        if(length(opt$information) > 1L && opt$information[2] == "default") {
            if(opt$se == "robust.cluster") {
                opt$information[2] <- "observed"
            } else {
                opt$information[2] <- "expected"
            }
        }
    }

    # multilevel
    # brute-force override (for now)
    if(opt$multilevel) {
        opt$meanstructure <- TRUE
        #opt$missing <- "listwise" # still needed for 0.6-8 (otherwise, we
                                  # we break tidySEM tests where they set
                                  # missing = "fiml" + multilevel

        # test
        if(length(opt$test) == 1L && opt$test == "default") {
            opt$test <- "standard"
        } else if(all(opt$test %in% c("none", "standard","yuan.bentler"))) {
            # nothing to do
        } else if(opt$se == "robust") {
            opt$test <- "yuan.bentler"
        } else {
            stop("lavaan ERROR: `test' argument must one of \"none\", \"standard\" or \"yuan.bentler\" in the multilevel case")
        }

        # se
        if(opt$se == "default") {
            opt$se <- "standard"
        } else if(opt$se %in% c("none", "standard", "robust.huber.white", "sandwich")) {
            # nothing to do
        } else if(opt$se == "robust.sem") {
            opt$se <- "robust.huber.white"
        } else {
            stop("lavaan ERROR: `se' argument must one of \"none\", \"standard\" or \"robust.huber.white\" in the multilevel case")
        }

        # information
        if(opt$information[1] == "default") {
            opt$information[1] <- "observed"
        }
        if(length(opt$information) > 1L && opt$information[2] == "default") {
            opt$information[2] <- "observed"
        }
    }


    # missing
    if(opt$missing == "default") {
        if(opt$mimic == "Mplus" && !opt$categorical &&
           opt$estimator %in% c("default", "ml", "mlr")) {
            # since version 5?
            opt$missing <- "ml"
            # check later if this is ok
        } else {
            opt$missing <- "listwise"
        }
    } else if(opt$missing %in% c("ml", "direct", "fiml")) {
        #if(opt$categorical && opt$estimator != "mml") {
        if(opt$categorical) {
            stop("lavaan ERROR: missing = ", dQuote(opt$missing),
                 " not available in the categorical setting")
        }
        opt$missing <- "ml"
        if(opt$estimator %in% c("mlm", "mlmv", "gls", "wls", "wlsm", "wlsmv",
                                "uls", "ulsm", "ulsmv", "pml")) {
            stop("lavaan ERROR: missing=\"ml\" is not allowed for estimator ",
                 dQuote(opt$estimator))
        }
    } else if(opt$missing %in% c("ml.x", "direct.x", "fiml.x")) {
        #if(opt$categorical && opt$estimator != "mml") {
        if(opt$categorical) {
            stop("lavaan ERROR: missing = ", dQuote(opt$missing),
                 " not available in the categorical setting")
        }
        opt$missing <- "ml.x"
        if(opt$estimator %in% c("mlm", "mlmv", "gls", "wls", "wlsm", "wlsmv",
                                "uls", "ulsm", "ulsmv", "pml")) {
            stop("lavaan ERROR: missing=\"ml\" is not allowed for estimator ",
                 dQuote(opt$estimator))
        }
    } else if(opt$missing %in% c("two.stage", "twostage", "two-stage",
                                 "two.step",  "twostep",  "two-step")) {
        opt$missing <- "two.stage"
        if(opt$categorical) {
            stop("lavaan ERROR: missing=\"two.stage\" not available in the categorical setting")
        }
        if(opt$estimator %in% c("mlm", "mlmv", "gls", "wls", "wlsm", "wlsmv",
                                "uls", "ulsm", "ulsmv", "pml", "mml")) {
            stop("lavaan ERROR: missing=\"two.stage\" is not allowed for estimator MLM, MLMV, GLS, ULS, ULSM, ULSMV, DWLS, WLS, WLSM, WLSMV, DLS, PML, MML")
        }
    } else if(opt$missing %in% c("robust.two.stage", "robust.twostage",
                                 "robust.two-stage", "robust-two-stage",
                                 "robust.two.step",  "robust.twostep",
                                 "robust-two-step")) {
        opt$missing <- "robust.two.stage"
        if(opt$categorical) {
            stop("lavaan ERROR: missing=\"robust.two.stage\" not available in the categorical setting")
        }
        if(opt$estimator %in% c("mlm", "mlmv", "gls", "wls", "wlsm", "wlsmv",
                                "uls", "ulsm", "ulsmv", "pml", "mml")) {
            stop("lavaan ERROR: missing=\"robust.two.stage\" is not allowed for estimator MLM, MLMV, GLS, ULS, ULSM, ULSMV, DWLS, WLS, WLSM, WLSMV, DLS, PML, MML")
        }
    } else if(opt$missing == "listwise") {
        # nothing to do
    } else if(opt$missing == "pairwise") {
        # nothing to do
    } else if(opt$missing == "available.cases") {
        # nothing to do, or warn if not categorical?
    } else if(opt$missing == "doubly.robust") {
        if(opt$estimator != "pml") {
            stop("lavaan ERROR: doubly.robust option only available for estimator PML")
        }
    } else if(opt$missing == "doubly_robust") {
        opt$missing <- "doubly.robust"
        if(opt$estimator != "pml") {
            stop("lavaan ERROR: doubly.robust option only available for estimator PML")
        }
    } else if(opt$missing == "available_cases") {
        opt$missing <- "available.cases"
    } else {
        stop("lavaan ERROR: unknown value for `missing' argument: ", opt$missing, "\n")
    }

    # default test statistic
    if(length(opt$test) == 1L && opt$test == "default") {
        if(opt$missing == "two.stage" ||
           opt$missing == "robust.two.stage") {
            opt$test <- "satorra.bentler"
        } else {
            opt$test <- "standard"
        }
    }

    # rename if needed
    if(length(target.idx <- which(opt$test %in%
        c("satorra", "sb", "satorra.bentler", "satorra-bentler",
          "m.adjusted", "m", "mean.adjusted", "mean-adjusted"))) > 0L) {
        opt$test[target.idx] <- "satorra.bentler"
    }
    if(length(target.idx <- which(opt$test %in%
        c("yuan", "yb", "yuan.bentler", "yuan-bentler"))) > 0L) {
        opt$test[target.idx] <- "yuan.bentler"
    }
    if(length(target.idx <- which(opt$test %in%
        c("yuan.bentler.mplus", "yuan-bentler.mplus",
          "yuan-bentler-mplus"))) > 0L) {
        opt$test[target.idx] <- "yuan.bentler.mplus"
    }
    if(length(target.idx <- which(opt$test %in%
        c("mean.var.adjusted", "mean-var-adjusted", "mv", "second.order",
          "satterthwaite", "mv.adjusted"))) > 0L) {
        opt$test[target.idx] <- "mean.var.adjusted"
    }
    if(length(target.idx <- which(opt$test %in%
        c("mplus6", "scale.shift", "scaled.shifted",
          "scaled-shifted"))) > 0L) {
        opt$test[target.idx] <- "scaled.shifted"
    }
    if(length(target.idx <- which(opt$test %in%
        c("bootstrap", "boot", "bollen.stine", "bollen-stine"))) > 0L) {
        opt$test[target.idx] <- "bollen.stine"
    }

    # check missing
    if(opt$missing %in% c("ml", "ml.x") && opt$se == "robust.sem") {
        warning("lavaan WARNING: missing will be set to ",
                    dQuote("listwise"), " for se = ",
                    dQuote(opt$se) )
        opt$missing <- "listwise"
    }
    if(opt$missing %in% c("ml", "ml.x") &&
       any(opt$test %in% c("satorra.bentler",
                       "mean.var.adjusted", "scaled.shifted"))) {
        warning("lavaan WARNING: missing will be set to ",
                    dQuote("listwise"), " for satorra.bentler style test")
        opt$missing <- "listwise"
    }

    # missing = "two.stage"
    if(opt$missing == "two.stage" ||
       opt$missing == "robust.two.stage") {
        opt$meanstructure <- TRUE
        # se
        if(opt$se == "default") {
            if(opt$missing == "two.stage") {
                opt$se <- "two.stage"
            } else {
                opt$se <- "robust.two.stage"
            }
        } else if(opt$missing == "two.stage" &&
                  opt$se      == "two.stage") {
            # nothing to do
        } else if(opt$missing == "robust.two.stage" &&
                  opt$se      == "robust.two.stage") {
            # nothing to do
        } else {
            warning("lavaan WARNING: se will be set to ",
                     dQuote(opt$missing), " if missing = ",
                     dQuote(opt$missing) )
            opt$se <- opt$missing
        }
        # information
        if(opt$information[1] == "default") {
            # for both two.stage and robust.two.stage
            opt$information[1] <- "observed"
        } else if(opt$information[1] == "first.order") {
            warning("lavaan WARNING: information will be set to ",
                     dQuote("observed"), " if missing = ",
                     dQuote(opt$missing) )
            opt$information[1] <- "observed"
        }
        # observed.information (ALWAYS "h1" for now)
        opt$observed.information[1] <- "h1"

        if(length(opt$information) > 1L && opt$information[2] == "default") {
            # for both two.stage and robust.two.stage
            opt$information[2] <- "observed"
        }
        # observed.information (ALWAYS "h1" for now)
        opt$observed.information[2] <- "h1"


        # test
        if(length(opt$test) > 1L) {
            warning("lavaan WARNING: test= argument can only contain a single element if missing = ", dQuote(opt$missing), " (taking the first)" )
            opt$test <- opt$test[1]
        }

        if(length(opt$test) == 1L && opt$test == "default") {
            opt$test <- "satorra.bentler"
        } else if(length(opt$test) == 1L && opt$test %in%
            c("satorra", "sb", "satorra.bentler", "satorra-bentler")) {
            opt$test <- "satorra.bentler"
        } else {
            warning("lavaan WARNING: test will be set to ",
                     dQuote("satorra.bentler"), " if missing = ",
                     dQuote(opt$missing) )
            opt$test <- "satorra.bentler"
        }
    }



    # meanstructure
    if(is.logical(opt$meanstructure)) {
        if(opt$meanstructure == FALSE) {
            # user explicitly wants meanstructure == FALSE
            # check for conflicting arguments
            #if(opt$estimator %in% c("mlm", "mlmv", "mlr", "mlf", "ulsm", "ulsmv", "wlsm", "wlsmv", "pml")) {
            #    warning("lavaan WARNING: estimator forces meanstructure = TRUE")
            #}
            if(opt$missing %in% c("ml", "ml.x", "two.stage")) {
                warning("lavaan WARNING: missing argument forces meanstructure = TRUE")
            }
        }
    } else if(opt$meanstructure == "default") {
        # by default: no meanstructure!
        if(opt$estimator == "pml") {
            opt$meanstructure <- TRUE
        } else if(opt$mimic == "Mplus") {
            opt$meanstructure <- TRUE
        } else {
            opt$meanstructure <- FALSE
        }
    } else {
        stop("lavaan ERROR: meanstructure must be TRUE, FALSE or \"default\"\n")
    }

    # estimator and se
    if(opt$se == "boot" || opt$se == "bootstrap") {
        opt$se <- "bootstrap"
        opt$information[1] <- "observed"
        if(length(opt$information) > 1L && opt$information[2] == "default") {
            opt$information[2] <- "observed"
        }
        opt$bootstrap <- as.integer(opt$bootstrap)
        stopifnot(opt$bootstrap > 0L)
    }

    # default estimator
    if(opt$estimator == "default") {
        if(opt$categorical) {
            opt$estimator <- "wlsmv"
        } else {
            opt$estimator <- "ml"
        }
    }

    # backwards compatibility (0.4 -> 0.5)
    if(opt$se == "robust.mlm") opt$se <- "robust.sem"
    if(opt$se == "robust.mlr") opt$se <- "robust.huber.white"

    if(opt$estimator == "ml") {
        opt$estimator <- "ML"
        if(opt$se == "default") {
            opt$se <- "standard"
        } else if(opt$se %in% c("bootstrap", "none",
                  "external", "standard", "robust.huber.white",
                  "robust.cluster", "robust.cluster.sem",
                  "two.stage", "robust.two.stage", "robust.sem")) {
            # nothing to do
        } else if(opt$se == "first.order") {
            # backwards compatibility
            opt$se <- "standard"
            opt$information[1] <- "first.order"
            if(length(opt$information) > 1L &&
               opt$information[2] == "default") {
                opt$information[2] <- "first.order"
            }
        } else if(opt$se == "observed") {
            opt$se <- "standard"
            opt$information[1] <- "observed"
            if(length(opt$information) > 1L &&
               opt$information[2] == "default") {
                opt$information[2] <- "observed"
            }
        } else if(opt$se == "expected") {
            opt$se <- "standard"
            opt$information[1] <- "expected"
            if(length(opt$information) > 1L &&
               opt$information[2] == "default") {
                opt$information[2] <- "expected"
            }
        } else if(opt$se == "robust") {
            if(opt$missing %in% c("ml", "ml.x")) {
                opt$se <- "robust.huber.white"
            } else if(opt$missing == "two.stage") {
                opt$se <- "robust.two.stage"
            } else {
                opt$se <- "robust.sem"
            }
        } else {
            stop("lavaan ERROR: unknown value for `se' argument when estimator is ML: ",
                 opt$se, "\n")
        }

    } else if(opt$estimator == "mlm"   ||
              opt$estimator == "mlmv"  ||
              opt$estimator == "mlmvs") {
        est.orig <- opt$estimator
        if(! (length(opt$test) == 1L && opt$test == "none") ) {
            if(opt$estimator == "mlm") {
                opt$test <- "satorra.bentler"
            } else if(opt$estimator == "mlmv") {
                opt$test <- "scaled.shifted"
            } else if(opt$estimator == "mlmvs") {
                opt$test <- "mean.var.adjusted"
            }
        }
        opt$estimator <- "ML"
        #opt$meanstructure <- TRUE
        if(opt$se == "bootstrap") {
            stop("lavaan ERROR: use ML estimator for bootstrap")
        }
        if(opt$se != "none" && opt$se != "external") {
            opt$se <- "robust.sem"
        }
        #opt$information[1] <- "expected"
        # in 0.6, we allow for information = "observed" as well
        opt$missing <- "listwise"
    } else if(opt$estimator == "mlf") {
        opt$estimator <- "ML"
        #opt$meanstructure <- TRUE
        if(opt$se == "bootstrap") {
            stop("lavaan ERROR: use ML estimator for bootstrap")
        }
        if(opt$se != "none" && opt$se != "external") {
            opt$se <- "standard"
            opt$information[1] <- "first.order"
            if(length(opt$information) > 1L &&
               opt$information[2] == "default") {
                opt$information[2] <- "first.order"
            }
        }
    } else if(opt$estimator == "mlr") {
        opt$estimator <- "ML"
        if(opt$se == "bootstrap") {
            stop("lavaan ERROR: use ML estimator for bootstrap")
        }
        if(opt$se != "none" && opt$se != "external") opt$se <- "robust.huber.white"
        if( !(length(opt$test) == 1L && opt$test == "none") &&
            opt$se != "external" ) {
            if(opt$mimic == "Mplus") {
                opt$test <- "yuan.bentler.mplus"
            } else if(opt$mimic == "EQS") {
                opt$test <- "yuan.bentler"
            } else {
                opt$test <- "yuan.bentler.mplus" # for now
            }
        }
        #opt$meanstructure <- TRUE
    } else if(opt$estimator == "gls") {
        opt$estimator <- "GLS"
        if(opt$se == "default" || opt$se == "standard") {
            opt$se <- "standard"
        } else if(opt$se == "none" ||
                  opt$se == "bootstrap" ||
                  opt$se == "external") {
            # nothing to do
        } else {
            stop("lavaan ERROR: invalid value for `se' argument when estimator is GLS: ",
                 opt$se, "\n")
        }
        if(!all(opt$test %in% c("standard","none"))) {
            stop("lavaan ERROR: invalid value for `test' argument when estimator is GLS: ",
                 opt$test, "\n")
        }
        opt$missing <- "listwise"
    } else if(opt$estimator == "ntrls") {
        opt$estimator <- "NTRLS"
        if(opt$se == "default" || opt$se == "standard") {
            opt$se <- "standard"
        } else if(opt$se == "none" ||
                  opt$se == "bootstrap" ||
                  opt$se == "external") {
            # nothing to do
        } else {
            stop("lavaan ERROR: invalid value for `se' argument when estimator is NTRLS: ",
                 opt$se, "\n")
        }
        if(!all(opt$test %in% c("standard","none"))) {
            stop("lavaan ERROR: invalid value for `test' argument when estimator is NTRLS: ",
                 opt$test, "\n")
        }
        opt$missing <- "listwise"
    } else if(opt$estimator == "wls") {
        opt$estimator <- "WLS"
        if(opt$se == "default" || opt$se == "standard") {
            opt$se <- "standard"
        } else if(opt$se == "none" ||
                  opt$se == "bootstrap" ||
                  opt$se == "external") {
            # nothing to do
        } else if(opt$se == "robust.sem") {
            # nothing to do
        } else if(opt$se == "robust") {
            opt$se <- "robust.sem"
        } else {
            stop("lavaan ERROR: invalid value for `se' argument when estimator is WLS: ",
                 opt$se, "\n")
        }
        if(!all(opt$test %in% c("standard","none"))) {
            stop("lavaan ERROR: invalid value for `test' argument when estimator is WLS: ",
                 opt$test, "\n")
        }
        #opt$missing <- "listwise"
    } else if(opt$estimator == "dls") {
        opt$optim.method <- "gn"
        opt$optim.gn.tol.x <- 1e-5
        opt$sample.cov.rescale <- TRUE
        opt$estimator <- "DLS"
        if(opt$se == "default") {
            opt$se <- "robust.sem"
        } else if(opt$se == "none" ||
                  opt$se == "bootstrap" ||
                  opt$se == "external") {
            # nothing to do
        } else if(opt$se == "robust.sem") {
            # nothing to do
        } else if(opt$se == "robust") {
            opt$se <- "robust.sem"
        } else {
            stop("lavaan ERROR: invalid value for `se' argument when estimator is DLS: ",
                 opt$se, "\n")
        }
        if(opt$test == "default" || opt$test == "standard") {
            opt$test <- "satorra.bentler"
        } else if(!all(opt$test %in% c("standard","none","satorra.bentler",
                            "mean.adjusted",
                            "mean.var.adjusted","scaled.shifted"))) {
            stop("lavaan ERROR: invalid value for `test' argument when estimator is DLS: ",
                 opt$test, "\n")
        }
        if(opt$missing %in% c("fiml", "ml", "direct")) {
            stop("lavaan ERROR: missing data is not supported if estimator is DLS; use missing = \"listwise\"")
        }
        opt$missing <- "listwise"

        # check estimator.args
        if(is.null(opt$estimator.args)) {
            opt$estimator.args <- list(dls.a = 1.0, dls.GammaNT = "sample")
        } else {
            if(is.null(opt$estimator.args$dls.a)) {
                opt$estimator.args$dls.a <- 1.0
            } else {
                stopifnot(is.numeric(opt$estimator.args$dls.a))
                if(opt$estimator.args$dls.a < 0.0 ||
                   opt$estimator.args$dls.a > 1.0) {
                    stop("lavaan ERROR: dls.a value in estimator.args must be between 0 and 1.")
                }
            }
            if(is.null(opt$estimator.args$dls.GammaNT)) {
                opt$estimator.args$dls.GammaNT <- "sample"
            } else {
                stopifnot(is.character(opt$estimator.args$dls.GammaNT))
                opt$estimator.args$dls.GammaNT <-
                    tolower(opt$estimator.args$dls.GammaNT)
                if(!opt$estimator.args$dls.GammaNT %in% c("sample", "model")) {
                    stop("lavaan ERROR: dls.GammaNT value in estimator.args must be either \"sample\" or \"model\".")
                }
            }
        }
    } else if(opt$estimator == "dwls") {
        opt$estimator <- "DWLS"
        if(opt$se == "default" || opt$se == "standard") {
            opt$se <- "standard"
        } else if(opt$se == "none" ||
                  opt$se == "bootstrap" ||
                  opt$se == "external") {
            # nothing to do
        } else if(opt$se == "robust.sem") {
            # nothing to do
        } else if(opt$se == "robust") {
            opt$se <- "robust.sem"
        } else {
            stop("lavaan ERROR: invalid value for `se' argument when estimator is DWLS: ",
                 opt$se, "\n")
        }
        if(!all(opt$test %in% c("standard","none","satorra.bentler",
                            "mean.adjusted",
                            "mean.var.adjusted","scaled.shifted"))) {
            stop("lavaan ERROR: invalid value for `test' argument when estimator is DWLS: ",
                 opt$test, "\n")
        }
        #opt$missing <- "listwise"
    } else if(opt$estimator == "wlsm") {
        opt$estimator <- "DWLS"
        if(opt$se == "bootstrap") {
            stop("lavaan ERROR: use (D)WLS estimator for bootstrap")
        }
        if(opt$se != "none" && opt$se != "external") opt$se <- "robust.sem"
        if(all(opt$test %in% c("mean.var.adjusted", "satorra.bentler",
                           "scaled.shifted"))) {
            # nothing to do
        } else if(! (length(opt$test) == 1L && opt$test == "none") ) {
            opt$test <- "satorra.bentler"
        }
        #opt$missing <- "listwise"
    } else if(opt$estimator == "wlsmv") {
        opt$estimator <- "DWLS"
        if(opt$se == "bootstrap") {
            stop("lavaan ERROR: use (D)WLS estimator for bootstrap")
        }
        if(opt$se != "none" && opt$se != "external") opt$se <- "robust.sem"
        if(all(opt$test %in% c("mean.var.adjusted", "satorra.bentler",
                           "scaled.shifted"))) {
            # nothing to do
        } else if(! (length(opt$test) == 1L && opt$test == "none") ) {
            opt$test <- "scaled.shifted"
        }
        #opt$missing <- "listwise"
    } else if(opt$estimator == "wlsmvs") {
        opt$estimator <- "DWLS"
        if(opt$se == "bootstrap") {
            stop("lavaan ERROR: use (D)WLS estimator for bootstrap")
        }
        if(opt$se != "none" && opt$se != "external") opt$se <- "robust.sem"
        if(all(opt$test %in% c("mean.var.adjusted", "satorra.bentler",
                           "scaled.shifted"))) {
            # nothing to do
        } else if(! (length(opt$test) == 1L && opt$test == "none") ) {
            opt$test <- "mean.var.adjusted"
        }
        #opt$missing <- "listwise"
    } else if(opt$estimator == "uls") {
        opt$estimator <- "ULS"
        if(opt$se == "default" || opt$se == "standard") {
            opt$se <- "standard"
        } else if(opt$se == "none" ||
                  opt$se == "bootstrap" ||
                  opt$se == "external") {
            # nothing to do
        } else if(opt$se == "robust.sem") {
            # nothing to do
        } else if(opt$se == "robust") {
            opt$se <- "robust.sem"
        } else {
            stop("lavaan ERROR: invalid value for `se' argument when estimator is ULS: ",
                 opt$se, "\n")
        }
        if(!all(opt$test %in% c("standard","none", "satorra.bentler",
                            "mean.adjusted",
                            "mean.var.adjusted","scaled.shifted"))) {
            stop("lavaan ERROR: invalid value for `test' argument when estimator is ULS: ",
                 opt$test, "\n")
        }
        #opt$missing <- "listwise"
    } else if(opt$estimator == "ulsm") {
        opt$estimator <- "ULS"
        if(opt$se == "bootstrap") {
            stop("lavaan ERROR: use ULS estimator for bootstrap")
        }
        if(opt$se != "none" && opt$se != "external") opt$se <- "robust.sem"
        if(all(opt$test %in% c("mean.var.adjusted", "satorra.bentler",
                           "scaled.shifted"))) {
            # nothing to do
        } else if(! (length(opt$test) == 1L && opt$test == "none") ) {
            opt$test <- "satorra.bentler"
        }
        #opt$missing <- "listwise"
    } else if(opt$estimator == "ulsmv") {
        opt$estimator <- "ULS"
        if(opt$se == "bootstrap") {
            stop("lavaan ERROR: use ULS estimator for bootstrap")
        }
        if(opt$se != "none" && opt$se != "external") opt$se <- "robust.sem"
        if(all(opt$test %in% c("mean.var.adjusted", "satorra.bentler",
                           "scaled.shifted"))) {
            # nothing to do
        } else if(! (length(opt$test) == 1L && opt$test == "none") ) {
            opt$test <- "scaled.shifted"
        }
        #opt$missing <- "listwise"
    } else if(opt$estimator == "ulsmvs") {
        opt$estimator <- "ULS"
        if(opt$se == "bootstrap") {
            stop("lavaan ERROR: use ULS estimator for bootstrap")
        }
        if(opt$se != "none" && opt$se != "external") opt$se <- "robust.sem"
        if(all(opt$test %in% c("mean.var.adjusted", "satorra.bentler",
                           "scaled.shifted"))) {
            # nothing to do
        } else if(! (length(opt$test) == 1L && opt$test == "none") ) {
            opt$test <- "mean.var.adjusted"
        }
        #opt$missing <- "listwise"
    } else if(opt$estimator == "pml") {
        opt$estimator <- "PML"
        opt$information[1] <- "observed"
        if(length(opt$information) > 1L &&
               opt$information[2] == "default") {
            opt$information[2] <- "observed"
        }
        if(opt$se == "default") {
            opt$se <- "robust.huber.white"
        }
        if(! (length(opt$test) == 1L && opt$test == "none") ) {
            opt$test <- "mean.var.adjusted"
        }
        #opt$missing <- "listwise"
    } else if(opt$estimator %in% c("fml","umn")) {
        opt$estimator <- "FML"
        opt$information[1] <- "observed"
        if(length(opt$information) > 1L &&
               opt$information[2] == "default") {
            opt$information[2] <- "observed"
        }
        if(opt$se == "default") {
            opt$se <- "standard"
        }
        if(! (length(opt$test) == 1L && opt$test == "none") ) {
            opt$test <- "standard"
        }
        #opt$missing <- "listwise"
    } else if(opt$estimator == "reml") {
        opt$estimator <- "REML"
        opt$information[1] <- "observed"
        if(length(opt$information) > 1L &&
               opt$information[2] == "default") {
            opt$information[2] <- "observed"
        }
        if(opt$se == "default") {
            opt$se <- "standard"
        }
        if(! (length(opt$test) == 1L && opt$test == "none") ) {
            opt$test <- "standard"
        }
        opt$missing <- "listwise"
    } else if(opt$estimator %in% c("mml")) {
        opt$estimator <- "MML"
        opt$information[1] <- "observed"
        opt$meanstructure <- TRUE
        if(length(opt$information) > 1L &&
               opt$information[2] == "default") {
            opt$information[2] <- "observed"
        }
        if(opt$se == "default") {
            opt$se <- "standard"
        }
        if(length(opt$test) == 1L && opt$test == "default") {
            opt$test <- "none"
        }
        #opt$missing <- "listwise"
        if(opt$link == "default") {
            #opt$link <- "logit"
            opt$link <- "probit"
        } else if(opt$link %in% c("logit","probit")) {
            # nothing to do
        } else {
            stop("lavaan ERROR: link must be `logit' or `probit'")
        }
        # check for parameterization
        if(opt$parameterization == "default") {
            opt$parameterization <- "mml"
        } else {
            stop("lavaan WARNING: parameterization argument is ignored if estimator = MML")
        }
    } else if(opt$estimator == "none") {
        if(opt$se == "default") {
            opt$se <- "none"
        }
        if(length(opt$test) == 1L && opt$test == "default") {
            opt$test <- "none"
        }
    } else {
        stop("lavaan ERROR: unknown value for `estimator' argument: ", opt$estimator, "\n")
    }


    # special stuff for categorical
    if(opt$categorical) {
        opt$meanstructure <- TRUE # Mplus style
        if(opt$estimator == "ML") {
            stop("lavaan ERROR: estimator ML for ordered data is not supported yet. Use WLSMV instead.")
        }
    }

    # link
    if(opt$link == "logit") {
        if(opt$estimator != "mml") {
             warning("lavaan WARNING: link will be set to ",
                    dQuote("probit"), " for estimator = ",
                    dQuote(opt$estimator) )
        }
    }

    # likelihood approach (wishart or normal) + sample.cov.rescale
    if(!opt$estimator %in% c("ML", "REML", "PML", "FML","NTRLS")) {
        if(opt$likelihood != "default") {
            stop("lavaan ERROR: likelihood argument is only relevant if estimator = ML")
        }
        if(opt$sample.cov.rescale == "default") {
            opt$sample.cov.rescale <- FALSE
        }# else {
        #    warning("sample.cov.rescale argument is only relevant if estimator = ML")
        #}
    } else { # ml and friends
        if(opt$estimator %in% c("PML", "FML")) {
            opt$likelihood <- "normal"
        } else if(opt$likelihood == "default") {
           opt$likelihood <- "normal"
            if(opt$mimic == "EQS"    ||
               opt$mimic == "LISREL" ||
               opt$mimic == "AMOS") {
                opt$likelihood <- "wishart"
            }
        } else if(opt$likelihood == "wishart" || opt$likelihood == "normal") {
            # nothing to do
        } else {
            stop("lavaan ERROR: invalid value for `likelihood' argument: ",
                 opt$likelihood, "\n")
        }

        if(opt$sample.cov.rescale == "default") {
            opt$sample.cov.rescale <- FALSE
            if(opt$likelihood == "normal") {
                opt$sample.cov.rescale <- TRUE
            }
        } else if(!is.logical(opt$sample.cov.rescale)) {
            stop("lavaan ERROR: sample.cov.rescale must be either \"default\", TRUE, or FALSE")
        } else {
            # nothing to do
        }
    }

    # se information
    if(opt$information[1] == "default") {
        if(opt$missing %in% c("ml", "ml.x") ||
           opt$se == "robust.huber.white"   ||
           opt$se == "first.order") {
           #nchar(opt$constraints) > 0L) {
            opt$information[1] <- "observed"
        } else {
            opt$information[1] <- "expected"
        }
    } else if(opt$information[1] %in%
              c("observed", "expected", "first.order")) {
        # nothing to do
    } else {
        stop("lavaan ERROR: information must be either \"expected\", \"observed\", or \"first.order\"\n")
    }

    # first.order information can not be used with robust
    if(opt$information[1] == "first.order" &&
       opt$se %in% c("robust.huber.white", "robust.sem")) {
        stop("lavaan ERROR: information must be either \"expected\" or \"observed\" if robust standard errors are requested.")
    }

    # test information
    if(length(opt$information) == 1L) {
        opt$information <- rep(opt$information, 2L)
    }
    if(opt$information[2] == "default") {
        if(opt$missing %in% c("ml", "ml.x") ||
           opt$se == "robust.huber.white"   ||
           opt$se == "first.order") {
           #nchar(opt$constraints) > 0L) {
            opt$information[2] <- "observed"
        } else {
            opt$information[2] <- "expected"
        }
    } else if(opt$information[2] %in%
              c("observed", "expected", "first.order")) {
        # nothing to do
    } else {
        stop("lavaan ERROR: information[2] must be either \"expected\", \"observed\", or \"first.order\"\n")
    }

    # first.order information can not be used with robust
    if(opt$information[2] == "first.order" &&
       opt$test %in% c("satorra.bentler",
                       "yuan.bentler",
                       "yuan.bentler.mplus",
                       "mean.var.adjusted",
                       "scaled.shifted")) {
        stop("lavaan ERROR: information must be either \"expected\" or \"observed\" if robust test statistics are requested.")
    }

    # information meat
    if(length(opt$information.meat) > 1L) {
        warning("lavaan WARNING: only first element of information.meat is used")
        opt$information.meat <- opt$information.meat[1]
    }
    if(opt$information.meat == "default") {
        opt$information.meat <- "first.order"
    } else if(opt$information.meat %in% c("first.order")) {
        # nothing to do
    } else {
        stop("lavaan ERROR: information.meat must be \"first.order\" (for now) \n")
    }

    if(opt$observed.information[1] == "hessian" ||
       opt$observed.information[1] == "h1") {
        # nothing to do
    } else {
        stop("lavaan ERROR: observed.information must be either \"hessian\", or \"h1\"\n")
    }

    if(length(opt$observed.information) == 1L) {
        opt$observed.information <- rep(opt$observed.information, 2L)
    }
    if(opt$observed.information[2] == "hessian" ||
       opt$observed.information[2] == "h1") {
        # do nothing
    } else if(opt$observed.information[2] == "default") {
        if(opt$test %in% c("satorra.bentler",
                           "yuan.bentler",
                           "yuan.bentler.mplus",
                           "mean.var.adjusted",
                           "scaled.shifted")) {
            if(opt$estimator == "PML" || opt$test == "yuan.bentler.mplus") {
                opt$observed.information[2] <- "hessian"
            } else {
                opt$observed.information[2] <- "h1" # CHANGED in 0.6-6!
            }
        } else {
            # default is "hessian"
            opt$observed.information[2] <- "hessian"
        }
    } else {
        stop("lavaan ERROR: observed.information[2] must be either \"hessian\", or \"h1\"\n")
    }

    if(opt$h1.information[1] == "structured" ||
       opt$h1.information[1] == "unstructured") {
        # nothing to do
    } else {
        stop("lavaan ERROR: h1.information must be either \"structured\" or \"unstructured\"\n")
    }

    if(length(opt$h1.information) == 1L) {
        opt$h1.information <- rep(opt$h1.information, 2L)
    }
    if(opt$h1.information[2] == "structured" ||
       opt$h1.information[2] == "unstructured") {
        # nothing to do
    } else {
        stop("lavaan ERROR: h1.information[2] must be either \"structured\" or \"unstructured\"\n")
    }

    if(length(opt$h1.information.meat) > 1L) {
        warning("lavaan WARNING: only first element of h1.information.meat is used")
        opt$h1.information.meat <- opt$h1.information.meat[1]
    }
    if(opt$h1.information.meat == "default") {
        opt$h1.information.meat <- opt$h1.information[1]
    } else if(opt$h1.information.meat == "structured" ||
              opt$h1.information.meat == "unstructured") {
        # nothing to do
    } else {
        stop("lavaan ERROR: h1.information.meat must be either \"structured\" or \"unstructured\"\n")
    }

    # check information if estimator is uls/wls and friends
    if(opt$estimator %in% c("ULS", "WLS", "DWLS")) {
        if(opt$information[1] != "expected") {
            warning("lavaan WARNING: information will be set to ",
                    dQuote("expected"), " for estimator = ", dQuote(opt$estimator))
            opt$information[1] <- "expected"
            opt$information[2] <- "expected"
        }
        opt$h1.information[1] <- "unstructured" # FIXME: allow option?
        opt$h1.information[2] <- "unstructured" # FIXME: allow option?
    }


    # omega information
    if(opt$omega.information == "default") {
        opt$omega.information <- opt$information[2] # test version!
    } else if(opt$omega.information %in% c("expected", "observed")) {
        # nothing to do
    } else {
        stop("lavaan ERROR: omega.information must be either \"expected\" or \"observed\"")
    }

    if(opt$omega.h1.information == "default") {
        #opt$omega.h1.information <- opt$h1.information[2] # test version!
        opt$omega.h1.information <- "unstructured"
    } else if(opt$omega.h1.information %in% c("structured", "unstructured")) {
        # nothing to do
    } else {
        stop("lavaan ERROR: omega.h1.information must be either \"structured\" or \"unstructured\"")
    }

    # omega information.meat
    if(opt$omega.information.meat == "default") {
        opt$omega.information.meat <- "first.order"
    } else if(opt$omega.information %in% c("first.order")) {
        # nothing to do
    } else {
        stop("lavaan ERROR: omega.information.meat must be \"first.order\"")
    }

    if(opt$omega.h1.information.meat == "default") {
        opt$omega.h1.information.meat <- opt$omega.h1.information
    } else if(opt$omega.h1.information.meat %in%
              c("structured", "unstructured")) {
        # nothing to do
    } else {
        stop("lavaan ERROR: omega.h1.information.meat must be either \"structured\" or \"unstructured\"")
    }




    # conditional.x
    if(is.logical(opt$conditional.x)) {
    } else if(opt$conditional.x == "default") {
        if(opt$estimator == "ML" && (opt$mimic == "Mplus" ||
                                     opt$mimic == "lavaan")) {
            opt$conditional.x <- FALSE
        } else if(opt$categorical) {
            opt$conditional.x <- TRUE
        } else {
            opt$conditional.x <- FALSE
        }
    } else {
        stop("lavaan ERROR: conditional.x must be TRUE, FALSE or \"default\"\n")
    }

    # if conditional.x, always use a meanstructure
    if(opt$conditional.x) {
        opt$meanstructure <- TRUE
    }

    # fixed.x
    if(is.logical(opt$fixed.x)) {
        if(opt$conditional.x && opt$fixed.x == FALSE) {
            stop("lavaan ERROR: fixed.x = FALSE is not supported when conditional.x = TRUE.")
        }
    } else if(opt$fixed.x == "default") {
        if(opt$estimator %in% c("MML", "ML") && (opt$mimic == "Mplus" ||
                                     opt$mimic == "lavaan")) {
            opt$fixed.x <- TRUE
        } else if(opt$conditional.x) {
            opt$fixed.x <- TRUE
        } else {
            opt$fixed.x <- FALSE
        }
    } else {
        stop("lavaan ERROR: fixed.x must be TRUE, FALSE or \"default\"\n")
    }


    # meanstructure again
    if(opt$missing %in% c("ml", "ml.x") || opt$model.type == "growth") {
        opt$meanstructure <- TRUE
    }
    if("intercepts" %in% opt$group.equal ||
       "means" %in% opt$group.equal) {
        opt$meanstructure <- TRUE
    }
    #if(opt$se == "robust.huber.white" ||
    #   opt$se == "robust.sem" ||
    #   opt$test == "satorra.bentler" ||
    #   opt$test == "mean.var.adjusted" ||
    #   opt$test == "scaled.shifted" ||
    #   opt$test == "yuan.bentler") {
    #    opt$meanstructure <- TRUE
    #}
    stopifnot(is.logical(opt$meanstructure))
    stopifnot(is.logical(opt$verbose))
    stopifnot(is.logical(opt$warn))

    if(opt$debug) {
        opt$verbose <- opt$warn <- TRUE
    }

    # zero cell frequencies
    if(is.character(opt$zero.add) && opt$zero.add == "default") {
        # default: c(0.5, 0.0)
        opt$zero.add <- c(0.5, 0.0)
        # FIXME: TODO: mimic EQS , LISREL (0.0, 0.0)
    } else if(is.numeric(opt$zero.add)) {
        if(length(opt$zero.add) == 1L) {
            opt$zero.add <- c(opt$zero.add, opt$zero.add)
        } else if(length(opt$zero.add) > 2L) {
            warning("lavaan WARNING: argument `zero.add' only uses the first two numbers")
            opt$zero.add <- opt$zero.add[1:2]
        }
    } else {
       stop("lavaan ERROR: argument `zero.add' must be numeric or \"default\"")
    }

    if(is.character(opt$zero.keep.margins) &&
       opt$zero.keep.margins == "default") {
        if(opt$mimic %in% c("lavaan", "Mplus")) {
            opt$zero.keep.margins <- TRUE
        } else {
            opt$zero.keep.margins <- FALSE
        }
    } else if(is.logical(opt$zero.keep.margins)) {
        # nothing to do
    } else {
        stop("lavaan ERROR: argument `zero.keep.margins' must be logical or \"default\"")
    }

    # parameterization
    if(opt$parameterization == "default") {
        # for now, default is always delta
        opt$parameterization <- "delta"
    } else if(opt$parameterization %in% c("delta", "theta", "mml")) {
        # nothing to do
    } else {
        stop("lavaan ERROR: argument `parameterization' should be `delta' or `theta'")
    }



    # std.lv vs auto.fix.first # new in 0.6-5 (used to be in sem/cfa/growth)
    if(opt$std.lv) {
        opt$auto.fix.first <- FALSE
    }

    # std.lv vs effect.coding # new in 0.6-4
    if(is.logical(opt$effect.coding)) {
        if(opt$effect.coding) {
            opt$effect.coding <- c("loadings", "intercepts")
        } else {
            opt$effect.coding <- ""
        }
    } else if(length(opt$effect.coding) == 0L) {
        # nothing to do
    } else if(all(opt$effect.coding %in% c("loadings", "intercepts",
                                           "mg.lv.efa.variances",
                                           "mg.lv.variances",
                                           "mg.lv.means",
                                           "mg.lv.intercepts"))) {
        # nothing to do
    } else {
        stop("lavaan ERROR: unknown value for ", sQuote("effect.coding"),
             " argument: ", opt$effect.coding, "\n")
    }

    # if we use effect coding for the factor loadings, we don't need/want
    # std.lv = TRUE
    if("loadings" %in% opt$effect.coding) {
        if(opt$std.lv) {
            stop("lavaan ERROR: std.lv is set to FALSE but effect.coding contains ", dQuote("loadings"))
            opt$std.lv <- FALSE
        }
        # shut off auto.fix.first
        opt$auto.fix.first <- FALSE
    }

    # test again

    # unless test = "none", always add test = "standard" as the
    # first entry
    # NO: this breaks lavaan.survey pval.pFsum, which has the following check:
    # if (!lavInspect(lavaan.fit, "options")$test %in% c("satorra.bentler",
    #    "mean.var.adjusted", "Satterthwaite")) {
    #    stop("Please refit the model with Satorra-Bentler (MLM) or Satterthwaite (MLMVS) adjustment.")
    #}
    #if(! (length(opt$test) == 1L && opt$test == "none") ) {
    #    opt$test <- c("standard", opt$test)
    #    opt$test <- unique(opt$test)
    #}

    # final check
    wrong.idx <- which(! opt$test %in% c("none", "standard", "satorra.bentler",
                            "yuan.bentler", "yuan.bentler.mplus",
                            "mean.var.adjusted", "scaled.shifted",
                            "bollen.stine"))
    if(length(wrong.idx) > 0L) {
        stop("lavaan ERROR: invalid option(s) for test argument: ",
             paste(dQuote(opt$test[wrong.idx]), collapse = " "),
             "\n\t\t",
             "Possible options are: \"none\", \"standard\",
             \"satorra.bentler\", \"yuan.bentler\", \"yuan.bentler.mplus\",
             \"mean.var.adjusted\", \"scaled.shifted\", or \"bollen.stine\"")
    }


    # optim.bounds
    if(!is.null(opt$optim.bounds) && length(opt$optim.bounds) > 0L) {
        # opt$bounds should be "default"
        if(is.null(opt$bounds) || opt$bounds == "none") {
            opt$bounds <- "user"
        } else {
            stop("lavaan ERROR: bounds and optim.bounds arguments can not be used together")
        }
    }

    # bounds
    if(is.null(opt$bounds)) {
        opt$bounds <- "none" # for now
    } else if(is.logical(opt$bounds)) {
        if(opt$bounds) {
            opt$bounds <- "default"
        } else {
            opt$bounds <- "none"
        }
    }

    # handle different 'profiles'
    if(opt$bounds == "none") {
        opt$optim.bounds <- list(lower = character(0L),
                                 upper = character(0L))
    } else if(opt$bounds == "user") {
        # nothing to do
    } else if(opt$bounds == "default") {
        opt$optim.bounds <- list(lower = c("ov.var", "lv.var", "loadings"),
                                 upper = c("ov.var", "lv.var", "loadings"),
                                 lower.factor = c(1.2, 1.0, 1.1),
                                 upper.factor = c(1.2, 1.3, 1.1),
                                 min.reliability.marker = 0.1,
                                 min.var.lv.exo = 0.0,
                                 min.var.lv.endo = 0.0)
    } else if(opt$bounds == "pos.var") {
        opt$optim.bounds <- list(lower = c("ov.var", "lv.var"),
                                 lower.factor = c(1, 1),
                                 min.reliability.marker = 0.0,
                                 min.var.lv.exo = 0.0,
                                 min.var.lv.endo = 0.0)
    } else if(opt$bounds == "pos.ov.var") {
        opt$optim.bounds <- list(lower = c("ov.var"),
                                 lower.factor = 1)
    } else if(opt$bounds == "pos.lv.var") {
        opt$optim.bounds <- list(lower = c("lv.var"),
                                 lower.factor = 1,
                                 min.reliability.marker = 0.0,
                                 min.var.lv.exo = 0.0,
                                 min.var.lv.endo = 0.0)
    } else {
        stop("lavaan ERROR: unknown `bounds' option: ", opt$bounds)
    }


    # rotation
    opt$rotation <- tolower(opt$rotation)
    if(opt$rotation %in% c("crawfer", "crawford.ferguson", "crawford-ferguson",
                           "crawfordferguson")) {
        opt$rotation <- "cf"
    }
    if(opt$rotation %in% c("varimax", "quartimax", "orthomax", "cf", "oblimin",
                     "quartimin", "geomin", "entropy", "mccammon", "infomax",
                     "tandem1", "tandem2", "none",
                     "oblimax", "bentler", "simplimax", "target", "pst")) {
        # nothing to do
    } else if(opt$rotation %in% c("cf-quartimax", "cf-varimax", "cf-equamax",
                            "cf-parsimax", "cf-facparsim")) {
        # nothing to do here; we need M/P to set cf.gamma
    } else {
        txt <- c("Rotation method ", dQuote(opt$rotation), " not supported. ",
        "Supported rotation methods are: varimax, quartimax, orthomax, cf, ",
        "oblimin, quartimin, geomin, entropy, mccammon, infomax,",
        "tandem1, tandem2, oblimax, bentler, simplimax, target, pst, ",
        "crawford-ferguson,  cf-quartimax,  cf-varimax, cf-equamax, ",
        "cf-parsimax, cf-facparsim")
        stop(lav_txt2message(txt, header = "lavaan ERROR:"))
    }

    # rotation.se
    if(!opt$rotation.se %in% c("delta", "bordered")) {
        stop("lavaan ERROR: rotation.se option must be either \"delta\" or \"bordered\".")
    }

    # rotations.args
    if(!is.list(opt$rotation.args)) {
        stop("lavaan ERROR: rotation.args should be be list.")
    }

    # if target, check target matrix
    if(opt$rotation == "target" || opt$rotation == "pst") {
        target <- opt$rotation.args$target
        if(is.null(target) || !is.matrix(target)) {
            stop("lavaan ERROR: ",
                 "rotation target matrix is NULL, or not a matrix")
        }
    }
    if(opt$rotation == "pst") {
        target.mask <- opt$rotation.args$target.mask
        if(is.null(target.mask) || !is.matrix(target.mask)) {
            stop("lavaan ERROR: ",
                 "rotation target.mask matrix is NULL, or not a matrix")
        }
    }
    # if NAs, force opt$rotation to be 'pst' and create target.mask
    if(opt$rotation == "target" && anyNA(target)) {
        opt$rotation <- "pst"
        target.mask <- matrix(1, nrow = nrow(target), ncol = ncol(target))
        target.mask[ is.na(target) ] <- 0
        opt$rotation.args$target.mask <- target.mask
    }

    # set row.weights
    opt$rotation.args$row.weights <- tolower(opt$rotation.args$row.weights)
    if(opt$rotation.args$row.weights == "default") {
        # the default is "none", except for varimax
        if(opt$rotation == "varimax") {
            opt$rotation.args$row.weights <- "kaiser"
        } else {
            opt$rotation.args$row.weights <- "none"
        }
    } else if(opt$rotation.args$row.weights %in% c("cureton-mulaik",
              "cureton.mulaik", "cm")) {
    } else if(opt$rotation.args$row.weights %in% c("kaiser", "none")) {
        # nothing to do
    } else {
        stop("lavaan ERROR: rotation.args$row.weights should be \"none\",",
             " \"kaiser\" or \"cureton-mulaik\".")
    }

    # check opt$rotation.args$algorithm
    opt$rotation.args$algorithm <- tolower(opt$rotation.args$algorithm)
    if(opt$rotation.args$algorithm %in% c("gpa", "pairwise")) {
        # nothing to do
    } else {
        stop("lavaan ERROR: opt$rotation.args$algorithm must be gpa or pairwise")
    }

    # order.lv.by
    opt$rotation.args$order.lv.by <- tolower(opt$rotation.args$order.lv.by)
    if(opt$rotation.args$order.lv.by %in% c("sumofsquares", "index", "none")) {
        # nothing to do
    } else {
        stop("lavaan ERROR: rotation.args$order.lv.by should be \"none\",",
             " \"index\" or \"sumofsquares\".")
    }



    # group.w.free
    #if(opt$group.w.free && opt$categorical) {
    #    stop("lavaan ERROR: group.w.free = TRUE is not supported (yet) in the categorical setting.")
    #}

    if(opt$debug) { cat("lavaan DEBUG: lavaanOptions OUT\n"); str(opt) }

    opt
}
