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
                correlation        = FALSE,     # correlation structure
                effect.coding      = FALSE,     # TRUE implies
                                                # c("loadings", "intercepts")
                ceq.simple         = FALSE,      # treat simple eq cons special?
                parameterization   = "default",
                #ov.order           = "model",  # how to 'order' ov's in the
                #                               # partable

                auto.fix.first     = FALSE,
                auto.fix.single    = FALSE,
                auto.var           = FALSE,
                auto.cov.lv.x      = FALSE,
                auto.cov.y         = FALSE,
                auto.th            = FALSE,
                auto.delta         = FALSE,
                auto.efa           = FALSE,

                # rotation
                rotation           = "geomin",
                rotation.se        = "bordered", # "bordered" or "delta"
                rotation.args      = list(orthogonal     = FALSE,
                                          row.weights    = "default",
                                          std.ov         = TRUE,
                                          geomin.epsilon = 0.001, # was 0.01 < 0.6-10
                                          orthomax.gamma = 1,
                                          cf.gamma       = 0,
                                          oblimin.gamma  = 0,
                                          promax.kappa   = 4,
                                          target         = matrix(0,0,0),
                                          target.mask    = matrix(0,0,0),
                                          rstarts        = 30L,
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
                sample.icov        = TRUE,
                ridge              = FALSE,
                ridge.constant     = "default",

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

                # default test statistic for scaling
                scaled.test            = "standard",

                # old approach trace.UGamma2
                ug2.old.approach       = FALSE,

                # bootstrap
                bootstrap              = 1000L,

                # gamma
                gamma.n.minus.one      = FALSE,
                gamma.unbiased         = FALSE,

                # optimization
                control                = list(),
                optim.method           = "default", # gn for DLS, nlminb rest
                optim.attempts         = 4L,
                optim.force.converged  = FALSE,
                optim.gradient         = "analytic",
                optim.init_nelder_mead = FALSE,
                optim.var.transform    = "none",
                optim.parscale         = "none",
                optim.partrace         = FALSE,
                optim.dx.tol           = 1e-03, # not too strict
                optim.bounds           = list(),
                em.iter.max            = 10000L,
                em.fx.tol              = 1e-08,
                em.dx.tol              = 1e-04,
                em.zerovar.offset      = 0.0001,
                em.h1.iter.max         = 500L,
                em.h1.tol              = 1e-05, # was 1e-06 < 0.6-9
                em.h1.warn             = TRUE,
                optim.gn.iter.max      = 200L,
                optim.gn.stephalf.max  = 10L,
                optim.gn.tol.x         = 1e-05,

                # numerical integration
                integration.ngh        = 21L,

                # parallel
                parallel               = "no",
                ncpus                  = parallel::detectCores() - 1L,
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
        opt$optim.partrace <- TRUE
    }

    # everything lowercase
    opt.old <- opt
    opt <- lapply(opt, function(x) { if(is.character(x)) tolower(x) else x})
    # except group,group.partial, which may contain capital letters
    opt$group.label <- opt.old$group.label
    opt$group.partial <- opt.old$group.partial
    rm(opt.old)

    # first of all: set estimator
    if(opt$estimator == "default") {
        if(opt$.categorical) {
            opt$estimator <- "wlsmv"
        } else {
            opt$estimator <- "ml"
        }
    }
    # store lower-case estimator name
    orig.estimator <- opt$estimator

    # rename names of test statistics if needed, and check for invalid values
    opt$test <- lav_test_rename(opt$test, check = TRUE)

    # same for scaled.test
    opt$scaled.test <- lav_test_rename(opt$scaled.test, check = TRUE)

    # rename names of se values, and check for invalid values
    # pass-through function: may change value of information
    # for backwards compatibility (eg if se = "expected")
    opt <- lav_options_check_se(opt)

    # do.fit implies se="none and test="none" (unless not default)
    if(!opt$do.fit) {
        if(opt$se == "default") {
            opt$se <- "none"
        }
        if(opt$test[1] == "default") {
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
    } else if(opt$mimic %in% c("lm", "LM", "regression")) {
        opt$mimic <- "lm"
    } else {
        stop("lavaan ERROR: mimic must be \"lavaan\", \"Mplus\" or \"EQS\" \n")
    }

    # group.equal and group.partial
    if(opt$group.equal[1] == "none") {
        opt$group.equal <- character(0)
    } else if(is.null(opt$group.equal) || all(nchar(opt$group.equal) == 0L)) {
        if(opt$mimic == "Mplus") {
            if(opt$.categorical) {
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
    if(opt$.categorical && "intercepts" %in% opt$group.equal) {
        opt$group.equal <- unique(c(opt$group.equal, "thresholds"))
    }
    if(opt$.categorical && "thresholds" %in% opt$group.equal) {
        opt$group.equal <- unique(c(opt$group.equal, "intercepts"))
    }

    # representation
    if(opt$representation == "default") {
        opt$representation <- "LISREL"
    } else if(opt$representation %in% c("lisrel", "LISREL")) {
        opt$representation <- "LISREL"
    #} else if(opt$representation %in% c("eqs", "EQS", "bentler-weeks")) {
    #    opt$representation <- "EQS"
    } else if(opt$representation %in% c("ram", "RAM")) {
        opt$representation <- "RAM"
    } else {
        stop("lavaan ERROR: representation must be \"LISREL\" or \"RAM\" \n")
    }

    # clustered
    # brute-force override (for now)
    if(opt$.clustered && !opt$.multilevel) {
        opt$meanstructure <- TRUE

        if(opt$estimator == "mlr") {
            opt$estimator <- "ml"
            opt$test <- "yuan.bentler.mplus"
            opt$se <- "robust.cluster"
        } else if(opt$estimator == "mlm") {
            opt$estimator <- "ml"
            opt$test <- "satorra.bentler"
            opt$se <- "robust.cluster.sem"
        } else if(opt$.categorical) {
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
            if(opt$se == "robust.cluster" && opt$estimator == "ml") {
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
    if(opt$.multilevel) {
        opt$meanstructure <- TRUE

        # test
        if(length(opt$test) == 1L && opt$test == "default") {
            # ok, will be set later
        } else if(all(opt$test %in% c("none", "standard","yuan.bentler"))) {
            # nothing to do
        } else {
            stop("lavaan ERROR: `test' argument must one of \"none\", \"standard\" or \"yuan.bentler\" in the multilevel case")
        }

        # se
        if(opt$se == "default") {
            # ok, will be set later
        } else if(opt$se %in% c("none", "standard", "robust.huber.white",
                                "sandwich")) {
            # nothing to do
        } else if(opt$se == "robust") {
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
        if(opt$mimic == "Mplus" && !opt$.categorical &&
           opt$estimator %in% c("default", "ml", "mlr")) {
            # since version 5?
            opt$missing <- "ml"
            # check later if this is ok
        } else {
            opt$missing <- "listwise"
        }
    } else if(opt$missing %in% c("ml", "direct", "fiml")) {
        if(opt$.categorical) {
            stop("lavaan ERROR: missing = ", dQuote(opt$missing),
                 " not available in the categorical setting")
        }
        opt$missing <- "ml"
        if(opt$estimator %in% c("mlm", "mlmv", "gls", "wls", "wlsm", "wlsmv",
                                "uls", "ulsm", "ulsmv", "pml", "dls")) {
            stop("lavaan ERROR: missing=\"ml\" is not allowed for estimator ",
                 dQuote(toupper(opt$estimator)))
        }
    } else if(opt$missing %in% c("ml.x", "direct.x", "fiml.x")) {
        if(opt$.categorical) {
            stop("lavaan ERROR: missing = ", dQuote(opt$missing),
                 " not available in the categorical setting")
        }
        opt$missing <- "ml.x"
        if(opt$estimator %in% c("mlm", "mlmv", "gls", "wls", "wlsm", "wlsmv",
                                "uls", "ulsm", "ulsmv", "pml", "dls")) {
            stop("lavaan ERROR: missing=\"ml\" is not allowed for estimator ",
                 dQuote(toupper(opt$estimator)))
        }
    } else if(opt$missing %in% c("two.stage", "twostage", "two-stage",
                                 "two.step",  "twostep",  "two-step")) {
        opt$missing <- "two.stage"
        if(opt$.categorical) {
            stop("lavaan ERROR: missing=\"two.stage\" not available in the categorical setting")
        }
        if(opt$estimator %in% c("mlm", "mlmv", "gls", "wls", "wlsm", "wlsmv",
                                "uls", "ulsm", "ulsmv", "pml", "mml", "dls")) {
            stop("lavaan ERROR: missing=\"two.stage\" is not allowed for estimator ",
                 dQuote(toupper(opt$estimator)))
        }
    } else if(opt$missing %in% c("robust.two.stage", "robust.twostage",
                                 "robust.two-stage", "robust-two-stage",
                                 "robust.two.step",  "robust.twostep",
                                 "robust-two-step",
                                 "two.stage.robust", "twostage.robust",
                                 "two-stage.robust", "two-stage-robust",
                                 "two.step.robust",  "twostep.robust",
                                 "two-step-robust")) {
        opt$missing <- "robust.two.stage"
        if(opt$.categorical) {
            stop("lavaan ERROR: missing=\"robust.two.stage\" not available in the categorical setting")
        }
        if(opt$estimator %in% c("mlm", "mlmv", "gls", "wls", "wlsm", "wlsmv",
                                "uls", "ulsm", "ulsmv", "pml", "mml", "dls")) {
            stop("lavaan ERROR: missing=\"robust.two.stage\" is not allowed for estimator ",
                 dQuote(toupper(opt$estimator)))
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

        # new in 0.6-9: ALWAS h1.information = "unstructured"
        opt$h1.information <- c("unstructured", "unstructured")

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

    # bootstrap
    if(opt$se == "bootstrap") {
        opt$information[1] <- "observed"
        if(length(opt$information) > 1L && opt$information[2] == "default") {
            opt$information[2] <- "observed"
        }
        opt$bootstrap <- as.integer(opt$bootstrap)
        stopifnot(opt$bootstrap > 0L)
    }





    ##################################################################
    # ML and friends: MLF, MLM, MLMV, MLMVS, MLR                     #
    ##################################################################
    if(opt$estimator %in% c("ml", "mlf", "mlm", "mlmv", "mlmvs", "mlr")) {

        # set estimator
        opt$estimator <- "ML"

        # se
        if(opt$se == "bootstrap" &&
           orig.estimator %in% c("mlf", "mlm", "mlmv", "mlmvs", "mlr")) {
            stop("lavaan ERROR: use ML estimator for bootstrap")
        } else if(opt$se == "default") {
            if(orig.estimator %in% c("ml", "mlf")) {
                opt$se <- "standard"
            } else if(orig.estimator %in% c("mlm", "mlmv", "mlmvs")) {
                opt$se <- "robust.sem"
            } else if(orig.estimator == "mlr") {
                opt$se <- "robust.huber.white"
            }
        } else if(opt$se == "robust") {
            if(opt$missing %in% c("ml", "ml.x")) {
                opt$se <- "robust.huber.white"
            } else if(opt$missing == "two.stage") { # needed?
                opt$se <- "two.stage"
            } else if(opt$missing == "robust.two.stage") { # needed?
                opt$se <- "robust.two.stage"
            } else {
                opt$se <- "robust.sem"
            }
        }

        # information
        if(orig.estimator == "mlf") {
            if(opt$information[1] == "default") {
                opt$information[1] <- "first.order"
            }
            if(length(opt$information) > 1L &&
               opt$information[2] == "default") {
                opt$information[2] <- "first.order"
            }
        }

        # test
        if( !opt$test[1] == "none" ) {
            if(orig.estimator %in% c("ml", "mlf")) {
                if(opt$test[1] == "default") {
                    opt$test <- "standard"
                } else {
                    opt$test <- union("standard", opt$test)
                }
            } else if(orig.estimator == "mlm") {
                if(opt$test[1] == "default") {
                    opt$test <- "satorra.bentler"
                } else {
                    opt$test <- union("satorra.bentler", opt$test)
                }
            } else if(orig.estimator == "mlmv") {
                if(opt$test[1] == "default") {
                    opt$test <- "scaled.shifted"
                } else {
                    opt$test <- union("scaled.shifted", opt$test)
                }
            } else if(orig.estimator == "mlmvs") {
                if(opt$test[1] == "default") {
                    opt$test <- "mean.var.adjusted"
                } else {
                    opt$test <- union("mean.var.adjusted", opt$test)
                }
            } else if(orig.estimator == "mlr") {
                if(opt$mimic == "EQS") {
                    mlr.test <- "yuan.bentler"
                } else if(opt$mimic == "Mplus") {
                    mlr.test <- "yuan.bentler.mplus"
                } else {
                    mlr.test <- "yuan.bentler.mplus" # for now
                }
                if(opt$test[1] == "default") {
                    opt$test <- mlr.test
                } else {
                    opt$test <- union(mlr.test, opt$test)
                }
            }
        }


    ##################################################################
    # GLS                                                            #
    ##################################################################
    } else if(opt$estimator == "gls") {

        # estimator
        opt$estimator <- "GLS"

        # FIXME: catch categorical, clustered, ...

        # se
        if(opt$se == "default") {
            opt$se <- "standard"
        }

        # test
        if(opt$test[1] == "default") {
            opt$test <- "standard"
        }
        bad.idx <- which(!opt$test %in% c("standard", "none",
                                          "browne.residual.nt", # == standard
                                          "browne.residual.nt.model",
                                          "browne.residual.adf",
                                          "browne.residual.adf.model"))
        if(length(bad.idx) > 0L) {
            stop("lavaan ERROR: invalid value(s) in test= argument when estimator is GLS:\n\t\t",
                 paste(opt$test[bad.idx], collapse = " "), "\n")
        }

        # missing
        opt$missing <- "listwise" # also pairwise?


    ##################################################################
    # NTRLS (experimental)                                           #
    ##################################################################
    } else if(opt$estimator == "ntrls") {

        # optim.gradient
        opt$optim.gradien <- "numerical"

        # estimator
        opt$estimator <- "NTRLS"

        # se
        if(opt$se == "default") {
            opt$se <- "standard"
        }

        # test
        if(opt$test[1] == "default") {
            opt$test <- "standard"
        }
        bad.idx <- which(!opt$test %in% c("standard", "none",
                                          "browne.residual.nt",
                                          "browne.residual.nt.model",
                                          "browne.residual.adf",
                                          "browne.residual.adf.model"))
        if(length(bad.idx) > 0L) {
            stop("lavaan ERROR: invalid value(s) in test= argument when estimator is NTRLS:\n\t\t",
                 paste(opt$test[bad.idx], collapse = " "), "\n")
        }

        # missing
        opt$missing <- "listwise"

    ##################################################################
    # catML (experimental)                                           #
    ##################################################################
    } else if(opt$estimator == "catml") {

        # optim.gradient
        opt$optim.gradient <- "numerical" # for now

        # estimator
        opt$estimator <- "catML"

        # force correlation = TRUE, and categorical = FALSE
        opt$correlation <- TRUE
        opt$.categorical <- FALSE # we 'pretend' to have continuous data!

        # se
        if(opt$se == "default") {
            opt$se <- "robust.sem" # for now
        }

        # test
        if(opt$test[1] == "default") {
            opt$test <- "satorra.bentler"
        }

        # missing
        if(opt$missing %in% c("listwise", "pairwise")) {
            # nothing to do
        } else if(opt$missing == "default") {
            opt$missing <- "listwise"
        } else {
            stop("lavaan ERROR: missing argument should be listwise or pairwise if estimator is catML")
        }


    ##################################################################
    # WLS                                                            #
    ##################################################################
    } else if(opt$estimator == "wls") {

        # estimator
        opt$estimator <- "WLS"

        # se
        if(opt$se == "default") {
            opt$se <- "standard"
        }

        # test
        if(opt$test[1] == "default") {
            opt$test <- "standard"
        }
        bad.idx <- which(!opt$test %in% c("standard", "none",
                                          "browne.residual.nt",
                                          "browne.residual.nt.model",
                                          "browne.residual.adf", # == standard
                                          "browne.residual.adf.model"))
        if(length(bad.idx) > 0L) {
            stop("lavaan ERROR: invalid value(s) in test= argument when estimator is WLS:\n\t\t",
                 paste(opt$test[bad.idx], collapse = " "), "\n")
        }

        # missing
        #opt$missing <- "listwise" (could be pairwise)


    ##################################################################
    # DLS                                                            #
    ##################################################################
    } else if(opt$estimator == "dls") {

        # sample.cov.rescale
        if(is.logical(opt$sample.cov.rescale)) {
            # nothing to do
        } else if(opt$sample.cov.rescale == "default") {
            opt$sample.cov.rescale <- TRUE
        } else {
            stop("lavaan ERROR: sample.cov.rescale value must be logical.")
        }

        # estimator
        opt$estimator <- "DLS"

        # se
        if(opt$se == "default") {
            opt$se <- "robust.sem"
        }

        # test
        if(opt$test[1] == "default") {
            opt$test <- "satorra.bentler"
        }
        bad.idx <- which(!opt$test %in% c("standard", "none",
                                          "satorra.bentler",
                                          "browne.residual.nt", # == standard
                                          "browne.residual.nt.model",
                                          "browne.residual.adf",
                                          "browne.residual.adf.model"))
        if(length(bad.idx) > 0L) {
            stop("lavaan ERROR: invalid value(s) in test= argument when estimator is DLS:\n\t\t",
                 paste(opt$test[bad.idx], collapse = " "), "\n")
        }

        # always include "satorra.bentler"
        if(opt$test[1] %in% c("browne.residual.nt", "browne.residual.adf",
                              "browne.residual.nt.model",
                              "browne.residual.adf.model")) {
            opt$test <- union("satorra.bentler", opt$test)
        }

        # missing
        opt$missing <- "listwise"

        # estimator.args
        if(is.null(opt$estimator.args)) {
            opt$estimator.args <- list(dls.a = 1.0, dls.GammaNT = "model",
                                       dls.FtimesNmin1 = FALSE)
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
                opt$estimator.args$dls.GammaNT <- "model"
            } else {
                stopifnot(is.character(opt$estimator.args$dls.GammaNT))
                opt$estimator.args$dls.GammaNT <-
                    tolower(opt$estimator.args$dls.GammaNT)
                if(!opt$estimator.args$dls.GammaNT %in% c("sample", "model")) {
                    stop("lavaan ERROR: dls.GammaNT value in estimator.args must be either \"sample\" or \"model\".")
                }
            }
            if(is.null(opt$estimator.args$dls.FtimesNminus1)) {
                opt$estimator.args$dls.FtimesNminus1 <- FALSE
            } else {
                stopifnot(is.logical(opt$estimator.args$dls.FtimesNminus1))
            }
        }

        if(opt$estimator.args$dls.GammaNT == "sample") {
            if(opt$optim.method %in% c("nlminb", "gn")) {
                # nothing to do
            } else if(opt$optim.method == "default") {
                opt$optim.method <- "gn"
            } else {
                stop("lavaan ERROR: optim.method must be either nlminb or gn if estimator is DLS.")
            }
        } else {
            if(opt$optim.method %in% c("gn")) {
                # nothing to do
            } else if(opt$optim.method == "default") {
                opt$optim.method <- "gn"
            } else if(opt$optim.method == "nlminb") {
                opt$optim.gradient = "numerical"
            } else {
                stop("lavaan ERROR: optim.method must be either nlminb or gn if estimator is DLS.")
            }
        }

    ##################################################################
    # DWLS, WLSM, WLSMV, WLSMVS                                      #
    ##################################################################
    } else if(opt$estimator  %in% c("dwls", "wlsm", "wlsmv", "wlsmvs")) {

        # estimator
        opt$estimator <- "DWLS"

        # se
        if(opt$se == "bootstrap" &&
            orig.estimator %in% c("wlsm", "wlsmv", "wlsmvs")) {
            stop("lavaan ERROR: use (D)WLS estimator for bootstrap")
        } else if(opt$se == "default") {
            if(orig.estimator == "dwls") {
                opt$se <- "standard"
            } else {
                opt$se <- "robust.sem"
            }
        } else if(opt$se == "robust") {
            opt$se <- "robust.sem"
        }

        # test
        if( !opt$test[1] == "none" ) {
            if(orig.estimator == "dwls") {
                if(opt$test[1] == "default") {
                    opt$test <- "standard"
                } else {
                    opt$test <- union("standard", opt$test)
                }
            } else if(orig.estimator == "wlsm") {
                if(opt$test[1] == "default") {
                    opt$test <- "satorra.bentler"
                } else {
                    opt$test <- union("satorra.bentler", opt$test)
                }
            } else if(orig.estimator == "wlsmv") {
                if(opt$test[1] == "default") {
                    opt$test <- "scaled.shifted"
                } else {
                    opt$test <- union("scaled.shifted", opt$test)
                }
            } else if(orig.estimator == "wlsmvs") {
                if(opt$test[1] == "default") {
                    opt$test <- "mean.var.adjusted"
                } else {
                    opt$test <- union("mean.var.adjusted", opt$test)
                }
            }
        }


    ##################################################################
    # ULS, ULSM, ULSMV, ULSMVS                                       #
    ##################################################################
    } else if(opt$estimator %in% c("uls", "ulsm", "ulsmv", "ulsmvs")) {

        # estimator
        opt$estimator <- "ULS"

        # se
        if(opt$se == "bootstrap" &&
            orig.estimator %in% c("ulsm", "ulsmv", "ulsmvs")) {
            stop("lavaan ERROR: use ULS estimator for bootstrap")
        } else if(opt$se == "default") {
            if(orig.estimator == "uls") {
                opt$se <- "standard"
            } else {
                opt$se <- "robust.sem"
            }
        } else if(opt$se == "robust") {
            opt$se <- "robust.sem"
        }

        # test
        if( !opt$test[1] == "none" ) {
            if(orig.estimator == "uls") {
                if(opt$test[1] == "default") {
                    opt$test <- "standard"
                } else {
                    opt$test <- union("standard", opt$test)
                }
            } else if(orig.estimator == "ulsm") {
                if(opt$test[1] == "default") {
                    opt$test <- "satorra.bentler"
                } else {
                    opt$test <- union("satorra.bentler", opt$test)
                }
            } else if(orig.estimator == "ulsmv") {
                if(opt$test[1] == "default") {
                    opt$test <- "scaled.shifted"
                } else {
                    opt$test <- union("scaled.shifted", opt$test)
                }
            } else if(orig.estimator == "ulsmvs") {
                if(opt$test[1] == "default") {
                    opt$test <- "mean.var.adjusted"
                } else {
                    opt$test <- union("mean.var.adjusted", opt$test)
                }
            }
        }


    ##################################################################
    # PML                                                            #
    ##################################################################
    } else if(opt$estimator == "pml") {

        # estimator
        opt$estimator <- "PML"

        # se
        if(opt$se == "default") {
            opt$se <- "robust.huber.white"
        }

        # information
        opt$information[1] <- "observed"
        if(length(opt$information) > 1L &&
               opt$information[2] == "default") {
            opt$information[2] <- "observed"
        }
        if(length(opt$observed.information) > 1L &&
               opt$observed.information[2] == "default") {
            opt$observed.information[2] <- "hessian"
        }

        # test
        if(length(opt$test) > 1L) {
            stop("lavaan ERROR: only a single test statistic is allow when estimator is PML,")
        }
        if(!opt$test[1] == "none") {
            opt$test <- "mean.var.adjusted"
        }


    ##################################################################
    # FML - UMN                                                      #
    ##################################################################
    } else if(opt$estimator %in% c("fml","umn")) {

        # estimator
        opt$estimator <- "FML"

        # optim.gradient
        opt$optim.gradient <- "numerical"

        # se
        if(opt$se == "default") {
            opt$se <- "standard"
        }

        # information
        opt$information[1] <- "observed"
        if(length(opt$information) > 1L &&
               opt$information[2] == "default") {
            opt$information[2] <- "observed"
        }

        # test
        if(!opt$test[1] == "none") {
            opt$test <- "standard"
        }

    ##################################################################
    # REML                                                           #
    ##################################################################
    } else if(opt$estimator == "reml") {

        # estimator
        opt$estimator <- "REML"

        # se
        if(opt$se == "default") {
            opt$se <- "standard"
        }

        # information
        opt$information[1] <- "observed"
        if(length(opt$information) > 1L &&
               opt$information[2] == "default") {
            opt$information[2] <- "observed"
        }

        # test
        if(!opt$test[1] == "none") {
            opt$test <- "standard"
        }

        # missing
        opt$missing <- "listwise"

    ##################################################################
    # MML                                                            #
    ##################################################################
    } else if(opt$estimator %in% c("mml")) {

        # estimator
        opt$estimator <- "MML"

        # se
        if(opt$se == "default") {
            opt$se <- "standard"
        }

        # information
        opt$information[1] <- "observed"
        opt$meanstructure <- TRUE
        if(length(opt$information) > 1L &&
               opt$information[2] == "default") {
            opt$information[2] <- "observed"
        }

        # test
        opt$test <- "none"

        # link
        if(opt$link == "default") {
            #opt$link <- "logit"
            opt$link <- "probit"
        } else if(opt$link %in% c("logit","probit")) {
            # nothing to do
        } else {
            stop("lavaan ERROR: link must be `logit' or `probit'")
        }

        # parameterization
        if(opt$parameterization == "default") {
            opt$parameterization <- "mml"
        } else {
            stop("lavaan WARNING: parameterization argument is ignored if estimator = MML")
        }


    ##################################################################
    # FABIN                                                          #
    ##################################################################
    } else if(opt$estimator %in% c("fabin", "fabin2", "fabin3")) {
        # experimental, for cfa or sam only

        # estimator
        if(opt$estimator == "fabin") {
            opt$estimator <- "FABIN2"
        } else {
            opt$estimator <- toupper(opt$estimator)
        }

        # se
        if(opt$se == "default") {
            opt$se <- "none"
        }

        # test
        opt$test <- "none" # for now

        # missing
        opt$missing <- "listwise" # for now (until we have two-stage working)

        # options
        if(is.null(opt$estimator.args)) {
            opt$estimator.args <- list(thetapsi.method = "ULS",
                                       lambda.cutoff = 1.05)
        } else {
            if(is.null(opt$estimator.args$thetapsi.method)) {
                opt$estimator.args$thetapsi.method <- "GLS"
            } else {
                opt$estimator.args$thetapsi.method <-
                    toupper(opt$estimator.args$thetapsi.method)
                if(opt$estimator.args$thetapsi.method %in% c("ULS",
                                                             "GLS", "WLS")) {
                    if(opt$estimator.args$thetapsi.method == "WLS") {
                        opt$estimator.args$thetapsi.method <- "GLS"
                    }
                } else {
                    stop("lavaan ERROR: unknown value for estimator.args$thetapsi.method option: ", opt$estimator.args$thetapsi.method)
                }
            }
        }

        # brute-force override
        opt$optim.method <- "noniter"
        opt$start <- "simple"

    ##################################################################
    # NONE                                                           #
    ##################################################################
    } else if(opt$estimator == "none") {

        # se
        if(opt$se == "default") {
            opt$se <- "none"
        }

        # test
        if(opt$test[1] == "default") {
            opt$test <- "none"
        }


    } else {
        stop("lavaan ERROR: unknown value for estimator= argument: ",
             opt$estimator, "\n")
    }








    # optim.method - if still "default" at this point -> set to "nlminb"
    if(opt$optim.method == "default") {
        opt$optim.method <- "nlminb"
    }


    # special stuff for categorical
    if(opt$.categorical) {
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
    if(!opt$estimator %in% c("ML", "REML", "PML", "FML","NTRLS","catML")) {
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
       any(opt$test %in% c("satorra.bentler",
                       "yuan.bentler",
                       "yuan.bentler.mplus",
                       "mean.var.adjusted",
                       "scaled.shifted"))) {
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
        if(any(opt$test %in% c("satorra.bentler",
                           "yuan.bentler",
                           "yuan.bentler.mplus",
                           "mean.var.adjusted",
                           "scaled.shifted"))) {
            if(length(opt$test) > 1L) {
                opt$observed.information[2] <- "h1" # CHANGED in 0.6-6!
                if(any(opt$test == "yuan.bentler.mplus")) {
                    warning("observed.information for ALL test statistics is set to h1.")
                }
            } else {
                if(opt$estimator == "PML" ||
                   opt$test[1] == "yuan.bentler.mplus") {
                    opt$observed.information[2] <- "hessian"
                } else {
                    opt$observed.information[2] <- "h1" # CHANGED in 0.6-6!
                }
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
        } else if(opt$.categorical) {
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
        #if(opt$conditional.x && opt$fixed.x == FALSE && !opt$.multilevel) {
        if(opt$conditional.x && opt$fixed.x == FALSE) {
            stop("lavaan ERROR: fixed.x = FALSE is not supported when conditional.x = TRUE.")
        }
        if(opt$fixed.x && is.character(opt$start) && opt$start == "simple") {
            warning("lavaan WARNING: start = \"simple\" implies fixed.x = FALSE")
            opt$fixed.x <- FALSE
        }
    } else if(opt$fixed.x == "default") {
        if(opt$estimator %in% c("MML", "ML") &&
           (opt$mimic == "Mplus" || opt$mimic == "lavaan") &&
           is.character(opt$start) && opt$start != "simple") { # new in 0.6-12
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

    # add scaled.test to test (if not already there)
    if(opt$scaled.test != "standard") {
        if(length(opt$test) == 1L && opt$test[1] == "standard") {
            opt$test <- unique(c(opt$test, opt$scaled.test))
        } else {
            opt$test <- unique(c(opt$scaled.test, opt$test))
        }

        # make sure "standard" comes first
        standard.idx <- which(opt$test == "standard")[1]
        if(length(standard.idx) > 0L && standard.idx != 1L) {
            opt$test <- c("standard", opt$test[-standard.idx])
        }
    }


    # final check
    wrong.idx <- which(! opt$test %in% c("none", "standard", "satorra.bentler",
                            "yuan.bentler", "yuan.bentler.mplus",
                            "mean.var.adjusted", "scaled.shifted",
                            "browne.residual.adf", "browne.residual.nt",
                            "browne.residual.nt.model",
                            "browne.residual.adf.model",
                            "bollen.stine"))
    if(length(wrong.idx) > 0L) {
        txt <- c("invalid option(s) for test argument: ",
                 paste(dQuote(opt$test[wrong.idx]), collapse = " "), ". ",
                 "Possible options are: \"none\", \"standard\",
                 \"browne.residual.adf\", \"browne.residual.nt\",
                 \"browne.residual.adf.model\", \"browne.residual.nt.model\",
                 \"satorra.bentler\", \"yuan.bentler\", \"yuan.bentler.mplus\",
                 \"mean.var.adjusted\",
                 \"scaled.shifted\", or \"bollen.stine\"")
        stop(lav_txt2message(txt, header = "lavaan ERROR:"))
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
    } else if(opt$bounds == "default" ||
              opt$bounds == "wide") {
        opt$optim.bounds <- list(lower = c("ov.var", "lv.var", "loadings",
                                           "covariances"),
                                 upper = c("ov.var", "lv.var", "loadings",
                                           "covariances"),
                                 lower.factor = c(1.05, 1.0, 1.1, 1.0),
                                 upper.factor = c(1.20, 1.3, 1.1, 1.0),
                                 min.reliability.marker = 0.1,
                                 min.var.lv.endo = 0.005)
    } else if(opt$bounds == "standard") {
        opt$optim.bounds <- list(lower = c("ov.var", "lv.var", "loadings",
                                           "covariances"),
                                 upper = c("ov.var", "lv.var", "loadings",
                                           "covariances"),
                                 lower.factor = c(1.0, 1.0, 1.0, 0.999),
                                 upper.factor = c(1.0, 1.0, 1.0, 0.999),
                                 min.reliability.marker = 0.1,
                                 min.var.lv.endo = 0.005)
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
                     "tandem1", "tandem2", "none", "promax",
                     "oblimax", "bentler", "simplimax", "target", "pst")) {
        # nothing to do
    } else if(opt$rotation %in% c("cf-quartimax", "cf-varimax", "cf-equamax",
                            "cf-parsimax", "cf-facparsim")) {
        # nothing to do here; we need M/P to set cf.gamma
    } else if(opt$rotation %in% c("bi-quartimin", "biquartimin")) {
        opt$rotation <- "biquartimin"
    } else if(opt$rotation %in% c("bi-geomin", "bigeomin")) {
        opt$rotation <- "bigeomin"
    } else {
        txt <- c("Rotation method ", dQuote(opt$rotation), " not supported. ",
        "Supported rotation methods are: varimax, quartimax, orthomax, cf, ",
        "oblimin, quartimin, geomin, entropy, mccammon, infomax,", "promax",
        "tandem1, tandem2, oblimax, bentler, simplimax, target, pst, ",
        "crawford-ferguson,  cf-quartimax,  cf-varimax, cf-equamax, ",
        "cf-parsimax, cf-facparsim", "biquartimin", "bigeomin")
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

    # force orthogonal for some rotation algorithms
    if(opt$rotation %in% c("varimax", "entropy", "mccammon",
                           "tandem1", "tandem2") ) {
        opt$rotation.args$orthogonal <- TRUE
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
        # the default is "none", except for varimax and promax
        if(opt$rotation %in% c("varimax", "promax")) {
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
    } else if(opt$rotation %in% c("bi-geomin", "bigeomin",
                                  "bi-quartimin", "biquartimin")) {
        opt$rotation.args$order.lv.by <- "none"
    } else {
        stop("lavaan ERROR: rotation.args$order.lv.by should be \"none\",",
             " \"index\" or \"sumofsquares\".")
    }

    # no standard errors for promax (for now)...
    if(tolower(opt$rotation) == "promax") {
        opt$se <- "none"
        opt$rotation.args$algorithm <- "promax"
        opt$rotation.args$rstarts <- 0L
    }




    # correlation
    if(opt$correlation) {
        if(opt$missing == "ml") {
            stop("lavaan ERROR: correlation structures only work for complete data (for now).")
        }
        if(opt$.multilevel) {
            stop("lavaan ERROR: correlation structures only work for single-level data.")
        }
        if(opt$conditional.x) {
            stop("lavaan ERROR: correlation structures only work for conditional.x = FALSE (for now).")
        }
        if(opt$representation == "RAM") {
            stop("lavaan ERROR: correlation structures only work for representation = \"LISREL\".")
        }
        if(opt$fixed.x) {
            # first fix eliminate.pstar.idx in lav_mvnorm_information_expected()
            stop("lavaan ERROR: correlation structures only work for fixed.x = FALSE (for now).")
        }
    }



    # group.w.free
    #if(opt$group.w.free && opt$.categorical) {
    #    stop("lavaan ERROR: group.w.free = TRUE is not supported (yet) in the categorical setting.")
    #}

    # in order not to break semTools and blavaan, we restore categorical:
    opt$categorical <- opt$.categorical

    if(opt$debug) { cat("lavaan DEBUG: lavaanOptions OUT\n"); str(opt) }

    opt
}


# rename names of se values, and check for invalid values
lav_options_check_se <- function(opt = NULL) {

    # se must be a character string
    if(!is.character(opt$se)) {
        opt$se <- "default"
    }

    # unlike test=, se= should be a single character string
    if(length(opt$se) > 1L) {
        warning("lavaan WARNING: se= argument should be a single character string;\n\t\t  ",
                "Only the first entry (", dQuote(opt$se[1]),
                ") is used.", sep = "")
        opt$se <- opt$se[1]
    }

    # backwards compatibility (0.4 -> 0.5)
    if(opt$se == "robust.mlm") {
        opt$se <- "robust.sem"
    } else if(opt$se == "robust.mlr") {
        opt$se <- "robust.huber.white"
    } else if(opt$se == "first.order") {
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
    }

    # convenience
    else if(opt$se == "sandwich") {
        # (since 0.6-6)
        opt$se <- "robust.huber.white"
    } else if(opt$se == "boot") {
        opt$se <- "bootstrap"
    }

    # handle generic 'robust' (except clustered/multilvel)
    #else if(opt$se == "robust" && !opt$.clustered && !opt$.multilevel) {
    #    if(opt$missing %in% c("ml", "ml.x")) {
    #        opt$se <- "robust.huber.white"
    #    } else if(opt$missing == "two.stage") {
    #        opt$se <- "two.stage"
    #    } else if(opt$missing == "robust.two.stage") {
    #        opt$se <- "robust.two.stage"
    #    } else {
    #        # depends on estimator!
    #        opt$se <- "robust.sem"
    #    }
    #}

    # check for invalid names
    if(!opt$se %in% c("default", "none", "standard",
                      "bootstrap", "external",
                      "robust", "robust.sem", "robust.huber.white",
                      "two.stage", "robust.two.stage",
                      "robust.cluster", "robust.cluster.sem")) {
        stop("lavaan ERROR: invalid value in se= argument:\n\t\t",
             dQuote(opt$se))
    }

    # check for invalid names per estimator
    orig.estimator <- tolower(opt$estimator)

    # GLS, NTRLS, FML, UMN
    ok.flag <- TRUE
    if(orig.estimator %in% c("gls", "ntrls", "fml", "umn")) {
        ok.flag <- opt$se %in% c("default", "none", "standard",
                                 "bootstrap", "external")
    }

    # WLS, DLS, DWLS, WLSM, WLSMV, WLSMVS, ULS, ULSM, ULSMV, ULSMVS
    else if(orig.estimator %in% c("wls", "dls",
                                  "dwls", "wlsm", "wlsmv", "wlsmvs",
                                  "uls", "ulsm", "ulsmv", "ulsmvs")) {
        ok.flag <- opt$se %in% c("default", "none", "standard",
                                 "bootstrap", "external",
                                 "robust", "robust.sem")
    }

    # PML
    else if(orig.estimator  == "pml") {
        ok.flag <- opt$se %in% c("default", "none", "standard",
                                 "bootstrap", "external",
                                 "robust.huber.white")
    }

    # FABIN
    else if(orig.estimator %in% c("fabin", "fabin2", "fabin3")) {
        ok.flag <- opt$se %in% c("default", "none", "bootstrap", "external")
    }

    # OTHERS
    else if(orig.estimator %in% c("fml", "umn", "mml", "reml")) {
        ok.flag <- opt$se %in% c("default", "none", "standard", "external")
    }

    if(!ok.flag) {
        stop("lavaan ERROR: invalid value in se= argument for estimator ",
             toupper(orig.estimator), ":\n\t\t", dQuote(opt$se), sep = "")
    }

    opt
}


