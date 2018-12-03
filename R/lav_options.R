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
                std.lv             = FALSE,
                parameterization   = "default",

                auto.fix.first     = FALSE,
                auto.fix.single    = FALSE,
                auto.var           = FALSE,
                auto.cov.lv.x      = FALSE,
                auto.cov.y         = FALSE,
                auto.th            = FALSE,
                auto.delta         = FALSE,

                # seat belts
                safe.ov.var.ub     = FALSE,
                save.ov.var.lb     = FALSE,

                # full data
                std.ov             = FALSE,
                missing            = "default",

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
                likelihood             = "default",
                link                   = "default",
                representation         = "default",
                do.fit                 = TRUE,

                # inference
                information            = "default",
                h1.information         = "structured",
                #h1.information.se      = "structured",
                #h1.information.test    = "structured",
                se                     = "default",
                test                   = "default",
                bootstrap              = 1000L,
                observed.information   = "hessian",
                gamma.n.minus.one      = FALSE,
                #gamma.unbiased         = FALSE,

                # optimization
                control                = list(),
                optim.method           = "nlminb",
                optim.method.cor       = "nlminb",
                optim.force.converged  = FALSE,
                optim.gradient         = "analytic",
                optim.init_nelder_mead = FALSE,
                optim.var.transform    = "none",
                optim.parscale         = "none",
                optim.dx.tol           = 1e-04,
                em.iter.max            = 10000L,
                em.fx.tol              = 1e-08,
                em.dx.tol              = 1e-04,
                em.zerovar.offset      = 0.0001,

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

                # more models/info
                h1                     = TRUE,
                baseline               = TRUE,
                baseline.conditional.x.free.slopes = TRUE,
                implied                = TRUE,
                loglik                 = TRUE,

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
    opt$group <- opt.old$group
    opt$group.label <- opt.old$group.label
    opt$group.partial <- opt.old$group.partial
    opt$cluster <- opt.old$cluster


    # do.fit implies se="none and test="none" (unless not default)
    if(!opt$do.fit) {
        if(opt$se == "default") {
            opt$se <- "none"
        }
        if(opt$test == "default") {
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
    } else if(is.null(opt$group.equal) || nchar(opt$group.equal) == 0L) {
        if(opt$mimic == "Mplus" && !is.null(opt$group)) {
            if(opt$categorical) {
                opt$group.equal <- c("loadings", "thresholds")
            } else {
                opt$group.equal <- c("loadings", "intercepts")
            }
        } else {
            opt$group.equal <- character(0)
        }
    } else if(length(opt$group.equal) == 0) {
        # nothing to do
    } else if(all(opt$group.equal %in% c("loadings", "intercepts", "means",
                                         "regressions", "residuals",
                                         "residual.covariances", "thresholds",
                                         "lv.variances", "lv.covariances"))) {
        # nothing to do
    } else {
        stop("lavaan ERROR: unknown value for `group.equal' argument: ",
             opt$group.equal, "\n")
    }
    if(is.null(opt$group.partial) || nchar(opt$group.partial) == 0L) {
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
        if(opt$test == "default") {
            opt$test <- "yuan.bentler.mplus"
        } else if(opt$test %in% c("none", "standard",
                                  "satorra.bentler",
                                  "yuan.bentler","yuan.bentler.mplus")) {
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
        if(opt$information == "default") {
            if(opt$se == "robust.cluster") {
                opt$information <- "observed"
            } else {
                opt$information <- "expected"
            }
        }
        #} else if(opt$information %in% c("observed", "first.order")) {
        #    # nothing to do
        #} else {
        #    stop("lavaan ERROR: `information' argument must be \"observed\" in the multilevel case (for now)")
        #}

        #opt$fixed.x = FALSE
        #opt$control <- list(gradient = "numerical")
    }

    # multilevel
    # brute-force override (for now)
    if(opt$multilevel) {
        opt$meanstructure <- TRUE
        opt$missing <- "listwise"

        # test
        if(opt$test == "default") {
            opt$test <- "standard"
        } else if(opt$test %in% c("none", "standard","yuan.bentler")) {
            # nothing to do
        } else if(opt$se == "robust") {
            opt$test <- "yuan.bentler"
        } else {
            stop("lavaan ERROR: `test' argument must one of \"none\", \"standard\" or \"yuan.bentler\" in the multilevel case")
        }

        # se
        if(opt$se == "default") {
            opt$se <- "standard"
        } else if(opt$se %in% c("none", "standard", "robust.huber.white")) {
            # nothing to do
        } else if(opt$se == "robust.sem") {
            opt$se <- "robust.huber.white"
        } else {
            stop("lavaan ERROR: `se' argument must one of \"none\", \"standard\" or \"robust.huber.white\" in the multilevel case")
        }

        # information
        if(opt$information == "default") {
            opt$information <- "observed"
        }
        #} else if(opt$information %in% c("observed", "first.order")) {
        #    # nothing to do
        #} else {
        #    stop("lavaan ERROR: `information' argument must be \"observed\" in the multilevel case (for now)")
        #}

        #opt$fixed.x = FALSE
        #opt$control <- list(gradient = "numerical")
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
            stop("lavaan ERROR: missing=\"two.stage\" is not allowed for estimator MLM, MLMV, GLS, ULS, ULSM, ULSMV, DWLS, WLS, WLSM, WLSMV, PML, MML")
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
            stop("lavaan ERROR: missing=\"robust.two.stage\" is not allowed for estimator MLM, MLMV, GLS, ULS, ULSM, ULSMV, DWLS, WLS, WLSM, WLSMV, PML, MML")
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
    if(opt$test == "default") {
        if(opt$missing == "two.stage" ||
           opt$missing == "robust.two.stage") {
            opt$test <- "satorra.bentler"
        } else {
            opt$test <- "standard"
        }
    } else if(opt$test %in% c("none", "standard")) {
        # nothing to do
    } else if(opt$test == "satorra" ||
              opt$test == "sb" ||
              opt$test == "SB" ||
              opt$test == "satorra.bentler" ||
              opt$test == "satorra-bentler") {
        opt$test <- "satorra.bentler"
    } else if(opt$test == "yuan" ||
              opt$test == "yb" ||
              opt$test == "YB" ||
              opt$test == "yuan.bentler" ||
              opt$test == "yuan-bentler") {
        opt$test <- "yuan.bentler"
    } else if(opt$test == "yuan.bentler.mplus" ||
              opt$test == "yuan-bentler.mplus" ||
              opt$test == "yuan-bentler-mplus") {
        opt$test <- "yuan.bentler.mplus"
    } else if(opt$test == "m.adjusted" ||
              opt$test == "m" ||
              opt$test == "mean.adjusted" ||
              opt$test == "mean-adjusted") {
        opt$test <- "satorra.bentler"
    } else if(opt$test == "mean.var.adjusted" ||
              opt$test == "mean-var-adjusted" ||
              opt$test == "mv" ||
              opt$test == "second.order" ||
              opt$test == "satterthwaite" ||
              opt$test == "Satterthwaite" ||
              opt$test == "mv.adjusted") {
        opt$test <- "mean.var.adjusted"
    } else if(opt$test == "mplus6" ||
              opt$test == "scale.shift" ||
              opt$test == "scaled.shifted") {
        opt$test <- "scaled.shifted"
    } else if(opt$test == "bootstrap" ||
              opt$test == "boot" ||
              opt$test == "bollen.stine" ||
              opt$test == "bollen-stine") {
        opt$test <- "bollen.stine"
    } else {
        stop("lavaan ERROR: `test' argument must one of \"none\", \"standard\",
            \"satorra.bentler\", \"yuan.bentler\", \"yuan.bentler.mplus\",
            \"mean.var.adjusted\", \"scaled.shifted\",
            \"bollen.stine\", or \"bootstrap\"")
    }

    # check missing
    if(opt$missing %in% c("ml", "ml.x") && opt$se == "robust.sem") {
        warning("lavaan WARNING: missing will be set to ",
                    dQuote("listwise"), " for se = ",
                    dQuote(opt$se) )
        opt$missing <- "listwise"
    }
    if(opt$missing %in% c("ml", "ml.x") &&
       opt$test %in% c("satorra.bentler",
                       "mean.var.adjusted", "scaled.shifted")) {
        warning("lavaan WARNING: missing will be set to ",
                    dQuote("listwise"), " for test = ",
                    dQuote(opt$test) )
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
        if(opt$information == "default") {
            # for both two.stage and robust.two.stage
            opt$information <- "observed"
        } else if(opt$information == "first.order") {
            warning("lavaan WARNING: information will be set to ",
                     dQuote("observed"), " if missing = ",
                     dQuote(opt$missing) )
            opt$information <- "observed"
        }
        # observed.information (ALWAYS "h1" for now)
        opt$observed.information <- "h1"
        # test
        if(opt$test == "default" ||
           opt$test == "satorra.bentler") {
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
        # unless there is a group argument? (added since 0.4-10)
        # if(!is.null(opt$group)) opt$meanstructure <- TRUE
    } else {
        stop("lavaan ERROR: meanstructure must be TRUE, FALSE or \"default\"\n")
    }

    # estimator and se
    if(opt$se == "boot" || opt$se == "bootstrap") {
        opt$se <- "bootstrap"
        opt$information <- "observed"
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
            opt$information <- "first.order"
        } else if(opt$se == "observed") {
            opt$se <- "standard"
            opt$information <- "observed"
        } else if(opt$se == "expected") {
            opt$se <- "standard"
            opt$information <- "expected"
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
        if(opt$test != "none") {
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
        if(opt$se != "none" && opt$se != "external") opt$se <- "robust.sem"
        #if(!(opt$information %in% c("expected", "default"))) {
        #    warning("lavaan WARNING: information will be set to ",
        #            dQuote("expected"), " for estimator = ",
        #            dQuote(toupper(est.orig)) )
        #}
        #opt$information <- "expected"
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
            opt$information <- "first.order"
        }
    } else if(opt$estimator == "mlr") {
        opt$estimator <- "ML"
        if(opt$se == "bootstrap") {
            stop("lavaan ERROR: use ML estimator for bootstrap")
        }
        if(opt$se != "none" && opt$se != "external") opt$se <- "robust.huber.white"
        if(opt$test != "none" && opt$se != "external") {
            if(opt$mimic == "Mplus" || opt$test == "yuan.bentler.mplus") {
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
        if(!opt$test %in% c("standard","none")) {
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
        if(!opt$test %in% c("standard","none")) {
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
        if(!opt$test %in% c("standard","none")) {
            stop("lavaan ERROR: invalid value for `test' argument when estimator is WLS: ",
                 opt$test, "\n")
        }
        #opt$missing <- "listwise"
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
        if(!opt$test %in% c("standard","none","satorra.bentler",
                            "mean.adjusted",
                            "mean.var.adjusted","scaled.shifted")) {
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
        if(opt$test %in% c("mean.var.adjusted", "satorra.bentler",
                           "scaled.shifted")) {
            # nothing to do
        } else if(opt$test != "none") {
            opt$test <- "satorra.bentler"
        }
        #opt$missing <- "listwise"
    } else if(opt$estimator == "wlsmv") {
        opt$estimator <- "DWLS"
        if(opt$se == "bootstrap") {
            stop("lavaan ERROR: use (D)WLS estimator for bootstrap")
        }
        if(opt$se != "none" && opt$se != "external") opt$se <- "robust.sem"
        if(opt$test %in% c("mean.var.adjusted", "satorra.bentler",
                           "scaled.shifted")) {
            # nothing to do
        } else if(opt$test != "none") {
            opt$test <- "scaled.shifted"
        }
        #opt$missing <- "listwise"
    } else if(opt$estimator == "wlsmvs") {
        opt$estimator <- "DWLS"
        if(opt$se == "bootstrap") {
            stop("lavaan ERROR: use (D)WLS estimator for bootstrap")
        }
        if(opt$se != "none" && opt$se != "external") opt$se <- "robust.sem"
        if(opt$test %in% c("mean.var.adjusted", "satorra.bentler",
                           "scaled.shifted")) {
            # nothing to do
        } else if(opt$test != "none") {
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
        if(!opt$test %in% c("standard","none", "satorra.bentler",
                            "mean.adjusted",
                            "mean.var.adjusted","scaled.shifted")) {
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
        if(opt$test %in% c("mean.var.adjusted", "satorra.bentler",
                           "scaled.shifted")) {
            # nothing to do
        } else if(opt$test != "none") {
            opt$test <- "satorra.bentler"
        }
        #opt$missing <- "listwise"
    } else if(opt$estimator == "ulsmv") {
        opt$estimator <- "ULS"
        if(opt$se == "bootstrap") {
            stop("lavaan ERROR: use ULS estimator for bootstrap")
        }
        if(opt$se != "none" && opt$se != "external") opt$se <- "robust.sem"
        if(opt$test %in% c("mean.var.adjusted", "satorra.bentler",
                           "scaled.shifted")) {
            # nothing to do
        } else if(opt$test != "none") {
            opt$test <- "scaled.shifted"
        }
        #opt$missing <- "listwise"
    } else if(opt$estimator == "ulsmvs") {
        opt$estimator <- "ULS"
        if(opt$se == "bootstrap") {
            stop("lavaan ERROR: use ULS estimator for bootstrap")
        }
        if(opt$se != "none" && opt$se != "external") opt$se <- "robust.sem"
        if(opt$test %in% c("mean.var.adjusted", "satorra.bentler",
                           "scaled.shifted")) {
            # nothing to do
        } else if(opt$test != "none") {
            opt$test <- "mean.var.adjusted"
        }
        #opt$missing <- "listwise"
    } else if(opt$estimator == "pml") {
        opt$estimator <- "PML"
        opt$information <- "observed"
        if(opt$se == "default")
            opt$se <- "robust.huber.white"
        if(opt$test != "none") opt$test <- "mean.var.adjusted"
        #opt$missing <- "listwise"
    } else if(opt$estimator %in% c("fml","umn")) {
        opt$estimator <- "FML"
        opt$information <- "observed"
        if(opt$se == "default")
            opt$se <- "standard"
        if(opt$test != "none") opt$test <- "standard"
        #opt$missing <- "listwise"
    } else if(opt$estimator == "reml") {
        opt$estimator <- "REML"
        opt$information <- "observed"
        if(opt$se == "default")
            opt$se <- "standard"
        if(opt$test != "none") opt$test <- "standard"
        opt$missing <- "listwise"
    } else if(opt$estimator %in% c("mml")) {
        opt$estimator <- "MML"
        opt$information <- "observed"
        if(opt$se == "default")
            opt$se <- "standard"
        if(opt$test == "default")
            opt$test <- "none"
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
        if(opt$test == "default") {
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
        } else {
            warning("sample.cov.rescale argument is only relevant if estimator = ML")
        }
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

    # information
    if(opt$information == "default") {
        if(opt$missing %in% c("ml", "ml.x") ||
           opt$se == "robust.huber.white"   ||
           opt$se == "first.order") {
           #nchar(opt$constraints) > 0L) {
            opt$information <- "observed"
        } else {
            opt$information <- "expected"
        }
    } else if(opt$information %in% c("observed", "expected", "first.order")) {
        # nothing to do
    } else {
        stop("lavaan ERROR: information must be either \"expected\", \"observed\", or \"first.order\"\n")
    }

    if(opt$h1.information == "structured" ||
       opt$h1.information == "unstructured") {
        # nothing to do
    } else {
        stop("lavaan ERROR: h1.information must be either \"structured\" or \"unstructured\"\n")
    }
    #if(opt$h1.information.test == "structured" ||
    #   opt$h1.information.test == "unstructured") {
    #    # nothing to do
    #} else {
    #    stop("lavaan ERROR: h1.information.se must be either \"structured\" or \"unstructured\"\n")
    #}

    # check information if estimator is uls/wls and friends
    if(opt$estimator %in% c("ULS", "WLS", "DWLS")) {
        if(opt$information != "expected") {
            warning("lavaan WARNING: information will be set to ",
                    dQuote("expected"), " for estimator = ", dQuote(opt$estimator))
            opt$information <- "expected"
        }
        #if(opt$h1.information != "structured") { # default value
        #    warning("lavaan WARNING: h1.information is not used if estimator = ", dQuote(opt$estimator))
        opt$h1.information <- "unstructured"
        #}
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

    if(opt$debug) { cat("lavaan DEBUG: lavaanOptions OUT\n"); str(opt) }

    opt
}
