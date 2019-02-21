# model estimation
lav_model_estimate <- function(lavmodel       = NULL,
                               lavpartable    = NULL, # for parscale = "stand"
                               lavh1          = NULL, # for multilevel + parsc
                               lavsamplestats = NULL,
                               lavdata        = NULL,
                               lavoptions     = NULL,
                               lavcache       = list(),
                               do.fit         = TRUE) {

    estimator     <- lavoptions$estimator
    verbose       <- lavoptions$verbose
    debug         <- lavoptions$debug
    ngroups       <- lavsamplestats@ngroups

    if(lavsamplestats@missing.flag || estimator == "PML") {
        group.weight <- FALSE
    } else {
        group.weight <- TRUE
    }

    # temp test
    if(lavoptions$partrace) {
        # fx + parameter values
        PENV <- new.env()
        PENV$PARTRACE <- matrix(NA, nrow=0, ncol=lavmodel@nx.free + 1L)
    }

    # starting values (ignoring equality constraints)
    x.unpack <- lav_model_get_parameters(lavmodel)

    # 1. parameter scaling (to handle data scaling, not parameter scaling)
    parscale <- rep(1.0, length(x.unpack))

    # for < 0.6 compatibility
    if(is.null(lavoptions$optim.parscale)) {
        lavoptions$optim.parscale <- "none"
    }

    if(lavoptions$optim.parscale == "none") {
        # do nothing, but still set SCALE, as before
    } else if(lavoptions$optim.parscale %in% c("stand", "st", "standardize",
                                               "standarized", "stand.all")) {
        # rescale parameters as if the data was standardized
        # new in 0.6-2
        #
        # FIXME: this works well, as long as the variances of the
        #        latent variables (which we do not know) are more or less
        #        equal to 1.0 (eg std.lv = TRUE)
        #
        #        Once we have better estimates of those variances, we could
        #        use them to set the scale
        #

        if(lavdata@nlevels > 1L) {
            if(length(lavh1) > 0L) {
                OV.VAR <- lapply(lavh1$implied$cov, diag)
            } else {
                OV.VAR <- lapply(do.call(c, lapply(lavdata@Lp, "[[", "ov.idx")),
                                     function(x) rep(1, length(x) ))
            }
        } else {
            if(lavoptions$conditional.x) {
                OV.VAR <- lavsamplestats@res.var
            } else {
                OV.VAR <- lavsamplestats@var
            }
        }

        if(lavoptions$std.lv) {
            parscale <- lav_standardize_all(lavobject = NULL,
                            est = rep(1, length(lavpartable$lhs)),
                            est.std = rep(1, length(lavpartable$lhs)),
                            cov.std = FALSE, ov.var = OV.VAR,
                            lavmodel = lavmodel, lavpartable = lavpartable,
                            cov.x = lavsamplestats@cov.x)
        } else {
            # needs good estimates for lv variances!
            parscale <- lav_standardize_all(lavobject = NULL,
                            est = rep(1, length(lavpartable$lhs)),
                            #est.std = rep(1, length(lavpartable$lhs)),
                            # here, we use whatever the starting values are
                            # for the latent variances...
                            cov.std = FALSE, ov.var = OV.VAR,
                            lavmodel = lavmodel, lavpartable = lavpartable,
                            cov.x = lavsamplestats@cov.x)
        }

        # in addition, take sqrt for variance parameters
        var.idx <- which(lavpartable$op == "~~" &
                         lavpartable$lhs == lavpartable$rhs)
        if(length(var.idx) > 0L) {
            parscale[var.idx] <- sqrt(abs(parscale[var.idx]))
        }
        parscale <- parscale[ lavpartable$free > 0 ]
    }
    # parscale should obey the equality constraints
    if(lavmodel@eq.constraints && lavoptions$optim.parscale != "none") {
        # pack
        p.pack <- as.numeric( (parscale - lavmodel@eq.constraints.k0) %*%
                               lavmodel@eq.constraints.K )
        # unpack
        parscale <- as.numeric(lavmodel@eq.constraints.K %*% p.pack) +
                               lavmodel@eq.constraints.k0
    }
    if(debug) {
        cat("parscale = ", parscale, "\n")
    }
    z.unpack <- x.unpack * parscale

    # 2. pack (apply equality constraints)
    if(lavmodel@eq.constraints) {
        z.pack <- as.numeric( (z.unpack - lavmodel@eq.constraints.k0) %*%
                              lavmodel@eq.constraints.K )
    } else {
        z.pack <- z.unpack
    }

    # 3. transform (already constrained) variances to standard deviations?
    # TODO
    #if(lavoptions$optim.var.transform == "sqrt" &&
    #       length(lavmodel@x.free.var.idx) > 0L) {
    #    # transforming variances using atan (or another sigmoid function?)
    #    # FIXME: better approach?
    #    #start.x[lavmodel@x.free.var.idx] <-
    #    #    atan(start.x[lavmodel@x.free.var.idx])
    #    start.x[lavmodel@x.free.var.idx] <-
    #        sqrt(start.x[lavmodel@x.free.var.idx]) # assuming positive var
    #}

    # final starting values for optimizer
    start.x <- z.pack
    if(debug) {
        cat("start.x = ", start.x, "\n")
    }


    # bounds? (new in 0.6-2)
    if(is.null(lavpartable$lower)) {
        lower <- -Inf
    } else {
        lower <- lavpartable$lower[ lavpartable$free > 0 ]
        if(lavmodel@eq.constraints) {
            # collect non -Inf elements
            not.inf.val <- lower[is.finite(lower)]

            lowerb <- lower
            lowerb[lower == -Inf] <- -999999

            # pack
            l.pack <- as.numeric( (lowerb - lavmodel@eq.constraints.k0) %*%
                                  lavmodel@eq.constraints.K )
            lower <- l.pack
            lower[!(lower %in% not.inf.val)] <- -Inf
        }
    }
    if(is.null(lavpartable$upper)) {
        upper <- +Inf
    } else {
        upper <- lavpartable$upper[ lavpartable$free > 0 ]
        if(lavmodel@eq.constraints) {
            # collect non -Inf elements
            not.inf.val <- upper[is.finite(upper)]

            upperb <- upper
            upperb[upper == Inf] <- 999999

            # pack
            l.pack <- as.numeric( (upperb - lavmodel@eq.constraints.k0) %*%
                                  lavmodel@eq.constraints.K )
            upper <- l.pack
            upper[!(upper %in% not.inf.val)] <- Inf
        }
    }




    # function to be minimized
    objective_function <- function(x, verbose = FALSE, infToMax = FALSE) {

        # 3. standard deviations to variances
        # WARNING: x is still packed here!
        #if(lavoptions$optim.var.transform == "sqrt" &&
        #   length(lavmodel@x.free.var.idx) > 0L) {
        #    #x[lavmodel@x.free.var.idx] <- tan(x[lavmodel@x.free.var.idx])
        #    x.var <- x[lavmodel@x.free.var.idx]
        #    x.var.sign <- sign(x.var)
        #    x[lavmodel@x.free.var.idx] <- x.var.sign * (x.var * x.var) # square!
        #}

        # 2. unpack
        if(lavmodel@eq.constraints) {
            x <- as.numeric(lavmodel@eq.constraints.K %*% x) +
                            lavmodel@eq.constraints.k0
        }

        # 1. unscale
        x <- x / parscale

        # update GLIST (change `state') and make a COPY!
        GLIST <- lav_model_x2GLIST(lavmodel, x = x)

        fx <- lav_model_objective(lavmodel       = lavmodel,
                                  GLIST          = GLIST,
                                  lavsamplestats = lavsamplestats,
                                  lavdata        = lavdata,
                                  lavcache       = lavcache,
                                  verbose        = verbose,
                                  forcePD        = FALSE)

        # only for PML: divide by N (to speed up convergence)
        if(estimator == "PML") {
            fx <- fx / lavsamplestats@ntotal
        }



        if(debug || verbose) {
            cat("Objective function  = ", sprintf("%18.16f", fx), "\n", sep="")
        }
        if(debug) {
            #cat("Current unconstrained parameter values =\n")
            #tmp.x <- lav_model_get_parameters(lavmodel, GLIST=GLIST, type="unco")
            #print(tmp.x); cat("\n")
            cat("Current free parameter values =\n"); print(x); cat("\n")
        }

        if(lavoptions$partrace) {
            PENV$PARTRACE <- rbind(PENV$PARTRACE, c(fx, x))
        }

        # for L-BFGS-B
        #if(infToMax && is.infinite(fx)) fx <- 1e20
        if(!is.finite(fx)) {
            fx <- 1e20
        }

        fx
    }

    gradient_function <- function(x, verbose = FALSE, infToMax = FALSE) {

        # transform variances back
        #if(lavoptions$optim.var.transform == "sqrt" &&
        #   length(lavmodel@x.free.var.idx) > 0L) {
        #    #x[lavmodel@x.free.var.idx] <- tan(x[lavmodel@x.free.var.idx])
        #    x.var <- x[lavmodel@x.free.var.idx]
        #    x.var.sign <- sign(x.var)
        #    x[lavmodel@x.free.var.idx] <- x.var.sign * (x.var * x.var) # square!
        #}

        # 2. unpack
        if(lavmodel@eq.constraints) {
            x <- as.numeric(lavmodel@eq.constraints.K %*% x) +
                            lavmodel@eq.constraints.k0
        }

        # 1. unscale
        x <- x / parscale

        # update GLIST (change `state') and make a COPY!
        GLIST <- lav_model_x2GLIST(lavmodel, x = x)

        dx <- lav_model_gradient(lavmodel       = lavmodel,
                                 GLIST          = GLIST,
                                 lavsamplestats = lavsamplestats,
                                 lavdata        = lavdata,
                                 lavcache       = lavcache,
                                 type           = "free",
                                 group.weight   = group.weight, ### check me!!
                                 verbose        = verbose,
                                 forcePD        = TRUE)

        if(debug) {
            cat("Gradient function (analytical) =\n"); print(dx); cat("\n")
        }

        # 1. scale (note: divide, not multiply!)
        dx <- dx / parscale

        # 2. pack
        if(lavmodel@eq.constraints) {
            dx <- as.numeric( dx %*% lavmodel@eq.constraints.K )
        }

        # 3. transform variances back
        #if(lavoptions$optim.var.transform == "sqrt" &&
        #   length(lavmodel@x.free.var.idx) > 0L) {
        #    x.var <- x[lavmodel@x.free.var.idx] # here in 'var' metric
        #    x.var.sign <- sign(x.var)
        #    x.var <- abs(x.var)
        #    x.sd <- sqrt(x.var)
        #    dx[lavmodel@x.free.var.idx] <-
        #        ( 2 * x.var.sign * dx[lavmodel@x.free.var.idx] * x.sd )
        #}

        # only for PML: divide by N (to speed up convergence)
        if(estimator == "PML") {
            dx <- dx / lavsamplestats@ntotal
        }

        if(debug) {
            cat("Gradient function (analytical, after eq.constraints.K) =\n"); print(dx); cat("\n")
        }

        dx
    }

    gradient_function_numerical <- function(x, verbose=FALSE) {

        # NOTE: no need to 'tranform' anything here (var/eq)
        # this is done anyway in objective_function

        # numerical approximation using the Richardson method
        npar <- length(x)
        h <- 10e-6
        dx <- numeric( npar )

        ## FIXME: call lav_model_objective directly!!
        for(i in 1:npar) {
            x.left <- x.left2 <- x.right <- x.right2 <- x
            x.left[i]  <- x[i] - h; x.left2[i]  <- x[i] - 2*h
            x.right[i] <- x[i] + h; x.right2[i] <- x[i] + 2*h
            fx.left   <- objective_function(x.left)
            fx.left2  <- objective_function(x.left2)
            fx.right  <- objective_function(x.right)
            fx.right2 <- objective_function(x.right2)
            dx[i] <- (fx.left2 - 8*fx.left + 8*fx.right - fx.right2)/(12*h)
        }

        #dx <- lavGradientC(func=objective_function, x=x)
        # does not work if pnorm is involved... (eg PML)

        if(debug) {
            cat("Gradient function (numerical) =\n"); print(dx); cat("\n")
        }

        dx
    }


    # check if the initial values produce a positive definite Sigma
    # to begin with -- but only for estimator="ML"
    #if(estimator %in% c("ML","PML","FML","MML")) {
    if(estimator %in% c("ML","FML","MML")) {
        Sigma.hat <- computeSigmaHat(lavmodel, extra=TRUE, debug=lavoptions$debug)
        for(g in 1:ngroups) {
            if(!attr(Sigma.hat[[g]], "po")) {
                group.txt <- ifelse(ngroups > 1,
                                    paste(" in group ",g,".",sep=""), ".")
                if(debug) print(Sigma.hat[[g]])
                stop("lavaan ERROR: initial model-implied matrix (Sigma) is not positive definite;\n  check your model and/or starting parameters", group.txt)
                x <- start.x
                fx <- as.numeric(NA)
                attr(fx, "fx.group") <- rep(as.numeric(NA), ngroups)
                attr(x, "converged")  <- FALSE
                attr(x, "iterations") <- 0L
                attr(x, "control")    <- lavoptions@control
                attr(x, "fx")         <- fx
                return(x)
            }
        }
    }


    # parameter scaling
    # FIXME: what is the best way to set the scale??
    # current strategy: if startx > 1.0, we rescale by using
    # 1/startx
    SCALE <- rep(1.0, length(start.x))
    if(lavoptions$optim.parscale == "none") {
        idx <- which(abs(start.x) > 1.0)
        if(length(idx) > 0L) {
            SCALE[idx] <- abs(1.0/start.x[idx])
        }
    }
    if(debug) {
        cat("SCALE = ", SCALE, "\n")
    }


    # first some nelder mead steps? (default = FALSE)
    INIT_NELDER_MEAD <- lavoptions$optim.init_nelder_mead

    # gradient: analytic, numerical or NULL?
    if(is.character(lavoptions$optim.gradient)) {
        if(lavoptions$optim.gradient %in% c("analytic","analytical")) {
            GRADIENT <- gradient_function
        } else if(lavoptions$optim.gradient %in% c("numerical", "numeric")) {
            GRADIENT <- gradient_function_numerical
        } else if(lavoptions$optim.gradient %in% c("NULL", "null")) {
            GRADIENT <- NULL
        } else {
            warning("lavaan WARNING: gradient should be analytic, numerical or NULL")
        }
    } else if(is.logical(lavoptions$optim.gradient)) {
        if(lavoptions$optim.gradient) {
            GRADIENT <- gradient_function
        } else {
            GRADIENT <- NULL
        }
    } else if(is.null(lavoptions$optim.gradient)) {
        GRADIENT <- gradient_function
    }


    # default optimizer
    if(length(lavmodel@ceq.nonlinear.idx) == 0L &&
       length(lavmodel@cin.linear.idx)    == 0L &&
       length(lavmodel@cin.nonlinear.idx) == 0L) {
        if(is.null(lavoptions$optim.method)) {
            OPTIMIZER <- "NLMINB"
            #OPTIMIZER <- "BFGS"  # slightly slower, no bounds; better scaling!
            #OPTIMIZER <- "L-BFGS-B"  # trouble with Inf values for fx!
        } else {
            OPTIMIZER <- toupper(lavoptions$optim.method)
            stopifnot(OPTIMIZER %in% c("NLMINB0", "NLMINB1", "NLMINB2",
                      "NLMINB", "BFGS", "L-BFGS-B", "NONE"))
            if(OPTIMIZER == "NLMINB1") {
                OPTIMIZER <- "NLMINB"
            }
        }
    } else {
        OPTIMIZER <- "NLMINB.CONSTR"
    }

    if(INIT_NELDER_MEAD) {
        if(verbose) cat("Initial Nelder-Mead step:\n")
        trace <- 0L; if(verbose) trace <- 1L
        optim.out <- optim(par=start.x,
                           fn=objective_function,
                           method="Nelder-Mead",
                           #control=list(maxit=10L,
                           #             parscale=SCALE,
                           #             trace=trace),
                           hessian=FALSE,
                           verbose=verbose)
        cat("\n")
        start.x <- optim.out$par
    }



    if(OPTIMIZER == "NLMINB0") {
        if(verbose) cat("Quasi-Newton steps using NLMINB0 (no analytic gradient):\n")
        #if(debug) control$trace <- 1L;
        control.nlminb <- list(eval.max=20000L,
                               iter.max=10000L,
                               trace=0L,
                               #abs.tol=1e-20, ### important!! fx never negative
                               abs.tol=(.Machine$double.eps * 10),
                               rel.tol=1e-10,
                               #step.min=2.2e-14, # in =< 0.5-12
                               step.min=1.0, # 1.0 in < 0.5-21
                               step.max=1.0,
                               x.tol=1.5e-8,
                               xf.tol=2.2e-14)
        control.nlminb <- modifyList(control.nlminb, lavoptions$control)
        control <- control.nlminb[c("eval.max", "iter.max", "trace",
                                    "step.min", "step.max",
                                    "abs.tol", "rel.tol", "x.tol", "xf.tol")]
        #cat("DEBUG: control = "); print(str(control.nlminb)); cat("\n")
        optim.out <- nlminb(start=start.x,
                            objective=objective_function,
                            gradient=NULL,
                            lower=lower,
                            upper=upper,
                            control=control,
                            scale=SCALE,
                            verbose=verbose)
        if(verbose) {
            cat("convergence status (0=ok): ", optim.out$convergence, "\n")
            cat("nlminb message says: ", optim.out$message, "\n")
            cat("number of iterations: ", optim.out$iterations, "\n")
            cat("number of function evaluations [objective, gradient]: ",
                optim.out$evaluations, "\n")
        }

        iterations <- optim.out$iterations
        x          <- optim.out$par
        if(optim.out$convergence == 0) {
            # optimizer is happy, but do we have a zero gradient
            # check 'unscaled' gradient (for now)
            grad <- gradient_function_numerical(x)
            if( all(abs(grad) < lavoptions$optim.dx.tol) ) {
                converged <- TRUE
            } else {
                # warning here?
                converged <- FALSE
            }
        } else {
            converged <- FALSE
        }

    } else if(OPTIMIZER == "NLMINB") {
        if(verbose) cat("Quasi-Newton steps using NLMINB:\n")
        #if(debug) control$trace <- 1L;
        control.nlminb <- list(eval.max=20000L,
                               iter.max=10000L,
                               trace=0L,
                               #abs.tol=1e-20, ### important!! fx never negative
                               abs.tol=(.Machine$double.eps * 10),
                               rel.tol=1e-10,
                               #step.min=2.2e-14, # in =< 0.5-12
                               step.min=1.0, # 1.0 in < 0.5-21
                               step.max=1.0,
                               x.tol=1.5e-8,
                               xf.tol=2.2e-14)
        control.nlminb <- modifyList(control.nlminb, lavoptions$control)
        control <- control.nlminb[c("eval.max", "iter.max", "trace",
                                    "step.min", "step.max",
                                    "abs.tol", "rel.tol", "x.tol", "xf.tol")]
        #cat("DEBUG: control = "); print(str(control.nlminb)); cat("\n")
        optim.out <- nlminb(start=start.x,
                            objective=objective_function,
                            gradient=GRADIENT,
                            lower=lower,
                            upper=upper,
                            control=control,
                            scale=SCALE,
                            verbose=verbose)
        if(verbose) {
            cat("convergence status (0=ok): ", optim.out$convergence, "\n")
            cat("nlminb message says: ", optim.out$message, "\n")
            cat("number of iterations: ", optim.out$iterations, "\n")
            cat("number of function evaluations [objective, gradient]: ",
                optim.out$evaluations, "\n")
        }

        iterations <- optim.out$iterations
        x          <- optim.out$par
        if(optim.out$convergence == 0) {
            # optimizer is happy, but do we have a zero gradient
            # check 'unscaled' gradient (for now)

            # should we treat bounds differently?


            grad <- gradient_function(x)
            if( all(abs(grad) < lavoptions$optim.dx.tol) ) {
                converged <- TRUE
            } else {
                # warning here?
                converged <- FALSE
            }
        } else {
            converged <- FALSE
        }
    } else if(OPTIMIZER == "BFGS") {

        # warning: Bollen example with estimator=GLS does NOT converge!
        # (but WLS works!)
        # - BB.ML works too

        control.bfgs <- list(trace=0L, fnscale=1,
                             parscale=SCALE, ## or not?
                             ndeps=1e-3,
                             maxit=10000,
                             abstol=1e-20,
                             reltol=1e-10,
                             REPORT=1L)
        control.bfgs <- modifyList(control.bfgs, lavoptions$control)
        control <- control.bfgs[c("trace", "fnscale", "parscale", "ndeps",
                                  "maxit", "abstol", "reltol", "REPORT")]
        #trace <- 0L; if(verbose) trace <- 1L
        optim.out <- optim(par=start.x,
                           fn=objective_function,
                           gr=GRADIENT,
                           method="BFGS",
                           control=control,
                           hessian=FALSE,
                           verbose=verbose)
        if(verbose) {
            cat("convergence status (0=ok): ", optim.out$convergence, "\n")
            cat("optim BFGS message says: ", optim.out$message, "\n")
            #cat("number of iterations: ", optim.out$iterations, "\n")
            cat("number of function evaluations [objective, gradient]: ",
                optim.out$counts, "\n")
        }

        #iterations <- optim.out$iterations
        iterations <- optim.out$counts[1]
        x          <- optim.out$par
        if(optim.out$convergence == 0L) {
            converged <- TRUE
        } else {
            converged <- FALSE
        }
    } else if(OPTIMIZER == "L-BFGS-B") {

        # warning, does not cope with Inf values!!

        control.lbfgsb <- list(trace=0L, fnscale=1,
                               parscale=SCALE, ## or not?
                               ndeps=1e-3,
                               maxit=10000,
                               REPORT=1L,
                               lmm=5L,
                               factr=1e7,
                               pgtol=0)
        control.lbfgsb <- modifyList(control.lbfgsb, lavoptions$control)
        control <- control.lbfgsb[c("trace", "fnscale", "parscale",
                                    "ndeps", "maxit", "REPORT", "lmm",
                                    "factr", "pgtol")]
        optim.out <- optim(par=start.x,
                           fn=objective_function,
                           gr=GRADIENT,
                           method="L-BFGS-B",
                           lower=lower,
                           upper=upper,
                           control=control,
                           hessian=FALSE,
                           verbose=verbose,
                           infToMax=TRUE)
        if(verbose) {
            cat("convergence status (0=ok): ", optim.out$convergence, "\n")
            cat("optim L-BFGS-B message says: ", optim.out$message, "\n")
            #cat("number of iterations: ", optim.out$iterations, "\n")
            cat("number of function evaluations [objective, gradient]: ",
                optim.out$counts, "\n")
        }

        #iterations <- optim.out$iterations
        iterations <- optim.out$counts[1]
        x          <- optim.out$par
        if(optim.out$convergence == 0L) {
            converged <- TRUE
        } else {
            converged <- FALSE
        }
    } else if(OPTIMIZER == "NLMINB.CONSTR") {

        ocontrol <- list(verbose=verbose)
        if(!is.null(lavoptions$control$control.outer)) {
            ocontrol <- c(lavoptions$control$control.outer, verbose=verbose)
        }
        control.nlminb <- list(eval.max=20000L,
                               iter.max=10000L,
                               trace=0L,
                               #abs.tol=1e-20,
                               abs.tol=(.Machine$double.eps * 10),
                               rel.tol=1e-9) # 1e-10 seems 'too strict'
        control.nlminb <- modifyList(control.nlminb, lavoptions$control)
        control <- control.nlminb[c("eval.max", "iter.max", "trace",
                                    "abs.tol", "rel.tol")]
        cin <- cin.jac <- ceq <- ceq.jac <- NULL
        if(!is.null(body(lavmodel@cin.function))) cin     <- lavmodel@cin.function
        if(!is.null(body(lavmodel@cin.jacobian))) cin.jac <- lavmodel@cin.jacobian
        if(!is.null(body(lavmodel@ceq.function))) ceq     <- lavmodel@ceq.function
        if(!is.null(body(lavmodel@ceq.jacobian))) ceq.jac <- lavmodel@ceq.jacobian
        trace <- FALSE; if(verbose) trace <- TRUE
        optim.out <- nlminb.constr(start = start.x,
                                   objective=objective_function,
                                   gradient=GRADIENT,
                                   control=control,
                                   scale=SCALE,
                                   verbose=verbose,
                                   cin = cin, cin.jac = cin.jac,
                                   ceq = ceq, ceq.jac = ceq.jac,
                                   control.outer = ocontrol
                                  )
        if(verbose) {
            cat("convergence status (0=ok): ", optim.out$convergence, "\n")
            cat("nlminb.constr message says: ", optim.out$message, "\n")
            cat("number of outer iterations: ", optim.out$outer.iterations, "\n")
            cat("number of inner iterations: ", optim.out$iterations, "\n")
            cat("number of function evaluations [objective, gradient]: ",
                optim.out$evaluations, "\n")
        }

        iterations <- optim.out$iterations
        x          <- optim.out$par
        if(optim.out$convergence == 0) {
            converged <- TRUE
        } else {
            converged <- FALSE
        }
    } else if(OPTIMIZER == "NONE") {
        x <- start.x
        iterations <- 0L
        converged <- TRUE
        control <- list()
        optim.out <- list()
    }

    fx <- objective_function(x) # to get "fx.group" attribute

    # transform back
    # 3.
    #if(lavoptions$optim.var.transform == "sqrt" &&
    #       length(lavmodel@x.free.var.idx) > 0L) {
    #    #x[lavmodel@x.free.var.idx] <- tan(x[lavmodel@x.free.var.idx])
    #    x.var <- x[lavmodel@x.free.var.idx]
    #    x.var.sign <- sign(x.var)
    #    x[lavmodel@x.free.var.idx] <- x.var.sign * (x.var * x.var) # square!
    #}

    # 2. unpack
    if(lavmodel@eq.constraints) {
        x <- as.numeric(lavmodel@eq.constraints.K %*% x) +
                        lavmodel@eq.constraints.k0
    }

    # 1. unscale
    x <- x / parscale

    attr(x, "converged")  <- converged
    attr(x, "iterations") <- iterations
    attr(x, "control")    <- control
    attr(x, "fx")         <- fx
    attr(x, "parscale")   <- parscale
    if(!is.null(optim.out$con.jac)) attr(x, "con.jac")    <- optim.out$con.jac
    if(!is.null(optim.out$lambda))  attr(x, "con.lambda") <- optim.out$lambda
    if(lavoptions$partrace) {
        attr(x, "partrace") <- PENV$PARTRACE
    }

    x
}

# backwards compatibility
# estimateModel <- lav_model_estimate

