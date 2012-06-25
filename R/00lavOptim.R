# eventually, this file will contain all routines related to
# optimization -- YR 21 june 2012


# super class -- virtual statistical model that needs to be optimized
lavRefOptim <- setRefClass("lavOptim",

# fields
fields = list(
    npar            = "integer",     # number of free model parameters
    theta           = "numeric",     # the model parameters (free only)
    theta.start     = "numeric",     # starting values
    optim.method    = "character",   # optimization method
    optim.control   = "list",        # control parameter for optimization method
    optim.out       = "list"         # optimization results
),

# methods
methods = list(

show = function() {
    cat(class(.self), "model parameters (theta):\n")
    print(theta)
},

objective = function(x) {
    cat("this is a dummy function\n")
    return(as.numeric(NA))
},

gradient = function(x) {
    cat("this is dummy function\n")
    return(rep(as.numeric(NA), npar))
},

hessian = function(x) {
    cat("this is dummy function\n")
    return(matrix(as.numeric(NA), npar, npar))
}, 

start = function() {
    rep(0, npar)    
},

optimize = function(method = "nlminb", control = list(), verbose = FALSE,
                    start = NULL) {
    method <- tolower(method)
    hessian <- FALSE
    if( method == "none" ) {
        optim.method <<- "none"
    } else if( method %in% c("nlminb", "quasi-newton", "quasi.newton", 
                      "nlminb.hessian") ) {
        optim.method <<- "nlminb"
        if(verbose)
            control$trace <- 1L
        if(method == "nlminb.hessian")
            hessian <- TRUE
        control.nlminb <- list(eval.max=20000L,
                               iter.max=10000L,
                               trace=0L,
                               #abs.tol=1e-20, ### important!! fx never negative
                               abs.tol=(.Machine$double.eps * 10),
                               rel.tol=1e-10,
                               x.tol=1.5e-8,
                               step.min=2.2e-14)
        control.nlminb <- modifyList(control.nlminb, control)
        optim.control <<- control.nlminb[c("eval.max", "iter.max", "trace",
                                           "abs.tol", "rel.tol", "x.tol",
                                           "step.min")]
        
    } else if( method %in% c("newton", "newton-raphson", "newton.raphson") ) {
        optim.method <<- "newton"
        if(verbose)
            control$verbose <- TRUE
        control.nr <- list(grad.tol = 1e-6,
                           iter.max = 200L,
                           inner.max = 20L,
                           verbose   = FALSE)
        control.nr <- modifyList(control.nr, control)
        optim.control <<- control.nr[c("grad.tol", "iter.max", "inner.max",
                                       "verbose")]
    } else {
        stop("unknown optim method: ", optim.method)
    }

    # starting values
    if(is.null(start)) {
        theta.start <<- .self$start()
    } else {
        stopifnot(length(start) == npar)
        theta.start <<- start
    }

    # run objective function to intialize (and see if starting values
    # are valid
    tmp <- objective(theta.start)

    if(optim.method == "newton") {
        out <- lavOptimNewtonRaphson(object=.self, control = optim.control)
        optim.out <<- out
    } else if(optim.method == "nlminb") {
        if(!hessian) {
            out <- nlminb(start = theta, objective = .self$objective,
                          gradient = .self$gradient, control = optim.control)
        } else {
            out <- nlminb(start = theta, objective = .self$objective,
                          gradient = .self$gradient, 
                          hessian = .self$hessian,
                          control = optim.control)
        }
        # FIXME: use generic fields
        optim.out <<- out
    }
}

))


# this is a simple/naive Newton Raphson implementation
# - minimization only
# - it assumes that the hessian is always positive definite (no check!)
# - it may do some backstepping, but there is no guarantee that it will
#   converge
# this function is NOT for general-purpose optimization, but should only be 
# used or simple (convex!) problems (eg. estimating polychoric/polyserial 
# correlations, probit regressions, ...)
# 
lavOptimNewtonRaphson <- function(object,
                                  control = list(iter.max  = 100L,
                                                 grad.tol  = 1e-6,
                                                 inner.max = 20L,
                                                 verbose   = FALSE)) {

    # housekeeping
    converged <- FALSE; message <- character(0); inner <- 0L

    # we start with classic newton: step length alpha_k = 1
    alpha_k <- 1

    # current estimates
    fx.old <- object$objective()
    gradient <- object$gradient()
    norm.grad <- sqrt( crossprod(gradient) )
    max.grad <- max(abs(gradient))

    # start NR steps
    for(i in seq_len(control$iter.max)) {

        if(control$verbose) {
            cat("NR step ", sprintf("%2d", (i-1L)), ": max.grad = ", 
                sprintf("%12.9f", max.grad), " norm.grad = ", 
                sprintf("%12.9f", norm.grad), "\n", sep="")
        }

        # simple convergence test
        if(max.grad < control$grad.tol) {
            converged <- TRUE
            message <- paste("converged due to max.grad < tol (",
                             control$grad.tol, ")", sep="")
            if(control$verbose) cat("NR ", message, "\n", sep="")
            break
        }

        # compute search direction 'p_k'
        hessian <- object$hessian()
        p_k <- as.numeric( solve(hessian, gradient) )

        # update theta and fx
        theta <- object$theta - alpha_k * p_k
        fx.new <- object$objective(theta)

        # check if we minimize
        if(fx.new > fx.old) {
            # simple backstepping
            for(j in seq_len(control$inner.max)) {
                inner <- inner + 1L
                alpha_k <- alpha_k/2
                theta <<- object$theta + alpha_k * p_k
                fx.new <- object$objective(theta)
                if(fx.new < fx.old)
                    break
            }
            if(fx.new > fx.old) { # it didn't work... bail out
                message <- paste("backstepping failed after ",
                                 control$inner.max, "iterations")
                return(list(fx=fx.new, converged=FALSE, message=message,
                            max.grad=max.grad,iterations=i, backsteps=inner))
            }
            # reset alpha
            alpha_k <- 1
        }

        # update
        gradient <- object$gradient()
        max.grad <- max(abs(gradient))
        norm.grad <- sqrt( crossprod(gradient) )
        fx.old <- fx.new
    }

    list(fx          = fx.old,
         converged   = converged,
         message     = message,
         max.grad    = max.grad,
         norm.grad   = norm.grad,
         iterations  = i-1L,
         backsteps   = inner)
}



