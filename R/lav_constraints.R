lav_constraints_parse <- function(partable = NULL, constraints = NULL,
                                  theta = NULL,
                                  debug = FALSE) {

    # just in case we do not have a $free column in partable
    if(is.null(partable$free)) {
        partable$free <- rep(1L, length(partable$lhs))
    }

    # from the partable: free parameters
    if(!is.null(theta)) {
        # nothing to do
    } else if(!is.null(partable$est)) {
        theta <- partable$est[ partable$free > 0L ]
    } else if(!is.null(partable$start)) {
        theta <- partable$start[ partable$free > 0L ]
    } else {
        theta <- rep(0, length(partable$lhs))
    }

    # number of free parameters
    npar <- length(theta)

    # parse the constraints
    if(is.null(constraints)) {
        LIST <- NULL
    } else if(!is.character(constraints)) {
        stop("lavaan ERROR: constraints should be a string")
    } else {
        FLAT <- lavParseModelString( constraints )
        CON <- attr(FLAT, "constraints")
        LIST <- list()
        if(length(CON) > 0L) {
            lhs = unlist(lapply(CON, "[[", "lhs"))
             op = unlist(lapply(CON, "[[",  "op"))
            rhs = unlist(lapply(CON, "[[", "rhs"))
            LIST$lhs        <- c(LIST$lhs,        lhs)
            LIST$op         <- c(LIST$op,         op)
            LIST$rhs        <- c(LIST$rhs,        rhs)
        } else {
            stop("lavaan ERROR: no constraints found in constraints argument")
        }
    }

    # variable definitions
    def.function <- lav_partable_constraints_def(partable, con = LIST,
                                                 debug = debug)

    # construct ceq/ciq functions
    ceq.function <- lav_partable_constraints_ceq(partable, con = LIST,
                                                 debug = debug)
    # linear or nonlinear?
    ceq.linear.idx <- lav_constraints_linear_idx(func = ceq.function,
                                                 npar = npar)
    ceq.nonlinear.idx <- lav_constraints_nonlinear_idx(func = ceq.function,
                                                       npar = npar)

    # inequalities
    cin.function <- lav_partable_constraints_ciq(partable, con = LIST,
                                                 debug = debug)

    # linear or nonlinear?
    cin.linear.idx <- lav_constraints_linear_idx(func = cin.function,
                                                 npar = npar)
    cin.nonlinear.idx <- lav_constraints_nonlinear_idx(func = cin.function,
                                                       npar = npar)

    # Jacobians
    if(!is.null(body(ceq.function))) {
        ceq.JAC <- try(lav_func_jacobian_complex(func = ceq.function,
                                                 x = theta), silent=TRUE)
        if(inherits(ceq.JAC, "try-error")) { # eg. pnorm()
            ceq.JAC <- lav_func_jacobian_simple(func = ceq.function, x = theta)
        }

        # constants
        # do we have a non-zero 'rhs' elements? FIXME!!! is this reliable??
        ceq.rhs <- -1 * ceq.function( numeric(npar) )
    } else {
        ceq.JAC <- matrix(0, nrow = 0L, ncol = npar)
        ceq.rhs <- numeric(0L)
    }

    if(!is.null(body(cin.function))) {
        cin.JAC <- try(lav_func_jacobian_complex(func = cin.function,
                                                 x = theta), silent=TRUE)
        if(inherits(cin.JAC, "try-error")) { # eg. pnorm()
            cin.JAC <- lav_func_jacobian_simple(func = cin.function, x = theta)
        }

        # constants
        # do we have a non-zero 'rhs' elements? FIXME!!! is this reliable??
        cin.rhs <- -1 * cin.function( numeric(npar) )
    } else {
        cin.JAC <- matrix(0, nrow = 0L, ncol = npar)
        cin.rhs <- numeric(0L)
    }

    OUT <- list(def.function      = def.function,

                ceq.function      = ceq.function,
                ceq.JAC           = ceq.JAC,
                ceq.rhs           = ceq.rhs,
                ceq.linear.idx    = ceq.linear.idx,
                ceq.nonlinear.idx = ceq.nonlinear.idx,

                cin.function      = cin.function,
                cin.JAC           = cin.JAC,
                cin.rhs           = cin.rhs,
                cin.linear.idx    = cin.linear.idx,
                cin.nonlinear.idx = cin.nonlinear.idx)

    OUT
}

lav_constraints_linear_idx <- function(func = NULL, npar = NULL) { 

    if(is.null(func) || is.null(body(func))) return(integer(0L))
    
    # seed 1: rnorm
    A0 <- lav_func_jacobian_complex(func = func, x = rnorm(npar))

    # seed 2: rnorm
    A1 <- lav_func_jacobian_complex(func = func, x = rnorm(npar))

    A0minA1 <- A0 - A1
    linear <- apply(A0minA1, 1,  function(x) all(x == 0))
    which(linear)
}

lav_constraints_nonlinear_idx <- function(func = NULL, npar = NULL) {   
  
    if(is.null(func) || is.null(body(func))) return(integer(0L))

    # seed 1: rnorm
    A0 <- lav_func_jacobian_complex(func = func, x = rnorm(npar))

    # seed 2: rnorm
    A1 <- lav_func_jacobian_complex(func = func, x = rnorm(npar))

    A0minA1 <- A0 - A1
    linear <- apply(A0minA1, 1,  function(x) all(x == 0))
    which(!linear)
}


# FIXME: is there a more elegant/robust way to do this??
lav_constraints_check_linear <- function(model) {

     # seed 1: rnorm
     A.ceq <- A.cin <- matrix(0, model@nx.free, 0)
     if(!is.null(body(model@ceq.function)))
         A.ceq <- t(lav_func_jacobian_complex(func=model@ceq.function, x=rnorm(model@nx.free)))
     if(!is.null(body(model@cin.function)))
         A.cin <- t(lav_func_jacobian_complex(func=model@cin.function, x=rnorm(model@nx.free)))
     A0 <- cbind(A.ceq, A.cin)    

     # seed 2: rnorm
     A.ceq <- A.cin <- matrix(0, model@nx.free, 0)
     if(!is.null(body(model@ceq.function)))
         A.ceq <- t(lav_func_jacobian_complex(func=model@ceq.function, x=rnorm(model@nx.free)))
     if(!is.null(body(model@cin.function)))
         A.cin <- t(lav_func_jacobian_complex(func=model@cin.function, x=rnorm(model@nx.free)))
     A1 <- cbind(A.ceq, A.cin)

     A0minA1 <- all.equal(A0, A1)
     if(is.logical(A0minA1) && A0minA1 == TRUE)
         return(TRUE)
     else
         return(FALSE)
}

