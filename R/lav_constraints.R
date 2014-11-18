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

