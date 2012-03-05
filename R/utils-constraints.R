# FIXME: is there a more elegant/robust way to do this??
checkLinearConstraints <- function(model) {

     # seed 1: rnorm
     A.ceq <- A.cin <- matrix(0, model@nx.free, 0)
     if(!is.null(body(model@ceq.function)))
         A.ceq <- t(jacobian(func=model@ceq.function, x=rnorm(model@nx.free)))
     if(!is.null(body(model@cin.function)))
         A.cin <- t(jacobian(func=model@cin.function, x=rnorm(model@nx.free)))
     A0 <- cbind(A.ceq, A.cin)    

     # seed 2: rnorm
     A.ceq <- A.cin <- matrix(0, model@nx.free, 0)
     if(!is.null(body(model@ceq.function)))
         A.ceq <- t(jacobian(func=model@ceq.function, x=rnorm(model@nx.free)))
     if(!is.null(body(model@cin.function)))
         A.cin <- t(jacobian(func=model@cin.function, x=rnorm(model@nx.free)))
     A1 <- cbind(A.ceq, A.cin)

     A0minA1 <- all.equal(A0, A1)
     if(is.logical(A0minA1) && A0minA1 == TRUE)
         return(TRUE)
     else
         return(FALSE)
}
