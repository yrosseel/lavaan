# FIXME: is there a more elegant/robust way to do this??
checkLinearConstraints <- function(model) {

     # seed = zero vector
     A.ceq <- t(jacobian(func=model@ceq.function, x=rep(0,model@nx.free)))
     A.cin <- t(jacobian(func=model@cin.function, x=rep(0,model@nx.free)))
     A0 <- cbind(A.ceq, A.cin)    

     # seed = unit vector
     A.ceq <- t(jacobian(func=model@ceq.function, x=rep(1,model@nx.free)))
     A.cin <- t(jacobian(func=model@cin.function, x=rep(1,model@nx.free)))
     A1 <- cbind(A.ceq, A.cin)

     if(sum(A0 - A1) == 0)
         return(TRUE)
     else
         return(FALSE)
}
