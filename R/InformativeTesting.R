InformativeTesting <- function(model = NULL, data, constraints = NULL, 
                               R = 1000L, 
                               type = "bollen.stine", return.LRT = TRUE, 
                               calibrate = FALSE, calibrate.R = 500L, 
                               calibrate.alpha = 0.05, 
                               parallel = c("no", "multicore", "snow"), 
                               ncpus = 1L, cl = NULL,
                               verbose = FALSE, 
                               stoptest = 1, 
                               conclusion = FALSE, 
                               ...){
    
    #checks
    stopifnot(stoptest >= 0 & stoptest <= 1)
    
    #fit unconstrained (free) model and H1
    fit.free <- sem(model, data = data, se = "none", test = "standard", ...) 
    fit.ineq <- sem(model, data = data, se = "none", test = "standard", 
                    constraints = constraints, ...) 

    #generate H0 and H2 models
    con.idx  <- (max(fit.free@User$id) + 1L):max(fit.ineq@User$id)
 
    user.ineq  <- fit.ineq@User
    user.equal <- user.ineq
    user.equal$op[con.idx] <- "=="

    #fit H0
    fit.equal <- sem(user.equal, data = data, se="none", 
                     test = "standard", ...)
 
    #prepare
    fitA1 <- fit.equal
    fitA2 <- fit.ineq
    fitB1 <- fit.ineq
    fitB2 <- fit.free
 
    #Bootstrap LRT values voor Type A (h0 vs. h1)
    lrt.bootA <- lavaan:::bootstrapLRT(fitA1, fitA2, 
                                       R = R, 
                                       type = type, 
                                       verbose = verbose, 
                                       return.LRT = return.LRT, 
                                       calibrate = calibrate, 
                                       calibrate.R = calibrate.R, 
                                       calibrate.alpha = calibrate.alpha,
                                       parallel = parallel, 
                                       ncpus = ncpus, 
                                       cl = cl)

    #check whether computation needs to be continued
    if(lrt.bootA > stoptest){
     cat("p.value A = ", lrt.bootA, "\n")
     stop("Stoptest: Plug-in p-value for test Type A is larger 
                     than the pre-specified value. Further 
                     computation is not necessary.")
    }
 
    #Bootstrap LRT values voor Type B (h1 vs. h2)
    lrt.bootB <- lavaan:::bootstrapLRT(fitB1, fitB2, 
                                       R = R, 
                                       type = type, 
                                       verbose = verbose, 
                                       return.LRT = return.LRT, 
                                       calibrate = calibrate, 
                                       calibrate.R = calibrate.R, 
                                       calibrate.alpha = calibrate.alpha,
                                       parallel = parallel, 
                                       ncpus = ncpus, 
                                       cl = cl)

    #Conclusions
    if(conclusion){
        
      conclusion0 <-
"  NB. Results for Type A and Type B are inconclusive since 
  no calibration is performed.

  In order to obtain meaningfull results set argument 'calibrate=TRUE'. "
  
      conclusion1 <- 
"  Hypothesis test Type A is significant. 
  Hypothesis test Type B is significant.

  The null-hypothesis is rejected in favor of the order 
  constrained hypothesis. 

  The order constrained hypothesis is rejected in favor 
  of the unconstrained hypothesis."
 
      conclusion2 <- 
"  Hypothesis test Type A is significant. 
  Hypothesis test Type B is non-significant. 

  The null-hypothesis is rejected in favor of the order 
  constrained hypothesis.

  The order constrained hypothesis cannot be rejected 
  in favor of the unconstrained hypothesis."

      conclusion3 <- 
"  Hypothesis test Type A is non-significant.   
  Hypothesis test Type B is significant.

  The null-hypothesis cannot be rejected in favor of the 
  order constrained hypothesis."

      conclusion4 <-
"  Hypothesis test Type A is non-significant.   
  Hypothesis test Type B is non-significant. 

  The null-hypothesis cannot be rejected in favor of the 
  order constrained hypothesis."

      conclusion5 <-
"  No conclusion can be drawn."

          if(calibrate==FALSE){
            concl <- conclusion0
          }
          else if(calibrate & 
            lrt.bootA <= attr(lrt.bootA,"adj.alpha") & 
            lrt.bootB <= attr(lrt.bootB,"adj.alpha")){
              concl <- conclusion1
          }
          else if(calibrate & 
            lrt.bootA <= attr(lrt.bootA,"adj.alpha") & 
            lrt.bootB > attr(lrt.bootB,"adj.alpha")){
              concl <- conclusion2
          }
          else if(calibrate & 
            lrt.bootA > attr(lrt.bootA,"adj.alpha") & 
            lrt.bootB <= attr(lrt.bootB,"adj.alpha")){
              concl <- conclusion3
          }
          else if(calibrate & 
            lrt.bootA > attr(lrt.bootA,"adj.alpha") & 
            lrt.bootB > attr(lrt.bootB,"adj.alpha")){
              concl <- conclusion4
          }
          else if(calibrate){
            concl <- conclusion5
          }
    }
    
    output <- list(fitA1 = fitA1, fitA2 = fitA2, fitB1 = fitB1, fitB2 = fitB2,
                   lrt.bootA = lrt.bootA, lrt.bootB = lrt.bootB, 
                   calibrate.alpha = calibrate.alpha, 
                   conclusion = conclusion, 
                   calibrate.R = calibrate.R,
                   calibrate = calibrate, 
                   type = type, 
                   R = R)
      
    if(conclusion){ 
      output$concl = concl
    }        
 
    class(output) <- "InformativeTesting"
          
      return(output)
}

figures <- function(lrt.A = 1L,
                    lrt.B = 1L,
                    plugin.pvalues.A = 1L,
                    plugin.pvalues.B = 1L,
                    lrt.obs.A = 1L,
                    lrt.obs.B = 1L,
                    pvalue.A = 1L,
                    pvalue.B = 1L,
                    calibrate = calibrate,
                    return.LRT = return.LRT){

       if(calibrate & return.LRT){
         par(mfrow = c(2,2))
       }else{
         par(mfrow = c(1,2))
       }

       #LRT distribution Type A  
       p       <- sprintf("  %15s  %1.3f", "Plug-in p-value:", pvalue.A)
       obs.lrt <- sprintf("  %8s %4.2f", "Observed LRT:", lrt.obs.A)
       obj     <- c(obs.lrt, p)
       lrt     <- hist(lrt.A, nclass=17, plot=FALSE)
          plot(lrt, freq = TRUE, col = "grey", border = par("fg"),
               main = paste("LRT-values Distribution - Type A"),
               xlab = "LRT-values", ylab="Freqency", xlim = range(lrt.A), axes = TRUE)
            abline(v = lrt.obs.A, col = "red", lty = 1, lwd = 2)
            legend("topright", obj, col = c("red", "white"), lty = 1, lwd = 2, cex = 0.75, bty="n")



       #LRT distribution Type B
       p       <- sprintf("  %15s  %1.3f", "Plug-in p-value:", pvalue.B)
       obs.lrt <- sprintf("  %8s %4.2f", "Observed LRT:", lrt.obs.B)
       obj     <- c(obs.lrt, p)
       lrt     <- hist(lrt.B, nclass=17, plot=FALSE)
          plot(lrt, freq = TRUE, col = "grey", border = par("fg"),
               main = paste("LRT-values Distribution - Type B"),
               xlab = "LRT-values", ylab="Freqency", xlim = range(lrt.B), axes = TRUE)
             abline(v = lrt.obs.B, col = "red", lty=1, lwd = 2)
             legend("topright", obj, col = c("red", "white"), lty = 1, lwd = 2, cex = 0.75, bty="n")

       if(calibrate){

           #plug-in p-values distribution Type A
           pvalues.A <- hist(plugin.pvalues.A, nclass=17, plot=FALSE)
             plot(pvalues.A, freq = TRUE, col = "grey", border = par("fg"),
                  main = paste("Plug-in p-values Distribution - Type A"),
                  xlab = "Plug-in p-values", ylab="Freqency", xlim = range(plugin.pvalues.A),
                  axes = TRUE, labels = FALSE, add = FALSE, ann = TRUE)

           #plug-in p-values distribution Type B
           pvalues.B <- hist(plugin.pvalues.B, nclass=17, plot=FALSE)
             plot(pvalues.B, freq = TRUE, col = "grey", border = par("fg"),
                  main = paste("Plug-in p-values Distribution - Type B"),
                  xlab = "Plug-in p-values", ylab="Freqency", xlim = range(plugin.pvalues.B),
                  axes = TRUE, labels = FALSE, add = FALSE, ann = TRUE)
       }
}


print.InformativeTesting <- function(x, ...) {
 object <- x
 cat("\n")
 cat("  Order Constrained Hypothesis Testing:\n\n\n")
  if(object$calibrate){
        h1.txt <- sprintf(" %-8s  %12s  %16s  %15s", "", "Adjusted",
                          object$type, "Significant/")
        h2.txt <- sprintf(" %-8s  %12s  %16s  %15s", "", "alpha levels",
                          "Plug-in p-values", "Non-significant")
        cat(h1.txt, "\n")
        cat(h2.txt, "\n")
        cat("-----------------------------------------------------------\n")
        if(object$lrt.bootA <= attr(object$lrt.bootA,"adj.alpha")){
         sig <- "Significant"
        }else{
          sig <- "Non-significant"
        }
        t0.txt <- sprintf("  Type A:")
        t1.txt <- sprintf("  %6s %1.3f", "", attr(object$lrt.bootA,"adj.alpha"))
        t2.txt <- sprintf("  %10s %5.3f", "", object$lrt.bootA[1])
        t3.txt <- sprintf("  %15s", sig)
        cat(t0.txt, t1.txt, t2.txt, t3.txt, "\n", sep="")
        if(object$lrt.bootB <= attr(object$lrt.bootB,"adj.alpha")){
          sig <- "Significant"
        }else{
          sig <- "Non-significant"
         }
        t0.txt <- sprintf("  Type B:")
        t1.txt <- sprintf("  %6s %1.3f", "", attr(object$lrt.bootB,"adj.alpha"))
        t2.txt <- sprintf("  %10s %5.3f", "", object$lrt.bootB[1])
        t3.txt <- sprintf("  %15s", sig)
        cat(t0.txt, t1.txt, t2.txt, t3.txt, "\n\n", sep="")

        if(object$conclusion){
          cat("  Conclusion:\n\n")
          t0.txt <- sprintf(object$concl, "\n\n")
          cat(t0.txt, "\n\n", sep="")
        }
        # cat("  For more details see summary --> summary(object)")
  }else{
        h1.txt <- sprintf(" %-8s  %12s  %16s  %15s", "", "",
                      object$type, "Significant/")
        h2.txt <- sprintf(" %-8s  %12s  %16s  %15s", "", "alpha levels",
                          "Plug-in p-values", "Non-significant")
        cat(h1.txt, "\n")
        cat(h2.txt, "\n")
        cat("---------------------------------------------------------\n")
        if(object$lrt.bootA <= object$calibrate.alpha){
          sig <- "Significant"
        }else{
          sig <- "Non-significant"
         }
        t0.txt <- sprintf("  Type A:")
        t1.txt <- sprintf("  %6s %1.3f", "", object$calibrate.alpha)
        t2.txt <- sprintf("  %10s %5.3f", "", object$lrt.bootA[1])
        t3.txt <- sprintf("  %15s", sig)
        cat(t0.txt, t1.txt, t2.txt, t3.txt, "\n", sep="")
        if(object$lrt.bootB <= object$calibrate.alpha){
          sig <- "Significant"
        }else{
          sig <- "Non-significant"
         }
        t0.txt <- sprintf("  Type B:")
        t1.txt <- sprintf("  %6s %1.3f", "", object$calibrate.alpha)
        t2.txt <- sprintf("  %10s %5.3f", "", object$lrt.bootB[1])
        t3.txt <- sprintf("  %15s", sig)
        cat(t0.txt, t1.txt, t2.txt, t3.txt, "\n\n\n", sep="")
        if(object$conclusion & object$calibrate==FALSE){
          cat("  Conclusion:\n\n")
          t0.txt <- sprintf(object$concl, "\n\n")
          cat(t0.txt, "\n\n", sep="")
        }
        #cat("  For more details see summary --> summary(object)")
  }
  #cat("\n\n")
}

summary.InformativeTesting <- function(object, ...){
     cat("\n")
     cat("  Variable names in model       :", unlist(object$fitA1@Sample@ov.names), "\n")
     cat("  Number of variables           :", object$fitA1@Model@nvar, "\n")
     cat("  Number of groups              :", object$fitA1@Sample@ngroups, "\n")
     cat("  Used sample size per group    :", unlist(object$fitA1@Sample@nobs), "\n")
     cat("  Used sample size              :", sum(unlist(object$fitA1@Sample@nobs)), "\n")
     cat("  Total sample size             :", sum(unlist(object$fitA1@Sample@norig)), "\n\n")
     cat("  Estimator                     :", object$fitA1@Options$estimator, "\n")
     cat("  Missing data                  :", object$fitA1@Options$missing, "\n\n")
     cat("  Bootstrap method              :", object$type, "\n")
     if(object$calibrate){
     cat("  Calibrated alpha Type A       :", attr(object$lrt.bootA,"adj.alpha"), "\n")
     cat("  Calibrated alpha Type B       :", attr(object$lrt.bootB,"adj.alpha"), "\n\n\n")
     }else{
       cat("  ","\n\n\n")
     }
     cat("  Summary models observed data:\n\n")
     h.txt <- sprintf("  %-6s  %-9s  %10s  %2s  %6s", "Type", "Model ",
                      "Chi-square", "df", "p-value")
     cat(h.txt, "\n")
     cat("--------------------------------------------\n")
     t0.txt <- sprintf("  %-6s","Type A")
     t1.txt <- sprintf("  %-9s", "H0 (==)")
     t2.txt <- sprintf("  %10.3f",fitMeasures(object$fitA1, "chisq"))
     t3.txt <- sprintf("  %2.0f",fitMeasures(object$fitA1, "df"))
     t4.txt <- sprintf("    %1.3f",fitMeasures(object$fitA1, "pvalue"))
     cat(t0.txt, t1.txt, t2.txt, t3.txt, t4.txt, "\n", sep="")
     t0.txt <- sprintf("  %6s", "")
     t1.txt <- sprintf("  H1 (<,>)")
     t2.txt <- sprintf("  %11.3f",fitMeasures(object$fitA2, "chisq"))
     t3.txt <- sprintf("  %2.0f",fitMeasures(object$fitA2, "df"))
     t4.txt <- sprintf("    %1.3f",fitMeasures(object$fitA2, "pvalue"))
     cat(t0.txt, t1.txt, t2.txt, t3.txt, t4.txt, "\n\n", sep="")
     t0.txt <- sprintf("  %-6s","Type B")
     t1.txt <- sprintf("  %-9s","H1 (<,>)")
     t2.txt <- sprintf("  %10.3f",fitMeasures(object$fitB1, "chisq"))
     t3.txt <- sprintf("  %2.0f",fitMeasures(object$fitB1, "df"))
     t4.txt <- sprintf("    %1.3f",fitMeasures(object$fitB1, "pvalue"))
     cat(t0.txt, t1.txt, t2.txt, t3.txt, t4.txt, "\n", sep="")
     t0.txt <- sprintf("  %6s", "")
     t1.txt <- sprintf("  %-8s","Hu")
     t2.txt <- sprintf("  %11.3f",fitMeasures(object$fitB2, "chisq"))
     t3.txt <- sprintf("  %2.0f",fitMeasures(object$fitB2, "df"))
     t4.txt <- sprintf("    %1.3f",fitMeasures(object$fitB2, "pvalue"))
     cat(t0.txt, t1.txt, t2.txt, t3.txt, t4.txt, "\n", sep="")
     cat("  (u=unconstrained)")

     cat("\n\n")
     print(object)
     
}





   



