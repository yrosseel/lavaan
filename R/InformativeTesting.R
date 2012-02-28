# This code is contributed by Leonard Vanbrabant <L.G.F.Vanbrabant@uu.nl>
InformativeTesting <- function(model = NULL, data, constraints = NULL, R = 1000L, 
                               type = "bollen.stine", return.LRT = TRUE, 
                               calibrate = FALSE, calibrate.R = 500L, 
                               calibrate.alpha = 0.05, 
                               parallel = c("no", "multicore", "snow"), 
                               ncpus = 1L, cl = NULL, verbose = FALSE, 
                               stoptest = NULL, conclusion = FALSE, ...){
  
  #checks
  if(!is.null(stoptest)){
        stopifnot(stoptest >=0 & stoptest <=1)
  }
  
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
  if(!is.null(stoptest)){      
    if(lrt.bootA[1] > stoptest){
      cat("p.value A = ", lrt.bootA, "\n")
      stop("Stoptest: Plug-in p-value for test Type A is larger 
                       than the pre-specified value. Further 
           computation is not necessary.")
      }
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
      "  NB. Results for Type A and Type B are inconclusive since no calibration is 
  performed.
    
    In order to obtain meaningfull results set argument 'calibrate=TRUE'. "
  
    conclusion1 <- 
      "  Hypothesis test Type A is significant. 
  Hypothesis test Type B is significant.
    
    The null-hypothesis is rejected in favor of the order constrained hypothesis. 
    
    The order constrained hypothesis is rejected in favor of the unconstrained 
    hypothesis."
 
    conclusion2 <- 
      "  Hypothesis test Type A is significant. 
  Hypothesis test Type B is non-significant. 
    
    The null-hypothesis is rejected in favor of the order constrained hypothesis.
    
    The order constrained hypothesis cannot be rejected in favor of the 
    unconstrained hypothesis."

    conclusion3 <- 
      "  Hypothesis test Type A is non-significant.   
  Hypothesis test Type B is significant.
    
    The null-hypothesis cannot be rejected in favor of the order constrained 
    hypothesis."

    conclusion4 <-
      "  Hypothesis test Type A is non-significant.   
  Hypothesis test Type B is non-significant. 
    
    The null-hypothesis cannot be rejected in favor of the order constrained 
    hypothesis."

    conclusion5 <-
      "  No conclusion can be drawn."
    
    if(calibrate == FALSE){
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
                 calibrate.alpha = calibrate.alpha, calibrate.R = calibrate.R, 
                 calibrate = calibrate, return.LRT = return.LRT, 
                 type = type, conclusion = conclusion, R = R)
  
  if(conclusion){ 
    output$concl = concl
  }        
  
  class(output) <- "InformativeTesting"
  
  return(output)
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
      cat("\n")
      cat("  Conclusion:\n\n")
      t0.txt <- sprintf(object$concl, "\n\n")
      cat(t0.txt, "\n\n", sep="")  
    }  
    #   cat("  For more details see summary --> summary(object)")
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
    cat("\n")
    # cat("  For more details see summary --> summary(object)")
  }
}



summary.InformativeTesting <- function(object, ...){
  cat("\n")
  cat("  Variable names in model       :", unlist(object$fitA1@Data@ov.names), "\n")  
  cat("  Number of variables           :", object$fitA1@Model@nvar, "\n")  
  cat("  Number of groups              :", object$fitA1@Data@ngroups, "\n")  
  cat("  Used sample size per group    :", unlist(object$fitA1@Sample@nobs), "\n")
  cat("  Used sample size              :", sum(unlist(object$fitA1@Data@nobs)), "\n")
  cat("  Total sample size             :", sum(unlist(object$fitA1@Data@norig)), "\n\n")
  cat("  Estimator                     :", object$fitA1@Options$estimator, "\n")
  cat("  Missing data                  :", object$fitA1@Options$missing, "\n\n")
  cat("  Bootstrap method              :", object$type, "\n")
  cat("  Calibrated alpha              :", object$calibrate, "\n")
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
  cat("\n\n\n")  
  
  print(object)
}


plot.InformativeTesting <- function(x, ..., 
                                    type = c("all", 
                                             "LRT.A","LRT.B", 
                                             "plugin.pvalues.A", 
                                             "plugin.pvalues.B"),
                                    main = "main title(s)",
                                    xlab = "xlabel",
                                    ylab = "ylabel",
                                    freq = TRUE,
                                    cex.main = 1,
                                    cex.lab = NULL,
                                    nclass = NULL,
                                    col = "grey",
                                    border = par("fg"),
                                    axes = TRUE,
                                    vline.LRT = FALSE, 
                                    vline.col = "red", 
                                    lty = NULL,
                                    lwd = NULL,
                                    legend = FALSE,
                                    cex.legend = 0.75,
                                    loc.legend = "topright"
                                    ){

  object <- x
  
  #checks
  return.LRT <- object$return.LRT
  calibrate <- object$calibrate
  
  stopifnot(return.LRT | calibrate)
  stopifnot(type %in% c("all", "LRT.A", "LRT.B", 
                        "plugin.pvalues.A", 
                        "plugin.pvalues.B"))
  
  if(type == "plugin.pvalues.A" | type == "plugin.pvalues.B"){
    stopifnot(calibrate)
  }
  if(type == "LRT.A" | type == "LRT.B"){
    stopifnot(return.LRT)
  }
  
  if(missing(type))
    type <- "all"
  
  #multiple images
  if(type == "all" & (return.LRT == FALSE | calibrate == FALSE)){
    par(mfrow = c(1, 2))
    
  }else if(type == "all" & return.LRT & calibrate){
    par(mfrow = c(2, 2))
    
  }else if(length(type) == 1){
    par(mfrow = c(1, 1))
  }
  
  
  #prepare for loop
  if(return.LRT | calibrate){
    if(type == "all"){
      a <- 1L
      b <- 2L
    }else if(type == "LRT.A" | type == "plugin.pvalues.A"){
      a <- 1L
      b <- 1L        
    }else if(type == "LRT.B" | type == "plugin.pvalues.B"){
      a <- 2L
      b <- 2L
    }
  }
  
  
  if(type == "LRT.A" | type == "LRT.B"){
    calibrate <- FALSE
  }else if(type == "plugin.pvalues.A" | type == "plugin.pvalues.B"){
    return.LRT <- FALSE  
  }
  
  #prepare  
  if(return.LRT){
    lrt     <- matrix(as.numeric(NA), length(attr(object$lrt.bootA, 
                                                  "LRT")), 2)
    lrt.obs <- matrix(as.numeric(NA), 1, 2) 
    p.obs   <- matrix(as.numeric(NA), 1, 2)
    
    lrt[, 1]     <- attr(object$lrt.bootA, "LRT")
    lrt.obs[, 1]  <- attr(object$lrt.bootA, "LRT.original")
    p.obs[, 1]     <- object$lrt.bootA[1]
    
    lrt[, 2]     <- attr(object$lrt.bootB, "LRT")
    lrt.obs[, 2]  <- attr(object$lrt.bootB, "LRT.original") 
    p.obs[, 2]   <- object$lrt.bootB[1]
    
  }
  
  if(calibrate){
    pp <- matrix(as.numeric(NA), length(attr(object$lrt.bootA, 
                                             "plugin.pvalues")), 2)
    
    pp[, 1] <- attr(object$lrt.bootA, "plugin.pvalues")
    pp[, 2] <- attr(object$lrt.bootB, "plugin.pvalues")
    
  } 
  
  #prepare titles
  
  if(main == "main title(s)" & type == "all"){
    main <- c("Bootrapped LRT values - type A", 
              "Bootrapped LRT values - type B",
              "Bootrapped plug-in p-values - type A", 
              "Bootrapped plug-in p-values - type B")
    
    title <- as.matrix(main)
    
    
  }else if(main == "main title(s)" & type != "all"){
    main <- "NA"
    title <- main
  }else if(main != "main title(s)" & type == "all"){
    title <- matrix("NA", 4, 1)
    main <- as.matrix(main)
    
    for(k in 1:length(main)){
      title[k, ] <- as.matrix(main[k, ])
    }
  }else if(main != "main title(s)" & type != "all"){
    title <- as.matrix(main)
  }
  
  #make titles
  if(return.LRT){
    
    for(i in a:b){
      main <- title
      
      if(xlab == "xlabel")
        xlab <- "LRT values"
      if(ylab == "ylabel")
        ylab <- "frequency"
      
      
      if(type == "all" & return.LRT & calibrate){
        main <- main[i, ] 
      }else if(type == "all" & return.LRT & calibrate == FALSE){
        main <- main[i, ]
      }else if(type == "LRT.A" | type == "LRT.B"){
        main <- main          
      }
      
      #make histogram
      plot <- hist(lrt[, i], nclass = nclass, xlim = range(lrt[, i]), 
                   plot = FALSE) 
      
      plot(plot, freq = freq, main = main, xlab = xlab, ylab = ylab, 
           cex.main = cex.main, cex.lab = cex.lab, nclass = nclass, 
           col = col, border = border, axes = axes)
      
      #Add vertical line obs. LRT value
      if(vline.LRT){
        abline(v = lrt.obs[, i], col = vline.col, lty = lty, lwd = lwd)
      }        
      #Legend attributes
      if(legend){
        
        plugin.pvalue  <- sprintf("  %15s  %1.3f", "Plug-in p-value:", 
                                  p.obs[i]) 
        obs.lrt        <- sprintf("  %15s  %4.2f", "Observed LRT:", 
                                  lrt.obs[, i])
        legend.obj <- c(obs.lrt, plugin.pvalue)
        
        legend(loc.legend, legend.obj, lty = lty, 
               lwd = lwd, cex = cex.legend, bty="n")
      }
      
    }
  }
  
  if(calibrate){ 
    
    for(i in a:b){
      main <- title
      
      if(xlab == "xlabel")
        xlab <- "Plug-in p-values"
      if(ylab == "ylabel")
        ylab <- "frequency"
      
      
      if(type == "all" & return.LRT & calibrate){
        main <- main[i + 2, ] 
      }else if(type =="all" & main != "main title(s)" & return.LRT == FALSE){
        main <- main[i, ] 
      }else if(type == "all" & return.LRT == FALSE){
        main <- main[i + 2, ]
      }else if(type == "plugin.pvalues.A" | type == "plugin.pvalues.B"){
        main <- main
      }
      
      plot <- hist(pp[, i], nclass = nclass, xlim = range(0:1), 
                   plot = FALSE) 
      
      plot(plot, freq = TRUE, main = main, xlab = xlab, ylab = ylab, 
           cex.main = cex.main, cex.lab = cex.lab, nclass = nclass, 
           col = col, border = border, axes = axes)
    }
    
  }
  
  
}




