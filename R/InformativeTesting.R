# This code is contributed by Leonard Vanbrabant <L.G.F.Vanbrabant@uu.nl>
InformativeTesting <- function(model = NULL, data, constraints = NULL, 
                               R = NULL, type = "bollen.stine",
                               return.LRT = TRUE, 
                               double.bootstrap = "FDB",
                               double.bootstrap.R = 500L, 
                               double.bootstrap.alpha = 0.05,
                               parallel = c("no", "multicore", "snow"), 
                               ncpus = 1L, cl = NULL, verbose = FALSE, ...){
  #prepare
  if (missing(R) & double.bootstrap == "standard") R <- 1000L
  if (missing(R) & double.bootstrap == "FDB") R <- 2000L
  
  #fit unconstrained (free) model and H1
  fit.free <- sem(model, data = data, se = "none", test = "standard", ...) 
  fit.ineq <- sem(model, data = data, se = "none", test = "standard", 
                  constraints = constraints, ...) 
  
  #generate H0 and H2 models
  con.idx  <- (max(fit.free@User$id) + 1L):max(fit.ineq@User$id)
  
  user.equal  <- fit.ineq@User
  user.equal$op[con.idx] <- "=="
  
  #fit H0
  fit.equal <- sem(user.equal, data = data, se="none", 
                   test = "standard", ...)
  
  
  #Bootstrap LRT values voor Type A (h0 vs. h1)
  lrt.bootA <- bootstrapLRT(fit.equal, fit.ineq, 
                            R = R, type = type, verbose = verbose,
                            return.LRT = return.LRT, 
                            double.bootstrap = double.bootstrap,
                            double.bootstrap.R = double.bootstrap.R, 
                            double.bootstrap.alpha = double.bootstrap.alpha,
                            parallel = parallel, ncpus = ncpus, cl = cl)
  
  
  #Bootstrap LRT values voor Type B (h1 vs. h2)
  lrt.bootB <- bootstrapLRT(fit.ineq, fit.free, 
                            R = R, type = type, verbose = verbose, 
                            return.LRT = return.LRT, 
                            double.bootstrap = double.bootstrap,
                            double.bootstrap.R = double.bootstrap.R, 
                            double.bootstrap.alpha = double.bootstrap.alpha,
                            parallel = parallel, ncpus = ncpus, cl = cl)
  
  
  output <- list(fitA1 = fit.equal, fitA2 = fit.ineq, fitB2 = fit.free,
                 lrt.bootA = lrt.bootA, lrt.bootB = lrt.bootB,
                 double.bootstrap = double.bootstrap, 
                 return.LRT = return.LRT, 
                 type = type)
  
  class(output) <- "InformativeTesting"
  
  return(output)
}



print.InformativeTesting <- function(x, ...) {
  object <- x
  cat("\n")
  cat("Order Constrained Hypothesis Testing:\n\n")
  
  p.A <- object$lrt.bootA[1]
  p.B <- object$lrt.bootB[1]
  sig.A <- sig.B <- "Non-significant"
  
  if (object$double.bootstrap == "FDB") {
    p.A <- attr(object$lrt.bootA, "adj.pvalue")
    p.B <- attr(object$lrt.bootB, "adj.pvalue")
  }
      
  h2.txt <- sprintf(" %-8s  %12s  %16s  %15s", "", "alpha levels", 
                    "Plug-in p-values", "Non-significant")
  
  if (object$double.bootstrap == "standard") {
    h1.txt <- sprintf(" %-8s  %12s  %16s  %15s", "", "Adjusted", 
                      object$type, "Significant/")
    alpha.A <- attr(object$lrt.bootA, "adj.alpha") ##Adjusted p-value????
    alpha.B <- attr(object$lrt.bootB, "adj.alpha")
    if (object$lrt.bootA[1] <= attr(object$lrt.bootA, "adj.alpha")) { 
      sig.A <- "Significant"
    }
    if (object$lrt.bootB[1] <= attr(object$lrt.bootB, "adj.alpha")) { 
      sig.B <- "Significant" 
    }
    cat(h1.txt, "\n")
    cat(h2.txt, "\n")
  }
  else if (object$double.bootstrap == "FDB") {
    h0.txt <- sprintf(" %-8s  %12s  %16s  %15s", "", "", 
                      "FDB", "")
    h1.txt <- sprintf(" %-8s  %12s  %16s  %15s", "", "", 
                      object$type, "Significant/")
    p.A <- attr(object$lrt.bootA, "adj.pvalue")
    p.B <- attr(object$lrt.bootB, "adj.pvalue")
    alpha.A <- alpha.B <- 0.05
    if (attr(object$lrt.bootA, "adj.pvalue") <= 0.05) { 
      sig.A <- "Significant"
    }
    if (attr(object$lrt.bootB, "adj.pvalue") <= 0.05) { 
      sig.B <- "Significant"
    }
    cat(h0.txt, "\n")
    cat(h1.txt, "\n")
    cat(h2.txt, "\n")
  }
  else if (object$double.bootstrap == "no") {
    h1.txt <- sprintf(" %-8s  %12s  %16s  %15s", "", "", 
                      object$type, "Significant/")
    alpha.A <- alpha.B <- 0.05
    sig.A <- sig.B <- "Inconclusive*"
    cat(h1.txt, "\n")
    cat(h2.txt, "\n")
  }  
    cat("-----------------------------------------------------------\n") 
    t0.txt <- sprintf("  Type A:")
    t1.txt <- sprintf("  %6s %1.3f", "", alpha.A) 
    t2.txt <- sprintf("  %10s %5.3f", "", p.A)
    t3.txt <- sprintf("  %15s", sig.A)
    cat(t0.txt, t1.txt, t2.txt, t3.txt, "\n", sep="")
    t0.txt <- sprintf("  Type B:")
    t1.txt <- sprintf("  %6s %1.3f", "", alpha.B) 
    t2.txt <- sprintf("  %10s %5.3f", "", p.B)
    t3.txt <- sprintf("  %15s", sig.B)
    cat(t0.txt, t1.txt, t2.txt, t3.txt, "\n\n", sep="")
    if (sig.A == "Inconclusive*"){
      cat("  *For meaningfull results set double.bootstrap")
    }
}


summary.InformativeTesting <- function(object, 
                                       fit = c("h0","hi","hu"), ...)
{
  stopifnot (fit %in% c("h0","hi","hu"))
  if (!missing(fit)){
    if (fit == "h0") {
      summary(object$fitA1)
    }
    else if (fit == "hi") {
      summary(object$fitA2)
    }
    else if (fit == "hu") {
      summary(object$fitB2)
    }
  }
  if (missing(fit)) {
    cat("  \n")
    cat("  Variable names in model       :", unlist(object$fitA1@Data@ov.names[1]), "\n")  
    cat("  Number of variables           :", object$fitA1@Model@nvar[1], "\n")  
    cat("  Number of groups              :", object$fitA1@Data@ngroups, "\n")  
    cat("  Used sample size per group    :", unlist(object$fitA1@Data@nobs), "\n")
    cat("  Used sample size              :", sum(unlist(object$fitA1@Data@nobs)), "\n")
    cat("  Total sample size             :", sum(unlist(object$fitA1@Data@norig)), "\n\n")
    cat("  Estimator                     :", object$fitA1@Options$estimator, "\n")
    cat("  Missing data                  :", object$fitA1@Options$missing, "\n")
    cat("  Bootstrap method              :", object$type, "\n")
    cat("  Double bootstrap method       :", object$double.bootstrap, "\n")
    cat(  "\n\nFitted models:\n\n")
    h.txt <- sprintf("  %-6s  %-9s  %10s  %2s  %6s", "Type", "Model ", 
                     "Chi-square", "df", "p-value")
    cat(h.txt, "\n")
    cat(  "--------------------------------------------\n")
    t0.txt <- sprintf("  %-6s","Type A")
    t1.txt <- sprintf("  %-9s", "H0 (==)")
    t2.txt <- sprintf("  %10.3f",fitMeasures(object$fitA1, "chisq")) 
    t3.txt <- sprintf("  %2.0f",fitMeasures(object$fitA1, "df")) 
    t4.txt <- sprintf("    %1.3f",fitMeasures(object$fitA1, "pvalue"))
    cat(t0.txt, t1.txt, t2.txt, t3.txt, t4.txt, "\n", sep="")
    t0.txt <- sprintf("  %6s", "")
    t1.txt <- sprintf("  Hi (<,>)")
    t2.txt <- sprintf("  %11.3f",fitMeasures(object$fitA2, "chisq")) 
    t3.txt <- sprintf("  %2.0f",fitMeasures(object$fitA2, "df")) 
    t4.txt <- sprintf("    %1.3f",fitMeasures(object$fitA2, "pvalue"))
    cat(t0.txt, t1.txt, t2.txt, t3.txt, t4.txt, "\n\n", sep="")
    t0.txt <- sprintf("  %-6s","Type B")
    t1.txt <- sprintf("  %-9s","Hi (<,>)")
    t2.txt <- sprintf("  %10.3f",fitMeasures(object$fitA2, "chisq")) 
    t3.txt <- sprintf("  %2.0f",fitMeasures(object$fitA2, "df")) 
    t4.txt <- sprintf("    %1.3f",fitMeasures(object$fitA2, "pvalue"))
    cat(t0.txt, t1.txt, t2.txt, t3.txt, t4.txt, "\n", sep="")
    t0.txt <- sprintf("  %6s", "")
    t1.txt <- sprintf("  %-8s","Hu")
    t2.txt <- sprintf("  %11.3f",fitMeasures(object$fitB2, "chisq")) 
    t3.txt <- sprintf("  %2.0f",fitMeasures(object$fitB2, "df")) 
    t4.txt <- sprintf("    %1.3f",fitMeasures(object$fitB2, "pvalue"))
    cat(t0.txt, t1.txt, t2.txt, t3.txt, t4.txt, "\n", sep="")
    cat("\n  (i = informative, u = unconstrained)")
    cat("\n\n\n")  
    
    print(object)
  }
}


plot.InformativeTesting <- function(x, ..., 
                                    type = "all",
                                    main = "main",
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
                                    loc.legend = "topright") 
 {
  
  #prepare
  object <- x
  return.LRT <- object$return.LRT
  double.bootstrap <- object$double.bootstrap
  
  #checks
  stopifnot(type %in% c("all", "LRT.A", "LRT.B", 
                        "ppvalues.A", "ppvalues.B"))
  if (type == "ppvalues.A" | type == "ppvalues.B") {
    stopifnot(double.bootstrap == "standard")
  }
  if (type == "LRT.A" | type == "LRT.B") {
    stopifnot(return.LRT)
  }
  if (type == "all" & !return.LRT) { 
    stopifnot (double.bootstrap != "FDB")
  }
  
  #prepare data
  pvalue <- rep(as.numeric(NA), 2) 
   pvalue[1]   <- object$lrt.bootA[1]
   pvalue[2]   <- object$lrt.bootB[1]
  
  if (ylab == "ylabel") {
    y.lab  <- "Frequency"
  } else {
    y.lab <- ylab
  }
  
  if (return.LRT) {
    lrt.obs <- rep(as.numeric(NA), 2) 
      lrt.obs[1]  <- attr(object$lrt.bootA, "LRT.original")
      lrt.obs[2]  <- attr(object$lrt.bootB, "LRT.original")
        
      lrt.A <- attr(object$lrt.bootA, "LRT")
      lrt.B <- attr(object$lrt.bootB, "LRT")
    
    lrt <- list(data.frame(lrt.A), data.frame(lrt.B))
    lrt <- do.call(rbind.fill, lrt)
            
    if (xlab == "xlabel") { 
      x.lrt <- c("Bootstrapped LRT values")
    }
    else {
      x.lrt <- xlab
    }
    if (main == "main") {
      m.lrt <- c("Distribution of bootstrapped LRT values - type A", 
                 "Distribution of bootstrapped LRT values - type B")
    }
    else {
      m.lrt <- main
    }
  } 
  if (double.bootstrap == "FDB") {
      lrt.q <- rep(as.numeric(NA), 2)
        lrt.q[1] <- attr(object$lrt.bootA, "lrt.q") 
        lrt.q[2] <- attr(object$lrt.bootB, "lrt.q")
      adj.pvalue <- rep(as.numeric(NA), 2)
        adj.pvalue[1] <- attr(object$lrt.bootA, "adj.pvalue")
        adj.pvalue[2] <- attr(object$lrt.bootB, "adj.pvalue")
      
      if (xlab == "xlabel") {
        x.lrt <- c("Bootstrapped LRT values")
      }
      else {
        x.lrt <- xlab
      }
      if (main == "main") {
        m.lrt <- c("Distribution of first level LRT values - type A", 
                   "Distribution of first level LRT values - type B")
      }
      else {
        m.lrt <- main
      }
    }
  
  if (double.bootstrap == "standard") {
    ppvalue.A <- attr(object$lrt.bootA, "plugin.pvalues")
    ppvalue.B <- attr(object$lrt.bootB, "plugin.pvalues")
    
    ppv <- list(data.frame(ppvalue.A), data.frame(ppvalue.B))
    ppv <- do.call(rbind.fill, ppv)
    
    if (xlab == "xlabel") {
      x.ppv  <- c("Bootstrapped plug-in p values")
    }
    else {
      x.ppv <- xlab    
    }
    if (main == "main") {
      m.ppv  <- c("Distribution of plug-in p values - type A", 
                  "Distribution of plug-in p values - type B")
    }else {
      m.ppv <- main
    }  
  }
  
  if (return.LRT & type == "all" & double.bootstrap != "standard") {
    par(mfrow = c(1, 2))
  }
  if (return.LRT & type == "all" & double.bootstrap == "standard") {
    par(mfrow = c(2, 2))
  }
  if (!return.LRT & (double.bootstrap == "standard" | 
       double.bootstrap == "FDB")) {
    par(mfrow = c(1, 2))
  }
  if (type != "all") {
    par(mfrow = c(1, 1))
  }
    
  #make histograms
  if (type == "LRT.A" | type == "LRT.B") {
    double.bootstrap = "no"
  }
  if (type == "ppvalues.A" | type == "ppvalues.B") {
    return.LRT <- FALSE
  }
  if ((type == "LRT.A" & return.LRT) |
      (type == "ppvalues.A" & double.bootstrap == "standard")) { 
    a = 1L
    b = 1L
  }
  if ((type == "LRT.B" & return.LRT) |
      (type == "ppvalues.B" & double.bootstrap == "standard")) {
    a = 2L
    b = 2L
  }
  if (type == "all") {
    a = 1L
    b = 2L
  }  
  
  for (i in a:b) {
    if (return.LRT) {
      plot <- hist(lrt[, i], nclass = nclass, xlim = range(lrt[, i]), plot = FALSE) 
      plot(plot, freq = freq, main = m.lrt[i], xlab = x.lrt, ylab = y.lab, 
           cex.main = cex.main, cex.lab = cex.lab, nclass = nclass, 
           col = col, border = border, axes = axes) 
      
      if (vline.LRT & double.bootstrap != "FDB") {
        abline(v = lrt.obs[i], col = vline.col, lty = lty, lwd = lwd)
      }
      else if (vline.LRT & double.bootstrap == "FDB") {
        abline(v = lrt.q[i], col = vline.col, lty = lty, lwd = lwd)
      }
      
      if (legend & double.bootstrap != "FDB") {
          ppval    <- sprintf("  %15s  %1.3f", "plug-in p value", pvalue[i]) 
          obs.lrt  <- sprintf("  %15s  %4.2f", "obs LRT", lrt.obs[i])
            legend.obj <- c(obs.lrt, ppval)
          legend(loc.legend, legend.obj, lty = lty, 
                 lwd = lwd, cex = cex.legend, bty="n")
      }
      if (legend & double.bootstrap == "FDB") {
        ppval      <- sprintf("  %15s  %1.3f", "FDB p-value", adj.pvalue[i]) 
        obs.lrt.q  <- sprintf("  %15s  %4.2f", "obs LRT.q", lrt.q[i])
        legend.obj <- c(obs.lrt.q, ppval)
        legend(loc.legend, legend.obj, lty = lty, 
               lwd = lwd, cex = cex.legend, bty="n")
      }
    }
    if (double.bootstrap == "standard") {      
        plot <- hist(ppv[, i], nclass = nclass, xlim = range(0:1), 
                     plot = FALSE) 
        plot(plot, freq = TRUE, main = m.ppv[i], xlab = x.ppv, ylab = y.lab, 
             cex.main = cex.main, cex.lab = cex.lab, nclass = nclass, 
             col = col, border = border, axes = axes) 
      }
  }
}




