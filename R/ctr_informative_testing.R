# This code is contributed by Leonard Vanbrabant <L.G.F.Vanbrabant@hotmail.com>
InformativeTesting <- function(model = NULL, data, constraints = NULL, 
                               R = 1000L, type = "bollen.stine",
                               return.LRT = TRUE, 
                               double.bootstrap = "standard",
                               double.bootstrap.R = 500L,
                               double.bootstrap.alpha = 0.05,
                               parallel = c("no", "multicore", "snow"), 
                               ncpus = 1L, cl = NULL, verbose = FALSE, ...){
  
  fit.B1 <- sem(model, data = data, se = "none", test = "standard", ...) 
  
  fit.B0 <- fit.A1 <- sem(model, data = data, se = "none", test = "standard", 
                          constraints = constraints, ...) 
  
  con.idx  <- (max(fit.B1@ParTable$id) + 1L):max(fit.A1@ParTable$id)
  
  user.equal  <- fit.A1@ParTable
  user.equal$op[con.idx] <- "=="
  
  fit.A0 <- sem(user.equal, data = data, se="none", test = "standard", ...)
  
  lrt.bootA <- bootstrapLRT(fit.A0, fit.A1, 
                            R = R, type = type, verbose = verbose,
                            return.LRT = return.LRT, 
                            double.bootstrap = double.bootstrap,
                            double.bootstrap.R = double.bootstrap.R, 
                            double.bootstrap.alpha = double.bootstrap.alpha,
                            parallel = parallel, ncpus = ncpus, cl = cl)
  
  lrt.bootB <- bootstrapLRT(fit.B0, fit.B1, 
                            R = R, type = type, verbose = verbose, 
                            return.LRT = return.LRT, 
                            double.bootstrap = double.bootstrap,
                            double.bootstrap.R = double.bootstrap.R, 
                            double.bootstrap.alpha = double.bootstrap.alpha,
                            parallel = parallel, ncpus = ncpus, cl = cl)
  
  output <- list(fit.A0 = fit.A0, fit.A1 = fit.A1, fit.B1 = fit.B1,
                 lrt.bootA = lrt.bootA, lrt.bootB = lrt.bootB,
                 double.bootstrap = double.bootstrap,
                 double.bootstrap.alpha = double.bootstrap.alpha,
                 return.LRT = return.LRT, type = type)
  
  class(output) <- "InformativeTesting"
  
  return(output)
}



print.InformativeTesting <- function(x, ...) {
  object <- x
  
  cat("  \n")
  cat("  Variable names in model       :", unlist(object$fit.A1@Data@ov.names[1]), "\n")  
  cat("  Number of variables           :", object$fit.A1@Model@nvar[1], "\n")  
  cat("  Number of groups              :", object$fit.A1@Data@ngroups, "\n")  
  cat("  Used sample size per group    :", unlist(object$fit.A1@Data@nobs), "\n")
  cat("  Used sample size              :", sum(unlist(object$fit.A1@Data@nobs)), "\n")
  cat("  Total sample size             :", sum(unlist(object$fit.A1@Data@norig)), "\n\n")
  cat("  Estimator                     :", object$fit.A1@Options$estimator, "\n")
  cat("  Missing data                  :", object$fit.A1@Options$missing, "\n")
  cat("  Bootstrap method              :", object$type, "\n")
  cat("  Double bootstrap method       :", object$double.bootstrap, "\n\n\n")
  
  p.A <- object$lrt.bootA[1]
  p.B <- object$lrt.bootB[1]
  sig.A <- sig.B <- "Non-significant"
  a.A <- a.B <- object$double.bootstrap.alpha
  
  if (object$double.bootstrap == "FDB" | object$double.bootstrap == "standard"){
    if (object$double.bootstrap == "FDB") {
      adj.p.A <- attr(object$lrt.bootA, "adj.pvalue")
      adj.p.B <- attr(object$lrt.bootB, "adj.pvalue")
    }
    if (object$double.bootstrap == "standard") {
      adj.p.A <- attr(object$lrt.bootA, "adj.pvalue")
      adj.p.B <- attr(object$lrt.bootB, "adj.pvalue")
      adj.a.A <- attr(object$lrt.bootA, "adj.alpha")
      adj.a.B <- attr(object$lrt.bootB, "adj.alpha")
    }
    if (attr(object$lrt.bootA, "adj.pvalue") < object$double.bootstrap.alpha) 
      sig.A <- "Significant"
    if (attr(object$lrt.bootB, "adj.pvalue") < object$double.bootstrap.alpha) 
      sig.B <- "Significant"
  }
  else {
    sig.A <- sig.B <- "Inconclusive*"
  }
  
  cat("Order Constrained Hypothesis Testing:\n\n")
  h0.txt <- sprintf("  %-22s  %8s  %12s", "", "Type A", "Type B") 
  cat(h0.txt, "\n")
  cat("  -------------------------------------------------\n")
  t00.txt <- sprintf("  %-22s", "LR statistic")
  t01.txt <- sprintf("  %s %7.2f", "", attr(object$lrt.bootA, "LRT.original"))
  t02.txt <- sprintf("  %s %11.2f", "", attr(object$lrt.bootB, "LRT.original"))
  
  if (object$double.bootstrap == "standard"){
    t10.txt <- sprintf("  %-22s", "Adjusted alpha")
    t11.txt <- sprintf("  %s %7.2f", "", adj.a.A)
    t12.txt <- sprintf("  %s %11.2f", "", adj.a.B)
  }
  if (object$double.bootstrap == "FDB" | object$double.bootstrap == "standard"){
    t20.txt <- sprintf("  %-22s", "Adjusted p-value")
    t21.txt <- sprintf("  %s %7.2f", "", adj.p.A)
    t22.txt <- sprintf("  %s %11.2f", "", adj.p.B)
  }
  t30.txt <- sprintf("  %-22s", "P-value")
  t31.txt <- sprintf("  %s %7.2f", "", p.A)
  t32.txt <- sprintf("  %s %11.2f", "", p.B)
  t40.txt <- sprintf("  %-22s", "Alpha")
  t41.txt <- sprintf("  %s %7.2f", "", a.A)
  t42.txt <- sprintf("  %s %11.2f", "", a.B)
  t50.txt <- sprintf("  %10s", "Significance")
  t51.txt <- sprintf("  %18s", sig.A)
  t52.txt <- sprintf("  %15s", sig.B)
  
  cat(t00.txt, t01.txt, t02.txt, "\n\n", sep="")
  
  if (object$double.bootstrap == "FDB"){
    cat(t30.txt, t31.txt, t32.txt, "\n", sep="")
    cat(t20.txt, t21.txt, t22.txt, "\n\n", sep="")
    cat(t40.txt, t41.txt, t42.txt, "\n\n", sep="") 
    cat(t50.txt, t51.txt, t52.txt, "\n", sep="")
  } 
  else if(object$double.bootstrap == "standard"){
    cat(t30.txt, t31.txt, t32.txt, "\n", sep="")
    cat(t10.txt, t11.txt, t12.txt, "\n\n", sep="")
    cat(t20.txt, t21.txt, t22.txt, "\n", sep="")
    cat(t40.txt, t41.txt, t42.txt, "\n\n", sep="") 
    cat(t50.txt, t51.txt, t52.txt, "\n", sep="")
  }
  else if(object$double.bootstrap == "no"){
    cat(t30.txt, t31.txt, t32.txt, "\n", sep="")
    cat(t40.txt, t41.txt, t42.txt, "\n\n", sep="") 
    cat("  *For meaningfull results set double.bootstrap")
  }
  
}


plot.InformativeTesting <- function(x, ..., 
                                    type = "all",
                                    main = "main",
                                    xlab = "xlabel",
                                    ylab = "Frequency",
                                    freq = TRUE,
                                    breaks = 15,
                                    cex.main = 1,
                                    cex.lab = NULL,
                                    cex.axis = NULL,
                                    col = "grey",
                                    border = par("fg"),
                                    axes = TRUE,
                                    vline = TRUE, 
                                    vline.col = c("red", "blue"), 
                                    lty = c(1,2),
                                    lwd = 1,
                                    legend = TRUE,
                                    bty = "o",
                                    cex.legend = 0.75,
                                    loc.legend = "topright") 
{
  object <- x
  return.LRT <- object$return.LRT
  double.bootstrap <- object$double.bootstrap
  double.bootstrap.alpha <- object$double.bootstrap.alpha
  
  stopifnot(type %in% c("all", "LRT.A", "LRT.B", 
                        "ppvalues.A", "ppvalues.B"))
  if (type == "ppvalues.A" || type == "ppvalues.B") 
    stopifnot(double.bootstrap == "standard")
  if (type == "LRT.A" || type == "LRT.B") stopifnot(return.LRT)
  if (type == "all" & !return.LRT) stopifnot (double.bootstrap != "FDB")
  
  pvalue <- rep(as.numeric(NA), 2) 
  pvalue[1]   <- object$lrt.bootA[1]
  pvalue[2]   <- object$lrt.bootB[1]
  
  y.lab <- ylab
  
  if (return.LRT) {
    lrt.obs <- rep(as.numeric(NA), 2) 
    lrt.obs[1]  <- attr(object$lrt.bootA, "LRT.original")
    lrt.obs[2]  <- attr(object$lrt.bootB, "LRT.original")
    
    lrt.A <- attr(object$lrt.bootA, "LRT")
    lrt.B <- attr(object$lrt.bootB, "LRT")
    
    if (length(lrt.A) - length(lrt.B) < 0L) {
      lrt <- cbind(c(lrt.A, rep(NA, length(lrt.B) - length(lrt.A))), lrt.B)
    }
    else { 
      lrt <- cbind(lrt.A, c(lrt.B, rep(NA, length(lrt.A) - length(lrt.B))))
    }
    
    if (xlab == "xlabel") { 
      x.lrt <- c("Bootstrapped LR statistic values")
    }
    else {
      x.lrt <- xlab
    }
    
    if (main == "main") {
      m.lrt <- c("Distribution of LR statistic values - Type A", 
                 "Distribution of LR statistic values - Type B")
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
      x.lrt <- c("Bootstrapped LR statistic values")
    }
    else {
      x.lrt <- xlab
    }
    
    if (main == "main") {
      m.lrt <- c("Distribution of LR statistic values - Type A", 
                 "Distribution of LR statistic values - Type B")
    }
    else {
      m.lrt <- main
    }
  }
  
  if (double.bootstrap == "standard") {
    ppvalue.A <- attr(object$lrt.bootA, "plugin.pvalues")
    ppvalue.B <- attr(object$lrt.bootB, "plugin.pvalues")
    adj.a <- rep(as.numeric(NA), 2)
    adj.a[1] <- quantile(ppvalue.A, double.bootstrap.alpha)
    adj.a[2] <- quantile(ppvalue.B, double.bootstrap.alpha)
    adj.ppv <- rep(as.numeric(NA), 2)
    adj.ppv[1] <- attr(object$lrt.bootA, "adj.pvalue")
    adj.ppv[2] <- attr(object$lrt.bootB, "adj.pvalue")
    
    if (length(ppvalue.A) - length(ppvalue.B) < 0L) {
      ppv <- cbind(c(ppvalue.A, rep(NA, length(ppvalue.B) - 
                                      length(ppvalue.A))), ppvalue.B)
    }
    else { 
      ppv <- cbind(ppvalue.A, c(ppvalue.B, rep(NA, length(ppvalue.A) - 
                                                 length(ppvalue.B))))
    }
    
    if (xlab == "xlabel") {
      x.ppv  <- c("Bootstrapped plug-in p-values")
    }
    else {
      x.ppv <- xlab    
    }
    
    if (main == "main") {
      m.ppv  <- c("Distribution of plug-in p-values - Type A", 
                  "Distribution of plug-in p-values - Type B")
    }
    else {
      m.ppv <- main
    }  
  }
  
  if (return.LRT & type == "all" & double.bootstrap != "standard") {
    par(mfrow = c(1, 2))
  }
  else if (return.LRT & type == "all" & double.bootstrap == "standard") {
    par(mfrow = c(2, 2))
  }
  else if (!return.LRT & (double.bootstrap == "standard" | 
                            double.bootstrap == "FDB")) {
    par(mfrow = c(1, 2))
  }
  else if (type != "all") {
    par(mfrow = c(1, 1))
  }
  
  if (double.bootstrap == "standard") {
    if (type == "LRT.A" | type == "LRT.B") double.bootstrap = "no"
  }
  
  if (type == "ppvalues.A" | type == "ppvalues.B") return.LRT <- FALSE
  
  if ((type == "LRT.A" & return.LRT) |
        (type == "ppvalues.A" & double.bootstrap == "standard")) { 
    a = 1L
    b = 1L
  }
  else if ((type == "LRT.B" & return.LRT) |
             (type == "ppvalues.B" & double.bootstrap == "standard")) {
    a = 2L
    b = 2L
  }
  else if (type == "all") {
    a = 1L
    b = 2L
  }  
  
  for (i in a:b) {
    if (return.LRT) {
      plot <- hist(lrt[,i], plot = FALSE, breaks=breaks) 
      plot(plot, freq = freq, main = m.lrt[i], xlab = x.lrt, ylab = y.lab, 
           cex.axis = cex.axis, cex.main = cex.main, cex.lab = cex.lab, 
           col = col, border = border, axes = axes, ...) 
      
      if (vline) abline(v = lrt.obs[i], col = vline.col[1], lty = lty[1], 
                        lwd = lwd)
      if (vline & double.bootstrap == "FDB") abline(v = lrt.q[i], 
                                                    col = vline.col[2], 
                                                    lty = lty[2], lwd = lwd)
      
      if (legend & double.bootstrap != "FDB") {
        ppval    <- sprintf("  %-15s  %.2f", "plug-in p value", pvalue[i]) 
        obs.lrt  <- sprintf("  %-15s  %05.2f", "observed LRT", lrt.obs[i])
        legend.obj <- c(obs.lrt, ppval)
        
        if (!vline){
          legend(loc.legend, legend.obj, lty = c(0, 0),   
                 lwd = lwd, cex = cex.legend, bty = bty)
        }
        else {
          legend(loc.legend, legend.obj, lty = c(lty[1], 0), col = vline.col[1],  
                 lwd = lwd, cex = cex.legend, bty = bty)
        }
      }
      else if (legend & double.bootstrap == "FDB") {
        obs.lrt    <- sprintf("  %-8s  %05.2f", "observed LRT", lrt.obs[i])
        obs.lrt.q  <- sprintf("  %-8s  %05.2f", "LRT.Q", lrt.q[i])
        ppval      <- sprintf("  %-8s  %.2f", "FDB p-value", adj.pvalue[i]) 
        legend.obj <- c(obs.lrt, obs.lrt.q, ppval)
        
        if (!vline) {
          legend(loc.legend, legend.obj, lty = c(0,0,0), 
                 col = c(vline.col, "blue"), lwd = lwd, cex = cex.legend, 
                 bty = bty)
        }
        else {
          legend(loc.legend, legend.obj, lty = c(lty[1],lty[2],0), 
                 col = c(vline.col[1], vline.col[2]), lwd = lwd, 
                 cex = cex.legend, bty = bty)
        }
      }
    }
    
    if (double.bootstrap == "standard") {      
      plot <- hist(ppv[, i], plot = FALSE, breaks=breaks) 
      plot(plot, freq = TRUE, main = m.ppv[i], xlab = x.ppv, ylab = y.lab, 
           cex.axis = cex.axis, cex.main = cex.main, cex.lab = cex.lab, 
           xlim = range(0:1), col = col, border = border, axes = axes, ...) 
      if (vline) {
        abline(v = adj.a[i], col = vline.col[1], lty = lty[1], lwd = lwd)
        abline(v = adj.ppv[i], col = vline.col[2], lty = lty[2], lwd = lwd)
      }
      
      if (legend) {
        adja    <- sprintf("  %-10s  %.2f", "Adj.alpha", adj.a[i]) 
        adjp    <- sprintf("  %-10s  %.2f", "Adj.p-value", adj.ppv[i])
        legend.obj <- c(adja, adjp)
        if (!vline) {
          legend(loc.legend, legend.obj, lty = 0, col = vline.col,  
                 lwd = lwd, cex = cex.legend, bty = bty)
        }
        else {
          legend(loc.legend, legend.obj, lty = lty, col = vline.col,  
                 lwd = lwd, cex = cex.legend, bty = bty)
        }
      }
    }
  }
}
