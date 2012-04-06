# This code is contributed by Leonard Vanbrabant <L.G.F.Vanbrabant@hotmail.com>
InformativeTesting <- function(model = NULL, data, constraints = NULL, 
                               R = NULL, type = "bollen.stine",
                               return.LRT = TRUE, 
                               double.bootstrap = "FDB",
                               double.bootstrap.R = 500L,
                               double.bootstrap.alpha = 0.05,
                               parallel = c("no", "multicore", "snow"), 
                               ncpus = 1L, cl = NULL, verbose = FALSE, ...){
  
  if (missing(R) & double.bootstrap == "standard" |
      double.bootstrap == "no") { R <- 1000L 
  }
  else if (missing(R) & double.bootstrap == "FDB") { R <- 2000L 
  }
  
  fit.free <- sem(model, data = data, se = "none", test = "standard", ...) 
  fit.ineq <- sem(model, data = data, se = "none", test = "standard", 
                  constraints = constraints, ...) 
  
  con.idx  <- (max(fit.free@User$id) + 1L):max(fit.ineq@User$id)
  
  user.equal  <- fit.ineq@User
  user.equal$op[con.idx] <- "=="
  
  fit.equal <- sem(user.equal, data = data, se="none", 
                   test = "standard", ...)
   
  lrt.bootA <- bootstrapLRT(fit.equal, fit.ineq, 
                            R = R, type = type, verbose = verbose,
                            return.LRT = return.LRT, 
                            double.bootstrap = double.bootstrap,
                            double.bootstrap.R = double.bootstrap.R, 
                            double.bootstrap.alpha = double.bootstrap.alpha,
                            parallel = parallel, ncpus = ncpus, cl = cl)
  
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
                 double.bootstrap.alpha = double.bootstrap.alpha,
                 return.LRT = return.LRT, 
                 type = type)
  
  class(output) <- "InformativeTesting"
  
  return(output)
}




print.InformativeTesting <- function(x, ...) {
  object <- x
  cat("\n")
  p.A <- object$lrt.bootA[1]
  p.B <- object$lrt.bootB[1]
  sig.A <- sig.B <- "Non-significant"
 
  if (object$double.bootstrap == "FDB") {
    p.A <- attr(object$lrt.bootA, "adj.pvalue")
    p.B <- attr(object$lrt.bootB, "adj.pvalue")
  }
  
  if (object$double.bootstrap == "standard") {
    p.adj.A <- attr(object$lrt.bootA, "adj.pvalue")
    p.adj.B <- attr(object$lrt.bootB, "adj.pvalue")
    alpha.adj.A <- attr(object$lrt.bootA, "adj.alpha")
    alpha.adj.B <- attr(object$lrt.bootB, "adj.alpha")
    alpha.A <- alpha.B <- object$double.bootstrap.alpha
    
    if (attr(object$lrt.bootA, "adj.pvalue") < object$double.bootstrap.alpha) 
      sig.A <- "Significant"
    if (attr(object$lrt.bootB, "adj.pvalue") < object$double.bootstrap.alpha) 
      sig.B <- "Significant"
    
    cat("Order Constrained Hypothesis Testing:\n\n\n")
    hd.txt <- sprintf(" %-16s %1s %-15s", "Double bootstrap", "=", 
                      object$double.bootstrap)
    ht.txt <- sprintf(" %-16s %1s %-15s", "Bootstrap type", "=", 
                      object$type)
    cat(ht.txt,"\n")
    cat(hd.txt,"\n\n\n")
    h0.txt <- sprintf("  %-22s  %8s  %15s", "", "Type A", "Type B") 
    cat(h0.txt, "\n")
    cat("  -------------------------------------------------\n")
    t0.txt <- sprintf("  %-22s", "LR statistic")
    t1.txt <- sprintf("  %2s %05.2f", "", attr(object$lrt.bootA, "LRT.original"))
    t2.txt <- sprintf("  %9s %05.2f", "", attr(object$lrt.bootB, "LRT.original"))
    cat(t0.txt, t1.txt, t2.txt, "\n\n", sep="")
    t0.txt <- sprintf("  %-22s", "Adjusted alpha")
    t1.txt <- sprintf("  %3s %.2f", "", alpha.adj.A)
    t2.txt <- sprintf("  %10s %.2f", "", alpha.adj.B)
    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
    t0.txt <- sprintf("  %-22s", "P-value")
    t1.txt <- sprintf("  %3s %.2f", "", p.A)
    t2.txt <- sprintf("  %10s %.2f", "", p.B)
    cat(t0.txt, t1.txt, t2.txt, "\n\n", sep="")
    t0.txt <- sprintf("  %-22s", "Alpha")
    t1.txt <- sprintf("  %3s %.2f", "", alpha.A)
    t2.txt <- sprintf("  %10s %.2f", "", alpha.B)
    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
    t0.txt <- sprintf("  %-22s", "Adjusted p-value")
    t1.txt <- sprintf("  %3s %.2f", "", p.adj.A)
    t2.txt <- sprintf("  %10s %.2f", "", p.adj.B)
    cat(t0.txt, t1.txt, t2.txt, "\n\n", sep="")
    t0.txt <- sprintf("  %10s", "Significance")
    t1.txt <- sprintf("  %18s", sig.A)
    t2.txt <- sprintf("  %15s", sig.B)
    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
  }
  else if (object$double.bootstrap == "FDB") {
    cat("Order Constrained Hypothesis Testing:\n\n\n")
    hd.txt <- sprintf(" %-16s %1s %-15s", "Double bootstrap", "=", 
                      object$double.bootstrap)
    ht.txt <- sprintf(" %-16s %1s %-15s", "Bootstrap type", "=", object$type)
    cat(ht.txt,"\n")
    cat(hd.txt,"\n\n")
    p.A <- attr(object$lrt.bootA, "adj.pvalue")
    p.B <- attr(object$lrt.bootB, "adj.pvalue")
    alpha.A <- alpha.B <- object$double.bootstrap.alpha
    
    if (attr(object$lrt.bootA, "adj.pvalue") < object$double.bootstrap.alpha) { 
      sig.A <- "Significant"
    }
    if (attr(object$lrt.bootB, "adj.pvalue") < object$double.bootstrap.alpha) { 
      sig.B <- "Significant"
    }
    h0.txt <- sprintf("  %-22s  %8s  %15s", "", "Type A", "Type B") 
    cat(h0.txt, "\n")
    cat("  -------------------------------------------------\n")
    t0.txt <- sprintf("  %-22s", "LR statistic")
    t1.txt <- sprintf("  %2s %05.2f", "", attr(object$lrt.bootA, "LRT.original"))
    t2.txt <- sprintf("  %9s %05.2f", "", attr(object$lrt.bootB, "LRT.original"))
    cat(t0.txt, t1.txt, t2.txt, "\n\n", sep="")
    t0.txt <- sprintf("  %-22s", "Alpha")
    t1.txt <- sprintf("  %3s %.2f", "", alpha.A)
    t2.txt <- sprintf("  %10s %.2f", "", alpha.B)
    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
    t0.txt <- sprintf("  %-22s", "Adjusted p-value")
    t1.txt <- sprintf("  %3s %.2f", "", p.A)
    t2.txt <- sprintf("  %10s %.2f", "", p.B)
    cat(t0.txt, t1.txt, t2.txt, "\n\n", sep="")
    t0.txt <- sprintf("  %10s", "Significance")
    t1.txt <- sprintf("  %18s", sig.A)
    t2.txt <- sprintf("  %15s", sig.B)
    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
  }
  else if (object$double.bootstrap == "no") {
    alpha.A <- alpha.B <- object$double.bootstrap.alpha
    sig.A <- sig.B <- "Inconclusive*"
    cat("Order Constrained Hypothesis Testing:\n\n\n")
    hd.txt <- sprintf(" %-16s %1s %-15s", "Double bootstrap", "=", 
                      object$double.bootstrap)
    ht.txt <- sprintf(" %-16s %1s %-15s", "Bootstrap type", "=", object$type)
    cat(ht.txt,"\n")
    cat(hd.txt,"\n\n\n")
    h0.txt <- sprintf("  %-22s  %8s  %15s", "", "Type A", "Type B") 
    cat(h0.txt, "\n")
    cat("  -------------------------------------------------\n")
    t0.txt <- sprintf("  %-22s", "LR statistic")
    t1.txt <- sprintf("  %2s %05.2f", "", attr(object$lrt.bootA, "LRT.original"))
    t2.txt <- sprintf("  %9s %05.2f", "", attr(object$lrt.bootB, "LRT.original"))
    cat(t0.txt, t1.txt, t2.txt, "\n\n", sep="")
    t0.txt <- sprintf("  %-22s", "Alpha")
    t1.txt <- sprintf("  %3s %.2f", "", alpha.A)
    t2.txt <- sprintf("  %10s %.2f", "", alpha.B)
    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
    t0.txt <- sprintf("  %-22s", "P-value")
    t1.txt <- sprintf("  %3s %.2f", "", p.A)
    t2.txt <- sprintf("  %10s %.2f", "", p.B)
    cat(t0.txt, t1.txt, t2.txt, "\n\n", sep="")
    t0.txt <- sprintf("  %10s", "Significance")
    t1.txt <- sprintf("  %18s", sig.A)
    t2.txt <- sprintf("  %15s", sig.B)
    cat(t0.txt, t1.txt, t2.txt, "\n\n", sep="")
    cat("  *For meaningfull results set double.bootstrap")
  }  
}


summary.InformativeTesting <- function(object, 
                                       fit = c("H0","Hi","Hu"), 
                                       fitmeasures = FALSE,
                                       ...)
{
  if (!missing(fit)) {
    fit <- tolower(fit) 
    stopifnot (fit %in% c("h0","hi","hu"))
    
    if ("h0" %in% fit) {
      fit.model <- sem(object$fitA1@User, 
                       group = object$fitA1@Options$group, ...)
    }
    else if ("hi" %in% fit) {
      fit.model <- sem(object$fitA2@User, 
                       group = object$fitA2@Options$group, ...)
    }
    else if ("hu" %in% fit) {
      fit.model <- sem(object$fitB2@User, 
                       group = object$fitB2@Options$group, ...)
    }
    
    if ("hu" %in% fit)     
      fit.m <- fitmeasures(fit.model, c("chisq", "df", "pvalue", "cfi",   
                                        "tli", "logl", "aic", "bic", "rmsea"))
    
    summary(fit.model)
    
    if (fitmeasures & "hu" %in% fit) { 
      cat("\n")
      cat("Fit Measures \n\n")
      print(fit.m)
    }
    else if (fitmeasures & ("h0" %in% fit | "hi" %in% fit)) {
      cat("\n")
      cat("Fit Measures only available for the unconstrained model Hu \n\n")
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
    cat("  Double bootstrap method       :", object$double.bootstrap, "\n\n")
    
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
                                    cex.axis = NULL,
                                    nclass = NULL,
                                    col = "grey",
                                    border = par("fg"),
                                    axes = TRUE,
                                    vline = FALSE, 
                                    vline.col = c("red", "blue"), 
                                    lty = c(1,2),
                                    lwd = 1,
                                    legend = FALSE,
                                    cex.legend = 0.75,
                                    loc.legend = "topright") 
{
  object <- x
  return.LRT <- object$return.LRT
  double.bootstrap <- object$double.bootstrap
  double.bootstrap.alpha <- object$double.bootstrap.alpha
  
  stopifnot(type %in% c("all", "LRT.A", "LRT.B", 
                        "ppvalues.A", "ppvalues.B"))
  if (type == "ppvalues.A" | type == "ppvalues.B") 
    stopifnot(double.bootstrap == "standard")
  if (type == "LRT.A" | type == "LRT.B") stopifnot(return.LRT)
  if (type == "all" & !return.LRT) stopifnot (double.bootstrap != "FDB")
  
  pvalue <- rep(as.numeric(NA), 2) 
  pvalue[1]   <- object$lrt.bootA[1]
  pvalue[2]   <- object$lrt.bootB[1]
  
  if (ylab == "ylabel") {
    y.lab  <- "Frequency"
  }
  else {
    y.lab <- ylab
  }
  
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
      x.ppv  <- c("Bootstrapped plug-in p values")
    }
    else {
      x.ppv <- xlab    
    }
    
    if (main == "main") {
      m.ppv  <- c("Distribution of plug-in p values - type A", 
                  "Distribution of plug-in p values - type B")
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
      plot <- hist(lrt[,i], nclass = nclass, plot = FALSE) 
      plot(plot, freq = freq, main = m.lrt[i], xlab = x.lrt, ylab = y.lab, 
           cex.axis = cex.axis, cex.main = cex.main, cex.lab = cex.lab, 
           col = col, border = border, axes = axes) 
      
      if (vline) abline(v = lrt.obs[i], col = vline.col[1], lty = lty[1], 
                        lwd = lwd)
      if (vline & double.bootstrap == "FDB") abline(v = lrt.q[i], 
                                                    col = vline.col[2], 
                                                    lty = lty[2], lwd = lwd)
      
      if (legend & double.bootstrap != "FDB") {
        ppval    <- sprintf("  %-15s  %.2f", "plug-in p value", pvalue[i]) 
        obs.lrt  <- sprintf("  %-15s  %05.2f", "obs LRT", lrt.obs[i])
        legend.obj <- c(obs.lrt, ppval)
        
        if (!vline){
          legend(loc.legend, legend.obj, lty = c(0, 0),   
                 lwd = lwd, cex = cex.legend, bty="n")
        }
        else {
          legend(loc.legend, legend.obj, lty = c(lty[1], 0), col = vline.col[1],  
                 lwd = lwd, cex = cex.legend, bty="n")
        }
      }
      else if (legend & double.bootstrap == "FDB") {
        obs.lrt    <- sprintf("  %-8s  %05.2f", "obs LRT", lrt.obs[i])
        obs.lrt.q  <- sprintf("  %-8s  %05.2f", "obs LRT.Q", lrt.q[i])
        ppval      <- sprintf("  %-8s  %.2f", "FDB p-value", adj.pvalue[i]) 
        legend.obj <- c(obs.lrt, obs.lrt.q, ppval)
        
        if (!vline) {
          legend(loc.legend, legend.obj, lty = c(0,0,0), 
                 col = c(vline.col, "blue"), lwd = lwd, cex = cex.legend, 
                 bty="n")
        }
        else {
          legend(loc.legend, legend.obj, lty = c(lty[1],lty[2],0), 
                 col = c(vline.col[1], vline.col[2]), lwd = lwd, 
                 cex = cex.legend, bty="n")
        }
      }
    }
    
    if (double.bootstrap == "standard") {      
      plot <- hist(ppv[, i], nclass = nclass, plot = FALSE) 
      plot(plot, freq = TRUE, main = m.ppv[i], xlab = x.ppv, ylab = y.lab, 
           cex.axis = cex.axis, cex.main = cex.main, cex.lab = cex.lab, 
           xlim = range(0:1), col = col, border = border, axes = axes) 
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
                 lwd = lwd, cex = cex.legend, bty="n")
        }
        else {
          legend(loc.legend, legend.obj, lty = lty, col = vline.col,  
                 lwd = lwd, cex = cex.legend, bty="n")
        }
      }
    }
  }
}
