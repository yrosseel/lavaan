# This code is contributed by Leonard Vanbrabant <L.G.F.Vanbrabant@gmail.com>
InformativeTesting <- function(model = NULL, data, constraints = NULL,
                               R = 1000L, type = "bollen.stine",
                               return.LRT = TRUE,
                               double.bootstrap = "standard",
                               double.bootstrap.R = 249,
                               double.bootstrap.alpha = 0.05,
                               parallel = c("no", "multicore", "snow"),
                               ncpus = 1L, cl = NULL, verbose = FALSE, ...) {
  fit.B1 <- sem(model, ...,
    data = data,
    se   = "none",
    test = "standard"
  )

  fit.B0 <- fit.A1 <- sem(model, ...,
    data        = data,
    se          = "none",
    test        = "standard",
    constraints = constraints
  )

  # con.idx  <- (max(fit.B1@ParTable$id) + 1L):max(fit.A1@ParTable$id)
  #
  # user.equal  <- fit.A1@ParTable
  # user.equal$op[con.idx] <- "=="

  user.equal <- fit.A1@ParTable
  CON <- attr(
    lavParseModelString(constraints, parser = fit.B1@Options$parser),
    "constraints"
  )
  for (con in 1:length(CON)) {
    if (CON[[con]]$op %in% c("<", ">")) {
      this.lhs <- CON[[con]]$lhs
      this.op <- CON[[con]]$op
      this.rhs <- CON[[con]]$rhs

      # find this line in user.equal@ParTable
      idx <- which(
        user.equal$lhs == this.lhs,
        user.equal$op == this.op,
        user.equal$rhs == this.rhs
      )
      if (length(idx) == 0L) { # not found, give warning?
        lav_msg_stop(gettext("no inequality constraints (<, >) found."))
      }

      # change op to ==
      user.equal$op[idx] <- "=="
    }
  }

  fit.A0 <- sem(user.equal, ...,
    data = data,
    se   = "none",
    test = "standard"
  )

  lrt.bootA <- bootstrapLRT(fit.A0, fit.A1,
    R                      = R,
    type                   = type,
    verbose                = verbose,
    return.LRT             = return.LRT,
    double.bootstrap       = double.bootstrap,
    double.bootstrap.R     = double.bootstrap.R,
    double.bootstrap.alpha = double.bootstrap.alpha,
    parallel               = parallel,
    ncpus                  = ncpus,
    cl                     = cl
  )

  lrt.bootB <- bootstrapLRT(fit.B0, fit.B1,
    R                      = R,
    type                   = type,
    verbose                = verbose,
    return.LRT             = return.LRT,
    double.bootstrap       = double.bootstrap,
    double.bootstrap.R     = double.bootstrap.R,
    double.bootstrap.alpha = double.bootstrap.alpha,
    parallel               = parallel,
    ncpus                  = ncpus,
    cl                     = cl
  )

  output <- list(
    fit.A0 = fit.A0, fit.A1 = fit.A1, fit.B1 = fit.B1,
    lrt.bootA = lrt.bootA, lrt.bootB = lrt.bootB,
    double.bootstrap = double.bootstrap,
    double.bootstrap.alpha = double.bootstrap.alpha,
    return.LRT = return.LRT, type = type
  )

  class(output) <- "InformativeTesting"

  return(output)
}


print.InformativeTesting <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  object <- x
  cat("\nInformativeTesting: Order/Inequality Constrained Hypothesis Testing:\n\n")
  cat("  Variable names in model         :", unlist(object$fit.A1@Data@ov.names[1]), "\n")
  cat("  Number of variables             :", object$fit.A1@Model@nvar[1], "\n")
  cat("  Number of groups                :", object$fit.A1@Data@ngroups, "\n")
  cat("  Used sample size per group      :", unlist(object$fit.A1@Data@nobs), "\n")
  cat("  Used sample size                :", sum(unlist(object$fit.A1@Data@nobs)), "\n")
  cat("  Total sample size               :", sum(unlist(object$fit.A1@Data@norig)), "\n\n")
  cat("  Estimator                       :", object$fit.A1@Options$estimator, "\n")
  cat("  Missing data                    :", object$fit.A1@Options$missing, "\n")
  cat("  Bootstrap method                :", object$type, "\n")
  cat("  Double bootstrap method         :", object$double.bootstrap, "\n")

  dbtype <- object$double.bootstrap
  # original LRT for hypothesis test Type A
  TsA <- attr(object$lrt.bootA, "LRT.original")
  # original LRT for hypothesis test Type B
  TsB <- attr(object$lrt.bootB, "LRT.original")
  # unadjusted pvalues for Ts
  pvalueA <- object$lrt.bootA[1]
  pvalueB <- object$lrt.bootB[1]
  alpha <- object$double.bootstrap.alpha

  ###
  if (dbtype == "no") {
    cat(
      "\n\n  Type A test: H0: all restriktions active (=)", "\n",
      "          vs. H1: at least one restriktion strictly true (>)", "\n"
    )
    cat("         Test statistic: ", format(round(TsA, digits), nsmall = digits), ", unadjusted p-value: ",
      if (pvalueA < 1e-04) {
        "<0.0001"
      } else {
        format(round(pvalueA, digits), nsmall = digits)
      }, " (alpha = ", alpha, ") ", "\n\n",
      sep = ""
    )
    cat(
      "  Type B test: H0: all restriktions true", "\n",
      "          vs. H1: at least one restriktion false", "\n"
    )
    cat("         Test statistic: ", format(round(TsB, digits), nsmall = digits), ", unadjusted p-value: ",
      if (pvalueB < 1e-04) {
        "<0.0001"
      } else {
        format(round(pvalueB, digits), nsmall = digits)
      }, " (alpha = ", alpha, ") ", "\n",
      sep = ""
    )
  } else if (dbtype == "FDB") {
    # adjusted pvalues for Ts
    adj.pvalueA <- attr(object$lrt.bootA, "adj.pvalue")
    adj.pvalueB <- attr(object$lrt.bootB, "adj.pvalue")
    cat(
      "\n\n  Type A test: H0: all restriktions active (=)", "\n",
      "          vs. H1: at least one restriktion strictly true (>)", "\n"
    )
    cat("         Test statistic: ", format(round(TsA, digits), nsmall = digits), ", adjusted p-value: ",
      if (adj.pvalueA < 1e-04) {
        "<0.0001"
      } else {
        format(round(adj.pvalueA, digits), nsmall = digits)
      }, " (alpha = ",
      alpha, ") ", "\n\n",
      sep = ""
    )
    cat(
      "  Type B test: H0: all restriktions true", "\n",
      "          vs. H1: at least one restriktion false", "\n"
    )
    cat("         Test statistic: ", format(round(TsB, digits),
      nsmall = digits
    ), ", adjusted p-value: ",
    if (adj.pvalueB < 1e-04) {
      "<0.0001"
    } else {
      format(round(adj.pvalueB, digits), nsmall = digits)
    }, " (alpha = ",
    alpha, ") ", "\n",
    sep = ""
    )
  } else if (dbtype == "standard") {
    # adjusted nominal levels
    adj.alphaA <- attr(object$lrt.bootA, "adj.alpha")
    adj.alphaB <- attr(object$lrt.bootB, "adj.alpha")
    # adjusted pvalues for Ts
    adj.pvalueA <- attr(object$lrt.bootA, "adj.pvalue")
    adj.pvalueB <- attr(object$lrt.bootB, "adj.pvalue")
    cat(
      "\n\n  Type A test: H0: all restriktions active (=)", "\n",
      "          vs. H1: at least one restriktion strictly true (>)", "\n"
    )
    cat("         Test statistic: ", format(round(TsA, digits),
      nsmall = digits
    ), ", adjusted p-value: ",
    if (adj.pvalueA < 1e-04) {
      "<0.0001"
    } else {
      format(round(adj.pvalueA, digits), nsmall = digits)
    }, " (alpha = ", alpha, ") ", "\n",
    sep = ""
    )
    cat("                                ", "unadjusted p-value: ",
      if (pvalueA < 1e-04) {
        "<0.0001"
      } else {
        format(round(pvalueA, digits), nsmall = digits)
      }, " (alpha = ",
      format(round(adj.alphaA, digits), nsmall = digits), ") ", "\n\n",
      sep = ""
    )
    cat(
      "  Type B test: H0: all restriktions true", "\n",
      "          vs. H1: at least one restriktion false", "\n"
    )
    cat("         Test statistic: ", format(round(TsB, digits), nsmall = digits), ", adjusted p-value: ",
      if (adj.pvalueB < 1e-04) {
        "<0.0001"
      } else {
        format(round(adj.pvalueB, digits), nsmall = digits)
      }, " (alpha = ", alpha, ") ", "\n",
      sep = ""
    )
    cat("                               ", "unadjusted p-value: ",
      if (pvalueB < 1e-04) {
        "<0.0001"
      } else {
        format(round(pvalueB, digits), nsmall = digits)
      }, " (alpha = ",
      format(round(adj.alphaB, digits), nsmall = digits), ") ", "\n\n",
      sep = ""
    )
  }

  if (dbtype == "no") {
    cat("\n  No double bootstrap method is set. The results may be spurious.\n\n")
  }
}


plot.InformativeTesting <- function(x, ...,
                                    type = c("lr", "ppv"),
                                    main = "main",
                                    xlab = "xlabel",
                                    ylab = "Frequency",
                                    freq = TRUE,
                                    breaks = 15,
                                    cex.main = 1,
                                    cex.lab = 1,
                                    cex.axis = 1,
                                    col = "grey",
                                    border = par("fg"),
                                    vline = TRUE,
                                    vline.col = c("red", "blue"),
                                    lty = c(1, 2),
                                    lwd = 1,
                                    legend = TRUE,
                                    bty = "o",
                                    cex.legend = 1,
                                    loc.legend = "topright") {
  object <- x
  return.LRT <- object$return.LRT
  double.bootstrap <- object$double.bootstrap
  double.bootstrap.alpha <- object$double.bootstrap.alpha
  pvalue <- c(object$lrt.bootA[1], object$lrt.bootB[1])

  par(mfrow = c(1, 2))
  if (length(type) == 2) {
    par(mfrow = c(2, 2))
  }

  if (return.LRT && (type == "lr" || length(type) == 2)) {
    lrt.obs <- c(
      attr(object$lrt.bootA, "LRT.original"),
      attr(object$lrt.bootB, "LRT.original")
    )
    lrt.A <- attr(object$lrt.bootA, "LRT")
    lrt.B <- attr(object$lrt.bootB, "LRT")
    if (length(lrt.A) - length(lrt.B) < 0L) {
      lrt <- as.data.frame(cbind(c(lrt.A, rep(as.numeric(NA), length(lrt.B) -
        length(lrt.A))), lrt.B))
    } else {
      lrt <- as.data.frame(cbind(lrt.A, c(lrt.B, rep(
        as.numeric(NA),
        length(lrt.A) -
          length(lrt.B)
      ))))
    }
    names(lrt) <- c("lrt.A", " lrt.B")

    if (xlab == "xlabel") {
      xlab.lrt <- c("Bootstrapped LR values")
    }
    if (main == "main") {
      main.lrt <- c(
        "Distr. of LR values - Type A",
        "Distr. of LR values - Type B"
      )
    }

    for (i in 1:2) {
      plot <- hist(lrt[, i], plot = FALSE, breaks = breaks)
      plot(plot, ...,
        freq     = freq,
        main     = main.lrt[i],
        xlab     = xlab.lrt,
        ylab     = ylab,
        cex.axis = cex.axis,
        cex.main = cex.main,
        cex.lab  = cex.lab,
        col      = col,
        border   = border,
        axes     = FALSE,
        xaxt     = "n"
      )

      axis(side = 1)
      axis(side = 2)
      box(lty = 1, col = "black")

      if (vline) {
        abline(
          v = lrt.obs[i],
          col = vline.col[1],
          lty = lty[1],
          lwd = lwd
        )
      }
      if (legend) {
        ppvalue <- sprintf("%.2f", pvalue[i])
        obs.lrt <- sprintf("%.2f", lrt.obs[i])
        ppval <- paste0("plug-in p value = ", ppvalue)
        obs.lrt <- paste0("observed LR = ", obs.lrt)
        legend.obj <- c(obs.lrt, ppval)
        if (!vline) {
          legend(loc.legend, legend.obj,
            lty = c(0, 0),
            lwd = lwd,
            cex = cex.legend,
            bty = bty
          )
        } else {
          legend(loc.legend, legend.obj,
            lty = c(lty[1], 0),
            col = vline.col[1],
            lwd = lwd,
            cex = cex.legend,
            bty = bty
          )
        }
      }
    }
  }

  if (double.bootstrap == "standard" && (type == "ppv" || length(type) == 2)) {
    ppvalue.A <- attr(object$lrt.bootA, "plugin.pvalues")
    ppvalue.B <- attr(object$lrt.bootB, "plugin.pvalues")
    adj.a <- c(
      quantile(ppvalue.A, double.bootstrap.alpha),
      quantile(ppvalue.B, double.bootstrap.alpha)
    )
    adj.ppv <- c(
      attr(object$lrt.bootA, "adj.pvalue"),
      attr(object$lrt.bootB, "adj.pvalue")
    )
    if (length(ppvalue.A) - length(ppvalue.B) < 0L) {
      ppv <- as.data.frame(cbind(c(ppvalue.A, rep(NA, length(ppvalue.B) -
        length(ppvalue.A))), ppvalue.B))
    } else {
      ppv <- as.data.frame(cbind(ppvalue.A, c(ppvalue.B, rep(NA, length(ppvalue.A) -
        length(ppvalue.B)))))
    }
    names(ppv) <- c("ppA", "ppB")

    if (xlab == "xlabel") {
      xlab.ppv <- c("Bootstrapped plug-in p-values")
    }
    if (main == "main") {
      main.ppv <- c(
        "Distr. of plug-in p-values - Type A",
        "Distr. of plug-in p-values - Type B"
      )
    }

    for (i in 1:2) {
      plot <- hist(ppv[, i], plot = FALSE, breaks = breaks)
      plot(plot, ...,
        freq     = freq,
        main     = main.ppv[i],
        xlab     = xlab.ppv,
        ylab     = ylab,
        cex.axis = cex.axis,
        cex.main = cex.main,
        cex.lab  = cex.lab,
        col      = col,
        border   = border,
        axes     = FALSE,
        xaxt     = "n"
      )

      axis(side = 1, at = seq(0, 1, 0.1))
      axis(side = 2)
      box(lty = 1, col = "black")
      if (vline) {
        abline(
          v = adj.a[i],
          col = vline.col[1],
          lty = lty[1],
          lwd = lwd
        )
        abline(
          v = adj.ppv[i],
          col = vline.col[2],
          lty = lty[2],
          lwd = lwd
        )
      }
      if (legend) {
        adj.alpha <- sprintf("%.2f", adj.a[i])
        adj.pval <- sprintf("%.2f", adj.ppv[i])
        adja <- paste0("Adjusted alpha = ", adj.alpha)
        adjp <- paste0("Adjusted p-value = ", adj.pval)
        legend.obj <- c(adja, adjp)
        if (!vline) {
          legend(loc.legend, legend.obj,
            lty = 0,
            col = vline.col,
            lwd = lwd,
            cex = cex.legend,
            bty = bty
          )
        } else {
          legend(loc.legend, legend.obj,
            lty = lty,
            col = vline.col,
            lwd = lwd,
            cex = cex.legend,
            bty = bty
          )
        }
      }
    }
  }
}
