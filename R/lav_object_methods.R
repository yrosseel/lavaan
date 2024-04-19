# `methods' for fitted lavaan objects
#
# standard (S4) methods:
# - show()
# - summary()
# - coef()
# - fitted.values() + fitted()
# - vcov()
# - logLik()
# - nobs()
# - update()
# - anova()

# lavaan-specific methods:
#
# - parameterEstimates()
# - standardizedSolution()
# - parameterTable()
# - varTable()


setMethod(
  "show", "lavaan",
  function(object) {
    # efa?
    efa.flag <- object@Options$model.type == "efa"

    # show only basic information
    res <- lav_object_summary(object,
      fit.measures = FALSE,
      estimates = FALSE,
      modindices = FALSE,
      efa = efa.flag
    )
    if (efa.flag) {
      # print (standardized) loadings only
      class(res) <- c("lavaan.efa", "list")
      print(res)
    } else {
      # print lavaan header
      print(res)
    }
    invisible(res)
  }
)

setMethod(
  "summary", "lavaan",
  function(object, header = TRUE,
           fit.measures = FALSE,
           estimates = TRUE,
           ci = FALSE,
           fmi = FALSE,
           standardized = FALSE,
           std = standardized,
           std.nox = FALSE, # TODO: remove deprecated argument in early 2025
           remove.step1 = TRUE,
           remove.unused = TRUE,
           cov.std = TRUE,
           rsquare = FALSE,
           fm.args = list(
             standard.test = "default",
             scaled.test = "default",
             rmsea.ci.level = 0.90,
             rmsea.h0.closefit = 0.05,
             rmsea.h0.notclosefit = 0.08,
             robust = TRUE,
             cat.check.pd = TRUE
           ),
           modindices = FALSE,
           nd = 3L, cutoff = 0.3, dot.cutoff = 0.1) {
    # efa?
    efa.flag <- object@Options$model.type == "efa"

    res <- lav_object_summary(
      object = object, header = header,
      fit.measures = fit.measures, estimates = estimates,
      ci = ci, fmi = fmi, std = std, standardized = standardized,
      remove.step1 = remove.step1, remove.unused = remove.unused,
      cov.std = cov.std,
      rsquare = rsquare, efa = efa.flag,
      fm.args = fm.args, modindices = modindices
    )
    # res has class c("lavaan.summary", "list")

    # what about nd? only used if we actually print; save as attribute
    attr(res, "nd") <- nd

    # if efa, add cutoff and dot.cutoff, and change class
    if (efa.flag) {
      # class(res) <- c("lavaan.summary.efa", "list")
      attr(res, "cutoff") <- cutoff
      attr(res, "dot.cutoff") <- dot.cutoff
    }

    res
  }
)


setMethod(
  "coef", "lavaan",
  function(object, type = "free", labels = TRUE) {
    lav_object_inspect_coef(
      object = object, type = type,
      add.labels = labels, add.class = TRUE
    )
  }
)

standardizedSolution <- # nolint
  standardizedsolution <- function(object, # nolint
                                   type = "std.all",
                                   se = TRUE,
                                   zstat = TRUE,
                                   pvalue = TRUE,
                                   ci = TRUE,
                                   level = 0.95,
                                   cov.std = TRUE,
                                   remove.eq = TRUE,
                                   remove.ineq = TRUE,
                                   remove.def = FALSE,
                                   partable = NULL,
                                   GLIST = NULL, # nolint
                                   est = NULL,
                                   output = "data.frame") {
    stopifnot(type %in% c("std.all", "std.lv", "std.nox"))

    # check output= argument
    output <- tolower(output)
    if (output %in% c("data.frame", "table")) {
      output <- "data.frame"
    } else if (output %in% c("text", "pretty")) {
      output <- "text"
    } else {
      lav_msg_stop(gettextf(
        "output must be %s or %s", sQuote("data.frame"), sQuote("text"))
      )
    }

    # no zstat + pvalue if estimator is Bayes
    if (object@Options$estimator == "Bayes") {
      zstat <- pvalue <- FALSE
    }

    # no se if class is not lavaan
    # using class() -- can't use inherits(), as this includes blavaan
    if (class(object)[1L] != "lavaan") {
      if (missing(se) || !se) {
        se <- FALSE
        zstat <- FALSE
        pvalue <- FALSE
      }
    }

    if (is.null(partable)) {
      tmp.partable <- inspect(object, "list")
    } else {
      tmp.partable <- partable
    }
    tmp.list <- tmp.partable[, c("lhs", "op", "rhs", "exo")]
    if (!is.null(tmp.partable$group)) {
      tmp.list$group <- tmp.partable$group
    }
    if (!is.null(tmp.partable$block)) {
      tmp.list$block <- tmp.partable$block
    }
    if (sum(nchar(tmp.partable$label)) != 0L) {
      tmp.list$label <- tmp.partable$label
    }

    # add std and std.all columns
    if (type == "std.lv") {
      tmp.list$est.std <- lav_standardize_lv(object,
        est = est, GLIST = GLIST,
        partable = partable, cov.std = cov.std
      )
    } else if (type == "std.all") {
      tmp.list$est.std <- lav_standardize_all(object,
        est = est, GLIST = GLIST,
        partable = partable, cov.std = cov.std
      )
    } else if (type == "std.nox") {
      tmp.list$est.std <- lav_standardize_all_nox(object,
        est = est, GLIST = GLIST,
        partable = partable, cov.std = cov.std
      )
    }

    if (object@Options$se != "none" && se) {
      # add 'se' for standardized parameters
      tmp.vcov <- try(lav_object_inspect_vcov(object,
        standardized = TRUE,
        type = type, free.only = FALSE,
        add.labels = FALSE,
        add.class = FALSE
      ))
      if (inherits(tmp.vcov, "try-error") || is.null(tmp.vcov)) {
        tmp.list$se <- rep(NA, length(tmp.list$lhs))
        if (zstat) {
          tmp.list$z <- rep(NA, length(tmp.list$lhs))
        }
        if (pvalue) {
          tmp.list$pvalue <- rep(NA, length(tmp.list$lhs))
        }
      } else {
        tmp <- diag(tmp.vcov)
        # catch negative values
        min.idx <- which(tmp < 0)
        if (length(min.idx) > 0L) {
          tmp[min.idx] <- as.numeric(NA)
        }
        # now, we can safely take the square root
        tmp <- sqrt(tmp)

        # catch near-zero SEs
        zero.idx <- which(tmp < .Machine$double.eps^(1 / 4)) # was 1/2 < 0.6
        # was 1/3 < 0.6-9
        if (length(zero.idx) > 0L) {
          tmp[zero.idx] <- 0.0
        }
        tmp.list$se <- tmp

        # add 'z' column
        if (zstat) {
          tmp.se <- ifelse(tmp.list$se == 0.0, NA, tmp.list$se)
          tmp.list$z <- tmp.list$est.std / tmp.se
        }
        if (zstat && pvalue) {
          tmp.list$pvalue <- 2 * (1 - pnorm(abs(tmp.list$z)))
        }
      }
    }

    # simple symmetric confidence interval
    if (se && object@Options$se != "none" && ci) {
      # next three lines based on confint.lm
      a <- (1 - level) / 2
      a <- c(a, 1 - a)
      fac <- qnorm(a)
      # if(object@Options$se != "bootstrap") {
      ci <- tmp.list$est.std + tmp.list$se %o% fac
      # } else {
      #    ci <- rep(as.numeric(NA), length(tmp.list$est.std)) +
      #              tmp.list$se %o% fac
      # }

      tmp.list$ci.lower <- ci[, 1]
      tmp.list$ci.upper <- ci[, 2]
    }


    # if single group, remove group column
    if (object@Data@ngroups == 1L) tmp.list$group <- NULL

    # remove == rows?
    if (remove.eq) {
      eq.idx <- which(tmp.list$op == "==")
      if (length(eq.idx) > 0L) {
        tmp.list <- tmp.list[-eq.idx, ]
      }
    }
    # remove <> rows?
    if (remove.ineq) {
      ineq.idx <- which(tmp.list$op %in% c("<", ">"))
      if (length(ineq.idx) > 0L) {
        tmp.list <- tmp.list[-ineq.idx, ]
      }
    }
    # remove := rows?
    if (remove.def) {
      def.idx <- which(tmp.list$op == ":=")
      if (length(def.idx) > 0L) {
        tmp.list <- tmp.list[-def.idx, ]
      }
    }

    # remove attribute for data order
    attr(tmp.list, "ovda") <- NULL

    if (output == "text") {
      class(tmp.list) <- c(
        "lavaan.parameterEstimates", "lavaan.data.frame",
        "data.frame"
      )
      # tmp.list$exo is needed for printing, don't remove it
      attr(tmp.list, "group.label") <- object@Data@group.label
      attr(tmp.list, "level.label") <- object@Data@level.label
      # attr(tmp.list, "header") <- FALSE
    } else {
      tmp.list$exo <- NULL
      tmp.list$block <- NULL
      class(tmp.list) <- c("lavaan.data.frame", "data.frame")
    }

    tmp.list
  }

parameterEstimates <- # nolint
  parameterestimates <- function(object,
                                 # select columns
                                 se = TRUE,
                                 zstat = TRUE,
                                 pvalue = TRUE,
                                 ci = TRUE,
                                 standardized = FALSE,
                                 fmi = FALSE,
                                 # control
                                 level = 0.95,
                                 boot.ci.type = "perc",
                                 cov.std = TRUE,
                                 fmi.options = list(),
                                 # add rows
                                 rsquare = FALSE,
                                 # remove rows
                                 remove.system.eq = TRUE,
                                 remove.eq = TRUE,
                                 remove.ineq = TRUE,
                                 remove.def = FALSE,
                                 remove.nonfree = FALSE,
                                 remove.step1 = TRUE,
                                 remove.unused = FALSE,
                                 # output
                                 add.attributes = FALSE,
                                 output = "data.frame",
                                 header = FALSE) {
    if (inherits(object, "lavaan.fsr")) {
      return(object$PE)
    }

    # deprecated add.attributes (for psycho/blavaan)
    if (add.attributes) {
      output <- "text"
    }


    # no se if class is not lavaan
    # can't use inherits(), as this would return TRUE if object is from blavaan
    if (class(object)[1L] != "lavaan") {
      if (missing(se) || !se) {
        se <- FALSE
        zstat <- FALSE
        pvalue <- FALSE
      }
    }

    # check output= argument
    output <- tolower(output)
    if (output %in% c("data.frame", "table")) {
      output <- "data.frame"
      header <- FALSE
    } else if (output %in% c("text", "pretty")) {
      output <- "text"
    } else {
      lav_msg_stop(gettextf(
        "output must be %s or %s", sQuote("data.frame"), sQuote("text"))
      )
    }

    # check fmi
    if (fmi) {
      if (inherits(object, "lavaanList")) {
        lav_msg_warn(gettext(
          "fmi not available for object of class \"lavaanList\""))
        fmi <- FALSE
      }
      if (object@Options$se != "standard") {
        lav_msg_warn(gettext(
          "fmi only available if se = \"standard\""))
        fmi <- FALSE
      }
      if (object@Options$estimator != "ML") {
        lav_msg_warn(gettext(
          "fmi only available if estimator = \"ML\""))
        fmi <- FALSE
      }
      if (!object@SampleStats@missing.flag) {
        lav_msg_warn(gettext(
          "fmi only available if missing = \"(fi)ml\""))
        fmi <- FALSE
      }
      if (!object@optim$converged) {
        lav_msg_warn(gettext(
          "fmi not available; model did not converge"))
        fmi <- FALSE
      }
    }

    # no zstat + pvalue if estimator is Bayes
    if (object@Options$estimator == "Bayes") {
      zstat <- pvalue <- FALSE
    }

    tmp.partable <- as.data.frame(object@ParTable, stringsAsFactors = FALSE)
    tmp.list <- tmp.partable[, c("lhs", "op", "rhs", "free")]
    if (!is.null(tmp.partable$user)) {
      tmp.list$user <- tmp.partable$user
    }
    if (!is.null(tmp.partable$block)) {
      tmp.list$block <- tmp.partable$block
    } else {
      tmp.list$block <- rep(1L, length(tmp.list$lhs))
    }
    if (!is.null(tmp.partable$level)) {
      tmp.list$level <- tmp.partable$level
    } else {
      tmp.list$level <- rep(1L, length(tmp.list$lhs))
    }
    if (!is.null(tmp.partable$group)) {
      tmp.list$group <- tmp.partable$group
    } else {
      tmp.list$group <- rep(1L, length(tmp.list$lhs))
    }
    if (!is.null(tmp.partable$step)) {
      tmp.list$step <- tmp.partable$step
    }
    if (!is.null(tmp.partable$efa)) {
      tmp.list$efa <- tmp.partable$efa
    }
    if (!is.null(tmp.partable$label)) {
      tmp.list$label <- tmp.partable$label
    } else {
      tmp.list$label <- rep("", length(tmp.list$lhs))
    }
    if (!is.null(tmp.partable$exo)) {
      tmp.list$exo <- tmp.partable$exo
    } else {
      tmp.list$exo <- rep(0L, length(tmp.list$lhs))
    }
    if (inherits(object, "lavaanList")) {
      # per default: nothing!
      # if("partable" %in% object@meta$store.slots) {
      #    COF <- sapply(object@ParTableList, "[[", "est")
      #    tmp.list$est <- rowMeans(COF)
      # }
      tmp.list$est <- NULL
    } else if (!is.null(tmp.partable$est)) {
      tmp.list$est <- tmp.partable$est
    } else {
      tmp.list$est <- lav_model_get_parameters(object@Model,
        type = "user",
        extra = TRUE
      )
    }
    if (!is.null(tmp.partable$lower)) {
      tmp.list$lower <- tmp.partable$lower
    }
    if (!is.null(tmp.partable$upper)) {
      tmp.list$upper <- tmp.partable$upper
    }


    # add se, zstat, pvalue
    if (se && object@Options$se != "none") {
      tmp.list$se <- lav_object_inspect_se(object)
      # handle tiny SEs
      tmp.list$se <- ifelse(tmp.list$se < sqrt(.Machine$double.eps),
        0, tmp.list$se
      )
      tmp.se <- ifelse(tmp.list$se < sqrt(.Machine$double.eps), NA, tmp.list$se)
      if (zstat) {
        tmp.list$z <- tmp.list$est / tmp.se
        if (pvalue) {
          tmp.list$pvalue <- 2 * (1 - pnorm(abs(tmp.list$z)))
          # remove p-value if bounds have been used
          if (!is.null(tmp.partable$lower)) {
            b.idx <- which(abs(tmp.partable$lower - tmp.partable$est) <
              sqrt(.Machine$double.eps) &
              tmp.partable$free > 0L)
            if (length(b.idx) > 0L) {
              tmp.list$pvalue[b.idx] <- as.numeric(NA)
            }
          }
          if (!is.null(tmp.partable$upper)) {
            b.idx <- which(abs(tmp.partable$upper - tmp.partable$est) <
              sqrt(.Machine$double.eps) &
              tmp.partable$free > 0L)
            if (length(b.idx) > 0L) {
              tmp.list$pvalue[b.idx] <- as.numeric(NA)
            }
          }
        }
      }
    }

    # extract bootstrap data (if any)
    if (object@Options$se == "bootstrap" ||
      "bootstrap" %in% object@Options$test ||
      "bollen.stine" %in% object@Options$test) {
      tmp.boot <- lav_object_inspect_boot(object)
      bootstrap.seed <- attr(tmp.boot, "seed") # for bca
      error.idx <- attr(tmp.boot, "error.idx")
      if (length(error.idx) > 0L) {
        tmp.boot <- tmp.boot[-error.idx, , drop = FALSE] # drops attributes
      }
    } else {
      tmp.boot <- NULL
    }

    bootstrap.successful <- NROW(tmp.boot) # should be zero if NULL

    # confidence interval
    if (se && object@Options$se != "none" && ci) {
      # next three lines based on confint.lm
      a <- (1 - level) / 2
      a <- c(a, 1 - a)
      if (object@Options$se != "bootstrap") {
        fac <- qnorm(a)
        ci <- tmp.list$est + tmp.list$se %o% fac
      } else if (object@Options$se == "bootstrap") {
        # local copy of 'norm.inter' from boot package (not exported!)
        norm.inter <- function(t, alpha) {
          t <- t[is.finite(t)]
          tmp.r <- length(t)
          rk <- (tmp.r + 1) * alpha
          if (!all(rk > 1 & rk < tmp.r)) {
            lav_msg_warn(gettext("extreme order statistics used as endpoints"))
          }
          k <- trunc(rk)
          inds <- seq_along(k)
          out <- inds
          kvs <- k[k > 0 & k < tmp.r]
          tstar <- sort(t, partial = sort(union(c(1, tmp.r), c(kvs, kvs + 1))))
          ints <- (k == rk)
          if (any(ints)) out[inds[ints]] <- tstar[k[inds[ints]]]
          out[k == 0] <- tstar[1L]
          out[k == tmp.r] <- tstar[tmp.r]
          not <- function(v) xor(rep(TRUE, length(v)), v)
          temp <- inds[not(ints) & k != 0 & k != tmp.r]
          temp1 <- qnorm(alpha[temp])
          temp2 <- qnorm(k[temp] / (tmp.r + 1))
          temp3 <- qnorm((k[temp] + 1) / (tmp.r + 1))
          tk <- tstar[k[temp]]
          tk1 <- tstar[k[temp] + 1L]
          out[temp] <- tk + (temp1 - temp2) / (temp3 - temp2) * (tk1 - tk)
          cbind(round(rk, 2), out)
        }

        stopifnot(!is.null(tmp.boot))
        stopifnot(boot.ci.type %in% c(
          "norm", "basic", "perc",
          "bca.simple", "bca"
        ))
        if (boot.ci.type == "norm") {
          fac <- qnorm(a)
          boot.x <- colMeans(tmp.boot, na.rm = TRUE)
          boot.est <-
            lav_model_get_parameters(object@Model,
              GLIST = lav_model_x2GLIST(object@Model, boot.x),
              type = "user", extra = TRUE
            )
          bias.est <- (boot.est - tmp.list$est)
          ci <- (tmp.list$est - bias.est) + tmp.list$se %o% fac
        } else if (boot.ci.type == "basic") {
          ci <- cbind(tmp.list$est, tmp.list$est)
          alpha <- (1 + c(level, -level)) / 2

          # free.idx only
          qq <- apply(tmp.boot, 2, norm.inter, alpha)
          free.idx <- which(object@ParTable$free &
            !duplicated(object@ParTable$free))
          ci[free.idx, ] <- 2 * ci[free.idx, ] - t(qq[c(3, 4), ])

          # def.idx
          def.idx <- which(object@ParTable$op == ":=")
          if (length(def.idx) > 0L) {
            boot.def <- apply(tmp.boot, 1, object@Model@def.function)
            if (length(def.idx) == 1L) {
              boot.def <- as.matrix(boot.def)
            } else {
              boot.def <- t(boot.def)
            }
            qq <- apply(boot.def, 2, norm.inter, alpha)
            ci[def.idx, ] <- 2 * ci[def.idx, ] - t(qq[c(3, 4), ])
          }

          # TODO: add cin/ceq?
        } else if (boot.ci.type == "perc") {
          ci <- cbind(tmp.list$est, tmp.list$est)
          alpha <- (1 + c(-level, level)) / 2

          # free.idx only
          qq <- apply(tmp.boot, 2, norm.inter, alpha)
          free.idx <- which(object@ParTable$free &
            !duplicated(object@ParTable$free))
          ci[free.idx, ] <- t(qq[c(3, 4), ])

          # def.idx
          def.idx <- which(object@ParTable$op == ":=")
          if (length(def.idx) > 0L) {
            boot.def <- apply(tmp.boot, 1, object@Model@def.function)
            if (length(def.idx) == 1L) {
              boot.def <- as.matrix(boot.def)
            } else {
              boot.def <- t(boot.def)
            }
            qq <- apply(boot.def, 2, norm.inter, alpha)
            def.idx <- which(object@ParTable$op == ":=")
            ci[def.idx, ] <- t(qq[c(3, 4), ])
          }

          # TODO:  add cin/ceq?
        } else if (boot.ci.type == "bca.simple") {
          # no adjustment for scale!! only bias!!
          alpha <- (1 + c(-level, level)) / 2
          zalpha <- qnorm(alpha)
          ci <- cbind(tmp.list$est, tmp.list$est)

          # free.idx only
          free.idx <- which(object@ParTable$free &
            !duplicated(object@ParTable$free))
          x <- tmp.list$est[free.idx]
          for (i in seq_along(free.idx)) {
            t <- tmp.boot[, i]
            t <- t[is.finite(t)]
            t0 <- x[i]
            # check if we have variance (perhaps constrained to 0?)
            # new in 0.6-3
            if (var(t) == 0) {
              next
            }
            w <- qnorm(sum(t < t0) / length(t))
            a <- 0.0 #### !!! ####
            adj.alpha <- pnorm(w + (w + zalpha) / (1 - a * (w + zalpha)))
            qq <- norm.inter(t, adj.alpha)
            ci[free.idx[i], ] <- qq[, 2]
          }

          # def.idx
          def.idx <- which(object@ParTable$op == ":=")
          if (length(def.idx) > 0L) {
            x.def <- object@Model@def.function(x)
            boot.def <- apply(tmp.boot, 1, object@Model@def.function)
            if (length(def.idx) == 1L) {
              boot.def <- as.matrix(boot.def)
            } else {
              boot.def <- t(boot.def)
            }
            for (i in seq_along(def.idx)) {
              t <- boot.def[, i]
              t <- t[is.finite(t)]
              t0 <- x.def[i]
              w <- qnorm(sum(t < t0) / length(t))
              a <- 0.0 #### !!! ####
              adj.alpha <- pnorm(w + (w + zalpha) / (1 - a * (w + zalpha)))
              qq <- norm.inter(t, adj.alpha)
              ci[def.idx[i], ] <- qq[, 2]
            }
          }

          # TODO:
          # - add cin/ceq
        } else if (boot.ci.type == "bca") { # new in 0.6-12
          # we assume that the 'ordinary' (nonparametric) was used

          lavoptions <- object@Options
          ngroups <- object@Data@ngroups
          nobs <- object@SampleStats@nobs
          ntotal <- object@SampleStats@ntotal

          # we need enough bootstrap runs
          if (nrow(tmp.boot) < ntotal) {
            lav_msg_stop(gettextf(
            "BCa confidence intervals require more (successful) bootstrap runs
            (%1$s) than the number of observations (%2$s).",
            nrow(tmp.boot), ntotal))
          }

          # does not work with sampling weights (yet)
          if (!is.null(object@Data@weights[[1]])) {
            lav_msg_stop(
              gettext("BCa confidence intervals not available in
                      the presence of sampling weights."))
          }

          # check if we have a seed
          if (is.null(bootstrap.seed)) {
            lav_msg_stop(gettext("seed not available in tmp.boot object."))
          }

          # compute 'X' matrix with frequency indices (to compute
          # the empirical influence values using regression)
          tmp.freq <- lav_utils_bootstrap_indices(
            R = lavoptions$bootstrap,
            nobs = nobs, parallel = lavoptions$parallel[1],
            ncpus = lavoptions$ncpus, cl = lavoptions[["cl"]],
            iseed = bootstrap.seed, return.freq = TRUE,
            merge.groups = TRUE
          )
          if (length(error.idx) > 0L) {
            tmp.freq <- tmp.freq[-error.idx, , drop = FALSE]
          }
          stopifnot(nrow(tmp.freq) == nrow(tmp.boot))

          # compute empirical influence values (using regression)
          # remove first column per group
          first.idx <- sapply(object@Data@case.idx, "[[", 1L)
          tmp.lm <- lm.fit(x = cbind(1, tmp.freq[, -first.idx]), y = tmp.boot)
          tmp.beta <- unname(tmp.lm$coefficients)[-1, , drop = FALSE]
          tmp.ll <- rbind(0, tmp.beta)

          # compute 'a' for all parameters at once
          tmp.aa <- apply(tmp.ll, 2L, function(x) {
            tmp.l <- x - mean(x)
            sum(tmp.l^3) / (6 * sum(tmp.l^2)^1.5)
          })

          # adjustment for both bias AND scale
          alpha <- (1 + c(-level, level)) / 2
          zalpha <- qnorm(alpha)
          ci <- cbind(tmp.list$est, tmp.list$est)

          # free.idx only
          free.idx <- which(object@ParTable$free &
            !duplicated(object@ParTable$free))
          stopifnot(length(free.idx) == ncol(tmp.boot))
          x <- tmp.list$est[free.idx]
          for (i in seq_along(free.idx)) {
            t <- tmp.boot[, i]
            t <- t[is.finite(t)]
            t0 <- x[i]
            # check if we have variance (perhaps constrained to 0?)
            # new in 0.6-3
            if (var(t) == 0) {
              next
            }
            w <- qnorm(sum(t < t0) / length(t))
            a <- tmp.aa[i]
            adj.alpha <- pnorm(w + (w + zalpha) / (1 - a * (w + zalpha)))
            qq <- norm.inter(t, adj.alpha)
            ci[free.idx[i], ] <- qq[, 2]
          }

          # def.idx
          def.idx <- which(object@ParTable$op == ":=")
          if (length(def.idx) > 0L) {
            x.def <- object@Model@def.function(x)
            boot.def <- apply(tmp.boot, 1, object@Model@def.function)
            if (length(def.idx) == 1L) {
              boot.def <- as.matrix(boot.def)
            } else {
              boot.def <- t(boot.def)
            }

            # recompute empirical influence values
            tmp.lm <- lm.fit(x = cbind(1, tmp.freq[, -1]), y = boot.def)
            tmp.beta <- unname(tmp.lm$coefficients)[-1, , drop = FALSE]
            tmp.ll <- rbind(0, tmp.beta)

            # compute 'a' values for all def.idx parameters
            tmp.aa <- apply(tmp.ll, 2L, function(x) {
              tmp.l <- x - mean(x)
              sum(tmp.l^3) / (6 * sum(tmp.l^2)^1.5)
            })

            # compute bca ci
            for (i in seq_along(def.idx)) {
              t <- boot.def[, i]
              t <- t[is.finite(t)]
              t0 <- x.def[i]
              w <- qnorm(sum(t < t0) / length(t))
              a <- tmp.aa[i]
              adj.alpha <- pnorm(w + (w + zalpha) / (1 - a * (w + zalpha)))
              qq <- norm.inter(t, adj.alpha)
              ci[def.idx[i], ] <- qq[, 2]
            }
          }

          # TODO:
          # - add cin/ceq
        }
      }

      tmp.list$ci.lower <- ci[, 1]
      tmp.list$ci.upper <- ci[, 2]
    }

    # standardized estimates?
    # 28 March 2024: TDJ adds option to select specific types
    if (is.logical(standardized)) {
      if (standardized) {
        standardized <- c("std.lv", "std.all")
        if (length(lavNames(object, "ov.x")) && object@Options$fixed.x) {
          standardized <- c(standardized, "std.nox")
        }
      } else {
        standardized <- character(0)
      } # corresponds to standardized=FALSE
    } else {
      # !is.logical(standardized)
      standardized <- tolower(as.character(standardized))
      if ("std.nox" %in% standardized) {
        # sanity checks
        if (length(lavNames(object, "ov.x")) == 0) {
          lav_msg_note(gettext(
            "`std.nox' unavailable without fixed exogenous predictors"))
          standardized <- setdiff(standardized, "std.nox")
        } else if (!object@Options$fixed.x) {
          lav_msg_note(gettext("`std.nox' unavailable when fixed.x = FALSE"))
          standardized <- setdiff(standardized, "std.nox")
        }
      }
    }
    # Then add each requested type
    # (original source code, but now independently conditional)
    if ("std.lv" %in% standardized) {
      tmp.list$std.lv <- lav_standardize_lv(object, cov.std = cov.std)
    }
    if ("std.all" %in% standardized) {
      tmp.list$std.all <- lav_standardize_all(object,
        est.std = tmp.list$est.std,
        cov.std = cov.std
      )
    }
    if ("std.nox" %in% standardized) {
      tmp.list$std.nox <- lav_standardize_all_nox(object,
        est.std = tmp.list$est.std,
        cov.std = cov.std
      )
    }

    # rsquare?
    if (rsquare) {
      r2 <- lavTech(object, "rsquare", add.labels = TRUE)
      tmp.names <- unlist(lapply(r2, names))
      nel <- length(tmp.names)
      if (nel == 0L) {
        lav_msg_warn(
          gettext("rsquare = TRUE, but there are no dependent variables"))
      } else {
        if (lav_partable_nlevels(tmp.list) == 1L) {
          block <- rep(seq_along(r2), sapply(r2, length))
          first.block.idx <- which(!duplicated(tmp.list$block) &
            tmp.list$block > 0L)
          gval <- tmp.list$group[first.block.idx]
          if (length(gval) > 0L) {
            group <- rep(gval, sapply(r2, length))
          } else {
            # single block, single group
            group <- rep(1L, length(block))
          }
          r2 <- data.frame(
            lhs = tmp.names, op = rep("r2", nel),
            rhs = tmp.names, block = block, group = group,
            est = unlist(r2), stringsAsFactors = FALSE
          )
        } else {
          # add level column
          block <- rep(seq_along(r2), sapply(r2, length))
          first.block.idx <- which(!duplicated(tmp.list$block) &
            tmp.list$block > 0L)
          # always at least two blocks
          gval <- tmp.list$group[first.block.idx]
          group <- rep(gval, sapply(r2, length))
          lval <- tmp.list$level[first.block.idx]
          level <- rep(lval, sapply(r2, length))
          r2 <- data.frame(
            lhs = tmp.names, op = rep("r2", nel),
            rhs = tmp.names,
            block = block, group = group,
            level = level,
            est = unlist(r2), stringsAsFactors = FALSE
          )
        }
        # add step column if needed
        if (!is.null(tmp.list$step)) {
          r2$step <- 2L # per default
          # simplification: we assume that only the
          # observed indicators of latent variables are step 1
          ov.ind <- unlist(object@pta$vnames$ov.ind)
          step1.idx <- which(r2$lhs %in% ov.ind)
          r2$step[step1.idx] <- 1L
        }
        tmp.list <- lav_partable_merge(pt1 = tmp.list, pt2 = r2, warn = FALSE)
      }
    }

    # fractional missing information (if estimator="fiml")
    if (fmi) {
      se.orig <- tmp.list$se

      # new in 0.6-6, use 'EM' based (unstructured) sample statistics
      # otherwise, it would be as if we use expected info, while the
      # original use observed, producing crazy results
      if (object@Data@ngroups > 1L) {
        em.cov <- lapply(lavInspect(object, "sampstat.h1"), "[[", "cov")
        em.mean <- lapply(lavInspect(object, "sampstat.h1"), "[[", "mean")
      } else {
        em.cov <- lavInspect(object, "sampstat.h1")$cov
        em.mean <- lavInspect(object, "sampstat.h1")$mean
      }

      tmp.pt <- parTable(object)
      tmp.pt$ustart <- tmp.pt$est
      tmp.pt$start <- tmp.pt$est <- NULL

      this.options <- object@Options
      if (!is.null(fmi.options) && is.list(fmi.options)) {
        # modify original options
        this.options <- modifyList(this.options, fmi.options)
      }
      # override
      this.options$optim.method <- "none"
      this.options$sample.cov.rescale <- FALSE
      this.options$check.gradient <- FALSE
      this.options$baseline <- FALSE
      this.options$h1 <- FALSE
      this.options$test <- FALSE

      fit.complete <- lavaan(
        model = tmp.pt,
        sample.cov = em.cov,
        sample.mean = em.mean,
        sample.nobs = lavInspect(object, "nobs"),
        slotOptions = this.options
      )

      se.comp <- parameterEstimates(fit.complete,
        ci = FALSE, fmi = FALSE,
        zstat = FALSE, pvalue = FALSE, remove.system.eq = FALSE,
        remove.eq = FALSE, remove.ineq = FALSE,
        remove.def = FALSE, remove.nonfree = FALSE, remove.unused = FALSE,
        rsquare = rsquare, add.attributes = FALSE
      )$se

      se.comp <- ifelse(se.comp == 0.0, as.numeric(NA), se.comp)
      tmp.list$fmi <- 1 - (se.comp * se.comp) / (se.orig * se.orig)
    }

    # if single level, remove level column
    if (object@Data@nlevels == 1L) tmp.list$level <- NULL

    # if single group, remove group column
    if (object@Data@ngroups == 1L) tmp.list$group <- NULL

    # if single everything, remove block column
    if (object@Data@nlevels == 1L &&
      object@Data@ngroups == 1L) {
      tmp.list$block <- NULL
    }

    # if no user-defined labels, remove label column
    if (sum(nchar(object@ParTable$label)) == 0L) {
      tmp.list$label <- NULL
    }

    # remove non-free parameters? (but keep ==, >, < and :=)
    if (remove.nonfree) {
      nonfree.idx <- which(tmp.list$free == 0L &
        !tmp.list$op %in% c("==", ">", "<", ":="))
      if (length(nonfree.idx) > 0L) {
        tmp.list <- tmp.list[-nonfree.idx, ]
      }
    }

    # remove 'unused' parameters
    # these are parameters that are automatically added (user == 0),
    # but with their final (est) values fixed to their default values
    # (typically 1 or 0).
    # currently only intercepts and scaling-factors (for now)
    # should we also remove fixed-to-1 variances? (parameterization = theta)?
    if (remove.unused) {
      # intercepts
      int.idx <- which(tmp.list$op == "~1" &
        tmp.list$user == 0L &
        tmp.list$free == 0L &
        tmp.list$est == 0)
      if (length(int.idx) > 0L) {
        tmp.list <- tmp.list[-int.idx, ]
      }

      # scaling factors
      scaling.idx <- which(tmp.list$op == "~*~" &
        tmp.list$user == 0L &
        tmp.list$free == 0L &
        tmp.list$est == 1)
      if (length(scaling.idx) > 0L) {
        tmp.list <- tmp.list[-scaling.idx, ]
      }
    }


    # remove 'free' column
    tmp.list$free <- NULL

    # remove == rows?
    if (remove.eq) {
      eq.idx <- which(tmp.list$op == "==" & tmp.list$user == 1L)
      if (length(eq.idx) > 0L) {
        tmp.list <- tmp.list[-eq.idx, ]
      }
    }
    if (remove.system.eq) {
      eq.idx <- which(tmp.list$op == "==" & tmp.list$user != 1L)
      if (length(eq.idx) > 0L) {
        tmp.list <- tmp.list[-eq.idx, ]
      }
    }
    # remove <> rows?
    if (remove.ineq) {
      ineq.idx <- which(tmp.list$op %in% c("<", ">"))
      if (length(ineq.idx) > 0L) {
        tmp.list <- tmp.list[-ineq.idx, ]
      }
    }
    # remove := rows?
    if (remove.def) {
      def.idx <- which(tmp.list$op == ":=")
      if (length(def.idx) > 0L) {
        tmp.list <- tmp.list[-def.idx, ]
      }
    }

    # remove step 1 rows?
    if (remove.step1 && !is.null(tmp.list$step)) {
      step1.idx <- which(tmp.list$step == 1L)
      if (length(step1.idx) > 0L) {
        tmp.list <- tmp.list[-step1.idx, ]
      }
      # remove step column
      tmp.list$step <- NULL
    }

    # remove attribute for data order
    attr(tmp.list, "ovda") <- NULL

    # remove tmp.list$user
    tmp.list$user <- NULL

    if (output == "text") {
      class(tmp.list) <- c(
        "lavaan.parameterEstimates", "lavaan.data.frame",
        "data.frame"
      )
      if (header) {
        attr(tmp.list, "categorical") <- object@Model@categorical
        attr(tmp.list, "parameterization") <- object@Model@parameterization
        attr(tmp.list, "information") <- object@Options$information[1]
        attr(tmp.list, "information.meat") <- object@Options$information.meat
        attr(tmp.list, "se") <- object@Options$se
        attr(tmp.list, "group.label") <- object@Data@group.label
        attr(tmp.list, "level.label") <- object@Data@level.label
        attr(tmp.list, "bootstrap") <- object@Options$bootstrap
        attr(tmp.list, "bootstrap.successful") <- bootstrap.successful
        attr(tmp.list, "missing") <- object@Options$missing
        attr(tmp.list, "observed.information") <-
          object@Options$observed.information[1]
        attr(tmp.list, "h1.information") <- object@Options$h1.information[1]
        attr(tmp.list, "h1.information.meat") <-
          object@Options$h1.information.meat
        attr(tmp.list, "header") <- header
        # FIXME: add more!!
      }
    } else {
      tmp.list$exo <- NULL
      tmp.list$lower <- tmp.list$upper <- NULL
      class(tmp.list) <- c("lavaan.data.frame", "data.frame")
    }

    tmp.list
  }

parameterTable <- parametertable <- parTable <- partable <- # nolint
  function(object) {
    # convert to data.frame
    out <- as.data.frame(object@ParTable, stringsAsFactors = FALSE)

    class(out) <- c("lavaan.data.frame", "data.frame")
    out
  }

varTable <- vartable <- function(object, ov.names = names(object), # nolint
                                 ov.names.x = NULL,
                                 ordered = NULL, factor = NULL,
                                 as.data.frame. = TRUE) { # nolint

  if (inherits(object, "lavaan")) {
    tmp.var <- object@Data@ov
  } else if (inherits(object, "lavData")) {
    tmp.var <- object@ov
  } else if (inherits(object, "data.frame")) {
    tmp.var <- lav_dataframe_vartable(
      frame = object, ov.names = ov.names,
      ov.names.x = ov.names.x,
      ordered = ordered, factor = factor,
      as.data.frame. = FALSE
    )
  } else {
    lav_msg_stop(gettext("object must of class lavaan or a data.frame"))
  }

  if (as.data.frame.) {
    tmp.var <- as.data.frame(tmp.var,
      stringsAsFactors = FALSE,
      row.names = seq_along(tmp.var$name)
    )
    class(tmp.var) <- c("lavaan.data.frame", "data.frame")
  }

  tmp.var
}


setMethod(
  "fitted.values", "lavaan",
  function(object, type = "moments", labels = TRUE) {
    # lowercase type
    type <- tolower(type)

    # catch type="casewise"
    if (type %in% c("casewise", "case", "obs", "observations", "ov")) {
      return(lavPredict(object, type = "ov", label = labels))
    }

    lav_object_inspect_implied(object,
      add.labels = labels, add.class = TRUE,
      drop.list.single.group = TRUE
    )
  }
)


setMethod(
  "fitted", "lavaan",
  function(object, type = "moments", labels = TRUE) {
    fitted.values(object, type = type, labels = labels)
  }
)


setMethod(
  "vcov", "lavaan",
  function(object, type = "free", labels = TRUE, remove.duplicated = FALSE) {
    # check for convergence first!
    if (object@optim$npar > 0L && !object@optim$converged) {
      lav_msg_stop(gettext("model did not converge"))
    }

    if (object@Options$se == "none") {
      lav_msg_stop(gettext("vcov not available if se=\"none\""))
    }

    if (type == "user" || type == "joint" || type == "all" || type == "full" ||
      type == "complete") {
      if (remove.duplicated) {
        lav_msg_stop(gettext(
          "argument \"remove.duplicated\" not supported if type = \"user\""
        ))
      }
      tmp.varcov <- lav_object_inspect_vcov_def(object,
        joint = TRUE,
        add.labels = labels,
        add.class = TRUE
      )
    } else if (type == "free") {
      tmp.varcov <- lav_object_inspect_vcov(object,
        add.labels = labels,
        add.class = TRUE,
        remove.duplicated = remove.duplicated
      )
    } else {
      lav_msg_stop(gettext("type argument should be \"user\" or \"free\""))
    }

    tmp.varcov
  }
)


# logLik (so that we can use the default AIC/BIC functions from stats4(
setMethod(
  "logLik", "lavaan",
  function(object, ...) {
    if (object@Options$estimator != "ML") {
      lav_msg_warn(gettext("logLik only available if estimator is ML"))
    }
    if (object@optim$npar > 0L && !object@optim$converged) {
      lav_msg_warn(gettext("model did not converge"))
    }

    # new in 0.6-1: we use the @loglik slot (instead of fitMeasures)
    if (.hasSlot(object, "loglik")) {
      tmp.logl <- object@loglik
    } else {
      tmp.logl <- lav_model_loglik(
        lavdata = object@Data,
        lavsamplestats = object@SampleStats,
        lavimplied = object@implied,
        lavmodel = object@Model,
        lavoptions = object@Options
      )
    }

    logl <- tmp.logl$loglik
    attr(logl, "df") <- tmp.logl$npar ### note: must be npar, not df!!
    attr(logl, "nobs") <- tmp.logl$ntotal
    class(logl) <- "logLik"
    logl
  }
)

# nobs
if (!exists("nobs", envir = asNamespace("stats4"))) {
  setGeneric("nobs", function(object, ...) standardGeneric("nobs"))
}
setMethod(
  "nobs", signature(object = "lavaan"),
  function(object, ...) {
    object@SampleStats@ntotal
  }
)

# see: src/library/stats/R/update.R
setMethod(
  "update", signature(object = "lavaan"),
  function(object, model, add, ..., evaluate = TRUE) {
    call <- object@call
    if (is.null(call)) {
      lav_msg_stop(gettext("need an object with call slot"))
    }

    extras <- match.call(expand.dots = FALSE)$...

    if (!missing(model)) {
      # call$formula <- update.formula(formula(object), formula.)
      call$model <- model
    } else if (exists(as.character(call$model))) {
      call$model <- eval(call$model, parent.frame())
    } else if (is.character(call$model)) {
      ## do nothing
      ## call$model <- call$model
    } else {
      call$model <- parTable(object)
      call$model$est <- NULL
      call$model$se <- NULL
    }
    if (!is.null(call$slotParTable) && is.list(call$model)) {
      call$slotParTable <- call$model
    }

    if (length(extras) > 0) {
      ## check for call$slotOptions conflicts
      if (!is.null(call$slotOptions)) {
        same.names <- intersect(names(lavOptions()), names(extras))
        for (i in same.names) {
          call$slotOptions[[i]] <- extras[[i]]
          extras[i] <- NULL # not needed if they are in slotOptions
        }
      }
      existing <- !is.na(match(names(extras), names(call)))
      for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
      if (any(!existing)) {
        call <- c(as.list(call), extras[!existing])
        call <- as.call(call)
      }
    }

    if (missing(add) && !evaluate) {
      return(call)
    }
    ## for any of the other 3 scenarios, we need the updated fit

    ## Check if "add" and "model" are both strings; combine them
    if (missing(add)) {
      add.allready.in.partable <- TRUE # because nothing to add
    } else {
      if (is.character(add) && is.character(call$model)) {
        call$model <- c(call$model, add)
        add.allready.in.partable <- TRUE
      } else {
        add.allready.in.partable <- FALSE
      }
    }
    newfit <- eval(call, parent.frame())
    if (add.allready.in.partable && evaluate) {
      return(newfit)
    }

    ## only remaining situations: "add" exists, but either "add" or "model"
    ## is a parameter table, so update the parameter table in the call
    if (!(mode(add) %in% c("list", "character"))) {
      lav_msg_stop(
        gettext("'add' argument must be model syntax or parameter table.
                See ?lavaanify help page.")
      )
    }
    tmp.pt <- lav_object_extended(newfit, add = add)@ParTable
    tmp.pt$user <- NULL # get rid of "10" category used in lavTestScore()
    ## group == 0L in new rows
    tmp.pt$group[tmp.pt$group == 0L] <- tmp.pt$block[tmp.pt$group == 0L]
    # tmp.pt$plabel == "" in new rows.  Consequences?
    tmp.pt$est <- NULL
    tmp.pt$se <- NULL
    call$model <- tmp.pt

    if (evaluate) {
      eval(call, parent.frame())
    } else {
      call
    }
  }
)




setMethod(
  "anova", signature(object = "lavaan"),
  function(object, ...) {
    # NOTE: if we add additional arguments, it is not the same generic
    # anova() function anymore, and match.call will be screwed up

    # NOTE: we need to extract the names of the models from match.call here,
    #       otherwise, we loose them in the call stack

    mcall <- match.call(expand.dots = TRUE)
    dots <- list(...)

    # catch SB.classic and SB.H0
    # SB.classic <- TRUE; SB.H0 <- FALSE

    # arg.names <- names(dots)
    # arg.idx <- which(nchar(arg.names) > 0L)
    # if(length(arg.idx) > 0L) {
    #    if(!is.null(dots$SB.classic))
    #        SB.classic <- dots$SB.classic
    #    if(!is.null(dots$SB.H0))
    #        SB.H0 <- dots$SB.H0
    #    dots <- dots[-arg.idx]
    # }

    modp <- if (length(dots)) {
      sapply(dots, inherits, "lavaan")
    } else {
      logical(0)
    }
    tmp.names <- sapply(as.list(mcall)[c(FALSE, TRUE, modp)], deparse)

    # use do.call to handle changed dots
    # ans <- do.call("lavTestLRT", c(list(object = object,
    #               SB.classic = SB.classic, SB.H0 = SB.H0,
    #               model.names = tmp.names), dots))

    # ans
    lavTestLRT(object = object, ..., model.names = tmp.names)
  }
)
