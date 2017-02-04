#
# initial version: YR 25/03/2009

short.summary <- function(object) {

    # catch FAKE run
    FAKE <- FALSE
    if(object@Options$optim.method == "none") {
        FAKE <- TRUE
    }

    # Convergence or not?
    if(FAKE) {
        cat(sprintf("lavaan (%s) -- DRY RUN with 0 iterations\n",
                    packageDescription("lavaan", fields="Version")))
    } else if(object@optim$iterations > 0) {
        if(object@optim$converged) {
	    cat(sprintf("lavaan (%s) converged normally after %3i iterations\n",
                    packageDescription("lavaan", fields="Version"),
                    object@optim$iterations))
        } else {
            cat(sprintf("** WARNING ** lavaan (%s) did NOT converge after %i iterations\n", 
                packageDescription("lavaan", fields="Version"),
                object@optim$iterations))
            cat("** WARNING ** Estimates below are most likely unreliable\n")
        }
    } else {
        cat(sprintf("** WARNING ** lavaan (%s) model has NOT been fitted\n",
                    packageDescription("lavaan", fields="Version")))
        cat("** WARNING ** Estimates below are simply the starting values\n")
    }
    cat("\n")

    # number of free parameters
    #t0.txt <- sprintf("  %-40s", "Number of free parameters")
    #t1.txt <- sprintf("  %10i", object@optim$npar)
    #t2.txt <- ""
    #cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
    #cat("\n")
   
    # listwise deletion?
    listwise <- FALSE
    for(g in 1:object@Data@ngroups) {
       if(object@Data@nobs[[1L]] != object@Data@norig[[1L]]) {
           listwise <- TRUE
           break
       }
    }


    if(object@Data@ngroups == 1L) {
        if(listwise) {
            cat(sprintf("  %-40s", ""), sprintf("  %10s", "Used"), 
                                        sprintf("  %10s", "Total"),
                "\n", sep="")
        }
        t0.txt <- sprintf("  %-40s", "Number of observations")
        t1.txt <- sprintf("  %10i", object@Data@nobs[[1L]])
        t2.txt <- ifelse(listwise,
                  sprintf("  %10i", object@Data@norig[[1L]]), "")
        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
    } else {
        if(listwise) {
            cat(sprintf("  %-40s", ""), sprintf("  %10s", "Used"),  
                                        sprintf("  %10s", "Total"),
                "\n", sep="")
        }
        t0.txt <- sprintf("  %-40s", "Number of observations per group")
        cat(t0.txt, "\n")
        for(g in 1:object@Data@ngroups) {
            t.txt <- sprintf("  %-40s  %10i", object@Data@group.label[[g]],
                                              object@Data@nobs[[g]])
            t2.txt <- ifelse(listwise,
                      sprintf("  %10i", object@Data@norig[[g]]), "")
            cat(t.txt, t2.txt, "\n", sep="")
        }
    }
    cat("\n")

    # missing patterns?
    if(object@SampleStats@missing.flag) {
        if(object@Data@ngroups == 1L) {
            t0.txt <- sprintf("  %-40s", "Number of missing patterns")
            t1.txt <- sprintf("  %10i", 
                              object@Data@Mp[[1L]]$npatterns)
            cat(t0.txt, t1.txt, "\n\n", sep="")
        } else {
            t0.txt <- sprintf("  %-40s", "Number of missing patterns per group")
            cat(t0.txt, "\n")
            for(g in 1:object@Data@ngroups) {
                t.txt <- sprintf("  %-40s  %10i", object@Data@group.label[[g]],
                                 object@Data@Mp[[g]]$npatterns)
                cat(t.txt, "\n", sep="")
            }
            cat("\n")
        }
    }

    # Print Chi-square value for the user-specified (full/h0) model

    # robust/scaled statistics?
    if(object@Options$test %in% c("satorra.bentler", "yuan.bentler",
                                  "mean.var.adjusted",
                                  "scaled.shifted") &&
       length(object@test) > 1L) {
        scaled <- TRUE
        if(object@Options$test == "scaled.shifted")
            shifted <- TRUE
        else
            shifted <- FALSE
    } else {
        scaled <- FALSE
        shifted <- FALSE
    }

    # 0. heading
    #h.txt <- sprintf("\nChi-square test user model (h0)",
    #                 object@Options$estimator)
    t0.txt <- sprintf("  %-40s", "Estimator")
    t1.txt <- sprintf("  %10s", object@Options$estimator)
    t2.txt <- ifelse(scaled, 
              sprintf("  %10s", "Robust"), "")
    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

    # check if test == "none"
    if(object@Options$test != "none" && object@Options$estimator != "MML") {

        # 1. chi-square values
        t0.txt <- sprintf("  %-40s", "Minimum Function Test Statistic")  
        t1.txt <- sprintf("  %10.3f", object@test[[1]]$stat)
        t2.txt <- ifelse(scaled, 
                  sprintf("  %10.3f", object@test[[2]]$stat), "")
        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

        # 2. degrees of freedom
        t0.txt <- sprintf("  %-40s", "Degrees of freedom")
        t1.txt <- sprintf("  %10i",   object@test[[1]]$df)
        t2.txt <- ifelse(scaled, 
                         ifelse(round(object@test[[2]]$df) == 
                                object@test[[2]]$df,
                                sprintf("  %10i",   object@test[[2]]$df),
                                sprintf("  %10.3f", object@test[[2]]$df)),
                         "")
        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

        # 3. P-value
        if(is.na(object@test[[1]]$df)) {
            t0.txt <- sprintf("  %-40s", "P-value")
            t1.txt <- sprintf("  %10.3f", object@test[[1]]$pvalue)
            t2.txt <- ifelse(scaled,
                      sprintf("  %10.3f", object@test[[2]]$pvalue), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        } else if(object@test[[1]]$df > 0) {
            if(object@test[[1]]$refdistr == "chisq") {
                t0.txt <- sprintf("  %-40s", "P-value (Chi-square)")
            } else if(length(object@test) == 1L &&
                      object@test[[1]]$refdistr == "unknown") {
                t0.txt <- sprintf("  %-40s", "P-value (Unknown)")
            } else {
                t0.txt <- sprintf("  %-40s", "P-value")
            }
            t1.txt <- sprintf("  %10.3f", object@test[[1]]$pvalue)
            t2.txt <- ifelse(scaled,
                      sprintf("  %10.3f", object@test[[2]]$pvalue), "")
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        } else {
            # FIXME: should we do this? To warn that exact 0.0 was not obtained?
            if(object@optim$fx > 0) {
                t0.txt <- sprintf("  %-35s", "Minimum Function Value")
                t1.txt <- sprintf("  %15.13f", object@optim$fx)
                t2.txt <- ""
                cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
            }
        }

        # 3b. Do we have a Bollen-Stine p-value?
        if(object@Options$test == "bollen.stine") {
            t0.txt <- sprintf("  %-40s", "P-value (Bollen-Stine Bootstrap)")
            t1.txt <- sprintf("  %10.3f", object@test[[2]]$pvalue)
            cat(t0.txt, t1.txt, "\n", sep="")
        }

        # 4. Scaling correction factor
        if(scaled) {
            t0.txt <- sprintf("  %-40s", "Scaling correction factor")
            t1.txt <- sprintf("  %10s", "")
            t2.txt <- sprintf("  %10.3f", object@test[[2]]$scaling.factor)
            cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
            if(object@Options$test == "yuan.bentler") {
                if(object@Options$mimic == "Mplus") {
                    cat("    for the Yuan-Bentler correction (Mplus variant)\n")
                } else {
                    cat("    for the Yuan-Bentler correction\n")
                }
            } else if(object@Options$test == "satorra.bentler") {
                if(object@Options$mimic == "Mplus" &&
                   object@Options$estimator == "ML") {
                    cat("    for the Satorra-Bentler correction (Mplus variant)\n")
                } else if(object@Options$mimic == "Mplus" &&
                          object@Options$estimator == "DWLS") {
                    cat("    for the Satorra-Bentler correction (WLSM)\n")
                } else if(object@Options$mimic == "Mplus" &&
                          object@Options$estimator == "ULS") {
                    cat("    for the Satorra-Bentler correction (ULSM)\n")
                } else {
                    cat("    for the Satorra-Bentler correction\n")
                }
            } else if(object@Options$test == "mean.var.adjusted") {
                if(object@Options$mimic == "Mplus" &&
                   object@Options$estimator == "ML") {
                    cat("    for the mean and variance adjusted correction (MLMV)\n")
                } else if(object@Options$mimic == "Mplus" &&
                          object@Options$estimator == "DWLS") {
                    cat("    for the mean and variance adjusted correction (WLSMV)\n")
                } else if(object@Options$mimic == "Mplus" &&
                          object@Options$estimator == "ULS") {
                    cat("    for the mean and variance adjusted correction (ULSMV)\n")
                } else {
                    cat("    for the mean and variance adjusted correction\n")
                }
            }
        }

        # 4b. Shift parameter?
        if(shifted) {
            if(object@Data@ngroups == 1L) {
                t0.txt <- sprintf("  %-40s", "Shift parameter")
                t1.txt <- sprintf("  %10s", "")
                t2.txt <- sprintf("  %10.3f", 
                                  object@test[[2]]$shift.parameter)
                cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
            } else { # multiple groups, multiple shift values!
                cat("  Shift parameter for each group:\n")
                for(g in 1:object@Data@ngroups) {
                    t0.txt <- sprintf("    %-38s", object@Data@group.label[[g]])
                    t1.txt <- sprintf("  %10s", "")
                    t2.txt <- sprintf("  %10.3f",
                                     object@test[[2]]$shift.parameter[g])
                    cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
                }
            }
            if(object@Options$mimic == "Mplus" &&
               object@Options$estimator == "DWLS") {
                cat("    for simple second-order correction (WLSMV)\n")
            } else {
                cat("    for simple second-order correction (Mplus variant)\n")
            }
        }

        if(object@Data@ngroups > 1L) {
            cat("\n")
            cat("Chi-square for each group:\n\n")
            for(g in 1:object@Data@ngroups) {
                t0.txt <- sprintf("  %-40s", object@Data@group.label[[g]])
                t1.txt <- sprintf("  %10.3f", object@test[[1]]$stat.group[g])
                t2.txt <- ifelse(scaled, sprintf("  %10.3f", 
                                 object@test[[2]]$stat.group[g]), "")
                cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
            }
        } 
    } # test != none

    if(object@Options$estimator == "MML") {
        fm <- fitMeasures(object, c("logl", "npar", "aic", "bic", "bic2"))
        print.fit.measures(fm)
    }

    #cat("\n")
}

setMethod("show", "lavaan",
function(object) {

    # show only basic information
    short.summary(object)

})

setMethod("summary", "lavaan",
function(object, header       = TRUE,
                 fit.measures = FALSE, 
                 estimates    = TRUE,
                 ci           = FALSE, 
                 fmi          = FALSE,
                 standardized = FALSE,
                 rsquare      = FALSE, 
                 std.nox      = FALSE, 
                 modindices   = FALSE, 
                 nd = 3L) {

    if(std.nox) standardized <- TRUE

    # print the 'short' summary
    if(header) {
        short.summary(object)
    }

    # only if requested, the fit measures
    if(fit.measures) {
        if(object@Options$test == "none") {
            warning("lavaan WARNING: fit measures not available if test = \"none\"\n\n")
        } else if(object@optim$npar > 0L && !object@optim$converged) {
            warning("lavaan WARNING: fit measures not available if model did not converge\n\n")
        } else {
            print.fit.measures( fitMeasures(object, fit.measures="default") )
        }
    }

    if(estimates) {
        PE <- parameterEstimates(object, ci = ci, standardized = standardized,
                                 rsquare = rsquare, fmi = fmi,
                                 remove.eq = FALSE, remove.system.eq = TRUE,
                                 remove.ineq = FALSE, remove.def = FALSE, 
                                 add.attributes = TRUE)
        if(standardized && std.nox) {
            PE$std.all <- PE$std.nox
        }
        print(PE, nd = nd)
    }

    # modification indices?
    if(modindices) {
        cat("Modification Indices:\n\n")
        print( modificationIndices(object, standardized=TRUE) )
    }
})

# old summary (<0.5-19)
summary2 <- function(object, estimates=TRUE, fit.measures=FALSE, 
                     standardized=FALSE, 
                     rsquare=FALSE, std.nox=FALSE, modindices=FALSE) {

    if(std.nox) standardized <- TRUE

    # always print the 'short' summary
    short.summary(object)

    # only if requested, the fit measures
    if(fit.measures) {
        if(object@Options$test == "none") {
            warning("lavaan WARNING: fit measures not available if test = \"none\"\n\n")
        } else if(object@optim$npar > 0L && !object@optim$converged) {
            warning("lavaan WARNING: fit measures not available if model did not converge\n\n")
        } else {
            print.fit.measures( fitMeasures(object, fit.measures="default") )
        }
    }


    if(estimates) {

    # main part: parameter estimates
    cat("\nParameter estimates:\n\n")
    t0.txt <- sprintf("  %-40s", "Information")
    tmp.txt <- object@Options$information
    t1.txt <- sprintf("  %10s", paste(toupper(substring(tmp.txt,1,1)), 
			 	     substring(tmp.txt,2), sep=""))
    cat(t0.txt, t1.txt, "\n", sep="")
    t0.txt <- sprintf("  %-31s", "Standard Errors")
    tmp.txt <- object@Options$se
    t1.txt <- sprintf("  %19s", paste(toupper(substring(tmp.txt,1,1)),  
                                      substring(tmp.txt,2), sep=""))
    cat(t0.txt, t1.txt, "\n", sep="")
    if(object@Options$se == "bootstrap") {
        t0.txt <- sprintf("  %-40s", "Number of requested bootstrap draws")
        t1.txt <- sprintf("  %10i", object@Options$bootstrap)
        cat(t0.txt, t1.txt, "\n", sep="")
        t0.txt <- sprintf("  %-40s", "Number of successful bootstrap draws")
        t1.txt <- sprintf("  %10i", NROW(object@boot$coef))
        cat(t0.txt, t1.txt, "\n", sep="")
    }
    cat("\n")

    # local print function
    print.estimate <- function(name="ERROR", i=1, z.stat=TRUE) {
       
        # cut name if (still) too long
        name <- strtrim(name, width=13L)

        if(!standardized) {
            if(is.na(se[i])) {
                txt <- sprintf("    %-13s %9.3f %8.3f\n", name, est[i], se[i])
            } else if(se[i] == 0) {
                txt <- sprintf("    %-13s %9.3f\n", name, est[i])
            } else if(est[i]/se[i] > 9999.999) {
                txt <- sprintf("    %-13s %9.3f %8.3f\n", name, est[i], se[i])
            } else if(!z.stat) {
                txt <- sprintf("    %-13s %9.3f %8.3f\n", name, est[i], se[i])
            } else {
                z <- est[i]/se[i]
                pval <- 2 * (1 - pnorm( abs(z) ))
                txt <- sprintf("    %-13s %9.3f %8.3f %8.3f %8.3f\n",
                               name, est[i], se[i], z, pval)
            }
        } else {
            if(is.na(se[i])) {
                txt <- sprintf("    %-13s %9.3f %8.3f                   %8.3f %8.3f\n", name, est[i], se[i], est.std[i], est.std.all[i])
            } else if(se[i] == 0) {
                txt <- sprintf("    %-13s %9.3f                            %8.3f %8.3f\n", name, est[i], est.std[i], est.std.all[i])
            } else if(est[i]/se[i] > 9999.999) {
                txt <- sprintf("    %-13s %9.3f %8.3f                   %8.3f %8.3f\n", name, est[i], se[i], est.std[i], est.std.all[i])
            } else if(!z.stat) {
                txt <- sprintf("    %-13s %9.3f %8.3f                   %8.3f %8.3f\n", name, est[i], se[i], est.std[i], est.std.all[i])
            } else {
                z <- est[i]/se[i]
                pval <- 2 * (1 - pnorm( abs(z) ))
                txt <- sprintf("    %-13s %9.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n",
                               name, est[i], se[i], z, pval, est.std[i], est.std.all[i])
            }
        }
        cat(txt)
    }

    est <- lav_object_inspect_est(object)
    se  <- lav_object_inspect_se(object)
    if(rsquare || standardized) {
        est.std    <- standardize.est.lv(object)
        if(std.nox) {
            est.std.all <- standardize.est.all.nox(object, est.std=est.std)
        } else {
            est.std.all <- standardize.est.all(object, est.std=est.std)
        }
    } 

    for(g in 1:object@Data@ngroups) {
        ov.names <- vnames(object@ParTable, "ov", group=g)
        lv.names <- vnames(object@ParTable, "lv", group=g)

        # group header
        if(object@Data@ngroups > 1) {
            if(g > 1) cat("\n\n")
            cat("Group ", g, 
                " [", object@Data@group.label[[g]], "]:\n\n", sep="")
        }

        # estimates header
        if(!standardized) {
            cat("                   Estimate  Std.err  Z-value  P(>|z|)\n")
        } else {
            if(std.nox) {
                cat("                   Estimate  Std.err  Z-value  P(>|z|)   Std.lv  Std.nox\n")
            }
            else {
                cat("                   Estimate  Std.err  Z-value  P(>|z|)   Std.lv  Std.all\n")
            }
        }

        makeNames <- function(NAMES, LABELS) {
            multiB <- FALSE
            if(any(nchar(NAMES) != nchar(NAMES, "bytes")))
                multiB <- TRUE
            if(any(nchar(LABELS) != nchar(LABELS, "bytes")))
                multiB <- TRUE
            # labels?
            l.idx <- which(nchar(LABELS) > 0L)
            if(length(l.idx) > 0L) {
                if(!multiB) {
                    LABELS <- abbreviate(LABELS, 4)
                    LABELS[l.idx] <- paste(" (", LABELS[l.idx], ")", sep="")
                    MAX.L <- max(nchar(LABELS))
                    NAMES <- abbreviate(NAMES, minlength = (13 - MAX.L), 
                                        strict = TRUE)
                } else {
                    # do not abbreviate anything (eg in multi-byte locales)
                    MAX.L <- 4L
                }
                NAMES <- sprintf(paste("%-", (13 - MAX.L), "s%", MAX.L, "s",
                                       sep=""), NAMES, LABELS)
            } else {
                if(!multiB) {
                    NAMES <- abbreviate(NAMES, minlength = 13, strict = TRUE)
                } else {
                    NAMES <- sprintf(paste("%-", 13, "s", sep=""), NAMES)
                }
            }

            NAMES
        }

        NAMES <- object@ParTable$rhs

        # 1a. indicators ("=~") (we do show dummy indicators)
        mm.idx <- which( object@ParTable$op == "=~" & 
                        !object@ParTable$lhs %in% ov.names &
                         object@ParTable$group == g)
        if(length(mm.idx)) {
            cat("Latent variables:\n")
            lhs.old <- ""
            NAMES[mm.idx] <- makeNames(  object@ParTable$rhs[mm.idx],
                                       object@ParTable$label[mm.idx])
            for(i in mm.idx) {
                lhs <- object@ParTable$lhs[i]
                if(lhs != lhs.old) cat("  ", lhs, " =~\n", sep="")
                print.estimate(name=NAMES[i], i)
                lhs.old <- lhs
            }
            cat("\n")
        }

        # 1b. formative/composites ("<~")
        fm.idx <- which( object@ParTable$op == "<~" &
                         object@ParTable$group == g)
        if(length(fm.idx)) {
            cat("Composites:\n")
            lhs.old <- ""
            NAMES[fm.idx] <- makeNames(  object@ParTable$rhs[fm.idx],
                                       object@ParTable$label[fm.idx])
            for(i in fm.idx) {
                lhs <- object@ParTable$lhs[i]
                if(lhs != lhs.old) cat("  ", lhs, " <~\n", sep="")
                print.estimate(name=NAMES[i], i)
                lhs.old <- lhs
            }
            cat("\n")
        }

        # 2. regressions
        eqs.idx <- which(object@ParTable$op == "~" & object@ParTable$group == g)
        if(length(eqs.idx) > 0) {
            cat("Regressions:\n")
            lhs.old <- ""
            NAMES[eqs.idx] <- makeNames(  object@ParTable$rhs[eqs.idx],
                                        object@ParTable$label[eqs.idx])
            for(i in eqs.idx) {
                lhs <- object@ParTable$lhs[i]
                if(lhs != lhs.old) cat("  ", lhs, " ~\n", sep="")
                print.estimate(name=NAMES[i], i)
                lhs.old <- lhs
            }
            cat("\n")
        }

        # 3. covariances
        cov.idx <- which(object@ParTable$op == "~~" & 
                         !object@ParTable$exo &
                         object@ParTable$lhs != object@ParTable$rhs &
                         object@ParTable$group == g)
        if(length(cov.idx) > 0) {
            cat("Covariances:\n")
            lhs.old <- ""
            NAMES[cov.idx] <- makeNames(  object@ParTable$rhs[cov.idx],
                                        object@ParTable$label[cov.idx])
            for(i in cov.idx) {
                lhs <- object@ParTable$lhs[i]
                if(lhs != lhs.old) cat("  ", lhs, " ~~\n", sep="")
                print.estimate(name=NAMES[i], i)
                lhs.old <- lhs
            }
            cat("\n")
        }

        # 4. intercepts/means
        #ord.names <- vnames(object@ParTable, type="ov.ord", group=g)
        int.idx <- which(object@ParTable$op == "~1" & 
                         #!object@ParTable$lhs %in% ord.names &
                         !object@ParTable$exo &
                         object@ParTable$group == g)
        if(length(int.idx) > 0) {
            cat("Intercepts:\n")
            NAMES[int.idx] <- makeNames(  object@ParTable$lhs[int.idx],
                                        object@ParTable$label[int.idx])
            for(i in int.idx) {
                print.estimate(name=NAMES[i], i)
            }
            cat("\n")
        }

        # 4b thresholds
        th.idx <- which(object@ParTable$op == "|" &
                        object@ParTable$group == g)
        if(length(th.idx) > 0) {
            cat("Thresholds:\n")
            NAMES[th.idx] <- makeNames(  paste(object@ParTable$lhs[th.idx],
                                               "|",
                                               object@ParTable$rhs[th.idx],
                                               sep=""),
                                         object@ParTable$label[th.idx])
            for(i in th.idx) {
                print.estimate(name=NAMES[i], i)
            }
            cat("\n")
        }

        # 5. (residual) variances
        var.idx <- which(object@ParTable$op == "~~" &
                         !object@ParTable$exo &
                         object@ParTable$lhs == object@ParTable$rhs &
                         object@ParTable$group == g)
        if(length(var.idx) > 0) {
            cat("Variances:\n")
            NAMES[var.idx] <- makeNames(  object@ParTable$rhs[var.idx],
                                        object@ParTable$label[var.idx])
            for(i in var.idx) {
                if(object@Options$mimic == "lavaan") {
                    print.estimate(name=NAMES[i], i, z.stat=FALSE)
                } else {
                    print.estimate(name=NAMES[i], i, z.stat=TRUE)
                }
            }
            cat("\n")
        }

        # 6. latent response scales
        delta.idx <- which(object@ParTable$op == "~*~" &
                         object@ParTable$group == g)
        if(length(delta.idx) > 0) {
            cat("Scales y*:\n")
            NAMES[delta.idx] <- makeNames(  object@ParTable$rhs[delta.idx],
                                            object@ParTable$label[delta.idx])
            for(i in delta.idx) {
                print.estimate(name=NAMES[i], i, z.stat=TRUE)
            }
            cat("\n")
        }

        # 7. group weight
        group.idx <- which(object@ParTable$lhs == "group" &
                           object@ParTable$op == "%" &
                           object@ParTable$group == g)
        if(length(group.idx) > 0) {
            cat("Group weight:\n")
            NAMES[group.idx] <- makeNames(  object@ParTable$rhs[group.idx],
                                            object@ParTable$label[group.idx])
            for(i in group.idx) {
                print.estimate(name=NAMES[i], i, z.stat=TRUE)
            }
            cat("\n")
        }

    } # ngroups

    # 6. variable definitions
    def.idx <- which(object@ParTable$op == ":=")
    if(length(def.idx) > 0) {
        if(object@Data@ngroups > 1) cat("\n")
        cat("Defined parameters:\n")
        NAMES[def.idx] <- makeNames(  object@ParTable$lhs[def.idx], "")
        for(i in def.idx) {
            print.estimate(name=NAMES[i], i)
        }
        cat("\n")
    }

    # 7. constraints
    cin.idx <- which((object@ParTable$op == "<" | 
                      object@ParTable$op == ">"))
    ceq.idx <- which(object@ParTable$op == "==" & object@ParTable$user == 1L)
    if(length(cin.idx) > 0L || length(ceq.idx) > 0L) {
        # set small negative values to zero, to avoid printing " -0.000"
        slack <- ifelse(abs(est) < 1e-5, 0, est)
        #slack[cin.idx] <- object@Model@cin.function(object@optim$x)
        #slack[ceq.idx] <- object@Model@ceq.function(object@optim$x)
       
        if(object@Data@ngroups > 1 && length(def.idx) == 0L) cat("\n")
        cat("Constraints:                               Slack (>=0)\n")
        for(i in c(cin.idx,ceq.idx)) {
            lhs <- object@ParTable$lhs[i]
             op <- object@ParTable$op[i]
            rhs <- object@ParTable$rhs[i]
            if(rhs == "0" && op == ">") {
                con.string <- paste(lhs, " - 0", sep="")
            } else if(rhs == "0" && op == "<") {
                con.string <- paste(rhs, " - (", lhs, ")", sep="")
            } else if(rhs != "0" && op == ">") {
                con.string <- paste(lhs, " - (", rhs, ")", sep="")
            } else if(rhs != "0" && op == "<") {
                con.string <- paste(rhs, " - (", lhs, ")", sep="")
            } else if(rhs == "0" && op == "==") {
                con.string <- paste(lhs, " - 0", sep="")
            } else if(rhs != "0" && op == "==") {
                con.string <- paste(lhs, " - (", rhs, ")", sep="")
            }
            con.string <- abbreviate(con.string, 41, strict = TRUE)
            txt <- sprintf("    %-41s %8.3f\n", 
                           con.string, slack[i])
            cat(txt)
        }   
        cat("\n")
    }

    } # parameter estimates


    # R-square?
    if(rsquare) {
        r2 <- lav_object_inspect_rsquare(object, est.std.all=est.std.all,
                  drop.list.single.group = FALSE, add.labels = TRUE,
                  add.class = FALSE)
        for(g in 1:object@Data@ngroups) {
            if(object@Data@ngroups > 1) {
                cat("R-Square Group ", g, " [", 
                    object@Data@group.label[[g]], "]",
                    ":\n\n", sep="")
            } else {
                cat("R-Square:\n\n")
            } 
            for(i in 1:length(r2[[g]])) {
                t1.txt <- sprintf("    %-13s %9.3f\n", names(r2[[g]])[i], 
                                  r2[[g]][i])
                cat(t1.txt)
            }
            if(g < object@Data@ngroups) cat("\n")
        }
    }

    # modification indices?
    if(modindices) {
        cat("Modification Indices:\n\n")
        print( modificationIndices(object, standardized=TRUE) )
    }

}



setMethod("coef", "lavaan",
function(object, type="free", labels=TRUE) {

    if(type == "user" || type == "all") {
        type <- "user"
        idx <- 1:length( object@ParTable$lhs )
    } else if(type == "free") {
        idx <- which(object@ParTable$free > 0L & !duplicated(object@ParTable$free))
    } else {
        stop("lavaan ERROR: argument `type' must be one of free or user")
    }
    EST <- lav_object_inspect_est(object)
    cof <- EST[idx]
  
    # labels?
    if(labels) names(cof) <- lav_partable_labels(object@ParTable, type=type)

    # class
    class(cof) <- c("lavaan.vector", "numeric")

    cof
})

standardizedSolution <- standardizedsolution <- function(object,
                                                         type = "std.all",
                                                         se = TRUE,
                                                         zstat = TRUE,
                                                         pvalue = TRUE,
                                                         remove.eq = TRUE,
                                                         remove.ineq = TRUE,
                                                         remove.def = FALSE) {

    stopifnot(type %in% c("std.all", "std.lv", "std.nox"))

    # no zstat + pvalue if estimator is Bayes
    if(object@Options$estimator == "Bayes") {
        zstat <- pvalue <- FALSE
    }

    # no se if class is not lavaan
    if(class(object) != "lavaan") {
        if(missing(se) || !se) {
            se <- FALSE
            zstat <- FALSE
            pvalue <- FALSE
        }
    }

    PARTABLE <- inspect(object, "list")
    free.idx <- which(PARTABLE$free > 0L)
    LIST <- PARTABLE[,c("lhs", "op", "rhs")]
    if(!is.null(PARTABLE$group)) {
        LIST$group <- PARTABLE$group
    }

    # add std and std.all columns
    if(type == "std.lv") {
        LIST$est.std     <- standardize.est.lv(object)
    } else if(type == "std.all") {
        LIST$est.std <- standardize.est.all(object)
    } else if(type == "std.nox") {
        LIST$est.std <- standardize.est.all.nox(object)
    }

    if(object@Options$se != "none" && se) {
        # add 'se' for standardized parameters
        VCOV <- try(lav_object_inspect_vcov(object, standardized = TRUE,
                                            type = type, free.only = FALSE,
                                            add.labels = FALSE, 
                                            add.class = FALSE))
        if(inherits(VCOV, "try-error")) {
            LIST$se <- rep(NA, length(LIST$lhs))
            if(zstat) {
                LIST$z  <- rep(NA, length(LIST$lhs))
            }
            if(pvalue) {
                 LIST$pvalue <- rep(NA, length(LIST$lhs))
            }
        } else {
            tmp <- diag(VCOV)
            # catch negative values
            min.idx <- which(tmp < 0)
            if(length(min.idx) > 0L) {
                tmp[min.idx] <- as.numeric(NA)
            }
            # now, we can safely take the square root
            tmp <- sqrt(tmp)
            # catch near-zero SEs
            zero.idx <- which(tmp < sqrt(.Machine$double.eps))
            if(length(zero.idx) > 0L) {
                tmp[zero.idx] <- 0.0
            }
            LIST$se <- tmp
 
            # add 'z' column
            if(zstat) {
                 tmp.se <- ifelse( LIST$se == 0.0, NA, LIST$se)
                 LIST$z <- LIST$est.std / tmp.se
            }
            if(zstat && pvalue) {
                 LIST$pvalue <- 2 * (1 - pnorm( abs(LIST$z) ))
            }
        }
    }

    # if single group, remove group column
    if(object@Data@ngroups == 1L) LIST$group <- NULL

    # remove == rows?
    if(remove.eq) {
        eq.idx <- which(LIST$op == "==")
        if(length(eq.idx) > 0L) {
            LIST <- LIST[-eq.idx,]
        }
    }
    # remove <> rows?
    if(remove.ineq) {
        ineq.idx <- which(LIST$op == "<" || LIST$op == ">")
        if(length(ineq.idx) > 0L) {
            LIST <- LIST[-ineq.idx,]
        }
    }
    # remove := rows?
    if(remove.def) {
        def.idx <- which(LIST$op == ":=")
        if(length(def.idx) > 0L) {
            LIST <- LIST[-def.idx,]
        }
    }

    # always add attributes (for now)
    class(LIST) <- c("lavaan.data.frame", "data.frame")
    LIST
}

parameterEstimates <- parameterestimates <- function(object,
                                                     se    = TRUE,
                                                     zstat = TRUE,
                                                     pvalue = TRUE,
                                                     ci = TRUE, 
                                                     level = 0.95, 
                                                     boot.ci.type = "perc",
                                                     standardized = FALSE, 
                                                     fmi = FALSE,
                                                     remove.system.eq = TRUE,
                                                     remove.eq = TRUE,
                                                     remove.ineq = TRUE,
                                                     remove.def = FALSE,
                                                     rsquare = FALSE,
                                                     add.attributes = FALSE) {

    # no se if class is not lavaan
    if(class(object) != "lavaan") {
        if(missing(se) || !se) {
            se <- FALSE
            zstat <- FALSE
            pvalue <- FALSE
        }
    }

    # check fmi
    if(fmi) {
        if(inherits(object, "lavaanList")) {
            warning("lavaan WARNING: fmi not available for object of class \"lavaanList\"")
            fmi <- FALSE
        }
        if(object@Options$se != "standard") {
            warning("lavaan WARNING: fmi only available if se = \"standard\"")
            fmi <- FALSE
        }
        if(object@Options$estimator != "ML") {
            warning("lavaan WARNING: fmi only available if estimator = \"ML\"")
            fmi <- FALSE
        }
        if(!object@SampleStats@missing.flag) {
            warning("lavaan WARNING: fmi only available if missing = \"(fi)ml\"")
            fmi <- FALSE
        }
        if(!object@optim$converged) {
            warning("lavaan WARNING: fmi not available; model did not converge")
            fmi <- FALSE
        }
    }

    # no zstat + pvalue if estimator is Bayes
    if(object@Options$estimator == "Bayes") {
        zstat <- pvalue <- FALSE
    }

    PARTABLE <- as.data.frame(object@ParTable, stringsAsFactors = FALSE)
    LIST <- PARTABLE[,c("lhs", "op", "rhs")]
    if(!is.null(PARTABLE$user)) {
        LIST$user <- PARTABLE$user
    }
    if(!is.null(PARTABLE$group)) {
        LIST$group <- PARTABLE$group
    } else {
        LIST$group <- rep(1L, length(LIST$lhs))
    }
    if(!is.null(PARTABLE$label)) {
        LIST$label <- PARTABLE$label
    } else {
        LIST$label <- rep("", length(LIST$lhs))
    }
    if(!is.null(PARTABLE$exo)) {
        LIST$exo <- PARTABLE$exo
    } else {
        LIST$exo <- rep(0L, length(LIST$lhs))
    }
    if(inherits(object, "lavaanList")) {
        # per default: nothing!
        #if("partable" %in% object@meta$store.slots) {
        #    COF <- sapply(object@ParTableList, "[[", "est")
        #    LIST$est <- rowMeans(COF)
        #} 
        LIST$est <- NULL
    } else if(!is.null(PARTABLE$est)) {
        LIST$est <- PARTABLE$est
    } else {
        LIST$est <- lav_model_get_parameters(object@Model, type = "user",
                                             extra = TRUE)    
    }


    # add se, zstat, pvalue
    if(se && object@Options$se != "none") {
        LIST$se <- lav_object_inspect_se(object)
        tmp.se <- ifelse(LIST$se == 0.0, NA, LIST$se)
        if(zstat) { 
            LIST$z <- LIST$est / tmp.se
            if(pvalue) {
                LIST$pvalue <- 2 * (1 - pnorm( abs(LIST$z) ))
            }
        }
    }

    # extract bootstrap data (if any)
    BOOT <- lav_object_inspect_boot(object)
    bootstrap.successful <- NROW(BOOT) # should be zero if NULL

    # confidence interval
    if(se && object@Options$se != "none" && ci) {
        # next three lines based on confint.lm
        a <- (1 - level)/2; a <- c(a, 1 - a)
        if(object@Options$se != "bootstrap") {
            fac <- qnorm(a)
            ci <- LIST$est + LIST$se %o% fac
        } else if(object@Options$se == "bootstrap") {

            # local copy of 'norm.inter' from boot package (not exported!)
            norm.inter <- function(t, alpha)  {
                t <- t[is.finite(t)]; R <- length(t); rk <- (R + 1) * alpha
                if (!all(rk > 1 & rk < R)) 
                     warning("extreme order statistics used as endpoints")
                k <- trunc(rk); inds <- seq_along(k)
                out <- inds; kvs <- k[k > 0 & k < R]
                tstar <- sort(t, partial = sort(union(c(1, R), c(kvs, kvs+1))))
                ints <- (k == rk)
                if (any(ints)) out[inds[ints]] <- tstar[k[inds[ints]]]
                out[k == 0] <- tstar[1L]
                out[k == R] <- tstar[R]
                not <- function(v) xor(rep(TRUE,length(v)),v)
                temp <- inds[not(ints) & k != 0 & k != R]
                temp1 <- qnorm(alpha[temp])
                temp2 <- qnorm(k[temp]/(R+1))
                temp3 <- qnorm((k[temp]+1)/(R+1))
                tk <- tstar[k[temp]]
                tk1 <- tstar[k[temp]+1L]
                out[temp] <- tk + (temp1-temp2)/(temp3-temp2)*(tk1 - tk)
                cbind(round(rk, 2), out)
            }

            stopifnot(!is.null(BOOT))
            stopifnot(boot.ci.type %in% c("norm","basic","perc","bca.simple"))
            if(boot.ci.type == "norm") {
                fac <- qnorm(a)
                boot.x <- colMeans(BOOT)
                boot.est <- 
                    lav_model_get_parameters(object@Model, 
                                       GLIST=lav_model_x2GLIST(object@Model, boot.x), 
                                       type="user", extra=TRUE)
                bias.est <- (boot.est - LIST$est)
                ci <- (LIST$est - bias.est) + LIST$se %o% fac
            } else if(boot.ci.type == "basic") {
                ci <- cbind(LIST$est, LIST$est)
                alpha <- (1 + c(level, -level))/2

                # free.idx only
                qq <- apply(BOOT, 2, norm.inter, alpha)
                free.idx <- which(object@ParTable$free & 
                                  !duplicated(object@ParTable$free))
                ci[free.idx,] <- 2*ci[free.idx,] - t(qq[c(3,4),])

                # def.idx
                def.idx <- which(object@ParTable$op == ":=")
                if(length(def.idx) > 0L) {
                    BOOT.def <- apply(BOOT, 1, object@Model@def.function)
                    if(length(def.idx) == 1L) {
                        BOOT.def <- as.matrix(BOOT.def)
                    } else {
                        BOOT.def <- t(BOOT.def)
                    }
                    qq <- apply(BOOT.def, 2, norm.inter, alpha)
                    ci[def.idx,] <- 2*ci[def.idx,] - t(qq[c(3,4),])
                }

                # TODO: add cin/ceq?
               
            } else if(boot.ci.type == "perc") {
                ci <- cbind(LIST$est, LIST$est)
                alpha <- (1 + c(-level, level))/2

                # free.idx only
                qq <- apply(BOOT, 2, norm.inter, alpha)
                free.idx <- which(object@ParTable$free & 
                                  !duplicated(object@ParTable$free))
                ci[free.idx,] <- t(qq[c(3,4),])

                # def.idx
                def.idx <- which(object@ParTable$op == ":=")
                if(length(def.idx) > 0L) {
                    BOOT.def <- apply(BOOT, 1, object@Model@def.function)
                    if(length(def.idx) == 1L) {
                        BOOT.def <- as.matrix(BOOT.def)
                    } else {
                        BOOT.def <- t(BOOT.def)
                    }
                    qq <- apply(BOOT.def, 2, norm.inter, alpha)
                    def.idx <- which(object@ParTable$op == ":=")
                    ci[def.idx,] <- t(qq[c(3,4),])
                }

                # TODO:  add cin/ceq?

            } else if(boot.ci.type == "bca.simple") {
               # no adjustment for scale!! only bias!!
               alpha <- (1 + c(-level, level))/2
               zalpha <- qnorm(alpha)
               ci <- cbind(LIST$est, LIST$est)

               # free.idx only
               free.idx <- which(object@ParTable$free & 
                                 !duplicated(object@ParTable$free))
               x <- LIST$est[free.idx]
               for(i in 1:length(free.idx)) {
                   t <- BOOT[,i]; t <- t[is.finite(t)]; t0 <- x[i]
                   w <- qnorm(sum(t < t0)/length(t))
                   a <- 0.0 #### !!! ####
                   adj.alpha <- pnorm(w + (w + zalpha)/(1 - a*(w + zalpha)))
                   qq <- norm.inter(t, adj.alpha)
                   ci[free.idx[i],] <- qq[,2]
               }

               # def.idx
               def.idx <- which(object@ParTable$op == ":=")
               if(length(def.idx) > 0L) {
                   x.def <- object@Model@def.function(x)
                   BOOT.def <- apply(BOOT, 1, object@Model@def.function)
                   if(length(def.idx) == 1L) {
                       BOOT.def <- as.matrix(BOOT.def)
                   } else {
                       BOOT.def <- t(BOOT.def)
                   }
                   for(i in 1:length(def.idx)) {
                       t <- BOOT.def[,i]; t <- t[is.finite(t)]; t0 <- x.def[i]
                       w <- qnorm(sum(t < t0)/length(t))
                       a <- 0.0 #### !!! ####
                       adj.alpha <- pnorm(w + (w + zalpha)/(1 - a*(w + zalpha)))
                       qq <- norm.inter(t, adj.alpha)
                       ci[def.idx[i],] <- qq[,2]
                   }
               }
               
               # TODO:
               # - add cin/ceq
            }
        }

        LIST$ci.lower <- ci[,1]; LIST$ci.upper <- ci[,2]    
    }

    # standardized estimates?
    if(standardized) {
        LIST$std.lv  <- standardize.est.lv(object)
        LIST$std.all <- standardize.est.all(object, est.std=LIST$est.std)
        LIST$std.nox <- standardize.est.all.nox(object, est.std=LIST$est.std)
    }

    # rsquare?
    if(rsquare) {
        r2 <- lavTech(object, "rsquare", add.labels = TRUE)
        NAMES <- unlist(lapply(r2, names)); nel <- length(NAMES)
        R2 <- data.frame( lhs = NAMES, op = rep("r2", nel), rhs = NAMES,
                          group = rep(1:length(r2), sapply(r2, length)),
                          est = unlist(r2), stringsAsFactors = FALSE )
        LIST <- lav_partable_merge(pt1 = LIST, pt2 = R2, warn = FALSE)
    }

    # fractional missing information (if estimator="fiml")
    if(fmi) {
        SE.orig <- LIST$se
        lavmodel <- object@Model; implied <- object@implied
        COV <- if(lavmodel@conditional.x) implied$res.cov else implied$cov
        MEAN <- if(lavmodel@conditional.x) implied$res.int else implied$mean

        # provide rownames
        for(g in 1:object@Data@ngroups)
            rownames(COV[[g]]) <- object@Data@ov.names[[g]]

        # if estimator="ML" and likelihood="normal" --> rescale
        if(object@Options$estimator == "ML" &&
           object@Options$likelihood == "normal") {
            for(g in 1:object@Data@ngroups) {
                N <- object@Data@nobs[[g]]
                COV[[g]] <- (N+1)/N * COV[[g]]
            }
        }

        # fit another model, using the model-implied moments as input data
        step2 <- lavaan(slotOptions  = object@Options,
                        slotParTable = object@ParTable,
                        sample.cov   = COV,
                        sample.mean  = MEAN,
                        sample.nobs  = object@Data@nobs)
        SE2 <- lav_object_inspect_se(step2)
        SE.step2 <- ifelse(SE2 == 0.0, as.numeric(NA), SE2)
        if(rsquare) {
            # add additional elements, since LIST$se is now longer
            r2.idx <- which(LIST$op == "r2")
            if(length(r2.idx) > 0L) {
                SE.step2 <- c(SE.step2, rep(as.numeric(NA), length(r2.idx)))
            }
        }
        LIST$fmi <- 1-(SE.step2*SE.step2/(SE.orig*SE.orig))
    }

    # if single group, remove group column
    if(object@Data@ngroups == 1L) LIST$group <- NULL

    # if no user-defined labels, remove label column
    if(sum(nchar(object@ParTable$label)) == 0L) LIST$label <- NULL

    # remove == rows?
    if(remove.eq) {
        eq.idx <- which(LIST$op == "==" & LIST$user == 1L)
        if(length(eq.idx) > 0L) {
            LIST <- LIST[-eq.idx,]
        }
    }
    if(remove.system.eq) {
        eq.idx <- which(LIST$op == "==" & LIST$user != 1L)
        if(length(eq.idx) > 0L) {
            LIST <- LIST[-eq.idx,]
        }
    }
    # remove <> rows?
    if(remove.ineq) {
        ineq.idx <- which(LIST$op == "<" || LIST$op == ">")
        if(length(ineq.idx) > 0L) {
            LIST <- LIST[-ineq.idx,]
        }
    }
    # remove := rows?
    if(remove.def) {
        def.idx <- which(LIST$op == ":=")
        if(length(def.idx) > 0L) {
            LIST <- LIST[-def.idx,]
        }
    }

    # remove LIST$user
    LIST$user <- NULL

    if(add.attributes) {
        class(LIST) <- c("lavaan.parameterEstimates", "lavaan.data.frame",
                         "data.frame")
        attr(LIST, "information") <- object@Options$information
        attr(LIST, "se") <- object@Options$se
        attr(LIST, "group.label") <- object@Data@group.label
        attr(LIST, "bootstrap") <- object@Options$bootstrap
        attr(LIST, "bootstrap.successful") <- bootstrap.successful
        # FIXME: add more!!
    } else {
        LIST$exo <- NULL
        class(LIST) <- c("lavaan.data.frame", "data.frame")
    }

    LIST
}

parameterTable <- parametertable <- parTable <- partable <-
        function(object) {

    # convert to data.frame
    out <- as.data.frame(object@ParTable, stringsAsFactors = FALSE)

    class(out) <- c("lavaan.data.frame", "data.frame")
    out
}

varTable <- vartable <- function(object, ov.names=names(object), 
                                 ov.names.x=NULL, 
                                 ordered = NULL, factor = NULL,
                                 as.data.frame.=TRUE) {

    if(inherits(object, "lavaan")) {
        VAR <- object@Data@ov
    } else if(inherits(object, "lavData")) {
        VAR <- object@ov
    } else if(inherits(object, "data.frame")) {
        VAR <- lav_dataframe_vartable(frame = object, ov.names = ov.names, 
                                      ov.names.x = ov.names.x, 
                                      ordered = ordered, factor = factor,
                                      as.data.frame. = FALSE)
    } else {
        stop("object must of class lavaan or a data.frame")
    } 

    if(as.data.frame.) {
        VAR <- as.data.frame(VAR, stringsAsFactors=FALSE,
                             row.names=1:length(VAR$name))
        class(VAR) <- c("lavaan.data.frame", "data.frame")
    }

    VAR
}


setMethod("fitted.values", "lavaan",
function(object, type = "moments", labels=TRUE) {

    # lowercase type
    type <- tolower(type)

    # catch type="casewise"
    if(type %in% c("casewise","case","obs","observations","ov")) {
        return( lavPredict(object, type = "ov", label = labels) )
    }

    lav_object_inspect_implied(object,
               add.labels = labels, add.class = TRUE,
               drop.list.single.group = TRUE)
})


setMethod("fitted", "lavaan",
function(object, type = "moments", labels=TRUE) {
     fitted.values(object, type = type, labels = labels)
})


setMethod("vcov", "lavaan",
function(object, labels = TRUE, remove.duplicated = FALSE) {

    # check for convergence first!
    if(object@optim$npar > 0L && !object@optim$converged)
        stop("lavaan ERROR: model did not converge")

    if(object@Options$se == "none") {
        stop("lavaan ERROR: vcov not available if se=\"none\"")
    }

    VarCov <- lav_object_inspect_vcov(object,
                                      add.labels = labels,
                                      add.class = TRUE,
                                      remove.duplicated = remove.duplicated)

    VarCov
})


# logLik (so that we can use the default AIC/BIC functions from stats4(
setMethod("logLik", "lavaan",
function(object, ...) {
    if(object@Options$estimator != "ML") {
        warning("lavaan WARNING: logLik only available if estimator is ML")
    }
    if(object@optim$npar > 0L && !object@optim$converged) {
        warning("lavaan WARNING: model did not converge")
    }
    
    logl.df <- fitMeasures(object, c("logl", "npar", "ntotal"))
    names(logl.df) <- NULL
    logl <- logl.df[1]
    attr(logl, "df") <- logl.df[2]    ### note: must be npar, not df!!
    attr(logl, "nobs") <- logl.df[3]
    class(logl) <- "logLik"
    logl
})

# nobs
if(!exists("nobs", envir=asNamespace("stats4"))) {
    setGeneric("nobs", function(object, ...) standardGeneric("nobs"))
}
setMethod("nobs", signature(object = "lavaan"),
function(object, ...) {
    object@SampleStats@ntotal
})

# see: src/library/stats/R/update.R
setMethod("update", signature(object = "lavaan"),
function(object, model, ..., evaluate = TRUE) {

    call <- object@call
    if(is.null(call))
        stop("need an object with call slot")

    extras <- match.call(expand.dots = FALSE)$...

    if(!missing(model))
        #call$formula <- update.formula(formula(object), formula.)
        call$model <- model

    if(length(extras) > 0) {
        existing <- !is.na(match(names(extras), names(call)))
        for(a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if(any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    if (evaluate) {
        eval(call, parent.frame())
    }
    else call
})


setMethod("anova", signature(object = "lavaan"),
function(object, ...) {

    # NOTE: if we add additional arguments, it is not the same generic
    # anova() function anymore, and match.call will be screwed up

    # NOTE: we need to extract the names of the models from match.call here,
    #       otherwise, we loose them in the call stack

    mcall <- match.call(expand.dots = TRUE)
    dots <- list(...)

    # catch SB.classic and SB.H0
    SB.classic <- TRUE; SB.H0 <- FALSE

    arg.names <- names(dots)
    arg.idx <- which(nchar(arg.names) > 0L)
    if(length(arg.idx) > 0L) {
        if(!is.null(dots$SB.classic))
            SB.classic <- dots$SB.classic
        if(!is.null(dots$SB.H0))
            SB.H0 <- dots$SB.H0           
        dots <- dots[-arg.idx]
    }

    modp <- if(length(dots))
        sapply(dots, is, "lavaan") else logical(0)
    mods <- c(list(object), dots[modp])
    NAMES <- sapply(as.list(mcall)[c(FALSE, TRUE, modp)], deparse)

    # use do.call to handle changed dots
    ans <- do.call("lavTestLRT", c(list(object = object, 
                   SB.classic = SB.classic, SB.H0 = SB.H0, 
                   model.names = NAMES), dots))

    ans
})


