#
# initial version: YR 25/03/2009

short.summary <- function(object) {

    # Convergence or not?
    if(object@Fit@iterations > 0) {
        if(object@Fit@converged) {
	    cat(sprintf("lavaan (%s) converged normally after %3i iterations\n",
                    packageDescription("lavaan", fields="Version"),
                    object@Fit@iterations))
        } else {
            cat(sprintf("** WARNING ** lavaan (%s) did NOT converge after %i iterations\n", 
                packageDescription("lavaan", fields="Version"),
                object@Fit@iterations))
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
    #t1.txt <- sprintf("  %10i", object@Fit@npar)
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
       length(object@Fit@test) > 1L) {
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
    if(object@Options$test != "none") {

        # 1. chi-square values
        t0.txt <- sprintf("  %-40s", "Minimum Function Test Statistic")  
        t1.txt <- sprintf("  %10.3f", object@Fit@test[[1]]$stat)
        t2.txt <- ifelse(scaled, 
                  sprintf("  %10.3f", object@Fit@test[[2]]$stat), "")
        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

        # 2. degrees of freedom
        t0.txt <- sprintf("  %-40s", "Degrees of freedom")
        t1.txt <- sprintf("  %10i",   object@Fit@test[[1]]$df)
        t2.txt <- ifelse(scaled, 
                         ifelse(round(object@Fit@test[[2]]$df) == 
                                object@Fit@test[[2]]$df,
                                sprintf("  %10i",   object@Fit@test[[2]]$df),
                                sprintf("  %10.3f", object@Fit@test[[2]]$df)),
                         "")
        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

        # 3. P-value
        if(object@Fit@test[[1]]$refdistr == "chisq") {
            t0.txt <- sprintf("  %-40s", "P-value (Chi-square)")
        } else if(length(object@Fit@test) == 1L &&
                  object@Fit@test[[1]]$refdistr == "unknown") {
            t0.txt <- sprintf("  %-40s", "P-value (Unknown)")
        } else {
            t0.txt <- sprintf("  %-40s", "P-value")
        }
        t1.txt <- sprintf("  %10.3f", object@Fit@test[[1]]$pvalue)
        t2.txt <- ifelse(scaled,
                  sprintf("  %10.3f", object@Fit@test[[2]]$pvalue), "")
        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

        # 3b. Do we have a Bollen-Stine p-value?
        if(object@Options$test == "bollen.stine") {
            t0.txt <- sprintf("  %-40s", "P-value (Bollen-Stine Bootstrap)")
            t1.txt <- sprintf("  %10.3f", object@Fit@test[[2]]$pvalue)
            cat(t0.txt, t1.txt, "\n", sep="")
        }

        # 4. Scaling correction factor
        if(scaled) {
            t0.txt <- sprintf("  %-40s", "Scaling correction factor")
            t1.txt <- sprintf("  %10s", "")
            t2.txt <- sprintf("  %10.3f", object@Fit@test[[2]]$scaling.factor)
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
                                  object@Fit@test[[2]]$shift.parameter)
                cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
            } else { # multiple groups, multiple shift values!
                cat("  Shift parameter for each group:\n")
                for(g in 1:object@Data@ngroups) {
                    t0.txt <- sprintf("    %-38s", object@Data@group.label[[g]])
                    t1.txt <- sprintf("  %10s", "")
                    t2.txt <- sprintf("  %10.3f",
                                     object@Fit@test[[2]]$shift.parameter[g])
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
                t1.txt <- sprintf("  %10.3f", object@Fit@test[[1]]$stat.group[g])
                t2.txt <- ifelse(scaled, sprintf("  %10.3f", 
                                 object@Fit@test[[2]]$stat.group[g]), "")
                cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
            }
        } 
    } # test != none

    cat("\n")
}

setMethod("show", "lavaan",
function(object) {

    # show only basic information
    short.summary(object)

})

setMethod("summary", "lavaan",
function(object, estimates=TRUE, fit.measures=FALSE, standardized=FALSE, 
         rsquare=FALSE, std.nox=FALSE, modindices=FALSE) {

    if(std.nox) standardized <- TRUE

    # always print the 'short' summary
    short.summary(object)

    # only if requested, the fit measures
    if(fit.measures) {
        if(object@Options$test == "none") {
            warning("lavaan WARNING: fit measures not available if test = \"none\"\n\n")
        } else if(object@Fit@npar > 0L && !object@Fit@converged) {
            warning("lavaan WARNING: fit measures not available if model did not converge\n\n")
        } else {
            print.fit.measures( fitMeasures(object, fit.measures="default") )
        }
    }


    if(estimates) {

    # main part: parameter estimates
    cat("Parameter estimates:\n\n")
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
        t1.txt <- sprintf("  %10i", object@Options$boot)
        cat(t0.txt, t1.txt, "\n", sep="")
        t0.txt <- sprintf("  %-40s", "Number of successful bootstrap draws")
        t1.txt <- sprintf("  %10i", nrow(attr(object@Fit@est, "BOOT.COEF")))
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

    est <- object@Fit@est
    se  <- object@Fit@se
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
        ord.names <- vnames(object@ParTable, type="ov.ord", group=g)
        int.idx <- which(object@ParTable$op == "~1" & 
                         !object@ParTable$lhs %in% ord.names &
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
    ceq.idx <- which(object@ParTable$op == "==")
    if(length(cin.idx) > 0L || length(ceq.idx) > 0L) {
        # set small negative values to zero, to avoid printing " -0.000"
        slack <- ifelse(abs(est) < 1e-5, 0, est)
        #slack[cin.idx] <- object@Model@cin.function(object@Fit@x)
        #slack[ceq.idx] <- object@Model@ceq.function(object@Fit@x)
       
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
        r2 <- rsquare(object, est.std.all=est.std.all)
        if(object@Data@ngroups == 1L) r2 <- list(r2)
        for(g in 1:object@Data@ngroups) {
            if(object@Data@ngroups > 1) {
                cat("R-Square Group ", object@Data@group.label[[g]], 
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

})


free.parameters <- function(object, type="free") {

    GLIST <- object@Model@GLIST

    for(mm in 1:length(GLIST)) {
        # erase everything
        GLIST[[mm]][,] <- 0.0

        dimnames(GLIST[[mm]]) <- object@Model@dimNames[[mm]]

        # fill in free/unco parameter counts
        if(type == "free") {
            m.el.idx <- object@Model@m.free.idx[[mm]]
            x.el.idx <- object@Model@x.free.idx[[mm]]
        } else if(type == "unco") {
            m.el.idx <- object@Model@m.unco.idx[[mm]]
            x.el.idx <- object@Model@x.unco.idx[[mm]]
        } else if(type == "user") {
            m.el.idx <- object@Model@m.user.idx[[mm]]
            x.el.idx <- object@Model@x.user.idx[[mm]]
        } else {
            stop("lavaan ERROR: unknown type argument:", type, )
        }
        GLIST[[mm]][m.el.idx] <- x.el.idx

        # class
        if(object@Model@isSymmetric[mm]) {
            class(GLIST[[mm]]) <- c("lavaan.matrix.symmetric", "matrix")
        } else {
            class(GLIST[[mm]]) <- c("lavaan.matrix", "matrix")
        }
    }

    GLIST
}

standard.errors <- function(object) {

    GLIST <- object@Model@GLIST

    for(mm in 1:length(GLIST)) {
        # labels
        dimnames(GLIST[[mm]]) <- object@Model@dimNames[[mm]]

        # fill in starting values
        m.user.idx <- object@Model@m.user.idx[[mm]]
        x.user.idx <- object@Model@x.user.idx[[mm]]
        GLIST[[mm]][m.user.idx] <- object@Fit@se[x.user.idx]

        # class
        if(object@Model@isSymmetric[mm]) {
            class(GLIST[[mm]]) <- c("lavaan.matrix.symmetric", "matrix")
        } else {
            class(GLIST[[mm]]) <- c("lavaan.matrix", "matrix")
        }
    }

    GLIST
}

starting.values <- function(object) {

    GLIST <- object@Model@GLIST

    for(mm in 1:length(GLIST)) {
        # labels
        dimnames(GLIST[[mm]]) <- object@Model@dimNames[[mm]]

        # fill in starting values
        m.user.idx <- object@Model@m.user.idx[[mm]]
        x.user.idx <- object@Model@x.user.idx[[mm]]
        GLIST[[mm]][m.user.idx] <- object@Fit@start[x.user.idx]

        # class
        if(object@Model@isSymmetric[mm]) {
            class(GLIST[[mm]]) <- c("lavaan.matrix.symmetric", "matrix")
        } else {
            class(GLIST[[mm]]) <- c("lavaan.matrix", "matrix")
        }
    }

    GLIST
}

parameter.values <- function(object) {

    GLIST <- object@Model@GLIST

    for(mm in 1:length(GLIST)) {
        # labels
        dimnames(GLIST[[mm]]) <- object@Model@dimNames[[mm]]

        # fill in starting values
        m.user.idx <- object@Model@m.user.idx[[mm]]
        x.user.idx <- object@Model@x.user.idx[[mm]]
        GLIST[[mm]][m.user.idx] <- object@Fit@est[x.user.idx]

        # class
        if(object@Model@isSymmetric[mm]) {
            class(GLIST[[mm]]) <- c("lavaan.matrix.symmetric", "matrix")
        } else {
            class(GLIST[[mm]]) <- c("lavaan.matrix", "matrix")
        }
    }

    GLIST
}

parameter.list <- function(object) {

    LIST <- as.data.frame(object@ParTable, stringsAsFactors = FALSE)
    LIST
}


derivatives <- function(object) {
 
    GLIST <- computeGradient(object@Model, GLIST=NULL, 
                             samplestats=object@SampleStats, type="allofthem",
                             estimator=object@Options$estimator, 
                             verbose=FALSE, forcePD=TRUE,
                             group.weight=TRUE, constraints=FALSE)
    names(GLIST) <- names(object@Model@GLIST)

    for(mm in 1:length(GLIST)) {
        # labels
        dimnames(GLIST[[mm]]) <- object@Model@dimNames[[mm]]

        # class
        if(object@Model@isSymmetric[mm]) {
            class(GLIST[[mm]]) <- c("lavaan.matrix.symmetric", "matrix")
        } else {
            class(GLIST[[mm]]) <- c("lavaan.matrix", "matrix")
        }
    }

    GLIST
}

setMethod("coef", "lavaan",
function(object, type="free", labels=TRUE) {

    if(type == "user" || type == "all") {
        type <- "user"
        idx <- 1:length( object@ParTable$lhs )
    } else if(type == "free") {
        idx <- which(object@ParTable$free > 0L & !duplicated(object@ParTable$free))
    } else if(type == "unco") {
        idx <- which(object@ParTable$unco > 0L & 
                     !duplicated(object@ParTable$unco))
    } else {
        stop("argument `type' must be one of free, unco, or user")
    }
    cof <- object@Fit@est[idx]
  
    # labels?
    if(labels) names(cof) <- getParameterLabels(object@ParTable, type=type)

    # class
    class(cof) <- c("lavaan.vector", "numeric")

    cof
})

standardizedSolution <- standardizedsolution <- function(object, type="std.all") {

    stopifnot(type %in% c("std.all", "std.lv", "std.nox"))

    LIST <- inspect(object, "list")
    LIST <- LIST[,c("lhs", "op", "rhs", "group")]

    # add std and std.all columns
    if(type == "std.lv") {
        LIST$est.std     <- standardize.est.lv(object)
    } else if(type == "std.all") {
        LIST$est.std <- standardize.est.all(object)
    } else if(type == "std.nox") {
        LIST$est.std <- standardize.est.all.nox(object)
    }

    # add 'se' for standardized parameters
    # TODO!!
    LIST$se <- rep(NA, length(LIST$lhs))

    # add 'z' column
    tmp.se <- ifelse( LIST$se == 0.0, NA, LIST$se)
    LIST$z <- LIST$est.std / tmp.se
    LIST$pvalue <- 2 * (1 - pnorm( abs(LIST$z) ))

    # if single group, remove group column
    if(object@Data@ngroups == 1L) LIST$group <- NULL

    class(LIST) <- c("lavaan.data.frame", "data.frame")
    LIST
}

parameterEstimates <- parameterestimates <- 
    function(object, 
             ci = TRUE, level = 0.95, boot.ci.type = "perc",
             standardized = FALSE, 
             fmi = "default") {

    # fmi or not
    FMI <- fmi
    if(fmi == "default") {
        if(object@SampleStats@missing.flag &&
           object@Fit@converged &&
           object@Options$estimator == "ML" &&
           object@Options$se == "standard")
            FMI <- TRUE
        else
            FMI <- FALSE
    }
    if(FMI && object@Options$se != "standard") {
        warning("lavaan WARNING: fmi only available if se=\"standard\"")
        FMI <- FALSE
    }

    LIST <- inspect(object, "list")
    LIST <- LIST[,c("lhs", "op", "rhs", "group", "label")]
    # add est and se column
    est <- object@Fit@est
    BOOT <- attr(est, "BOOT.COEF")
    attributes(est) <- NULL
    LIST$est <- est

    if(object@Options$se != "none") {
        LIST$se  <- object@Fit@se
        tmp.se <- ifelse( LIST$se == 0.0, NA, LIST$se)
        LIST$z <- LIST$est / tmp.se
        LIST$pvalue <- 2 * (1 - pnorm( abs(LIST$z) ))
    }

    # confidence interval
    if(object@Options$se != "none" && ci) {
        # next three lines based on confint.lm
        a <- (1 - level)/2; a <- c(a, 1 - a)
        if(object@Options$se != "bootstrap") {
            fac <- qnorm(a)
            ci <- LIST$est + LIST$se %o% fac
        } else if(object@Options$se == "bootstrap") {
            stopifnot(!is.null(BOOT))
            stopifnot(boot.ci.type %in% c("norm","basic","perc","bca.simple"))
            if(boot.ci.type == "norm") {
                fac <- qnorm(a)
                boot.x <- colMeans(BOOT)
                boot.est <- 
                    getModelParameters(object@Model, 
                                       GLIST=x2GLIST(object@Model, boot.x), 
                                       type="user", extra=TRUE)
                bias.est <- (boot.est - LIST$est)
                ci <- (LIST$est - bias.est) + LIST$se %o% fac
            } else if(boot.ci.type == "basic") {
                ci <- cbind(LIST$est, LIST$est)
                alpha <- (1 + c(level, -level))/2

                # free.idx only
                qq <- apply(BOOT, 2, boot:::norm.inter, alpha)
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
                    qq <- apply(BOOT.def, 2, boot:::norm.inter, alpha)
                    ci[def.idx,] <- 2*ci[def.idx,] - t(qq[c(3,4),])
                }

                # TODO: add cin/ceq?
               
            } else if(boot.ci.type == "perc") {
                ci <- cbind(LIST$est, LIST$est)
                alpha <- (1 + c(-level, level))/2

                # free.idx only
                qq <- apply(BOOT, 2, boot:::norm.inter, alpha)
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
                    qq <- apply(BOOT.def, 2, boot:::norm.inter, alpha)
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
                   qq <- boot:::norm.inter(t, adj.alpha)
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
                       qq <- boot:::norm.inter(t, adj.alpha)
                       ci[def.idx[i],] <- qq[,2]
                   }
               }
               
               # TODO:
               # - add cin/ceq
            }
        }

        LIST$ci.lower <- ci[,1]; LIST$ci.upper <- ci[,2]    
    }

    # add std and std.all columns
    if(standardized) {
        LIST$std.lv  <- standardize.est.lv(object)
        LIST$std.all <- standardize.est.all(object, est.std=LIST$est.std)
        LIST$std.nox <- standardize.est.all.nox(object, est.std=LIST$est.std)
    }

    # fractional missing information (if estimator="fiml")
    if(FMI) {
        SE.orig <- LIST$se
        COV <- object@Fit@Sigma.hat
        MEAN <- object@Fit@Mu.hat

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
        SE.step2 <- ifelse(step2@Fit@se == 0.0, as.numeric(NA), step2@Fit@se)
        LIST$fmi <- 1-(SE.step2^2/SE.orig^2)
    }

    # if single group, remove group column
    if(object@Data@ngroups == 1L) LIST$group <- NULL

    # if no user-defined labels, remove label column
    if(sum(nchar(object@ParTable$label)) == 0L) LIST$label <- NULL


    class(LIST) <- c("lavaan.data.frame", "data.frame")
    LIST
}

parameterTable <- parametertable <- parTable <- partable <-
        function(object) {
    inspect(object, "list")            
}

varTable <- vartable <- function(object, ov.names=names(object), 
                                 ov.names.x=NULL, as.data.frame.=TRUE) {

    if(class(object) == "lavaan") {
        VAR <- object@Data@ov
    } else if(class(object) == "data.frame") {
        OV <- lapply(object[,unique(unlist(c(ov.names,ov.names.x))),drop=FALSE],
                     function(x)
                  list(nobs=sum(!is.na(x)),
                       type=class(x)[1],
                       mean=ifelse(class(x)[1] == "numeric",
                                   mean(x, na.rm=TRUE), as.numeric(NA)),
                       var=ifelse(class(x)[1] == "numeric",
                                  var(x, na.rm=TRUE), as.numeric(NA)),
                       nlevels=nlevels(x),
                       lnames=paste(levels(x),collapse="|")
                      ))
        VAR <- list()
        VAR$name    <- names(OV)
        VAR$idx  <- match(VAR$name, names(object))
        VAR$nobs <- unname(sapply(OV, "[[", "nobs"))
        VAR$type <- unname(sapply(OV, "[[", "type"))
        VAR$mean <- unname(sapply(OV, "[[", "mean"))
        VAR$var  <- unname(sapply(OV, "[[", "var"))
        VAR$nlev <- unname(sapply(OV, "[[", "nlevels"))
        VAR$lnam <- unname(sapply(OV, "[[", "lnames"))
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



setMethod("inspect", "lavaan",
function(object, what="free") {

    if(length(what) > 1) {
        stop("`what' arguments contains multiple arguments; only one is allowed")
    }

    # be case insensitive
    what <- tolower(what)

    if(what == "free") {
        free.parameters(object, type="free")
    } else if(what == "unco") {
        free.parameters(object, type="unco")
    } else if(what == "user") {
        free.parameters(object, type="user")
    } else if(what == "se" ||
              what == "std.err" ||
              what == "standard.errors") {
        standard.errors(object)
    } else if(what == "start" ||
              what == "starting.values") {
        starting.values(object)
    } else if(what == "coef" ||
              what == "coefficients" ||
              what == "parameters" ||
              what == "parameter.estimates" ||
              what == "parameter.values" ||
              what == "estimates" ||
              what == "x" ||
              what == "est") {
        parameter.values(object)
    } else if(what == "list") {
        parameter.list(object)
    } else if(what == "dx" ||
              what == "gradient" ||
              what == "derivatives") {
        derivatives(object)
    } else if(what == "fit" ||
              what == "fitmeasures" ||
              what == "fit.measures" ||
              what == "fit.indices") {
        fitMeasures(object)
    } else if(what == "std" ||
              what == "std.coef" ||
              what == "standardized" ||
              what == "standardizedsolution" ||
              what == "standardized.solution") {
        standardizedSolution(object)
    } else if(what == "mi" ||
              what == "modindices" ||
              what == "modification" ||
              what == "modificationindices" ||
              what == "modification.indices") {
         modificationIndices(object)
    } else if(what == "samp" ||
              what == "sample" ||
              what == "samplestatistics" ||
              what == "sampstat") {
        sampStat(object)
    } else if(what == "rsquare" || 
              what == "r-square" ||
              what == "r2") {
         rsquare(object)
    } else if(what == "wls.v") {
        getWLS.V(object, drop.list.single.group=TRUE)
    } else if(what == "nacov") {
        getSampleStatsNACOV(object)    
    } else if(what == "modelcovlv"  ||
              what == "modelcov.lv" ||
              what == "cov.lv") {
        getModelCovLV(object, correlation.metric=FALSE, labels=TRUE)
    } else if(what == "modelcorlv"  ||
              what == "modelcor.lv" ||
              what == "cor.lv") {
        getModelCorLV(object, labels=TRUE)
    } else if(what == "modelcovall"  ||
              what == "modelcov.all" ||
              what == "cov.all") {
        getModelCov(object, correlation.metric=FALSE, labels=TRUE)
    } else if(what == "modelcorall"  ||
              what == "modelcor.all" ||
              what == "cor.all") {
        getModelCor(object, labels=TRUE)
    } else if(what == "converged") {
        object@Fit@converged
    } else {
        stop("unknown `what' argument in inspect function: `", what, "'")
    }

})

# fixme, should we export this function?
sampStat <- function(object, labels=TRUE) {

    G <- object@Data@ngroups
    ov.names <- object@Data@ov.names

    OUT <- vector("list", length=G)
    for(g in 1:G) {
        OUT[[g]]$cov  <- object@SampleStats@cov[[g]]
        if(labels) 
            rownames(OUT[[g]]$cov) <- colnames(OUT[[g]]$cov) <- ov.names[[g]]
        class(OUT[[g]]$cov) <- c("lavaan.matrix.symmetric", "matrix")

        #if(object@Model@meanstructure) {
            OUT[[g]]$mean <- as.numeric(object@SampleStats@mean[[g]])
            if(labels) names(OUT[[g]]$mean) <- ov.names[[g]]
            class(OUT[[g]]$mean) <- c("lavaan.vector", "numeric")
        #}

        if(object@Model@categorical) {
            OUT[[g]]$th <- as.numeric(object@SampleStats@th[[g]])
            if(length(object@Model@num.idx[[g]]) > 0L) {
                OUT[[g]]$th <- OUT[[g]]$th[-object@Model@num.idx[[g]]]
            }
            if(labels) {
                names(OUT[[g]]$th) <- 
                    vnames(object@ParTable, type="th", group=g)
            }
            class(OUT[[g]]$th) <- c("lavaan.vector", "numeric")
        }

        if(object@Model@categorical &&
           object@Model@nexo > 0L) {
            OUT[[g]]$slopes  <- object@SampleStats@slopes[[g]]
            if(labels) {
                rownames(OUT[[g]]$slopes) <- ov.names[[g]]
                colnames(OUT[[g]]$slopes) <- 
                    vnames(object@ParTable, type="ov.x", group=g)
                class(OUT[[g]]$slopes) <- c("lavaan.matrix", "matrix")
            }
        }
    }

    if(G == 1) {
        OUT <- OUT[[1]]
    } else {
        names(OUT) <- unlist(object@Data@group.label)
    }

    OUT
}


setMethod("fitted.values", "lavaan",
function(object, labels=TRUE) {

    G <- object@Data@ngroups
    ov.names <- object@Data@ov.names

    OUT <- vector("list", length=G)
    for(g in 1:G) {
        OUT[[g]]$cov  <- object@Fit@Sigma.hat[[g]]
        if(labels) 
            rownames(OUT[[g]]$cov) <- colnames(OUT[[g]]$cov) <- ov.names[[g]]
        class(OUT[[g]]$cov) <- c("lavaan.matrix.symmetric", "matrix")

        #if(object@Model@meanstructure) {
            OUT[[g]]$mean <- as.numeric(object@Fit@Mu.hat[[g]])
            if(labels) names(OUT[[g]]$mean) <- ov.names[[g]]
            class(OUT[[g]]$mean) <- c("lavaan.vector", "numeric")
        #}

        if(object@Model@categorical) {
            OUT[[g]]$th <- as.numeric(object@Fit@TH[[g]])
            if(length(object@Model@num.idx[[g]]) > 0L) {
                OUT[[g]]$th <- OUT[[g]]$th[-object@Model@num.idx[[g]]]
            }
            if(labels) {
                names(OUT[[g]]$th) <-
                    vnames(object@ParTable, type="th", group=g)
            }
            class(OUT[[g]]$th) <- c("lavaan.vector", "numeric")
        }
    }

    if(G == 1) {
        OUT <- OUT[[1]]
    } else {
        names(OUT) <- unlist(object@Data@group.label)
    }

    OUT
})


setMethod("fitted", "lavaan",
function(object, labels=TRUE) {
     fitted.values(object, labels=labels)
})


setMethod("vcov", "lavaan",
function(object, labels=TRUE) {

    # check for convergence first!
    if(object@Fit@npar > 0L && !object@Fit@converged)
        stop("lavaan ERROR: model did not converge")

    if(object@Fit@npar == 0) {
        VarCov <- matrix(0,0,0)
    } else {
        VarCov <- estimateVCOV(object@Model, samplestats=object@SampleStats, 
                               options=object@Options,
                               #data=eval(object@call[["data"]], 
                               #          parent.frame()) 
                               data=object@Data
                              )
    }

    if(labels) {
        colnames(VarCov) <- rownames(VarCov) <- 
            getParameterLabels(object@ParTable, type="free")
    }

    class(VarCov) <- c("lavaan.matrix.symmetric", "matrix")
    VarCov
})


rsquare <- function(object, est.std.all=NULL) {

    if(is.null(est.std.all)) est.std.all <- standardize.est.all(object)
    ngroups <- object@Data@ngroups
    partable <- object@ParTable
    partable$rsquare <- 1.0 - est.std.all
    # no values > 1.0
    partable$rsquare[partable$rsquare > 1.0] <- as.numeric(NA)
    r2 <- vector("list", length=ngroups)

    for(g in 1:ngroups) {
        ind.names   <- partable$rhs[ which(partable$op == "=~" & partable$group == g) ]
        eqs.y.names <- partable$lhs[ which(partable$op == "~"  & partable$group == g) ]
        y.names <- unique( c(ind.names, eqs.y.names) )

        idx <- which(partable$op == "~~" & partable$lhs %in% y.names & 
                     partable$rhs == partable$lhs & partable$group == g)
        tmp <- partable$rsquare[idx]; names(tmp) <- partable$lhs[idx]
        r2[[g]] <- tmp
    }

    if(ngroups == 1) {
        r2 <- r2[[1]]
    } else {
        names(r2) <- unlist(object@Data@group.label)
    }

    r2
}

# logLik (so that we can use the default AIC/BIC functions from stats4(
setMethod("logLik", "lavaan",
function(object, ...) {
    if(object@Options$estimator != "ML") {
        stop("lavaan ERROR: logLik only available if estimator is ML")
    }
    if(object@Fit@npar > 0L && !object@Fit@converged)
        stop("lavaan ERROR: model did not converge")
    
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



# this is based on the anova function in the lmer package
setMethod("anova", signature(object = "lavaan"),
function(object, ...) {

    if(object@Fit@npar > 0L && !object@Fit@converged)
        stop("lavaan ERROR: model did not converge")

    # NOTE: if we add additional arguments, it is not the same generic
    # anova() function anymore, and match.call will be screwed up

    # default arguments
    SB.classic <- FALSE
    SB.H0      <- FALSE

    mcall <- match.call(expand.dots = TRUE)
    dots <- list(...)
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

    # some general properties (taken from the first model)
    estimator <- object@Options$estimator
    likelihood <- object@Options$likelihood
    ngroups <- object@Data@ngroups
    nobs <- object@SampleStats@nobs
    ntotal <- object@SampleStats@ntotal
    npar <- object@Fit@npar

    # shortcut for single argument (just plain LRT)
    if(!any(modp)) {
        aic <- bic <- c(NA, NA)
        if(estimator == "ML") {
            aic <- c(NA, AIC(object))
            bic <- c(NA, BIC(object))
        }

        val <- data.frame(Df = c(0, object@Fit@test[[1L]]$df),
                          AIC = aic,
                          BIC = bic,
                          Chisq = c(0, object@Fit@test[[1L]]$stat),
                          "Chisq diff" = c(NA, object@Fit@test[[1L]]$stat),
                          "Df diff" = c(NA, object@Fit@test[[1L]]$df),
                          "Pr(>Chisq)" = c(NA, object@Fit@test[[1L]]$pvalue),
                          row.names = c("Saturated", "Model"),
                          check.names = FALSE)
        attr(val, "heading") <- "Chi Square Test Statistic (unscaled)\n"
        class(val) <- c("anova", class(val))
        return(val)
    }

    # list of models
    mods <- c(list(object), dots[modp])
    names(mods) <- sapply(as.list(mcall)[c(FALSE, TRUE, modp)], as.character)

    # put them in order (using number of free parameters)
    nfreepar <- sapply(mods, function(x) x@Fit@npar)
    if(any(duplicated(nfreepar))) { ## FIXME: what to do here?
        # what, same number of free parameters?
        # maybe, we need to count number of constraints
        ncon <- sapply(mods, function(x) { nrow(x@Model@con.jac) })
        nfreepar <- nfreepar - ncon
    }

    mods <- mods[order(nfreepar, decreasing = TRUE)]

    # here come the checks
    if(TRUE) {
        # 1. same data? we check the covariance matrices of the first group
        ## wow FIXME: we may need to reorder the rows/columns first!!
        #COVS <- lapply(mods, function(x) slot(slot(x, "Sample"), "cov")[[1]])
        #if(!all(sapply(COVS, all.equal, COVS[[1]]))) {
        #    stop("lavaan ERROR: models must be fit to the same data")
        #}
        # 2. nested models? *different* npars?
     
        # TODO!
        
        # 3. all meanstructure?
    }

    # 
    mods.scaled <- unlist( lapply(mods, function(x) {
        any(c("satorra.bentler", "yuan.bentler", 
              "mean.var.adjusted", "scaled.shifted") %in% 
            unlist(sapply(slot(slot(x, "Fit"), "test"), "[", "test")) ) }))

    if(all(mods.scaled)) {
        scaled <- TRUE
        # which type?
        TEST <- object@Fit@test[[2]]$test
    } else if(!all(mods.scaled)) {
        scaled <- FALSE
        TEST <- "standard"
    } else {
        stop("lavaan ERROR: some models (but not all) have scaled test statistics")
    }

    # which models have used a MEANSTRUCTURE?
    mods.meanstructure <- sapply(mods, function(x) { 
                                 unlist(slot(slot(x, "Model"), 
                                                     "meanstructure"))})
    if(all(mods.meanstructure)) {
        meanstructure <- "ok"
    } else if(sum(mods.meanstructure) == 0) {
        meanstructure <- "ok"
    } else {
        stop("lavaan ERROR: some models (but not all) have a meanstructure")
    }

    # collect statistics for each model
    Df    <- unlist(lapply(mods, function(x) slot(slot(x, "Fit"), 
                        "test")[[1]]$df))
    Chisq <- unlist(lapply(mods, function(x) slot(slot(x, "Fit"), 
                        "test")[[1]]$stat))

    # difference statistics
    Chisq.delta  <- c(NA, diff(Chisq)) 
    Df.delta     <- c(NA, diff(Df))
    
    # correction for scaled test statistics
    if(scaled) {
        if(SB.classic && TEST %in% c("satorra.bentler", "yuan.bentler")) {
            # use formula from Satorra & Bentler 2001
            scaling.factor <- unlist(lapply(mods, 
                function(x) slot(slot(x, "Fit"), "test")[[2]]$scaling.factor))
            cd1 <- diff(scaling.factor * Df)/diff(Df)

            # check for negative scaling factors
            if(any(cd1 < 0)) {
                warning("lavaan WARNING: some scaling factors are negative: [",
                        paste(round(cd1, 3), collapse=" "),"]; rerun with SB.classic=FALSE")
                cd1[cd1 < 0] <- NA
            }
            cd <- c(NA, cd1)
            Chisq.delta <- Chisq.delta/cd

            # extract scaled Chisq for each model
            Chisq <- unlist(lapply(mods, function(x) slot(slot(x, "Fit"),
                            "test")[[2]]$stat))
        } else {
            # see Mplus Web Note 10 (2006)
            for(m in seq_len(length(mods) - 1L)) {

                if(mods[[m]]@Fit@test[[1]]$df == mods[[m+1]]@Fit@test[[1]]$df) {
                    warnings("lavaan WARNING: some models have the same number of free parameters")
                    next
                }

                if(SB.H0) {
                    # evaluate under H0
                    stop("SB.H0 has not been implemented yet. FIXME!")
                }

                # original M (Satorra)
                Delta1 <- computeDelta(mods[[m]]@Model)
                WLS.V <- getWLS.V(object)
                Gamma <- getSampleStatsNACOV(object)

                # weight WLS.V
                for(g in 1:ngroups) {
                    WLS.V[[g]] <- nobs[[g]]/ntotal * WLS.V[[g]]
                }

                # information matrix
                P1 <- matrix(0, nrow=npar, ncol=npar)
                for(g in 1:ngroups) {
                     P1 <- P1 + (t(Delta1[[g]]) %*% WLS.V[[g]] %*% Delta1[[g]])
                }
                P1.inv <- solve(P1)

                # compute A for these two nested models
                p1 <- mods[[m   ]]@ParTable # partable h1
                p0 <- mods[[m+1L]]@ParTable # partable h0
                af <- getConstraintsFunction(p1,p0)
                A <- lavJacobianC(func=af, x=mods[[m   ]]@Fit@x)

                trace.UGamma  <- numeric( ngroups )
                trace.UGamma2 <- numeric( ngroups )
                for(g in 1:ngroups) {
                    U <- WLS.V[[g]]  %*% Delta1[[g]] %*%
                         (P1.inv %*% t(A) %*% solve(A %*% P1.inv %*% t(A)) %*% A %*% P1.inv) %*% t(Delta1[[g]]) %*% WLS.V[[g]]
                    UG <- U %*% Gamma[[g]]; tUG <- t(UG)
                    trace.UGamma[g]  <- ntotal/nobs[[g]] * sum( U * Gamma[[g]] )
                    trace.UGamma2[g] <- ntotal/nobs[[g]] * sum( UG * tUG )
                }

                tr.M  <- sum(trace.UGamma)
                tr2.M <- sum(trace.UGamma2)

                # adjust Df.delta?
                if(TEST == "mean.var.adjusted") {
                    # NOT needed for scaled.shifted; see
                    # see 'Simple Second Order Chi-Square Correction' 2010 paper
                    # on www.statmodel.com, section 4
                    # Df.delta[m+1L] <- floor((tr.M1^2 / tr2.M1) + 0.5)
                    Df.delta[m+1L] <- tr.M^2 / tr2.M
                } 

                scaling.factor <- tr.M / Df.delta[m+1L]
                if(scaling.factor < 0) scaling.factor <- as.numeric(NA)
                Chisq.delta[m+1L] <- Chisq.delta[m+1L]/scaling.factor
            }
        } 
    }

    # Pvalue
    Pvalue.delta <- pchisq(Chisq.delta, Df.delta, lower.tail = FALSE)

    aic <- bic <- rep(NA, length(mods))
    if(estimator == "ML") {
        aic <- sapply(mods, FUN=AIC)
        bic <- sapply(mods, FUN=BIC)
    }

    val <- data.frame(Df = Df,
                      AIC = aic,
                      BIC = bic,
                      Chisq = Chisq,
                      "Chisq diff" = Chisq.delta,
                      "Df diff" = Df.delta,
                      "Pr(>Chisq)" = Pvalue.delta,
                      row.names = names(mods), 
                      check.names = FALSE)

    if(scaled) {
        attr(val, "heading") <- 
            paste("Scaled Chi Square Difference Test (test = ",
                  TEST, ")\n", sep="")
    } else {
        attr(val, "heading") <- "Chi Square Difference Test\n"
    }
    class(val) <- c("anova", class(val))

    return(val)

})


getWLS.est <- function(object, drop.list.single.group=FALSE) {

    # shortcuts
    samplestats   = object@SampleStats
    estimator     = object@Options$estimator
    G             = object@Data@ngroups
    categorical   = object@Model@categorical
    meanstructure = object@Model@meanstructure
    fixed.x       = object@Model@fixed.x

    # compute moments for all groups
    Sigma.hat <- computeSigmaHat(object@Model)
    if(meanstructure && !categorical) {
        Mu.hat <- computeMuHat(object@Model)
    } else if(categorical) {
        TH <- computeTH(object@Model)
        if(fixed.x)
            PI <- computePI(object@Model)
    }

    WLS.est <- vector("list", length=samplestats@ngroups)
    for(g in 1:samplestats@ngroups) {
        if(categorical) {
            if(fixed.x) {
                WLS.est[[g]] <- c(TH[[g]],vec(PI[[g]]),
                                  diag(Sigma.hat[[g]])[num.idx[[g]]],
                                  vech(Sigma.hat[[g]], diagonal=FALSE))
            } else {
                WLS.est[[g]] <- c(TH[[g]],diag(Sigma.hat[[g]])[num.idx[[g]]],
                                  vech(Sigma.hat[[g]], diagonal=FALSE))
            }
        } else if(meanstructure) {
            WLS.est[[g]] <- c(Mu.hat[[g]], vech(Sigma.hat[[g]]))
        } else {
            WLS.est[[g]] <- vech(Sigma.hat[[g]])
        }
    }

    OUT <- WLS.est
    if(G == 1 && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(G > 1) names(OUT) <- unlist(object@Data@group.label)
    }

    OUT
}

getWLS.V <- function(object, Delta=computeDelta(object@Model), 
                     drop.list.single.group=FALSE) {

    # shortcuts
    samplestats = object@SampleStats
    estimator   = object@Options$estimator
    G           = object@Data@ngroups

    WLS.V       <- vector("list", length=samplestats@ngroups)
    if(estimator == "GLS"  ||
       estimator == "WLS"  ||
       estimator == "DWLS" ||
       estimator == "ULS") {
        # for GLS, the WLS.V22 part is: 0.5 * t(D) %*% [S.inv %x% S.inv] %*% D
        # for WLS, the WLS.V22 part is: Gamma
        WLS.V <- samplestats@WLS.V
    } else if(estimator == "ML") {
        Sigma.hat <- computeSigmaHat(object@Model)
        if(object@Model@meanstructure) Mu.hat <- computeMuHat(object@Model)
        for(g in 1:samplestats@ngroups) {
            if(samplestats@missing.flag) {
                WLS.V[[g]] <- compute.Abeta(Sigma.hat=Sigma.hat[[g]],
                                            Mu.hat=Mu.hat[[g]],
                                            samplestats=samplestats,
                                            data=object@data, group=g,
                                            information="expected")
            } else {
                # WLS.V22 = 0.5*t(D) %*% [Sigma.hat.inv %x% Sigma.hat.inv]%*% D
                WLS.V[[g]] <-
                    compute.Abeta.complete(Sigma.hat=Sigma.hat[[g]],
                                           meanstructure=object@Model@meanstructure)
            }
        }
    } # ML


    OUT <- WLS.V
    if(G == 1 && drop.list.single.group) {
        OUT <- OUT[[1]]
    } else {
        if(G > 1) names(OUT) <- unlist(object@Data@group.label)
    }

    OUT
}

getSampleStatsNACOV <- function(object) {

    if(object@Options$se == "robust.mlr")
        stop("not done yet; FIX THIS!")

    # shortcuts
    samplestats = object@SampleStats
    estimator   = object@Options$estimator

    NACOV       <- vector("list", length=samplestats@ngroups)

    if(estimator == "GLS"  ||
       estimator == "WLS"  ||
       estimator == "DWLS" ||
       estimator == "ULS") {
        NACOV <- samplestats@NACOV
    } else if(estimator == "ML") {
        for(g in 1:samplestats@ngroups) {
            NACOV[[g]] <- 
                compute.Gamma(object@Data@X[[g]], 
                              meanstructure=object@Options$meanstructure)
            if(object@Options$mimic == "Mplus") {
                G11 <- ( samplestats@cov[[g]] * (samplestats@nobs[[g]]-1)/samplestats@nobs[[g]] )
                NACOV[[g]][1:nrow(G11), 1:nrow(G11)] <- G11
            }
        }
    }

    NACOV
}

getHessian <- function(object) {
    # lazy approach: take -1 the observed information
    E <- computeObservedInformation(object@Model, 
                                    samplestats=object@SampleStats,
                                    X=object@Data@X,
                                    type="free",
                                    estimator=object@Options$estimator,
                                    group.weight=TRUE)

    -E
}

getVariability <- function(object) {
    # lazy approach: get it from Nvcov.first.order
    NACOV <- Nvcov.first.order(object@Model,
                               samplestats=object@SampleStats,
                               data=object@Data,
                               estimator=object@Options$estimator)

    B0 <- attr(NACOV, "B0")

    if(object@Options$estimator == "PML") {
        B0 <- B0 * object@SampleStats@ntotal
    }
    
    B0
}

getModelCorLV <- function(object, labels=TRUE) {
    getModelCovLV(object, correlation.metric=TRUE, labels=labels)
}

getModelCovLV <- function(object, correlation.metric=FALSE, labels=TRUE) {

    G <- object@Data@ngroups

    # compute lv covar
    OUT <- lavaan:::computeETA(object@Model, samplestats=fit@SampleStats)

    # correlation?
    if(correlation.metric) {
        OUT <- lapply(OUT, cov2cor)
    }
    
    # we need psi matrix for labels
    psi.group <- which(names(object@Model@GLIST) == "psi")
    if(labels) {
        for(g in 1:G) {
            psi.idx <- psi.group[g]
            NAMES <- object@Model@dimNames[[psi.idx]][[1L]]
            colnames(OUT[[g]]) <- rownames(OUT[[g]]) <- NAMES
            class(OUT[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
    }

    if(G == 1) {
        OUT <- OUT[[1]]
    } else {
        names(OUT) <- unlist(object@Data@group.label)
    }

    OUT
}

getModelCor <- function(object, labels=TRUE) {
    getModelCov(object, correlation.metric=TRUE, labels=labels)
}

getModelCov <- function(object, correlation.metric=FALSE, labels=TRUE) {

    G <- object@Data@ngroups

    # compute extended model implied covariance matrix (both ov and lv)
    OUT <- lavaan:::computeCOV(object@Model, samplestats=fit@SampleStats)

    # correlation?
    if(correlation.metric) {
        OUT <- lapply(OUT, cov2cor)
    }

    # we need lambda + psi matrix for labels
    lambda.group <- which(names(object@Model@GLIST) == "lambda")
    psi.group <- which(names(object@Model@GLIST) == "psi")
    if(labels) {
        for(g in 1:G) {
            lambda.idx <- lambda.group[g]
            psi.idx <- psi.group[g]
            NAMES <- c(object@Model@dimNames[[lambda.idx]][[1L]],
                       object@Model@dimNames[[psi.idx]][[1L]])
            colnames(OUT[[g]]) <- rownames(OUT[[g]]) <- NAMES
            class(OUT[[g]]) <- c("lavaan.matrix.symmetric", "matrix")
        }
    }

    if(G == 1) {
        OUT <- OUT[[1]]
    } else {
        names(OUT) <- unlist(object@Data@group.label)
    }

    OUT
}

