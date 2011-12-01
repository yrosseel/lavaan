#
# initial version: YR 25/03/2009

short.summary <- function(object) {

    # Convergence or not?
    if(object@Fit@iterations > 0) {
        if(object@Fit@converged) {
	    cat(sprintf("lavaan (%s) converged normally after %i iterations\n",
                        packageDescription("lavaan", fields="Version"),
                        object@Fit@iterations))
        } else {
            cat("\n")
            cat(sprintf("** WARNING ** Lavaan (%s) did NOT converge after %i iterations\n", 
                packageDescription("lavaan", fields="Version"),
                object@Fit@iterations))
            cat("** WARNING ** Estimates below are most likely unreliable\n")
        }
    } else {
        cat("\n")
        cat(sprintf("Lavaan (%s)\n",
                    packageDescription("lavaan", fields="Version")))
        cat("** WARNING ** Model has not been fitted\n")
        cat("** WARNING ** Estimates below are simply the starting values\n")
    }
    cat("\n")
   
    # Short Model description + number of observations
    #h.txt <- sprintf("%-40s%11s  %9s", "", "Independent", "Dependent")
    #cat(h.txt,"\n")
    #t0.txt <- sprintf("  %-38s", "Number of observed variables")
    #t1.txt <- sprintf("  %9i", length( object@User@ov.x.idx ))
    #t2.txt <- sprintf("  %9i", length( object@User@ov.mm.idx ))
    #cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

    #t0.txt <- sprintf("  %-38s", "Number of continuous latent variables")
    #t1.txt <- sprintf("  %9i", length( object@User@lv.x.idx ))
    #t2.txt <- sprintf("  %9i", length( object@User@lv.y.idx ))
    #cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
    #cat("\n")

    #t0.txt <- sprintf("  %-38s", "Number of free parameters")
    #t1.txt <- sprintf("  %9i", object@User@npar )
    #cat(t0.txt, t1.txt, "\n", sep="")

    #t0.txt <- sprintf("  %-38s", "Number of groups")
    #t1.txt <- sprintf("  %9i", object@Sample@ngroups)
    #cat(t0.txt, t1.txt, "\n", sep="")

    # listwise deletion?
    listwise <- FALSE
    if(!any(unlist(lapply(object@Sample@missing, "[[", "flag"))) &&
       (object@Sample@nobs[[1L]] < object@Sample@norig[[1L]])) {
        listwise <- TRUE
    }

    if(object@Sample@ngroups == 1L) {
        if(listwise) {
            cat(sprintf("  %-40s", ""), sprintf("  %10s", "Used"), 
                                        sprintf("  %10s", "Total"),
                "\n", sep="")
        }
        t0.txt <- sprintf("  %-40s", "Number of observations")
        t1.txt <- sprintf("  %10i", object@Sample@nobs[[1L]])
        t2.txt <- ifelse(listwise,
                  sprintf("  %10i", object@Sample@norig[[1L]]), "")
        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
    } else {
        if(listwise) {
            cat(sprintf("  %-40s", ""), sprintf("  %10s", "Used"),  
                                        sprintf("  %10s", "Total"),
                "\n", sep="")
        }
        t0.txt <- sprintf("  %-40s", "Number of observations per group")
        cat(t0.txt, "\n")
        for(g in 1:object@Sample@ngroups) {
            t.txt <- sprintf("  %-40s  %10i", object@Sample@group.label[[g]],
                                                object@Sample@nobs[[g]])
            t2.txt <- ifelse(listwise,
                      sprintf("  %10i", object@Sample@norig[[g]]), "")
            cat(t.txt, t2.txt, "\n", sep="")
        }
    }
    cat("\n")

    # missing patterns?
    if(object@Sample@missing[[1L]]$flag && object@Sample@ngroups == 1L) {
        t0.txt <- sprintf("  %-40s", "Number of missing patterns")
        t1.txt <- sprintf("  %10i", object@Sample@missing[[1]]$npatterns)
        cat(t0.txt, t1.txt, "\n\n", sep="")
    }
    if(any(unlist(lapply(object@Sample@missing, "[[", "flag"))) && 
       object@Sample@ngroups > 1L) {
        t0.txt <- sprintf("  %-40s", "Number of missing patterns per group")
        cat(t0.txt, "\n")
        for(g in 1:object@Sample@ngroups) {
            t.txt <- sprintf("  %-40s  %10i", object@Sample@group.label[[g]],
                             object@Sample@missing[[g]]$npatterns)
            cat(t.txt, "\n", sep="")
        }
        cat("\n")
    }

    # Print Chi-square value for the user-specified (full/h0) model

    # robust/scaled statistics?
    if(object@Options$test %in% c("satorra.bentler", "yuan.bentler")) {
        scaled <- TRUE
    } else {
        scaled <- FALSE
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
        t0.txt <- sprintf("  %-40s", "Minimum Function Chi-square")  
        t1.txt <- sprintf("  %10.3f", object@Fit@test[[1]]$stat)
        t2.txt <- ifelse(scaled, 
                  sprintf("  %10.3f", object@Fit@test[[2]]$stat), "")
        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

        # 2. degrees of freedom
        t0.txt <- sprintf("  %-40s", "Degrees of freedom")
        t1.txt <- sprintf("  %10i", object@Fit@test[[1]]$df); 
        t2.txt <- ifelse(scaled, 
                  sprintf("  %10i", object@Fit@test[[2]]$df), "")
        cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

        # 3. P-value
        t0.txt <- sprintf("  %-40s", "P-value")
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
                if(object@Options$mimic == "Mplus") {
                    cat("    for the Satorra-Bentler correction (Mplus variant)\n")
                } else {
                    cat("    for the Satorra-Bentler correction\n")
                }
            }
        } 

        if(object@Sample@ngroups > 1L) {
            cat("\n")
            cat("Chi-square for each group:\n\n")
            for(g in 1:object@Sample@ngroups) {
                t0.txt <- sprintf("  %-40s", object@Sample@group.label[[g]])
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
            cat("lavaan WARNING: fit measures not available if test = \"none\"\n\n")
        } else {
            print.fit.measures( fitMeasures(object) )
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
    t0.txt <- sprintf("  %-38s", "Standard Errors")
    tmp.txt <- object@Options$se
    t1.txt <- sprintf("  %12s", paste(toupper(substring(tmp.txt,1,1)),  
                                      substring(tmp.txt,2), sep=""))
    cat(t0.txt, t1.txt, "\n", sep="")
    if(object@Options$se == "bootstrap") {
        t0.txt <- sprintf("  %-40s", "Number of requested bootstrap draws")
        t1.txt <- sprintf("  %10i", object@Options$boot)
        cat(t0.txt, t1.txt, "\n", sep="")
        t0.txt <- sprintf("  %-40s", "Number of successful bootstrap draws")
        t1.txt <- sprintf("  %10i", nrow(attr(object@Fit@est, "BOOT")))
        cat(t0.txt, t1.txt, "\n", sep="")
    }
    cat("\n")

    # local print function
    print.estimate <- function(name="ERROR", i=1, z.stat=TRUE) {
       
        # cut name if (still) too long
        name <- substr(name, 1, 13)

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

    for(g in 1:object@Sample@ngroups) {
        ov.names <- vnames(object@User, "ov", group=g)
        lv.names <- vnames(object@User, "lv", group=g)

        # group header
        if(object@Sample@ngroups > 1) {
            if(g > 1) cat("\n\n")
            cat("Group ", g, 
                " [", object@Sample@group.label[[g]], "]:\n\n", sep="")
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
            # labels?
            l.idx <- which(nchar(LABELS) > 0L)
            if(length(l.idx) > 0L) {
                LABELS <- abbreviate(LABELS, 4)
                LABELS[l.idx] <- paste(" (", LABELS[l.idx], ")", sep="")
                MAX.L <- max(nchar(LABELS))
                NAMES <- abbreviate(NAMES, minlength = (13 - MAX.L), 
                                    strict = TRUE)
                NAMES <- sprintf(paste("%-", (13 - MAX.L), "s%", MAX.L, "s",
                                       sep=""), NAMES, LABELS)
            } else {
                NAMES <- abbreviate(NAMES, minlength = 13, strict = TRUE)
            }
        }

        NAMES <- object@User$rhs

        # 1. indicators ("=~") (we do show dummy indicators)
        mm.idx <- which( object@User$op == "=~" & 
                        !object@User$lhs %in% ov.names &
                         object@User$group == g)
        if(length(mm.idx)) {
            cat("Latent variables:\n")
            lhs.old <- ""
            NAMES[mm.idx] <- makeNames(  object@User$rhs[mm.idx],
                                       object@User$label[mm.idx])
            for(i in mm.idx) {
                lhs <- object@User$lhs[i]
                if(lhs != lhs.old) cat("  ", lhs, " =~\n", sep="")
                print.estimate(name=NAMES[i], i)
                lhs.old <- lhs
            }
            cat("\n")
        }

        # 2. regressions
        eqs.idx <- which(object@User$op == "~" & object@User$group == g)
        if(length(eqs.idx) > 0) {
            cat("Regressions:\n")
            lhs.old <- ""
            NAMES[eqs.idx] <- makeNames(  object@User$rhs[eqs.idx],
                                        object@User$label[eqs.idx])
            for(i in eqs.idx) {
                lhs <- object@User$lhs[i]
                if(lhs != lhs.old) cat("  ", lhs, " ~\n", sep="")
                print.estimate(name=NAMES[i], i)
                lhs.old <- lhs
            }
            cat("\n")
        }

        # 3. covariances
        cov.idx <- which(object@User$op == "~~" & 
                         !object@User$fixed.x &
                         object@User$lhs != object@User$rhs &
                         object@User$group == g)
        if(length(cov.idx) > 0) {
            cat("Covariances:\n")
            lhs.old <- ""
            NAMES[cov.idx] <- makeNames(  object@User$rhs[cov.idx],
                                        object@User$label[cov.idx])
            for(i in cov.idx) {
                lhs <- object@User$lhs[i]
                if(lhs != lhs.old) cat("  ", lhs, " ~~\n", sep="")
                print.estimate(name=NAMES[i], i)
                lhs.old <- lhs
            }
            cat("\n")
        }

        # 4. intercepts/means
        int.idx <- which(object@User$op == "~1" & 
                         !object@User$fixed.x &
                         object@User$group == g)
        if(length(int.idx) > 0) {
            cat("Intercepts:\n")
            NAMES[int.idx] <- makeNames(  object@User$lhs[int.idx],
                                        object@User$label[int.idx])
            for(i in int.idx) {
                print.estimate(name=NAMES[i], i)
            }
            cat("\n")
        }

        # 5. (residual) variances
        var.idx <- which(object@User$op == "~~" &
                         !object@User$fixed.x &
                         object@User$lhs == object@User$rhs &
                         object@User$group == g)
        if(length(var.idx) > 0) {
            cat("Variances:\n")
            NAMES[var.idx] <- makeNames(  object@User$rhs[var.idx],
                                        object@User$label[var.idx])
            for(i in var.idx) {
                if(object@Options$mimic == "lavaan") {
                    print.estimate(name=NAMES[i], i, z.stat=FALSE)
                } else {
                    print.estimate(name=NAMES[i], i, z.stat=TRUE)
                }
            }
            cat("\n")
        }

    } # ngroups

    # 6. variable definitions
    def.idx <- which(object@User$op == ":=")
    if(length(def.idx) > 0) {
        if(object@Sample@ngroups > 1) cat("\n")
        cat("Defined parameters:\n")
        NAMES[def.idx] <- makeNames(  object@User$lhs[def.idx], "")
        for(i in def.idx) {
            print.estimate(name=NAMES[i], i)
        }
        cat("\n")
    }

    # 7. constraints
    cin.idx <- which((object@User$op == "<" | 
                      object@User$op == ">"))
    ceq.idx <- which(object@User$op == "==")
    if(length(cin.idx) > 0L || length(ceq.idx) > 0L) {
        # set small negative values to zero, to avoid printing " -0.000"
        slack <- ifelse(abs(est) < 1e-5, 0, est)
        #slack[cin.idx] <- object@Model@cin.function(object@Fit@x)
        #slack[ceq.idx] <- object@Model@ceq.function(object@Fit@x)
       
        if(object@Sample@ngroups > 1 && length(def.idx) == 0L) cat("\n")
        cat("Constraints:                               Slack (>=0)\n")
        for(i in c(cin.idx,ceq.idx)) {
            lhs <- object@User$lhs[i]
             op <- object@User$op[i]
            rhs <- object@User$rhs[i]
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
        if(object@Sample@ngroups == 1L) r2 <- list(r2)
        for(g in 1:object@Sample@ngroups) {
            if(object@Sample@ngroups > 1) {
                cat("R-Square Group ", object@Sample@group.label[[g]], 
                    ":\n\n", sep="")       
            } else {
                cat("R-Square:\n\n")
            } 
            for(i in 1:length(r2[[g]])) {
                t1.txt <- sprintf("    %-13s %9.3f\n", names(r2[[g]])[i], 
                                  r2[[g]][i])
                cat(t1.txt)
            }
            if(g < object@Sample@ngroups) cat("\n")
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

    LIST <- as.data.frame(object@User, stringsAsFactors = FALSE)
    LIST
}


derivatives <- function(object) {
 
    GLIST <- computeGradient(object@Model, GLIST=NULL, 
                             sample=object@Sample, type="allofthem",
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
        idx <- 1:length( object@User$lhs )
    } else if(type == "free") {
        idx <- which(object@User$free > 0L & !duplicated(object@User$free))
    } else if(type == "unco") {
        idx <- which(object@User$free.uncon > 0L & 
                     !duplicated(object@User$free.uncon))
    } else {
        stop("argument `type' must be one of free, unco, or user")
    }
    cof <- object@Fit@est[idx]
  
    # labels?
    if(labels) names(cof) <- getParameterLabels(object@User, type=type)

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
    if(object@Sample@ngroups == 1L) LIST$group <- NULL

    class(LIST) <- c("lavaan.data.frame", "data.frame")
    LIST
}

parameterEstimates <- parameterestimates <- 
    function(object, 
             ci = TRUE, level = 0.95, boot.ci.type = "perc",
             standardized = FALSE) {

    LIST <- inspect(object, "list")
    LIST <- LIST[,c("lhs", "op", "rhs", "group", "label")]
    # add est and se column
    est <- object@Fit@est
    BOOT <- attr(est, "BOOT")
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
                free.idx <- which(object@User$free & 
                                  !duplicated(object@User$free))
                ci[free.idx,] <- 2*ci[free.idx,] - t(qq[c(3,4),])

                # def.idx
                def.idx <- which(object@User$op == ":=")
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
                free.idx <- which(object@User$free & 
                                  !duplicated(object@User$free))
                ci[free.idx,] <- t(qq[c(3,4),])

                # def.idx
                def.idx <- which(object@User$op == ":=")
                if(length(def.idx) > 0L) {
                    BOOT.def <- apply(BOOT, 1, object@Model@def.function)
                    if(length(def.idx) == 1L) {
                        BOOT.def <- as.matrix(BOOT.def)
                    } else {
                        BOOT.def <- t(BOOT.def)
                    }
                    qq <- apply(BOOT.def, 2, boot:::norm.inter, alpha)
                    def.idx <- which(object@User$op == ":=")
                    ci[def.idx,] <- t(qq[c(3,4),])
                }

                # TODO:  add cin/ceq?

            } else if(boot.ci.type == "bca.simple") {
               # no adjustment for scale!! only bias!!
               alpha <- (1 + c(-level, level))/2
               zalpha <- qnorm(alpha)
               ci <- cbind(LIST$est, LIST$est)

               # free.idx only
               free.idx <- which(object@User$free & 
                                 !duplicated(object@User$free))
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
               def.idx <- which(object@User$op == ":=")
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

    # if single group, remove group column
    if(object@Sample@ngroups == 1L) LIST$group <- NULL

    # if no user-defined labels, remove label column
    if(sum(nchar(object@User$label)) == 0L) LIST$label <- NULL


    class(LIST) <- c("lavaan.data.frame", "data.frame")
    LIST
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
              what == "sampstat") {
        sampStat(object)
    } else if(what == "rsquare" || 
              what == "r-square" ||
              what == "r2") {
         rsquare(object)
    } else if(what == "converged") {
        object@Fit@converged
    } else {
        stop("unknown `what' argument in inspect function: `", what, "'")
    }

})

# fixme, should we export this function?
sampStat <- function(object, labels=TRUE) {

    G <- object@Sample@ngroups
    ov.names <- object@Sample@ov.names

    OUT <- vector("list", length=G)
    for(g in 1:G) {
        OUT[[g]]$cov  <- object@Sample@cov[[g]]
        if(labels) rownames(OUT[[g]]$cov) <- colnames(OUT[[g]]$cov) <- ov.names
        class(OUT[[g]]$cov) <- c("lavaan.matrix.symmetric", "matrix")

        #if(object@Model@meanstructure) {
            OUT[[g]]$mean <- as.numeric(object@Sample@mean[[g]])
            if(labels) names(OUT[[g]]$mean) <- ov.names
            class(OUT[[g]]$mean) <- c("lavaan.vector", "numeric")
        #}
    }

    if(G == 1) {
        OUT <- OUT[[1]]
    } else {
        names(OUT) <- unlist(object@Sample@group.label)
    }

    OUT
}


setMethod("fitted.values", "lavaan",
function(object, labels=TRUE) {

    G <- object@Sample@ngroups
    ov.names <- object@Sample@ov.names

    OUT <- vector("list", length=G)
    for(g in 1:G) {
        OUT[[g]]$cov  <- object@Fit@Sigma.hat[[g]]
        if(labels) rownames(OUT[[g]]$cov) <- colnames(OUT[[g]]$cov) <- ov.names
        class(OUT[[g]]$cov) <- c("lavaan.matrix.symmetric", "matrix")

        #if(object@Model@meanstructure) {
            OUT[[g]]$mean <- as.numeric(object@Fit@Mu.hat[[g]])
            if(labels) names(OUT[[g]]$mean) <- ov.names
            class(OUT[[g]]$mean) <- c("lavaan.vector", "numeric")
        #}
    }

    if(G == 1) {
        OUT <- OUT[[1]]
    } else {
        names(OUT) <- unlist(object@Sample@group.label)
    }

    OUT
})


setMethod("fitted", "lavaan",
function(object, labels=TRUE) {
     fitted.values(object, labels=labels)
})


setMethod("vcov", "lavaan",
function(object, labels=TRUE) {

    if(object@Fit@npar == 0) {
        VarCov <- matrix(0,0,0)
    } else {
        VarCov <- estimateVCOV(object@Model, sample=object@Sample, 
                               options=object@Options,
                               data=eval(object@call[["data"]], 
                                         parent.frame()) )
    }

    if(labels) {
        colnames(VarCov) <- rownames(VarCov) <- 
            getParameterLabels(object@User, type="free")
    }

    class(VarCov) <- c("lavaan.matrix.symmetric", "matrix")
    VarCov
})


rsquare <- function(object, est.std.all=NULL) {

    if(is.null(est.std.all)) est.std.all <- standardize.est.all(object)
    ngroups <- object@Sample@ngroups
    user <- object@User
    user$rsquare <- 1.0 - est.std.all
    # no values > 1.0
    user$rsquare[user$rsquare > 1.0] <- as.numeric(NA)
    r2 <- vector("list", length=ngroups)

    for(g in 1:ngroups) {
        ind.names   <- user$rhs[ which(user$op == "=~" & user$group == g) ]
        eqs.y.names <- user$lhs[ which(user$op == "~"  & user$group == g) ]
        y.names <- unique( c(ind.names, eqs.y.names) )

        idx <- which(user$op == "~~" & user$lhs %in% y.names & 
                     user$rhs == user$lhs & user$group == g)
        tmp <- user$rsquare[idx]; names(tmp) <- user$lhs[idx]
        r2[[g]] <- tmp
    }

    if(ngroups == 1) {
        r2 <- r2[[1]]
    } else {
        names(r2) <- unlist(object@Sample@group.label)
    }

    r2
}

# logLik (so that we can use the default AIC/BIC functions from stats4(
setMethod("logLik", "lavaan",
function(object, ...) {
    if(object@Options$estimator != "ML") {
        stop("lavaan ERROR: logLik only available if estimator is ML")
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
    object@Sample@ntotal
})

# see: src/library/stats/R/update.R
setMethod("update", signature(object = "lavaan"),
function(object, model.syntax., ..., evaluate = TRUE) {

    call <- object@call
    if(is.null(call))
        stop("need an object with call slot")

    extras <- match.call(expand.dots = FALSE)$...

    if(!missing(model.syntax.))
        #call$formula <- update.formula(formula(object), formula.)
        call$model.syntax <- model.syntax.

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

    mcall <- match.call(expand.dots = TRUE)
    dots <- list(...)
    modp <- if(length(dots))
        sapply(dots, is, "lavaan") else logical(0)

    # single argument version is not supported (what should be display?)
    if(!any(modp)) stop("lavaan ERROR: need at least two models to compare")

    # list of models
    mods <- c(list(object), dots[modp])
    names(mods) <- sapply(as.list(mcall)[c(FALSE, TRUE, modp)], as.character)

    # put them in order (using number of free parameters)
    mods <- mods[order(sapply(lapply(mods, logLik), attr, "df"), 
                 decreasing = TRUE)]

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

    # which models have used a `scaled' test statistic?
    mods.scaled <- unlist( lapply(mods, function(x) {
        any(c("satorra.bentler", "yuan.bentler") %in% 
            unlist(sapply(slot(slot(x, "Fit"), "test"), "[", "test")) ) }))

    if(all(mods.scaled)) {
        scaled <- TRUE
    } else if(!all(mods.scaled)) {
        scaled <- FALSE
    } else {
        warning("lavaan WARNING: some models (but not all) have scaled test statistics")
        scaled <- FALSE
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
        # use formula from mplus web note (www.statmodel.com)
        scaling.factor <- unlist(lapply(mods, function(x) slot(slot(x, "Fit"),
                                     "test")[[2]]$scaling.factor))
        cd <- c(NA, diff(scaling.factor * Df)/diff(Df))
        Chisq.delta <- Chisq.delta/cd

        # print out scaled Chisq for each model
        Chisq <- unlist(lapply(mods, function(x) slot(slot(x, "Fit"),
                        "test")[[2]]$stat))
    }

    # Pvalue
    Pvalue.delta <- pchisq(Chisq.delta, Df.delta, lower = FALSE)

    val <- data.frame(Df = Df,
                      AIC = sapply(mods, AIC),
                      BIC = sapply(mods, BIC),
                      Chisq = Chisq,
                      "Chisq diff" = Chisq.delta,
                      "Df diff" = Df.delta,
                      "Pr(>Chisq)" = Pvalue.delta,
                      row.names = names(mods), 
                      check.names = FALSE)

    if(scaled) {
        attr(val, "heading") <- "Scaled Chi Square Difference Test\n"
    } else {
        attr(val, "heading") <- "Chi Square Difference Test\n"
    }
    class(val) <- c("anova", class(val))

    return(val)

})



