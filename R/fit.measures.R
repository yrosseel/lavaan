fitMeasures <- fitmeasures <- function(object, fit.measures="all") {

    # do we have a test statistic?
    if(object@Fit@test[[1]]$test == "none") {
        stop("lavaan ERROR: please refit the model with test=\"standard\"")
    }

    if("all" %in% fit.measures) {
       class.flag <- TRUE
    } else {
       class.flag <- FALSE
    }

    # reference: Muthen technical Appendix 5

    # collect info from the lavaan slots
    GLIST <- object@Model@GLIST
    N <- object@SampleStats@ntotal
    #q <- length(vnames(object@ParTable, "ov.x"))
    #p <- nvar - q
    npar <- object@Fit@npar
    fx <- object@Fit@fx
    fx.group <- object@Fit@fx.group
    meanstructure <- object@Model@meanstructure
    categorical   <- object@Model@categorical
    multigroup    <- object@Data@ngroups > 1L
    estimator     <- object@Options$estimator
    test          <- object@Options$test
    G <- object@Data@ngroups  # number of groups
    X2 <- object@Fit@test[[1]]$stat
    df <- object@Fit@test[[1]]$df
   
    if(test %in% c("satorra.bentler", "yuan.bentler", "mean.adjusted",
                   "mean.var.adjusted", "scaled.shifted")) {
        scaled <- TRUE
    } else {
        scaled <- FALSE
    }

    # scaled X2
    if(scaled) {
        X2.scaled <- object@Fit@test[[2]]$stat
        df.scaled <- object@Fit@test[[2]]$df
    }


    # select 'default' fit measures
    fit.measures <- tolower(fit.measures)
    if("all" %in% fit.measures) {
        if(scaled) {
            chisq <- c("chisq", "df", "pvalue", "chisq.scaled", "df.scaled",
                       "pvalue.scaled", "chisq.scaling.factor")
        } else {
            chisq <- c("chisq", "df", "pvalue")
        }
        if(scaled) {
            baseline <- c("baseline.chisq", "baseline.df", "baseline.pvalue",
                          "baseline.chisq.scaled", "baseline.df.scaled", 
                          "baseline.pvalue.scaled", 
                          "baseline.chisq.scaling.factor")
        } else {
            baseline <- c("baseline.chisq", "baseline.df", "baseline.pvalue")
        }
        if(scaled) {
            cfi.tli <- c("cfi", "tli", "cfi.scaled", "tli.scaled")
        } else {
            cfi.tli <- c("cfi", "tli")
        }
        if(scaled && object@Options$test == "yuan.bentler") {
            logl <- c("logl", "unrestricted.logl", "npar", "aic", "bic",
                      "scaling.factor.h1", "scaling.factor.h0", "ntotal")
        } else {
            logl <- c("logl", "unrestricted.logl", "npar", "aic", "bic", "ntotal")
        }
        #if(object@Options$mimic == "Mplus") {
            logl <- c(logl, "bic2")
        #}
        if(scaled) {
            rmsea <- c("rmsea", "rmsea.scaled")
            rmsea.ci <- c("rmsea", "rmsea.ci.lower", "rmsea.ci.upper",
                          "rmsea.scaled", "rmsea.ci.lower.scaled", 
                           "rmsea.ci.upper.scaled")
            rmsea.full <- c("rmsea", "rmsea.ci.lower", "rmsea.ci.upper", 
                            "rmsea.pvalue",
                            "rmsea.scaled", "rmsea.ci.lower.scaled", 
                            "rmsea.ci.upper.scaled",
                            "rmsea.pvalue.scaled")
        } else {
            rmsea <- c("rmsea")
            rmsea.ci <- c("rmsea", "rmsea.ci.lower", "rmsea.ci.upper")
            rmsea.full <- c("rmsea", 
                            "rmsea.ci.lower", "rmsea.ci.upper", "rmsea.pvalue")
        }

        if(categorical) {
            srmr <- character(0)
        } else {
            srmr <- c("srmr", "srmr_nomean")
        }
  

        if(estimator == "ML") {
            fit.measures <- c(chisq, baseline, cfi.tli, logl, rmsea.full, srmr)
        } else {
            if(multigroup) {
                fit.measures <- c(chisq, baseline, cfi.tli, 
                                  rmsea.ci, srmr) 
            } else {
                fit.measures <- c(chisq, baseline, cfi.tli, 
                                  rmsea.full, srmr)
            }
        }
    }
    indices <- list()

    # Chi-square value estimated model (H0)
    if(any(c("chisq", "chisq.scaled", 
             "chisq.scaling.factor") %in% fit.measures)) {
	indices["chisq"] <- X2
        if(scaled) {
            indices["chisq.scaled"] <- X2.scaled
            indices["chisq.scaling.factor"] <- object@Fit@test[[2]]$scaling.factor
        } 
    }
    if(any(c("df", "df.scaled") %in% fit.measures)) {
        indices["df"] <- df
        if(scaled) {
            indices["df.scaled"] <- df.scaled
        }
    }
    if(any(c("pvalue", "pvalue.scaled") %in% fit.measures)) {
        indices["pvalue"] <- object@Fit@test[[1]]$pvalue
        if(scaled) {
            indices["pvalue.scaled"] <- object@Fit@test[[2]]$pvalue
        }
    }


    if(any(c("cfi", "cfi.scaled", "tli", "tli.scaled",
             "baseline.chisq", "baseline.chisq.scaled",
             "baseline.pvalue", "baseline.pvalue.scaled") %in% fit.measures)) {

        # call explicitly independence model
        # this is not strictly needed for ML, but it is for
        # GLS and WLS
        # and MLM and MLR to get the scaling factor(s)!
        #if(estimator == "ML") {
        #    if(object@SampleStats@missing.flag) {
        #        do.fit <- TRUE
        #    } else {
        #        do.fit <- FALSE
        #    }
        #} else {
        #    do.fit <- TRUE
        #}

        #OV.X <- character(0L)
        #if(object@Options$mimic == "Mplus") 
        #    OV.X <- vnames(object@ParTable, type="ov.x")

        #indep.syntax <- 
        #    syntax.independence.model(ov.names   = object@Data@ov.names,
        #                              ov.names.x = OV.X,
        #                              sample.cov = object@SampleStats@cov)
        #fit.indep <- update(object, model = indep.syntax, 
        #                    se = "none", do.fit=TRUE, 
        #                    constraints = "",
        #                    verbose = FALSE, warn = FALSE)
        #OCALL <- as.list(object@call); OCALL$env <- NULL; OCALL[[1]] <- NULL
        #NCALL <- list(model = indep.syntax, se = "none", do.fit = TRUE, 
        #              constraints = "", verbose = FALSE, warn = FALSE)
        #CALL  <- modifyList(OCALL, NCALL)
        #cat("DEBUG!\n"); print(as.list(object@call$env)); cat("*******\n")
        #fit.indep <- do.call("lavaan", args=CALL, envir=object@call$env)
        #fit.indep <- do.call("lavaan", args=CALL)

        fit.indep <- independence.model.fit(object)
                            
        X2.null <- fit.indep@Fit@test[[1]]$stat
        df.null <- fit.indep@Fit@test[[1]]$df
        if(scaled) {
            X2.null.scaled <- fit.indep@Fit@test[[2]]$stat
            df.null.scaled <- fit.indep@Fit@test[[2]]$df
        } 

        if("baseline.chisq" %in% fit.measures) {
            indices["baseline.chisq"] <- X2.null
            if(scaled) {
                indices["baseline.chisq.scaled"] <- X2.null.scaled
            } 
        }
        if("baseline.df" %in% fit.measures) {
            indices["baseline.df"] <- df.null
            if(scaled) {
                indices["baseline.df.scaled"] <- df.null.scaled
            }
        }
        if("baseline.pvalue" %in% fit.measures) {
            indices["baseline.pvalue"] <- fit.indep@Fit@test[[1]]$pvalue
            if(scaled) {
                indices["baseline.pvalue.scaled"] <- 
                    fit.indep@Fit@test[[2]]$pvalue
            }
        }
        if("baseline.chisq.scaling.factor" %in% fit.measures) {
            indices["baseline.chisq.scaling.factor"] <-
                fit.indep@Fit@test[[2]]$scaling.factor
        }

        # CFI 
        if("cfi" %in% fit.measures) {
            indices["cfi"] <- ( 1 - max(c(X2 - df,0)) / 
                                    max( c(X2-df, X2.null-df.null, 0) ) 
                              )
        }
        if("cfi.scaled" %in% fit.measures) {
            indices["cfi.scaled"] <- 
                ( 1 - max( c(X2.scaled - df.scaled,0)) / 
                      max( c(X2.scaled - df.scaled, 
                             X2.null.scaled - df.null.scaled, 0) ) 
                )
        }

        # TLI 
        if("tli" %in% fit.measures) {
            if(df > 0) {
                TLI <- (X2.null/df.null - X2/df)/(X2.null/df.null - 1)
            } else {
                TLI <- 1
            }
            indices["tli"] <- TLI
        }
        if("tli.scaled" %in% fit.measures) {
            if(df > 0) {
                TLI <- (X2.null.scaled/df.null.scaled - 
                        X2.scaled/df.scaled) / 
                       (X2.null.scaled/df.null.scaled - 1)
            } else {
                TLI <- 1
            }
            indices["tli.scaled"] <- TLI
        }
    }

    if(estimator == "ML") {
        if("logl" %in% fit.measures ||
           "unrestricted.logl" %in% fit.measures ||
           "npar" %in% fit.measures ||
           "aic" %in% fit.measures ||
           "bic" %in% fit.measures) {

            # logl H1 -- unrestricted (aka saturated) model
            logl.H1.group <- numeric(G)
            for(g in 1:G) {
                nvar <- ncol(object@SampleStats@cov[[g]])
                if(!object@SampleStats@missing.flag) {
                    Ng <- object@SampleStats@nobs[[g]]
                    c <- Ng*nvar/2 * log(2 * pi)
                    logl.H1.group[g] <- ( -c -(Ng/2) *
                                          object@SampleStats@cov.log.det[[g]]
                                          - (Ng/2)*nvar )
                } else { # missing patterns case
                    pat <- object@Data@Mp[[g]]$pat
                    Ng <- object@Data@nobs[[g]]
                    ni <- as.numeric(apply(pat, 1, sum) %*% 
                                     as.integer(rownames(pat)))
                    fx.full <- object@SampleStats@missing.h1[[g]]$h1
                    logl.H1.group[g] <- - (ni/2 * log(2 * pi)) - 
                                              (Ng/2 * fx.full)
                }
            }
            if(G > 1) {
                logl.H1 <- sum(logl.H1.group)
            } else {
                logl.H1 <- logl.H1.group[1]
            }

            if("unrestricted.logl" %in% fit.measures) {
                indices["unrestricted.logl"] <- logl.H1
            }

            # logl H0
            logl.H0.group <- numeric(G)
            for(g in 1:G) {
                Ng <- object@SampleStats@nobs[[g]]
                logl.H0.group[g] <- -Ng * (fx.group[g] - logl.H1.group[g]/Ng)
            }
            if(G > 1) {
                logl.H0 <- sum(logl.H0.group)
            } else {
                logl.H0 <- logl.H0.group[1]
            }
           
            if("logl" %in% fit.measures) {
                indices["logl"] <- logl.H0
            }

            # Number of free parameters
            if("npar" %in% fit.measures) {
                indices["npar"] <- npar
            }

            # AIC
            AIC <-  -2*logl.H0 + 2*npar
            if("aic" %in% fit.measures) {
                indices["aic"] <- AIC
            }

            # BIC
            if("bic" %in% fit.measures) {
                BIC <- -2*logl.H0 + npar*log(N)
                indices["bic"] <- BIC

                # add sample-size adjusted bic
                N.star <- (N + 2) / 24
                BIC2 <- -2*logl.H0 + npar*log(N.star)
                indices["bic2"] <- BIC2
            }
           

            # scaling factor for MLR
            if(object@Options$test == "yuan.bentler") {
                indices["scaling.factor.h1"] <- 
                    object@Fit@test[[2]]$scaling.factor.h1
                indices["scaling.factor.h0"] <- 
                    object@Fit@test[[2]]$scaling.factor.h0
            }
        }
    }

    N.RMSEA <- max(N, X2*2) # FIXME: good strategy??
    if(any(c("rmsea","rmsea.scaled") %in% fit.measures)) {
        # RMSEA
        if(df > 0) {
            if(scaled) {
                d <- sum(object@Fit@test[[2]]$trace.UGamma)
            } 
            if(object@Options$mimic %in% c("Mplus", "lavaan")) {
                GG <- 0
                RMSEA <- sqrt( max( c((X2/N)/df - 1/(N-GG), 0) ) ) * sqrt(G)
                if(scaled) {
                    RMSEA.scaled <- 
                         sqrt( max( c((X2/N)/d - 1/(N-GG), 0) ) ) * sqrt(G)
                }
            } else {
                RMSEA <- sqrt( max( c((X2/N)/df - 1/N, 0) ) )
                if(scaled) {
                    RMSEA.scaled <- sqrt( max( c((X2/N)/d - 1/N, 0) ) )
                }
            }
        } else {
            RMSEA <- RMSEA.scaled <- 0
        }
        indices["rmsea"] <- RMSEA
        if(scaled) {
            indices["rmsea.scaled"] <- RMSEA.scaled
        }
    }

    if("rmsea.ci.lower" %in% fit.measures) {
        lower.lambda <- function(lambda) {
            (pchisq(X2, df=df, ncp=lambda) - 0.95)
        }
        if(df < 1 || lower.lambda(0) < 0.0) {
            indices["rmsea.ci.lower"] <- 0
        } else {
            lambda.l <- try(uniroot(f=lower.lambda, lower=0, upper=X2)$root)
            if(inherits(lambda.l, "try-error")) { lambda.l <- NA }
            if(object@Options$mimic %in% c("lavaan", "Mplus")) {
                GG <- 0
                indices["rmsea.ci.lower"] <- 
                    sqrt( lambda.l/((N-GG)*df) ) * sqrt(G)
            } else {
                indices["rmsea.ci.lower"] <- sqrt( lambda.l/(N*df) )
            }
        }
    }

    if("rmsea.ci.lower.scaled" %in% fit.measures) {
        df2 <- sum(object@Fit@test[[2]]$trace.UGamma)
        lower.lambda <- function(lambda) {
            (pchisq(X2, df=df2, ncp=lambda) - 0.95)
        }
        if(df < 1 || df2 < 1 || lower.lambda(0) < 0.0) {
            indices["rmsea.ci.lower.scaled"] <- 0
        } else {
            lambda.l <- try(uniroot(f=lower.lambda, lower=0, upper=X2)$root)
            if(inherits(lambda.l, "try-error")) { lambda.l <- NA }
            if(object@Options$mimic %in% c("lavaan", "Mplus")) {
                GG <- 0
                indices["rmsea.ci.lower.scaled"] <- 
                    sqrt( lambda.l/((N-GG)*df2) ) * sqrt(G)
            } else {
                indices["rmsea.ci.lower.scaled"] <- sqrt( lambda.l/(N*df2) )
            }
        }
    }

    if("rmsea.ci.upper" %in% fit.measures) {
        upper.lambda <- function(lambda) {
            (pchisq(X2, df=df, ncp=lambda) - 0.05)
        }
        if(df < 1 || upper.lambda(N.RMSEA) > 0 || upper.lambda(0) < 0) {
            indices["rmsea.ci.upper"] <- 0
        } else {
            lambda.u <- try(uniroot(f=upper.lambda, lower=0, upper=N.RMSEA)$root)
            if(inherits(lambda.u, "try-error")) { lambda.u <- NA }
            if(object@Options$mimic %in% c("lavaan", "Mplus")) {
                GG <- 0
                indices["rmsea.ci.upper"] <- 
                    sqrt( lambda.u/((N-GG)*df) ) * sqrt(G)
            } else {
                indices["rmsea.ci.upper"] <- sqrt( lambda.u/(N*df) )
            }
        }
    }

    if("rmsea.ci.upper.scaled" %in% fit.measures) {
        df2 <- sum(object@Fit@test[[2]]$trace.UGamma)
        upper.lambda <- function(lambda) {
            (pchisq(X2, df=df2, ncp=lambda) - 0.05)
        }
        if(df < 1 || df2 < 1 || upper.lambda(N.RMSEA) > 0) {
            indices["rmsea.ci.upper.scaled"] <- 0
        } else {
            lambda.u <- try(uniroot(f=upper.lambda, lower=0, upper=N.RMSEA)$root)
            if(inherits(lambda.u, "try-error")) { lambda.u <- NA }
            if(object@Options$mimic %in% c("lavaan", "Mplus")) {
                GG <- 0
                indices["rmsea.ci.upper.scaled"] <- 
                    sqrt( lambda.u/((N-GG)*df2) ) * sqrt(G)
            } else {
                indices["rmsea.ci.upper.scaled"] <- 
                    sqrt( lambda.u/(N*df2) )
            }
        }
    }    

    if("rmsea.pvalue" %in% fit.measures) {
        if(df > 0) {
            if(object@Options$mimic %in% c("lavaan","Mplus")) {
                GG <- 0
                ncp <- (N-GG)*df*0.05^2/G
                indices["rmsea.pvalue"] <- 
                    1 - pchisq(X2, df=df, ncp=ncp)
            } else {
                indices["rmsea.pvalue"] <- 
                    1 - pchisq(X2, df=df, ncp=(N*df*0.05^2))
            }
        } else {
            indices["rmsea.pvalue"] <- 1
        }
    }

    if("rmsea.pvalue.scaled" %in% fit.measures) {
        df2 <- sum(object@Fit@test[[2]]$trace.UGamma)
        if(df > 0) {
            if(object@Options$mimic %in% c("lavaan", "Mplus")) {
                GG <- 0
                ncp <- (N-GG)*df2*0.05^2/G
                indices["rmsea.pvalue.scaled"] <- 
                    1 - pchisq(X2, df=df2, ncp=ncp)
            } else {
                indices["rmsea.pvalue.scaled"] <-
                    1 - pchisq(X2, df=df2, ncp=(N*df2*0.05^2))
            }
        } else {
            indices["rmsea.pvalue.scaled"] <- 1
        }
    }

    if("srmr" %in% fit.measures) {
        # SRMR
        srmr.group <- numeric(G)
        srmr_nomean.group <- numeric(G)
        for(g in 1:G) {
            # observed
            if(!object@SampleStats@missing.flag) {
                S <- object@SampleStats@cov[[g]]
                M <- object@SampleStats@mean[[g]]
            } else {
                # EM estimates
                S <- object@SampleStats@missing.h1[[g]]$sigma
                M <- object@SampleStats@missing.h1[[g]]$mu
            }
            nvar <- ncol(S)

            # estimated
            Sigma.hat <- object@Fit@Sigma.hat[[g]]
            Mu.hat    <- object@Fit@Mu.hat[[g]]

            # standardized residual covariance matrix
            # this is the Hu and Bentler definition, not the Bollen one!
            sqrt.d <- 1/sqrt(diag(S))
            D <- diag(sqrt.d, ncol=length(sqrt.d))
            R <- D %*% (S - Sigma.hat) %*% D

            # this is what the Mplus documentation suggest, 
            # but is not what is used!
            #sqrt.d2 <- 1/sqrt(diag(Sigma.hat))
            #D2 <- diag(sqrt.d2, ncol=length(sqrt.d2))
            #R <- D %*% S %*% D   - D2 %*% Sigma.hat %*% D2

            if(meanstructure) {
                # standardized residual mean vector
                R.mean <- D %*% (M - Mu.hat)
                #R.mean <-  D %*% M - D2 %*% Mu.hat
                e <- nvar*(nvar+1)/2 + nvar
                srmr.group[g] <- sqrt( (sum(R[lower.tri(R, diag=TRUE)]^2) +
                                        sum(R.mean^2))/ e )
                e <- nvar*(nvar+1)/2
                srmr_nomean.group[g] <- sqrt( sum(R[lower.tri(R, diag=TRUE)]^2) / e )
                srmr_nomean.group[g] 
            } else {
                e <- nvar*(nvar+1)/2
                srmr_nomean.group[g] <- srmr.group[g] <- sqrt( sum(R[lower.tri(R, diag=TRUE)]^2) / e )
                
            }
        }
        
        if(G > 1) {
            ## FIXME: get the scaling right
            SRMR <- as.numeric( (unlist(object@SampleStats@nobs) %*% srmr.group) / object@SampleStats@ntotal )
            SRMR_NOMEAN <- as.numeric( (unlist(object@SampleStats@nobs) %*% srmr_nomean.group) / object@SampleStats@ntotal )
        } else {
            SRMR <- srmr.group[1]
            SRMR_NOMEAN <- srmr_nomean.group[1]
        }

        indices["srmr"] <- SRMR
        indices["srmr_nomean"] <- SRMR_NOMEAN
    }

    if("ntotal" %in% fit.measures) {
        indices["ntotal"] <- object@SampleStats@ntotal
    }

    # do we have everything that we requested?
    idx.missing <- which(is.na(match(fit.measures, names(indices))))
    if(length(idx.missing) > 0L) {
        cat("lavaan WARNING: some requested fit measure(s) are not available for this model:\n")
        print( fit.measures[ idx.missing ] )
        cat("\n")
    }
    
    out <- unlist(indices[fit.measures])

    if(length(out) > 0L) {
        class(out) <- c("lavaan.vector", "numeric")
    } else {
        return( invisible(numeric(0)) )
    }

    out
}

# print a nice summary of the fit measures
print.fit.measures <- function(x) {

   # scaled?
   scaled <- "chisq.scaled" %in% names(x)

   # independence model
   if("baseline.chisq" %in% names(x)) {
       cat("Chi-square test baseline model:\n\n")
       t0.txt <- sprintf("  %-40s", "Minimum Function Chi-square")
       t1.txt <- sprintf("  %10.3f", x["baseline.chisq"])
       t2.txt <- ifelse(scaled,
                 sprintf("  %10.3f", x["baseline.chisq.scaled"]), "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

       t0.txt <- sprintf("  %-40s", "Degrees of freedom")
       t1.txt <- sprintf("  %10i", x["baseline.df"])
       t2.txt <- ifelse(scaled,
                 sprintf("  %10i", x["baseline.df.scaled"]), "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

       t0.txt <- sprintf("  %-40s", "P-value")
       t1.txt <- sprintf("  %10.3f", x["baseline.pvalue"])
       t2.txt <- ifelse(scaled,
                 sprintf("  %10.3f", x["baseline.pvalue.scaled"]), "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
    }

   # cfi/tli
   if("cfi" %in% names(x)) {
       cat("\nFull model versus baseline model:\n\n")
       t0.txt <- sprintf("  %-40s", "Comparative Fit Index (CFI)")
       t1.txt <- sprintf("  %10.3f", x["cfi"])
       t2.txt <- ifelse(scaled,
                 sprintf("  %10.3f", x["cfi.scaled"]), "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

       t0.txt <- sprintf("  %-40s", "Tucker-Lewis Index (TLI)")
       t1.txt <- sprintf("  %10.3f", x["tli"])
       t2.txt <- ifelse(scaled,
                 sprintf("  %10.3f", x["tli.scaled"]), "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
   }

   # likelihood
   if("logl" %in% names(x)) {
       cat("\nLoglikelihood and Information Criteria:\n\n")
       t0.txt <- sprintf("  %-40s", "Loglikelihood user model (H0)")
       t1.txt <- sprintf("  %10.3f", x["logl"])
       t2.txt <- ifelse(scaled,
                 sprintf("  %10.3f", x["logl"]), "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       #cat(t0.txt, t1.txt, "\n", sep="")
       if(!is.na(x["scaling.factor.h0"])) {
           t0.txt <- sprintf("  %-40s", "Scaling correction factor")
           t1.txt <- sprintf("  %10s", "")
           t2.txt <- sprintf("  %10.3f", x["scaling.factor.h0"])
           cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
           cat("    for the MLR correction\n")
       }

       t0.txt <- sprintf("  %-40s", "Loglikelihood unrestricted model (H1)")
       t1.txt <- sprintf("  %10.3f", x["unrestricted.logl"])
       t2.txt <- ifelse(scaled,
                 sprintf("  %10.3f", x["unrestricted.logl"]), "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       #cat(t0.txt, t1.txt, "\n", sep="")
       if(!is.na(x["scaling.factor.h1"])) {
           t0.txt <- sprintf("  %-40s", "Scaling correction factor")
           t1.txt <- sprintf("  %10s", "")
           t2.txt <- sprintf("  %10.3f", x["scaling.factor.h1"])
           cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
           cat("    for the MLR correction\n")
       }

       cat("\n")
       t0.txt <- sprintf("  %-40s", "Number of free parameters")
       t1.txt <- sprintf("  %10i", x["npar"])
       t2.txt <- ifelse(scaled,
                 sprintf("  %10i", x["npar"]), "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")

       t0.txt <- sprintf("  %-40s", "Akaike (AIC)")
       t1.txt <- sprintf("  %10.3f", x["aic"])
       t2.txt <- ifelse(scaled,
                 sprintf("  %10.3f", x["aic"]), "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       #cat(t0.txt, t1.txt, "\n", sep="")
       t0.txt <- sprintf("  %-40s", "Bayesian (BIC)")
       t1.txt <- sprintf("  %10.3f", x["bic"])
       t2.txt <- ifelse(scaled,
                 sprintf("  %10.3f", x["bic"]), "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       #cat(t0.txt, t1.txt, "\n", sep="")
       if(!is.na(x["bic2"])) {
           t0.txt <- sprintf("  %-40s", "Sample-size adjusted Bayesian (BIC)")
           t1.txt <- sprintf("  %10.3f", x["bic2"])
           t2.txt <- ifelse(scaled,
                     sprintf("  %10.3f", x["bic2"]), "")
           cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
        }
   }

   # RMSEA
   if("rmsea" %in% names(x)) {
       cat("\nRoot Mean Square Error of Approximation:\n\n")
       t0.txt <- sprintf("  %-40s", "RMSEA")
       t1.txt <- sprintf("  %10.3f", x["rmsea"])
       t2.txt <- ifelse(scaled,
                 sprintf("  %10.3f", x["rmsea.scaled"]), "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       if("rmsea.ci.lower" %in% names(x)) {
           t0.txt <- sprintf("  %-38s", "90 Percent Confidence Interval")
           t1.txt <- sprintf("  %5.3f", x["rmsea.ci.lower"])
           t2.txt <- sprintf("  %5.3f", x["rmsea.ci.upper"])
           t3.txt <- ifelse(scaled,
                     sprintf("       %5.3f  %5.3f", x["rmsea.ci.lower.scaled"],
                                               x["rmsea.ci.upper.scaled"]), "")
           cat(t0.txt, t1.txt, t2.txt, t3.txt, "\n", sep="")
       }
       if("rmsea.pvalue" %in% names(x)) {
           t0.txt <- sprintf("  %-40s", "P-value RMSEA <= 0.05")
           t1.txt <- sprintf("  %10.3f", x["rmsea.pvalue"])
           t2.txt <- ifelse(scaled,
                 sprintf("  %10.3f", x["rmsea.pvalue.scaled"]), "")
           cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       }
   }

   # SRMR
   if("srmr" %in% names(x)) {
       cat("\nStandardized Root Mean Square Residual:\n\n")
       t0.txt <- sprintf("  %-40s", "SRMR")
       t1.txt <- sprintf("  %10.3f", x["srmr"])
       t2.txt <- ifelse(scaled,
                 sprintf("  %10.3f", x["srmr"]), "")
       cat(t0.txt, t1.txt, t2.txt, "\n", sep="")
       #cat(t0.txt, t1.txt, "\n", sep="")
   }

   cat("\n")
}


