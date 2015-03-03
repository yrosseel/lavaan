ctr_pml_plrt <- function(lavobject = NULL, lavmodel = NULL, lavdata = NULL,
                         lavsamplestats = NULL, lavpartable = NULL,
                         lavoptions = NULL, x = NULL, VCOV = NULL,
                         lavcache = NULL) {

    if(!is.null(lavobject)) {
        lavmodel <- lavobject@Model
        lavdata <- lavobject@Data
        lavoptions <- lavobject@Options
        lavsamplestats <- lavobject@SampleStats
        lavcache <- lavobject@Cache
        lavpartable <- lavobject@ParTable
    }


    if(is.null(x)) {
        # compute 'fx' = objective function value 
        # (NOTE: since 0.5-18, NOT divided by N!!)
        fx <- lav_model_objective(lavmodel       = lavmodel,
                                  lavsamplestats = lavsamplestats,
                                  lavdata        = lavdata,
                                  lavcache       = lavcache,
                                  estimator      = "PML")
        H0.fx <- as.numeric(fx)
        H0.fx.group <- attr(fx, "fx.group")
    } else {
        H0.fx <- attr(attr(x, "fx"), "fx.pml")
        H0.fx.group <- attr(attr(x, "fx"), "fx.group")
    }

    # fit a saturated model 'fittedSat'
    ModelSat <- lav_partable_unrestricted(lavobject      = NULL,
                                          lavdata        = lavdata,
                                          lavoptions     = lavoptions,
                                          lavsamplestats = lavsamplestats)

    # FIXME: se="none", test="none"??
    Options <- lavoptions
    Options$verbose <- FALSE
    Options$se <- "none"
    Options$test <- "none"
    fittedSat <- lavaan(ModelSat, slotOptions = Options,
                        slotSampleStats = lavsamplestats,
                        slotData = lavdata, slotCache = lavcache)
    fx <- lav_model_objective(lavmodel = fittedSat@Model,
                              lavsamplestats = fittedSat@SampleStats,
                              lavdata = fittedSat@Data,
                              lavcache = fittedSat@Cache,
                              estimator = "PML")
    SAT.fx <- as.numeric(fx)
    SAT.fx.group <- attr(fx, "fx.group")

    # we also need a `saturated model', but where the moments are based
    # on the model-implied sample statistics under H0
    ModelSat2 <- 
        lav_partable_unrestricted(lavobject      = NULL,
                                  lavdata        = lavdata,
                                  lavoptions     = lavoptions,
                                  lavsamplestats = NULL,
                                  sample.cov     = computeSigmaHat(lavmodel),
                                  sample.mean    = computeMuHat(lavmodel), 
                                  sample.th      = computeTH(lavmodel))

    #Options$se <- lavoptions$se
    fittedSat2 <- lavaan(ModelSat2, 
                        control=list(optim.method          = "none",
                                     optim.force.converged = TRUE) ,
                        #slotOptions = lavoptions,
                        slotOptions = Options,
                        slotSampleStats = lavsamplestats,
                        slotData = lavdata, slotCache = lavcache)

    # the code below was contributed by Myrsini Katsikatsou (Jan 2015)

# for now, only a single group is supported:    
g = 1L


########################### The code for PLRT for overall goodness of fit

# First define the number of non-redundant elements of the (fitted)
# covariance/correlation matrix of the underlying variables.
#nvar <- lavmodel@nvar[[g]] 
#dSat <- nvar*(nvar-1)/2
#if(length(lavmodel@num.idx[[g]]) > 0L) {
#    dSat <- dSat + length(lavmodel@num.idx[[g]])
#}

# select `free' parameters (excluding thresholds) from fittedSat2 model
PT.Sat2 <- fittedSat2@ParTable
dSat.idx <- PT.Sat2$free[ PT.Sat2$free > 0L & PT.Sat2$op != "|" ] # remove thresholds

# Secondly, we need to specify the indices of the rows/columns of vcov(), hessian, and
# variability matrix that refer to all SEM parameters except thresholds.
PT <- lavpartable
index.par <- PT$free[PT$free > 0L & PT$op != "|"]

# Thirdly, specify the sample size.
# nsize <- lavdata@nobs[[g]]
nsize <- lavsamplestats@ntotal

# Now we can proceed to the computation of the quantities needed for PLRT.
# Briefly, to say that PLRT is equal to the difference of two quadratic forms.
# To compute the first and second moment adjusted PLRT we should compute
# the asymptotic mean and variance of each quadratic quantity as well as
# their asymptotic covariance.

##### Section 1. Compute the asymptotic mean and variance of the first quadratic quantity
# Below I assume that lavobject is the output of lavaan function. I guess
# vcov(lavobject) can be substituted by VCOV object insed lavaan function
# defined at lines 703 -708. But what is the object inside lavaan function
# for getHessian(lavobject)?
if(is.null(VCOV)) {
    VCOV <- lav_model_vcov(lavmodel       = lavmodel,
                           lavsamplestats = lavsamplestats,
                           lavoptions     = lavoptions,
                           lavdata        = lavdata,
                           lavpartable    = lavpartable,
                           lavcache       = lavcache)
}
InvG_to_psipsi_attheta0 <- (1 * VCOV )[index.par, index.par]  #G^psipsi(theta0)
#below the lavaan function getHessian is used
#Hattheta0 <- (-1) * H0.Hessian
#Hattheta0 <- H0.Hessian
#InvHattheta0 <- solve(Hattheta0)
InvHattheta0 <- attr(VCOV, "E.inv")
InvH_to_psipsi_attheta0 <- InvHattheta0[index.par, index.par]   #H^psipsi(theta0)
if(lavmodel@eq.constraints) {
    IN <- InvH_to_psipsi_attheta0
    IN.npar <- ncol(IN)

    # create `bordered' matrix
    if(nrow(lavmodel@con.jac) > 0L) {
        H <- lavmodel@con.jac[, index.par]
        inactive.idx <- attr(H, "inactive.idx")
        lambda <- lavmodel@con.lambda # lagrangean coefs
        if(length(inactive.idx) > 0L) {
            H <- H[-inactive.idx,,drop=FALSE]
            lambda <- lambda[-inactive.idx]
        }
        if(nrow(H) > 0L) {
            H0 <- matrix(0,nrow(H),nrow(H))
            H10 <- matrix(0, ncol(IN), nrow(H))
            DL <- 2*diag(lambda, nrow(H), nrow(H))
            # FIXME: better include inactive + slacks??
            E3 <- rbind( cbind(     IN,  H10, t(H)),
                         cbind( t(H10),   DL,  H0),
                         cbind(      H,   H0,  H0)  )
            Inv_of_InvH_to_psipsi_attheta0 <-
                MASS::ginv(IN)[1:IN.npar, 1:IN.npar, drop = FALSE]
        } else {
            Inv_of_InvH_to_psipsi_attheta0 <- solve(IN)
        }
    }
} else {
    Inv_of_InvH_to_psipsi_attheta0 <- 
        solve(InvH_to_psipsi_attheta0) #[H^psipsi(theta0)]^(-1)
}

H0tmp_prod1 <- Inv_of_InvH_to_psipsi_attheta0 %*% InvG_to_psipsi_attheta0
H0tmp_prod2 <- H0tmp_prod1 %*% H0tmp_prod1
E_tww <- sum(diag(H0tmp_prod1)) #expected mean of the first quadratic quantity
var_tww <- 2* sum(diag(H0tmp_prod2)) #variance of the first quadratic quantity

##### Section 2: Compute the asymptotic mean and variance of the second quadratic quantity.
# Now we need to evaluate the fitted (polychoric) correlation/ covariance matrix
# using the estimates of SEM parameters derived under the fitted model
# which is the model of the null hypothesis. We also need to compute the
# vcov matrix of these estimates (estimates of polychoric correlations)
# as well as the related hessian and variability matrix. 
tmp.options <- fittedSat2@Options
tmp.options$se <- lavoptions$se
VCOV.Sat2 <- lav_model_vcov(lavmodel       = fittedSat2@Model,
                            lavsamplestats = fittedSat2@SampleStats,
                            lavoptions     = tmp.options,
                            lavdata        = fittedSat2@Data,
                            lavpartable    = fittedSat2@ParTable,
                            lavcache       = fittedSat2@Cache)
InvG_to_sigmasigma_attheta0 <- VCOV.Sat2[dSat.idx, dSat.idx, drop = FALSE]  #G^sigmasigma(theta0)
#Hattheta0 <- (-1)* getHessian(fittedSat2)
#Hattheta0 <- getHessian(fittedSat2)
#InvHattheta0 <- solve(Hattheta0)
InvHattheta0 <- attr(VCOV.Sat2, "E.inv")
InvH_to_sigmasigma_attheta0 <- InvHattheta0[dSat.idx, dSat.idx, drop = FALSE] #H^sigmasigma(theta0)
Inv_of_InvH_to_sigmasigma_attheta0 <- solve(InvH_to_sigmasigma_attheta0) #[H^sigmasigma(theta0)]^(-1)
H1tmp_prod1 <- Inv_of_InvH_to_sigmasigma_attheta0 %*% InvG_to_sigmasigma_attheta0
H1tmp_prod2 <- H1tmp_prod1 %*% H1tmp_prod1
E_tzz <- sum(diag(H1tmp_prod1))     #expected mean of the second quadratic quantity
var_tzz <- 2* sum(diag(H1tmp_prod2))#variance of the second quadratic quantity



##### Section 3: Compute the asymptotic covariance of the two quadratic quantities
deltamat <- computeDelta(lavmodel)[[g]] # [[1]] to be substituted by g?
# The above gives the derivatives of thresholds and polychoric correlations
# with respect to SEM param (including thresholds) evaluated under H0.
# From deltamat we need to exclude the rows and columns referring to thresholds.
# For this:
# free_TH_indices <- lavmodel@x.free.idx[[6]]
free_TH_indices <- PT$free[PT$free > 0L & PT$op == "|"]
noTH <- length(free_TH_indices)
#nrowsDelta <- noTH + dSat
nrowsDelta <- noTH + length(dSat.idx)
drhodpsi_mat <- deltamat[(noTH+1):nrowsDelta , index.par ]
tmp_prod <- t(drhodpsi_mat) %*% Inv_of_InvH_to_sigmasigma_attheta0 %*%
            drhodpsi_mat %*%
            InvG_to_psipsi_attheta0 %*% H0tmp_prod1
cov_tzztww <- 2* sum(diag(tmp_prod))

##### Section 4: compute the adjusted PLRT and its p-value
# PLRTH0Sat <- 2*nsize*(lavfit@fx - fittedSat@Fit@fx)
PLRTH0Sat <- 2*(H0.fx - SAT.fx)
asym_mean_PLRTH0Sat <- E_tzz - E_tww
asym_var_PLRTH0Sat  <- var_tzz + var_tww -2*cov_tzztww
scaling.factor <- (asym_mean_PLRTH0Sat / (asym_var_PLRTH0Sat/2) )
FSA_PLRT_SEM <- (asym_mean_PLRTH0Sat / (asym_var_PLRTH0Sat/2) )* PLRTH0Sat
adjusted_df  <- (asym_mean_PLRTH0Sat^2) / (asym_var_PLRTH0Sat/2)
# In some very few cases (simulations show very few cases in small sample sizes)
# the adjusted_df is a negative number, we should then
# print a warning like: "The adjusted df is computed to be a negative number
# and for this the first and second moment adjusted PLRT is not computed." .
pvalue <- 1-pchisq(FSA_PLRT_SEM, df=adjusted_df )

list(PLRTH0Sat = PLRTH0Sat, stat = FSA_PLRT_SEM, df = adjusted_df, p.value = pvalue, scaling.factor = scaling.factor)
}
############################################################################


ctr_pml_aic_bic <- function(lavobject) {

    lavmodel <- lavobject@Model
    lavdata <- lavobject@Data
    lavoptions <- lavobject@Options
    lavsamplestats <- lavobject@SampleStats
    lavcache <- lavobject@Cache
    lavpartable <- lavobject@ParTable

    nsize <- lavobject@SampleStats@ntotal

########################## The code for PL version fo AIC and BIC
# The following should be done because it is not the pl log-likelihood
# that is maximized but a fit function that should be minimized. So, we
# should find the value of log-PL at the estimated parameters through the
# value of the fitted function.
# The following may need to be updated if we change the fit function
# so that it is correct for the case of missing values as well.
# Below lavcache as defined in the code of lavaan function and
# [[1]] may need to be substituted by g.
#nsize <- lavdata@nobs[[1]]
#ObsProp <- lavcache[[1]]$bifreq / nsize
#zero.idx <- which(prop == 0.0)
#prop <- prop[-zero.idx]
#logPL <-  (-1)* nsize * ( lavfit@fx - sum( prop*log(prop) )  )
#fmin <- lav_model_objective(lavmodel = lavobject@Model, lavsamplestats = lavobject@SampleStats, lavdata = lavobject@Data, lavcache = lavobject@Cache, estimator = "PML")
#logPL <- sum(attr(fmin, "logl.group"))
    logPL <- lavobject@Fit@logl

# Find the right dimension of the (theta) parameter.
# lavobject is the output of function lavaan when we fit the model of interest.
# I don't know how else to get the Hessian through the code of lavaan function.
# That's why I'm using below the function getHessian(lavobject).
#Hes <- (-1)* getHessian(lavobject)
#Hes <- getHessian(lavobject)
#InvG <- nsize * vcov(lavobject)
#InvG <- 1 * vcov(lavobject)
# I guess vcov(lavobject) can be substituted by object VCOV computed inside
# the lavaan function.
#dimTheta <- sum(diag(Hes %*% InvG))

    VCOV <- lav_model_vcov(lavmodel       = lavmodel,
                           lavsamplestats = lavsamplestats,
                           lavoptions     = lavoptions,
                           lavdata        = lavdata,
                           lavpartable    = lavpartable,
                           lavcache       = lavcache)

    H.inv <- attr(VCOV, "E.inv")
    B0.group <- attr(VCOV, "B0.group") # B0 = J(theta) per group

    if(lavsamplestats@ngroups > 1L) {
        # groups weights
        B0 <- (lavsamplestats@nobs[[1]]/lavsamplestats@ntotal) * B0.group[[1]]
        for(g in 2:lavsamplestats@ngroups) {
            B0 <- B0 + (lavsamplestats@nobs[[g]]/lavsamplestats@ntotal) * B0.group[[g]]
        }
    } else {
        B0 <- B0.group[[1]]
    }

    dimTheta <- sum(B0 * H.inv)


# computations of PL versions of AIC and BIC
PL_AIC <- (-2)*logPL + 2*dimTheta
PL_BIC <- (-2)*logPL + dimTheta *log(nsize)

list(logPL = logPL, PL_AIC = PL_AIC, PL_BIC = PL_BIC)
}
