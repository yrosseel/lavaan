# How does lavaan work?

library(lavaan)

# PART 1: from syntax to matrices

model <- 'f =~ x1 + a*x2 + 3*x3'

# parsing the syntax
lavParseModelString(model, as.data.frame. = TRUE)
# creates a 'FLAT' initial parameter table
FLAT <- lavParseModelString(model)
lavNames(FLAT)
lavNames(FLAT, "lv")

# lavaanify()
# - first creates FLAT
# - needs information from the data (eg number of groups)
# - builds the initial parameter table
# - the parameter table is the core representation of the model in lavaan
lavaanify(model)

# if sem/cfa is used, more parameters are set free
lavaanify(model, auto = TRUE)

# example: equality constraints using labels
model <- 'f =~ x1 + a*x2 + a*x3'
lavaanify(model, auto = TRUE)

# alternative for 'simple' equality constraints
# (will become the default soon)
lavaanify(model, auto = TRUE, ceq.simple = TRUE)

# explicit equality constraints
model <- 'f =~ x1 + a*x2 + b*x3; a == b'
lavaanify(model, auto = TRUE, ceq.simple = TRUE)

# multiple groups/blocks
model <- 'f =~ x1 + c(a1,a2)*x2 + c(b1, b2)*x3'
lavaanify(model, auto = TRUE, ngroups = 2)

# matrix representation: LISREL-all-y
model <- 'f =~ x1 + x2 + x3'
PT <- lavaanify(model, auto = TRUE, as.data.frame. = TRUE)
PT
# map every parameter to a matrix element
MAT <- as.data.frame(lavaan:::lav_lisrel(PT))
cbind(PT, MAT)

# alternative matrix representation: RAM
MAT <- as.data.frame(lavaan:::lav_ram(PT))
cbind(PT, MAT)

# first create default lavoptions
fit <- sem(model = 'f =~x1 + x2 + x3', data = HolzingerSwineford1939)
lavoptions <- fit@Options

# create 'Model' object (S4 class)
# (often called lavmodel in internal functions)
Model <- lavaan:::lav_model(PT, lavoptions = lavoptions)
# another representation of the model, more suitable for computations
Model@nmat
Model@isSymmetric
Model@dimNames

# the matrix representation is in GLIST (group list)
# GLIST may contain the matrices of multiple groups, levels, blocks, ...
MLIST <- Model@GLIST

# MLIST is the matrix representation of a single block
# MLIST is used to compute model-based statistics, e.g.,
# the model-implied variance-covariance matrix
lavaan:::lav_lisrel_sigma(MLIST)




# PART 2: lavaan() workflow:
# 18 steps in total (many of them can be skipped)

# 1. lavParseModelString + ov.names + ngroups + nlevels ...

# 2. set/check options (eg meanstructure = TRUE)
#   lav_options_set() # lav_options.R

# 3. check the (raw) data, or the sample statistics
#   lavData()
lavdata <- lavaan:::lavData(data = HolzingerSwineford1939,
                            ov.names = c("x1", "x2", "x3"))
slotNames(lavdata)
lavdata@ov.names
lavdata@ngroups
lavdata@X[[1]] # the raw data in group/block 1

# 4. lavpartable: create the parameter table
# needs to know: how many groups, how many categories, ...
lavpartable <- lavaanify(model = FLAT, auto = TRUE)

# 4b. compute and store partable attributes (ov.names, ov.names.x, ...)
lavpta <- lav_partable_attributes(lavpartable)
lavpta$vnames$ov

# 5. lavsamplestats
# compute sample statistics (cov, mean, Gamma, ...)
lavsamplestats <- lavaan:::lav_samplestats_from_data(lavdata,
                                                     lavoptions = lavoptions)
slotNames(lavsamplestats)
lavsamplestats@cov[[1]] # observed covariance matrix first group/block
lavsamplestats@icov[[1]] # inverse observed covariance matrix
lavsamplestats@cov.log.det[[1]] # log determinant covariance matrix

# 6. lavh1
# summary statistics of the 'unrestricted' (h1) model
# with complete data, this is trivial (cov, mean)
# when data is missing, we need to estimate cov/mean using EM
lavh1 <- lavaan:::lav_h1_implied_logl(lavdata = lavdata,
                             lavsamplestats   = lavsamplestats,
                             lavoptions       = lavoptions)
lavh1$implied

# 7. parameter bounds (needs lavh1$implied)
# only used if bounded estimation is requested
lavoptions$bounds <- "standard"
lavoptions$optim.bounds <- list(lower = c("ov.var", "loadings"),
                                upper = c("ov.var", "loadings"),
                                min.reliability.marker = 0.1)
lavpartable <- lavaan:::lav_partable_add_bounds(partable = lavpartable,
                lavh1 = lavh1, lavdata = lavdata, lavsamplestats = lavsamplestats,
                lavoptions = lavoptions)
# remove bounds again to save space
lavpartable$lower <- NULL
lavpartable$upper <- NULL

# 8. compute starting values
lavpartable$start <- lavaan:::lav_start(start.method = lavoptions$start,
                               lavpartable     = lavpartable,
                               lavsamplestats  = lavsamplestats)
lavpartable

# 9. lavmodel: create internal model representation (with GLIST)
lavmodel <- lavaan:::lav_model(lavpartable      = lavpartable,
                               lavoptions       = lavoptions)
lavmodel@GLIST

# 10. lavcache: compute some additional summary statistis
# only used when estimator = "PML", and "MML" (for now)
lavcache <- list()

# 11. estimation
# - default: lav_model_estimate() + nlminb() (quasi-Newton optimization)
# - lav_optim_gn(): Gauss-Newton optimization
# - lav_optim_noniter(): non-iterative procedures
# - lav_mvnorm_cluster_em_h0: EM for multilevel models
lavaan:::lav_verbose(TRUE) # be verbose
x <- try(lavaan:::lav_model_estimate(lavmodel        = lavmodel,
                            lavpartable     = lavpartable,
                            lavsamplestats  = lavsamplestats,
                            lavdata         = lavdata,
                            lavoptions      = lavoptions,
                            lavcache        = lavcache))
lavaan:::lav_verbose(FALSE) # switch off verbose

# store parameters in lavmodel
lavmodel <- lav_model_set_parameters(lavmodel, x = as.numeric(x))

# store parameters in @ParTable$est
lavpartable$est <- lav_model_get_parameters(lavmodel = lavmodel,
                                            type = "user", extra = TRUE)
lavpartable

# 12. lavimplied + lavloglik
# store model implied statistics in @implied
# if likelihood-based method, store also loglik in @loglik
lavimplied <- lav_model_implied(lavmodel)
lavloglik <- lavaan:::lav_model_loglik(lavdata        = lavdata,
                              lavsamplestats = lavsamplestats,
                              lavimplied     = lavimplied,
                              lavmodel       = lavmodel,
                              lavoptions     = lavoptions)

# 13. compute 'vcov': the variance matrix of the (free) parameters
# this is needed to compute standard errors
VCOV <- lavaan:::lav_model_vcov(lavmodel        = lavmodel,
                       lavsamplestats  = lavsamplestats,
                       lavoptions      = lavoptions,
                       lavdata         = lavdata,
                       lavpartable     = lavpartable,
                       lavcache        = lavcache,
                       lavimplied      = lavimplied,
                       lavh1           = lavh1)
VCOV

# prepare lavvcov slot
lavvcov <- list(se = lavoptions$se, information = lavoptions$information,
                vcov = VCOV)

# store standard errors in parameter table
lavpartable$se <- lavaan:::lav_model_vcov_se(lavmodel = lavmodel,
                                    lavpartable = lavpartable,
                                    VCOV = VCOV)

# 14. compute global test statistic (chi-square)
# trivial for standard test (=N * F_ML)
# more work for 'robust' test statistics (eg test = "satorra-bentler")
TEST <- lavaan:::lav_model_test(lavmodel            = lavmodel,
                       lavpartable         = lavpartable,
                       lavsamplestats      = lavsamplestats,
                       lavimplied          = lavimplied,
                       lavh1               = lavh1,
                       lavoptions          = lavoptions,
                       x                   = x,
                       VCOV                = VCOV,
                       lavdata             = lavdata,
                       lavcache            = lavcache,
                       lavloglik           = lavloglik)
lavtest <- TEST

# 14bis. lavfit
# store 'fit'information
# no longer used, but if I remove it, a dozen (old) packages break...

# 15. fit baseline model
fit.indep <- try(lavaan:::lav_object_independence(object = NULL,
                                         lavsamplestats = lavsamplestats,
                                         lavdata        = lavdata,
                                         lavcache       = lavcache,
                                         lavoptions     = lavoptions,
                                         lavh1          = lavh1), silent = TRUE)

# 16. rotation
# only needed if efa() blocks are defined
# lavmodel <- lav_model_efa_rotate(lavmodel   = lavmodel,
#                                  lavoptions = lavoptions)

# 17. create lavaan object

# don't run, as some pieces have not been created...
# lavaan <- new("lavaan",
#               version      = as.character(packageVersion("lavaan")),
#               call         = mc,                  # match.call
#               timing       = timing,              # list
#               Options      = lavoptions,          # list
#               ParTable     = lavpartable,         # list
#               pta          = lavpta,              # list
#               Data         = lavdata,             # S4 class
#               SampleStats  = lavsamplestats,      # S4 class
#               Model        = lavmodel,            # S4 class
#               Cache        = lavcache,            # list
#               Fit          = lavfit,              # S4 class
#               boot         = lavboot,             # list
#               optim        = lavoptim,            # list
#               implied      = lavimplied,          # list
#               loglik       = lavloglik,           # list
#               vcov         = lavvcov,             # list
#               test         = lavtest,             # list
#               h1           = lavh1,               # list
#               baseline     = lavbaseline,         # list
#               internal     = list(),              # empty list
#               external     = list()               # empty list
# )

# 18. post-fitting check of parameters
# lavInspect(lavaan, "post.check")

# the sem/cfa/growth function just set some
# options to user-friendly settings:

# default options for sem/cfa call
# mc$int.ov.free     = TRUE
# mc$int.lv.free     = FALSE
# #mc$auto.fix.first  = !std.lv
# mc$auto.fix.first  = TRUE # (re)set in lav_options_set
# mc$auto.fix.single = TRUE
# mc$auto.var        = TRUE
# mc$auto.cov.lv.x   = TRUE
# mc$auto.cov.y      = TRUE
# mc$auto.th         = TRUE
# mc$auto.delta      = TRUE
# mc$auto.efa        = TRUE


# PART 3: extractor functions
fit <- sem(model = 'f =~ x1 + x2 + x3 + x4', data = HolzingerSwineford1939)

parameterEstimates(fit)
# = subset of parTable(fit), but with additional columsn (z, pvalues, ...)

parameterEstimates(fit, output = "text")
# this is a big part of the summary() output

# summary()
# first creates summary output as a list:
out <- summary(fit)
names(out)
class(out)
# with a specific print function
out

# fit indices
fitMeasures(fit)
fitMeasures(fit, output = "matrix")
fitMeasures(fit, output = "text")

# lavInspect/lavTech
lavInspect(fit, "est")
lavTech(fit, "est")


# PART 4: lavaan slots

fit <- sem(model = 'f =~ x1 + x2 + x3 + x4', data = HolzingerSwineford1939)
class(fit)
slotNames(fit)

# 1. lavaan version used to create this object
fit@version

# 2. user-specified call
fit@call

# 3. timings of several substeps
unlist(fit@timing)

# 4. options used for this object
unlist(fit@Options)

# 5. the parameter table
fit@ParTable  # list
parTable(fit) # return as data.frame

# 6. the parameter table attributes
names(fit@pta)
fit@pta$vnames$ov
fit@pta$vidx$ov
fit@pta$nfac
fit@pta$nblocks

# 7. Data slot (S4)
fit@Data # has its own print function
slotNames(fit@Data)
as.data.frame(fit@Data@ov)
str(fit@Data)

# 8. SampleStats (S4)
fit@SampleStats
slotNames(fit@SampleStats)
fit@SampleStats@cov[[1]] # list with element per group

# 9. Model (S4)
slotNames(fit@Model)
fit@Model@x.free.idx # parameter index in parameter table
fit@Model@m.free.idx # parameter index in model matrix

# 10. Cache (list)
# cached information, only used for estimator PML and MML (for now)

# 11. Fit
# deprecated, only kept to avoid breaking some (old) packages (eg rsem)

# 12. boot (list)
# only used if bootstrapping was used
fitb <- sem(model = 'f =~ x1 + x2 + x3 + x4', data = HolzingerSwineford1939,
            se = "bootstrap", bootstrap = 100L, verbose = TRUE)
head(fitb@boot$coef)

# 13. optim -- info about (iterative) optimization process
str(fit@optim)

# 14. loglik -- loglikelihood information (ML estimation only)
unlist(fit@loglik)

# 15. implied -- implied sample statistics (per group)
fit@implied

# 16. vcov -- variance-covariance matrix of the (free) model parameters
fit@vcov

# 17. test -- global test statistics
fit2 <- sem(model = 'f =~ x1 + x2 + x3 + x4', data = HolzingerSwineford1939,
            test = "satorra-bentler")
names(fit2@test)
fit2@test$satorra.bentler

# 18. h1 -- sample statistics + logl unrestricted/saturated model
# often similar to info in @SampleStats, but not if missing data,
# multilevel data, ...
fit@h1$implied # this is what is used for lavInspect(fit, "implied")
fit@h1$logl

# 19. baseline -- information about the baseline model (needed for CFI/TLI)
names(fit@baseline)
as.data.frame(fit@baseline$partable)
fit@baseline$test$standard

# 20. internal -- list to store specials flags/info, for internal use only

# 21. external -- list to store specials flags/info if you are an external
# developer


# PART 5: source code structure

# 5 main functions:

# 1) xxx_lavaan.R
# 2) xxx_lavaanList.R
# 3) xxx_efa.R (new in 0.6-13)
# 4) xxx_sam.R
# 5) xxx_fsr.R (only for repro reasons; access via lavaan:::fsr())

# mixture of new, old, and very old code
# very old code: functions do not start with lav_ prefix
# for example:
lavaan:::lav_model_sigma

# files that start with ctr_  contain contributed code
# written by others
# for example: ctr_pml_plrt.R (written by Myrsini Katsikatsou)
# (with only minor edits by YR)

# 00class.R contains S4 class definitions
# 00generic.R defines 'generic' functions (that can be specialized)
# fitMeasures(), lavInspect(), lavTech()

# zzz_OLDNAMES.R contains aliases for (very old) function names
# that are still used by some packages
# eg computeExpectedInformation <- lav_model_information_expected

# zzz.R traditionally contains the (in)famous startup message

# most files start with the lav_ prefix
# the second term often refers the type of object for which the file
# contains functions, for example
# lav_matrix.R
# lav_partable_subset.R
# lav_model_estimate.R
# lav_object_post_check.R

# but sometimes, it refers to what is created, or what is done
# lav_test.R  # creates @test slot
# lav_print.R # prints various objects

# for ML, an important set of functions are:
# lav_mvnorm.R
# lav_mvnorm_h1.R
# lav_mvnorm_missing.R
# lav_mvnorm_missing_h1.R
# lav_mvnorm_cluster.R
# lav_mvnorm_cluster_missing.R

# the standard ML discrepancy function is found at the top of
# lav_objective.R

# sometimes, lavaan needs to do trivial things (like regression), but
# needs a bit more info than what we get from standard functions (like lm):
# lav_uvreg.R          # univariate regression
# lav_uvord.R          # univariate probit regression
# lav_mvreg.R          # multivariate regression
# lav_mvreg_cluster.R  # multivariate twolevel regression

# during iterative estimation, we need to compute the value of the
# objective function (i.e., the discrepancy function) many times

# lav_model_estimate.R # runs the iterations
# and defines the following objective function:

   # function to be minimized
    objective_function <- function(x, verbose = FALSE, infToMax = FALSE,
                                   debug = FALSE) {
        # 2. unpack
        if(lavmodel@eq.constraints) {
            x <- as.numeric(lavmodel@eq.constraints.K %*% x) +
                            lavmodel@eq.constraints.k0
        }

        # 1. unscale
        x <- x / parscale

        # update GLIST (change `state') and make a COPY!
        GLIST <- lav_model_x2GLIST(lavmodel, x = x)

        fx <- lav_model_objective(lavmodel       = lavmodel,
                                  GLIST          = GLIST,
                                  lavsamplestats = lavsamplestats,
                                  lavdata        = lavdata,
                                  lavcache       = lavcache,
                                  verbose        = verbose)

        # only for PML: divide by N (to speed up convergence)
        if(estimator == "PML") {
            fx <- fx / lavsamplestats@ntotal
        }

        if(debug || verbose) {
            cat("  objective function  = ",
                sprintf("%18.16f", fx), "\n", sep="")
        }
        if(debug) {
            cat("Current free parameter values =\n"); print(x); cat("\n")
        }

        if(lavoptions$optim.partrace) {
            PENV$PARTRACE <- rbind(PENV$PARTRACE, c(fx, x))
        }

        # for L-BFGS-B
        #if(infToMax && is.infinite(fx)) fx <- 1e20
        if(!is.finite(fx)) {
            fx.group <- attr(fx, "fx.group")
            fx <- 1e20
            attr(fx, "fx.group") <- fx.group # only for lav_model_fit()
        }

        fx
    }

# lav_model_objective() can be found in lav_model_objective.R
# 1) compute model implied summary statistics (for each group)
#    using eg lav_model_sigma()
# 2) compute value discrepancy function
#    eg estimator.GLS() or estimator.ML() # see lav_objective.R
# 3) if multiple groups, combine the values using group weights
# 4) return value (fx)











