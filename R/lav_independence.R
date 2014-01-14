# fit an 'independence' model 
# useful for CFI/TLI fit measures
# for now, we only compute the 'default' independence model where
# the diagonal elements of Sigma are freely estimable

# note that for ML (and ULS and DWLS), the 'estimates' of the 
# independence model are simply the observed variances
# but for GLS and WLS, this is not the case!!

# starting from a model.syntax
independence.model <- function(model.syntax = '', ...) {

    # process
    no.fit <- cfa(model=model.syntax, ..., data=NULL, sample.cov=NULL)

    # reconstruct model.syntax...
    OV.X <- character(0L)
    if(no.fit@Options$mimic %in% c("lavaan", "Mplus"))
        OV.X <- vnames(no.fit@ParTable, type="ov.x", group=1L)
    model.syntax <- 
        syntax.independence.model(ov.names   = no.fit@Data@ov.names[[1L]],
                                  ov.names.x = OV.X,
                                  sample.cov = no.fit@SampleStats@cov)
    # refit
    lavaan <- cfa(model=model.syntax, ...)

    lavaan
}


# starting from a fitted lavaan object
independence.model.fit2 <- function(object) {

    # construct syntax for independence model
    OV.X <- character(0L)
    if(object@Options$mimic %in% c("lavaan", "Mplus"))
        OV.X <- vnames(object@ParTable, type="ov.x", group=1L)
    model.syntax <- 
        syntax.independence.model(ov.names   = object@Data@ov.names[[1L]],
                                  ov.names.x = OV.X,
                                  sample.cov = object@SampleStats@cov)

    # refit
    lavaan <- update(object, model = model.syntax)

    lavaan
}

#### FIXME!!!!! we should get rid of this function.... !
#### FIXME: will not work if two different models in multiple groups
independence.model.fit <- function(object) {

    mc <- match.call()
    timing <- list()
    categorical <- object@Model@categorical

    # construct parameter Table for independence model
    #if(!categorical) {
        OV.X <- lapply(as.list(1:object@Data@ngroups),
                       function(x) vnames(object@ParTable, type="ov.x", x))
    #} else {
    #    OV.X <- NULL
    #}

    # what with fixed.x?
    if(object@Options$mimic %in% c("lavaan", "Mplus")) {
        FIXED.X = object@Model@fixed.x
    } else if(object@Options$mimic == "EQS") {
        # always ignore fixed.x
        OV.X = NULL
        FIXED.X = FALSE
    } else if(object@Options$mimic == "LISREL") {
        # always ignore fixed.x??? CHECKME!!
        OV.X = NULL
        FIXED.X = FALSE
    }

    # construct
    lavaanParTable <- 
        lav_partable_independence(ov.names      = object@Data@ov.names,
                                  ov            = object@Data@ov,
                                  ov.names.x    = OV.X,
                                  sample.cov    = object@SampleStats@cov,
                                  meanstructure = object@Model@meanstructure,
                                  sample.mean   = object@SampleStats@mean,
                                  sample.th     = object@SampleStats@th,
                                  parameterization = object@Options$parameterization,
                                  fixed.x       = FIXED.X)
   
    # fit?
    do.fit <- TRUE

    # 1. lavaanOptions
    lavaanOptions <- object@Options
    lavaanOptions$se      <- "none"
    lavaanOptions$do.fit  <- do.fit
    lavaanOptions$verbose <- FALSE
    lavaanOptions$warn    <- FALSE

    # 2b. change meanstructure flag?
    if(any(lavaanParTable$op == "~1")) lavaanOptions$meanstructure <- TRUE

    # 3. 
    lavaanData             <- object@Data
    lavaanSampleStats      <- object@SampleStats

    # 4. 
    lavaanStart <-
        lav_start(partable    = lavaanParTable,
                  samplestats = lavaanSampleStats,
                  model.type  = lavaanOptions$model.type,
                  debug       = lavaanOptions$debug)
    lavaanParTable$start <- lavaanStart

    # 5. 
    lavaanModel <-
        lav_model(partable         = lavaanParTable,
                  representation   = lavaanOptions$representation,
                  th.idx           = lavaanSampleStats@th.idx,
                  parameterization = lavaanOptions$parameterization,
                  debug            = lavaanOptions$debug)

    # cache
    lavaanCache <- object@Cache

    # 6.
    x <- VCOV <- TEST <- NULL
    if(do.fit) {
        x <- lav_model_estimate(lavaanModel,
                           samplestats  = lavaanSampleStats,
                           X            = object@Data@X,
                           cache        = lavaanCache,
                           options      = lavaanOptions)
                           # control???
        lavaanModel <- lav_model_set_parameters(lavaanModel, x = x,
                          estimator=lavaanOptions$estimator)
        if(!is.null(attr(x, "con.jac")))
            lavaanModel@con.jac <- attr(x, "con.jac")
    }

    # 7.
    
    # 8.
    # NOTE: Mplus 6 BUG??
    # - if estimator = WLSMV, baseline model is NOT using
    #   scaled.shifted, but mean.(var.)adusted!!
    test.options <- lavaanOptions
    #if(test.options$test == "scaled.shifted")
    #    test.options$test <- "mean.var.adjusted"
    TEST <- computeTestStatistic(lavaanModel,
                                 partable      = lavaanParTable,
                                 samplestats   = lavaanSampleStats,
                                 options       = test.options,
                                 x             = x,
                                 VCOV          = VCOV,
                                 cache         = lavaanCache,
                                 data          = lavaanData)

    # 9. collect information about model fit (S4)
    lavaanFit <- Fit(partable = lavaanParTable,
                     model    = lavaanModel,
                     x        = x,
                     VCOV     = VCOV,
                     TEST     = TEST)

    # 10. construct lavaan object
    lavaan <- new("lavaan",
                  call        = mc,                     # match.call
                  timing      = timing,                 # list
                  Options     = lavaanOptions,          # list
                  ParTable    = lavaanParTable,         # list
                  Data        = lavaanData,             # S3 class
                  SampleStats = lavaanSampleStats,      # S4 class
                  Model       = lavaanModel,            # S4 class
                  Cache       = lavaanCache,
                  Fit         = lavaanFit               # S4 class
                 )

    lavaan
}

