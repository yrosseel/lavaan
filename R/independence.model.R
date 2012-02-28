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
        OV.X <- vnames(no.fit@User, type="ov.x", group=1L)
    model.syntax <- 
        syntax.independence.model(ov.names   = no.fit@Data@ov.names[[1L]],
                                  ov.names.x = OV.X,
                                  sample.cov = no.fit@Sample@cov)
    # refit
    lavaan <- cfa(model=model.syntax, ...)

    lavaan
}


# starting from a fitted lavaan object
independence.model.fit2 <- function(object) {

    # construct syntax for independence model
    OV.X <- character(0L)
    if(object@Options$mimic %in% c("lavaan", "Mplus"))
        OV.X <- vnames(object@User, type="ov.x", group=1L)
    model.syntax <- 
        syntax.independence.model(ov.names   = object@Data@ov.names[[1L]],
                                  ov.names.x = OV.X,
                                  sample.cov = object@Sample@cov)

    # refit
    lavaan <- update(object, model = model.syntax)

    lavaan
}

#### FIXME!!!!! we should get rid of this function.... !
#### FIXME: will not work if two different models in multiple groups
independence.model.fit <- function(object) {

    mc <- match.call()
    timing <- list()

    # construct parameter Table for independence model
    OV.X <- lapply(as.list(1:object@Sample@ngroups),
                   function(x) vnames(object@User, type="ov.x", x))

    # construct
    lavaanUser <- independenceModel(ov.names   = object@Data@ov.names,
                                    ov.names.x = OV.X,
                                    sample.cov = object@Sample@cov,
                                    meanstructure = object@Model@meanstructure,
                                    sample.mean = object@Sample@mean,
                                    fixed.x    = object@Model@fixed.x)
    # fit?
    do.fit <- TRUE

    # 1. lavaanOptions
    lavaanOptions <- object@Options
    lavaanOptions$se      <- "none"
    lavaanOptions$do.fit  <- do.fit
    lavaanOptions$verbose <- FALSE
    lavaanOptions$warn    <- FALSE

    # 2b. change meanstructure flag?
    if(any(lavaanUser$op == "~1")) lavaanOptions$meanstructure <- TRUE

    # 3. 
    lavaanData             <- object@Data
    lavaanSampleStats      <- object@Sample

    # 4. 
    lavaanStart <-
        StartingValues(user       = lavaanUser,
                       sample     = lavaanSampleStats,
                       model.type = lavaanOptions$model.type,
                       debug      = lavaanOptions$debug)

    # 5. 
    lavaanModel <-
        Model(user           = lavaanUser,
              start          = lavaanStart,
              representation = lavaanOptions$representation,
              debug          = lavaanOptions$debug)

    # 6.
    x <- VCOV <- TEST <- NULL
    if(do.fit) {
        x <- estimateModel(lavaanModel,
                           sample  = lavaanSampleStats,
                           options = lavaanOptions)
                           # control???
        lavaanModel <- setModelParameters(lavaanModel, x = x)
        if(!is.null(attr(x, "con.jac")))
            lavaanModel@con.jac <- attr(x, "con.jac")
    }

    # 7.
    
    # 8.
    TEST <- computeTestStatistic(lavaanModel,
                                 user    = lavaanUser,
                                 sample  = lavaanSampleStats,
                                 options = lavaanOptions,
                                 x       = x,
                                 VCOV    = VCOV,
                                 data    = lavaanData)

    # 9. collect information about model fit (S4)
    lavaanFit <- Fit(user  = lavaanUser,
                     start = lavaanStart,
                     model = lavaanModel,
                     x     = x,
                     VCOV  = VCOV,
                     TEST  = TEST)

    # 10. construct lavaan object
    lavaan <- new("lavaan",
                  call    = mc,                     # match.call
                  timing  = timing,                 # list
                  Options = lavaanOptions,          # list
                  User    = lavaanUser,             # list
                  Data    = lavaanData,             # S3 class
                  Sample  = lavaanSampleStats,      # S4 class
                  Model   = lavaanModel,            # S4 class
                  Fit     = lavaanFit               # S4 class
                 )

    lavaan
}

