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

#### FIXME!!!!! we should get rid of this function.... !
#### FIXME: will not work if two different models in multiple groups
independence.model.fit <- function(object) {

    OV.X <- lapply(as.list(1:object@Data@ngroups),
                       function(x) vnames(object@ParTable, type="ov.x", x))

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
    lavpartable <- 
        lav_partable_independence(ov.names      = object@Data@ov.names,
                                  ov            = object@Data@ov,
                                  ov.names.x    = OV.X,
                                  sample.cov    = object@SampleStats@cov,
                                  meanstructure = object@Model@meanstructure,
                                  sample.mean   = object@SampleStats@mean,
                                  sample.th     = object@SampleStats@th,
                                  parameterization = object@Options$parameterization,
                                  fixed.x       = FIXED.X)

    # adapt options
    lavoptions <- object@Options
    lavoptions$se      <- "none" ## FIXME: if test = scaled, we need it anyway?
    lavoptions$do.fit  <- TRUE
    lavoptions$verbose <- FALSE
    lavoptions$warn    <- FALSE
    if(any(lavpartable$op == "~1")) lavoptions$meanstructure <- TRUE

    FIT <- lavaan(lavpartable,  
                  slotOptions     = lavoptions,
                  slotSampleStats = object@SampleStats,
                  slotData        = object@Data,
                  slotCache       = object@Cache)

    FIT
}

