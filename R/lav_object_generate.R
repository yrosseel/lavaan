# here, we generate new models based on the original model in lavobject
# 1. the independence model
# 2. the unrestricted model
# 3. ...


# 1. fit an 'independence' model 
#    note that for ML (and ULS and DWLS), the 'estimates' of the 
#    independence model are simply the observed variances
#    but for GLS and WLS, this is not the case!!
lav_object_independence <- function(object, se = FALSE, verbose = FALSE,
                                    warn = FALSE) {

    # construct parameter table for independence model
    lavpartable <- lav_partable_independence(object)

    # adapt options
    lavoptions <- object@Options

    # se
    if(se) {
        if(lavoptions$se == "none") {
            lavoptions$se <- "standard"
        }
    } else {
        ## FIXME: if test = scaled, we need it anyway?
        lavoptions$se <- "none" 
    }
   
    # ALWAYS do.fit
    lavoptions$do.fit  <- TRUE

    # verbose?
    if(verbose) {
        lavoptions$verbose <- TRUE
    } else {
        lavoptions$verbose <- FALSE
    }

    # warn?
    if(warn) {
        lavoptions$warn <- TRUE
    } else {
        lavoptions$warn <- FALSE
    }


    # needed?
    if(any(lavpartable$op == "~1")) lavoptions$meanstructure <- TRUE

    FIT <- lavaan(lavpartable,  
                  slotOptions     = lavoptions,
                  slotSampleStats = object@SampleStats,
                  slotData        = object@Data,
                  slotCache       = object@Cache)

    FIT
}


# 2. unrestricted model
lav_object_unrestricted <- function(object, se = FALSE, verbose = FALSE,
                                    warn = FALSE) {

    # construct parameter table for unrestricted model
    lavpartable <- lav_partable_unrestricted(object)

    # adapt options
    lavoptions <- object@Options

    # se
    if(se) {
        if(lavoptions$se == "none") {
            lavoptions$se <- "standard"
        }
    } else {
        ## FIXME: if test = scaled, we need it anyway?
        lavoptions$se <- "none" 
    }
   
    # ALWAYS do.fit
    lavoptions$do.fit  <- TRUE

    # verbose?
    if(verbose) {
        lavoptions$verbose <- TRUE
    } else {
        lavoptions$verbose <- FALSE
    }

    # warn?
    if(warn) {
        lavoptions$warn <- TRUE
    } else {
        lavoptions$warn <- FALSE
    }

    # needed?
    if(any(lavpartable$op == "~1")) lavoptions$meanstructure <- TRUE

    FIT <- lavaan(lavpartable,
                  slotOptions     = lavoptions,
                  slotSampleStats = object@SampleStats,
                  slotData        = object@Data,
                  slotCache       = object@Cache)

    FIT
}
