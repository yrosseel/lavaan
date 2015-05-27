# here, we generate new models based on the original model in lavobject
# 1. the independence model
# 2. the unrestricted model
# 3. model + extra parameters (for modindices/lavTestScore)


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


# 3. extended model
lav_object_extended <- function(object, add = NULL,
                                remove.duplicated = TRUE,
                                all.free = FALSE,
                                verbose = FALSE, warn = FALSE, 
                                do.fit = FALSE) {

    # partable original model
    partable <- object@ParTable[c("lhs","op","rhs","group","free",
                                  "exo","label","plabel")] # do we need 'exo'?
    if(all.free) {
        partable$user <- rep(1L, length(partable$lhs))
        non.free.idx <- which(partable$free == 0L)
        partable$free[ non.free.idx ] <- 1L
        partable$user[ non.free.idx ] <- 10L
    }
 
    # replace 'start' column, since lav_model will fill these in in GLIST
    partable$start <- parameterEstimates(object,
                          remove.eq = FALSE, remove.ineq = FALSE)$est

    # add new parameters, extend model
    if(is.list(add)) {
        stopifnot(!is.null(add$lhs),
                  !is.null(add$op),
                  !is.null(add$rhs))
        ADD <- add
    } else if(is.character(add)) {
        ADD <- lavaanify(add, ngroups = object@Data@ngroups)
        ADD <- ADD[,c("lhs","op","rhs","group","user","label")]
        remove.idx <- which(ADD$user == 0)
        if(length(remove.idx) > 0L) {
            ADD <- ADD[-remove.idx,]
        }
        ADD$start <- rep( 0, nrow(ADD))
        ADD$free  <- rep( 1, nrow(ADD))
        ADD$user  <- rep(10, nrow(ADD))
    }

    # merge
    LIST <- lav_partable_merge(partable, ADD,
                               remove.duplicated = remove.duplicated,
                               warn = FALSE)

    # remove nonpar?
    #if(remove.nonpar) {
    #    nonpar.idx <- which(LIST$op %in% c("==", ":=", "<", ">"))
    #    if(length(nonpar.idx) > 0L) {
    #        LIST <- LIST[-nonpar.idx,]
    #    }
    #}

    # redo 'free'
    free.idx <- which(LIST$free > 0)
    LIST$free[free.idx] <- 1:length(free.idx)
    
    # adapt options
    lavoptions <- object@Options

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
    if(any(LIST$op == "~1")) lavoptions$meanstructure <- TRUE

    FIT <- lavaan(LIST,
                  do.fit          = do.fit,
                  slotOptions     = lavoptions,
                  slotSampleStats = object@SampleStats,
                  slotData        = object@Data,
                  slotCache       = object@Cache)

    FIT
}
