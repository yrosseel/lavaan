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

    # set baseline/h1 to FALSE
    lavoptions$h1 <- FALSE
    lavoptions$baseline <- FALSE
    lavoptions$loglik <- TRUE # eg for multilevel
    lavoptions$implied <- TRUE #, needed for loglik
    lavoptions$check.start <- FALSE
    lavoptions$check.gradient <- FALSE
    lavoptions$check.post <- FALSE
    lavoptions$check.vcov <- FALSE

    # ALWAYS do.fit
    lavoptions$do.fit  <- TRUE

    # verbose?
    lavoptions$verbose <- verbose

    # warn?
    lavoptions$warn <- warn

    # needed?
    if(any(lavpartable$op == "~1")) lavoptions$meanstructure <- TRUE

    # FIXME: it is crucial that the order of the ov's, as returned by
    # lavNames() remains the same
    # so lavNames(object) should equal lavNames(lavpartable)
    # otherwise, we will use the wrong sample statistics!!!
    #
    # this seems ok now, because we first generate the covariances in
    # lavpartable, and they should be in the right order (unlike the
    # intercepts)

    if(.hasSlot(object, "h1"))  {
        lavh1 <- object@h1
    } else {
        lavh1 <- lav_h1_logl(lavdata = object@Data,
                             lavsamplestats = object@SampleStats,
                             lavoptions = object@Options)
    }

    FIT <- lavaan(lavpartable,
                  slotOptions     = lavoptions,
                  slotSampleStats = object@SampleStats,
                  slotData        = object@Data,
                  slotCache       = object@Cache,
                  sloth1          = lavh1)

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

    if(.hasSlot(object, "h1"))  {
        lavh1 <- object@h1
    } else {
        lavh1 <- lav_h1_logl(lavdata = object@Data,
                             lavsamplestats = object@SampleStats,
                             lavoptions = object@Options)
    }

    FIT <- lavaan(lavpartable,
                  slotOptions     = lavoptions,
                  slotSampleStats = object@SampleStats,
                  slotData        = object@Data,
                  slotCache       = object@Cache,
                  sloth1          = lavh1)

    FIT
}


# 3. extended model
lav_object_extended <- function(object, add = NULL,
                                remove.duplicated = TRUE,
                                all.free = FALSE,
                                verbose = FALSE, warn = FALSE,
                                do.fit = FALSE) {

    # partable original model
    partable <- object@ParTable[c("lhs", "op", "rhs", "free", "exo", "label",
                                 "plabel")]

    # new in 0.6-3: check for non-parameters
    nonpar.idx <- which(partable$op %in% c("==", ":=", "<", ">"))

    # always add block/group/level
    if(!is.null(object@ParTable$group)) {
        partable$group <- object@ParTable$group
    } else {
        partable$group <- rep(1L, length(partable$lhs))
        if(length(nonpar.idx) > 0L) {
            partable$group[nonpar.idx] <- 0L
        }
    }
    if(!is.null(object@ParTable$level)) {
        partable$level <- object@ParTable$level
    } else {
        partable$level <- rep(1L, length(partable$lhs))
        if(length(nonpar.idx) > 0L) {
            partable$level[nonpar.idx] <- 0L
        }
    }
    if(!is.null(object@ParTable$block)) {
        partable$block <- object@ParTable$block
    } else {
        partable$block <- rep(1L, length(partable$lhs))
        if(length(nonpar.idx) > 0L) {
            partable$block[nonpar.idx] <- 0L
        }
    }

    # TDJ: Added to prevent error when lav_partable_merge() is called below.
    #      Problematic if object@ParTable is missing one of the requested slots,
    #      which returns a NULL slot with a missing <NA> name.  For example:
    #        example(cfa)
    #        lav_partable_independence(lavdata = fit@Data, lavpta = fit@pta,
    #                                  lavoptions = lavInspect(fit, "options"))
    #     Has no "label" or "plabel" elements.
    empties <- which(sapply(partable, is.null))
    if(length(empties)) {
        partable[empties] <- NULL
    }

    if(all.free) {
        partable$user <- rep(1L, length(partable$lhs))
        non.free.idx <- which(partable$free == 0L & partable$op != "==" &
                              partable$op != ":=" & partable$op != "<" &
                              partable$op != ">")
        partable$free[ non.free.idx ] <- 1L
        partable$user[ non.free.idx ] <- 10L
    }

    # replace 'start' column, since lav_model will fill these in in GLIST
    partable$start <- parameterEstimates(object, remove.system.eq = FALSE,
                                         remove.def = FALSE,
                                         remove.eq = FALSE,
                                         remove.ineq = FALSE)$est

    # add new parameters, extend model
    if(is.list(add)) {
        stopifnot(!is.null(add$lhs),
                  !is.null(add$op),
                  !is.null(add$rhs))
        ADD <- add
    } else if(is.character(add)) {
        ngroups <- lav_partable_ngroups(partable)
        ADD.orig <- lavaanify(add, ngroups = ngroups)
        ADD <- ADD.orig[,c("lhs","op","rhs","user","label")] # minimum

        # always add block/group/level
        if(!is.null(ADD.orig$group)) {
            ADD$group <- ADD.orig$group
        } else {
            ADD$group <- rep(1L, length(ADD$lhs))
        }
        if(!is.null(ADD.orig$level)) {
            ADD$level <- ADD.orig$level
        } else {
            ADD$level <- rep(1L, length(ADD$lhs))
        }
        if(!is.null(ADD.orig$block)) {
            ADD$block <- ADD.orig$block
        } else {
            ADD$block <- rep(1L, length(ADD$lhs))
        }

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
    lavoptions$verbose <- verbose

    # warn?
    lavoptions$warn <- warn

    # do.fit?
    lavoptions$do.fit <- do.fit

    # needed?
    if(any(LIST$op == "~1")) lavoptions$meanstructure <- TRUE

    if(.hasSlot(object, "h1"))  {
        lavh1 <- object@h1
    } else {
        # old object -- for example 'usemmodelfit' in package 'pompom'

        # add a few fields
        lavoptions$h1 <- FALSE
        lavoptions$implied <- FALSE
        lavoptions$baseline <- FALSE
        lavoptions$loglik <- FALSE

        # add a few slots
        object@Data@weights <- vector("list", object@Data@ngroups)
        object@Model@estimator <- object@Options$estimator

        lavh1 <- lav_h1_logl(lavdata = object@Data,
                             lavsamplestats = object@SampleStats,
                             lavoptions = object@Options)
    }

    FIT <- lavaan(LIST,
                  slotOptions     = lavoptions,
                  slotSampleStats = object@SampleStats,
                  slotData        = object@Data,
                  slotCache       = object@Cache,
                  sloth1          = lavh1)

    FIT
}
