# classic score test (= Langrage Multiplier test)
#
# 'extra' contains extra parameters that are currently fixed (to zero),
# but should be released
#
# NOTE: - this is just a 'test' implementation
#       - it can not handle non-identified parameters,
#         redundant parameters, etc...
#
# YR 1 dec 2014

lavTestScore <- function(object, extra = NULL, release = NULL, 
                         verbose = FALSE) {

    if(object@Fit@npar > 0L && !object@Fit@converged)
        stop("lavaan ERROR: model did not converge")

    if(is.null(extra) || nchar(extra) == 0L) {
        stop("lavaan ERROR: extra parameter syntax is empty")
    }

    partable <- parTable(object)
    partable$start <- parameterEstimates(fit, remove.eq = FALSE, 
                                         remove.ineq = FALSE)$est
    nid <- max(partable$id)
    nfree <- max(partable$free)

    # parse extra syntax
    FLAT <- lavParseModelString(extra, as.data.frame = TRUE)
    FLAT$mod.idx <- FLAT$fixed <- FLAT$start <- FLAT$label <- NULL
    nflat <- nrow(FLAT)
    if(is.null(FLAT$group)) {
        FLAT$group <- rep(1L, nflat)
    }

    # check for duplicated elements
    TMP <- rbind(FLAT[,    c("lhs","op","rhs","group")], 
                 partable[,c("lhs","op","rhs","group")])
    idx <- which(duplicated(TMP, fromLast=TRUE)) # idx should be in FLAT
    if(length(idx)) {
        warning("lavaan WARNING: parameters already in the model are ignored:\n", 
        paste(apply(FLAT[idx, c("lhs","op","rhs")], 1, 
              paste, collapse=" "), collapse="\n"))
        FLAT <- FLAT[-idx,]
    }

    if(nrow(FLAT) == 0L) {
        stop("lavaan ERROR: extra parameters table is empty")
    }

    # what about group?
    nflat <- nrow(FLAT)
    FLAT$id <- (nid+1L):(nid+nflat)
    FLAT$start <- 0
    FLAT$free <- (nfree + 1L):(nfree + nflat)
    FLAT$label <- rep("", nflat)
    FLAT$plabel <- paste("p", FLAT$id, "__", sep = "")
    FLAT$user <- rep(0L, nflat)
    #FLAT$unco <- (max(partable$unco) + 1L):(max(partable$unco) + nflat)

    # build new 'extended' model
    NEW <- base::merge(partable, FLAT, all = TRUE)


    # fit model, without any iterations
    fit <- lavaan(NEW, 
                  start = NEW,
                  slotOptions = object@Options,
                  slotSampleStats = object@SampleStats,
                  slotData = object@Data,
                  do.fit = FALSE)

    dx <- lav_model_gradient(lavmodel = fit@Model, 
                             lavsamplestats = fit@SampleStats,
                             lavdata = fit@Data,
                             lavcache = fit@Cache,
                             type = "free",
                             estimator = fit@Options$estimator)

    lavoptions <- fit@Options
    lavoptions$se <- object@Options$se

    J <- lav_model_vcov(lavmodel = fit@Model, 
                        lavoptions = lavoptions,
                        lavpartable = fit@ParTable,
                        # control??
                        lavsamplestats = fit@SampleStats,
                        lavdata = fit@Data,
                        lavcache = fit@Cache)

    J <- J * nobs(fit)
    stat <- nobs(fit) * as.numeric(t(dx) %*% J %*% dx)
    df <- nrow(FLAT)
    pvalue <- 1 - pchisq(stat, df=df)

    list(stat = stat, df = df, p.value = pvalue, se = lavoptions$se)
}
