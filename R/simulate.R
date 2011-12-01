# simulate data starting from a user-specified model
#
# initial version: YR 24 jan 2011
# revision for 0.4-11: YR 21 okt 2011
simulateData <- function(
                         # user-specified model syntax
                         model.syntax = '',
                         model.type      = "sem",

                         # model modifiers
                         meanstructure   = "default",
                         int.ov.free     = TRUE,
                         int.lv.free     = FALSE,
                         fixed.x         = "default", # or FALSE?
                         orthogonal      = FALSE,
                         std.lv          = FALSE,

                         auto.fix.first  = !std.lv,
                         auto.fix.single = TRUE,
                         auto.var        = TRUE,
                         auto.cov.lv.x   = TRUE,
                         auto.cov.y      = TRUE,
                         ...,
                         sample.nobs = 500L,
                         group.label = paste("G", 1:ngroups, sep=""),
                         seed = NULL,
                         empirical = FALSE,
                         return.type = "data.frame"
                        )
{
    if(!is.null(seed)) set.seed(seed)
    if(!exists(".Random.seed", envir = .GlobalEnv))
        runif(1)               # initialize the RNG if necessary
    RNGstate <- .Random.seed

    # run lavaan to set up the model matrices
    fit <- lavaan(model.syntax, 
                  meanstructure=meanstructure, 
                  int.ov.free=int.ov.free, int.lv.free=int.lv.free,
                  fixed.x=fixed.x,
                  orthogonal=orthogonal,
                  std.lv=std.lv,
                  auto.fix.first=auto.fix.first,
                  auto.fix.single=auto.fix.single,
                  auto.var=auto.var,
                  auto.cov.lv.x=auto.cov.lv.x,
                  auto.cov.y=auto.cov.y,
                  sample.nobs=sample.nobs,
                  ...)

    # the model-implied moments for the population
    Sigma.hat <- computeSigmaHat(fit@Model)
       Mu.hat <- computeMuHat(fit@Model)

    # ngroups
    ngroups <- length(sample.nobs)

    # prepare 
    X <- vector("list", length=ngroups)
    out <- vector("list", length=ngroups)

    for(g in 1:ngroups) {
        # FIXME: change to rmvnorm once we include the library?
        X[[g]] <- MASS.mvrnorm(n = sample.nobs[g],
                               mu = Mu.hat[[g]],
                               Sigma = Sigma.hat[[g]],
                               empirical = empirical)

        if(return.type == "data.frame") X[[g]] <- as.data.frame(X[[g]])
    }

    if(return.type == "matrix") {
        if(ngroups == 1L) {
            return(X[[1L]])
        } else {
            return(X)
        }

    } else if (return.type == "data.frame") {
        Data <- X[[1L]]

        # if multiple groups, add group column
        if(ngroups > 1L) {
            for(g in 2:ngroups) {
                Data <- rbind(Data, X[[g]])
            }
            Data$group <- rep(1:ngroups, times=sample.nobs)
        }
        var.names <- vnames(fit@User, type="ov", group=1L)
        if(ngroups > 1L) var.names <- c(var.names, "group")
        names(Data) <- var.names
        return(Data)

    } else if (return.type == "cov") {
        if(ngroups == 1L) {
            return(cov(X[[1L]]))
        } else {
            cov.list <- lapply(X, cov)
            return(cov.list)
        }
    }
}
