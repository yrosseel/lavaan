# This code is written by YR (using lavaan components), but based on 
# research code written by Mariska Barendse (Groningen/Amsterdam, NL)
#
# September 2013
#
# Three fit indices for the PML estimator (if all categorical, no exo)
# - Cp(max)
# - CF
# - CM

# FIXME: how to handle multiple groups??

# Mariska Barendse Cp statistic
#lav_tables_fit_Cp <- function(object, alpha = 0.05) {
#
#    out <- lavTablesFit(object, statistic = "LR", p.value = TRUE)
#
#    # Bonferonni adjusted p-value
#    ntests <- length(out$lhs)
#    out$alpha.adj <- alpha / ntests
#    #out$pval <- pchisq(out$LR, df=out$df, lower.tail = FALSE)
#
#    # remove LR.h0.pval
#    #out$LR.h0.pval <- NULL
#
#    out
#}

lav_tables_fit_Cp <- function(object, alpha = 0.05) {

    if(!all(object@Data@ov$type == "ordered")) {
        return(list(LR=as.numeric(NA), df=as.numeric(NA), 
               p.value=as.numeric(NA), p.value.Bonferroni=as.numeric(NA)))
    }

    TF <- lavTablesFit(object, statistic = "LR", p.value = TRUE)

    # Bonferonni adjusted p-value
    ntests <- length(TF$lhs)
    TF$alpha.adj <- alpha / ntests

    out <- subset(TF, TF$LR.pval < TF$alpha.adj)

    # find largest LR
    max.idx <- which(TF$LR == max(TF$LR))

    extra <- list(LR=unname(TF$LR[max.idx]), df=unname(TF$df[max.idx]), 
                  lhs=TF$lhs[max.idx],
                  rhs=TF$rhs[max.idx],
                  group=TF$group[max.idx],
                  p.value=unname(TF$LR.pval[max.idx]),
                  ntests=ntests,
                  p.value.Bonferroni=unname(TF$LR.pval[max.idx]*length(TF$lhs)))

    attr(out, "CpMax") <- extra

    class(out) <- c("lavaan.tables.fit.Cp", "lavaan.data.frame", "data.frame")

    out
}

print.lavaan.tables.fit.Cp <- function(x, ...) {
   cat("Hier komt de header:\n")
   tmp <- x
   class(tmp) <- c("lavaan.data.frame", "data.frame")
   print(tmp)
   cat("blah blah blah\n")
}

# Mariska Barendse CF statistic
lav_tables_fit_CF <- function(object, est = "h0") {

    # check object class
    if(!class(object) %in% c("lavaan")) {
        stop("lavaan ERROR: object must be an object of class lavaan")
    }
    ngroups <- length( object@Data@X )
    CF.group <- rep(as.numeric(NA), ngroups)
    DF.group <- rep(as.numeric(NA), ngroups)

    # check if all ordered
    if(!all(object@Data@ov$type == "ordered")) {
        CF <- as.numeric(NA)
        attr(CF, "CF.group") <- CF.group
        attr(CF, "DF.group") <- DF.group
        return(CF)
    }

    # ord var in this group
    ov.ord <- unique(unlist(object@pta$vnames$ov.ord))
    ov.idx <- which(ov.ord %in% object@Data@ov$name)
    ov.nlev <- object@Data@ov$nlev[ ov.idx ]

    # h0 or h1?
    if(est == "h0") {
        Sigma.hat <- object@Fit@Sigma.hat
        TH        <- object@Fit@TH
        DF        <- prod(ov.nlev) - object@Fit@npar - 1L
    } else {
        ## FIXME: npar should be extract from 'saturated' fit
        Sigma.hat <- object@SampleStats@cov
        TH        <- object@SampleStats@th
        nvar      <- length(ov.ord)
        npar      <- (sum(ov.nlev - 1) + nvar*(nvar-1)/2)*object@Data@ngroups
        DF        <- prod(ov.nlev) - npar - 1L
    }

    for(g in 1:ngroups) {
        F.group <- estimator.FML(Sigma.hat = Sigma.hat[[g]],
                                      TH        = TH[[g]],
                                      th.idx    = object@Model@th.idx[[g]],
                                      num.idx   = object@Model@num.idx[[g]],
                                      X         = object@Model@X[[g]],
                                      cache     = object@Cache[[g]])

        CF.group[g] <- 2*object@Data@nobs[[g]]*F.group
    }

    # check for negative values
    CF.group[CF.group < 0] <- 0.0

    # global test statistic
    CF <- sum(CF.group)

    attr(CF, "CF.group") <- CF.group
    attr(CF, "DF") <- DF
    attr(CF, "rpat.observed") <- sapply(object@Data@Rp, "[[", "npatterns")
    attr(CF, "rpat.total") <- sapply(object@Data@Rp, "[[", "total.patterns")
    attr(CF, "rpat.empty") <- sapply(object@Data@Rp, "[[", "empty.patterns") 

    class(CF) <- c("lavaan.tables.fit.CF", "numeric")

    CF
}

print.lavaan.tables.fit.CF <- function(x, ...) {
    cat("This is CF:", x, "\n")
}

lav_tables_fit_CM <- function(object) {

   CF.h0 <- lav_tables_fit_CF(object, est = "h0")
   CF.h1 <- lav_tables_fit_CF(object, est = "h1")

   CF.h0.group <- attr(CF.h0, "CF.group")
   CF.h1.group <- attr(CF.h1, "CF.group")
   DF.h0 <- attr(CF.h0, "DF")
   DF.h1 <- attr(CF.h1, "DF")

   attributes(CF.h0) <- NULL
   attributes(CF.h1) <- NULL

   CM <- CF.h0 - CF.h1
   attr(CM, "CM.group") <- CF.h0.group - CF.h1.group
   attr(CM, "DF") <- DF.h0 - DF.h1

   CM
}
