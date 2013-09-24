# This code is written by YR (using lavaan components), but based on 
# research code written by Mariska Barendse (Groningen/Amsterdam, NL)
#
# September 2013
#
# Three fit indices for the PML estimator (if all categorical, no exo)
# - Cp(max)
# - CF
# - CM

# Mariska Barendse Cp statistic
lav_tables_fit_Cp <- function(object) {

    out <- lavTablesFit(object, stat = "LR", p.value = TRUE)

    # Bonferonni adjusted p-value
    ntests <- length(out$lhs)
    out$pval.adj <- pchisq(out$LR, df=out$df, lower.tail = FALSE) * ntests

    # remove LR.h0.pval
    #out$LR.h0.pval <- NULL

    out
}

lav_tables_fit_CpMax <- function(object) {
    out <- lav_tables_fit_Cp(object = object)

    # find largest LR
    max.idx <- which(out$LR == max(out$LR))

    list(LR=unname(out$LR[max.idx]), df=unname(out$df[max.idx]), 
         p.value=unname(out$LR.pval[max.idx]),
         p.value.Bonferroni=unname(out$LR.pval[max.idx]*length(out$lhs)))
}

# Mariska Barendse CF statistic
lav_tables_fit_CF <- function(object, est = "h0") {

    # check object class
    if(!class(object) %in% c("lavaan")) {
        stop("lavaan ERROR: object must be an object of class lavaan")
    }
    ngroups <- length( object@Data@X )

    CF.group <- numeric(ngroups)
    DF.group <- numeric(ngroups)

    # h0 or h1?
    if(est == "h0") {
        Sigma.hat <- object@Fit@Sigma.hat
        TH        <- object@Fit@TH
    } else {
        Sigma.hat <- object@SampleStats@cov
        TH        <- object@Fit@TH
    }

    for(g in 1:ngroups) {
        logLik.group <- estimator.FML(Sigma.hat = Sigma.hat[[g]],
                                      TH        = TH[[g]],
                                      th.idx    = object@Model@th.idx[[g]],
                                      num.idx   = object@Model@num.idx[[g]],
                                      X         = object@Model@X[[g]],
                                      cache     = object@Cache[[g]])

        freq <- as.numeric( rownames(object@Data@Rp[[g]]$pat) )
        CF.group[g] <- 2*logLik.group + 2*sum(freq*log(freq/sum(freq)))


        # ord var in this group
        ov.ord <- object@pta$vnames$ov.ord[[g]]
        ov.idx <- which(ov.ord %in% object@Data@ov$name)
        ov.nlev <- object@Data@ov$nlev[ ov.idx ]

        DF.group[g] <- prod(ov.nlev) - object@Fit@npar - 1L
    }

    # check for negative values
    CF.group[CF.group < 0] <- 0.0

    # global test statistic
    CF <- sum(CF.group)

    attr(CF, "CF.group") <- CF.group
    attr(CF, "DF.group") <- DF.group

    CF
}

lav_tables_fit_CM <- function(object) {

   CF.h0 <- lav_tables_fit_CF(object, est = "h0")
   CF.h1 <- lav_tables_fit_CF(object, est = "h1")

   CF.h0.group <- attr(CF.h0, "CF.group")
   CF.h1.group <- attr(CF.h1, "CF.group")
   DF.h0.group <- attr(CF.h0, "DF.group")
   DF.h1.group <- attr(CF.h1, "DF.group")

   attributes(CF.h0) <- NULL
   attributes(CF.h1) <- NULL

   CM <- CF.h0 - CF.h1
   attr(CM, "CM.group") <- CF.h0.group - CF.h1.group
   attr(CM, "DF.group") <- DF.h0.group - DF.h1.group

   CM
}
