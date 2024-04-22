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
# lav_tables_fit_Cp <- function(object, alpha = 0.05) {
#
#    out <- lavTablesFit(object, statistic = "G2", p.value = TRUE)
#
#    # Bonferonni adjusted p-value
#    ntests <- length(out$lhs)
#    out$alpha.adj <- alpha / ntests
#    #out$pval <- pchisq(out$G2, df=out$df, lower.tail = FALSE)
#
#    # remove G2.h0.pval
#    #out$G2.h0.pval <- NULL
#
#    out
# }

lavTablesFitCp <- function(object, alpha = 0.05) {
  lavdata <- object@Data

  if (!any(lavdata@ov$type == "ordered")) {
    return(list(
      G2 = as.numeric(NA), df = as.numeric(NA),
      p.value = as.numeric(NA), p.value.Bonferroni = as.numeric(NA)
    ))
  }

  TF <- lavTables(object,
    dimension = 2L, type = "table",
    statistic = "G2", p.value = TRUE
  )

  # Bonferonni adjusted p-value
  ntests <- length(TF$lhs)
  TF$alpha.adj <- alpha / ntests

  out <- subset(TF, TF$G2.pval < TF$alpha.adj)

  # find largest G2
  max.idx <- which(TF$G2 == max(TF$G2))

  extra <- list(
    G2 = unname(TF$G2[max.idx]), df = unname(TF$df[max.idx]),
    lhs = TF$lhs[max.idx],
    rhs = TF$rhs[max.idx],
    group = TF$group[max.idx],
    p.value = unname(TF$G2.pval[max.idx]),
    ntests = ntests,
    p.value.Bonferroni = unname(TF$G2.pval[max.idx] * length(TF$lhs))
  )

  attr(out, "CpMax") <- extra

  class(out) <- c("lavaan.tables.fit.Cp", "lavaan.data.frame", "data.frame")

  out
}

print.lavaan.tables.fit.Cp <- function(x, ...) {
  cat("CP-values that are significant at a Bonferroni adjusted level of significance\n")
  tmp <- x
  class(tmp) <- c("lavaan.data.frame", "data.frame")
  print(tmp)
}

# Mariska Barendse CF statistic
lavTablesFitCf <- function(object) {
  # check object class
  if (!inherits(object, "lavaan")) {
    lav_msg_stop(gettext("object must be an object of class lavaan"))
  }
  lavdata <- object@Data
  lavpta <- object@pta
  lavmodel <- object@Model
  lavcache <- object@Cache
  implied <- object@implied

  CF.group <- rep(as.numeric(NA), lavdata@ngroups)
  DF.group <- rep(as.numeric(NA), lavdata@ngroups)

  # check if all ordered
  if (!any(lavdata@ov$type == "ordered")) {
    CF <- as.numeric(NA)
    attr(CF, "CF.group") <- CF.group
    attr(CF, "DF.group") <- DF.group
    return(CF)
  }

  # ord var in this group
  ov.ord <- unique(unlist(lavpta$vnames$ov.ord))
  ov.idx <- which(ov.ord %in% lavdata@ov$name)
  ov.nlev <- lavdata@ov$nlev[ov.idx]

  Sigma.hat <- if (lavmodel@conditional.x) implied$res.cov else implied$cov
  TH <- if (lavmodel@conditional.x) implied$res.th else implied$th
  DF <- prod(ov.nlev) - object@optim$npar - 1L

  for (g in seq_len(lavdata@ngroups)) {
    F.group <- estimator.FML(
      Sigma.hat = Sigma.hat[[g]],
      TH = TH[[g]],
      th.idx = lavmodel@th.idx[[g]],
      num.idx = lavmodel@num.idx[[g]],
      X = lavdata@X[[g]],
      lavcache = lavcache[[g]]
    )
    CF.group[g] <- 2 * lavdata@nobs[[g]] * F.group
  }

  # check for negative values
  CF.group[CF.group < 0] <- 0.0

  # global test statistic
  CF <- sum(CF.group)

  attr(CF, "CF.group") <- CF.group
  attr(CF, "DF") <- DF
  attr(CF, "rpat.observed") <- sapply(lavdata@Rp, "[[", "npatterns")
  attr(CF, "rpat.total") <- sapply(lavdata@Rp, "[[", "total.patterns")
  attr(CF, "rpat.empty") <- sapply(lavdata@Rp, "[[", "empty.patterns")

  class(CF) <- c("lavaan.tables.fit.Cf", "numeric")

  CF
}

print.lavaan.tables.fit.Cf <- function(x, ...) {
  cat("Total response patterns: ", attr(x, "rpat.total"), "\n")
  cat("Observed response patterns: ", attr(x, "rpat.observed"), "\n")
  cat("Empty response patterns: ", attr(x, "rpat.empty"), "\n")
  cat("Cf results may be biased because of large numbers of empty cells in the multivariate contingency table\n")
  cat("Cf-value, overall:\n")
  CF <- unclass(x)
  attributes(CF) <- NULL
  print(CF)
  CF.group <- attr(x, "CF.group")
  if (length(CF.group) > 1L) {
    cat("Cf-value, per group:\n")
    print(CF.group)
  }
  cat("Degrees of freedom\n")
  print(attr(x, "DF"))
}

lavTablesFitCm <- function(object) {
  lavdata <- object@Data
  lavoptions <- object@Options

  CF.h0 <- lavTablesFitCf(object)

  # fit unrestricted model
  h1 <- lavCor(lavdata,
    estimator = lavoptions$estimator,
    se = "none", test = "none", output = "lavaan"
  )
  CF.h1 <- lavTablesFitCf(h1)

  CF.h0.group <- attr(CF.h0, "CF.group")
  CF.h1.group <- attr(CF.h1, "CF.group")
  DF.h0 <- attr(CF.h0, "DF")
  DF.h1 <- attr(CF.h1, "DF")

  attributes(CF.h0) <- NULL
  attributes(CF.h1) <- NULL

  CM <- CF.h0 - CF.h1
  attr(CM, "CM.group") <- CF.h0.group - CF.h1.group
  attr(CM, "DF") <- DF.h0 - DF.h1

  class(CM) <- c("lavaan.tables.fit.Cm", "numeric")

  CM
}


print.lavaan.tables.fit.Cm <- function(x, ...) {
  # cat("The percentage of empty cells\n") #weet niet goed want FML werkt niet
  # cat("CM results may be a little biased because of large numbers of empty cells in the multivariate contingency table\n")
  cat("Cm-value, overall:\n")
  CM <- unclass(x)
  attributes(CM) <- NULL
  print(CM)
  CM.group <- attr(x, "CM.group")
  if (length(CM.group) > 1L) {
    cat("Cm-value, per group:\n")
    print(CM.group)
  }
  cat("Degrees of freedom:\n")
  print(attr(x, "DF"))
}
