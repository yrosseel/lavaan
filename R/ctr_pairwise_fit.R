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

lavTablesFitCp <- function(object, alpha = 0.05) {   # nolint
  lavdata <- object@Data

  if (!any(lavdata@ov$type == "ordered")) {
    return(list(
      G2 = as.numeric(NA), df = as.numeric(NA),
      p.value = as.numeric(NA), p.value.Bonferroni = as.numeric(NA)
    ))
  }

  tf <- lavTables(object,
    dimension = 2L, type = "table",
    statistic = "G2", p.value = TRUE
  )

  # Bonferonni adjusted p-value
  ntests <- length(tf$lhs)
  tf$alpha.adj <- alpha / ntests

  out <- subset(tf, tf$G2.pval < tf$alpha.adj)

  # find largest G2
  max_idx <- which(tf$G2 == max(tf$G2))

  extra <- list(
    G2 = unname(tf$G2[max_idx]), df = unname(tf$df[max_idx]),
    lhs = tf$lhs[max_idx],
    rhs = tf$rhs[max_idx],
    group = tf$group[max_idx],
    p.value = unname(tf$G2.pval[max_idx]),
    ntests = ntests,
    p.value.Bonferroni = unname(tf$G2.pval[max_idx] * length(tf$lhs))
  )

  attr(out, "CpMax") <- extra

  class(out) <- c("lavaan.tables.fit.Cp", "lavaan.data.frame", "data.frame")

  out
}

lav_tables_fit_cp_print <- function(x, ...) {
  cat("CP-values that are significant at a Bonferroni",
      "adjusted level of significance\n")
  tmp <- x
  class(tmp) <- c("lavaan.data.frame", "data.frame")
  print(tmp)
}

# Mariska Barendse CF statistic
lavTablesFitCf <- function(object) {    # nolint
  # check object class
  if (!inherits(object, "lavaan")) {
    lav_msg_stop(gettext("object must be an object of class lavaan"))
  }
  lavdata <- object@Data
  lavpta <- object@pta
  lavmodel <- object@Model
  lavcache <- object@Cache
  implied <- object@implied

  cf_group <- rep(as.numeric(NA), lavdata@ngroups)
  df_group <- rep(as.numeric(NA), lavdata@ngroups)

  # check if all ordered
  if (!any(lavdata@ov$type == "ordered")) {
    cf <- as.numeric(NA)
    attr(cf, "CF.group") <- cf_group
    attr(cf, "DF.group") <- df_group
    return(cf)
  }

  # ord var in this group
  ov_ord <- unique(unlist(lavpta$vnames$ov.ord))
  ov_idx <- which(ov_ord %in% lavdata@ov$name)
  ov_nlev <- lavdata@ov$nlev[ov_idx]

  sigma_hat <- if (lavmodel@conditional.x) implied$res.cov else implied$cov
  th <- if (lavmodel@conditional.x) implied$res.th else implied$th
  df_1 <- prod(ov_nlev) - object@optim$npar - 1L

  for (g in seq_len(lavdata@ngroups)) {
    f_group <- lav_model_objective_fml(
      Sigma.hat = sigma_hat[[g]],
      TH = th[[g]],
      th.idx = lavmodel@th.idx[[g]],
      num.idx = lavmodel@num.idx[[g]],
      X = lavdata@X[[g]],
      lavcache = lavcache[[g]]
    )
    cf_group[g] <- 2 * lavdata@nobs[[g]] * f_group
  }

  # check for negative values
  cf_group[cf_group < 0] <- 0.0

  # global test statistic
  cf <- sum(cf_group)

  attr(cf, "CF.group") <- cf_group
  attr(cf, "DF") <- df_1
  attr(cf, "rpat.observed") <- sapply(lavdata@Rp, "[[", "npatterns")
  attr(cf, "rpat.total") <- sapply(lavdata@Rp, "[[", "total.patterns")
  attr(cf, "rpat.empty") <- sapply(lavdata@Rp, "[[", "empty.patterns")

  class(cf) <- c("lavaan.tables.fit.Cf", "numeric")

  cf
}

lav_tables_fit_cf_print <- function(x, ...) {
  cat("Total response patterns: ", attr(x, "rpat.total"), "\n")
  cat("Observed response patterns: ", attr(x, "rpat.observed"), "\n")
  cat("Empty response patterns: ", attr(x, "rpat.empty"), "\n")
  cat("Cf results may be biased because of large numbers of empty",
      "cells in the multivariate contingency table\n")
  cat("Cf-value, overall:\n")
  cf <- unclass(x)
  attributes(cf) <- NULL
  print(cf)
  cf_group <- attr(x, "CF.group")
  if (length(cf_group) > 1L) {
    cat("Cf-value, per group:\n")
    print(cf_group)
  }
  cat("Degrees of freedom\n")
  print(attr(x, "DF"))
}

lavTablesFitCm <- function(object) {  # nolint
  lavdata <- object@Data
  lavoptions <- object@Options

  cf_h0 <- lavTablesFitCf(object)

  # fit unrestricted model
  h1 <- lav_object_cor(lavdata,
    estimator = lavoptions$estimator,
    se = "none", test = "none", output = "lavaan"
  )
  cf_h1 <- lavTablesFitCf(h1)

  cf_h0_group <- attr(cf_h0, "CF.group")
  cf_h1_group <- attr(cf_h1, "CF.group")
  df_h0 <- attr(cf_h0, "DF")
  df_h1 <- attr(cf_h1, "DF")

  attributes(cf_h0) <- NULL
  attributes(cf_h1) <- NULL

  cm_1 <- cf_h0 - cf_h1
  attr(cm_1, "CM.group") <- cf_h0_group - cf_h1_group
  attr(cm_1, "DF") <- df_h0 - df_h1

  class(cm_1) <- c("lavaan.tables.fit.Cm", "numeric")

  cm_1
}


lav_tables_fit_cm_print <- function(x, ...) {
  # cat("The percentage of empty cells\n") #weet niet goed want FML werkt niet
  # cat("CM results may be a little biased because of large numbers of empty ",
  #      "cells in the multivariate contingency table\n")
  cat("Cm-value, overall:\n")
  cm_1 <- unclass(x)
  attributes(cm_1) <- NULL
  print(cm_1)
  cm_group <- attr(x, "CM.group")
  if (length(cm_group) > 1L) {
    cat("Cm-value, per group:\n")
    print(cm_group)
  }
  cat("Degrees of freedom:\n")
  print(attr(x, "DF"))
}
