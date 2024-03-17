# compute 'Omega' == A1^{-1} B1 A1^{-1}
# where A1 is the expected/observed information matrix of the unrestricted (h1)
# model, and B1 is the first-order information matrix of the unrestricted (h1)
# model
#
# but the exact result will depend on the options:
# for 'A':
# - omega.information ("expected" or "observed")
# - omega.h1.information ("structured" or "unstructured")
# for 'B':
# - omega.information.meat ("first-order")
# - omega.h1.information.meat ("structured" or "unstructured")
#
# special case: if data is complete, A is expected/unstructured, and B is
# unstructured, we get (sample-based) 'Gamma'
#
# YR 28 Oct 2020

lav_model_h1_omega <- function(lavobject = NULL,
                               lavmodel = NULL,
                               lavsamplestats = NULL,
                               lavdata = NULL,
                               lavimplied = NULL,
                               lavh1 = NULL,
                               lavcache = NULL,
                               lavoptions = NULL) {
  if (!is.null(lavobject) && inherits(lavobject, "lavaan")) {
    lavmodel <- lavobject@Model
    lavsamplestats <- lavobject@SampleStats
    lavdata <- lavobject@Data
    lavimplied <- lavobject@implied
    if (.hasSlot(lavobject, "h1")) {
      lavh1 <- lavobject@h1
    } else {
      lavh1 <- lav_h1_implied_logl(
        lavdata = lavobject@Data,
        lavsamplestats = lavobject@SampleStats,
        lavoptions = lavobject@Options
      )
    }
    lavcache <- lavobject@Cache
    lavoptions <- lavobject@Options
  }

  # sanity check
  if (length(lavh1) == 0L) {
    lavh1 <- lav_h1_implied_logl(
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      lavoptions = lavoptions
    )
  }
  if (length(lavimplied) == 0L) {
    lavimplied <- lav_model_implied(lavmodel = lavmodel)
  }

  # set options for A
  A1.options <- lavoptions
  A1.options$information <- lavoptions$omega.information
  A1.options$h1.information <- lavoptions$omega.h1.information

  B1.options <- lavoptions
  B1.options$information <- lavoptions$omega.information.meat # unused
  B1.options$h1.information <- lavoptions$omega.h1.information.meat

  # information
  information <- lavoptions$omega.information

  # compute A1 (per group)
  if (information == "observed") {
    A1 <- lav_model_h1_information_observed(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = A1.options
    )
  } else if (information == "expected") {
    A1 <- lav_model_h1_information_expected(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = A1.options
    )
  } else if (information == "first.order") { # not needed?
    A1 <- lav_model_h1_information_firstorder(
      lavmodel = lavmodel,
      lavsamplestats = lavsamplestats, lavdata = lavdata,
      lavimplied = lavimplied, lavh1 = lavh1,
      lavcache = lavcache, lavoptions = A1.options
    )
  }

  # compute B1 (per group)
  B1 <- lav_model_h1_information_firstorder(
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats, lavdata = lavdata,
    lavimplied = lavimplied, lavh1 = lavh1,
    lavcache = lavcache, lavoptions = B1.options
  )

  # return Omega per group
  Omega <- vector("list", length = lavdata@ngroups)
  trace.h1 <- numeric(lavdata@ngroups)
  h1.ndat <- numeric(lavdata@ngroups)
  for (g in seq_len(lavdata@ngroups)) {
    A1.g <- A1[[g]]
    B1.g <- B1[[g]]

    # mask independent 'fixed-x' variables
    zero.idx <- which(diag(A1.g) == 0)
    if (length(zero.idx) > 0L) {
      A1.inv <- matrix(0, nrow(A1.g), ncol(A1.g))
      a1 <- A1.g[-zero.idx, -zero.idx, drop = FALSE]
      a1.inv <- solve(a1)
      A1.inv[-zero.idx, -zero.idx] <- a1.inv
    } else {
      A1.inv <- solve(A1.g)
    }
    trace.h1[g] <- sum(B1.g * t(A1.inv))
    h1.ndat[g] <- ncol(A1.g) - length(zero.idx)

    Omega[[g]] <- A1.inv %*% B1.g %*% A1.inv
  }

  # store trace.h1 as an attribute (to be used in yuan-bentler)
  attr(Omega, "trace.h1") <- trace.h1
  attr(Omega, "h1.ndat") <- h1.ndat
  attr(Omega, "A.information") <- paste(A1.options$information,
    A1.options$h1.information,
    sep = "."
  )
  attr(Omega, "B.information") <- paste(B1.options$information,
    B1.options$h1.information,
    sep = "."
  )

  Omega
}
