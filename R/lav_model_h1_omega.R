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
  lav_model_h1_unpack(lavobject, need_h1 = TRUE, need_implied = TRUE)

  # set options for A
  a1_options <- lavoptions
  a1_options$information <- lavoptions$omega.information
  a1_options$h1.information <- lavoptions$omega.h1.information

  b1_options <- lavoptions
  b1_options$information <- lavoptions$omega.information.meat # unused
  b1_options$h1.information <- lavoptions$omega.h1.information.meat

  # compute A1 (per group)
  a1_1 <- lav_model_h1_info(
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats, lavdata = lavdata,
    lavimplied = lavimplied, lavh1 = lavh1,
    lavcache = lavcache, lavoptions = a1_options
  )

  # compute B1 (per group)
  b1 <- lav_model_h1_info_firstorder(
    lavmodel = lavmodel,
    lavsamplestats = lavsamplestats, lavdata = lavdata,
    lavimplied = lavimplied, lavh1 = lavh1,
    lavcache = lavcache, lavoptions = b1_options
  )

  # return Omega per group
  omega <- vector("list", length = lavdata@ngroups)
  trace_h1 <- numeric(lavdata@ngroups)
  h1_ndat <- numeric(lavdata@ngroups)
  for (g in seq_len(lavdata@ngroups)) {
    a1_g <- a1_1[[g]]
    b1_g <- b1[[g]]

    # mask independent 'fixed-x' variables
    zero_idx <- which(diag(a1_g) == 0)
    if (length(zero_idx) > 0L) {
      a1_inv_1 <- matrix(0, nrow(a1_g), ncol(a1_g))
      a1 <- a1_g[-zero_idx, -zero_idx, drop = FALSE]
      a1_inv <- solve(a1)
      a1_inv_1[-zero_idx, -zero_idx] <- a1_inv
    } else {
      a1_inv_1 <- solve(a1_g)
    }
    trace_h1[g] <- sum(b1_g * t(a1_inv_1))
    h1_ndat[g] <- ncol(a1_g) - length(zero_idx)

    omega[[g]] <- a1_inv_1 %*% b1_g %*% a1_inv_1
  }

  # store trace.h1 as an attribute (to be used in yuan-bentler)
  attr(omega, "trace.h1") <- trace_h1
  attr(omega, "h1.ndat") <- h1_ndat
  attr(omega, "A.information") <- paste(a1_options$information,
    a1_options$h1.information,
    sep = "."
  )
  attr(omega, "B.information") <- paste(b1_options$information,
    b1_options$h1.information,
    sep = "."
  )

  omega
}
