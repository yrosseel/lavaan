# read in information from Mplus difftest output, return as list
#
# line 1: test statistic (unscaled)
# line 2: number of groups
# line 3: number of sample statistics (ndat)
# line 4: number of free parameters (npar)
# delta (ndat x npar)
# P1 (E.inv) lav_matrix_vechr(npar x npar)
# V1 (NVarCov) lav_matrix_vechr(npar x npar)
lavutils_mplus_readdifftest <- function(file = "deriv.dat") {
  ### FIXME: does not work for multiple groups yet!!!

  raw <- scan(file, quiet = TRUE)
  T1 <- raw[1] # function value (usually T1 * 2 * nobs to get X2)
  ngroups <- as.integer(raw[2])
  ndat <- as.integer(raw[3])
  npar <- as.integer(raw[4])
  pstar <- npar * (npar + 1) / 2

  # delta
  offset <- 4L
  delta_raw <- raw[offset + seq_len(npar * ndat)]
  Delta <- matrix(delta_raw, nrow = ndat, ncol = npar, byrow = TRUE)

  # P1
  offset <- 4L + npar * ndat
  p1_raw <- raw[offset + seq_len(pstar)]
  P1 <- lav_matrix_lower2full(p1_raw)

  # (robust) NACOV npar
  offset <- 4L + npar * ndat + pstar
  nacov_raw <- raw[offset + seq_len(pstar)]
  V1 <- lav_matrix_lower2full(nacov_raw)

  # just for fun, M1
  # M1 <- (P1 - P1 %*% H %*% solve(t(H) %*% P1 %*% H) %*% t(H) %*% P1) %*% V1

  list(
    T1 = T1, ngroups = ngroups, ndat = ndat, npar = npar, pstar = pstar,
    Delta = Delta, P1 = P1, V1 = V1
  )
}
