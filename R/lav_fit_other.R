# various fit measures

# - lav_fit_cn
# - lav_fit_wrmr
# - lav_fit_mfi
# - lav_fit_ecvi

# Y.R. 21 July 2022

# Hoelter Critical N (CN)
lav_fit_cn <- function(X2 = NULL, df = NULL, N = NULL, alpha = 0.05) {
  # catch df=0, X2=0
  if (df == 0 && X2 < .Machine$double.eps) {
    CN <- as.numeric(NA)
  } else {
    CN <- qchisq(p = (1 - alpha), df = df) / (X2 / N) + 1
  }

  CN
}

# WRMR
# we use the definition: wrmr = sqrt ( 2*N*F / nel )
# Note: when multiple groups, 'nel' only seems to correspond to the
# first group???
lav_fit_wrmr <- function(X2 = NULL, nel = NULL) {
  if (nel > 0) {
    WRMR <- sqrt(X2 / nel)
  } else {
    WRMR <- as.numeric(NA)
  }

  WRMR
}

# MFI - McDonald Fit Index (McDonald, 1989)
lav_fit_mfi <- function(X2 = NULL, df = NULL, N = NULL) {
  MFI <- exp(-0.5 * (X2 - df) / N)
  MFI
}

# ECVI - cross-validation index (Brown & Cudeck, 1989, eq 5)
# "In the special case where F = F_ML, Equation 5 [=ECVI] is the
# rescaled AIC employed by Cudeck and Browne (1983, Equation 5.1). This
# result is concordant with a finding of Stone (1977). He showed under general
# conditions that if the "leaving one out at a time" method of cross-validation
# (Stone, 1974; Geisser, 1975) is employed, a log-likelihood measure of
# predictive validity is asymptotically equivalent to the AIC." (p. 448)

# not defined for multiple groups and/or models with meanstructures
# TDJ: According to Dudgeon (2004, p. 317), "ECVI requires no adjustment
#      when a model is fitted simultaneously in multiple samples."
#      And I think the lack of mean structure in Brown & Cudeck (1989)
#      was a matter of habitual simplification back then, not necessity.
# YR: - why does Dudgeon eq 22 use (df + 2*npar) instead of (2*npar)??
lav_fit_ecvi <- function(X2 = NULL, npar = npar, N = N) {
  ECVI <- X2 / N + (2 * npar) / N
  ECVI
}
