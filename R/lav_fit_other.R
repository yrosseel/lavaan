# various fit measures

# - lav_fit_cn
# - lav_fit_wrmr
# - lav_fit_mfi
# - lav_fit_ecvi

# Y.R. 21 July 2022

# Hoelter Critical N (CN)
lav_fit_cn <- function(X2 = NULL, df = NULL, N = NULL, alpha = 0.05) {

    # catch df=0, X2=0
    if(df == 0 && X2 < .Machine$double.eps) {
        CN <- as.numeric(NA)
    } else {
        CN <- qchisq(p = (1 - alpha), df = df)/(X2/N) + 1
    }

    CN
}

# WRMR
# we use the definition: wrmr = sqrt ( 2*N*F / nel )
# Note: when multiple groups, 'nel' only seems to correspond to the
# first group???
lav_fit_wrmr <- function(X2 = NULL, nel = NULL) {

    if(nel > 0) {
        WRMR <- sqrt( X2 / nel )
    } else {
        WRMR <- as.numeric(NA)
    }

    WRMR
}

# MFI - McDonald Fit Index (McDonald, 1989)
lav_fit_mfi <- function(X2 = NULL, df = NULL, N = NULL) {
    MFI <- exp(-0.5 * (X2 - df)/N)
    MFI
}

# ECVI - cross-validation index (Brown & Cudeck, 1989)
# not defined for multiple groups and/or models with meanstructures
# TDJ: According to Dudgeon (2004, p. 317), "ECVI requires no adjustment
#      when a model is fitted simultaneously in multiple samples."
#      And I think the lack of mean structure in Brown & Cudeck (1989)
#      was a matter of habitual simplification back then, not necessity.
lav_fit_ecvi <- function(X2 = NULL, npar = npar, N = N) {
    ECVI <- X2/N + (2 * npar)/N
    ECVI
}


