# utility functions for the sam() function
# YR 4 April 2023

# construct 'mapping matrix' M using either "ML", "GLS" or "ULS" method
# optionally return MTM (for ML)
#
# by construction, M %*% LAMBDA = I (the identity matrix)
lav_sam_mapping_matrix <- function(LAMBDA = NULL, THETA = NULL,
                                   S = NULL, S.inv = NULL,
                                   method = "ML", warn = TRUE) {

  method <- toupper(method)

  # ULS
  # M == solve( t(LAMBDA) %*% LAMBDA ) %*% t(LAMBDA)
  #   == MASS:::ginv(LAMBDA)
  if (method == "ULS") {
    # M == solve( t(LAMBDA) %*% LAMBDA ) %*% t(LAMBDA)
    #   == MASS:::ginv(LAMBDA)
    M <- try(tcrossprod(solve(crossprod(LAMBDA)), LAMBDA),
      silent = TRUE
    )
    if (inherits(M, "try-error")) {
      if (warn) {
        lav_msg_warn(gettext(
          "cannot invert crossprod(LAMBDA); using generalized inverse"))
      }
      M <- MASS::ginv(LAMBDA)
    }


    # GLS
    # M == solve( t(LAMBDA) %*% S.inv %*% LAMBDA ) %*% t(LAMBDA) %*% S.inv
  } else if (method == "GLS") {
    if (is.null(S.inv)) {
      S.inv <- try(solve(S), silent = TRUE)
    }
    if (inherits(S.inv, "try-error")) {
      if (warn) {
        lav_msg_warn(gettext("S is not invertible; switching to ULS method"))
      }
      M <- lav_sam_mapping_matrix(LAMBDA = LAMBDA, method = "ULS")
    } else {
      tLSinv <- t(LAMBDA) %*% S.inv
      tLSinvL <- tLSinv %*% LAMBDA
      M <- try(solve(tLSinvL, tLSinv), silent = TRUE)
      if (inherits(M, "try-error")) {
        if (warn) {
          lav_msg_warn(gettext("problem contructing mapping matrix;
                               switching to generalized inverse"))
        }
        M <- MASS::ginv(tLSinvL) %*% tLSinv
      }
    }

    # ML
    # M == solve(t(LAMBDA) %*% THETA.inv %*% LAMBDA) %*% t(LAMBDA) %*% THETA.inv
  } else if (method == "ML") {
    # Problem: if THETA has zero elements on the diagonal, we cannot invert

    # As we do not have access to Sigma(.inv), we cannot
    # use the trick as in lavPredict() (where we replace THETA.inv by
    # Sigma.inv)

    # old method (<0.6-16): remove rows/cols with zero values
    # on the diagonal of THETA, invert the submatrix and
    # set the '0' diagonal elements to one. This resulted in (somewhat)
    # distorted results.

    # new in 0.6-16: we use the Wall & Amemiya (2000) method using
    # the so-called 'T' transformation

    # new in 0.6-18: if we cannot use the marker method (eg growth models),
    # use the 'old' method anyway: remove zero rows/cols and invert submatrix

    zero.theta.idx <- which(abs(diag(THETA)) < 1e-4) # be conservative
    if (length(zero.theta.idx) == 0L) {
      # ok, no zero diagonal elements: try to invert THETA
      if (lav_matrix_is_diagonal(THETA)) {
        THETA.inv <- diag(1 / diag(THETA), nrow = nrow(THETA))
      } else {
        THETA.inv <- try(solve(THETA), silent = TRUE)
        if (inherits(THETA, "try-error")) {
          THETA.inv <- NULL
        }
      }
    } else {
      # see if we can use marker method
      marker.idx <- lav_utils_get_marker(LAMBDA = LAMBDA, std.lv = TRUE)
      if (any(is.na(marker.idx))) {
        THETA.inv <- try(lav_matrix_symmetric_inverse(THETA), silent = TRUE)
        if (inherits(THETA.inv, "try-error")) {
          THETA.inv <- NULL
        } else {
          diag(THETA.inv)[zero.theta.idx] <- 1
        }
      } else {
        # try tmat method later on
        THETA.inv <- NULL
      }
    }

    # could we invert THETA?
    if (!is.null(THETA.inv)) {
      # ha, all is good; compute M the usual way
      tLTi <- t(LAMBDA) %*% THETA.inv
      tLTiL <- tLTi %*% LAMBDA
      M <- try(solve(tLTiL, tLTi), silent = TRUE)
      if (inherits(M, "try-error")) {
        if (warn) {
          lav_msg_warn(gettext(
            "problem contructing ML mapping matrix; switching to ULS"))
        }
        M <- lav_sam_mapping_matrix(LAMBDA = LAMBDA, method = "ULS")
      }
    } else {
      # use W&A2000's method using the 'T' transformation
      M <- try(lav_sam_mapping_matrix_tmat(
        LAMBDA = LAMBDA,
        THETA = THETA
      ), silent = TRUE)
      if (inherits(M, "try-error")) {
        if (warn) {
          lav_msg_warn(gettext(
            "problem contructing ML mapping matrix; switching to ULS"))
        }
        M <- lav_sam_mapping_matrix(LAMBDA = LAMBDA, method = "ULS")
      }
    }
  } # ML

  M
}

# use 'T' transformation to create the 'Bartlett/ML' mapping matrix
# see Wall & Amemiya (2000) eq (7)
# see also Fuller 1987 page 357 (where T is called H), and page 364

# although Fuller/W&A always assumed that THETA is diagonal,
# their method seems to work equally well for non-diagonal THETA
#
# in our implementation:
# - we do NOT reorder the rows of LAMBDA
# - if std.lv = TRUE, we first rescale to 'create' marker indicators
#   and then rescale back at the end
#
lav_sam_mapping_matrix_tmat <- function(LAMBDA = NULL,
                                        THETA = NULL,
                                        marker.idx = NULL,
                                        std.lv = NULL) {
  LAMBDA <- as.matrix.default(LAMBDA)
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)

  # do we have marker.idx?
  if (is.null(marker.idx)) {
    # 'marker' indicator has a single non-zero element in a row
    marker.idx <- lav_utils_get_marker(LAMBDA = LAMBDA, std.lv = TRUE)
    if (any(is.na(marker.idx))) {
      lav_msg_stop(gettext("no clear markers in LAMBDA matrix"))
    }
  }

  # std.lv TRUE or FALSE?
  if (is.null(std.lv)) {
    std.lv <- FALSE
    if (any(diag(LAMBDA[marker.idx, , drop = FALSE]) != 1)) {
      std.lv <- TRUE
    }
  }

  # if std.lv = TRUE, rescale
  if (std.lv) {
    MARKER <- LAMBDA[marker.idx, , drop = FALSE]
    marker.inv <- 1 / diag(MARKER)
    LAMBDA <- t(t(LAMBDA) * marker.inv)
  }

  # compute 'T' matrix
  TMAT <- lav_sam_tmat(
    LAMBDA = LAMBDA, THETA = THETA,
    marker.idx = marker.idx
  )

  # ML mapping matrix
  M <- TMAT[marker.idx, , drop = FALSE]

  if (std.lv) {
    M <- M * marker.inv
  }

  M
}

# create 'T' matrix (tmat) for T-transformation
#
# Notes: - here we assume that LAMBDA has unity markers (no std.lv = TRUE)
#        - TMAT is NOT symmetric!
#        - Yc %*% t(TMAT) transforms the data in such a way that we get:
#           1)  Bartlett factor scores in the marker columns
#           2) 'V' values in the non-marker columns, where:
#               V = Yc - Yc[,marker.idx] %*% t(LAMBDA)
#
lav_sam_tmat <- function(LAMBDA = NULL,
                         THETA = NULL,
                         marker.idx = NULL) {
  LAMBDA <- as.matrix.default(LAMBDA)
  nvar <- nrow(LAMBDA)
  nfac <- ncol(LAMBDA)

  # do we have marker.idx?
  if (is.null(marker.idx)) {
    # 'marker' indicator has a single 1 element in a row
    marker.idx <- lav_utils_get_marker(LAMBDA = LAMBDA, std.lv = FALSE)
    if (any(is.na(marker.idx))) {
      lav_msg_stop(gettext("no clear markers in LAMBDA matrix"))
    }
  }

  # construct 'C' matrix
  C2 <- diag(nvar)
  C2[, marker.idx] <- -1 * LAMBDA
  C <- C2[-marker.idx, , drop = FALSE]

  # compute Sigma.ve and Sigma.vv
  Sigma.ve <- C %*% THETA
  # Sigma.vv <- C %*% THETA %*% t(C)
  Sigma.vv <- Sigma.ve %*% t(C)

  # construct 'Gamma' (and Gamma2) matrix
  # Gamma <- (t(Sigma.ve) %*% solve(Sigma.vv))[marker.idx,, drop = FALSE]
  Gamma <- try(t(solve(Sigma.vv, Sigma.ve)[, marker.idx, drop = FALSE]),
    silent = TRUE
  )
  if (inherits(Gamma, "try-error")) {
    tmp <- t(Sigma.ve) %*% MASS::ginv(Sigma.vv)
    Gamma <- tmp[marker.idx, , drop = FALSE]
  }
  Gamma2 <- matrix(0, nfac, nvar)
  Gamma2[, -marker.idx] <- Gamma
  Gamma2[, marker.idx] <- diag(nfac)

  # transformation matrix 'T' (we call it here 'Tmat')
  Tmat <- matrix(0, nvar, nvar)
  Tmat[-marker.idx, ] <- C
  Tmat[marker.idx, ] <- -Gamma2 %*% C2

  Tmat
}


# compute VETA
# - if alpha.correction == 0     -> same as local SAM (or MOC)
# - if alpha.correction == (N-1) -> same as FSR+Bartlett
lav_sam_veta <- function(M = NULL, S = NULL, THETA = NULL,
                         alpha.correction = 0L, lambda.correction = TRUE,
                         N = 20L, dummy.lv.idx = integer(0L), extra = FALSE) {
  # MSM
  MSM <- M %*% S %*% t(M)

  # MTM
  MTM <- M %*% THETA %*% t(M)

  # new in 0.6-16: make sure MTM is pd
  # (otherwise lav_matrix_symmetric_diff_smallest_root will fail)
  if (length(dummy.lv.idx) > 0L) {
    MTM.nodummy <- MTM[-dummy.lv.idx, -dummy.lv.idx, drop = FALSE]
    MTM.nodummy <- zapsmall(lav_matrix_symmetric_force_pd(
      MTM.nodummy,
      tol = 1e-04
    ))
    MTM[-dummy.lv.idx, -dummy.lv.idx] <- MTM.nodummy
  } else {
    MTM <- zapsmall(lav_matrix_symmetric_force_pd(MTM, tol = 1e-04))
  }

  # apply small sample correction (if requested)
  if (alpha.correction > 0) {
    alpha.N1 <- alpha.correction / (N - 1)
    if (alpha.N1 > 1.0) {
      alpha.N1 <- 1.0
    } else if (alpha.N1 < 0.0) {
      alpha.N1 <- 0.0
    }
    MTM <- (1 - alpha.N1) * MTM
    alpha <- alpha.correction
  } else {
    alpha <- alpha.correction
  }

  if (lambda.correction) {
    # use Fuller (1987) approach to ensure VETA is positive
    lambda <- try(lav_matrix_symmetric_diff_smallest_root(MSM, MTM),
      silent = TRUE
    )
    if (inherits(lambda, "try-error")) {
      lav_msg_warn(gettext("failed to compute lambda"))
      VETA <- MSM - MTM # and hope for the best
    } else {
      cutoff <- 1 + 1 / (N - 1)
      if (lambda < cutoff) {
        lambda.star <- lambda - 1 / (N - 1)
        VETA <- MSM - lambda.star * MTM
      } else {
        VETA <- MSM - MTM
      }
    }
  } else {
    VETA <- MSM - MTM
  }

  # extra attributes?
  if (extra) {
    attr(VETA, "lambda") <- lambda
    attr(VETA, "alpha") <- alpha
  }

  VETA
}

# compute EETA = E(Eta) = M %*% [YBAR - NU]
lav_sam_eeta <- function(M = NULL, YBAR = NULL, NU = NULL) {
  EETA <- M %*% (YBAR - NU)
  EETA
}

# compute veta including quadratic/interaction terms
lav_sam_veta2 <- function(FS = NULL, M = NULL,
                          VETA = NULL, EETA = NULL, THETA = NULL,
                          lv.names = NULL,
                          lv.int.names = NULL,
                          dummy.lv.names = character(0L),
                          alpha.correction = 0L,
                          lambda.correction = TRUE,
                          extra = FALSE) {

  # small utility function: var() divided by N
  varn <- function(x, N) {
    var(x, use = "pairwise.complete.obs") * (N - 1) / N
  }

  if (length(lv.int.names) == 0L) {
    lav_msg_stop(gettext("lv.int.names is empty: no lv quadratic/interaction
                         terms are provided"))
  }

  if (is.null(lv.names)) {
    lv.names <- paste("eta", seq_len(ncol(FS)), sep = "")
  }

  # MTM
  MTM <- M %*% THETA %*% t(M)

  # new in 0.6-16: make sure MTM is pd
  # (otherwise lav_matrix_symmetric_diff_smallest_root will fail)
  dummy.lv.idx <- which(lv.names %in% dummy.lv.names)
  if (length(dummy.lv.idx) > 0L) {
    MTM.nodummy <- MTM[-dummy.lv.idx, -dummy.lv.idx, drop = FALSE]
    MTM.nodummy <- zapsmall(lav_matrix_symmetric_force_pd(
      MTM.nodummy,
      tol = 1e-04
    ))
    MTM[-dummy.lv.idx, -dummy.lv.idx] <- MTM.nodummy
  } else {
    MTM <- zapsmall(lav_matrix_symmetric_force_pd(MTM, tol = 1e-04))
  }

  # augment to include intercept
  FS <- cbind(1, FS)
  N <- nrow(FS)
  MTM <- lav_matrix_bdiag(0, MTM)
  VETA <- lav_matrix_bdiag(0, VETA)
  EETA <- c(1, EETA)
  lv.names <- c("..int..", lv.names)
  nfac <- ncol(FS)

  idx1 <- rep(seq_len(nfac), each = nfac)
  idx2 <- rep(seq_len(nfac), times = nfac)

  NAMES <- paste(lv.names[idx1], lv.names[idx2], sep = ":")
  NAMES[seq_len(nfac)] <- lv.names

  FS2 <- FS[, idx1] * FS[, idx2]

  K.nfac <- lav_matrix_commutation(nfac, nfac)
  IK <- diag(nfac * nfac) + K.nfac

  EETA <- as.matrix(drop(EETA))
  VETAkMTM <- VETA %x% MTM

  # normal version (for now):
  Gamma.ME22 <- IK %*% (MTM %x% MTM)

  # ingredients (normal ME case)
  Var.FS2 <- varn(FS2, N)
  Var.ETAkME <- (tcrossprod(EETA) %x% MTM + VETAkMTM)
  Var.MEkETA <- lav_matrix_commutation_pre_post(Var.ETAkME)
  Var.ME2 <- Gamma.ME22

  cov.ETAkME.MEkETA <- lav_matrix_commutation_post(Var.ETAkME)
  cov.MEkETA.ETAkME <- t(cov.ETAkME.MEkETA)

  Var.ERROR <- (Var.ETAkME + Var.MEkETA + cov.ETAkME.MEkETA
    + cov.MEkETA.ETAkME + Var.ME2)

  # select only what we need
  colnames(Var.FS2) <- rownames(Var.FS2) <- NAMES
  colnames(Var.ERROR) <- rownames(Var.ERROR) <- NAMES
  lv.keep <- c(lv.names[-1], lv.int.names)
  Var.FS2 <- Var.FS2[lv.keep, lv.keep]
  Var.ERROR <- Var.ERROR[lv.keep, lv.keep]

  # apply small sample correction (if requested)
  if (alpha.correction > 0) {
    alpha.N1 <- alpha.correction / (N - 1)
    if (alpha.N1 > 1.0) {
      alpha.N1 <- 1.0
    } else if (alpha.N1 < 0.0) {
      alpha.N1 <- 0.0
    }
    Var.ERROR <- (1 - alpha.N1) * Var.ERROR
    alpha <- alpha.correction
  } else {
    alpha <- alpha.correction
  }

  if (lambda.correction) {
    # use Fuller (1987) approach to ensure VETA2 is positive
    lambda <- try(lav_matrix_symmetric_diff_smallest_root(
      Var.FS2,
      Var.ERROR
    ), silent = TRUE)
    if (inherits(lambda, "try-error")) {
      lav_msg_warn(gettext("failed to compute lambda"))
      VETA2 <- Var.FS2 - Var.ERROR # and hope for the best
    } else {
      cutoff <- 1 + 1 / (N - 1)
      if (lambda < cutoff) {
        lambda.star <- lambda - 1 / (N - 1)
        VETA2 <- Var.FS2 - lambda.star * Var.ERROR
      } else {
        VETA2 <- Var.FS2 - Var.ERROR
      }
    }
  } else {
    VETA2 <- Var.FS2 - Var.ERROR
  }

  # extra attributes?
  if (extra) {
    attr(VETA2, "lambda") <- lambda
    attr(VETA2, "alpha") <- alpha
  }

  VETA2
}

lav_sam_eeta2 <- function(EETA = NULL, VETA = NULL, lv.names = NULL,
                          lv.int.names = NULL) {
  if (length(lv.int.names) == 0L) {
    lav_msg_stop(gettext("lv.int.names is empty: no lv quadratic/interaction
                         terms are provided"))
  }

  if (is.null(lv.names)) {
    lv.names <- paste("eta", seq_len(ncol(VETA)), sep = "")
  }

  nfac <- nrow(VETA)
  idx1 <- rep(seq_len(nfac), each = nfac)
  idx2 <- rep(seq_len(nfac), times = nfac)
  NAMES <- c(lv.names, paste(lv.names[idx1], lv.names[idx2], sep = ":"))

  # E(\eta %x% \eta)
  EETA2 <- lav_matrix_vec(VETA) + EETA %x% EETA

  # add 1st order
  EETA2.aug <- c(EETA, EETA2)

  # select only what we need
  names(EETA2.aug) <- NAMES
  lv.keep <- c(lv.names, lv.int.names)
  EETA2.aug <- EETA2.aug[lv.keep]

  EETA2.aug
}

# compute veta including quadratic/interaction terms
lav_sam_fs2 <- function(FS = NULL, lv.names = NULL, lv.int.names = NULL) {
  varn <- function(x, N) {
    var(x) * (N - 1) / N
  }

  if (length(lv.int.names) == 0L) {
    lav_msg_stop(gettext("lv.int.names is empty: no lv quadratic/interaction
                         terms are provided"))
  }

  if (is.null(lv.names)) {
    lv.names <- paste("eta", seq_len(ncol(FS)), sep = "")
  }

  # augment to include intercept
  FS <- cbind(1, FS)
  N <- nrow(FS)
  lv.names <- c("..int..", lv.names)
  nfac <- ncol(FS)

  idx1 <- rep(seq_len(nfac), each = nfac)
  idx2 <- rep(seq_len(nfac), times = nfac)

  NAMES <- paste(lv.names[idx1], lv.names[idx2], sep = ":")

  FS2 <- FS[, idx1] * FS[, idx2]
  Var.FS2 <- varn(FS2, N)

  # select only what we need
  colnames(Var.FS2) <- rownames(Var.FS2) <- NAMES
  lv.main <- paste(lv.names[-1], "..int..", sep = ":")
  lv.keep <- c(lv.main, lv.int.names)
  Var.FS2 <- Var.FS2[lv.keep, lv.keep]

  Var.FS2
}

# create consistent lavaan object, based on (filled in) PT
lav_sam_step3_joint <- function(FIT = NULL, PT = NULL, sam.method = "local") {
  lavoptions <- FIT@Options

  lavoptions.joint <- lavoptions
  lavoptions.joint$optim.method <- "none"
  lavoptions.joint$optim.force.converged <- TRUE
  lavoptions.joint$check.gradient <- FALSE
  lavoptions.joint$check.start <- FALSE
  lavoptions.joint$check.post <- FALSE
  lavoptions.joint$rotation <- "none"
  lavoptions.joint$se <- "none"
  lavoptions.joint$store.vcov <- FALSE # we do this manually
  lavoptions.joint$verbose <- FALSE

  if (sam.method %in% c("local", "fsr")) {
    lavoptions.joint$baseline <- FALSE
    lavoptions.joint$sample.icov <- FALSE
    lavoptions.joint$h1 <- FALSE
    lavoptions.joint$test <- "none"
    lavoptions.joint$estimator <- "none"
  } else {
    lavoptions.joint$test <- lavoptions$test
    lavoptions.joint$estimator <- lavoptions$estimator
  }

  # set ustart values
  PT$ustart <- PT$est # as this is used if optim.method == "none"

  JOINT <- lavaan::lavaan(PT,
    slotOptions = lavoptions.joint,
    slotSampleStats = FIT@SampleStats,
    slotData = FIT@Data
  )
  JOINT
}

lav_sam_table <- function(JOINT = NULL, STEP1 = NULL, FIT.PA = FIT.PA,
                          mm.args = list(), struc.args = list(),
                          sam.method = "local",
                          local.options = list(), global.options = list()) {
  MM.FIT <- STEP1$MM.FIT

  sam.mm.table <- data.frame(
    Block = seq_len(length(STEP1$mm.list)),
    Latent = sapply(MM.FIT, function(x) {
      paste(unique(unlist(x@pta$vnames$lv)), collapse = ",")
    }),
    Nind = sapply(MM.FIT, function(x) {
      length(unique(unlist(x@pta$vnames$ov)))
    }),
    # Estimator = sapply(MM.FIT, function(x) { x@Model@estimator} ),
    Chisq = sapply(MM.FIT, function(x) {
      x@test[[1]]$stat
    }),
    Df = sapply(MM.FIT, function(x) {
      x@test[[1]]$df
    })
  )
  # pvalue = sapply(MM.FIT, function(x) {x@test[[1]]$pvalue}) )
  class(sam.mm.table) <- c("lavaan.data.frame", "data.frame")


  # extra info for @internal slot
  if (sam.method %in% c("local", "fsr")) {
    sam.struc.fit <- try(
      fitMeasures(
        FIT.PA,
        c(
          "chisq", "df", # "pvalue",
          "cfi", "rmsea", "srmr"
        )
      ),
      silent = TRUE
    )
    if (inherits(sam.struc.fit, "try-error")) {
      sam.struc.fit <- "(unable to obtain fit measures)"
      names(sam.struc.fit) <- "warning"
    }
    sam.mm.rel <- STEP1$REL
  } else {
    sam.struc.fit <- "no local fit measures available for structural part if sam.method is global"
    names(sam.struc.fit) <- "warning"
    sam.mm.rel <- numeric(0L)
  }


  SAM <- list(
    sam.method = sam.method,
    sam.local.options = local.options,
    sam.global.options = global.options,
    sam.mm.list = STEP1$mm.list,
    sam.mm.estimator = MM.FIT[[1]]@Model@estimator,
    sam.mm.args = mm.args,
    sam.mm.ov.names = lapply(MM.FIT, function(x) {
      x@pta$vnames$ov
    }),
    sam.mm.table = sam.mm.table,
    sam.mm.rel = sam.mm.rel,
    sam.struc.estimator = FIT.PA@Model@estimator,
    sam.struc.args = struc.args,
    sam.struc.fit = sam.struc.fit
  )
  SAM
}
