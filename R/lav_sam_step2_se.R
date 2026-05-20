# compute two-step standard errors for SAM models
#
# several possibilities:
# 1) se = "twostep": classic (but global) two-step corrected SEs
#     - create 'global' model, only to get the 'joint' information matrix
#     - partition information matrix (step 1, step 2)
#     - apply two-step correction for second step
#     - 'insert' these corrected SEs (and vcov) in JOINT
# 2) se = "standard": using I.22.inv (but without correction term)
# 3) se = "naive": grab (naive) VCOV from FIT.PA
# 4) se = "local": grab (robust) VCOV from FIT.PA

lav_sam_step2_se <- function(fit = NULL, joint = NULL,
                             step1 = NULL, step2 = NULL,
                             local_options = list()) {
  out <- list()
  sigma_11 <- step1$Sigma.11
  step1_free_idx <- step1$step1.free.idx
  step2_free_idx <- step2$step2.free.idx
  lavoptions <- fit@Options
  nlevels <- fit@pta$nlevels
  fit_pa <- step2$FIT.PA
  extra_id <- step2$extra.id

  # catch empty step2.free.idx
  if (length(step2_free_idx) == 0L) {
    # no (free) structural parameters at all!
    out <- list(
      V1 = matrix(0, 0, 0), V2 = matrix(0, 0, 0),
      VCOV = matrix(0, 0, 0)
    )
    return(out)
  }

  if (!lavoptions$se %in%
    c("none", "standard", "naive", "twostep", "twostep.robust",
      "local", "local.nt")) {
    lav_msg_warn(gettextf(
      "unknown se= argument: %s. Switching to twostep.",
      lavoptions$se
    ))
  }

  if (lavoptions$se == "none") {
    return(out)
  }

  if (lav_verbose()) {
    cat("Computing ", lavoptions$se, " standard errors ... ", sep = "")
  }

  if (lavoptions$se %in% c("naive", "twostep", "twostep.robust")) {
    info <- lavInspect(joint, "information")
    i_12 <- info[step1_free_idx, step2_free_idx]
    i_22 <- info[step2_free_idx, step2_free_idx]
    i_21 <- info[step2_free_idx, step1_free_idx]
  }

  # V2
  if (nlevels > 1L) {
    # FIXME: not ok for multigroup multilevel
    n <- fit@Data@Lp[[1]]$nclusters[[2]] # first group only
  } else {
    n <- nobs(fit)
  }

  # total number of free parameters STRUC
  if (fit_pa@Model@ceq.simple.only) {
    npar <- fit_pa@Model@nx.unco
    pts_free <- fit_pa@ParTable$free
    pts_free[pts_free > 0] <- seq_len(npar)
  } else {
    npar <- fit_pa@Model@nx.free
    pts_free <- fit_pa@ParTable$free
  }

  # do we have 'extra' free parameters in FIT.PA that are not free in JOINT?
  step2_rm_idx <- integer(0L)
  if (length(extra_id) > 0L) {
    id_idx <- which(fit_pa@ParTable$id %in% extra_id &
      fit_pa@ParTable$free > 0L)
    step2_rm_idx <- pts_free[id_idx]
  }

  # Fix for EFA/ESEM: when rotation is used, FIT.PA@Model@con.jac includes
  # columns for rotation identification constraints that are not part of
  # step2.free.idx. These extra columns cause dimension mismatch in
  # lav_model_information_augment_invert(). Remove them via rm.idx.
  if (nrow(fit_pa@Model@con.jac) > 0L) {
    n_jac_cols <- ncol(fit_pa@Model@con.jac)
    n_step2 <- length(step2_free_idx)
    if (n_jac_cols > n_step2) {
      step2_rm_idx <- union(step2_rm_idx, (n_step2 + 1):n_jac_cols)
    }
  }

  # invert augmented information, for I.22 block only
  # new in 0.6-16 (otherwise, eq constraints in struc part are ignored)
  if (lavoptions$se %in% c("standard", "twostep", "twostep.robust")) {
    i_22_inv <-
      lav_model_information_augment_invert(
        lavmodel = fit_pa@Model,
        information = i_22,
        inverted = TRUE,
        use_ginv = FALSE, # if interaction, SEs end up smaller than naive...
        rm_idx = step2_rm_idx
      )
    if (inherits(i_22_inv, "try-error")) {
      # hm, not good
      if (lavoptions$se != "naive") {
        lav_msg_warn(gettext(
          "problem inverting information matrix (I.22); -> switching
            to naive standard errors!"
        ))
        lavoptions$se <- "naive"
      }
    }
  } # se needs I.22.inv

  # method below has the advantage that we can use a 'robust' vcov
  # for the joint model;
  # but does not work if we have equality constraints in the MM!
  # -> D will be singular
  # A <- JOINT@vcov$vcov[ step2.free.idx,  step2.free.idx]
  # B <- JOINT@vcov$vcov[ step2.free.idx, -step2.free.idx]
  # C <- JOINT@vcov$vcov[-step2.free.idx,  step2.free.idx]
  # D <- JOINT@vcov$vcov[-step2.free.idx, -step2.free.idx]
  # I.22.inv <- A - B %*% solve(D) %*% C

  # se = "standard"
  if (lavoptions$se == "standard") {
    vcov_1 <- 1 / n * i_22_inv
    out$VCOV <- vcov_1

  # se = "naive" or "local": grab VCOV directly from FIT.PA
  } else if (lavoptions$se %in% c("naive", "local", "local.nt")) {
    if (is.null(fit_pa@vcov$vcov)) {
      fit_pa@Options$se <- "standard"
      vcov_1 <- lavTech(fit_pa, "vcov")
    } else {
      vcov_1 <- fit_pa@vcov$vcov
    }
    if (length(step2_rm_idx) > 0L) {
      vcov_1 <- vcov_1[-step2_rm_idx, -step2_rm_idx]
    }
    # order rows/cols of VCOV, so that they correspond with the (step 2)
    # parameters of the JOINT model.
    # step2$pt.idx/pts.idx include rows for defined parameters (:=, <, >);
    # vcov_1 only has rows for free parameters, so filter those out first.
    is_free <- fit_pa@ParTable$free[step2$pts.idx] > 0L
    idx <- sort.int(step2$pt.idx[is_free], index.return = TRUE)$ix
    vcov_1 <- vcov_1[idx, idx]

    out$VCOV <- vcov_1

  # se = "twostep" or "twostep.robust"
  } else if (lavoptions$se  == "twostep" || lavoptions$se == "twostep.robust") {

    if (lavoptions$se  == "twostep") {
      v2 <- 1 / n * i_22_inv # not the same as FIT.PA@vcov$vcov!!
      v1 <- i_22_inv %*% i_21 %*% sigma_11 %*% i_12 %*% i_22_inv
    } else if (lavoptions$se == "twostep.robust") {
      # following Yuan & Chan 2002, eqs 4, 10, 11, 12, 13 and 14
      # but for V11, V12, V21, V22: we use index '1' for step1, and '2'
      # for step 2!!

      # a <- -1 * info[step2_free_idx, step2_free_idx, drop = FALSE]
      m_b <- -1 * info[step2_free_idx, step1_free_idx, drop = FALSE]

      # get P (for a single group!! for now)
      p <- lav_sam_step1_local_jac(step1 = step1, fit = fit, p_only = TRUE)

      # get V22
      if (is.null(joint@SampleStats@NACOV[[1]])) {
        joint@SampleStats@NACOV <- lavTech(joint, "gamma")
      }
      tmp <- lav_model_nvcov_robust_sem(
        lavmodel = joint@Model, lavsamplestats = joint@SampleStats,
        lavcache = joint@cache, lavdata = joint@Data,
        lavimplied = joint@implied, lavh1 = joint@h1,
        lavoptions = joint@Options, use_ginv = FALSE,
        attr_delta = TRUE, attr_t_dvgvd = TRUE, attr_e_inv = TRUE,
        attr_wls_v = TRUE)
      # nvar_cov <- tmp[, ] # remove attributes
      delta <- attr(tmp, "Delta")
      # e_inv <- attr(tmp, "E.inv")
      wls_v <- attr(tmp, "WLS.V")
      t_dvgvd <- attr(tmp, "tDVGVD")

      v22 <- t_dvgvd[step2_free_idx, step2_free_idx, drop = FALSE] # ok

      # FIXME: for a single group only:
      v11 <- p %*% lavTech(joint, "gamma")[[1]] %*% t(p)
      v21 <- t(delta[[1]][, step2_free_idx]) %*% wls_v[[1]] %*%
                                  lavTech(joint, "gamma")[[1]] %*% t(p)
      v12 <- t(v21)

      #V11 <- NVarCov[step1.free.idx, step1.free.idx, drop = FALSE]
      #V21 <- tmp2[   step2.free.idx, step1.free.idx, drop = FALSE]
      #PI <- V22 + B %*% V12 + V21 %*% t(B) + B %*% V11 %*% t(B)
      #A.inv <- solve(A)
      a_inv <-  -1 * i_22_inv
      #VCOV <- A.inv %*% PI %*% t(A.inv)
      v2 <- 1 / n * (a_inv %*% v22 %*% a_inv)
      v1 <- 1 / n * (a_inv %*%
        (m_b %*% v12 + v21 %*% t(m_b) + m_b %*% v11 %*% t(m_b)) %*% a_inv)
    }

    # V for second step
    if (!is.null(local_options$alpha.correction) &&
      local_options$alpha.correction > 0) {
      alpha_n1 <- local_options$alpha.correction / (n - 1)
      if (alpha_n1 > 1.0) {
        alpha_n1 <- 1.0
      } else if (alpha_n1 < 0.0) {
        alpha_n1 <- 0.0
      }
      if (is.null(fit_pa@vcov$vcov)) {
        fit_pa@Options$se <- "standard"
        vcov_naive <- lavTech(fit_pa, "vcov")
      } else {
        vcov_naive <- fit_pa@vcov$vcov
      }
      if (length(step2_rm_idx) > 0L) {
        vcov_naive <- vcov_naive[-step2_rm_idx, -step2_rm_idx]
      }
      vcov_corrected <- v2 + v1
      vcov_1 <- alpha_n1 * vcov_naive + (1 - alpha_n1) * vcov_corrected
    } else {
      # no alpha correction
      vcov_1 <- v2 + v1
    }

    # store in out
    out$V2 <- v2
    out$V1 <- v1
    out$VCOV <- vcov_1
  } # twostep

  # store se
  out$se <- lavoptions$se # in case it changed

  if (lav_verbose()) {
    cat("done.\n")
  }

  out
}
