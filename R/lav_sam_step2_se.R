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

# se = "bootstrap": bootstrap the whole two-step procedure, and compute
# the vcov of the (combined step 1 + step 2) coefficients
lav_sam_step2_se_bootstrap <- function(sam_object = NULL, bootstrap = list()) {
  default_args <- list(R = 1000L, type = "ordinary",
                       show.progress = FALSE,
                       check.post = TRUE, keep.idx = FALSE)
  this_args <- modifyList(default_args, bootstrap)
  coef_1 <- lav_bootstrap_internal(object = sam_object,
    r = this_args$R, show_progress = this_args$show.progress,
    type = this_args$type, fun = "coef",
    check_post = this_args$check.post, keep_idx = this_args$keep.idx)
  coef_orig <- coef_1
  error_idx <- attr(coef_1, "error.idx")
  nfailed <- length(error_idx) # zero if NULL
  if (nfailed > 0L) {
    lav_msg_warn(gettextf(
      "%s bootstrap runs failed or did not converge.", nfailed))
  }
  notok <- length(attr(coef_1, "nonadmissible")) # zero if NULL
  if (notok > 0L) {
    lav_msg_warn(gettextf(
      "%s bootstrap runs resulted in nonadmissible solutions.", notok))
  }
  if (length(error_idx) > 0L) {
    # new in 0.6-13: we must still remove them!
    coef_1 <- coef_1[-error_idx, , drop = FALSE]
    # this also drops the attributes
  }
  nboot <- nrow(coef_1)
  var_cov <- cov(coef_1) * (nboot - 1) / nboot

  list(VCOV = var_cov, boot.coef = coef_orig, R = this_args$R)
}

# Case B local SEs: across-group (equality) constraints in the measurement
# model couple vech(VETA_g) across groups, so the NACOV of the stacked
# structural statistics (step1$Gamma.eta.full) has nonzero cross-group blocks
# that lavaan's multigroup robust.sem (a per-group NACOV list) cannot use.
# We therefore compute the robust sandwich ourselves, reusing the structural
# model's Delta, WLS.V (weight) and the augmented inverse information E.inv:
#
#   VCOV = (1/N^2) E.inv [ WD' (Dn Cov.full Dn) WD ] E.inv
#
# with WD = rbind_g( WLS.V_g %*% Delta_g ), Dn = diag(n_g) (repeated over the
# statistics of group g), Cov.full = the full covariance of the stacked
# structural statistics. When Cov.full is block-diagonal (Case A) this reduces
# exactly to lavaan's E.inv %*% (sum_g (n_g/N) WD_g' Gamma_g WD_g) %*% E.inv / N.
lav_sam_step2_se_local_caseB <- function(fit_pa, step1) {
  ngroups <- fit_pa@Data@ngroups
  ntotal  <- fit_pa@SampleStats@ntotal
  nobs_g  <- unlist(fit_pa@SampleStats@nobs)

  # structural model Delta, WLS.V (weight) and augmented E.inv
  tmp <- lav_model_nvcov_robust_sem(
    lavmodel = fit_pa@Model, lavsamplestats = fit_pa@SampleStats,
    lavcache = fit_pa@cache, lavdata = fit_pa@Data,
    lavimplied = fit_pa@implied, lavh1 = fit_pa@h1,
    lavoptions = fit_pa@Options, use_ginv = FALSE,
    attr_delta = TRUE, attr_e_inv = TRUE, attr_wls_v = TRUE)
  e_inv <- attr(tmp, "E.inv")
  delta <- attr(tmp, "Delta")
  wls_v <- attr(tmp, "WLS.V")

  # WD = rbind_g( W_g Delta_g ); dn = n_g repeated over group g's statistics
  wd_list <- vector("list", ngroups)
  dn_list <- vector("list", ngroups)
  for (g in seq_len(ngroups)) {
    if (is.matrix(wls_v[[g]])) {
      wd <- wls_v[[g]] %*% delta[[g]]
    } else {
      # DWLS/ULS: WLS.V holds the diagonal of the weight matrix
      wd <- wls_v[[g]] * delta[[g]]
    }
    wd_list[[g]] <- wd
    dn_list[[g]] <- rep.int(nobs_g[g], nrow(wd))
  }
  wd_all <- do.call(rbind, wd_list)
  dn <- unlist(dn_list)

  cov_full <- step1$Gamma.eta.full
  if (nrow(cov_full) != length(dn)) {
    lav_msg_stop(gettext(
      "internal error: dimension mismatch while assembling the multigroup
       (Case B) local sandwich."))
  }

  # M = WD' (Dn Cov.full Dn) WD  ;  (Dn Cov.full Dn)[i,j] = n_i cov[i,j] n_j
  dcd <- cov_full * outer(dn, dn)
  m_cov <- crossprod(wd_all, dcd %*% wd_all)
  vcov_1 <- (1 / (ntotal * ntotal)) * (e_inv %*% m_cov %*% e_inv)
  vcov_1
}

# grab the 'naive' vcov from FIT.PA (computing it if absent), removing
# any rows/cols in step2_rm_idx
lav_sam_step2_se_vcov_pa <- function(fit_pa, step2_rm_idx = integer(0L)) {
  if (is.null(fit_pa@vcov$vcov)) {
    fit_pa@Options$se <- "standard"
    vcov_pa <- lavTech(fit_pa, "vcov")
  } else {
    vcov_pa <- fit_pa@vcov$vcov
  }
  if (length(step2_rm_idx) > 0L) {
    vcov_pa <- vcov_pa[-step2_rm_idx, -step2_rm_idx, drop = FALSE]
  }
  vcov_pa
}

# Standard errors for latent (residual) variances that are 'fixed' (free == 0)
# in the JOINT model but are in fact data-derived. Under std.lv = TRUE, SAM
# step 1 fixes the *total* variance of each factor to 1; for an *endogenous*
# factor the residual variance is then a derived quantity
#   r_j = total_j - explained_j   with   VETA = (I - B)^-1 Psi (I - B)^-T
# (total_j fixed to 1). Although r_j is fixed-by-label in the joint model (and
# therefore excluded from the joint vcov, giving se = 0), it is a function of
# the free structural parameters and has genuine sampling variability. We
# recover its SE by the delta method, using a numerical jacobian of r_j
# through the model-implied VETA and the (se-method-specific) joint vcov.
# Exogenous factor variances have no incoming paths, so their jacobian is zero
# and they keep se = 0 automatically.
#
# Returns a numeric vector aligned with the rows of joint@ParTable (0 for every
# row that is not such a derived latent variance).
lav_sam_step2_se_lv_var <- function(joint = NULL) {
  lavmodel    <- joint@Model
  lavpartable <- joint@ParTable
  VCOV        <- joint@vcov$vcov
  se_out      <- numeric(length(lavpartable$lhs))

  if (is.null(VCOV) || lavmodel@representation != "LISREL" || anyNA(VCOV)) {
    return(se_out)
  }

  # 'psi' matrix index in GLIST per block, plus its latent-variable names
  nmat      <- lavmodel@nmat
  mm_idx    <- lav_model_group_mm_indices(nmat)
  psi_glist <- which(names(lavmodel@GLIST) == "psi")
  if (length(psi_glist) == 0L) {
    return(se_out)
  }
  lv_block <- lapply(seq_len(lavmodel@nblocks), function(g) {
    psi_mm <- intersect(mm_idx[[g]], psi_glist)
    lavmodel@dimNames[[psi_mm[1]]][[1]]
  })

  # target rows: latent (residual) variances, fixed (free == 0), and not
  # explicitly fixed by the user (user != 1)
  target <- which(
    lavpartable$op == "~~" &
      lavpartable$lhs == lavpartable$rhs &
      lavpartable$free == 0L &
      lavpartable$user != 1L
  )
  target <- target[vapply(target, function(i) {
    lavpartable$lhs[i] %in% lv_block[[lavpartable$block[i]]]
  }, logical(1L))]
  if (length(target) == 0L) {
    return(se_out)
  }

  # parameter vector in the space matching VCOV (unco if simple eq constraints)
  ceq_simple <- lavmodel@ceq.simple.only &&
    nrow(lavmodel@ceq.simple.K) == nrow(VCOV)
  vtype <- if (ceq_simple) "unco" else "free"
  x0 <- lav_model_get_parameters(lavmodel) # nx.free
  if (ceq_simple) {
    x0 <- drop(lavmodel@ceq.simple.K %*% x0) # nx.unco
  }
  npar <- length(x0)
  if (npar != ncol(VCOV)) { # spaces do not line up: stay safe
    return(se_out)
  }

  # explained_j: model-implied total variance of factor j in block g, with its
  # own residual psi[j, j] set to zero
  explained <- function(x, g, j) {
    gl <- lav_model_x2glist(lavmodel, x = x, type = vtype)
    psi_mm <- intersect(mm_idx[[g]], psi_glist)[1]
    gl[[psi_mm]][j, j] <- 0
    lav_model_veta(lavmodel = lavmodel, glist = gl)[[g]][j, j]
  }

  eps <- 1e-6
  for (row in target) {
    g <- lavpartable$block[row]
    j <- match(lavpartable$lhs[row], lv_block[[g]])
    if (is.na(j)) next
    v <- try(
      {
        e0 <- explained(x0, g, j)
        jac <- numeric(npar)
        for (k in seq_len(npar)) {
          xk <- x0
          xk[k] <- xk[k] + eps
          jac[k] <- -(explained(xk, g, j) - e0) / eps
        }
        as.numeric(crossprod(jac, VCOV %*% jac))
      },
      silent = TRUE
    )
    if (!inherits(v, "try-error") && is.finite(v) && v >= 0) {
      se_out[row] <- sqrt(v)
    }
  }

  se_out
}

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

  # invert augmented information, for I.22 block only
  # new in 0.6-16 (otherwise, eq constraints in struc part are ignored)
  if (lavoptions$se %in% c("standard", "twostep", "twostep.robust")) {
    # Fix for EFA/ESEM: when rotation is used, FIT.PA@Model@con.jac includes
    # columns for rotation identification constraints that are not part of
    # step2.free.idx. These extra columns cause a dimension mismatch in
    # lav_model_info_augment_invert(). Remove them via rm.idx. NOTE: this
    # adjustment is only relevant for the augmented-inverse path below; it must
    # NOT be added to step2_rm_idx itself, because the naive/local path reads
    # the FIT.PA vcov directly and would then drop legitimate parameters (this
    # used to break multigroup models with across-group equality constraints,
    # where con.jac legitimately has more columns than step2.free.idx).
    aug_rm_idx <- step2_rm_idx
    if (nrow(fit_pa@Model@con.jac) > 0L) {
      n_jac_cols <- ncol(fit_pa@Model@con.jac)
      n_step2 <- length(step2_free_idx)
      if (n_jac_cols > n_step2) {
        aug_rm_idx <- union(step2_rm_idx, (n_step2 + 1):n_jac_cols)
      }
    }
    i_22_inv <-
      lav_model_info_augment_invert(
        lavmodel = fit_pa@Model,
        information = i_22,
        inverted = TRUE,
        use_ginv = FALSE, # if interaction, SEs end up smaller than naive...
        rm_idx = aug_rm_idx
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

  # note: computing I.22.inv as a Schur complement of JOINT@vcov$vcov
  # would allow a 'robust' vcov for the joint model, but does not work if
  # we have equality constraints in the MM (the (1,1) block is singular)

  # se = "standard"
  if (lavoptions$se == "standard") {
    vcov_1 <- 1 / n * i_22_inv
    out$VCOV <- vcov_1

  # se = "naive" or "local": grab VCOV directly from FIT.PA
  } else if (lavoptions$se %in% c("naive", "local", "local.nt")) {
    if (isTRUE(step1$caseB) && lavoptions$se %in% c("local", "local.nt")) {
      # across-group constraints in the measurement model: build the
      # cross-group sandwich (the per-group FIT.PA vcov is incomplete)
      vcov_1 <- lav_sam_step2_se_local_caseB(fit_pa = fit_pa, step1 = step1)
      if (length(step2_rm_idx) > 0L) {
        vcov_1 <- vcov_1[-step2_rm_idx, -step2_rm_idx, drop = FALSE]
      }
    } else {
      vcov_1 <- lav_sam_step2_se_vcov_pa(fit_pa, step2_rm_idx)
    }
    # order rows/cols of VCOV, so that they correspond with the (step 2)
    # parameters of the JOINT model, but remove := parameters first
    pt_idx <- step2$pt.idx
    nonpar_id <- fit@ParTable$id[fit@ParTable$op %in% c(":=", "==", "<", ">")]
    rm_idx <- which(step2$pt.idx %in% nonpar_id)
    if (length(rm_idx) > 0L) {
      pt_idx <- step2$pt.idx[-rm_idx]
    }
    idx <- sort.int(pt_idx, index.return = TRUE)$ix
    vcov_1 <- vcov_1[idx, idx, drop = FALSE]

    # drop parameters that are free in FIT.PA, but fixed in the JOINT
    # model (eg std.lv = TRUE: the lv (residual) variances are freed in
    # the structural model, but fixed (to 1) in the joint model); there is
    # no room for their sampling variability in the joint vcov, and
    # step2.free.idx does not include them
    pt_2 <- step2$PT
    keep_idx <- which(pt_2$free[sort.int(pt_idx)] > 0L)
    if (length(keep_idx) < length(pt_idx)) {
      vcov_1 <- vcov_1[keep_idx, keep_idx, drop = FALSE]
    }

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
      m_b <- -1 * info[step2_free_idx, step1_free_idx, drop = FALSE]

      # get P (for a single group!! for now)
      p <- lav_sam_step1_local_jac(step1 = step1, fit = fit, p_only = TRUE)
      # align the rows of P with step1.free.idx: the rows of P are in
      # ascending free-parameter order, while step1.free.idx (and hence
      # m_b below) is ordered per measurement block
      p_row_idx <- match(step1_free_idx, attr(p, "free.idx"))
      if (anyNA(p_row_idx)) {
        lav_msg_stop(gettext(
          "internal error: unable to align the step 1 jacobian (P) with
           the step 1 free parameters"))
      }
      p <- p[p_row_idx, , drop = FALSE]

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
      delta <- attr(tmp, "Delta")
      wls_v <- attr(tmp, "WLS.V")
      t_dvgvd <- attr(tmp, "tDVGVD")

      v22 <- t_dvgvd[step2_free_idx, step2_free_idx, drop = FALSE] # ok

      # FIXME: for a single group only:
      gamma_1 <- lavTech(joint, "gamma")[[1]]
      v11 <- p %*% gamma_1 %*% t(p)
      if (is.matrix(wls_v[[1]])) {
        wd <- wls_v[[1]] %*% delta[[1]][, step2_free_idx, drop = FALSE]
      } else {
        # DWLS/ULS: WLS.V holds the diagonal of the weight matrix
        wd <- wls_v[[1]] * delta[[1]][, step2_free_idx, drop = FALSE]
      }
      v21 <- t(wd) %*% gamma_1 %*% t(p)
      v12 <- t(v21)

      a_inv <- -1 * i_22_inv
      v2 <- 1 / n * (a_inv %*% v22 %*% a_inv)
      v1 <- 1 / n * (a_inv %*%
        (m_b %*% v12 + v21 %*% t(m_b) + m_b %*% v11 %*% t(m_b)) %*% a_inv)
    }

    # V for second step
    if (!is.null(local_options$alpha.correction) &&
      local_options$alpha.correction > 0) {
      alpha_n1 <- lav_sam_alpha_n1(local_options$alpha.correction, n)
      vcov_naive <- lav_sam_step2_se_vcov_pa(fit_pa, step2_rm_idx)
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
