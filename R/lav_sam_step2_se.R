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

# WLS.V_g %*% Delta_g (the 'WD' building block of the robust SEM sandwich),
# handling the DWLS/ULS case where lavaan stores the diagonal weight matrix as
# a plain vector. Optionally restrict Delta to a subset of columns.
lav_sam_wls_delta <- function(wls_v_g, delta_g, cols = NULL) {
  if (!is.null(cols)) {
    delta_g <- delta_g[, cols, drop = FALSE]
  }
  if (is.matrix(wls_v_g)) {
    wls_v_g %*% delta_g
  } else {
    # DWLS/ULS: WLS.V holds the diagonal of the weight matrix
    wls_v_g * delta_g
  }
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
    wd <- lav_sam_wls_delta(wls_v[[g]], delta[[g]])
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

# Casewise (huber.white) building blocks for the two-step sandwich
# (se = "twostep.huber.white"): the Yuan & Chan V22/V21/V11 pieces, but with
# casewise-SCORE crossproducts instead of Gamma-based expressions:
#
#   g2_i = casewise score of the JOINT model w.r.t. theta2, at (theta2.hat,
#          theta1.hat)
#   h_i  = casewise influence on theta1.hat: per measurement block b,
#          h_i^(b) = (N/n_b) I_b^{-1} g1_i^(b)   (unit block information I_b,
#          n_b = cases used by block b; zero rows for cases without data in b)
#
#   V22 = (1/N) sum g2_i g2_i',  V21 = (1/N) sum g2_i h_i',
#   V11 = (1/N) sum h_i h_i'
#
# Fed into the same A/B assembly as the twostep.robust branch, this yields the
# stacked M-estimation (Murphy-Topel) sandwich with first-order meat -- the
# two-step analogue of se = "robust.huber.white". Because the meat is built
# from casewise scores, it needs no Gamma: missing = "ml" (FIML scores per
# missing pattern) and multiple groups (zero scores outside a case's group)
# are covered by the same expressions.
#
# The scores are centered before the crossproducts: an exact no-op for
# sam.method = "global" (the joint theta2 score is solved: sum_i g2_i = 0, and
# sum_i g1_i = 0 per block), and the correct removal of the O(1) offset for
# the local methods (whose step-2 estimator does not solve the joint score
# exactly -- the same asymptotic-equivalence argument that underlies the
# twostep.robust linearization).
lav_sam_step2_se_hw_v <- function(joint = NULL, step1 = NULL,
                                  step2_free_idx = NULL, ntot = NULL) {
  npar <- max(joint@ParTable$free)
  # for the local sam methods, the joint object is assembled without fitting
  # and carries Options$estimator == "none"; the casewise scores are those of
  # the (ML) joint likelihood in Model@estimator
  if (!joint@Options$estimator %in% c("ML", "WLS", "GLS", "ULS")) {
    joint@Options$estimator <- joint@Model@estimator
  }
  # ignore_constraints: the sam models may carry synthesized INEQUALITY
  # constraints (the measurement blocks are fitted with bounds =
  # "wide.zerovar", which lavaan represents as linear cin rows in con.jac);
  # these are inactive at an interior solution, but lav_sc()'s constraint
  # projection would still project the scores onto them, corrupting the
  # casewise influences. Equality constraints are excluded from this se
  # flavour upstream, so skipping the projection is exact here.
  sc <- lav_sc(joint, remove_empty_cases = FALSE,
               ignore_constraints = TRUE)
  sc[is.na(sc)] <- 0 # empty cases contribute a zero score
  if (ncol(sc) != npar) {
    lav_msg_stop(gettext(
      "internal error: the joint casewise scores do not match the free
       parameters (twostep.huber.white)."))
  }
  g2 <- sc[, step2_free_idx, drop = FALSE]

  # step-1 casewise influence, in the joint free-parameter numbering
  h_full <- matrix(0, nrow(sc), npar)
  pt_free <- step1$PT.free
  for (mm in seq_along(step1$MM.FIT)) {
    fb <- step1$MM.FIT[[mm]]
    scb <- lav_sc(fb, remove_empty_cases = FALSE,
                  ignore_constraints = TRUE) # see note above
    scb[is.na(scb)] <- 0 # cases with no data in this block: zero score
    if (nrow(scb) != nrow(sc)) {
      lav_msg_stop(gettext(
        "internal error: case alignment failure between a measurement block
         and the joint model (twostep.huber.white)."))
    }
    ib <- lavTech(fb, "information") # unit information (honors information=)
    nb <- fb@SampleStats@ntotal
    hb <- (ntot / nb) * scb %*% solve(ib)
    # map the block's free parameters into the joint free numbering (the same
    # mapping used for the Sigma.11 assembly in lav_sam_step1())
    mm_idx <- step1$block.mm.idx[[mm]]
    ptm_idx <- step1$block.ptm.idx[[mm]]
    par_idx <- pt_free[mm_idx[ptm_idx]]
    keep_idx <- fb@ParTable$free[ptm_idx]
    ok <- par_idx > 0L & keep_idx > 0L # drop := etc.
    h_full[, par_idx[ok]] <- hb[, keep_idx[ok], drop = FALSE]
  }

  # center (see note above) and build the V pieces
  g2c <- t(t(g2) - colMeans(g2))
  h1c <- t(t(h_full) - colMeans(h_full))
  list(
    v22 = crossprod(g2c) / ntot,
    v21 = crossprod(g2c, h1c) / ntot, # cols = FULL joint numbering
    v11 = crossprod(h1c) / ntot # rows/cols = FULL joint numbering
  )
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

# the robust (sandwich) vcov of FIT.PA, using the NACOV (= Gamma.eta) it was
# fitted with -- FIT.PA itself was fitted with se = "standard" (so that the
# naive vcov is available for the alpha blending), so recompute here
lav_sam_step2_se_vcov_pa_robust <- function(fit_pa,
                                            step2_rm_idx = integer(0L)) {
  fit_pa@Options$se <- "robust.sem"
  fit_pa@vcov$vcov <- NULL # force recompute
  vcov_pa <- lavTech(fit_pa, "vcov")
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
      "twostep.huber.white", "local", "local.nt")) {
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

  # conditional.x + local sam method: the joint-score (Yuan & Chan)
  # linearization of the twostep.robust branch below describes the
  # joint-model-given-theta1 estimator, which under conditional.x (with
  # non-centered exogenous covariates) is NOT asymptotically equivalent to
  # the local step-2 estimator: the structural fit is saturated in the
  # latent-on-x slopes and cannot pool information across statistics the way
  # the joint estimator does, so the joint-score vcov underestimates the
  # sampling variability (badly so for binary covariates). Compute the
  # twostep.robust SEs as the structural-space sandwich instead (FIT.PA with
  # NACOV = Gamma.eta), which linearizes the actual two-step estimator and
  # makes them identical to se = "local" in this setting. (For sam.method =
  # "global" step 2 IS the joint estimator, and the branch below applies.)
  tsrobust_condx_flag <- lavoptions$se == "twostep.robust" &&
    fit@Model@conditional.x &&
    isTRUE(step1$sam.method %in% c("local", "fsr", "cfsr")) &&
    !is.null(step1$Gamma.eta) && !is.null(step1$Gamma.eta[[1]])

  # PML + twostep.huber.white: the assembled joint carries
  # Options$estimator == "none" and no PML cache (bifreq/long/uniweights),
  # while both the joint information (the numerical hessian of the PML
  # gradient) and the casewise scores need them. Patch a local copy of the
  # joint (R copy semantics: the returned object is not affected).
  if (lavoptions$se == "twostep.huber.white" &&
      joint@Model@estimator == "PML") {
    joint@Options$estimator <- "PML"
    if (is.null(joint@Cache[[1]]$long)) {
      joint@Cache <- lav_step10_cache(
        slot_cache = NULL,
        lavdata = joint@Data,
        lavmodel = joint@Model,
        lavpartable = joint@ParTable,
        lavoptions = joint@Options,
        sampling_weights = if (length(joint@Data@sampling.weights) > 0L) {
          joint@Data@sampling.weights
        } else {
          NULL
        }
      )
    }
  }

  if (lavoptions$se %in% c("naive", "twostep", "twostep.robust",
                           "twostep.huber.white")) {
    info <- lavInspect(joint, "information")
    i_12 <- info[step1_free_idx, step2_free_idx, drop = FALSE]
    i_22 <- info[step2_free_idx, step2_free_idx, drop = FALSE]
    i_21 <- info[step2_free_idx, step1_free_idx, drop = FALSE]
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
  if (lavoptions$se %in% c("standard", "twostep", "twostep.robust",
                           "twostep.huber.white")) {
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
      # 'extra' FIT.PA parameters that are not free in the JOINT model are
      # already in step2_rm_idx; only columns BEYOND n_step2 + those extras
      # are rotation-identification columns that must be removed in
      # addition. (Comparing n_jac_cols against n_step2 alone wrongly
      # trimmed legitimate constraint columns whenever FIT.PA had extra
      # parameters, eg multigroup group.equal = "regressions".)
      n_expected <- n_step2 + length(step2_rm_idx)
      if (n_jac_cols > n_expected) {
        aug_rm_idx <- union(step2_rm_idx, (n_expected + 1):n_jac_cols)
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
  # (also twostep.robust under conditional.x -- see tsrobust_condx_flag above)
  } else if (lavoptions$se %in% c("naive", "local", "local.nt") ||
             tsrobust_condx_flag) {
    if (isTRUE(step1$caseB) && lavoptions$se %in% c("local", "local.nt")) {
      # across-group constraints in the measurement model: build the
      # cross-group sandwich (the per-group FIT.PA vcov is incomplete)
      vcov_1 <- lav_sam_step2_se_local_caseB(fit_pa = fit_pa, step1 = step1)
      if (length(step2_rm_idx) > 0L) {
        vcov_1 <- vcov_1[-step2_rm_idx, -step2_rm_idx, drop = FALSE]
      }
    } else if (tsrobust_condx_flag) {
      vcov_1 <- lav_sam_step2_se_vcov_pa_robust(fit_pa, step2_rm_idx)
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

  # se = "twostep", "twostep.robust" or "twostep.huber.white"
  } else if (lavoptions$se %in% c("twostep", "twostep.robust",
                                  "twostep.huber.white")) {

    robust <- (lavoptions$se == "twostep.robust")
    huber_white <- (lavoptions$se == "twostep.huber.white")
    hw_v <- NULL
    if (huber_white) {
      # casewise-score V pieces; on failure fall back to plain twostep
      hw_v <- tryCatch(
        lav_sam_step2_se_hw_v(
          joint = joint, step1 = step1,
          step2_free_idx = step2_free_idx, ntot = n
        ),
        error = function(e) e
      )
      if (inherits(hw_v, "error")) {
        lav_msg_warn(gettextf(
          "huber.white standard errors (se = \"twostep.huber.white\") could
           not be computed (%s); twostep standard errors are reported
           instead.", conditionMessage(hw_v)))
        huber_white <- FALSE
        lavoptions$se <- "twostep"
      }
    }
    if (robust && fit@Model@conditional.x &&
        isTRUE(step1$sam.method %in% c("local", "fsr", "cfsr"))) {
      # only reached when Gamma.eta is unavailable (tsrobust_condx_flag was
      # FALSE): the joint-score correction below is not valid for the local
      # step-2 estimator under conditional.x -> plain twostep
      lav_msg_warn(gettext(
        "robust standard errors (se = \"twostep.robust\") need Gamma.eta
         under conditional.x = TRUE, which could not be computed; twostep
         standard errors are reported instead."))
      robust <- FALSE
      lavoptions$se <- "twostep"
    }
    if (robust) {
      # get P (the influence of the step 1 / measurement parameters on the
      # sample statistics) and align its rows with step1.free.idx (the rows of
      # P are in ascending free-parameter order). Single group: one matrix;
      # multigroup: the columns are the stacked statistics of all groups, with
      # attr 'stat.offset' giving the per-group column boundaries.
      p <- lav_sam_step1_local_jac(step1 = step1, fit = fit, p_only = TRUE)
      p_row_idx <- match(step1_free_idx, attr(p, "free.idx"))
      if (anyNA(p_row_idx)) {
        # P cannot be aligned with the step 1 free parameters: the robust
        # correction is not available -> fall back to (non-robust) twostep
        # SEs. (Should not happen: within-block and across-group equality
        # constraints are handled; this is a defensive guard.)
        lav_msg_warn(gettext(
          "robust standard errors (se = \"twostep.robust\") are not available
           for this model; twostep standard errors are reported instead."))
        robust <- FALSE
        lavoptions$se <- "twostep"
      }
    }
    if (huber_white) {
      # stacked M-estimation sandwich with casewise-score meat: same A/B
      # assembly as the twostep.robust branch below, with the V pieces
      # replaced by their casewise counterparts (see lav_sam_step2_se_hw_v)
      m_b <- -1 * info[step2_free_idx, step1_free_idx, drop = FALSE]
      v22 <- hw_v$v22
      v21 <- hw_v$v21[, step1_free_idx, drop = FALSE]
      v11 <- hw_v$v11[step1_free_idx, step1_free_idx, drop = FALSE]
      v12 <- t(v21)

      a_inv <- -1 * i_22_inv
      v2 <- 1 / n * (a_inv %*% v22 %*% a_inv)
      v1 <- 1 / n * (a_inv %*%
        (m_b %*% v12 + v21 %*% t(m_b) + m_b %*% v11 %*% t(m_b)) %*% a_inv)
    } else if (!robust) {
      v2 <- 1 / n * i_22_inv # not the same as FIT.PA@vcov$vcov!!
      v1 <- i_22_inv %*% i_21 %*% sigma_11 %*% i_12 %*% i_22_inv
    } else {
      # following Yuan & Chan 2002, eqs 4, 10, 11, 12, 13 and 14
      # but for V11, V12, V21, V22: we use index '1' for step1, and '2'
      # for step 2!!
      m_b <- -1 * info[step2_free_idx, step1_free_idx, drop = FALSE]
      stat_offset <- attr(p, "stat.offset") # NULL if single group
      p <- p[p_row_idx, , drop = FALSE]

      # get V22 (= the robust covariance of the step 2 estimating function,
      # already aggregated over groups by lav_model_nvcov_robust_sem())
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

      # V11 = Cov(sqrt(N) * step1 params) and V21 = Cov(sqrt(N) * step2 score,
      # sqrt(N) * step1 params). Because the groups are independent samples,
      # Gamma is block-diagonal across groups, so both aggregate as a sum over
      # groups. The group weights differ per block, however, because the step 2
      # estimating function G is the sample-size weighted average of the
      # per-group scores (G = sum_g f_g G_g), while the step 1 estimates gamma_g
      # are per-group. With Gamma_g = Cov(sqrt(n_g) s_g) and f_g = n_g/N:
      #   V22 = Cov(sqrt(N) G)         -> sum_g f_g   WD_g' Gamma_g WD_g  (= t_dvgvd)
      #   V21 = Cov(sqrt(N) G, gamma)  -> sum_g  1    WD_g' Gamma_g P_g'
      #   V11 = Cov(sqrt(N) gamma)     -> sum_g 1/f_g P_g  Gamma_g P_g'
      # (each replacement of a 'score' index by a 'parameter' index trades one
      # factor f_g for 1/f_g). This makes the three terms of V1 below scale
      # consistently, so for independent groups (Case A) the per-group result
      # equals a separate single-group fit. It also covers across-group
      # measurement constraints (Case B): the cross-group coupling is encoded
      # in the per-group columns of P. For a single group f_g = 1 and the loop
      # reduces to the original (unweighted) single-group expressions.
      gamma <- lavTech(joint, "gamma")
      ngroups <- joint@Data@ngroups
      ntot <- joint@SampleStats@ntotal
      v11 <- matrix(0, nrow(p), nrow(p))
      v21 <- matrix(0, length(step2_free_idx), nrow(p))
      for (g in seq_len(ngroups)) {
        fg <- joint@SampleStats@nobs[[g]] / ntot
        gamma_g <- gamma[[g]]
        if (is.null(stat_offset)) {
          p_g <- p
        } else {
          p_g <- p[, (stat_offset[g] + 1L):stat_offset[g + 1L], drop = FALSE]
        }
        wd_g <- lav_sam_wls_delta(wls_v[[g]], delta[[g]], step2_free_idx)
        v11 <- v11 + (1 / fg) * (p_g %*% gamma_g %*% t(p_g))
        v21 <- v21 + crossprod(wd_g, gamma_g %*% t(p_g))
      }
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
