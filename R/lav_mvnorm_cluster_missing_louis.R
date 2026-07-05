# observed information for the two-level model with missing data,
# using Louis's (1982) method
#
# YR + Claude, July 2026
#
# Louis's identity:
#   -H_obs(theta) = E[-H_complete | v.obs] - Cov[S_complete | v.obs]
# where C is the complete-data loglikelihood: treating the cluster random
# effects, the missing y-values and the missing z-values as 'complete',
# C is simply the sum of two independent Gaussian samples (a two-group
# model):
#   d_ji = y_ji - nu_j ~ N(Mu.W, Sigma.W)   (one term per level-1 unit)
#   u_j  = (nu_j[both], z_j) ~ N(Mu.B, Sigma.B) (one term per cluster)
# (cf. Asparouhov & Muthen, 2003)
#
# the first term only requires the expected complete-data moments (one
# E-step); the second term is computed in closed form below, using the
# posterior distribution of the cluster-level unknowns (the same
# quantities as the E-step of lav_mvn_cl_mi_em_engine)
#
# this replaces the numerical Hessian (4 x npar full-data gradient
# passes) by a single full-data pass

# ---------------------------------------------------------------------
# Cov[ S_complete | v.obs ] in the *moment* space, where S_complete is
# the complete-data score of (-2 * logl) with respect to
# (Mu.W, vech(Sigma.W), Mu.B, vech(Sigma.B)) -- the same packing (and
# scale) as lav_mvn_cl_mi_dlogl_2l_samp()
#
# notation, per cluster j (centered): w = (beta, ztilde[mis]) with
# posterior N(m1, C1); per unit i with missing pattern p:
#   e_ji = d_ji - Mu.W = a_i + F_p beta + eta_i,  eta_i ~ N(0, V_p)
# where a_i = G_p yc_i[obs], G_p = [I; B_p] (observed rows identity,
# missing rows the sweep regression), F_p = -G_p S_ob (S_ob selects the
# observed both-columns), V_p = C_p embedded in the missing block;
# and t_j = u_j - Mu.B = k_j + K_q w (k_j holds the observed ztilde)
#
# the law of total covariance then gives two parts:
# (1) E_w[ Cov(S | w) ]: only the within scores are random given w;
#     unit-wise Magnus-Neudecker moments, aggregated per pattern
# (2) Cov_w( E[S | w] ): E[S|w] is a quadratic polynomial in w; for
#     centered w ~ N(0, C1): Cov(b'w + w'Rw, ...) = b1'C b2 +
#     2 tr(R1 C R2 C)
lav_mvn_cl_mi_scov_louis <- function(y1 = NULL,
                                     y2 = NULL,
                                     lp = NULL,
                                     mp = NULL,
                                     mu_w = NULL, # implied Mu.W
                                     sigma_w = NULL, # implied Sigma.W
                                     mu_b = NULL, # implied Mu.B
                                     sigma_b = NULL, # implied Sigma.B
                                     extra = FALSE) {
  # map implied to 2l matrices (for the posterior machinery)
  out <- lav_mvn_cl_implied22l(
    lp = lp, mu_w = mu_w, mu_b = mu_b,
    sigma_w = sigma_w, sigma_b = sigma_b
  )
  mu_y <- out$mu.y
  mu_z <- out$mu.z
  sw <- out$sigma.w # == sigma_w (within/y1w space)
  sb <- out$sigma.b[lp$both.idx[[2]], lp$both.idx[[2]], drop = FALSE]
  sigma_zz <- out$sigma.zz
  syz <- out$sigma.yz[lp$both.idx[[2]], , drop = FALSE]

  # dimensions and indices
  nobs <- lp$nclusters[[1]]
  nclusters <- lp$nclusters[[2]]
  between_idx <- lp$between.idx[[2]]
  both_idx <- lp$both.idx[[2]]
  nz <- length(between_idx)
  nb <- length(both_idx)
  zp <- mp$Zp
  if (nz > 0L) {
    y1w <- y1[, -between_idx, drop = FALSE]
    z <- y2[, between_idx, drop = FALSE]
    zc <- t(t(z) - mu_z)
  } else {
    y1w <- y1
    z <- NULL
  }
  ny <- ncol(y1w)
  p2 <- nb + nz
  y1w_c <- t(t(y1w) - mu_y)

  # score-space dimensions and offsets:
  # [Mu.W (ny), vech Sigma.W, Mu.B (p2), vech Sigma.B]
  nvw <- ny * (ny + 1L) / 2L
  nvb <- p2 * (p2 + 1L) / 2L
  ntot <- ny + nvw + p2 + nvb
  i_muw <- seq_len(ny)
  i_sw <- ny + seq_len(nvw)
  i_mub <- ny + nvw + seq_len(p2)
  i_sb <- ny + nvw + p2 + seq_len(nvb)

  # constants
  ow <- solve.default(sw)
  ob <- solve.default(sigma_b) # implied Sigma.B == [[sb, syz],[szy, szz]]
  gow <- lav_mat_dup_pre(ow %x% ow) # Gdd_w (Omega x Omega): nvw x ny^2
  gob <- lav_mat_dup_pre(ob %x% ob) # nvb x p2^2
  kny <- lav_mat_com(ny, ny)

  # ---- posterior machinery (same as the E-step) --------------------
  # (shared with lav_mvn_cl_mi_estep_ranef)
  post <- lav_mvn_cl_mi_posterior(
    y1 = y1, y2 = y2, lp = lp, mp = mp,
    mu_w = mu_w, sigma_w = sigma_w, mu_b = mu_b, sigma_b = sigma_b
  )
  breg_p <- post$breg
  ccond_p <- post$ccond
  m1_list <- post$m1
  c1_list <- post$c1
  zmis_q <- post$zmis_q
  zpat2j <- post$zpat2j
  cluster_idx <- lp$cluster.idx[[2]]

  # ---- per-pattern constants + data aggregates ---------------------
  # G_p, F_p, V_p; per (pattern, cluster): n_pj and sum of yc[obs]
  np <- mp$npatterns
  g_list <- f_list <- v_list <- vector("list", np)
  ycsum_list <- jmap_list <- nvec_list <- vector("list", np)
  for (p in seq_len(np)) {
    o_idx <- which(mp$pat[p, ])
    na_idx <- which(!mp$pat[p, ])
    g_p <- matrix(0, ny, length(o_idx))
    g_p[o_idx, ] <- diag(length(o_idx))
    if (length(na_idx) > 0L) {
      g_p[na_idx, ] <- breg_p[[p]]
    }
    g_list[[p]] <- g_p
    # F_p = -G_p S_ob: columns for both-vars; zero if unobserved
    f_p <- matrix(0, ny, nb)
    obs_b <- match(both_idx, o_idx) # NA if this both-var is missing
    seen <- which(!is.na(obs_b))
    if (length(seen) > 0L) {
      f_p[, seen] <- -g_p[, obs_b[seen], drop = FALSE]
    }
    f_list[[p]] <- f_p
    v_p <- matrix(0, ny, ny)
    if (length(na_idx) > 0L) {
      v_p[na_idx, na_idx] <- ccond_p[[p]]
    }
    v_list[[p]] <- v_p
    # per-cluster aggregates for this pattern
    case_idx <- mp$case.idx[[p]]
    j_idx <- mp$j.idx[[p]]
    ycs <- rowsum.default(y1w_c[case_idx, o_idx, drop = FALSE],
      group = j_idx, reorder = FALSE
    )
    jmap_list[[p]] <- as.integer(rownames(ycs))
    ycsum_list[[p]] <- ycs
    nvec <- rowsum.default(rep(1, length(case_idx)),
      group = j_idx, reorder = FALSE
    )
    nvec_list[[p]] <- drop(nvec)
  }

  # ---- part 1: E_w[ Cov(S | w) ] -- per-pattern totals -------------
  # per pattern we need (summing over all units in the pattern):
  #   chat_tot = sum_i E[c_i]  and  cc_tot = sum_i E[c_i c_i']
  # with c_i = a_i + F_p beta, beta | v.obs ~ N(m1[1:nb], C1[bb])
  scov <- matrix(0, ntot, ntot)
  escore <- numeric(ntot) # E[S | v.obs] (Fisher identity check)
  se_tot <- matrix(0, ny, ny) # sum_i E[e e'] (for escore)
  chat_all <- numeric(ny) # sum_i E[c_i] over everything
  for (p in seq_len(np)) {
    g_p <- g_list[[p]]
    f_p <- f_list[[p]]
    v_p <- v_list[[p]]
    jm <- jmap_list[[p]]
    nvec <- nvec_list[[p]]
    ycs <- ycsum_list[[p]]
    o_idx <- which(mp$pat[p, ])
    freq <- mp$freq[p]

    # pattern totals
    m1b <- do.call(rbind, lapply(m1_list[jm], function(m) m[seq_len(nb)]))
    if (nb == 0L) {
      m1b <- matrix(0, length(jm), 0L)
    }
    # chat_tot = G_p sum_j ycsum_pj + F_p sum_j n_pj m1b_j
    chat_tot <- drop(g_p %*% colSums(ycs)) +
      drop(f_p %*% colSums(nvec * m1b))
    # cc_tot = G (sum yc yc') G' + cross terms + n-weighted quadratic
    yy <- crossprod(y1w_c[mp$case.idx[[p]], o_idx, drop = FALSE])
    m1s <- crossprod(ycs, m1b) # sum_j ycsum_pj m1b_j'
    w2 <- matrix(0, nb, nb) # sum_j n_pj (m1b m1b' + C1[bb])
    for (k in seq_along(jm)) {
      c1bb <- c1_list[[jm[k]]][seq_len(nb), seq_len(nb), drop = FALSE]
      w2 <- w2 + nvec[k] * (tcrossprod(m1b[k, ]) + c1bb)
    }
    cross <- g_p %*% m1s %*% t(f_p)
    cc_tot <- g_p %*% yy %*% t(g_p) + cross + t(cross) +
      f_p %*% w2 %*% t(f_p)

    # (mu.w, mu.w): 4 freq * Ow V Ow
    scov[i_muw, i_muw] <- scov[i_muw, i_muw] +
      4 * freq * (ow %*% v_p %*% ow)
    # (mu.w, sigma.w): 2 Ow [V x chat' + chat' x V] t(gow)
    m_x <- (v_p %x% t(chat_tot)) + (t(chat_tot) %x% v_p)
    scov[i_muw, i_sw] <- scov[i_muw, i_sw] + 2 * (ow %*% m_x %*% t(gow))
    # (sigma.w, sigma.w): gow (I+K)[freq VxV + V x cc + cc x V] t(gow)
    m_q <- freq * (v_p %x% v_p) + (v_p %x% cc_tot) + (cc_tot %x% v_p)
    m_q <- m_q + kny %*% m_q
    scov[i_sw, i_sw] <- scov[i_sw, i_sw] + gow %*% m_q %*% t(gow)

    # bookkeeping for the expected score
    se_tot <- se_tot + cc_tot + freq * v_p
    chat_all <- chat_all + chat_tot
  }
  # empty units: e ~ N(0, Sigma.W)
  n_empty <- length(mp$empty.idx)
  if (n_empty > 0L) {
    scov[i_muw, i_muw] <- scov[i_muw, i_muw] + 4 * n_empty * ow
    m_q <- sw %x% sw
    m_q <- m_q + kny %*% m_q
    scov[i_sw, i_sw] <- scov[i_sw, i_sw] + n_empty * (gow %*% m_q %*% t(gow))
    se_tot <- se_tot + n_empty * sw
  }

  # ---- part 2: Cov_w( E[S | w] ) -- per-cluster --------------------
  # E[S|w] = const + b'wtilde + wtilde' R wtilde, wtilde ~ N(0, C1)
  et_tot <- matrix(0, p2, p2) # sum_j E[t t'] (for escore)
  that_all <- numeric(p2)
  for (j in seq_len(nclusters)) {
    m1 <- m1_list[[j]]
    c1 <- c1_list[[j]]
    nw <- length(m1)
    nzm <- nw - nb
    # K_q: w -> t (p2 x nw)
    k_q <- matrix(0, p2, nw)
    if (nb > 0L) {
      k_q[seq_len(nb), seq_len(nb)] <- diag(nb)
    }
    if (nzm > 0L) {
      zm_idx <- zmis_q[[zpat2j[j]]]
      k_q[nb + zm_idx, nb + seq_len(nzm)] <- diag(nzm)
    }
    # that = k_j + K m1
    that <- numeric(p2)
    if (nz > 0L) {
      zrow <- zc[j, ]
      zrow[is.na(zrow)] <- 0
      that[nb + seq_len(nz)] <- zrow
    }
    that <- that + drop(k_q %*% m1)

    # B matrix (ntot x nw) and R stack (ntot x nw^2)
    bmat <- matrix(0, ntot, nw)
    rstk <- matrix(0, ntot, nw * nw)
    # within rows: sums over the patterns present in this cluster
    if (nb > 0L) {
      fsum <- matrix(0, ny, nb) # sum_p n_pj F_p
      usum <- matrix(0, ny * ny, nb) # sum_p (F x chat_pj + chat_pj x F)
      ffsum <- matrix(0, ny * ny, nb * nb) # sum_p n_pj (F x F)
      m1b <- m1[seq_len(nb)]
      for (p in seq_len(np)) {
        k <- match(j, jmap_list[[p]])
        if (is.na(k)) {
          next
        }
        n_pj <- nvec_list[[p]][k]
        f_p <- f_list[[p]]
        chat_pj <- drop(g_list[[p]] %*% ycsum_list[[p]][k, ]) +
          n_pj * drop(f_p %*% m1b)
        fsum <- fsum + n_pj * f_p
        usum <- usum + (f_p %x% chat_pj) + (chat_pj %x% f_p)
        ffsum <- ffsum + n_pj * (f_p %x% f_p)
      }
      bmat[i_muw, seq_len(nb)] <- -2 * (ow %*% fsum)
      bmat[i_sw, seq_len(nb)] <- -(gow %*% usum)
      # quadratic rows: R = reshape of -(gow ffsum), beta-block only;
      # embed the nb x nb block into the nw x nw vec-space
      rsw <- -(gow %*% ffsum) # nvw x nb^2
      bpos <- as.vector(outer(seq_len(nb), seq_len(nb),
        function(r, s) (s - 1L) * nw + r
      ))
      rstk[i_sw, bpos] <- rsw
    }
    # between rows
    bmat[i_mub, ] <- -2 * (ob %*% k_q)
    bmat[i_sb, ] <- -(gob %*% ((k_q %x% that) + (that %x% k_q)))
    rstk[i_sb, ] <- -(gob %*% (k_q %x% k_q))

    # symmetrize the quadratic part: R <- (R + R')/2
    com_idx <- as.vector(outer(seq_len(nw), seq_len(nw),
      function(r, s) (r - 1L) * nw + s
    ))
    rstk <- (rstk + rstk[, com_idx, drop = FALSE]) / 2

    # Cov contribution: B C1 B' + 2 R (C1 x C1) R'
    scov <- scov + bmat %*% c1 %*% t(bmat) +
      2 * (rstk %*% (c1 %x% c1) %*% t(rstk))

    # bookkeeping for the expected score
    et_tot <- et_tot + tcrossprod(that) + k_q %*% c1 %*% t(k_q)
    that_all <- that_all + that
  }
  # part 1 has no (mu.w x sigma.w) transpose yet
  scov[i_sw, i_muw] <- t(scov[i_muw, i_sw])

  # ---- expected score (Fisher identity: == observed-data score) ----
  escore[i_muw] <- -2 * drop(ow %*% chat_all)
  escore[i_sw] <- drop(nobs * lav_mat_vech_dd(ow) -
    gow %*% as.vector(se_tot))
  escore[i_mub] <- -2 * drop(ob %*% that_all)
  escore[i_sb] <- drop(nclusters * lav_mat_vech_dd(ob) -
    gob %*% as.vector(et_tot))

  res <- list(scov = scov, escore = escore)
  if (extra) {
    res$post <- list(
      m1 = m1_list, c1 = c1_list,
      breg = breg_p, ccond = ccond_p,
      mu_y = mu_y, mu_z = mu_z, both_idx = both_idx,
      between_idx = between_idx,
      zmis_q = if (nz > 0L) zmis_q else NULL,
      zpat2j = if (nz > 0L) zpat2j else NULL
    )
  }
  res
}

# ---------------------------------------------------------------------
# observed information (on the scale of lav_model_hessian(), i.e. the
# Hessian of the objective function -logl/N) for a two-level model with
# missing data, using Louis's method; single group only
lav_mvn_cl_mi_h_louis <- function(lavmodel = NULL,
                                  lavsamplestats = NULL,
                                  lavdata = NULL,
                                  ceq_simple = FALSE,
                                  h = 1e-05) {
  stopifnot(lavdata@ngroups == 1L, lavdata@nlevels == 2L)
  lp <- lavdata@Lp[[1]]
  mp <- lavdata@Mp[[1]]
  y1 <- lavdata@X[[1]]
  y2 <- lavsamplestats@YLp[[1]][[2]]$Y2
  ntotal <- lavsamplestats@ntotal
  nclusters <- lp$nclusters[[2]]
  # note: the complete-data model has one 'within' term per level-1 unit,
  # *including* units with all level-1 variables missing; nunits may
  # therefore exceed ntotal (which excludes the empty units)
  nunits <- lp$nclusters[[1]]

  # current parameter values and implied moments
  implied <- lav_model_implied(lavmodel)
  mu_w <- implied$mean[[1]]
  sigma_w <- implied$cov[[1]]
  mu_b <- implied$mean[[2]]
  sigma_b <- implied$cov[[2]]

  # 1. one E-step: expected complete-data moments (per unit/cluster)
  engine <- lav_mvn_cl_mi_em_engine(
    y1 = y1, y2 = y2, lp = lp, mp = mp,
    min_variance = .Machine$double.eps # no clamping here
  )
  out2 <- lav_mvn_cl_implied22l(
    lp = lp, mu_w = mu_w, mu_b = mu_b,
    sigma_w = sigma_w, sigma_b = sigma_b
  )
  theta <- engine$pack(
    out2$mu.w, out2$mu.b, out2$sigma.w, out2$sigma.b,
    out2$mu.z, out2$sigma.zz, out2$sigma.yz
  )
  estat <- engine$implied(engine$step(theta))

  # convention shift: the E-step machinery imputes nu centered at the
  # *total* mean (the '2l' convention: mu.b.2l = Mu.B + Mu.W[both]),
  # while the complete-data score below is written in the 'implied'
  # convention (d ~ N(Mu.W, Sigma.W), u ~ N(Mu.B, Sigma.B), with Mu.W
  # possibly nonzero for both-variables, e.g. with within-level
  # covariates); the score *vector* is identical under both conventions
  # (the residuals coincide), but the expected means differ by the fixed
  # shift s0 = Mu.W[both] (at the current parameter values):
  #   E[d.impl] = E[d.2l] + s0,   E[u.impl] = E[u.2l] - s0
  # (the central second moments are shift-invariant)
  both_idx <- lp$both.idx[[2]]
  nb <- length(both_idx)
  s0 <- drop(mu_w)[both_idx]
  estat_mu_w <- drop(estat$Mu.W)
  estat_mu_w[both_idx] <- estat_mu_w[both_idx] + s0
  estat_mu_b <- drop(estat$Mu.B)
  estat_mu_b[seq_len(nb)] <- estat_mu_b[seq_len(nb)] - s0

  # 2. Cov[ S_complete | v.obs ] in moment space (-2 logl scale)
  scov <- lav_mvn_cl_mi_scov_louis(
    y1 = y1, y2 = y2, lp = lp, mp = mp,
    mu_w = mu_w, sigma_w = sigma_w, mu_b = mu_b, sigma_b = sigma_b
  )$scov

  # 3. term A: jacobian of the complete-data gradient (in theta space),
  #    with the expected sufficient statistics held fixed; this is the
  #    gradient of the two-group model fitted to the expected moments,
  #    on the same scale as lav_model_grad() (objective = -logl/N);
  #    with simple equality constraints, we (like lav_model_delta) work
  #    in the unconstrained space, and reduce at the very end
  if (lavmodel@ceq.simple.only) {
    npar <- lavmodel@nx.unco
    type_glist <- "unco"
  } else {
    npar <- lavmodel@nx.free
    type_glist <- "free"
  }
  cgrad <- function(x) {
    glist2 <- lav_model_x2glist(lavmodel, type = type_glist, x = x)
    implied_x <- lav_model_implied(lavmodel, glist = glist2)
    delta_x <- lav_model_delta(lavmodel, glist = glist2)[[1]]
    # complete-data score (-2 logl scale) at the expected moments
    ow_x <- solve.default(implied_x$cov[[1]])
    resw <- estat_mu_w - drop(implied_x$mean[[1]])
    m_w <- nunits * (estat$Sigma.W + tcrossprod(resw))
    s_muw <- -2 * nunits * drop(ow_x %*% resw)
    s_sw <- lav_mat_vech_dd(nunits * ow_x - ow_x %*% m_w %*% ow_x)
    ob_x <- solve.default(implied_x$cov[[2]])
    resb <- estat_mu_b - drop(implied_x$mean[[2]])
    m_b <- nclusters * (estat$Sigma.B + tcrossprod(resb))
    s_mub <- -2 * nclusters * drop(ob_x %*% resb)
    s_sb <- lav_mat_vech_dd(nclusters * ob_x - ob_x %*% m_b %*% ob_x)
    drop(c(s_muw, s_sw, s_mub, s_sb) %*% delta_x) / (2 * ntotal)
  }
  x <- lav_model_get_parameters(lavmodel = lavmodel)
  if (lavmodel@ceq.simple.only) {
    # unpack to the unconstrained space
    x <- drop(x %*% t(lavmodel@ceq.simple.K))
  }
  term_a <- matrix(0, npar, npar)
  for (jj in seq_len(npar)) {
    x_r <- x_l <- x
    x_r[jj] <- x[jj] + h
    x_l[jj] <- x[jj] - h
    term_a[, jj] <- (cgrad(x_r) - cgrad(x_l)) / (2 * h)
  }
  term_a <- (term_a + t(term_a)) / 2

  # 4. term B: Delta' Cov Delta, on the objective scale;
  #    obj = (-2 logl)/(2N) and -2 H_O = E[-2 H_C] - Cov(S2)/2, so
  #    H_obj = term_a - Delta' Cov Delta / (4N)
  #    (lav_model_delta() columns are in the unconstrained space when
  #    ceq.simple.only, matching term_a)
  delta <- lav_model_delta(lavmodel = lavmodel)[[1]]
  term_b <- crossprod(delta, scov %*% delta) / (4 * ntotal)

  hessian <- term_a - term_b
  if (ceq_simple && lavmodel@ceq.simple.only) {
    # reduce to the free parameters
    k_mat <- lavmodel@ceq.simple.K
    hessian <- crossprod(k_mat, hessian %*% k_mat)
  }
  hessian
}
