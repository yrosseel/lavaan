# various ways to compute a (scaled) difference chi-square test statistic

# - 0.6-13: fix multiple-group UG^2 bug in Satorra.2000 (reported by
#           Gronneberg, Foldnes and Moss) when Satterthwaite = TRUE and
#           ngroups > 1L (use old_approach = TRUE to get the old result)

# All equality-constraint reductions below use lav_con_eq_basis(): one
# orthonormal basis (eq.constraints.K, orthonormalized ceq.simple.K, or the
# ceq.JAC null space when equality constraints coexist with inequality
# constraints/bounds -- where the @eq.constraints/@ceq.simple.only packing
# flags are both FALSE), so K / t(K) round-trips are consistent everywhere.

lav_test_diff_satorra2000 <- function(m1, m0, h1 = TRUE, a_method = "delta",
                                      m_a = NULL,
                                      satterthwaite = FALSE,
                                      scaled_shifted = FALSE,
                                      old_approach = FALSE) {
  if (scaled_shifted) {
    satterthwaite <- TRUE
  }

  # extract information from m1 and m2
  t1 <- m1@test[[1]]$stat
  r1 <- m1@test[[1]]$df

  t0 <- m0@test[[1]]$stat
  r0 <- m0@test[[1]]$df

  # m = difference between the df's
  m <- r0 - r1

  # check for identical df setting
  if (m == 0L) {
    return(list(
      t_delta = (t0 - t1), scaling.factor = as.numeric(NA),
      df_delta = m, a = as.numeric(NA), b = as.numeric(NA)
    ))
  }

  # check for (near) identical test statistics (despite m > 0)
  if (abs(t1 - t0) < sqrt(.Machine$double.eps)) {
    lav_msg_warn(gettext("the test statistic of the restricted model is (nearly)
                         identical to the test statistic of the full model;
                         check your models."))
  }


  # bail out here, if m == 0 (but we should catch this earlier)
  # if(m < 1L) {
  #    txt <- paste("Can not compute (scaled) difference test when ",
  #                 "the degrees of freedom (df) are the same for both ",
  #                 "models:\n",
  #                 "Df model 1 = ", r1, ", and Df model 2 = ", r0, "\n",
  #                 sep = "")
  #            stop(lav_txt2message(txt, header = "lavaan ERROR:"))
  # }

  m_gamma <- lavTech(m1, "Gamma") # the same for m1 and m0
  # check for NULL
  if (is.null(m_gamma)) {
    lav_msg_stop(gettext(
      "can not compute Gamma matrix; perhaps missing = \"ml\"?"))
  }

  if (h1) {
    wls_v <- lavTech(m1, "WLS.V")
    m_pi <- lav_model_delta(m1@Model)
    p <- lavTech(m1, "information")
    # needed? (yes, if h1 already has eq constraints)
    p_inv <- lav_model_info_augment_invert(m1@Model,
      information = p,
      inverted = TRUE
    )
    # compute 'A' matrix
    # NOTE: order of parameters may change between h1 and H0, so be
    # careful!
    if (is.null(m_a)) {
      m_a <- lav_test_diff_a(m1, m0, method = a_method, reference = "H1")
      # take into account equality constraints m1
      if (a_method == "delta") {
        eq_basis_m1 <- lav_con_eq_basis(m1@Model)
        if (!is.null(eq_basis_m1)) {
          m_a <- m_a %*% t(eq_basis_m1)
        }
      }
      if (lav_debug()) print(m_a)
    }
  } else {
    lav_msg_stop(gettext("not ready yet"))

    wls_v <- lavTech(m0, "WLS.V")
    m_pi <- lav_model_delta(m0@Model)
    p <- lavTech(m0, "information")
    # needed?
    p_inv <- lav_model_info_augment_invert(m0@Model,
      information = p,
      inverted = TRUE
    )

    # compute 'A' matrix
    # NOTE: order of parameters may change between H1 and H0, so be
    # careful!
    if (is.null(m_a)) {
      # m1, m0 OR m0, m1 (works for delta, but not for exact)
      m_a <- lav_test_diff_a(m1, m0, method = a_method, reference = "H0")
      # take into account equality constraints m1
      eq_basis_m0 <- lav_con_eq_basis(m0@Model)
      if (!is.null(eq_basis_m0)) {
        m_a <- m_a %*% t(eq_basis_m0)
      }
      if (lav_debug()) print(m_a)
    }
  }

  # compute tr UG per group
  ngroups <- m1@SampleStats@ngroups

  # safety check: m_a %*% p_inv %*% t(m_a) should NOT contain all-zero
  # rows/columns
  # FIXME: is this really needed? As we use ginv later on
  apa <- m_a %*% p_inv %*% t(m_a)
  c_sums <- colSums(apa)
  r_sums <- rowSums(apa)
  empty_idx <- which(abs(c_sums) < .Machine$double.eps^0.5 &
    abs(r_sums) < .Machine$double.eps^0.5)
  if (length(empty_idx) > 0) {
    m_a <- m_a[-empty_idx, , drop = FALSE]
  }
  if (nrow(m_a) == 0L) {
    # oops... abort!
    return(list(
      t_delta = (t0 - t1), scaling.factor = as.numeric(NA),
      df_delta = m, a = as.numeric(NA), b = as.numeric(NA)
    ))
  }

  # paapaap
  paapaap <- p_inv %*% t(m_a) %*% MASS::ginv(m_a %*%
                               p_inv %*% t(m_a)) %*% m_a %*% p_inv

  # compute scaling factor
  fg <- unlist(m1@SampleStats@nobs) / m1@SampleStats@ntotal

  # both traces live in the npar x npar space: with M = paapaap
  # (symmetric) and K_g = Pi_g' V_g Gamma_g V_g Pi_g, we have
  # U = V Pi M Pi' V, so that (Satorra 2000, eq. 23)
  #   tr(U Gamma)     = sum_g fg_g tr(M K_g)
  #   tr((U Gamma)^2) = tr((M K)^2), K = sum_g fg_g K_g
  # this replaces the stacked pstar x pstar computation; the fg * Gamma_g
  # with UNweighted V blocks convention is trace-equivalent to the
  # Gamma_g / fg + fg-weighted-V convention used elsewhere (see the
  # SCALING CONVENTIONS note in lav_samplestats_gamma.R)
  k_group <- vector("list", length = ngroups)
  trace_ugamma_group <- numeric(ngroups)
  for (g in 1:ngroups) {
    vp <- wls_v[[g]] %*% m_pi[[g]] # pstar x npar
    k_group[[g]] <- crossprod(vp, m_gamma[[g]] %*% vp)
    trace_ugamma_group[g] <- sum(paapaap * k_group[[g]])
  }
  trace_ugamma <- sum(fg * trace_ugamma_group)

  trace_ugamma2 <- as.numeric(NA)
  if (satterthwaite) {
    if (old_approach) {
      # this is what we did <0.6-13: also the second trace per group
      trace_ugamma2_group <- numeric(ngroups)
      for (g in 1:ngroups) {
        mk <- paapaap %*% k_group[[g]]
        trace_ugamma2_group[g] <- sum(mk * t(mk))
      }
      trace_ugamma2 <- sum(fg * trace_ugamma2_group)
    } else {
      # global approach (not group-specific)
      k_all <- fg[1] * k_group[[1]]
      for (g in seq_len(ngroups - 1L) + 1L) {
        k_all <- k_all + fg[g] * k_group[[g]]
      }
      mk <- paapaap %*% k_all
      trace_ugamma2 <- sum(mk * t(mk))
    }
  }

  if (satterthwaite && !scaled_shifted) {
    cd <- trace_ugamma2 / trace_ugamma
    df_delta <- trace_ugamma^2 / trace_ugamma2
    t_delta <- (t0 - t1) / cd
    a <- as.numeric(NA)
    b <- as.numeric(NA)
  } else if (satterthwaite && scaled_shifted) {
    a <- sqrt(m / trace_ugamma2)
    # b <- m - sqrt(m * trace_ugamma^2 / trace_ugamma2)
    b <- m - a * trace_ugamma
    df_delta <- m
    t_delta <- (t0 - t1) * a + b
    cd <- as.numeric(NA)
  } else {
    cd <- 1 / m * trace_ugamma
    df_delta <- m
    t_delta <- (t0 - t1) / cd
    a <- as.numeric(NA)
    b <- as.numeric(NA)
  }

  list(
    t_delta = t_delta, scaling.factor = cd, df_delta = df_delta,
    trace.ugamma = trace_ugamma, trace.ugamma2 = trace_ugamma2,
    a = a, b = b
  )
}

lav_test_diff_sb2001 <- function(m1, m0, test = 2) {
  # extract information from m1 and m2
  t1 <- m1@test[[1]]$stat
  r1 <- m1@test[[1]]$df
  c1 <- m1@test[[test]]$scaling.factor

  ## check for situations when scaling.factor would be NA
  if (r1 == 0) {
    ## saturated model
    c1 <- 1 # canceled out by 0 when calculating "cd"

  } else if (r1 > 0 && isTRUE(all.equal(t1, 0))) {
    ## perfect fit
    c1 <- 0 # cancels out r1 when calculating "cd"
  }

  t0 <- m0@test[[1]]$stat
  r0 <- m0@test[[1]]$df
  c0 <- m0@test[[test]]$scaling.factor

  # m = difference between the df's
  m <- r0 - r1

  # check for identical df setting
  if (m == 0L) {
    return(list(
      T.delta = (t0 - t1), scaling.factor = as.numeric(NA),
      df.delta = m
    ))
  }

  # compute c_d
  cd <- (r0 * c0 - r1 * c1) / m

  # warn if cd is negative
  if (cd < 0) {
    lav_msg_warn(gettext("scaling factor is negative"))
    cd <- as.numeric(NA)
  }

  # compute scaled difference test
  t_delta <- (t0 - t1) / cd

  list(T.delta = t_delta, scaling.factor = cd, df.delta = m)
}

lav_test_diff_sb2010 <- function(m1, m0, test = 2,
                                             h1 = FALSE) {
  ### FIXME: check if models are nested at the parameter level!!!

  # extract information from m1 and m2
  t1 <- m1@test[[1]]$stat
  r1 <- m1@test[[1]]$df
  # c1 <- m1@test[[test]]$scaling.factor
  # if (r1 == 0) { # saturated model
  #   c1 <- 1
  # }

  t0 <- m0@test[[1]]$stat
  r0 <- m0@test[[1]]$df
  c0 <- m0@test[[test]]$scaling.factor
  if (r0 == 0) { # should never happen
    c0 <- 1
  }

  # m = difference between the df's
  m <- r0 - r1

  # check for identical df setting
  if (m == 0L) {
    return(list(
      T.delta = (t0 - t1), scaling.factor = as.numeric(NA),
      df.delta = m
    ))
  }

  # generate `M10' model
  if (h1) {
    # M0 with M1 parameters
    m01 <- lav_test_diff_m10(m0, m1, test = TRUE)
    c01 <- m01@test[[test]]$scaling.factor

    # check if vcov is positive definite (new in 0.6)
    # if not, we may get negative values
    eigvals <- eigen(lavTech(m01, "information"),
      symmetric = TRUE, only.values = TRUE
    )$values
    if (any(eigvals < -1 * .Machine$double.eps^(3 / 4))) {
      lav_msg_warn(gettext(
        "information matrix of the M01 model is not positive definite."
      ))
      # "                  As a result, the scale-factor can not be computed.")
      # cd <- as.numeric(NA)
    } # else {
    # compute c_d
    # cd.01 <- (r0 * c01 - r1 * c0) / m ???
    cd <- (r0 * c0 - r1 * c01) / m
    # }
  } else {
    # M1 with M0 parameters (as in Satorra & Bentler 2010)
    m10 <- lav_test_diff_m10(m1, m0, test = TRUE)
    c10 <- m10@test[[test]]$scaling.factor

    # check if vcov is positive definite (new in 0.6)
    # if not, we may get negative values
    eigvals <- eigen(lavTech(m10, "information"),
      symmetric = TRUE, only.values = TRUE
    )$values
    if (any(eigvals < -1 * .Machine$double.eps^(3 / 4))) {
      lav_msg_warn(gettext(
        "information matrix of the M10 model is not positive definite."
      ))
      # "                  As a result, the scale-factor can not be computed.")
      # cd <- as.numeric(NA)
    } # else {
    # compute c_d
    cd <- (r0 * c0 - r1 * c10) / m
    # }
  }

  # compute scaled difference test
  t_delta <- (t0 - t1) / cd

  list(
    T.delta = t_delta, scaling.factor = cd, df.delta = m,
    T.delta.unscaled = (t0 - t1)
  )
}

# create a new model 'm10', where we use model 'm1', but we
# inject it with the values of 'm0'
lav_test_diff_m10 <- function(m1, m0, test = FALSE) {
  # switch of verbose/se/test
  options_1 <- m1@Options
  # switch of optim.gradient check
  options_1$check.gradient <- FALSE

  # should we compute se/test statistics?
  if (!test) {
    options_1$se <- "none"
    options_1$test <- "none"
  }

  pt_m0 <- lav_pt_set_cache(m0@ParTable, m0@pta)
  pt_m1 <- lav_pt_set_cache(m1@ParTable, m1@pta)

  # `extend' PT.M1 partable to include all `fixed-to-zero parameters'
  pt_m1_full <- lav_pt_full(
    partable = pt_m1,
    free = TRUE, start = TRUE
  )
  pt_m1_extended <- lav_pt_merge(pt_m1, pt_m1_full,
    remove_duplicated = TRUE, warn = FALSE
  )

  # remove most columns
  pt_m1_extended$start <- NULL # new in 0.6-4! (otherwise, they are used)
  pt_m1_extended$est <- NULL
  pt_m1_extended$se <- NULL

  # in addition, use 'NA' for free parameters in ustart column
  free_par_idx <- which(pt_m1_extended$free > 0L)
  pt_m1_extended$ustart[free_par_idx] <- as.numeric(NA)

  # `extend' PT.M0 partable to include all `fixed-to-zero parameters'
  pt_m0_full <- lav_pt_full(
    partable = pt_m0,
    free = TRUE, start = TRUE
  )
  pt_m0_extended <- lav_pt_merge(pt_m0, pt_m0_full,
    remove_duplicated = TRUE, warn = FALSE
  )
  # remove most columns, but not 'est'
  pt_m0_extended$ustart <- NULL
  pt_m0_extended$start <- NULL
  pt_m0_extended$se <- NULL


  # FIXME:
  # - check if H0 does not contain additional parameters...

  options_1$optim.method <- "none"
  options_1$optim.force.converged <- TRUE
  options_1$baseline <- FALSE
  options_1$h1 <- TRUE # needed after all (yuan.benter.mplus)
  options_1$start <- pt_m0_extended # new in 0.6!
  m10 <- lavaan(
    model = pt_m1_extended,
    slot_options = options_1,
    slot_sample_stats = m1@SampleStats,
    slot_data = m1@Data,
    slot_cache = m1@Cache,
    verbose = FALSE
  )

  m10
}

# compute the `A' matrix: the jacobian of the constraint function a(\delta)
# (see Satorra 2000)
#
#
#
lav_test_diff_a <- function(m1, m0, method = "delta", reference = "H1") {
  # FIXME!!!!

  if (method == "exact") {
    if (reference == "H1") {
      af <- lav_test_diff_af_h1(m1 = m1, m0 = m0)
      xx <- m1@optim$x
    } else { # evaluate under H0
      lav_msg_stop(gettext("not ready yet"))
      # af <- .test_compute_partable_A_diff_h0(m1 = m1, m0 = m0)
      xx <- m0@optim$x
    }
    m_a <- try(lav_func_jacobian_complex(func = af, x = xx), silent = TRUE)
    if (inherits(m_a, "try-error")) {
      m_a <- lav_func_jacobian_simple(func = af, x = xx)
    }
  } else if (method == "delta") {
    # use a numeric approximation of `A'
    delta1_list <- lav_model_delta(m1@Model)
    delta0_list <- lav_model_delta(m0@Model)
    delta1 <- do.call(rbind, delta1_list)
    delta0 <- do.call(rbind, delta0_list)

    # take into account equality constraints m0
    # note: delta has one column per non-collapsed free parameter; the basis
    # K maps it to the constrained space, so we post-multiply by K (NOT t(K))
    eq_basis_m0 <- lav_con_eq_basis(m0@Model)
    if (!is.null(eq_basis_m0)) {
      delta0 <- delta0 %*% eq_basis_m0
    }

    # take into account equality constraints m1
    eq_basis_m1 <- lav_con_eq_basis(m1@Model)
    if (!is.null(eq_basis_m1)) {
      delta1 <- delta1 %*% eq_basis_m1
    }

    # H <- solve(t(Delta1) %*% Delta1) %*% t(Delta1) %*% Delta0
    h <- MASS::ginv(delta1) %*% delta0
    m_a <- t(lav_mat_ortho_complement(h))
  }

  m_a
}


# for each parameter in H1 (m1), see if we have somehow constrained
# this parameter under H0 (m0)
#
# since we work 'under H0', we need to use the labels/constraints/def
# as they appear in H0. Unfortunately, the order of the parameters, and
# even the (p)labels may be different in the two models...
#
# Therefore, we will attempt to:
#   - change the 'order' of the 'free' column in m0, so that they map to
#     to the 'x' that we will provide from H1
#   - the plabels used in "==" constraints must be renamed, if necessary
#
lav_test_diff_af_h1 <- function(m1, m0) {
  pt_m0 <- lav_pt_set_cache(parTable(m0), m0@pta)
  pt_m1 <- lav_pt_set_cache(parTable(m1), m1@pta)

  # select .p*. parameters only
  m0_p_idx <- which(grepl("\\.p", pt_m0$plabel))
  np0 <- length(m0_p_idx)
  m1_p_idx <- which(grepl("\\.p", pt_m1$plabel))
  np1 <- length(m1_p_idx)

  # check if parameter space is the same
  if (np0 != np1) {
    lav_msg_stop(gettext(
      "unconstrained parameter set is not the same in m0 and m1"))
  }

  # split partable in 'parameter' and 'constraints' section
  pt_m0_part1 <- pt_m0[m0_p_idx, ]
  pt_m0_part2 <- pt_m0[-m0_p_idx, ]

  pt_m1_part1 <- pt_m1[m1_p_idx, ]
  pt_m1_part2 <- pt_m1[-m1_p_idx, ]

  # figure out relationship between m0 and m1
  p1_id <- lav_pt_map_id_p1_in_p2(pt_m0_part1, pt_m1_part1)
  p0_free_idx <- which(pt_m0_part1$free > 0)

  # change 'free' order in m0
  # NOTE: this only works all the free parameters in h0 are also free
  # in h1 (and if not, they will become fixed in h0)
  pt_m0_part1$free[p0_free_idx] <-
    pt_m1_part1$free[pt_m0_part1$id[p1_id][p0_free_idx]]

  # paste back
  pt_m0 <- rbind(pt_m0_part1, pt_m0_part2)
  pt_m1 <- rbind(pt_m1_part1, pt_m1_part2)

  # `extend' PT.M1 partable to include all `fixed-to-zero parameters'
  pt_m1_full <- lav_pt_full(
    partable = pt_m1,
    free = TRUE, start = TRUE
  )
  pt_m1_extended <- lav_pt_merge(pt_m1, pt_m1_full,
    remove_duplicated = TRUE, warn = FALSE
  )

  # `extend' PT.M0 partable to include all `fixed-to-zero parameters'
  pt_m0_full <- lav_pt_full(
    partable = pt_m0,
    free = TRUE, start = TRUE
  )
  pt_m0_extended <- lav_pt_merge(pt_m0, pt_m0_full,
    remove_duplicated = TRUE, warn = FALSE
  )

  p1 <- pt_m1_extended
  np1 <- length(p1$lhs)
  p0 <- pt_m0_extended
  np0 <- length(p0$lhs)

  con_function <- function() NULL
  formals(con_function) <- alist(.x. = , ... = )
  body_txt <- paste("{\nout <- numeric(0L)\n", sep = "")


  # first handle def + == constraints
  # but FIRST, remove == constraints that also appear in H1!!!

  # remove equivalent eq constraints from p0
  p0_1 <- p0

  p0_eq_idx <- which(p0$op == "==")
  p1_eq_idx <- which(p1$op == "==")
  p0_remove_idx <- integer(0L)
  if (length(p0_eq_idx) > 0L) {
    for (i in seq_along(p0_eq_idx)) {
      # e0 in p0
      e0 <- p0_eq_idx[i]
      lhs <- p0$lhs[e0]
      rhs <- p0$rhs[e0]

      # do we have an equivalent constraint in H1?
      # NOTE!! the (p)labels may differ

      # SO, we will use an 'empirical' approach: if we fill in (random)
      # values, and work out the constraint, do we get identical values?
      # if yes, constraint is equivalent, and we should NOT add it here

      if (length(p1_eq_idx) > 0) {
        # generate random parameter values
        xx1 <- rnorm(length(m1_p_idx))
        xx0 <- xx1[p1_id]

        con_h0_value <- m0@Model@ceq.function(xx0)[i]
        con_h1_values <- m1@Model@ceq.function(xx1)

        if (con_h0_value %in% con_h1_values) {
          p0_remove_idx <- c(p0_remove_idx, e0)
        }
      }
    }
  }
  if (length(p0_remove_idx) > 0L) {
    p0_1 <- p0_1[-p0_remove_idx, ]
  }

  # only for the UNIQUE equality constraints in H0, generate syntax
  defcon_txt <- lav_pt_con_ceq(p0_1, txt_only = TRUE)
  body_txt <- paste(body_txt, defcon_txt, "\n", sep = "")


  # for each parameter in p1, we 'check' is it is fixed to a constant in p0
  ncon <- length(which(p0_1$op == "=="))
  for (i in seq_len(np1)) {
    # p in p1
    lhs <- p1$lhs[i]
    op <- p1$op[i]
    rhs <- p1$rhs[i]
    group <- p1$group[i]

    # ignore '==', '<', '>' and ':=' for now
    if (op == "==" || op == ">" || op == "<" || op == ":=") next

    # search for corresponding parameter in p0
    p0_idx <- which(p0$lhs == lhs & p0$op == op & p0$rhs == rhs &
      p0$group == group)
    if (length(p0_idx) == 0L) {
      lav_msg_stop(
        gettextf("parameter in H1 not found in H0: %s",
        paste(lhs, op, rhs, "(group = ", group, ")", sep = " ")
      ))
    }

    # 4 possibilities: p is free/fixed in p1, p is free/fixed in p0
    if (p1$free[i] == 0L) {
      if (p0$free[p0_idx] == 0L) {
        # match, nothing to do
      } else {
        lav_msg_warn(
          gettextf("fixed parameter in H1 is free in H0: %s",
          paste("\"", lhs, " ", op, " ", rhs,
            "\" (group = ", group, ")",
            sep = ""
          )
        ))
      }
    } else {
      if (p0$free[p0_idx] == 0L) {
        # match, this is a constrained parameter in H0
        ncon <- ncon + 1L
        body_txt <- paste(body_txt,
          "out[", ncon, "] = .x.[", p1$free[i], "] - ",
          p0$ustart[p0_idx], "\n",
          sep = ""
        )
        next
      } else {
        # match, nothing to do
      }
    }
  }


  # wrap function
  body_txt <- paste(body_txt, "return(out)\n}\n", sep = "")
  body(con_function) <- parse(file = "", text = body_txt)

  con_function
}
