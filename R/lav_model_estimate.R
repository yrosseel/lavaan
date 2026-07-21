# overview NLMINB (default) versus CONSTR (=constrained optimization)

#                         | cin.simple | nonlinear | no cin
# ----------------------------------------------------------
# eq.constraints (linear) | CONSTR     | CONSTR    | NLMINB
# ceq.nonlinear           | CONSTR     | CONSTR    | CONSTR
# ceq.simple              | NLMINB     | CONSTR    | NLMINB
# no ceq                  | NLMINB     | CONSTR    | NLMINB



# model estimation
lav_model_est <- function(lavmodel = NULL,
                               lavpartable = NULL, # for parscale = "stand"
                               lavh1 = NULL, # for multilevel + parsc
                               lavsamplestats = NULL,
                               lavdata = NULL,
                               lavoptions = NULL,
                               lavcache = list(),
                               start = "model",
                               do_fit = TRUE) {
  lavpartable <- lav_pt_set_cache(lavpartable)
  estimator <- lavoptions$estimator
  verbose <- lav_verbose()
  debug <- lav_debug()
  ngroups <- lavsamplestats@ngroups

  if (lavsamplestats@missing.flag || estimator == "PML" ||
      lavdata@nlevels > 1L) {
    group_weight <- FALSE
  } else {
    group_weight <- TRUE
  }

  # backwards compatibility < 0.6-11
  if (is.null(lavoptions$optim.partrace)) {
    lavoptions$optim.partrace <- FALSE
  }

  if (lavoptions$optim.partrace) {
    # fx + parameter values
    penv <- new.env()
    penv$PARTRACE <- matrix(NA, nrow = 0, ncol = lavmodel@nx.free + 1L)
  }

  # starting values (ignoring equality constraints)
  x_unpack <- lav_model_get_parameters(lavmodel)

  # override? use simple instead? (new in 0.6-7)
  if (start == "simple") {
    start_1 <- numeric(length(lavpartable$lhs))
    # set loadings to 0.7
    loadings_idx <- which(lavpartable$free > 0L &
      lavpartable$op == "=~")
    if (length(loadings_idx) > 0L) {
      start_1[loadings_idx] <- 0.7
    }
    # set (only) variances to 1
    var_idx <- which(lavpartable$free > 0L &
      lavpartable$op == "~~" &
      lavpartable$lhs == lavpartable$rhs)
    if (length(var_idx) > 0L) {
      start_1[var_idx] <- 1
    }

    if (lavmodel@ceq.simple.only) {
      x_unpack <- start_1[lavpartable$free > 0L &
        !duplicated(lavpartable$free)]
    } else {
      x_unpack <- start_1[lavpartable$free > 0L]
    }

    # override? use random starting values instead? (new in 0.6-18)
  } else if (start == "random") {
    start_1 <- lav_pt_random(
      lavpartable = lavpartable,
      # needed if we still need to compute bounds:
      lavh1 = lavh1,
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      lavoptions = lavoptions
    )

    if (lavmodel@ceq.simple.only) {
      x_unpack <- start_1[lavpartable$free > 0L &
        !duplicated(lavpartable$free)]
    } else {
      x_unpack <- start_1[lavpartable$free > 0L]
    }
  }

  # new in 0.7-1: multilevel models where a whole block (level) is
  # saturated: the h1 estimates for this block are already available;
  # temporarily fix these parameters at their h1 values, so they do not
  # enter the optimization again as free parameters (this often helps
  # convergence); after optimization, they are free parameters again
  h1_sat_idx <- integer(0L)
  h1_sat_ok <- lavoptions$optim.fix.saturated
  if (is.null(h1_sat_ok)) { # backwards compatibility
    h1_sat_ok <- TRUE
  }
  if (h1_sat_ok && lavdata@nlevels > 1L && estimator == "ML" &&
      # only if the optimizer honours box constraints (NLMINB family)
      length(lavmodel@ceq.nonlinear.idx) == 0L &&
      (lavmodel@cin.simple.only ||
       (length(lavmodel@cin.linear.idx) == 0L &&
        length(lavmodel@cin.nonlinear.idx) == 0L)) &&
      (is.null(lavoptions$optim.method) ||
       toupper(lavoptions$optim.method) %in%
         c("NLMINB", "NLMINB0", "NLMINB1", "L.BFGS.B"))) {
    tmp <- lav_model_est_h1_saturated(
      lavmodel = lavmodel, lavpartable = lavpartable,
      lavdata = lavdata, lavh1 = lavh1
    )
    if (length(tmp$x.idx) > 0L) {
      h1_sat_idx <- tmp$x.idx
      x_unpack[h1_sat_idx] <- tmp$value
    }
  }

  # 1. parameter scaling (to handle data scaling, not parameter scaling)
  parscale <- rep(1.0, length(x_unpack))

  # for < 0.6 compatibility
  if (is.null(lavoptions$optim.parscale)) {
    lavoptions$optim.parscale <- "none"
  }

  if (lavoptions$optim.parscale == "none") {
    # do nothing, but still set SCALE, as before


    # 0.6-17:
    # only temporarily: 'keep' this mistake, and change it later:
    # (note the "standarized")
    # we only do this to avoid breaking a test in semlbci
    # } else if(lavoptions$optim.parscale %in% c("stand", "st", "standardize",
    #                                           "standarized", "stand.all")) {

    # this is what it should be:
  } else if (lavoptions$optim.parscale %in% c(
    "stand", "st", "standardize",
    "standardized", "stand.all"
  )) {
    # rescale parameters as if the data was standardized
    # new in 0.6-2
    #
    # FIXME: this works well, as long as the variances of the
    #        latent variables (which we do not know) are more or less
    #        equal to 1.0 (eg std.lv = TRUE)
    #
    #        Once we have better estimates of those variances, we could
    #        use them to set the scale
    #

    if (lavdata@nlevels > 1L) {
      if (length(lavh1) > 0L) {
        ov_var <- lapply(lavh1$implied$cov, diag)
      } else {
        ov_var <- lapply(
          do.call(c, lapply(lavdata@Lp, "[[", "ov.idx")),
          function(x) rep(1, length(x))
        )
      }
    } else {
      if (lavoptions$conditional.x) {
        ov_var <- lavsamplestats@res.var
      } else {
        ov_var <- lavsamplestats@var
      }
    }

    if (lavoptions$std.lv) {
      parscale <- lav_standardize_all(
        lavobject = NULL,
        est = rep(1, length(lavpartable$lhs)),
        est_std = rep(1, length(lavpartable$lhs)),
        cov_std = FALSE, ov_var = ov_var,
        lavmodel = lavmodel, lavpartable = lavpartable,
        cov_x = lavsamplestats@cov.x
      )
    } else {
      # needs good estimates for lv variances! under the marker
      # (fixed-1 loading) convention, a latent variable inherits the
      # SCALE of its marker indicator, so the marker's observed
      # variance is the right order of magnitude (an upper bound);
      # higher-order factors follow their marker chain; lv's without a
      # marker fall back to 1.0 (the old guess -- which used to be
      # applied to ALL lv's, making the scaling factors useless
      # whenever the marker variables were far from unit variance)
      lv_var <- vector("list", lavmodel@ngroups)
      group_values <- lav_pt_group_values(lavpartable)
      for (g in seq_len(lavmodel@ngroups)) {
        mm_in_group <- 1:lavmodel@nmat[g] + cumsum(c(0, lavmodel@nmat))[g]
        mlist <- lavmodel@GLIST[mm_in_group]
        mm_lambda <- mlist$lambda
        n_lv <- ncol(mm_lambda)
        lv_var[[g]] <- rep(1.0, n_lv)
        lv_names_g <- lav_pt_vnames(lavpartable, "lv",
                                    group = group_values[g])
        lv_names_g <- unique(unlist(lv_names_g))
        ov_names_g <- lav_pt_vnames(lavpartable, "ov",
                                    group = group_values[g])
        ov_names_g <- unique(unlist(ov_names_g))
        if (length(lv_names_g) > 0L && n_lv > 0L) {
          s2 <- setNames(rep(NA_real_, length(lv_names_g)), lv_names_g)
          # a few passes to resolve higher-order marker chains
          for (rep_i in 1:4) {
            for (l in lv_names_g) {
              if (!is.na(s2[[l]])) next
              m_idx <- which(lavpartable$op == "=~" &
                lavpartable$lhs == l &
                lavpartable$group == group_values[g] &
                lavpartable$free == 0L &
                !is.na(lavpartable$ustart) &
                lavpartable$ustart != 0)
              if (length(m_idx) == 0L) next
              mk <- lavpartable$rhs[m_idx[1]]
              if (mk %in% ov_names_g) {
                pos <- match(mk, ov_names_g)
                if (!is.na(pos) && pos <= length(ov_var[[g]])) {
                  s2[[l]] <- ov_var[[g]][pos]
                }
              } else if (mk %in% lv_names_g && !is.na(s2[[mk]])) {
                s2[[l]] <- s2[[mk]]
              }
            }
            if (!anyNA(s2)) break
          }
          s2[is.na(s2)] <- 1.0
          pos <- match(lv_names_g, colnames(mm_lambda))
          ok_pos <- which(!is.na(pos))
          if (length(ok_pos) > 0L) {
            lv_var[[g]][pos[ok_pos]] <- s2[ok_pos]
          }
        }
      }

      parscale <- lav_standardize_all(
        lavobject = NULL,
        est = rep(1, length(lavpartable$lhs)),
        # est.std = rep(1, length(lavpartable$lhs)),
        # here, we use whatever the starting values are
        # for the latent variances...
        cov_std = FALSE, ov_var = ov_var,
        lv_var = lv_var,
        lavmodel = lavmodel, lavpartable = lavpartable,
        cov_x = lavsamplestats@cov.x
      )
    }

    # note: up to 0.7-1, the variance rows took the square root of
    # their standardization factor -- an empirical damping for the
    # (then) crude lv-variance guess of 1.0; with the marker-based lv
    # scales above, the full factor is the correct one (the z-metric
    # variance then has the same order of magnitude as a standardized
    # variance), so the sqrt tweak is gone

    if (lavmodel@ceq.simple.only) {
      parscale <- parscale[lavpartable$free > 0 &
        !duplicated(lavpartable$free)]
    } else {
      parscale <- parscale[lavpartable$free > 0]
    }
  }
  # parscale should obey the equality constraints
  if (lavmodel@eq.constraints && lavoptions$optim.parscale != "none") {
    # pack
    p_pack <- as.numeric((parscale - lavmodel@eq.constraints.k0) %*%
      lavmodel@eq.constraints.K)
    # unpack
    parscale <- as.numeric(lavmodel@eq.constraints.K %*% p_pack) +
      lavmodel@eq.constraints.k0
  }
  if (debug) {
    cat("parscale = ", parscale, "\n")
  }
  z_unpack <- x_unpack * parscale

  # scale-aware repair of degenerate variance starts (only when
  # parameter scaling is active): some starting values are absolute
  # raw-metric constants (e.g., the 0.05 default for latent variances);
  # transported into the z metric they can collapse to ~0, which puts
  # the start on the boundary of the positive-definite region and can
  # derail the optimizer. give those a sane z-metric start (0.05, the
  # standardized-world default).
  if (lavoptions$optim.parscale != "none" &&
      length(z_unpack) == length(parscale)) {
    if (lavmodel@ceq.simple.only) {
      keep <- lavpartable$free > 0L & !duplicated(lavpartable$free)
    } else {
      keep <- lavpartable$free > 0L
    }
    is_var <- (lavpartable$op == "~~" &
      lavpartable$lhs == lavpartable$rhs)[keep]
    small_idx <- which(is_var & abs(z_unpack) < 0.01)
    if (length(small_idx) > 0L) {
      z_unpack[small_idx] <- 0.05
    }
  }

  # 2. pack (apply equality constraints)
  if (lavmodel@eq.constraints && ncol(lavmodel@eq.constraints.K) > 0L) {
    z_pack <- as.numeric((z_unpack - lavmodel@eq.constraints.k0) %*%
      lavmodel@eq.constraints.K)
  } else {
    z_pack <- z_unpack
  }

  # 3. transform (already constrained) variances to standard deviations?
  # TODO
  # if(lavoptions$optim.var.transform == "sqrt" &&
  #       length(lavmodel@x.free.var.idx) > 0L) {
  #    # transforming variances using atan (or another sigmoid function?)
  #    # FIXME: better approach?
  #    #start.x[lavmodel@x.free.var.idx] <-
  #    #    atan(start.x[lavmodel@x.free.var.idx])
  #    start.x[lavmodel@x.free.var.idx] <-
  #        sqrt(start.x[lavmodel@x.free.var.idx]) # assuming positive var
  # }

  # final starting values for optimizer
  start_x <- z_pack
  if (debug) {
    cat("start.x = ", start_x, "\n")
  }

  # user-specified bounds? (new in 0.6-2)
  if (is.null(lavpartable$lower)) {
    lower <- -Inf
  } else {
    if (lavmodel@ceq.simple.only) {
      free_idx <- which(lavpartable$free > 0L &
        !duplicated(lavpartable$free))
      lower <- lavpartable$lower[free_idx]
    } else if (lavmodel@eq.constraints) {
      # bounds have no effect any longer....
      # 0.6-19 -> we switch to constrained estimation
      #lav_msg_warn(gettext(
      #  "bounds have no effect in the presence of linear
      #          equality constraints"))
      lower <- -Inf
    } else {
      lower <- lavpartable$lower[lavpartable$free > 0L]
    }
  }
  if (is.null(lavpartable$upper)) {
    upper <- +Inf
  } else {
    if (lavmodel@ceq.simple.only) {
      free_idx <- which(lavpartable$free > 0L &
        !duplicated(lavpartable$free))
      upper <- lavpartable$upper[free_idx]
    } else if (lavmodel@eq.constraints) {
      # bounds have no effect any longer....
      if (is.null(lavpartable$lower)) {
        # bounds have no effect any longer....
        # 0.6-19 -> we switch to constrained estimation
        #lav_msg_warn(gettext(
        # "bounds have no effect in the presence of linear
        #          equality constraints"))
      }
      upper <- +Inf
    } else {
      upper <- lavpartable$upper[lavpartable$free > 0L]
    }
  }

  # parameter scaling: the optimizer iterates in the z = x * parscale
  # metric, so the box constraints must be transformed as well (up to
  # 0.7-1 they were applied to the z values as if they were raw
  # values, silently distorting any bounds when optim.parscale was
  # active)
  if (lavoptions$optim.parscale != "none" &&
      length(parscale) > 0L && all(is.finite(parscale)) &&
      all(parscale > 0)) {
    if (length(lower) %in% c(1L, length(parscale))) {
      lower <- lower * parscale
    }
    if (length(upper) %in% c(1L, length(parscale))) {
      upper <- upper * parscale
    }
  }

  # check for inconsistent lower/upper bounds
  # this may happen if we have equality constraints; qr() may switch
  # the sign...
  bad_idx <- which(lower > upper)
  if (length(bad_idx) > 0L) {
    # switch
    # tmp <- lower[bad.idx]
    # lower[bad.idx] <- upper[bad.idx]
    # upper[bad.idx] <- tmp
    lower[bad_idx] <- -Inf
    upper[bad_idx] <- +Inf
  }

  # new in 0.7-1: fix the parameters of saturated blocks at their h1
  # values during optimization, by forcing lower == upper == start
  if (length(h1_sat_idx) > 0L) {
    lower <- rep(lower, length.out = length(start_x))
    upper <- rep(upper, length.out = length(start_x))
    lower[h1_sat_idx] <- upper[h1_sat_idx] <- start_x[h1_sat_idx]
  }

  # function to be minimized
  objective_function <- function(x, verbose = FALSE, inf_to_max = FALSE,
                                 debug = FALSE) {
    # 3. standard deviations to variances
    # WARNING: x is still packed here!
    # if(lavoptions$optim.var.transform == "sqrt" &&
    #   length(lavmodel@x.free.var.idx) > 0L) {
    #    #x[lavmodel@x.free.var.idx] <- tan(x[lavmodel@x.free.var.idx])
    #    x.var <- x[lavmodel@x.free.var.idx]
    #    x.var.sign <- sign(x.var)
    #    x[lavmodel@x.free.var.idx] <- x.var.sign * (x.var * x.var) # square!
    # }

    # 2. unpack
    if (lavmodel@eq.constraints) {
      x <- as.numeric(lavmodel@eq.constraints.K %*% x) +
        lavmodel@eq.constraints.k0
    }

    # 1. unscale
    x <- x / parscale

    # update GLIST (change `state') and make a COPY!
    glist <- lav_model_x2glist(lavmodel, x = x)

    fx <- lav_model_objective(
      lavmodel = lavmodel,
      glist = glist,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavcache = lavcache
    )

    # only for PML: divide by N (to speed up convergence)
    if (estimator == "PML") {
      fx <- fx / lavsamplestats@ntotal
    }



    if (debug || verbose) {
      cat("  objective function  = ",
        sprintf("%18.16f", fx), "\n",
        sep = ""
      )
    }
    if (debug) {
      # cat("Current unconstrained parameter values =\n")
      # tmp.x <- lav_model_get_parameters(lavmodel, glist=GLIST, type="unco")
      # print(tmp.x); cat("\n")
      cat("Current free parameter values =\n")
      print(x)
      cat("\n")
    }

    if (lavoptions$optim.partrace) {
      penv$PARTRACE <- rbind(penv$PARTRACE, c(fx, x))
    }

    # for L-BFGS-B
    # if(infToMax && is.infinite(fx)) fx <- 1e20
    if (!is.finite(fx)) {
      fx_group <- attr(fx, "fx.group")
      fx <- 1e20
      attr(fx, "fx.group") <- fx_group # only for lav_model_fit()
    }

    fx
  }

  gradient_function <- function(x, verbose = FALSE, inf_to_max = FALSE,
                                debug = FALSE) {
    # transform variances back
    # if(lavoptions$optim.var.transform == "sqrt" &&
    #   length(lavmodel@x.free.var.idx) > 0L) {
    #    #x[lavmodel@x.free.var.idx] <- tan(x[lavmodel@x.free.var.idx])
    #    x.var <- x[lavmodel@x.free.var.idx]
    #    x.var.sign <- sign(x.var)
    #    x[lavmodel@x.free.var.idx] <- x.var.sign * (x.var * x.var) # square!
    # }

    # 2. unpack
    if (lavmodel@eq.constraints) {
      x <- as.numeric(lavmodel@eq.constraints.K %*% x) +
        lavmodel@eq.constraints.k0
    }

    # 1. unscale
    x <- x / parscale

    # update GLIST (change `state') and make a COPY!
    glist <- lav_model_x2glist(lavmodel, x = x)

    dx <- lav_model_grad(
      lavmodel = lavmodel,
      glist = glist,
      lavsamplestats = lavsamplestats,
      lavdata = lavdata,
      lavcache = lavcache,
      type = "free",
      group_weight = group_weight, ### check me!!
      ceq_simple = lavmodel@ceq.simple.only
    )

    if (debug) {
      cat("Gradient function (analytical) =\n")
      print(dx)
      cat("\n")
    }

    # 1. scale (note: divide, not multiply!)
    dx <- dx / parscale

    # 2. pack
    if (lavmodel@eq.constraints) {
      dx <- as.numeric(dx %*% lavmodel@eq.constraints.K)
    }

    # 3. transform variances back
    # if(lavoptions$optim.var.transform == "sqrt" &&
    #   length(lavmodel@x.free.var.idx) > 0L) {
    #    x.var <- x[lavmodel@x.free.var.idx] # here in 'var' metric
    #    x.var.sign <- sign(x.var)
    #    x.var <- abs(x.var)
    #    x.sd <- sqrt(x.var)
    #    dx[lavmodel@x.free.var.idx] <-
    #        ( 2 * x.var.sign * dx[lavmodel@x.free.var.idx] * x.sd )
    # }

    # only for PML: divide by N (to speed up convergence)
    if (estimator == "PML") {
      dx <- dx / lavsamplestats@ntotal
    }

    if (debug) {
      cat("Gradient function (analytical, after eq.constraints.K) =\n")
      print(dx)
      cat("\n")
    }

    dx
  }

  gradient_function_numerical <- function(x, verbose = FALSE, debug = FALSE) {
    # NOTE: no need to 'transform' anything here (var/eq)
    # this is done anyway in objective_function

    # numerical approximation using the Richardson method
    npar <- length(x)
    h <- 10e-6
    dx <- numeric(npar)

    ## FIXME: call lav_model_objective directly!!
    for (i in 1:npar) {
      x_left <- x_left2 <- x_right <- x_right2 <- x
      x_left[i] <- x[i] - h
      x_left2[i] <- x[i] - 2 * h
      x_right[i] <- x[i] + h
      x_right2[i] <- x[i] + 2 * h
      fx_left <- objective_function(x_left, verbose = FALSE, debug = FALSE)
      fx_left2 <- objective_function(x_left2, verbose = FALSE, debug = FALSE)
      fx_right <- objective_function(x_right, verbose = FALSE, debug = FALSE)
      fx_right2 <- objective_function(x_right2, verbose = FALSE, debug = FALSE)
      dx[i] <- (fx_left2 - 8 * fx_left + 8 * fx_right - fx_right2) / (12 * h)
    }

    # dx <- lavGradientC(func=objective_function, x=x)
    # does not work if pnorm is involved... (eg PML)

    if (debug) {
      cat("Gradient function (numerical) =\n")
      print(dx)
      cat("\n")
    }

    dx
  }

  gradient_function_numerical_complex <- function(x, verbose = FALSE, debug = FALSE) { # nolint
    dx <- Re(lav_func_grad_complex(
      func = objective_function, x = x,
      h = sqrt(.Machine$double.eps)
    ))
    # does not work if pnorm is involved... (eg PML)

    if (debug) {
      cat("Gradient function (numerical complex) =\n")
      print(dx)
      cat("\n")
    }

    dx
  }


  # check if the initial values produce a positive definite Sigma
  # to begin with -- but only for estimator="ML"
  if (estimator %in% c("ML", "FML", "MML")) {
    sigma_hat <- lav_model_sigma(lavmodel, extra = TRUE)
    for (g in 1:ngroups) {
      if (!attr(sigma_hat[[g]], "po")) {
        group_txt <-
          if (ngroups > 1) gettextf(" in group %s.", g) else "."
        if (debug) {
          print(sigma_hat[[g]][, ])
        }
        lav_msg_warn(gettext(
          "initial model-implied matrix (Sigma) is not positive definite;
          check your model and/or starting parameters"), group_txt)
        x <- start_x
        fx <- as.numeric(NA)
        attr(fx, "fx.group") <- rep(as.numeric(NA), ngroups)
        attr(x, "converged") <- FALSE
        attr(x, "iterations") <- 0L
        attr(x, "control") <- lavoptions$control
        attr(x, "fx") <- fx
        return(x)
      }
    }
  }


  # parameter scaling
  # FIXME: what is the best way to set the scale??
  # current strategy: if startx > 1.0, we rescale by using
  # 1/startx
  scale_1 <- rep(1.0, length(start_x))
  idx <- which(abs(start_x) > 1.0)
  if (length(idx) > 0L) {
    scale_1[idx] <- abs(1.0 / start_x[idx])
  }
  if (debug) {
    cat("SCALE = ", scale_1, "\n")
  }


  # first try: check if starting values return a finite value
  fx <- objective_function(start_x, verbose = verbose, debug = debug)
  if (!is.finite(fx)) {
    # emergency change of start.x
    if (length(h1_sat_idx) > 0L) {
      # do not touch the (fixed) h1 values of saturated blocks
      start_x[-h1_sat_idx] <- start_x[-h1_sat_idx] / 10
    } else {
      start_x <- start_x / 10
    }
  }



  # first some nelder mead steps? (default = FALSE)
  init_nelder_mead <- lavoptions$optim.init_nelder_mead

  # gradient: analytic, numerical or NULL?
  if (is.character(lavoptions$optim.gradient)) {
    if (lavoptions$optim.gradient %in% c("analytic", "analytical")) {
      gradient <- gradient_function
    } else if (lavoptions$optim.gradient %in% c("numerical", "numeric")) {
      gradient <- gradient_function_numerical
    } else if (lavoptions$optim.gradient %in% c("numeric.complex", "complex")) {
      gradient <- gradient_function_numerical_complex
    } else if (lavoptions$optim.gradient %in% c("NULL", "null")) {
      gradient <- NULL
    } else {
      lav_msg_warn(gettext("gradient should be analytic, numerical or NULL"))
    }
  } else if (is.logical(lavoptions$optim.gradient)) {
    if (lavoptions$optim.gradient) {
      gradient <- gradient_function
    } else {
      gradient <- NULL
    }
  } else if (is.null(lavoptions$optim.gradient)) {
    gradient <- gradient_function
  }


  # default optimizer
  if (length(lavmodel@ceq.nonlinear.idx) == 0L &&
      (lavmodel@cin.simple.only || (length(lavmodel@cin.linear.idx)    == 0L &&
                                    length(lavmodel@cin.nonlinear.idx) == 0L))
     ) {
    if (is.null(lavoptions$optim.method)) {
      optimizer <- "NLMINB"
      # OPTIMIZER <- "BFGS"  # slightly slower, no bounds; better scaling!
      # OPTIMIZER <- "L-BFGS-B"  # trouble with Inf values for fx!
    } else {
      optimizer <- toupper(lavoptions$optim.method)
      stopifnot(optimizer %in% c(
        "NLMINB0", "NLMINB1", "NLMINB2",
        "NLMINB", "BFGS", "L.BFGS.B", "NONE"
      ))
      if (optimizer == "NLMINB1") {
        optimizer <- "NLMINB"
      }
    }
  } else {
    if (is.null(lavoptions$optim.method)) {
      optimizer <- "NLMINB.CONSTR"
    } else {
      optimizer <- toupper(lavoptions$optim.method)
      stopifnot(optimizer %in% c("NLMINB.CONSTR", "NLMINB", "NONE"))
    }
    if (optimizer == "NLMINB") {
      optimizer <- "NLMINB.CONSTR"
    }
  }

  if (init_nelder_mead) {
    if (verbose) cat("  initial Nelder-Mead step:\n")
    # trace <- 0L
    # if (verbose) trace <- 1L
    optim_out <- optim(
      par = start_x,
      fn = objective_function,
      method = "Nelder-Mead",
      # control=list(maxit=10L,
      #             parscale=SCALE,
      #             trace=trace),
      hessian = FALSE,
      verbose = verbose, debug = debug
    )
    cat("\n")
    start_x <- optim_out$par
  }



  if (optimizer == "NLMINB0") {
    if (verbose)
      cat("  quasi-Newton steps using NLMINB0 (no analytic gradient):\n")
    # if(debug) control$trace <- 1L;
    control_nlminb <- list(
      eval.max = 20000L,
      iter.max = 10000L,
      trace = 0L,
      # abs.tol=1e-20, ### important!! fx never negative
      abs.tol = (.Machine$double.eps * 10),
      rel.tol = 1e-10,
      # step.min=2.2e-14, # in =< 0.5-12
      step.min = 1.0, # 1.0 in < 0.5-21
      step.max = 1.0,
      x.tol = 1.5e-8,
      xf.tol = 2.2e-14
    )
    control_nlminb <- modifyList(control_nlminb, lavoptions$control)
    control <- control_nlminb[c(
      "eval.max", "iter.max", "trace",
      "step.min", "step.max",
      "abs.tol", "rel.tol", "x.tol", "xf.tol"
    )]
    # cat("DEBUG: control = "); print(str(control.nlminb)); cat("\n")
    optim_out <- nlminb(
      start = start_x,
      objective = objective_function,
      gradient = NULL,
      lower = lower,
      upper = upper,
      control = control,
      scale = scale_1,
      verbose = verbose, debug = debug
    )
    if (verbose) {
      cat("  convergence status (0=ok): ", optim_out$convergence, "\n")
      cat("  nlminb message says: ", optim_out$message, "\n")
      cat("  number of iterations: ", optim_out$iterations, "\n")
      cat(
        "  number of function evaluations [objective, gradient]: ",
        optim_out$evaluations, "\n"
      )
    }

    # try again
    if (optim_out$convergence != 0L) {
      optim_out <- nlminb(
        start = start_x,
        objective = objective_function,
        gradient = NULL,
        lower = lower,
        upper = upper,
        control = control,
        scale = scale_1,
        verbose = verbose, debug = debug
      )
    }

    iterations <- optim_out$iterations
    x <- optim_out$par
    if (optim_out$convergence == 0L) {
      converged <- TRUE
    } else {
      converged <- FALSE
    }
  } else if (optimizer == "NLMINB") {
    if (verbose) cat("  quasi-Newton steps using NLMINB:\n")
    # if(debug) control$trace <- 1L;
    control_nlminb <- list(
      eval.max = 20000L,
      iter.max = 10000L,
      trace = 0L,
      # abs.tol=1e-20, ### important!! fx never negative
      abs.tol = (.Machine$double.eps * 10),
      rel.tol = 1e-10,
      # step.min=2.2e-14, # in =< 0.5-12
      step.min = 1.0, # 1.0 in < 0.5-21
      step.max = 1.0,
      x.tol = 1.5e-8,
      xf.tol = 2.2e-14
    )
    control_nlminb <- modifyList(control_nlminb, lavoptions$control)
    control <- control_nlminb[c(
      "eval.max", "iter.max", "trace",
      "step.min", "step.max",
      "abs.tol", "rel.tol", "x.tol", "xf.tol"
    )]
    # cat("DEBUG: control = "); print(str(control.nlminb)); cat("\n")
    optim_out <- nlminb(
      start = start_x,
      objective = objective_function,
      gradient = gradient,
      lower = lower,
      upper = upper,
      control = control,
      scale = scale_1,
      verbose = verbose, debug = debug
    )
    if (verbose) {
      cat("  convergence status (0=ok): ", optim_out$convergence, "\n")
      cat("  nlminb message says: ", optim_out$message, "\n")
      cat("  number of iterations: ", optim_out$iterations, "\n")
      cat(
        "  number of function evaluations [objective, gradient]: ",
        optim_out$evaluations, "\n"
      )
    }

    iterations <- optim_out$iterations
    x <- optim_out$par
    if (optim_out$convergence == 0L) {
      converged <- TRUE
    } else {
      converged <- FALSE
    }
  } else if (optimizer == "BFGS") {
    # warning: Bollen example with estimator=GLS does NOT converge!
    # (but WLS works!)
    # - BB.ML works too

    control_bfgs <- list(
      trace = 0L, fnscale = 1,
      parscale = scale_1, ## or not?
      ndeps = 1e-3,
      maxit = 10000,
      abstol = 1e-20,
      reltol = 1e-10,
      REPORT = 1L
    )
    control_bfgs <- modifyList(control_bfgs, lavoptions$control)
    control <- control_bfgs[c(
      "trace", "fnscale", "parscale", "ndeps",
      "maxit", "abstol", "reltol", "REPORT"
    )]
    # trace <- 0L; if(verbose) trace <- 1L
    optim_out <- optim(
      par = start_x,
      fn = objective_function,
      gr = gradient,
      method = "BFGS",
      control = control,
      hessian = FALSE,
      verbose = verbose, debug = debug
    )
    if (verbose) {
      cat("  convergence status (0=ok): ", optim_out$convergence, "\n")
      cat("  optim BFGS message says: ", optim_out$message, "\n")
      # cat("number of iterations: ", optim.out$iterations, "\n")
      cat(
        "  number of function evaluations [objective, gradient]: ",
        optim_out$counts, "\n"
      )
    }

    # iterations <- optim.out$iterations
    iterations <- optim_out$counts[1]
    x <- optim_out$par
    if (optim_out$convergence == 0L) {
      converged <- TRUE
    } else {
      converged <- FALSE
    }
  } else if (optimizer == "L.BFGS.B") {
    # warning, does not cope with Inf values!!

    control_lbfgsb <- list(
      trace = 0L, fnscale = 1,
      parscale = scale_1, ## or not?
      ndeps = 1e-3,
      maxit = 10000,
      REPORT = 1L,
      lmm = 5L,
      factr = 1e7,
      pgtol = 0
    )
    control_lbfgsb <- modifyList(control_lbfgsb, lavoptions$control)
    control <- control_lbfgsb[c(
      "trace", "fnscale", "parscale",
      "ndeps", "maxit", "REPORT", "lmm",
      "factr", "pgtol"
    )]
    optim_out <- optim(
      par = start_x,
      fn = objective_function,
      gr = gradient,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = control,
      hessian = FALSE,
      verbose = verbose, debug = debug,
      inf_to_max = TRUE
    )
    if (verbose) {
      cat("  convergence status (0=ok): ", optim_out$convergence, "\n")
      cat("  optim L-BFGS-B message says: ", optim_out$message, "\n")
      # cat("number of iterations: ", optim.out$iterations, "\n")
      cat(
        "  number of function evaluations [objective, gradient]: ",
        optim_out$counts, "\n"
      )
    }

    # iterations <- optim.out$iterations
    iterations <- optim_out$counts[1]
    x <- optim_out$par
    if (optim_out$convergence == 0L) {
      converged <- TRUE
    } else {
      converged <- FALSE
    }
  } else if (optimizer == "NLMINB.CONSTR") {
    ocontrol <- list(verbose = verbose)
    if (!is.null(lavoptions$control$control.outer)) {
      ocontrol <- c(lavoptions$control$control.outer, verbose = verbose)
    }
    control_nlminb <- list(
      eval.max = 20000L,
      iter.max = 10000L,
      trace = 0L,
      # abs.tol=1e-20,
      abs.tol = (.Machine$double.eps * 10),
      rel.tol = 1e-9, # 1e-10 seems 'too strict'
      step.min = 1.0, # 1.0 in < 0.5-21
      step.max = 1.0,
      x.tol = 1.5e-8,
      xf.tol = 2.2e-14
    )
    control_nlminb <- modifyList(control_nlminb, lavoptions$control)
    control <- control_nlminb[c(
      "eval.max", "iter.max", "trace",
      "abs.tol", "rel.tol"
    )]
    cin <- cin_jac <- ceq <- ceq_jac <- NULL
    if (!is.null(body(lavmodel@cin.function))) cin <- lavmodel@cin.function
    if (!is.null(body(lavmodel@cin.jacobian))) cin_jac <- lavmodel@cin.jacobian
    if (!is.null(body(lavmodel@ceq.function))) ceq <- lavmodel@ceq.function
    if (!is.null(body(lavmodel@ceq.jacobian))) ceq_jac <- lavmodel@ceq.jacobian
    trace <- FALSE
    if (verbose) trace <- TRUE
    optim_out <- nlminb_constr(
      start = start_x,
      objective = objective_function,
      gradient = gradient,
      control = control,
      scale = scale_1,
      verbose = verbose, debug = debug,
      lower = lower,
      upper = upper,
      cin = cin, cin_jac = cin_jac,
      ceq = ceq, ceq_jac = ceq_jac,
      control_outer = ocontrol
    )
    if (verbose) {
      cat("  convergence status (0=ok): ", optim_out$convergence, "\n")
      cat("  nlminb_constr message says: ", optim_out$message, "\n")
      cat("  number of outer iterations: ", optim_out$outer.iterations, "\n")
      cat("  number of inner iterations: ", optim_out$iterations, "\n")
      cat(
        "  number of function evaluations [objective, gradient]: ",
        optim_out$evaluations, "\n"
      )
    }

    iterations <- optim_out$iterations
    x <- optim_out$par
    if (optim_out$convergence == 0) {
      converged <- TRUE
    } else {
      converged <- FALSE
    }
  } else if (optimizer == "NONE") {
    x <- start_x
    iterations <- 0L
    converged <- TRUE
    control <- list()

    # if inequality constraints, add con.jac/lambda
    # needed for df!
    if (length(lavmodel@ceq.nonlinear.idx) == 0L &&
        (lavmodel@cin.simple.only ||
         (length(lavmodel@cin.linear.idx) == 0L &&
          length(lavmodel@cin.nonlinear.idx) == 0L))) {
      optim_out <- list()
    } else {
      # if inequality constraints, add con.jac/lambda
      # needed for df!

      optim_out <- list()
      if (is.null(body(lavmodel@ceq.function))) {
        ceq <- function(x, ...) {
          numeric(0)
        }
      } else {
        ceq <- lavmodel@ceq.function
      }
      if (is.null(body(lavmodel@cin.function))) {
        cin <- function(x, ...) {
          numeric(0)
        }
      } else {
        cin <- lavmodel@cin.function
      }
      ceq0 <- ceq(start_x)
      cin0 <- cin(start_x)
      con0 <- c(ceq0, cin0)
      jac <- rbind(
        numDeriv::jacobian(ceq, x = start_x),
        numDeriv::jacobian(cin, x = start_x)
      )
      nceq <- length(ceq(start_x))
      ncin <- length(cin(start_x))
      ncon <- nceq + ncin
      ceq_idx <- cin_idx <- integer(0)
      if (nceq > 0L) ceq_idx <- 1:nceq
      if (ncin > 0L) cin_idx <- nceq + 1:ncin
      cin_flag <- rep(FALSE, ncon)
      if (ncin > 0L) cin_flag[cin_idx] <- TRUE

      inactive_idx <- integer(0L)
      cin_idx <- which(cin_flag)
      if (ncin > 0L) {
        # optimizer == "NONE": no multiplier information is available for
        # the externally supplied solution, so a constraint that touches
        # the boundary is conservatively treated as binding (this only
        # feeds the df bookkeeping); no strict-complementarity test here
        slack <- 1e-05
        inactive_idx <- which(cin_flag & con0 > slack)
      }
      attr(jac, "inactive.idx") <- inactive_idx
      attr(jac, "cin.idx") <- cin_idx
      attr(jac, "ceq.idx") <- ceq_idx

      optim_out$con.jac <- jac
      optim_out$lambda <- rep(0, ncon)
    }
  }

  # new in 0.6-19
  # if NLMINB() + cin.simple.only, add con.jac and lambda to optim.out
  if (optimizer %in% c("NLMINB", "NLMINB0", "L.BFGS.B") &&
      (lavmodel@cin.simple.only || lavmodel@ceq.simple.only)) {

    if (lavmodel@cin.simple.only && nrow(lavmodel@cin.JAC) > 0L) {
      # JAC
      cin_jac_1 <- lavmodel@cin.JAC

      # lambda (post-hoc): the bound rows involve one parameter each, so
      # the rows of cin.JAC are mutually orthogonal and the least-squares
      # multiplier reduces to the plain inner product with the gradient.
      # NOTE: for ceq.simple models both the parameter vector and the
      # gradient live in the reduced (packed) space, while cin.function
      # and cin.JAC operate on the full (unco) space -- unpack BOTH
      # (evaluating cin.function on the packed x returned garbage
      # constraint values, silently corrupting the old classification)
      dx <- gradient(x)
      if (lavmodel@ceq.simple.only) {
        unpack_idx <- lavpartable$free[lavpartable$free > 0]
        x_unpack <- x[unpack_idx]
        dx_unpack <- dx[unpack_idx]
      } else if (lavmodel@eq.constraints) {
        # unreachable via the standard pipeline: eq.constraints is a
        # packing flag that is only TRUE when equality constraints are the
        # ONLY constraints, which contradicts cin.simple.only; kept as a
        # safety net
        x_unpack <- as.numeric(lavmodel@eq.constraints.K %*% x) +
          lavmodel@eq.constraints.k0
        dx_unpack <- as.numeric(lavmodel@eq.constraints.K %*% dx)
      } else {
        x_unpack <- x
        dx_unpack <- dx
      }
      con0 <- lavmodel@cin.function(x_unpack)
      cin_lambda <- drop(cin_jac_1 %*% dx_unpack)

      # strict complementarity: a bound is only ACTIVE (binding) when the
      # solution sits at the bound AND the multiplier exerts force on the
      # gradient (see lav_con_cin_inactive_idx()); an interior solution
      # that merely grazes a bound keeps its ordinary standard error
      inactive_idx <- lav_con_cin_inactive_idx(
        con0 = con0, lambda = cin_lambda, jac = cin_jac_1,
        cin_flag = rep(TRUE, nrow(cin_jac_1))
      )
      cin_lambda[inactive_idx] <- 0

      # remove all inactive rows
      #if (length(inactive.idx) > 0L) {
      #  cin.JAC <- cin.JAC[-inactive.idx, , drop = FALSE]
      #  cin.lambda <- cin.lambda[-inactive.idx]
      #  inactive.idx <- integer(0L)
      #}
    } else {
      npar <- length(lavpartable$free[lavpartable$free > 0])
      cin_jac_1 <- matrix(0, nrow = 0L, ncol = npar)
      inactive_idx <- integer(0L)
      cin_lambda <- numeric(0L)
    }

    if (lavmodel@ceq.simple.only && nrow(lavmodel@ceq.simple.K) > 0L) {
      ceq_jac_1 <- t(lav_mat_ortho_complement(lavmodel@ceq.simple.K))
      ceq_lambda <- numeric(nrow(ceq_jac_1))
    } else {
      npar <- length(lavpartable$free[lavpartable$free > 0])
      ceq_jac_1 <- matrix(0, nrow = 0L, ncol = npar)
      ceq_lambda <- numeric(0L)
    }

    # combine
    jac <- rbind(cin_jac_1, ceq_jac_1)
    attr(jac, "inactive.idx") <- inactive_idx
    attr(jac, "cin.idx") <- seq_len(nrow(cin_jac_1))
    attr(jac, "ceq.idx") <- nrow(cin_jac_1) + seq_len(nrow(ceq_jac_1))
    lambda <- c(cin_lambda, ceq_lambda)

    optim_out$con.jac <- jac
    optim_out$lambda <- lambda
  }


  fx <- objective_function(x) # to get "fx.group" attribute

  # check convergence
  warn_txt <- ""
  if (converged) {
    # check.gradient
    if (!is.null(gradient) &&
      optimizer %in% c("NLMINB", "BFGS", "L.BFGS.B")) {
      # gradient in the z (optimizer) metric: when parameter scaling is
      # active, this is the right metric for an absolute tolerance --
      # the z parameters live on a standardized-world scale, while the
      # RAW gradient components of (say) a variance of order 1e-6 are
      # amplified by the inverse scale and could never pass a fixed
      # cutoff, even at a machine-precision optimum
      dx <- gradient(x)

      if (converged && lavoptions$check.gradient &&
        any(abs(dx) > lavoptions$optim.dx.tol)) {
        # ok, identify the non-zero elements
        non_zero <- which(abs(dx) > lavoptions$optim.dx.tol)

        # which ones are 'boundary' points, defined by lower/upper?
        bound_idx <- integer(0L)
        if (!is.null(lavpartable$lower)) {
          bound_idx <- c(bound_idx, which(lower == x))
        }
        if (!is.null(lavpartable$upper)) {
          bound_idx <- c(bound_idx, which(upper == x))
        }
        # parameters fixed at their h1 values (saturated blocks) were
        # not part of the search
        if (length(h1_sat_idx) > 0L) {
          bound_idx <- c(bound_idx, h1_sat_idx)
        }
        if (length(bound_idx) > 0L) {
          non_zero <- non_zero[-which(non_zero %in% bound_idx)]
        }

        # this has many implications ... so should be careful to
        # avoid false alarm
        if (length(non_zero) > 0L) {
          converged <- FALSE
          warn_txt <- paste("the optimizer (", optimizer, ") ",
            "claimed the model converged,\n",
            "       but not all elements of the gradient are (near) zero;\n",
            "       the optimizer may not have found a local solution\n",
            "       use check.gradient = FALSE to skip this check.",
            sep = ""
          )
        }
      }
    } else {
      dx <- numeric(0L)
    }
  } else {
    dx <- numeric(0L)
    warn_txt <- "the optimizer warns that a solution has NOT been found!"
  }

  # transform back
  # 3.
  # if(lavoptions$optim.var.transform == "sqrt" &&
  #       length(lavmodel@x.free.var.idx) > 0L) {
  #    #x[lavmodel@x.free.var.idx] <- tan(x[lavmodel@x.free.var.idx])
  #    x.var <- x[lavmodel@x.free.var.idx]
  #    x.var.sign <- sign(x.var)
  #    x[lavmodel@x.free.var.idx] <- x.var.sign * (x.var * x.var) # square!
  # }

  # 2. unpack
  if (lavmodel@eq.constraints) {
    x <- as.numeric(lavmodel@eq.constraints.K %*% x) +
      lavmodel@eq.constraints.k0
  }

  # 1. unscale
  x <- x / parscale

  attr(x, "converged") <- converged
  attr(x, "start") <- start_x
  attr(x, "warn.txt") <- warn_txt
  attr(x, "iterations") <- iterations
  attr(x, "control") <- control
  attr(x, "fx") <- fx
  attr(x, "dx") <- dx
  attr(x, "parscale") <- parscale
  if (!is.null(optim_out$con.jac)) attr(x, "con.jac") <- optim_out$con.jac
  if (!is.null(optim_out$lambda)) attr(x, "con.lambda") <- optim_out$lambda
  if (lavoptions$optim.partrace) {
    attr(x, "partrace") <- penv$PARTRACE
  }

  x
}

# backwards compatibility
# estimateModel <- lav_model_estimate

# new in 0.7-1: identify free parameters that belong to a 'saturated'
# block of a multilevel model
#
# a block is considered to be saturated if it only contains (free,
# unconstrained) variances/covariances -- and means/intercepts -- of the
# observed variables of that block; in that case, the 'estimated' values of
# these parameters are already available in lavh1$implied, and they do not
# need to enter the optimization again as free parameters: they are
# (temporarily) fixed at their h1 values during optimization only; once
# model estimation is over, they are treated as free parameters again
# (eg to compute standard errors)
#
# returns a list with two elements:
#   - x.idx: positions of these parameters in the (packed) parameter vector
#   - value: the corresponding h1 estimates
lav_model_est_h1_saturated <- function(lavmodel = NULL,
                                       lavpartable = NULL,
                                       lavdata = NULL,
                                       lavh1 = NULL) {
  empty <- list(x.idx = integer(0L), value = numeric(0L))

  # multilevel only (for now), and we need the h1 estimates
  if (lavdata@nlevels == 1L || length(lavh1) == 0L ||
      is.null(lavh1$implied$cov)) {
    return(empty)
  }

  # general linear equality constraints: the optimizer operates in a
  # reduced (packed) parameter space; we cannot fix individual parameters
  if (lavmodel@eq.constraints) {
    return(empty)
  }

  nlevels <- lavdata@nlevels
  ngroups <- lavdata@ngroups

  # position of each free parameter in the (packed) parameter vector
  free_id <- lavpartable$free[lavpartable$free > 0L]
  uid <- unique(free_id) # ceq.simple: duplicated ids are packed out
  dup_id <- unique(free_id[duplicated(free_id)])

  # labels that show up in (in)equality constraints
  con_idx <- which(lavpartable$op %in% c("==", "<", ">"))
  con_ref <- unique(c(lavpartable$lhs[con_idx], lavpartable$rhs[con_idx]))

  x_idx <- integer(0L)
  value <- numeric(0L)

  for (g in seq_len(ngroups)) {
    for (l in seq_len(nlevels)) {
      b <- (g - 1L) * nlevels + l
      ovn <- lavdata@ov.names.l[[g]][[l]]
      p <- length(ovn)
      cov_h1 <- lavh1$implied$cov[[b]]
      mean_h1 <- lavh1$implied$mean[[b]]
      if (p == 0L || is.null(cov_h1) || nrow(cov_h1) != p ||
          length(mean_h1) != p) {
        next
      }

      row_idx <- which(lavpartable$block == b)
      if (length(row_idx) == 0L) {
        next
      }
      op <- lavpartable$op[row_idx]
      lhs <- lavpartable$lhs[row_idx]
      rhs <- lavpartable$rhs[row_idx]
      free <- lavpartable$free[row_idx]

      # only variances/covariances and intercepts of observed variables
      # (no latent variables, no regressions)
      if (any(!op %in% c("~~", "~1")) ||
          any(!lhs %in% ovn) ||
          any(op == "~~" & !rhs %in% ovn)) {
        next
      }

      # all p*(p+1)/2 variances/covariances must be present and free
      cov_flag <- op == "~~"
      if (any(free[cov_flag] == 0L)) {
        next
      }
      i1 <- match(lhs[cov_flag], ovn)
      i2 <- match(rhs[cov_flag], ovn)
      pair_id <- pmin(i1, i2) + p * pmax(i1, i2)
      if (length(unique(pair_id)) != p * (p + 1L) / 2L) {
        next
      }

      # intercepts/means: all present; either free, or fixed at their h1
      # value (eg the zero within-level means of shared variables)
      int_flag <- op == "~1"
      if (lavmodel@meanstructure) {
        if (sum(int_flag) != p || any(duplicated(lhs[int_flag]))) {
          next
        }
        int_fixed <- which(int_flag & free == 0L)
        if (length(int_fixed) > 0L) {
          fixed_val <- lavpartable$start[row_idx[int_fixed]]
          h1_val <- mean_h1[match(lhs[int_fixed], ovn)]
          if (any(abs(fixed_val - h1_val) > 1e-10)) {
            next
          }
        }
      } else if (any(int_flag)) {
        next
      }

      # none of the free parameters may be involved in equality or
      # inequality constraints
      block_id <- free[free > 0L]
      if (any(block_id %in% dup_id)) {
        next
      }
      lab <- c(lavpartable$label[row_idx], lavpartable$plabel[row_idx])
      lab <- lab[nzchar(lab)]
      if (length(con_ref) > 0L && any(lab %in% con_ref)) {
        next
      }

      # this block is saturated: collect positions + h1 values
      b_free <- which(free > 0L)
      b_value <- numeric(length(b_free))
      for (i in seq_along(b_free)) {
        r <- b_free[i]
        if (op[r] == "~~") {
          b_value[i] <- cov_h1[match(lhs[r], ovn), match(rhs[r], ovn)]
        } else {
          b_value[i] <- mean_h1[match(lhs[r], ovn)]
        }
      }
      x_idx <- c(x_idx, match(free[b_free], uid))
      value <- c(value, b_value)
    } # levels
  } # groups

  list(x.idx = x_idx, value = value)
}
