# shared infrastructure for the univariate/bivariate ('uv'/'bv') modules:
#
#   lav_uvord.R  univariate ordinal probit (thresholds + slopes)
#   lav_uvreg.R  univariate linear regression (incl. residual variance)
#   lav_bvord.R  bivariate ordinal (polychoric/tetrachoric)
#   lav_bvmix.R  bivariate linear/ordinal (polyserial/biserial)
#   lav_bvreg.R  bivariate linear (pearson)
#
# YR/LDW 2026: extracted from the five modules (categorical-code refactor)
#
# All five modules share the same design:
#
#  - a 'cache' environment created by lav_XX_init_cache() holding the data,
#    fixed quantities, and the current parameter vector cache$theta;
#  - compute functions lav_XX_{loglik,logl,lik,sc,grad,hessian}_cache()
#    that evaluate in the cache environment via with(cache, ...);
#    NOTE: assignments made inside with(cache, ...) PERSIST in the cache.
#    The compute functions form a chain: the gradient assumes the
#    (log)lik function has been evaluated at cache$theta first, and the
#    hessian additionally assumes the gradient has been.  The nlminb
#    wrappers below maintain this invariant; call them in order if you
#    use the _cache functions directly.
#  - nlminb wrapper functions ('min fns') created once per module by
#    lav_uvbv_min_fns() below.

# check/normalize case weights (wt)
lav_uvbv_wt_check <- function(y = NULL, wt = NULL) {
  if (is.null(wt)) {
    wt <- rep(1, length(y))
  } else {
    if (length(y) != length(wt)) {
      lav_msg_stop(gettext("length y is not the same as length wt"))
    }
    if (any(wt < 0)) {
      lav_msg_stop(gettext("all weights should be positive"))
    }
  }
  wt
}

# select which derivative functions nlminb receives, given optim_method:
#   nlminb/nlminb2 = objective + gradient + hessian
#   nlminb1        = objective + gradient
#   nlminb0        = objective only
lav_uvbv_optim_fns <- function(optim_method = "nlminb",
                               min_fns = NULL) {
  out <- min_fns
  if (optim_method == "nlminb" || optim_method == "nlminb2") {
    # all three
  } else if (optim_method == "nlminb0") {
    out$gradient <- out$hessian <- NULL
  } else if (optim_method == "nlminb1") {
    out$hessian <- NULL
  }
  out
}

# create the standard nlminb objective/gradient/hessian wrappers for a
# module, implementing the state protocol described above:
#  - the objective stores 'x' in cache$theta and evaluates the (log)lik
#    (populating the intermediate quantities in the cache);
#  - gradient/hessian re-evaluate the earlier links of the chain first
#    whenever they are called at an 'x' that differs from cache$theta.
# All three return -1 * (quantity) / cache$n (minimization, per case).
#
# pre_objective (optional): function(x, cache) returning a replacement
# objective value (e.g. +Inf for an inadmissible 'x') or NULL to proceed.
lav_uvbv_min_fns <- function(logl_fun = NULL,
                             grad_fun = NULL,
                             hessian_fun = NULL,
                             pre_objective = NULL) {
  force(logl_fun); force(grad_fun); force(hessian_fun); force(pre_objective)

  objective <- function(x, cache = NULL) {
    if (!is.null(pre_objective)) {
      bad <- pre_objective(x, cache)
      if (!is.null(bad)) {
        return(bad)
      }
    }
    cache$theta <- x
    -1 * logl_fun(cache = cache) / cache$n
  }

  gradient <- function(x, cache = NULL) {
    # check if x has changed
    if (!all(x == cache$theta)) {
      cache$theta <- x
      tmp <- logl_fun(cache = cache)
      rm(tmp)
    }
    -1 * grad_fun(cache = cache) / cache$n
  }

  hessian <- function(x, cache = NULL) {
    # check if x has changed
    if (!all(x == cache$theta)) {
      cache$theta <- x
      tmp <- logl_fun(cache = cache)
      tmp <- grad_fun(cache = cache)
      rm(tmp)
    }
    -1 * hessian_fun(cache = cache) / cache$n
  }

  # a module may have no analytic gradient/hessian
  if (is.null(grad_fun)) {
    gradient <- NULL
  }
  if (is.null(hessian_fun)) {
    hessian <- NULL
  }

  list(objective = objective, gradient = gradient, hessian = hessian)
}
