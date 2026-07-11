# lavPredict() contains a collection of `predict' methods
# the unifying theme is that they all rely on the (unknown, to be estimated)
# or (known, apriori specified) values for the latent variables
#
# lv: latent variables (aka `factor scores')
# ov: predict linear part of y_i
#
# - YR 11 June 2013: first version, in order to get factor scores for the
#                    categorical case
# - YR 12 Jan 2014: refactoring + lav_predict_fy (to be used by estimator MML)
#
# - YR 25 Mar 2025: transformed factor scores (correlation-preserving)

# overload standard R function `predict'
setMethod(
  "predict", "lavaan",
  function(object, newdata = NULL, ...) {
    dotdotdot <- list(...)
    if (length(dotdotdot) > 0L) {
      for (j in seq_along(dotdotdot)) {
        lav_msg_warn(gettextf(
          "Unknown argument %s for %s", sQuote(names(dotdotdot)[j]),
          sQuote("predict"))
        )
      }
    }
    lavPredict(
      object = object, newdata = newdata, type = "lv", method = "EBM",
      fsm = FALSE, rel = FALSE, optim_method = "bfgs"
    )
  }
)

# efaList version
lav_efalist_predict <- function(object, ...) {
  # kill object$loadings if present
  object[["loadings"]] <- NULL

  if (length(object) == 1L) {
    # unlist
    object <- object[[1]]
  } else {
    # use the 'last' one per default
    object <- object[[length(object)]]
  }

  predict(object, ...)
}

# public function
lavPredict <- function(object, newdata = NULL, # keep order of predict(), 0.6-7 # nolint
                       type = "lv", method = "EBM", transform = FALSE,
                       se = "none", acov = "none", label = TRUE, fsm = FALSE,
                       mdist = FALSE, rel = FALSE,
                       append_data = FALSE, assemble = FALSE, # or TRUE?
                       level = 1L, optim_method = "bfgs", eta = NULL,
                       parallel = c("auto", "no", "multicore", "snow"),
                       ncpus = NULL, cl = NULL,
                       drop_list_single_group = TRUE,
                       mdist_draws = 2000L,
                      ...) {
  dotdotdot <- list(...)
  lav_adapt_func(environment(), dotdotdot, NULL)
  # check object
  object <- lav_object_check_version(object)

  parallel <- match.arg(parallel)

  # catch efaList objects
  if (inherits(object, "efaList")) {
    # kill object$loadings if present
    object[["loadings"]] <- NULL
    if (length(object) == 1L) {
      # unlist
      object <- object[[1]]
    } else {
      # use the 'last' one per default
      object <- object[[length(object)]]
    }
  }

  stopifnot(inherits(object, "lavaan"))
  lavmodel <- object@Model
  lavdata <- object@Data
  lavsamplestats <- object@SampleStats
  lavh1 <- object@h1
  lavimplied <- object@implied
  lavpartable <- object@ParTable

  # random slopes? use the dedicated empirical Bayes route:
  # - level = 2: posterior means of the between latent variables
  #   (including the random slopes); se = "standard" gives the
  #   posterior standard deviations
  # - level = 1: within factor scores, conditional on the posterior
  #   means of the cluster-level random effects (exact posterior means)
  if (length(lavmodel@rv.ov) > 0L || length(lavmodel@rv.lv) > 0L) {
    if (tolower(type) != "lv") {
      lav_msg_stop(gettext(
        "for models with random slopes, only type = \"lv\" is
         supported (for now)."))
    }
    if (!tolower(method) %in% c("ebm", "regression")) {
      lav_msg_stop(gettext(
        "for models with random slopes, only method = \"EBM\" is
         supported (the posterior means)."))
    }
    if (!is.null(newdata)) {
      lav_msg_stop(gettext(
        "newdata is not supported (yet) for models with random
         slopes."))
    }
    if (fsm || mdist || rel || transform || append_data) {
      lav_msg_stop(gettext(
        "the fsm/mdist/rel/transform/append_data options are not
         supported (yet) for models with random slopes."))
    }
    se_flag <- !is.null(se) && se != "none"
    if (se_flag && level == 1L) {
      lav_msg_stop(gettext(
        "standard errors for the level-1 factor scores are not
         available (yet) for models with random slopes."))
    }
    eb <- lav_mvn_cl_rs_eb(
      lavmodel = lavmodel, lavdata = lavdata,
      lavcache = object@Cache, se = se_flag
    )
    if (level == 2L) {
      out1 <- eb$l2
      se1 <- eb$se2
    } else {
      out1 <- eb$l1
      se1 <- NULL
    }
    if (!label) {
      colnames(out1) <- NULL
    }
    class(out1) <- c("lavaan.matrix", "matrix")
    if (!drop_list_single_group) {
      out1 <- list(out1)
      if (!is.null(se1)) {
        se1 <- list(se1)
      }
    }
    if (se_flag) {
      attr(out1, "se") <- se1
    }
    return(out1)
  }

  # warn if the model does not contain any 'regular' latent variables
  if (length(lav_object_vnames(lavpartable, "lv.regular")) == 0L) {
    lav_msg_warn(gettextf("fitted model does not contain regular
    (i.e., measured) latent variables; the matrix of factor scores may
    contain no columns"))
  }

  res <- lav_predict_internal(
    lavmodel = lavmodel, lavdata = lavdata,
    lavsamplestats = lavsamplestats, lavimplied = lavimplied, lavh1 = lavh1,
    lavpartable = lavpartable, newdata = newdata, type = type, method = method,
    transform = transform, se = se, acov = acov, label = label,
    fsm = fsm, rel = rel,
    mdist = mdist, append_data = append_data, assemble = assemble,
    level = level, optim_method = optim_method, eta = eta,
    parallel = parallel, ncpus = ncpus, cl = cl,
    drop_list_single_group = drop_list_single_group,
    mdist_draws = mdist_draws
  )

  res
}

# internal version, to be used if lavobject does not exist yet
lav_predict_internal <- function(lavmodel = NULL,
                               lavdata = NULL,
                               lavsamplestats = NULL,
                               lavh1 = NULL,
                               lavimplied = NULL,
                               lavpartable = NULL,
                               # standard options
                               newdata = NULL, # keep order of predict(), 0.6-7
                               type = "lv", method = "EBM", transform = FALSE,
                               se = "none", acov = "none", label = TRUE,
                               fsm = FALSE, rel = FALSE,
                               mdist = FALSE,
                               append_data = FALSE, assemble = FALSE, # or TRUE?
                               level = 1L, optim_method = "bfgs", eta = NULL,
                               parallel = "no", ncpus = NULL, cl = NULL,
                               drop_list_single_group = TRUE,
                               mdist_draws = 2000L) {
  # type
  type <- tolower(type)
  lavpta <- lav_pt_attributes(lavpartable)
  if (type %in% c("latent", "lv", "factor", "factor.score", "factorscore")) {
    type <- "lv"
  } else if (type %in% c("ov", "yhat")) {
    type <- "yhat"
  } else if (type %in% c("residuals", "resid", "error")) {
    type <- "resid"
  }

  # if resid, not for categorical, unless mdist = TRUE: in that case the
  # expected residuals of the latent responses are computed by Monte Carlo
  # integration (Mansolf & Reise, 2017; see lav_predict_mdist_cat)
  if (type == "resid" && lavmodel@categorical && !mdist) {
    lav_msg_stop(gettext(
      "casewise residuals not available if data is categorical
      (unless mdist = TRUE)"))
  }

  # append_data? check level
  if (append_data && level > 1L) {
    lav_msg_warn(gettext("append_data not available if level > 1L"))
    append_data <- FALSE
  }

  # mdist? -> fsm = TRUE (continuous data only; for categorical data the
  # distances are computed by Monte Carlo integration, and no factor score
  # matrix is involved)
  mdist_cat_flag <- mdist && lavmodel@categorical
  if (mdist && !mdist_cat_flag) {
    fsm <- TRUE
  }

  # se?
  if (acov != "none") {
    se <- acov # acov implies se
  }
  if (se != "none") {
    if (is.logical(se) && se) {
      se <- "standard"
      if (acov != "none") {
        acov <- se # reverse-imply upstream
      }
    }
    if (type != "lv") {
      lav_msg_stop(gettext("standard errors only available if type = \"lv\""))
    }
    # NOTE: for categorical data, the standard errors are computed in
    # lav_predict_eta_ebm_ml() from the curvature of the (EBM/ML) objective at
    # the optimum, and so differ from one response pattern to the next.
    # if(lavdata@missing %in% c("ml", "ml.x")) {
    #    se <- acov <- "none"
    #    warning("lavaan WARNING: standard errors not
    #                 available (yet) for missing data + fiml")
    # }
  }

  # need full data set supplied
  if (is.null(newdata)) {
    # use internal copy:
    if (lavdata@data.type != "full") {
      lav_msg_stop(gettext(
        "sample statistics were used for fitting and newdata is empty"))
    } else if (is.null(lavdata@X[[1]])) {
      lav_msg_stop(gettext("no local copy of data; FIXME!"))
    } else {
      data_obs <- lavdata@X
      ov_names <- lavdata@ov.names
    }
    exo <- lavdata@eXo
  } else {
    ov <- lavdata@ov
    # for multilevel models, the cluster variable must be present in
    # newdata, so that a fresh cluster structure (Lp) can be built
    if (lavdata@nlevels > 1L &&
        !all(lavdata@cluster %in% names(newdata))) {
      lav_msg_stop(gettextf(
        "newdata does not contain the cluster variable(s): %s",
        lav_msg_view(lavdata@cluster, "none")))
    }
    new_data <- lav_lavdata(
      data = newdata,
      group = lavdata@group,
      cluster = lavdata@cluster,
      ov_names = lavdata@ov.names,
      ov_names_x = lavdata@ov.names.x,
      ov_names_l = lavdata@ov.names.l,
      ordered = ov$name[ov$type == "ordered"],
      lavoptions = list(
        std.ov = lavdata@std.ov,
        group.label = lavdata@group.label,
        level.label = lavdata@level.label,
        missing = lavdata@missing
      ), # was FALSE before?
      allow_single_case = TRUE
    )
    # if ordered, check if number of levels is still the same (new in 0.6-7)
    if (lavmodel@categorical) {
      orig_ordered_idx <- which(lavdata@ov$type == "ordered")
      orig_ordered_lev <- lavdata@ov$nlev[orig_ordered_idx]
      match_new_idx <- match(
        lavdata@ov$name[orig_ordered_idx],
        new_data@ov$name
      )
      new_ordered_lev <- new_data@ov$nlev[match_new_idx]
      if (any(orig_ordered_lev - new_ordered_lev != 0)) {
        lav_msg_stop(
          gettext("mismatch number of categories for some ordered variables
                  in newdata compared to original data.")
        )
      }
    }
    data_obs <- new_data@X
    exo <- new_data@eXo
    ov_names <- new_data@ov.names

    # for multilevel models, replace lavdata by the freshly-built
    # newdata object: the downstream code needs the newdata cluster
    # structure (Lp), missing patterns (Mp) and dimensions (new in
    # 0.7-1); lavsamplestats is kept as-is (only used for the
    # model-implied quantities)
    if (lavdata@nlevels > 1L) {
      lavdata <- new_data
    }
  }

  # categorical data + mdist: expected (squared) M-distances of the latent
  # responses, via Monte Carlo integration (Mansolf & Reise, 2017)
  if (mdist_cat_flag) {
    if (!type %in% c("lv", "resid", "yhat")) {
      lav_msg_stop(gettext(
        "mdist is only available if type is one of: lv resid yhat"))
    }
    if (level > 1L) {
      lav_msg_stop(gettext(
        "mdist with categorical data is not available if level > 1L"))
    }
    mdist_cat <- lav_predict_mdist_cat(
      lavmodel = lavmodel, lavdata = lavdata, lavimplied = lavimplied,
      data_obs = data_obs, exo = exo,
      type = type, method = method, ndraws = mdist_draws
    )
    # (squared metric, as in the continuous case)
    mdist_1 <- lapply(mdist_cat, "[[", "d2")
  }

  if (type == "lv") {
    if (!is.null(eta)) {
      lav_msg_warn(gettext("lvs will be predicted here;
                           supplying eta= has no effect"))
    }

    # post fit check (lv pd?)
    # ok <- lav_object_post_check(object)
    # if(!ok) {
    #    stop("lavaan ERROR: lavInspect(,\"post.check\") is not TRUE;
    #       factor scores can not be computed. See the WARNING message.")
    # }

    out <- lav_predict_eta(
      lavobject = NULL, lavmodel = lavmodel,
      lavdata = lavdata, lavsamplestats = lavsamplestats,
      lavimplied = lavimplied, lavpta = lavpta,
      se = se, acov = acov, level = level,
      data_obs = data_obs, exo = exo, method = method,
      fsm = fsm, rel = rel, transform = transform, optim_method = optim_method,
      parallel = parallel, ncpus = ncpus, cl = cl
    )

    # extract fsm here
    if (fsm) {
      fsm_1 <- attr(out, "fsm")
    }

  # extract rel here
  if (rel) {
    rel_1 <- attr(out, "rel")
  }

    # extract se here
    if (se != "none") {
      se_1 <- attr(out, "se")
      if (acov != "none") {
        acov_1 <- attr(out, "acov")
      }
    }

    # remove dummy lv? (removes attr!)
    out <- lapply(seq_len(lavdata@ngroups), function(g) {
      # determine block
      b <- lav_predict_block_idx(lavdata, g, level)
      lv_idx <- c(
        lavmodel@ov.y.dummy.lv.idx[[b]],
        lavmodel@ov.x.dummy.lv.idx[[b]]
      )
      ret <- out[[g]]
      if (length(lv_idx) > 0L) {
        ret <- out[[g]][, -lv_idx, drop = FALSE]
      }
      ret
    })

    # we need to remove the dummy's before we transform
    # (update 0.6-20: no longer needed... as transform happens internally)
    if (fsm) {
      fsm_1 <- lapply(seq_len(lavdata@ngroups), function(g) {
        # determine block
        b <- lav_predict_block_idx(lavdata, g, level)
        lv_idx <- c(
          lavmodel@ov.y.dummy.lv.idx[[b]],
          lavmodel@ov.x.dummy.lv.idx[[b]]
        )
        # ov.idx <- lavmodel@ov.x.dummy.ov.idx[[b]]
        # or should we use pta$vidx$ov.ind?
        ov_ind <- lavpta$vidx$ov.ind[[b]]
        ret <- fsm_1[[g]]
        if (length(lv_idx) > 0L) {
          if (is.matrix(fsm_1[[g]])) {
            ret <- fsm_1[[g]][-lv_idx, ov_ind, drop = FALSE]
          } else if (is.list(fsm_1[[g]])) {
            fsm_1[[g]] <- lapply(fsm_1[[g]], function(x) {
              ret <- x[-lv_idx, ov_ind, drop = FALSE]
              ret
            })
          }
        }
        ret
      })
    }

#   # new in 0.6-16
#   # we assume the dummy lv's have already been removed
#     if (transform) {
#       # VETA <- lav_model_veta(lavmodel = lavmodel, remove.dummy.lv = TRUE)
#       EETA <- lav_model_eeta(
#         lavmodel = lavmodel,
#         lavsamplestats = lavsamplestats, remove.dummy.lv = TRUE
#       )
#       # compute transformation matrix
#       if (tolower(method) %in% c("ebm", "regression")) {
#         tmat <- lav_predict_tmat_green(lavmodel = lavmodel,
#                                        lavimplied = lavimplied)
#       } else {
#         tmat <- lav_predict_tmat_det(lavmodel = lavmodel,
#                                      lavimplied = lavimplied)
#       }
#
#       # update fsm_1
#       if (fsm) {
#         fsm_1 <- lapply(seq_len(lavdata@ngroups), function(g) {
#           # determine block
#           if (lavdata@nlevels == 1L) {
#             b <- g
#           } else {
#             b <- (g - 1) * lavdata@nlevels + level
#           }
#           ret <- tmat[[b]] %*% fsm_1[[g]]
#           ret
#         })
#       }
#
#       out <- lapply(seq_len(lavdata@ngroups), function(g) {
#         # determine block
#         if (lavdata@nlevels == 1L) {
#           b <- g
#         } else {
#           b <- (g - 1) * lavdata@nlevels + level
#         }
#
#         FS.centered <- scale(out[[g]], center = TRUE, scale = FALSE)
#         #FS.cov <- crossprod(FS.centered) / nrow(FS.centered)
#         #FS.cov.inv <- try(solve(FS.cov), silent = TRUE)
#         #if (inherits(FS.cov.inv, "try-error")) {
#         #  lav_msg_warn(
#         #    gettext("could not invert (co)variance matrix of factor scores;
#         #            returning original factor scores."))
#         #  return(out[[g]])
#         #}
#         #fs.inv.sqrt <- lav_mat_sym_sqrt(FS.cov.inv)
#         #veta.sqrt <- lav_mat_sym_sqrt(VETA[[g]])
#         #tmp <- FS.centered %*% fs.inv.sqrt %*% veta.sqrt
#         tmp <- FS.centered %*% t(tmat[[b]])
#         ret <- t(t(tmp) + drop(EETA[[g]]))
#
#         ret
#       })
#     }

    # new in 0.6-17
    if (mdist && !mdist_cat_flag) {
      veta <- lav_model_veta(lavmodel = lavmodel, remove_dummy_lv = TRUE)
      eeta <- lav_model_eeta(
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats, remove_dummy_lv = TRUE
      )
      mdist_1 <- lapply(seq_len(lavdata@ngroups), function(g) {
        a <- fsm_1[[g]]
        sigma <- lavimplied$cov[[g]]
        if (transform) {
          fs_cov <- veta[[g]]
        } else {
          fs_cov <- a %*% sigma %*% t(a)
        }
        fs_cov_inv <- solve(fs_cov)
        # Mahalanobis distance
        fs_c <- t(t(out[[g]]) - eeta[[g]]) # center
        df_squared <- rowSums((fs_c %*% fs_cov_inv) * fs_c)
        ret <- df_squared # squared!
        ret
      })
    }

    # append original/new data? (also remove attr)
    if (append_data && level == 1L) {
      out <- lapply(seq_len(lavdata@ngroups), function(g) {
        ret <- cbind(out[[g]], data_obs[[g]])
        ret
      })
    }

    if (se != "none") {
      se_1 <- lapply(seq_len(lavdata@ngroups), function(g) {
        # determine block
        b <- lav_predict_block_idx(lavdata, g, level)
        lv_idx <- c(
          lavmodel@ov.y.dummy.lv.idx[[b]],
          lavmodel@ov.x.dummy.lv.idx[[b]]
        )
        ret <- se_1[[g]]
        if (length(lv_idx) > 0L) {
          ret <- se_1[[g]][, -lv_idx, drop = FALSE]
        }
        ret
      })
      if (acov != "none") {
        acov_1 <- lapply(seq_len(lavdata@ngroups), function(g) {
          # determine block
          b <- lav_predict_block_idx(lavdata, g, level)
          lv_idx <- c(
            lavmodel@ov.y.dummy.lv.idx[[b]],
            lavmodel@ov.x.dummy.lv.idx[[b]]
          )
          ret <- acov_1[[g]]
          if (length(lv_idx) > 0L) {
            if (is.matrix(acov_1[[g]])) {
              ret <- acov_1[[g]][-lv_idx, -lv_idx, drop = FALSE]
            } else if (is.list(acov_1[[g]])) {
              ret <- lapply(acov_1[[g]], function(x) {
                ret <- x[-lv_idx, -lv_idx, drop = FALSE]
                ret
              })
            }
          }
          ret
        })
      } # acov
    } # se

    # label?
    if (label) {
      for (g in seq_len(lavdata@ngroups)) {
        gg <- lav_predict_block_idx(lavdata, g, level)

        if (append_data) {
          colnames(out[[g]]) <- c(
            lavpta$vnames$lv[[gg]],
            ov_names[[g]]
          ) # !not gg
        } else {
          colnames(out[[g]]) <- lavpta$vnames$lv[[gg]]
        }

        if (fsm) {
          # column names = manifest indicators. For models with composites the
          # factor-score matrix also has columns for the composite indicators,
          # so use the full set of observed indicators (ov = ov.ind + ov.cind)
          # rather than ov.ind alone.
          fsm_ov <- if (lavmodel@composites) {
            lavpta$vnames$ov[[gg]]
          } else {
            lavpta$vnames$ov.ind[[gg]]
          }
          if (is.null(fsm_1[[g]])) {
            # skip
          } else if (is.matrix(fsm_1[[g]])) {
            if (ncol(fsm_1[[g]]) == length(fsm_ov)) {
              dimnames(fsm_1[[g]]) <- list(lavpta$vnames$lv[[gg]], fsm_ov)
            }
          } else if (is.list(fsm_1[[g]])) {
            fsm_1[[g]] <- lapply(fsm_1[[g]], function(x) {
              if (ncol(x) == length(fsm_ov)) {
                dimnames(x) <- list(lavpta$vnames$lv[[gg]], fsm_ov)
              }
              x
            })
          }
        }

        if (se != "none") {
          if (!is.null(se_1[[g]])) {
            colnames(se_1[[g]]) <- lavpta$vnames$lv[[gg]]
          }
        }

        if (rel) {
          if (!is.null(rel_1[[g]])) {
            names(rel_1[[g]]) <- lavpta$vnames$lv[[gg]]
          }
        }

        if (acov != "none") {
          if (is.null(acov_1[[g]])) {
            # skip
          } else if (is.matrix(acov_1[[g]])) {
            dimnames(acov_1[[g]]) <- list(
              lavpta$vnames$lv[[gg]],
              lavpta$vnames$lv[[gg]]
            )
          } else if (is.list(acov_1[[g]])) {
            acov_1[[g]] <- lapply(acov_1[[g]], function(x) {
              dimnames(x) <- list(
                lavpta$vnames$lv[[gg]],
                lavpta$vnames$lv[[gg]]
              )
              x
            })
          }
        }
      } # g

      # group.labels
      if (lavdata@ngroups > 1L) {
        names(out) <- lavdata@group.label
        if (se != "none") {
          names(se_1) <- lavdata@group.label
        }
        if (acov != "none") {
          names(acov_1) <- lavdata@group.label
        }
      }
    } # label

    # yhat: estimated value for the observed indicators, given (estimated)
    # factor scores
    # resid: y - yhat
  } else if (type %in% c("yhat", "resid")) {
    resid_flag <- type == "resid"
    if (resid_flag && mdist_cat_flag) {
      # categorical data: the expected casewise residuals of the latent
      # responses (a by-product of the Monte Carlo integration)
      out <- lapply(mdist_cat, "[[", "resid")
    } else {
      out <- lav_predict_yhat(
        lavobject = NULL, lavmodel = lavmodel,
        lavdata = lavdata, lavsamplestats = lavsamplestats,
        lavimplied = lavimplied,
        data_obs = data_obs, exo = exo,
        eta = eta, method = method, optim_method = optim_method,
        fsm = fsm,
        resid_flag = resid_flag
      )
      if (fsm) {
        fsm_1 <- attr(out, "fsm")
      }
    }

    # label?
    if (label) {
      for (g in seq_len(lavdata@ngroups)) {
        if (ncol(out[[g]]) == length(lavpta$vnames$ov[[g]])) {
          colnames(out[[g]]) <- lavpta$vnames$ov[[g]]
        } else {
          # conditional.x: the columns only cover the non-exogenous ov's
          colnames(out[[g]]) <- lavpta$vnames$ov.nox[[g]]
        }
      }
    }

    # mdist
    if (mdist && !mdist_cat_flag) {
      mm_lambda <- lav_model_lambda(
        lavmodel = lavmodel,
        remove_dummy_lv = FALSE,
        use_wmat = TRUE
      )
      mdist_1 <- lapply(seq_len(lavdata@ngroups), function(g) {
        sigma <- lavimplied$cov[[g]]
        la <- mm_lambda[[g]]
        if (type == "resid") {
          ila <- diag(ncol(sigma)) - la %*% fsm_1[[g]]
          omega_e <- ila %*% sigma %*% t(ila)
          eig <- eigen(omega_e, symmetric = TRUE)
          a <- eig$vectors[, seq_len(nrow(la) - ncol(la)),
            drop = FALSE
          ]
        } else if (type == "yhat") {
          laa <- la %*% fsm_1[[g]]
          omega_e <- laa %*% sigma %*% t(laa)
          eig <- eigen(omega_e, symmetric = TRUE)
          a <- eig$vectors[, seq_len(ncol(la)), drop = FALSE]
        }
        out_a <- apply(out[[g]], 1L, function(x) {
          colSums(a * x, na.rm = TRUE)
        })
        if (is.matrix(out_a)) {
          out_a <- t(out_a)
        } else {
          out_a <- as.matrix(out_a)
        }
        # if(lavmodel@meanstructure) {
        #     est.mean <- drop(t(lavimplied$mean[[g]]) %*% a)
        #     if(type == "resid") {
        #         obs.mean <- drop(lavh1$implied$mean[[g]] %*% a)
        #         est.mean <- drop(t(lavimplied$mean[[g]]) %*% a)
        #         outA.mean <- obs.mean - est.mean
        #     } else if(type == "yhat") {
        #         outA.mean <- est.mean
        #     }
        # } else {
        #     outA.mean <- colMeans(out_a)
        # }
        out_a_cov <- t(a) %*% omega_e %*% a
        out_a_cov_inv <- solve(out_a_cov)
        # Mahalanobis distance
        # out_a_c <- t( t(out_a) - outA.mean ) # center
        out_a_c <- out_a
        df_squared <- rowSums((out_a_c %*% out_a_cov_inv) * out_a_c)
        ret <- df_squared # squared!
        ret
      })
    }


    # density for each observed item, given (estimated) factor scores
  } else if (type == "fy") {
    out <- lav_predict_fy(
      lavobject = NULL, lavmodel = lavmodel,
      lavdata = lavdata, lavsamplestats = lavsamplestats,
      lavimplied = lavimplied,
      data_obs = data_obs, exo = exo,
      eta = eta, method = method, optim_method = optim_method
    )

    # label?
    if (label) {
      for (g in seq_len(lavdata@ngroups)) {
        colnames(out[[g]]) <- lavpta$vnames$ov[[g]]
      }
    }
  } else {
    lav_msg_stop(gettext("type must be one of: lv yhat fy"))
  }

  # lavaan.matrix
  out <- lapply(out, "class<-", c("lavaan.matrix", "matrix"))

  if (lavdata@ngroups == 1L && drop_list_single_group) {
    res <- out[[1L]]
  } else {
    res <- out
  }

  # assemble multiple groups into a single data.frame? (new in 0.6-4)
  if (lavdata@ngroups > 1L && assemble) {
    if (!is.null(newdata)) {
      lavdata <- new_data
    }
    data_1 <- matrix(as.numeric(NA),
      nrow = sum(unlist(lavdata@norig)),
      ncol = ncol(out[[1L]])
    ) # assume == per g
    colnames(data_1) <- colnames(out[[1L]])
    for (g in seq_len(lavdata@ngroups)) {
      data_1[lavdata@case.idx[[g]], ] <- out[[g]]
    }
    data_1 <- as.data.frame(data_1, stringsAsFactors = FALSE)

    if (!is.null(newdata)) {
      data_1[, lavdata@group] <- newdata[, lavdata@group]
    } else {
      # add group
      data_1[, lavdata@group] <- rep(as.character(NA), nrow(data_1))
      if (lavdata@missing == "listwise") {
        # we will lose the group label of omitted variables!
        data_1[unlist(lavdata@case.idx), lavdata@group] <-
          rep(lavdata@group.label, unlist(lavdata@nobs))
      } else {
        data_1[unlist(lavdata@case.idx), lavdata@group] <-
          rep(lavdata@group.label, unlist(lavdata@norig))
      }
    }

    res <- data_1
  }

  if (fsm && type == "lv") {
    attr(res, "fsm") <- fsm_1
  }

  if (rel && type == "lv") {
    attr(res, "rel") <- rel_1
  }

  if (mdist) {
    attr(res, "mdist") <- mdist_1
  }

  if (se != "none") {
    attr(res, "se") <- se_1
    # return full sampling covariance matrix?
    if (acov == "standard") {
      attr(res, "acov") <- acov_1
    }
  }

  res
}

# map a group index 'g' (and a level) to the corresponding block index 'b'
# in the per-block model matrices (GLIST, veta, eeta, ...). For single-level
# models the block is just the group; for multilevel models the blocks are
# interleaved per group as (g - 1) * nlevels + level.
lav_predict_block_idx <- function(lavdata, g, level = 1L) {
  if (lavdata@nlevels == 1L) {
    g
  } else {
    (g - 1L) * lavdata@nlevels + level
  }
}

# internal function
lav_predict_eta <- function(lavobject = NULL, # for convenience
                            # sub objects
                            lavmodel = NULL, lavdata = NULL,
                            lavsamplestats = NULL,
                            lavimplied = NULL, lavpta = NULL,
                            # new data
                            data_obs = NULL, exo = NULL,
                            # options
                            method = "EBM",
                            fsm = FALSE,
              rel = FALSE,
                            transform = FALSE,
                            se = "none", acov = "none",
                            level = 1L,
                            optim_method = "bfgs",
                            parallel = "no", ncpus = NULL, cl = NULL) {
  # full object?
  if (inherits(lavobject, "lavaan")) {
    lavdata <- lavobject@Data
  } else {
    stopifnot(!is.null(lavdata))
  }

  # method
  method <- tolower(method)

  # alias
  if (method == "regression") {
    method <- "ebm"
  } else if (method == "bartlett" || method == "bartlet") {
    method <- "ml"
  }

  # at this point, method is either "ebm" (= regression/EB) or "ml" (= Bartlett)
  if (!method %in% c("ebm", "ml")) {
    lav_msg_stop(gettextf("unknown method: %s.", method))
  }

  if (all(lavdata@ov$type == "numeric")) {
    # normal case: closed-form regression (EBM) or Bartlett (ML) scores
    out <- lav_predict_eta_normal(
      lavobject = lavobject,
      lavmodel = lavmodel, lavdata = lavdata,
      lavimplied = lavimplied, lavpta = lavpta, se = se, acov = acov,
      level = level, lavsamplestats = lavsamplestats,
      data_obs = data_obs, exo = exo, fsm = fsm, rel = rel,
      transform = transform, method = method
    )
  } else {
    # non-normal case: numerical optimisation of the latent scores
    out <- lav_predict_eta_ebm_ml(
      lavobject = lavobject,
      lavmodel = lavmodel, lavdata = lavdata,
      lavsamplestats = lavsamplestats, se = se, acov = acov,
      level = level, data_obs = data_obs, exo = exo,
      ml = (method == "ml"), optim_method = optim_method,
      parallel = parallel, ncpus = ncpus, cl = cl
    )
  }

  out
}


# (pseudo-)inverse used in the factor-score machinery.
# For a converged model the model-implied Sigma is positive definite, and
# (lambda' Sigma.inv lambda) is invertible in all but a few "by construction"
# singular settings (e.g. a single-indicator construct with a fixed zero
# residual variance). In particular, for models with composites both matrices
# are well-conditioned, so we use an exact (and fast) solve() here; we only fall
# back to the Moore-Penrose pseudo-inverse for the genuinely singular cases that
# originally motivated MASS::ginv(). This keeps the common path exact and avoids
# relying on ginv() for composites.
lav_predict_solve <- function(x) {
  out <- try(solve(x), silent = TRUE)
  if (inherits(out, "try-error")) {
    out <- MASS::ginv(x)
  }
  # drop dimnames so the result behaves exactly like MASS::ginv() downstream
  dimnames(out) <- NULL
  out
}

# factor scores - normal case
# NOTE: this is the classic 'regression' method; for the linear/continuous
#       case, this is equivalent to both EB and EBM
# factor scores - normal (continuous) case - closed-form
#
# Two methods are supported, differing only in the factor-score coefficient
# matrix 'C' (mapping centered data to latent scores) and in the corresponding
# conditional covariance used for standard errors:
#
#  - method = "regression" (aka EB/EBM/Thurstone): C = Psi Lambda' Sigma.inv
#       For the linear/continuous case this is equivalent to both EB and EBM.
#
#  - method = "bartlett" (aka ML): C = (Lambda' Sigma.inv Lambda).inv Lambda'
#                                      Sigma.inv
#       The usual formula uses Theta.inv, but to deal with a singular THETA
#       (zeroes on the diagonal) we use the 'GLS' version with Sigma.inv
#       instead (Bentler & Yuan, 1997, 'Optimal Conditionally Unbiased
#       Equivariant Factor Score Estimators', in Berkane (Ed) 'Latent variable
#       modeling with applications to causality', Springer-Verlag).
#
# NOTE: instead of solve(), we use lav_predict_solve() (a pseudo-inverse
#       fallback) for the special settings where -by construction-
#       (Lambda' Sigma.inv Lambda) is singular. For Bartlett scores this
#       destroys the conditionally-unbiased property.
lav_predict_eta_normal <- function(lavobject = NULL, # for convenience
                                   # sub objects
                                   lavmodel = NULL, lavdata = NULL,
                                   lavsamplestats = NULL,
                                   lavimplied = NULL, lavpta = NULL,
                                   # optional new data
                                   data_obs = NULL, exo = NULL,
                                   se = "none", acov = "none", level = 1L,
                                   fsm = FALSE, rel = FALSE,
                                   transform = FALSE,
                                   method = "regression") {
  # full object?
  if (inherits(lavobject, "lavaan")) {
    lavmodel <- lavobject@Model
    lavdata <- lavobject@Data
    lavsamplestats <- lavobject@SampleStats
    lavimplied <- lavobject@implied
    lavpta <- lavobject@pta
  } else {
    stopifnot(
      !is.null(lavmodel), !is.null(lavdata),
      !is.null(lavsamplestats), !is.null(lavimplied)
    )
  }

  # Bartlett (ML) or regression (EB/EBM) scores?
  bartlett <- tolower(method) %in% c("ml", "bartlett", "bartlet")

  # method-specific factor-score coefficient matrix 'C' (maps centered data to
  # latent scores), as a function of Lambda, Sigma.inv and Psi (= veta_g)
  #
  # higher-order factors (Bartlett only): a higher-order factor has no
  # observed indicators, so its column in Lambda is empty, and the
  # (pseudo-inverse based) Bartlett mapping returns all-zero scores for this
  # factor; if lambda_star (the 'collapsed' loadings, see below) is provided,
  # the rows for the higher-order factors are taken from the Bartlett mapping
  # based on lambda_star instead (the other rows are left untouched)
  compute_fsc <- function(lambda, sigma_inv, veta_g, lambda_star = NULL,
                          fixes = NULL) {
    if (bartlett) {
      fsc <- lav_predict_solve(t(lambda) %*% sigma_inv %*% lambda) %*%
        t(lambda) %*% sigma_inv
      if (!is.null(lambda_star)) {
        for (fx in fixes) {
          ls_t <- lambda_star[, fx$tcols, drop = FALSE]
          fsc_star <- lav_predict_solve(
            t(ls_t) %*% sigma_inv %*% ls_t) %*% t(ls_t) %*% sigma_inv
          fsc[fx$cols, ] <- fsc_star[match(fx$cols, fx$tcols), ,
                                     drop = FALSE]
        }
      }
      fsc
    } else {
      veta_g %*% t(lambda) %*% sigma_inv
    }
  }

  # method-specific conditional (co)variance of the factor scores, used for the
  # standard errors and the full sampling covariance matrix (acov)
  #
  # if 'fsc' (and 'sigma') are provided (higher-order factors, Bartlett only),
  # the general expression Var(fs - eta) = C Sigma C' - C Lambda Psi -
  # Psi Lambda' C' + Psi is used instead: the standard expression below
  # assumes C Lambda = I, which no longer holds for the higher-order rows
  compute_fscov <- function(lambda, sigma_inv, veta_g,
                            fsc = NULL, sigma = NULL) {
    if (bartlett) {
      # traditional formula uses solve(Lambda' Theta.inv Lambda); we use the
      # Sigma.inv version (minus Psi) to handle negative/zero variances, with
      # a pseudo-inverse fallback (lav_predict_solve)
      if (is.null(fsc)) {
        lav_predict_solve(t(lambda) %*% sigma_inv %*% lambda) - veta_g
      } else {
        cl_psi <- fsc %*% lambda %*% veta_g
        fsc %*% sigma %*% t(fsc) - cl_psi - t(cl_psi) + veta_g
      }
    } else {
      veta_g - veta_g %*% t(lambda) %*% sigma_inv %*% lambda %*% veta_g
    }
  }

  if (is.null(data_obs)) {
    data_obs <- lavdata@X
    newdata_flag <- FALSE
  } else {
    newdata_flag <- TRUE
  }
  # eXo not needed

  # missings? and missing = "ml"?
  if (lavdata@missing %in% c("ml", "ml.x")) {
    if (newdata_flag) {
      mp_1 <- vector("list", lavdata@ngroups)
      for (g in seq_len(lavdata@ngroups)) {
        mp_1[[g]] <- lav_data_mi_patterns(data_obs[[g]])
      }
    } else {
      mp_1 <- lavdata@Mp
    }
  }

  mm_lambda <- lav_model_lambda(lavmodel = lavmodel, remove_dummy_lv = FALSE,
                                use_wmat = TRUE)
  sigma_hat <- lavimplied$cov
  sigma_inv <- lapply(sigma_hat, lav_predict_solve)
  veta <- lav_model_veta(lavmodel = lavmodel)
  eeta <- lav_model_eeta(lavmodel = lavmodel, lavsamplestats = lavsamplestats)
  ey <- lav_model_ey(lavmodel = lavmodel, lavsamplestats = lavsamplestats)

  # Bartlett + higher-order factors (new in 0.7-1): collapse the measurement
  # chain, just like in sam(): LAMBDA.star = LAMBDA %*% solve(I - B), where B
  # only contains the factor loadings of the higher-order factors (the rows
  # of the latent indicators); for each higher-order factor, the columns of
  # LAMBDA.star for its 'layer' (the factor itself + all factors that are
  # neither an ancestor nor a descendant in the measurement chain, keeping
  # only the topmost ones) define a collapsed measurement model
  # y = LAMBDA.star ETA2 + epsilon.star, and the (Sigma.inv form of the)
  # Bartlett mapping based on these columns provides the mapping-matrix row
  # for that factor; factors of the same layer share the same mapping
  mm_lambda_star <- vector("list", length = length(mm_lambda))
  ho_fixes <- vector("list", length = length(mm_lambda))
  if (bartlett && !is.null(lavpta) && !is.null(lavpta$vnames$lv.ind)) {
    lambda_mm_idx <- which(names(lavmodel@GLIST) == "lambda")
    mm_idx_block <- lav_model_group_mm_indices(lavmodel@nmat)
    for (b in seq_along(mm_lambda)) {
      lv_ind_names <- lavpta$vnames$lv.ind[[b]]
      if (length(lv_ind_names) == 0L) next
      empty_idx <- which(apply(mm_lambda[[b]], 2L,
                               function(x) all(x == 0)))
      if (length(empty_idx) == 0L) next
      mlist <- lavmodel@GLIST[mm_idx_block[[b]]]
      if (is.null(mlist$beta)) next
      lv_names <- lavmodel@dimNames[[lambda_mm_idx[b]]][[2L]]
      # keep only the rows of the latent indicators; other (structural)
      # regressions must not enter the collapsed loadings
      this_beta <- mlist$beta
      this_beta[is.na(this_beta)] <- 0
      keep_row_idx <- match(lv_ind_names, lv_names)
      bstar <- matrix(0, nrow = nrow(this_beta), ncol = ncol(this_beta))
      bstar[keep_row_idx, ] <- this_beta[keep_row_idx, , drop = FALSE]
      ib_inv <- solve(diag(nrow(bstar)) - bstar)
      lambda_star_b <- mm_lambda[[b]] %*% ib_inv
      # reach[r, j] = TRUE: factor j influences factor r through the
      # measurement chain (r is a latent indicator of j, possibly indirectly)
      reach <- abs(ib_inv) > 1e-12
      diag(reach) <- FALSE
      # the factors we can fix: empty column in LAMBDA, but non-empty
      # column in LAMBDA.star
      fix_idx <- empty_idx[apply(lambda_star_b[, empty_idx, drop = FALSE],
                                 2L, function(x) any(x != 0))]
      if (length(fix_idx) == 0L) next
      fixes <- list()
      for (j in fix_idx) {
        # this factor's layer: neither ancestors nor descendants of j, ...
        allowed <- which(!reach[j, ] & !reach[, j])
        # ... keeping only the topmost factors (no ancestor within the set)
        tcols <- allowed[!apply(reach[allowed, allowed, drop = FALSE],
                                1L, any)]
        key <- paste(tcols, collapse = "-")
        if (is.null(fixes[[key]])) {
          fixes[[key]] <- list(tcols = tcols, cols = j)
        } else {
          fixes[[key]]$cols <- c(fixes[[key]]$cols, j)
        }
      }
      mm_lambda_star[[b]] <- lambda_star_b
      ho_fixes[[b]] <- fixes
    }
  }

  fs <- vector("list", length = lavdata@ngroups)
  if (fsm) {
    fsm_1 <- vector("list", length = lavdata@ngroups)
  }
  if (rel) {
    rel_1 <- vector("list", length = lavdata@ngroups)
  }
  if (transform) {
    # correlation-preserving transformation: Green (regression) or
    # Krijnen/McDonald (Bartlett) determinacy matrix
    if (bartlett) {
      tmat <- lav_predict_tmat_det(lavmodel = lavmodel, lavimplied = lavimplied)
    } else {
      tmat <- lav_predict_tmat_green(lavmodel = lavmodel,
                                     lavimplied = lavimplied)
    }
  }

  if (acov != "none") {
    se <- acov # acov implies SE
  }
  if (se != "none") {
    se_1 <- vector("list", length = lavdata@ngroups)
    # return full sampling covariance matrix?
    if (acov != "none") {
      acov_1 <- vector("list", length = lavdata@ngroups)
    }
  }

  for (g in 1:lavdata@ngroups) {
    # determine block
    b <- lav_predict_block_idx(lavdata, g, level)

    veta_g <- veta[[b]]
    eeta_g <- eeta[[b]]
    lambda_g <- mm_lambda[[b]]
    ey_g <- ey[[b]]
    sigma_inv_g <- sigma_inv[[b]]
    # higher-order factors (Bartlett only; NULL otherwise)
    lambda_star_g <- mm_lambda_star[[b]]
    ho_fixes_g <- ho_fixes[[b]]

    if (lavdata@nlevels > 1L) {
      lp <- lavdata@Lp[[g]]

      # cluster means (Y2); recompute from data_obs (rather than the
      # sample statistics) so that the same code path handles newdata,
      # where the cluster structure comes from lavdata (= the newdata
      # object; see lav_predict_internal)
      y2_g <- rowsum.default(data_obs[[g]],
        group = lp$cluster.idx[[2]], reorder = FALSE, na.rm = FALSE
      ) / lp$cluster.size[[2]]

      # implied for this group
      group_idx <- (g - 1) * lavdata@nlevels + seq_len(lavdata@nlevels)
      implied_group <- lapply(lavimplied, function(x) x[group_idx])

      # random effects (=random intercepts or cluster means)
      if (lavdata@missing %in% c("ml", "ml.x")) {
        # missing data: use the posterior means from the E-step,
        # conditioning on ALL the observed data (new in 0.7-1); the
        # missing values are replaced by their posterior means, so
        # that the (linear) factor score computation below yields the
        # exact posterior means E(eta | all observed data)
        mb_j <- lav_mvn_cl_mi_estep_ranef(
          y1 = data_obs[[g]], y2 = y2_g,
          lp = lp, mp = lavdata@Mp[[g]],
          mu_w = implied_group$mean[[1]],
          sigma_w = implied_group$cov[[1]],
          mu_b = implied_group$mean[[2]],
          sigma_b = implied_group$cov[[2]],
          se = FALSE, impute = TRUE
        )
      } else {
        out <- lav_mvn_cl_implied22l(
          lp = lp,
          implied = implied_group
        )
        mb_j <- lav_mvn_cl_em_estep_ranef(
          ylp = list(list(), list(Y2 = y2_g)), lp = lp,
          sigma_w = out$sigma.w, sigma_b = out$sigma.b,
          sigma_zz = out$sigma.zz, sigma_yz = out$sigma.yz,
          mu_z = out$mu.z, mu_w = out$mu.w, mu_b = out$mu.b,
          se = FALSE
        )
      }

      ov_idx <- lp$ov.idx

      if (level == 1L) {
        data_w <- data_obs[[g]][, ov_idx[[1]]]
        # with missing = "ml": use the (posterior-mean) completed data
        if (!is.null(attr(mb_j, "y1w.imputed"))) {
          data_w <- attr(mb_j, "y1w.imputed")
        }
        data_b <- mb_j[lp$cluster.idx[[2]], , drop = FALSE]

        # center
        data_obs_g <- data_w - data_b
      } else if (level == 2L) {
        data_b_1 <- matrix(0,
          nrow = nrow(mb_j),
          ncol = ncol(data_obs[[g]])
        )
        data_b_1[, ov_idx[[1]]] <- mb_j
        between_idx <- lp$between.idx[[2 * g]]
        if (length(between_idx) > 0L) {
          data_b_1[, between_idx] <- data_obs[[g]][
            !duplicated(lp$cluster.idx[[2]]),
            between_idx
          ]
          # with missing = "ml": use the (posterior-mean) completed
          # between-only variables
          if (!is.null(attr(mb_j, "z.imputed"))) {
            data_b_1[, between_idx] <- attr(mb_j, "z.imputed")
          }
        }
        data_obs_g <- data_b_1[, ov_idx[[2]]]
      } else {
        lav_msg_stop(gettext("only 2 levels are supported"))
      }
    } else {
      data_obs_g <- data_obs[[b]]
    }

    nfac <- ncol(veta_g)
    if (nfac == 0L) {
      fs[[g]] <- matrix(0, lavdata@nobs[[g]], nfac)
      next
    }

    # center data
    yc <- t(t(data_obs_g) - ey_g)

    # NOTE: do NOT scale yc by the sampling weights here. Factor scores are
    # per-case predictions; a case's score must not depend on its own
    # sampling weight. Weights already enter via the (weighted) parameter
    # estimates and the weighted mean (ey_g) used for centering above.

    # global factor score coefficient matrix 'C'
    fsc <- compute_fsc(lambda_g, sigma_inv_g, veta_g,
      lambda_star = lambda_star_g, fixes = ho_fixes_g
    )

    # transform?
    if (transform) {
      fsc <- tmat[[b]] %*% fsc
    }

    # store fsm?
    if (fsm) {
      fsm_g <- fsc
    }

    # reliability?
    if (rel) {
      var_f <- fsc %*% sigma_hat[[g]] %*% t(fsc)  # or S?
      cov_f_eta <- fsc %*% lambda_g %*% veta_g
      var_eta <- veta_g
      # FS.determinacy <- diag( diag(1/sqrt(diag(Var.f))) %*%
      #                         Cov.f.eta %*%
      #                         diag(1/sqrt(diag(Var.eta)))
      #                       )
      fs_determinacy <- (diag(cov_f_eta) /
                          (sqrt(diag(var_f)) * sqrt(diag(var_eta))))
      rel_g <- fs_determinacy * fs_determinacy
    }

    # compute factor scores
    # (for two-level models with missing = "ml", the data have been
    # completed by their posterior means above, so the complete-data
    # expressions apply)
    if (lavdata@missing %in% c("ml", "ml.x") && lavdata@nlevels == 1L) {
      # missing patterns for this group
      mp <- mp_1[[g]]

      # factor scores container
      fs_g <- matrix(as.numeric(NA), nrow(yc), ncol = length(eeta_g))

      if (se == "standard") {
        se_g <- matrix(as.numeric(NA), nrow(yc), ncol = length(eeta_g))
      }

      if (acov == "standard") {
        acov_g <- vector("list", length = mp$npatterns)
      }

      # compute FSC per pattern
      for (p in seq_len(mp$npatterns)) {
        var_idx <- mp$pat[p, ] # observed
        na_idx <- which(!var_idx) # missing

        # extract observed data for these (centered) cases
        oc <- yc[mp$case.idx[[p]], mp$pat[p, ], drop = FALSE]

        # invert sigma (Sigma_22, observed part only) for this pattern
        sigma_22_inv <- try(lav_mat_sym_inverse_update(
          s_inv = sigma_inv_g, rm_idx = na_idx, logdet = FALSE
        ), silent = TRUE)
        if (inherits(sigma_22_inv, "try-error")) {
          lav_msg_stop(gettext("Sigma_22.inv cannot be inverted"))
        }

        lambda <- lambda_g[var_idx, , drop = FALSE]
        lambda_star_p <- NULL
        if (!is.null(lambda_star_g)) {
          lambda_star_p <- lambda_star_g[var_idx, , drop = FALSE]
        }
        fsc <- compute_fsc(lambda, sigma_22_inv, veta_g,
          lambda_star = lambda_star_p, fixes = ho_fixes_g
        )

        # if FSC contains rows that are all-zero, replace by NA
        #
        # this happens eg if all the indicators of a single factor
        # are missing; then this column in lambda only contains zeroes
        # and therefore the corresponding row in FSC contains only
        # zeroes, leading to factor score 0
        #
        # showing 'NA' is better than getting 0
        #
        # (Note that this is not needed for the 'regression' method,
        #  only for Bartlett)
        if (bartlett) {
          fsc_nona <- fsc # keep a NA-free copy for the SEs below
          zero_idx <- which(apply(fsc, 1L, function(x) all(x == 0)))
          if (length(zero_idx) > 0L) {
            fsc[zero_idx, ] <- NA
          }
        }

        # factor scores for this pattern
        fs_g[mp$case.idx[[p]], ] <- t(fsc %*% t(oc) + eeta_g)

        # SE?
        if (se == "standard") {
          if (is.null(lambda_star_p)) {
            tmp <- compute_fscov(lambda, sigma_22_inv, veta_g)
          } else {
            tmp <- compute_fscov(lambda, sigma_22_inv, veta_g,
              fsc = fsc_nona, sigma = sigma_hat[[b]][var_idx, var_idx,
                                                     drop = FALSE]
            )
          }
          tmp_d <- diag(tmp)
          tmp_d[tmp_d < 1e-05] <- as.numeric(NA)

          # all cases in this pattern get the same SEs
          se_g[mp$case.idx[[p]], ] <- matrix(sqrt(tmp_d),
            nrow = length(mp$case.idx[[p]]),
            ncol = ncol(se_g), byrow = TRUE
          )
        }

        # acov?
        if (acov == "standard") {
          acov_g[[p]] <- tmp # for this pattern
        }
      } # p
    } else {
      # compute factor scores
      fs_g <- t(fsc %*% t(yc) + eeta_g)
    }

    # replace values in dummy lv's by their observed counterpart
    if (length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L && level == 1L) {
      fs_g[, lavmodel@ov.y.dummy.lv.idx[[g]]] <-
        data_obs_g[, lavmodel@ov.y.dummy.ov.idx[[g]], drop = FALSE]
    }
    if (length(lavmodel@ov.x.dummy.lv.idx[[g]]) > 0L && level == 1L) {
      fs_g[, lavmodel@ov.x.dummy.lv.idx[[g]]] <-
        data_obs_g[, lavmodel@ov.x.dummy.ov.idx[[g]], drop = FALSE]
    }

    fs[[g]] <- fs_g

    # FSM
    if (fsm) {
      fsm_1[[g]] <- fsm_g
    }

    # REL
    if (rel) {
      rel_1[[g]] <- rel_g
    }

    # standard error
    if (se == "standard" && !transform) { # for now
      if (lavdata@missing %in% c("ml", "ml.x") && lavdata@nlevels == 1L) {
        se_1[[g]] <- se_g
        if (acov == "standard") {
          acov_1[[g]] <- acov_g
        }
      } else { # complete data
        if (is.null(lambda_star_g)) {
          tmp <- compute_fscov(lambda_g, sigma_inv_g, veta_g)
        } else {
          tmp <- compute_fscov(lambda_g, sigma_inv_g, veta_g,
            fsc = fsc, sigma = sigma_hat[[b]]
          )
        }
        tmp_d <- diag(tmp)
        tmp_d[tmp_d < 1e-05] <- as.numeric(NA)
        se_1[[g]] <- matrix(sqrt(tmp_d), nrow = 1L)

        # return full sampling covariance matrix?
        if (acov == "standard") {
          acov_1[[g]] <- tmp
        }
      }
    } # se = "standard"
  } # g

  if (fsm) {
    attr(fs, "fsm") <- fsm_1
  }
  if (rel) {
    attr(fs, "rel") <- rel_1
  }
  if (se != "none") {
    attr(fs, "se") <- se_1
    # return full sampling covariance matrix?
    if (acov == "standard") {
      attr(fs, "acov") <- acov_1
    }
  }

  fs
}

# decide whether (and how) to parallelise the per-case factor-score
# optimisation for a given number of independent tasks 'n_tasks'.
#  - "no"        : always serial
#  - "multicore" : fork-based (unix only), as requested by the user
#  - "snow"      : PSOCK cluster, as requested by the user
#  - "auto"      : fork-based, but only when it is clearly worthwhile, i.e.
#                  enough independent optimisations to amortise the overhead
#                  (the dedup step often leaves only a handful of unique cases)
# returns one of "no", "multicore", "snow".
lav_predict_parallel_mode <- function(parallel = "no", ncpus = NULL,
                                      n_tasks = 1L) {
  # number of unique optimisations above which "auto" switches on forking
  auto_threshold <- 1000L

  if (is.null(ncpus) || ncpus < 2L || n_tasks < 2L || parallel == "no") {
    return("no")
  }
  unix <- .Platform$OS.type != "windows"
  if (parallel == "multicore") {
    if (unix) "multicore" else "no"
  } else if (parallel == "snow") {
    "snow"
  } else if (parallel == "auto") {
    # auto only uses fork, and only when there is enough work to gain from it
    if (unix && n_tasks >= auto_threshold) "multicore" else "no"
  } else {
    "no"
  }
}

# factor scores - EBM or ML
lav_predict_eta_ebm_ml <- function(lavobject = NULL, # for convenience
                                   # sub objects
                                   lavmodel = NULL, lavdata = NULL,
                                   lavsamplestats = NULL,
                                   # optional new data
                                   data_obs = NULL, exo = NULL,
                                   se = "none", acov = "none", level = 1L,
                                   ml = FALSE,
                                   optim_method = "bfgs",
                                   parallel = "no", ncpus = NULL, cl = NULL) {
  optim_method <- tolower(optim_method)

  stopifnot(optim_method %in% c("nlminb", "bfgs"))

  # parallelisation: the per-case optimisation below is embarrassingly parallel
  # and fully deterministic (no RNG handling needed). Resolve the number of
  # cores now; the actual serial/parallel decision is made per group, once the
  # number of unique cases to optimise is known (see lav_predict_parallel_mode).
  parallel <- parallel[1L]
  if (parallel != "no" && is.null(ncpus)) {
    ncpus <- max(1L, parallel::detectCores() - 2L, na.rm = TRUE)
  }
  if (parallel != "no") {
    loadNamespace("parallel")
  }

  ### FIXME: if all indicators of a factor are normal, can we not
  ###        just use the `classic' regression method??
  ###        (perhaps after whitening, to get uncorrelated factors...)

  # full object?
  if (inherits(lavobject, "lavaan")) {
    lavmodel <- lavobject@Model
    lavdata <- lavobject@Data
    lavsamplestats <- lavobject@SampleStats
  } else {
    stopifnot(
      !is.null(lavmodel), !is.null(lavdata),
      !is.null(lavsamplestats)
    )
  }

  # new data?
  if (is.null(data_obs)) {
    data_obs <- lavdata@X
  }
  if (is.null(exo)) {
    exo <- lavdata@eXo
  }

  # se?
  if (acov != "none") {
    se <- acov # acov implies SE
  }
  # For categorical data the factor scores are posterior-mode (EBM) or ML
  # estimates obtained by numerical optimisation. Their standard errors are
  # based on the curvature of the objective at the optimum: the (co)variance
  # matrix is the inverse of the Hessian of the negative log-posterior (EBM) or
  # negative log-likelihood (ML), i.e. a Laplace/observed-information
  # approximation. Unlike the linear-Gaussian case, this curvature depends on
  # the response pattern, so the standard errors differ from one (pattern of)
  # observation(s) to the next.
  se_flag <- (se != "none")
  acov_flag <- (acov != "none")
  if (se_flag) {
    se_1 <- vector("list", length = lavdata@ngroups)
    if (acov_flag) {
      acov_1 <- vector("list", length = lavdata@ngroups)
    }
  }

  vetax <- lav_model_vetax(lavmodel = lavmodel)
  vetax_inv <- vetax
  for (g in seq_len(lavdata@ngroups)) {
    if (nrow(vetax[[g]]) > 0L) {
      vetax_inv[[g]] <- solve(vetax[[g]])
    }
  }
  eetax <- lav_model_eetax(
    lavmodel = lavmodel, lavsamplestats = lavsamplestats,
    exo = exo, nobs = lapply(data_obs, NROW),
    remove_dummy_lv = TRUE
  ) ## FIXME?
  th_1 <- lav_model_th(lavmodel = lavmodel, delta = FALSE)
  mm_theta <- lav_model_theta(lavmodel = lavmodel)

  # check for zero entries in THETA (new in 0.6-4)
  for (g in seq_len(lavdata@ngroups)) {
    if (any(diag(mm_theta[[g]]) == 0)) {
      lav_msg_stop(gettext(
        "(residual) variance matrix THETA contains zero elements
        on the diagonal."))
    }
  }

  # Model matrices per group, used by the conditional-density evaluation. They
  # are loop-invariant within a group, so we extract them once here instead of
  # rebuilding them on every objective evaluation (the objective is called many
  # times per case by the optimiser).
  fy_mlist_g <- lapply(seq_len(lavdata@ngroups), function(gg) {
    mm_idx <- 1:lavmodel@nmat[gg] + cumsum(c(0, lavmodel@nmat))[gg]
    lavmodel@GLIST[mm_idx]
  })

  # local objective function: x = lv values. The per-case dummy ov.y values are
  # passed in (rather than read via a case index 'i') so the objective does not
  # depend on the enclosing loop variable -- this lets us evaluate cases via
  # lapply()/mclapply() instead of a for() loop.
  f_eta_i <- function(x, y_i, x_i, mu_i, dummy_y) {
    # add 'dummy' values (if any) for ov.y
    if (length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L) {
      x2 <- c(x - mu_i, dummy_y)
    } else {
      x2 <- x - mu_i
    }

    # conditional density of y, given eta.i(=x)
    log_fy <- lav_predict_fy_eta_i(
      lavmodel = lavmodel,
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      y_i = y_i,
      x_i = x_i,
      eta_i = matrix(x2, nrow = 1L), # <---- eta!
      theta_sd = theta_sd,
      th = th,
      th_idx = th_idx,
      g = g,
      log = TRUE,
      mlist = fy_mlist_g[[g]]
    )

    if (ml) {
      # NOTE: 'true' ML is simply  -1*sum(log.fy)
      #  - but there is no upper/lower bound for the extrema:
      #    a pattern of all (in)correct drives the 'theta' parameter
      #    towards +/- Inf
      # - therefore, we add a vague prior, just to stabilize
      #
      diff <- t(x) - mu_i
      v <- diag(length(x)) * 1e-05
      tmp <- as.numeric(0.5 * diff %*% v %*% t(diff))
      out <- 1 + tmp - sum(log_fy, na.rm = TRUE)
    } else {
      diff <- t(x) - mu_i
      v <- vetax_inv[[g]]
      tmp <- as.numeric(0.5 * diff %*% v %*% t(diff))
      out <- tmp - sum(log_fy, na.rm = TRUE)
    }
    out
  }

  fs <- vector("list", length = lavdata@ngroups)
  for (g in seq_len(lavdata@ngroups)) {
    nfac <- ncol(vetax[[g]])
    nfac2 <- nfac
    if (length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L) {
      nfac2 <- nfac2 + length(lavmodel@ov.y.dummy.lv.idx[[g]])
    }
    fs[[g]] <- matrix(as.numeric(NA), nrow(data_obs[[g]]), nfac2)

    # special case: no regular lv's
    if (nfac == 0) {
      # impute dummy ov.y (if any)
      fs[[g]][, lavmodel@ov.y.dummy.ov.idx[[g]]] <-
        data_obs[[g]][, lavmodel@ov.y.dummy.ov.idx[[g]], drop = FALSE]
      next
    }

    ## FIXME: factor scores not identical (but close) to Mplus
    #         if delta elements not equal to 1??
    # mm_in_group <- 1:lavmodel@nmat[g] + cumsum(c(0, lavmodel@nmat))[g]
    # mlist <- lavmodel@GLIST[mm_in_group]

    # check for negative values
    neg_var_idx <- which(diag(mm_theta[[g]]) < 0)
    if (length(neg_var_idx) > 0) {
      lav_msg_warn(
        gettext("factor scores could not be computed due to at least
                one negative (residual) variance"))
      next
    }

    # common values
    theta_sd <- sqrt(diag(mm_theta[[g]]))
    th <- th_1[[g]]
    th_idx <- lavmodel@th.idx[[g]]

    # Casewise optimisation. Two cases that share the same inputs -- the data
    # row y_i, the exogenous covariates x_i, and the conditional mean mu_i --
    # define exactly the same objective function and starting value; since the
    # optimiser is deterministic, they yield identical factor scores. We
    # therefore optimise once per unique combination and copy the result to all
    # matching cases. For categorical data many cases share the same response
    # pattern, so this can reduce the number of optimisations by an order of
    # magnitude.
    n <- nrow(data_obs[[g]])

    # key each case by everything the objective depends on, using the exact
    # (round-tripping) hexadecimal representation so that only bitwise-identical
    # rows are ever merged
    key_cols <- list(data_obs[[g]], eetax[[g]])
    if (!is.null(exo[[g]])) {
      key_cols <- c(key_cols, list(exo[[g]]))
    }
    key_mat <- do.call(cbind, key_cols)
    case_key <- do.call(paste, c(
      lapply(seq_len(ncol(key_mat)), function(j) sprintf("%a", key_mat[, j])),
      sep = "\r"
    ))
    rep_idx <- which(!duplicated(case_key)) # one representative case / pattern
    case_map <- rep_idx[match(case_key, case_key[rep_idx])] # case -> rep row

    # optimise the latent scores for a single (representative) case 'i'.
    # Returns a list with the score 'eta' (incl. any dummy ov.y values) and,
    # when requested, the standard errors 'se' and sampling covariance 'acov'
    # of the (regular) latent scores.
    solve_i <- function(i) {
      # eXo?
      if (!is.null(exo[[g]])) {
        x_i <- exo[[g]][i, , drop = FALSE]
      } else {
        x_i <- NULL
      }
      mu_i <- eetax[[g]][i, , drop = FALSE]
      y_i <- data_obs[[g]][i, , drop = FALSE]

      # dummy ov.y values for this case (passed to the objective and appended
      # to the score below)
      if (length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L) {
        dummy_y <- data_obs[[g]][i, lavmodel@ov.y.dummy.ov.idx[[g]],
          drop = FALSE
        ]
      } else {
        dummy_y <- NULL
      }

      start_1 <- numeric(nfac) # initial values for eta
      se_i <- rep(as.numeric(NA), nfac)
      acov_i <- matrix(as.numeric(NA), nfac, nfac)

      if (!all(is.na(y_i))) {
        # find best values for eta.i
        if (optim_method == "nlminb") {
          out <- nlminb(
            start = start_1, objective = f_eta_i,
            gradient = NULL, # for now
            control = list(rel.tol = 1e-8),
            y_i = y_i, x_i = x_i, mu_i = mu_i, dummy_y = dummy_y
          )
        } else if (optim_method == "bfgs") {
          out <- optim(
            par = start_1, fn = f_eta_i,
            gr = NULL,
            control = list(reltol = 1e-8, fnscale = 1.1),
            method = "BFGS",
            y_i = y_i, x_i = x_i, mu_i = mu_i, dummy_y = dummy_y
          )
        }
        if (out$convergence == 0L) {
          eta_reg <- out$par

          # standard errors from the curvature of the objective at the
          # optimum: f_eta_i is the negative log-posterior (EBM) or negative
          # log-likelihood (ML), so its Hessian at the mode is the inverse of
          # the (Laplace) sampling covariance of the latent scores.
          if (se_flag) {
            hess <- try(stats::optimHess(eta_reg, f_eta_i,
              y_i = y_i, x_i = x_i, mu_i = mu_i, dummy_y = dummy_y
            ), silent = TRUE)
            if (!inherits(hess, "try-error")) {
              tmp_acov <- try(solve(hess), silent = TRUE)
              if (!inherits(tmp_acov, "try-error")) {
                acov_i <- tmp_acov
                tmp_d <- diag(tmp_acov)
                tmp_d[tmp_d < 0] <- as.numeric(NA) # catch non-pd Hessians
                se_i <- sqrt(tmp_d)
              }
            }
          }
        } else {
          eta_reg <- rep(as.numeric(NA), nfac)
        }
      } else {
        eta_reg <- rep(as.numeric(NA), nfac)
      }

      # add dummy ov.y lv values
      eta_i <- eta_reg
      if (length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L) {
        eta_i <- c(eta_i, dummy_y)
      }

      list(eta = eta_i, se = se_i, acov = acov_i)
    }

    # optimise once per unique case, serially or in parallel
    run_mode <- lav_predict_parallel_mode(parallel, ncpus, length(rep_idx))
    if (run_mode == "multicore") {
      res_list <- parallel::mclapply(rep_idx, solve_i,
        mc.cores = min(ncpus, length(rep_idx))
      )
    } else if (run_mode == "snow") {
      if (is.null(cl)) {
        cl0 <- parallel::makePSOCKcluster(
          rep("localhost", min(ncpus, length(rep_idx)))
        )
        on.exit(parallel::stopCluster(cl0), add = TRUE)
        res_list <- parallel::parLapply(cl0, rep_idx, solve_i)
      } else {
        res_list <- parallel::parLapply(cl, rep_idx, solve_i)
      }
    } else {
      res_list <- lapply(rep_idx, solve_i)
    }
    for (u in seq_along(rep_idx)) {
      fs[[g]][rep_idx[u], ] <- res_list[[u]]$eta
    }

    # copy each representative's factor scores to all cases in its pattern
    fs[[g]] <- fs[[g]][case_map, , drop = FALSE]

    # expand the (per-pattern) standard errors / sampling covariances to all
    # cases. The se/acov refer to the nfac regular factors; they are placed in
    # the leading columns of an nfac2-wide structure (dummy ov.y columns left
    # as NA) so that the dummy-lv handling in lav_predict_internal() applies
    # identically to the scores and to their standard errors.
    if (se_flag) {
      u_of_case <- match(case_key, case_key[rep_idx]) # case -> unique index
      se_rep <- do.call(rbind, lapply(res_list, "[[", "se")) # n_unique x nfac
      se_g <- matrix(as.numeric(NA), n, nfac2)
      se_g[, seq_len(nfac)] <- se_rep[u_of_case, , drop = FALSE]
      se_1[[g]] <- se_g
      if (acov_flag) {
        acov_rep <- lapply(res_list, "[[", "acov") # list: nfac x nfac
        acov_1[[g]] <- lapply(u_of_case, function(u) {
          m <- matrix(as.numeric(NA), nfac2, nfac2)
          m[seq_len(nfac), seq_len(nfac)] <- acov_rep[[u]]
          m
        })
      }
    }
  }

  if (se_flag) {
    attr(fs, "se") <- se_1
    if (acov_flag) {
      attr(fs, "acov") <- acov_1
    }
  }

  fs
}

# predicted value for response y*_i, conditional on the predicted latent
# variable scores
# `measurement part':
#     y*_i = nu + lambda eta_i + K x_i + epsilon_i
#
#    where eta_i = latent variable value for i (either given or from predict)
#
# Two types: 1) nrow(eta) = nrow(X) (factor scores)
#            2) nrow(eta) = 1L (given values)
#
# in both cases, we return [nobs x nvar] matrix per group
lav_predict_yhat <- function(lavobject = NULL, # for convience
                             # sub objects
                             lavmodel = NULL, lavdata = NULL,
                             lavsamplestats = NULL,
                             lavimplied = NULL,
                             # new data
                             data_obs = NULL, exo = NULL,
                             # eta values
                             eta = NULL,
                             # options
                             method = "EBM",
                             duplicate = FALSE,
                             optim_method = "bfgs",
                             fsm = FALSE,
                             resid_flag = FALSE) {
  # full object?
  if (inherits(lavobject, "lavaan")) {
    lavmodel <- lavobject@Model
    lavdata <- lavobject@Data
    lavsamplestats <- lavobject@SampleStats
    lavimplied <- lavobject@implied
  } else {
    stopifnot(
      !is.null(lavmodel), !is.null(lavdata),
      !is.null(lavsamplestats), !is.null(lavimplied)
    )
  }

  # new data?
  if (is.null(data_obs)) {
    data_obs <- lavdata@X
  }
  if (is.null(exo)) {
    exo <- lavdata@eXo
  }

  # do we get values for eta? If not, use `predict' to get plausible values
  if (is.null(eta)) {
    eta <- lav_predict_eta(
      lavobject = NULL, lavmodel = lavmodel,
      lavdata = lavdata, lavsamplestats = lavsamplestats,
      lavimplied = lavimplied,
      data_obs = data_obs, exo = exo, method = method,
      optim_method = optim_method, fsm = fsm
    )
    fsm_1 <- attr(eta, "fsm")
  } else {
    # matrix
    if (is.matrix(eta)) { # user-specified?
      if (nrow(eta) == 1L) {
        tmp <- matrix(eta, lavsamplestats@ntotal, length(eta),
          byrow = TRUE
        )
      } else if (nrow(eta) != lavsamplestats@ntotal) {
        lav_msg_stop(gettext("nrow(eta) != lavsamplestats@ntotal"))
      } else {
        tmp <- eta
      }
      eta <- lapply(1:lavdata@ngroups, function(i) tmp[lavdata@case.idx[[i]], ])
      # vector: just 1 row of factor-scores
    } else if (is.numeric(eta)) {
      # convert to matrix
      tmp <- matrix(eta, lavsamplestats@ntotal, length(eta), byrow = TRUE)
      eta <- lapply(1:lavdata@ngroups, function(i) tmp[lavdata@case.idx[[i]], ])
    } else if (is.list(eta)) {
      stopifnot(lavdata@ngroups == length(eta))
    }
  }

  yhat <- lav_model_yhat(
    lavmodel = lavmodel, glist = NULL,
    lavsamplestats = lavsamplestats, exo = exo,
    nobs = lapply(data_obs, NROW),
    eta = eta, duplicate = duplicate
  )

  # if conditional.x, paste exo
  if (lavmodel@categorical && !is.null(exo)) {
    yhat <- lapply(seq_len(lavdata@ngroups), function(g) {
      ret <- cbind(yhat[[g]], exo[[g]])
      ret
    })
  }

  # residuals? compute y - yhat
  if (resid_flag) {
    res <- lapply(seq_len(lavdata@ngroups), function(g) {
      ret <- data_obs[[g]] - yhat[[g]]
      ret
    })
  } else {
    res <- yhat
  }

  # fsm?
  if (fsm) {
    attr(res, "fsm") <- fsm_1
  }

  res
}

# conditional density y -- assuming independence!!
# f(y_i | eta_i, x_i) for EACH item
#
lav_predict_fy <- function(lavobject = NULL, # for convience
                           # sub objects
                           lavmodel = NULL, lavdata = NULL,
                           lavsamplestats = NULL,
                           lavimplied = NULL,
                           # new data
                           data_obs = NULL, exo = NULL,
                           # eta values
                           eta = NULL,
                           # options
                           method = "EBM",
                           log_1 = FALSE,
                           optim_method = "bfgs") {
  # full object?
  if (inherits(lavobject, "lavaan")) {
    lavmodel <- lavobject@Model
    lavdata <- lavobject@Data
    lavsamplestats <- lavobject@SampleStats
    lavimplied <- lavobject@implied
  } else {
    stopifnot(
      !is.null(lavmodel), !is.null(lavdata),
      !is.null(lavsamplestats), !is.null(lavimplied)
    )
  }

  # new data?
  if (is.null(data_obs)) {
    data_obs <- lavdata@X
  }
  if (is.null(exo)) {
    exo <- lavdata@eXo
  }

  # we need the YHATs (per group)
  yhat <- lav_predict_yhat(
    lavobject = NULL, lavmodel = lavmodel,
    lavdata = lavdata, lavsamplestats = lavsamplestats,
    lavimplied = lavimplied,
    data_obs = data_obs, exo = exo, eta = eta, method = method,
    duplicate = FALSE, optim_method = optim_method
  )

  mm_theta <- lav_model_theta(lavmodel = lavmodel)
  th <- lav_model_th(lavmodel = lavmodel, delta = FALSE)

  fy <- vector("list", length = lavdata@ngroups)
  for (g in seq_len(lavdata@ngroups)) {
    fy[[g]] <- lav_predict_fy_internal(
      x = data_obs[[g]], yhat = yhat[[g]],
      th = th[[g]], mm_theta = mm_theta[[g]],
      num_idx = lavmodel@num.idx[[g]],
      th_idx = lavmodel@th.idx[[g]],
      link = lavmodel@link, log_1 = log_1
    )
  }

  fy
}


# single group, internal function
lav_predict_fy_internal <- function(x = NULL, yhat = NULL,
                                    th = NULL, mm_theta = NULL,
                                    num_idx = NULL, th_idx = NULL,
                                    link = NULL, log_1 = FALSE) {
  # shortcuts
  theta_var <- diag(mm_theta)

  # check size yhat (either 1L or Nobs rows)
  if (!(nrow(yhat) == 1L || nrow(yhat) == nrow(x))) {
    lav_msg_stop(gettext("nrow(yhat[[g]]) not 1L and not nrow(X))"))
  }

  fy_group <- matrix(0, nrow(x), ncol(x))
  # if(NORMAL) {
  #    if(nrow(yhat) == nrow(X)) {
  #        tmp <- (X - yhat)^2
  #    } else {
  #        tmp <- sweep(X, MARGIN=2, STATS=yhat, FUN="-")^2
  #    }
  #    tmp1 <- sweep(tmp, MARGIN=2, theta.var, "/")
  #    tmp2 <- exp( -0.5 * tmp1 )
  #    tmp3 <- sweep(tmp2, MARGIN=2, sqrt(2*pi*theta.var), "/")
  #    if(log.) {
  #        FY.group <- log(tmp3)
  #    } else {
  #        FY.group <- tmp3
  #    }
  # } else {
  # mixed items

  ord_idx <- unique(th_idx[th_idx > 0L])

  # first, NUMERIC variables
  if (length(num_idx) > 0L) {
    for (v in num_idx) {
      fy_group[, v] <- dnorm(x[, v],
        # yhat may change or not per case
        mean = yhat[, v],
        sd   = sqrt(theta_var[v]),
        log  = log_1
      )
    }
  }

  # second, ORDERED variables
  for (v in ord_idx) {
    th_y <- th[th_idx == v]
    th_y_1 <- c(-Inf, th_y, Inf)
    ncat <- length(th_y) + 1L
    fy <- numeric(ncat)
    theta_v <- sqrt(theta_var[v])
    yhat_v <- yhat[, v]

    # two cases: yhat.v is a scalar, or has length = nobs
    fy <- matrix(0, nrow = length(yhat_v), ncol = ncat)

    # for each category
    for (k in seq_len(ncat)) {
      if (link == "probit") {
        fy[, k] <- pnorm((th_y_1[k + 1] - yhat_v) / theta_v) -
          pnorm((th_y_1[k] - yhat_v) / theta_v)
      } else if (link == "logit") {
        fy[, k] <- plogis((th_y_1[k + 1] - yhat_v) / theta_v) -
          plogis((th_y_1[k] - yhat_v) / theta_v)
      } else {
        lav_msg_stop(gettext("link must be probit or logit"))
      }
    }

    # underflow
    idx <- which(fy < .Machine$double.eps)
    if (length(idx) > 0L) {
      fy[idx] <- .Machine$double.eps
    }

    # log?
    if (log_1) {
      fy <- log(fy)
    }

    # case-wise expansion/selection
    if (length(yhat_v) == 1L) {
      # expand category probabilities for all observations
      fy_group[, v] <- fy[1L, x[, v]]
    } else {
      # select correct category probability per observation
      fy_group[, v] <- fy[cbind(seq_len(nrow(fy)), x[, v])]
    }
  } # ord

  fy_group
}



# conditional density y -- assuming independence!!
# f(y_i | eta_i, x_i)
#
# but for a SINGLE observation y_i (and x_i), for given values of eta_i
#
lav_predict_fy_eta_i <- function(lavmodel = NULL, lavdata = NULL,
                                 lavsamplestats = NULL,
                                 y_i = NULL, x_i = NULL,
                                 eta_i = NULL, theta_sd = NULL, g = 1L,
                                 th = NULL, th_idx = NULL, log = TRUE,
                                 mlist = NULL) {
  # model matrices for this group; loop-invariant, so callers in a tight loop
  # may pass a precomputed 'mlist' to avoid rebuilding it on every call
  if (is.null(mlist)) {
    mm_in_group <- 1:lavmodel@nmat[g] + cumsum(c(0, lavmodel@nmat))[g]
    mlist <- lavmodel@GLIST[mm_in_group]
  }

  # linear predictor for all items
  yhat <-
    lav_lisrel_eyetax(
      mlist = mlist,
      exo = x_i,
      eta = eta_i,
      sample_mean = lavsamplestats@mean[[g]],
      ov_y_dummy_ov_idx = lavmodel@ov.y.dummy.ov.idx[[g]],
      ov_x_dummy_ov_idx = lavmodel@ov.x.dummy.ov.idx[[g]],
      ov_y_dummy_lv_idx = lavmodel@ov.y.dummy.lv.idx[[g]],
      ov_x_dummy_lv_idx = lavmodel@ov.x.dummy.lv.idx[[g]],
      delta = FALSE
    )

  # P(y_i | eta_i, x_i) for all items
  if (all(lavdata@ov$type == "numeric")) {
    # NORMAL case
    fy <- dnorm(y_i, mean = yhat, sd = theta_sd, log = log)
  } else {
    nvar_g <- lavmodel@nvar[g]
    var_type <- lavdata@ov$type[seq_len(nvar_g)]
    fy <- numeric(nvar_g)

    # numeric variables: vectorised normal density
    num_idx <- which(var_type == "numeric")
    if (length(num_idx) > 0L) {
      fy[num_idx] <- dnorm(y_i[num_idx],
        mean = yhat[num_idx], sd = theta_sd[num_idx], log = log
      )
    }

    # ordered variables: P(y = k) = Phi((tau_k - yhat)/sd) -
    # Phi((tau_{k-1} - yhat)/sd), with tau_0 = -Inf and tau_{ncat} = +Inf.
    # Computed for all ordered items at once (two vectorised pnorm calls
    # instead of a per-item loop with two scalar pnorm calls each).
    ord_idx <- which(var_type == "ordered")
    if (length(ord_idx) > 0L) {
      # missing responses -> NA density
      na_mask <- is.na(y_i[ord_idx])
      if (any(na_mask)) {
        fy[ord_idx[na_mask]] <- as.numeric(NA)
      }
      obs <- ord_idx[!na_mask]
      if (length(obs) > 0L) {
        k <- y_i[obs] # observed category (1-based)
        # number of categories per variable, and the 0-based offset of each
        # variable's thresholds within 'th' (match() locates the first
        # threshold, which is robust to placeholder entries for the
        # numeric variables, i.e. th.idx == 0)
        ncat <- tabulate(th_idx, nbins = nvar_g)[obs] + 1L
        th_off <- match(obs, th_idx) - 1L

        # upper threshold tau_k (Inf for the top category)
        upper <- th[th_off + k]
        upper[k == ncat] <- Inf
        # lower threshold tau_{k-1} (-Inf for the first category)
        lower <- rep(-Inf, length(obs))
        nz <- k > 1L
        lower[nz] <- th[th_off[nz] + k[nz] - 1L]

        sd_obs <- theta_sd[obs]
        yhat_obs <- yhat[obs]
        prob <- pnorm((upper - yhat_obs) / sd_obs) -
          pnorm((lower - yhat_obs) / sd_obs)
        prob[prob < .Machine$double.eps] <- .Machine$double.eps
        fy[obs] <- if (log) log(prob) else prob
      }
    }

    # guard against unexpected variable types
    bad_idx <- which(!var_type %in% c("numeric", "ordered"))
    if (length(bad_idx) > 0L) {
      v <- bad_idx[1L]
      lav_msg_stop(gettextf("unknown type: `%1$s' for variable: %2$s",
        lavdata@ov$type[v], lavdata@ov$name[v])
      )
    }
  }

  fy
}

# gather the (unconditional) per-block ingredients shared by both
# correlation-preserving transformations: model-implied Sigma, latent
# covariance Psi (veta) and factor loadings Lambda
lav_predict_tmat_blocks <- function(lavobject = NULL,
                                    lavmodel = NULL, lavimplied = NULL) {
  if (!is.null(lavobject)) {
    lavmodel <- lavobject@Model
    lavimplied <- lavobject@implied
  }

  if (is.null(lavimplied) || length(lavimplied) == 0L) {
    lavimplied <- lav_model_implied(lavmodel)
  }
  if (lavmodel@conditional.x) {
    lavimplied <- lav_model_implied_cond2uncond(lavimplied)
  }

  list(
    sigma = lavimplied$cov,
    veta = lav_model_veta(lavmodel = lavmodel, remove_dummy_lv = FALSE),
    lambda = lav_model_lambda(lavmodel, remove_dummy_lv = FALSE,
                              use_wmat = TRUE),
    nblocks = lavmodel@nblocks
  )
}

# compute `transformation' matrix to convert regression factor scores
# to (Green's) correlation-preserving factor scores
lav_predict_tmat_green <- function(lavobject = NULL,
                                   lavmodel = NULL, lavimplied = NULL) {
  mm <- lav_predict_tmat_blocks(lavobject, lavmodel, lavimplied)

  lapply(seq_len(mm$nblocks), function(b) {
    lav_predict_tmat_green_internal(mm$sigma[[b]], mm$veta[[b]], mm$lambda[[b]])
  })
}

# compute `transformation' matrix to convert Bartlett factor scores
# to (Krijnen/McDonald) correlation-preserving factor scores
lav_predict_tmat_det <- function(lavobject = NULL,
                                 lavmodel = NULL, lavimplied = NULL) {
  mm <- lav_predict_tmat_blocks(lavobject, lavmodel, lavimplied)

  lapply(seq_len(mm$nblocks), function(b) {
    lav_predict_tmat_det_internal(mm$sigma[[b]], mm$veta[[b]], mm$lambda[[b]])
  })
}

# single block only (Green/regression)
lav_predict_tmat_green_internal <- function(sigma_1 = NULL, veta = NULL,
                                            lambda = NULL) {
  sigma_inv <- solve(sigma_1)
  veta_sqrt <- lav_mat_sym_sqrt(veta)
  veta32 <- veta %*% veta_sqrt
  tmp <- veta32 %*% t(lambda) %*% sigma_inv %*% lambda %*% veta32
  tmp_inv_sqrt <- lav_mat_sym_sqrt(solve(tmp))
  veta_sqrt %*% tmp_inv_sqrt %*% veta_sqrt
}

# single block only (Bartlett/determinacy), for internal use in
# lav_sam_step1_local()
lav_predict_tmat_det_internal <- function(sigma_1 = NULL, veta = NULL,
                                          lambda = NULL) {
    sigma_inv <- solve(sigma_1)
    veta_sqrt <- lav_mat_sym_sqrt(veta)
    veta_inv_sqrt <- lav_mat_sym_sqrt(solve(veta))
    tmp <- veta_sqrt %*% t(lambda) %*% sigma_inv %*% lambda %*% veta_sqrt
    tmp_sqrt <- lav_mat_sym_sqrt(tmp)
    tmat <- veta_sqrt %*% tmp_sqrt %*% veta_inv_sqrt
    tmat
}
