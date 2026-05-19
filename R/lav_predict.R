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
      fsm = FALSE, rel = FALSE, optim.method = "bfgs"
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
lavPredict <- function(object, newdata = NULL, # keep order of predict(), 0.6-7 # nolint start
                       type = "lv", method = "EBM", transform = FALSE,
                       se = "none", acov = "none", label = TRUE, fsm = FALSE,
                       mdist = FALSE, rel = FALSE,
                       append.data = FALSE, assemble = FALSE, # or TRUE?
                       level = 1L, optim.method = "bfgs", ETA = NULL,
                       drop.list.single.group = TRUE) {                         # nolint end
  # check object
  object <- lav_object_check_version(object)

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
    mdist = mdist, append_data = append.data, assemble = assemble,
    level = level, optim_method = optim.method, eta = ETA,
    drop_list_single_group = drop.list.single.group
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
                               drop_list_single_group = TRUE) {
  # type
  type <- tolower(type)
  lavpta <- lav_partable_attributes(lavpartable)
  if (type %in% c("latent", "lv", "factor", "factor.score", "factorscore")) {
    type <- "lv"
  } else if (type %in% c("ov", "yhat")) {
    type <- "yhat"
  } else if (type %in% c("residuals", "resid", "error")) {
    type <- "resid"
  }

  # if resid, not for categorical
  if (type == "resid" && lavmodel@categorical) {
    lav_msg_stop(gettext(
      "casewise residuals not available if data is categorical"))
  }

  # append.data? check level
  if (append_data && level > 1L) {
    lav_msg_warn(gettext("append.data not available if level > 1L"))
    append_data <- FALSE
  }

  # mdist? -> fsm = TRUE
  if (mdist) {
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
    if (lavmodel@categorical) {
      se <- acov <- "none"
      lav_msg_warn(gettext(
        "standard errors not available (yet) for non-normal data"))
    }
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
    new_data <- lav_lavdata(
      data = newdata,
      group = lavdata@group,
      ov_names = lavdata@ov.names,
      ov_names_x = lavdata@ov.names.x,
      ordered = ov$name[ov$type == "ordered"],
      lavoptions = list(
        std.ov = lavdata@std.ov,
        group.label = lavdata@group.label,
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
  }

  if (type == "lv") {
    if (!is.null(eta)) {
      lav_msg_warn(gettext("lvs will be predicted here;
                           supplying ETA has no effect"))
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
      lavimplied = lavimplied, se = se, acov = acov, level = level,
      data_obs = data_obs, exo = exo, method = method,
      fsm = fsm, rel = rel, transform = transform, optim_method = optim_method
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
      if (lavdata@nlevels == 1L) {
        b <- g
      } else {
        b <- (g - 1) * lavdata@nlevels + level
      }
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
        if (lavdata@nlevels == 1L) {
          b <- g
        } else {
          b <- (g - 1) * lavdata@nlevels + level
        }
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
#         #fs.inv.sqrt <- lav_matrix_symmetric_sqrt(FS.cov.inv)
#         #veta.sqrt <- lav_matrix_symmetric_sqrt(VETA[[g]])
#         #tmp <- FS.centered %*% fs.inv.sqrt %*% veta.sqrt
#         tmp <- FS.centered %*% t(tmat[[b]])
#         ret <- t(t(tmp) + drop(EETA[[g]]))
#
#         ret
#       })
#     }

    # new in 0.6-17
    if (mdist) {
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
        if (lavdata@nlevels == 1L) {
          b <- g
        } else {
          b <- (g - 1) * lavdata@nlevels + level
        }
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
          if (lavdata@nlevels == 1L) {
            b <- g
          } else {
            b <- (g - 1) * lavdata@nlevels + level
          }
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
        if (lavdata@nlevels > 1L) {
          gg <- (g - 1) * lavdata@nlevels + level
        } else {
          gg <- g
        }

        if (append_data) {
          colnames(out[[g]]) <- c(
            lavpta$vnames$lv[[gg]],
            ov_names[[g]]
          ) # !not gg
        } else {
          colnames(out[[g]]) <- lavpta$vnames$lv[[gg]]
        }

        if (fsm) {
          if (is.null(fsm_1[[g]])) {
            # skip
          } else if (is.matrix(fsm_1[[g]])) {
            dimnames(fsm_1[[g]]) <- list(
              lavpta$vnames$lv[[gg]],
              # ov.names[[g]]) # !not gg
              lavpta$vnames$ov.ind[[gg]]
            )
          } else if (is.list(fsm_1[[g]])) {
            fsm_1[[g]] <- lapply(fsm_1[[g]], function(x) {
              dimnames(x) <- list(
                lavpta$vnames$lv[[gg]],
                # ov.names[[g]]) # !not gg
                lavpta$vnames$ov.ind[[gg]]
              )
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

    # label?
    if (label) {
      for (g in seq_len(lavdata@ngroups)) {
        colnames(out[[g]]) <- lavpta$vnames$ov[[g]]
      }
    }

    # mdist
    if (mdist) {
      mm_lambda <- lav_model_lambda(
        lavmodel = lavmodel,
        remove_dummy_lv = FALSE
      )
      mdist_1 <- lapply(seq_len(lavdata@ngroups), function(g) {
        sigma <- lavimplied$cov[[g]]
        la <- mm_lambda[[g]]
        if (type == "resid") {
          ila <- diag(ncol(sigma)) - la %*% fsm[[g]]
          omega_e <- ila %*% sigma %*% t(ila)
          eig <- eigen(omega_e, symmetric = TRUE)
          a <- eig$vectors[, seq_len(nrow(la) - ncol(la)),
            drop = FALSE
          ]
        } else if (type == "yhat") {
          laa <- la %*% fsm[[g]]
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

# internal function
lav_predict_eta <- function(lavobject = NULL, # for convenience
                            # sub objects
                            lavmodel = NULL, lavdata = NULL,
                            lavsamplestats = NULL,
                            lavimplied = NULL,
                            # new data
                            data_obs = NULL, exo = NULL,
                            # options
                            method = "EBM",
                            fsm = FALSE,
              rel = FALSE,
                            transform = FALSE,
                            se = "none", acov = "none",
                            level = 1L,
                            optim_method = "bfgs") {
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

  # normal case?
  if (all(lavdata@ov$type == "numeric")) {
    if (method == "ebm") {
      out <- lav_predict_eta_normal(
        lavobject = lavobject,
        lavmodel = lavmodel, lavdata = lavdata,
        lavimplied = lavimplied, se = se, acov = acov,
        level = level, lavsamplestats = lavsamplestats,
        data_obs = data_obs, exo = exo, fsm = fsm, rel = rel,
        transform = transform
      )
    } else if (method == "ml") {
      out <- lav_predict_eta_bartlett(
        lavobject = lavobject,
        lavmodel = lavmodel, lavdata = lavdata,
        lavimplied = lavimplied, se = se, acov = acov,
        level = level, lavsamplestats = lavsamplestats,
        data_obs = data_obs, exo = exo, fsm = fsm, rel = rel,
        transform = transform
      )
    } else {
      lav_msg_stop(gettextf("unknown method: %s.", method))
    }
  } else {
    if (method == "ebm") {
      out <- lav_predict_eta_ebm_ml(
        lavobject = lavobject,
        lavmodel = lavmodel, lavdata = lavdata,
        lavsamplestats = lavsamplestats, se = se, acov = acov,
        level = level, data_obs = data_obs, exo = exo,
        ml = FALSE, optim_method = optim_method
      )
    } else if (method == "ml") {
      out <- lav_predict_eta_ebm_ml(
        lavobject = lavobject,
        lavmodel = lavmodel, lavdata = lavdata,
        lavsamplestats = lavsamplestats, se = se, acov = acov,
        level = level, data_obs = data_obs, exo = exo,
        ml = TRUE, optim_method = optim_method
      )
    } else {
      lav_msg_stop(gettextf("unknown method: %s.", method))
    }
  }

  out
}


# factor scores - normal case
# NOTE: this is the classic 'regression' method; for the linear/continuous
#       case, this is equivalent to both EB and EBM
lav_predict_eta_normal <- function(lavobject = NULL, # for convenience
                                   # sub objects
                                   lavmodel = NULL, lavdata = NULL,
                                   lavsamplestats = NULL,
                                   lavimplied = NULL,
                                   # optional new data
                                   data_obs = NULL, exo = NULL,
                                   se = "none", acov = "none", level = 1L,
                                   fsm = FALSE, rel = FALSE,
                                   transform = FALSE) {
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
        mp_1[[g]] <- lav_data_missing_patterns(data_obs[[g]])
      }
    } else {
      mp_1 <- lavdata@Mp
    }
  }

  mm_lambda <- lav_model_lambda(lavmodel = lavmodel, remove_dummy_lv = FALSE)
  sigma_hat <- lavimplied$cov
  sigma_inv <- lapply(sigma_hat, MASS::ginv)
  veta <- lav_model_veta(lavmodel = lavmodel)
  eeta <- lav_model_eeta(lavmodel = lavmodel, lavsamplestats = lavsamplestats)
  ey <- lav_model_ey(lavmodel = lavmodel, lavsamplestats = lavsamplestats)

  fs <- vector("list", length = lavdata@ngroups)
  if (fsm) {
    fsm_1 <- vector("list", length = lavdata@ngroups)
  }
  if (rel) {
    rel_1 <- vector("list", length = lavdata@ngroups)
  }
  if (transform) {
    tmat <- lav_predict_tmat_green(lavmodel = lavmodel,
                                   lavimplied = lavimplied)
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
    if (lavdata@nlevels == 1L) {
      b <- g
    } else {
      b <- (g - 1) * lavdata@nlevels + level
    }

    veta_g <- veta[[b]]
    eeta_g <- eeta[[b]]
    lambda_g <- mm_lambda[[b]]
    ey_g <- ey[[b]]
    sigma_inv_g <- sigma_inv[[b]]

    if (lavdata@nlevels > 1L) {
      lp <- lavdata@Lp[[g]]
      ylp <- lavsamplestats@YLp[[g]]

      # implied for this group
      group_idx <- (g - 1) * lavdata@nlevels + seq_len(lavdata@nlevels)
      implied_group <- lapply(lavimplied, function(x) x[group_idx])

      # random effects (=random intercepts or cluster means)
      out <- lav_mvnorm_cluster_implied22l(
        lp = lp,
        implied = implied_group
      )
      mb_j <- lav_mvnorm_cluster_em_estep_ranef(
        ylp = ylp, lp = lp,
        sigma_w = out$sigma.w, sigma_b = out$sigma.b,
        sigma_zz = out$sigma.zz, sigma_yz = out$sigma.yz,
        mu_z = out$mu.z, mu_w = out$mu.w, mu_b = out$mu.b,
        se = FALSE
      )

      ov_idx <- lp$ov.idx

      if (level == 1L) {
        data_w <- data_obs[[g]][, ov_idx[[1]]]
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

    # sampling weights? -- CHECKME: needed??
    if (!is.null(lavdata@weights[[g]]) && level == 1L) {
      # EY.g is already weighted
      # use sampling.weights.normalization == "group"
      wt <- lavdata@weights[[g]]
      wt2 <- wt / sum(wt) * lavdata@nobs[[g]]
      yc <- yc * sqrt(wt2)
    }

    # global factor score coefficient matrix 'C'
    fsc <- veta_g %*% t(lambda_g) %*% sigma_inv_g

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
    if (lavdata@missing %in% c("ml", "ml.x")) {
      # missing patterns for this group
      mp <- mp_1[[g]]

      # factor scores container
      fs_g <- matrix(as.numeric(NA), nrow(yc), ncol = length(eeta_g))

      # if(fsm) {
      #    FSM.g <- vector("list", length = Mp$npatterns)
      # }

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
        sigma_22_inv <- try(lav_matrix_symmetric_inverse_update(
          s_inv = sigma_inv_g, rm_idx = na_idx, logdet = FALSE
        ), silent = TRUE)
        if (inherits(sigma_22_inv, "try-error")) {
          lav_msg_stop(gettext("Sigma_22.inv cannot be inverted"))
        }

        lambda <- lambda_g[var_idx, , drop = FALSE]
        fsc <- veta_g %*% t(lambda) %*% sigma_22_inv

        # FSM?
        # if(fsm) {
        #    tmp <- matrix(as.numeric(NA), nrow = ncol(lambda),
        #                  ncol = ncol(Yc))
        #    tmp[,var.idx] <- FSC
        #    FSM.g[[p]] <- tmp
        # }

        # factor score for this pattern
        fs_g[mp$case.idx[[p]], ] <- t(fsc %*% t(oc) + eeta_g)

        # SE?
        if (se == "standard") {
          tmp <- (veta_g - veta_g %*% t(lambda) %*%
            sigma_22_inv %*% lambda %*% veta_g)
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
      if (lavdata@missing %in% c("ml", "ml.x")) {
        se_1[[g]] <- se_g
        if (acov == "standard") {
          acov_1[[g]] <- acov_g
        }
      } else { # complete data
        tmp <- (veta_g - veta_g %*% t(lambda_g) %*%
          sigma_inv_g %*% lambda_g %*% veta_g)
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

# factor scores - normal case - Bartlett method
# NOTES: 1) this is the classic 'Bartlett' method; for the linear/continuous
#           case, this is equivalent to 'ML'
#        2) the usual formula is:
#               FSC = solve(lambda' theta.inv lambda) (lambda' theta.inv)
#           BUT to deal with singular THETA (with zeroes on the diagonal),
#           we use the 'GLS' version instead:
#               FSC = solve(lambda' sigma.inv lambda) (lambda' sigma.inv)
#           Reference: Bentler & Yuan (1997) 'Optimal Conditionally Unbiased
#                      Equivariant Factor Score Estimators'
#                      in Berkane (Ed) 'Latent variable modeling with
#                      applications to causality' (Springer-Verlag)
#        3) instead of solve(), we use MASS::ginv, for special settings where
#           -by construction- (lambda' sigma.inv lambda) is singular
#           note: this will destroy the conditionally unbiased property
#                 of Bartlett scores!!
lav_predict_eta_bartlett <- function(lavobject = NULL, # for convenience
                                     # sub objects
                                     lavmodel = NULL, lavdata = NULL,
                                     lavsamplestats = NULL,
                                     lavimplied = NULL,
                                     # optional new data
                                     data_obs = NULL, exo = NULL,
                                     se = "none", acov = "none", level = 1L,
                                     fsm = FALSE, rel = FALSE,
                                     transform = FALSE) {
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
        mp_1[[g]] <- lav_data_missing_patterns(data_obs[[g]])
      }
    } else {
      mp_1 <- lavdata@Mp
    }
  }

  mm_lambda <- lav_model_lambda(lavmodel = lavmodel, remove_dummy_lv = FALSE)
  sigma_hat <- lavimplied$cov
  sigma_inv <- lapply(lavimplied$cov, MASS::ginv)
  veta <- lav_model_veta(lavmodel = lavmodel) # for se only
  eeta <- lav_model_eeta(lavmodel = lavmodel, lavsamplestats = lavsamplestats)
  ey <- lav_model_ey(lavmodel = lavmodel, lavsamplestats = lavsamplestats)

  fs <- vector("list", length = lavdata@ngroups)
  if (fsm) {
    fsm_1 <- vector("list", length = lavdata@ngroups)
  }
  if (rel) {
    rel_1 <- vector("list", length = lavdata@ngroups)
  }
  if (transform) {
    tmat <- lav_predict_tmat_det(lavmodel = lavmodel,
                                 lavimplied = lavimplied)
  }

  if (acov != "none") se <- acov # acov implies SE
  if (se != "none") {
    se_1 <- vector("list", length = lavdata@ngroups)
    # return full sampling covariance matrix?
    if (acov != "none") {
      acov_1 <- vector("list", length = lavdata@ngroups)
    }
  }

  for (g in 1:lavdata@ngroups) {
    # determine block
    if (lavdata@nlevels == 1L) {
      b <- g
    } else {
      b <- (g - 1) * lavdata@nlevels + level
    }

    veta_g <- veta[[b]]
    eeta_g <- eeta[[b]]
    lambda_g <- mm_lambda[[b]]
    ey_g <- ey[[b]]
    sigma_inv_g <- sigma_inv[[b]]

    if (lavdata@nlevels > 1L) {
      lp <- lavdata@Lp[[g]]
      ylp <- lavsamplestats@YLp[[g]]

      # implied for this group
      group_idx <- (g - 1) * lavdata@nlevels + seq_len(lavdata@nlevels)
      implied_group <- lapply(lavimplied, function(x) x[group_idx])

      # random effects (=random intercepts or cluster means)
      # NOTE: is the 'ML' way not simply using the observed cluster
      #       means?
      out <- lav_mvnorm_cluster_implied22l(
        lp = lp,
        implied = implied_group
      )
      mb_j <- lav_mvnorm_cluster_em_estep_ranef(
        ylp = ylp, lp = lp,
        sigma_w = out$sigma.w, sigma_b = out$sigma.b,
        sigma_zz = out$sigma.zz, sigma_yz = out$sigma.yz,
        mu_z = out$mu.z, mu_w = out$mu.w, mu_b = out$mu.b,
        se = FALSE
      )

      ov_idx <- lp$ov.idx

      if (level == 1L) {
        data_w <- data_obs[[g]][, ov_idx[[1]]]
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

    # sampling weights? CHECKME: needed??
    if (!is.null(lavdata@weights[[g]]) && level == 1L) {
      # EY.g is already weighted
      # use sampling.weights.normalization == "group"
      wt <- lavdata@weights[[g]]
      wt2 <- wt / sum(wt) * lavdata@nobs[[g]]
      yc <- yc * sqrt(wt2)
    }

    # global factor score coefficient matrix 'C'
    fsc <- (MASS::ginv(t(lambda_g) %*% sigma_inv_g %*% lambda_g)
    %*% t(lambda_g) %*% sigma_inv_g)

    # transform?
    if (transform) {
      fsc <- tmat[[b]] %*% fsc
    }

    # store fsm?
    if (fsm) {
      # store fsm?
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
    if (lavdata@missing %in% c("ml", "ml.x")) {
      # missing patterns for this group
      mp <- mp_1[[g]]

      # factor scores container
      fs_g <- matrix(as.numeric(NA), nrow(yc), ncol = length(eeta_g))

      # if(fsm) {
      #    FSM.g <- vector("list", length = Mp$npatterns)
      # }

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
        sigma_22_inv <- try(lav_matrix_symmetric_inverse_update(
          s_inv = sigma_inv_g, rm_idx = na_idx, logdet = FALSE
        ), silent = TRUE)
        if (inherits(sigma_22_inv, "try-error")) {
          lav_msg_stop(gettext("Sigma_22.inv cannot be inverted"))
        }

        lambda <- lambda_g[var_idx, , drop = FALSE]
        fsc <- (MASS::ginv(t(lambda) %*% sigma_22_inv %*% lambda)
        %*% t(lambda) %*% sigma_22_inv)

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
        #
        zero_idx <- which(apply(fsc, 1L, function(x) all(x == 0)))
        if (length(zero_idx) > 0L) {
          fsc[zero_idx, ] <- NA
        }

        # FSM?
        # if(fsm) {
        #    tmp <- matrix(as.numeric(NA), nrow = ncol(lambda),
        #                  ncol = ncol(Yc))
        #    tmp[,var.idx] <- FSC
        #    FSM.g[[p]] <- tmp
        # }

        # factor scores for this pattern
        fs_g[mp$case.idx[[p]], ] <- t(fsc %*% t(oc) + eeta_g)

        # SE?
        if (se == "standard") {
          tmp <- (MASS::ginv(t(lambda) %*% sigma_22_inv %*% lambda)
          - veta_g)
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
      }

      # what about FSM? There is no single one, but as many as patterns
      # if(fsm) {
      #    # use 'global' version (just like in complete case)
      #    FSM[[g]] <- ( MASS::ginv(t(LAMBDA.g) %*% Sigma.inv.g %*%
      #                    LAMBDA.g) %*% t(LAMBDA.g) %*% Sigma.inv.g )
      # }
    } else {
      # compute factor scores
      fs_g <- t(fsc %*% t(yc) + eeta_g)
    }

    # replace values in dummy lv's by their observed counterpart
    if (length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L && level == 1L) {
      fs_g[, lavmodel@ov.y.dummy.lv.idx[[g]]] <-
        data_obs[[g]][, lavmodel@ov.y.dummy.ov.idx[[g]], drop = FALSE]
    }
    if (length(lavmodel@ov.x.dummy.lv.idx[[g]]) > 0L && level == 1L) {
      fs_g[, lavmodel@ov.x.dummy.lv.idx[[g]]] <-
        data_obs[[g]][, lavmodel@ov.x.dummy.ov.idx[[g]], drop = FALSE]
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
    if (se == "standard" && !transform) {
      if (lavdata@missing %in% c("ml", "ml.x")) {
        se_1[[g]] <- se_g
        if (acov == "standard") {
          acov_1[[g]] <- acov_g
        }
      } else { # complete data

        # the traditional formula is:
        #     solve(t(lambda) %*% solve(theta) %*% lambda)
        # but we replace it by
        #     solve( t(lambda) %*% solve(sigma) %*% lambda ) - psi
        # to handle negative variances
        # in addition, we use ginv
        tmp <- (MASS::ginv(t(lambda_g) %*% sigma_inv_g %*% lambda_g)
        - veta_g)
        tmp_d <- diag(tmp)
        tmp_d[tmp_d < 1e-05] <- as.numeric(NA)
        se_1[[g]] <- matrix(sqrt(tmp_d), nrow = 1L)

        # return full sampling covariance matrix?
        if (acov == "standard") {
          acov_1[[g]] <- tmp
        }
      }
    } # se
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

# factor scores - EBM or ML
lav_predict_eta_ebm_ml <- function(lavobject = NULL, # for convenience
                                   # sub objects
                                   lavmodel = NULL, lavdata = NULL,
                                   lavsamplestats = NULL,
                                   # optional new data
                                   data_obs = NULL, exo = NULL,
                                   se = "none", acov = "none", level = 1L,
                                   ml = FALSE,
                                   optim_method = "bfgs") {
  optim_method <- tolower(optim_method)

  stopifnot(optim_method %in% c("nlminb", "bfgs"))

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
  # if(se != "none") {
  #    warning("lavaan WARNING: standard errors are
  #               not available (yet) for the non-normal case")
  # }

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

  # local objective function: x = lv values
  f_eta_i <- function(x, y_i, x_i, mu_i) {
    # add 'dummy' values (if any) for ov.y
    if (length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L) {
      x2 <- c(x - mu_i, data_obs[[g]][i,
        lavmodel@ov.y.dummy.ov.idx[[g]],
        drop = FALSE
      ])
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
      log = TRUE
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

    # casewise for now
    n <- nrow(data_obs[[g]])
    for (i in 1:n) {
      # eXo?
      if (!is.null(exo[[g]])) {
        x_i <- exo[[g]][i, , drop = FALSE]
      } else {
        x_i <- NULL
      }
      mu_i <- eetax[[g]][i, , drop = FALSE]
      y_i <- data_obs[[g]][i, , drop = FALSE]

      ### DEBUG ONLY:
      # cat("i = ", i, "mu.i = ", mu.i, "\n")

      start_1 <- numeric(nfac) # initial values for eta

      if (!all(is.na(y_i))) {
        # find best values for eta.i
        if (optim_method == "nlminb") {
          out <- nlminb(
            start = start_1, objective = f_eta_i,
            gradient = NULL, # for now
            control = list(rel.tol = 1e-8),
            y_i = y_i, x_i = x_i, mu_i = mu_i
          )
        } else if (optim_method == "bfgs") {
          out <- optim(
            par = start_1, fn = f_eta_i,
            gr = NULL,
            control = list(reltol = 1e-8, fnscale = 1.1),
            method = "BFGS",
            y_i = y_i, x_i = x_i, mu_i = mu_i
          )
        }
        if (out$convergence == 0L) {
          eta_i <- out$par
        } else {
          eta_i <- rep(as.numeric(NA), nfac)
        }
      } else {
        eta_i <- rep(as.numeric(NA), nfac)
      }

      # add dummy ov.y lv values
      if (length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L) {
        eta_i <- c(eta_i, data_obs[[g]][i,
          lavmodel@ov.y.dummy.ov.idx[[g]],
          drop = FALSE
        ])
      }

      fs[[g]][i, ] <- eta_i
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
# Two types: 1) nrow(ETA) = nrow(X) (factor scores)
#            2) nrow(ETA) = 1L (given values)
#
# in both cases, we return [nobs x nvar] matrix per group
lav_predict_yhat <- function(lavobject = NULL, # for convience
                             # sub objects
                             lavmodel = NULL, lavdata = NULL,
                             lavsamplestats = NULL,
                             lavimplied = NULL,
                             # new data
                             data_obs = NULL, exo = NULL,
                             # ETA values
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

  # do we get values for ETA? If not, use `predict' to get plausible values
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
        lav_msg_stop(gettext("nrow(ETA) != lavsamplestats@ntotal"))
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
                           # ETA values
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
                                 th = NULL, th_idx = NULL, log = TRUE) {
  mm_in_group <- 1:lavmodel@nmat[g] + cumsum(c(0, lavmodel@nmat))[g]
  mlist <- lavmodel@GLIST[mm_in_group]

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
    fy <- numeric(lavmodel@nvar[g])
    for (v in seq_len(lavmodel@nvar[g])) {
      if (lavdata@ov$type[v] == "numeric") {
        ### FIXME!!! we can do all numeric vars at once!!
        fy[v] <- dnorm(y_i[v],
          mean = yhat[v], sd = theta_sd[v],
          log = log
        )
      } else if (lavdata@ov$type[v] == "ordered") {
        # handle missing value
        if (is.na(y_i[v])) {
          fy[v] <- as.numeric(NA)
        } else {
          th_y <- th[th_idx == v]
          th_y_1 <- c(-Inf, th_y, Inf)
          k <- y_i[v]
          p1 <- pnorm((th_y_1[k + 1] - yhat[v]) / theta_sd[v])
          p2 <- pnorm((th_y_1[k] - yhat[v]) / theta_sd[v])
          prob <- (p1 - p2)
          if (prob < .Machine$double.eps) {
            prob <- .Machine$double.eps
          }
          if (log) {
            fy[v] <- log(prob)
          } else {
            fy[v] <- prob
          }
        }
      } else {
        lav_msg_stop(gettextf("unknown type: `%1$s' for variable: %2$s",
          lavdata@ov$type[v], lavdata@ov$name[v])
        )
      }
    }
  }

  fy
}

# compute `transformation' matrix to convert regression factor scores
# to (Green's) correlation-preserving factor scores
lav_predict_tmat_green <- function(lavobject = NULL,
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
  sigma_1 <- lavimplied$cov
  veta_1 <- lav_model_veta(lavmodel = lavmodel, remove_dummy_lv = FALSE)
  mm_lambda <- lav_model_lambda(lavmodel, remove_dummy_lv = FALSE)

  nblocks <- lavmodel@nblocks
  tmat <- vector("list", length = nblocks)

  # compute tmat per block
  for (b in seq_len(nblocks)) {
    sigma_b <- sigma_1[[b]]
    sigma_b_inv <- solve(sigma_b)
    veta <- veta_1[[b]]
    veta_sqrt <- lav_matrix_symmetric_sqrt(veta)
    veta32 <- veta %*% veta_sqrt
    lambda <- mm_lambda[[b]]
    tmp <- veta32 %*% t(lambda) %*% sigma_b_inv %*% lambda %*% veta32
    tmp_inv_sqrt <- lav_matrix_symmetric_sqrt(solve(tmp))
    tmat[[b]] <- veta_sqrt %*% tmp_inv_sqrt %*% veta_sqrt
  }

  tmat
}

# compute `transformation' matrix to convert Bartlett factor scores
# to (Krijnen/McDonald) correlation-preserving factor scores
lav_predict_tmat_det <- function(lavobject = NULL,
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
  sigma_1 <- lavimplied$cov
  veta_1 <- lav_model_veta(lavmodel = lavmodel, remove_dummy_lv = FALSE)
  mm_lambda <- lav_model_lambda(lavmodel, remove_dummy_lv = FALSE)

  nblocks <- lavmodel@nblocks
  tmat <- vector("list", length = nblocks)

  # compute tmat per block
  for (b in seq_len(nblocks)) {
    sigma_b <- sigma_1[[b]]
    sigma_b_inv <- solve(sigma_b)
    veta <- veta_1[[b]]
    veta_sqrt <- lav_matrix_symmetric_sqrt(veta)
    veta_inv_sqrt <- lav_matrix_symmetric_sqrt(solve(veta))
    lambda <- mm_lambda[[b]]
    tmp <- veta_sqrt %*% t(lambda) %*% sigma_b_inv %*% lambda %*% veta_sqrt
    tmp_sqrt <- lav_matrix_symmetric_sqrt(tmp)
    tmat[[b]] <- veta_sqrt %*% tmp_sqrt %*% veta_inv_sqrt
  }

  tmat
}

# single block only, for internal use in lav_sam_step1_local()
lav_predict_tmat_det_internal <- function(sigma_1 = NULL, veta = NULL,
                                          lambda = NULL) {
    sigma_inv <- solve(sigma_1)
    veta_sqrt <- lav_matrix_symmetric_sqrt(veta)
    veta_inv_sqrt <- lav_matrix_symmetric_sqrt(solve(veta))
    tmp <- veta_sqrt %*% t(lambda) %*% sigma_inv %*% lambda %*% veta_sqrt
    tmp_sqrt <- lav_matrix_symmetric_sqrt(tmp)
    tmat <- veta_sqrt %*% tmp_sqrt %*% veta_inv_sqrt
    tmat
}
