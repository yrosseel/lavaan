# lavPredict() contains a collection of `predict' methods
# the unifying theme is that they all rely on the (unknown, to be estimated)
# or (known, apriori specified) values for the latent variables
#
# lv: lavtent variables (aka `factor scores')
# ov: predict linear part of y_i
#
# - YR 11 June 2013: first version, in order to get factor scores for the
#                    categorical case
# - YR 12 Jan 2014: refactoring + lav_predict_fy (to be used by estimator MML)
#

# overload standard R function `predict'
setMethod(
  "predict", "lavaan",
  function(object, newdata = NULL) {
    lavPredict(
      object = object, newdata = newdata, type = "lv", method = "EBM",
      fsm = FALSE, optim.method = "bfgs"
    )
  }
)

# efaList version
predict.efaList <- function(object, ...) {
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
lavPredict <- function(object, newdata = NULL, # keep order of predict(), 0.6-7
                       type = "lv", method = "EBM", transform = FALSE,
                       se = "none", acov = "none", label = TRUE, fsm = FALSE,
                       mdist = FALSE,
                       append.data = FALSE, assemble = FALSE, # or TRUE?
                       level = 1L, optim.method = "bfgs", ETA = NULL,
                       drop.list.single.group = TRUE) {
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
  # backward compatibility
  if (.hasSlot(object, "h1")) {
    lavh1 <- object@h1
  } else {
    lavh1 <- lav_h1_implied_logl(
      lavdata = object@Data,
      lavsamplestats = object@SampleStats,
      lavoptions = object@Options
    )
  }
  lavimplied <- object@implied

  res <- lav_predict_internal(
    lavmodel = lavmodel, lavdata = lavdata,
    lavsamplestats = lavsamplestats, lavimplied = lavimplied, lavh1 = lavh1,
    lavpartable = object@ParTable, newdata = newdata, type = type, method = method,
    transform = transform, se = se, acov = acov, label = label, fsm = fsm,
    mdist = mdist, append.data = append.data, assemble = assemble,
    level = level, optim.method = optim.method, ETA = ETA,
    drop.list.single.group = drop.list.single.group
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
                                 se = "none", acov = "none", label = TRUE, fsm = FALSE,
                                 mdist = FALSE,
                                 append.data = FALSE, assemble = FALSE, # or TRUE?
                                 level = 1L, optim.method = "bfgs", ETA = NULL,
                                 drop.list.single.group = TRUE) {
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
    lav_msg_stop(gettext("casewise residuals not available if data is categorical"))
  }

  # append.data? check level
  if (append.data && level > 1L) {
    lav_msg_warn(gettext("append.data not available if level > 1L"))
    append.data <- FALSE
  }

  # mdist? -> fsm = TRUE
  if (mdist) {
    fsm <- TRUE
  }

  # se?
  if (acov != "none") {
    se <- acov # ACOV implies SE
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
    #    warning("lavaan WARNING: standard errors not available (yet) for missing data + fiml")
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
      data.obs <- lavdata@X
      ov.names <- lavdata@ov.names
    }
    eXo <- lavdata@eXo
  } else {
    OV <- lavdata@ov
    newData <- lavData(
      data = newdata,
      group = lavdata@group,
      ov.names = lavdata@ov.names,
      ov.names.x = lavdata@ov.names.x,
      ordered = OV$name[OV$type == "ordered"],
      lavoptions = list(
        std.ov = lavdata@std.ov,
        group.label = lavdata@group.label,
        missing = lavdata@missing,
        warn = TRUE
      ), # was FALSE before?
      allow.single.case = TRUE
    )
    # if ordered, check if number of levels is till the same (new in 0.6-7)
    if (lavmodel@categorical) {
      orig.ordered.idx <- which(lavdata@ov$type == "ordered")
      orig.ordered.lev <- lavdata@ov$nlev[orig.ordered.idx]
      match.new.idx <- match(
        lavdata@ov$name[orig.ordered.idx],
        newData@ov$name
      )
      new.ordered.lev <- newData@ov$nlev[match.new.idx]
      if (any(orig.ordered.lev - new.ordered.lev != 0)) {
        lav_msg_stop(
          gettext("mismatch number of categories for some ordered variables
                  in newdata compared to original data.")
        )
      }
    }
    data.obs <- newData@X
    eXo <- newData@eXo
    ov.names <- newData@ov.names
  }

  if (type == "lv") {
    if (!is.null(ETA)) {
      lav_msg_warn(gettext("lvs will be predicted here;
                           supplying ETA has no effect"))
    }

    # post fit check (lv pd?)
    # ok <- lav_object_post_check(object)
    # if(!ok) {
    #    stop("lavaan ERROR: lavInspect(,\"post.check\") is not TRUE; factor scores can not be computed. See the WARNING message.")
    # }

    out <- lav_predict_eta(
      lavobject = NULL, lavmodel = lavmodel,
      lavdata = lavdata, lavsamplestats = lavsamplestats,
      lavimplied = lavimplied, se = se, acov = acov, level = level,
      data.obs = data.obs, eXo = eXo, method = method,
      fsm = fsm, optim.method = optim.method
    )

    # extract fsm here
    if (fsm) {
      FSM <- attr(out, "fsm")
    }

    # extract se here
    if (se != "none") {
      SE <- attr(out, "se")
      if (acov != "none") {
        ACOV <- attr(out, "acov")
      }
    }

    # remove dummy lv? (removes attr!)
    out <- lapply(seq_len(lavdata@ngroups), function(g) {
      # determine block
      if (lavdata@nlevels == 1L) {
        bb <- g
      } else {
        bb <- (g - 1) * lavdata@nlevels + level
      }
      lv.idx <- c(
        lavmodel@ov.y.dummy.lv.idx[[bb]],
        lavmodel@ov.x.dummy.lv.idx[[bb]]
      )
      ret <- out[[g]]
      if (length(lv.idx) > 0L) {
        ret <- out[[g]][, -lv.idx, drop = FALSE]
      }
      ret
    })

    # we need to remove the dummy's before we transform
    if (fsm) {
      FSM <- lapply(seq_len(lavdata@ngroups), function(g) {
        # determine block
        if (lavdata@nlevels == 1L) {
          bb <- g
        } else {
          bb <- (g - 1) * lavdata@nlevels + level
        }
        lv.idx <- c(
          lavmodel@ov.y.dummy.lv.idx[[bb]],
          lavmodel@ov.x.dummy.lv.idx[[bb]]
        )
        # ov.idx <- lavmodel@ov.x.dummy.ov.idx[[bb]]
        # or should we use pta$vidx$ov.ind?
        ov.ind <- lavpta$vidx$ov.ind[[bb]]
        ret <- FSM[[g]]
        if (length(lv.idx) > 0L) {
          if (is.matrix(FSM[[g]])) {
            ret <- FSM[[g]][-lv.idx, ov.ind, drop = FALSE]
          } else if (is.list(FSM[[g]])) {
            FSM[[g]] <- lapply(FSM[[g]], function(x) {
              ret <- x[-lv.idx, ov.ind, drop = FALSE]
              ret
            })
          }
        }
        ret
      })
    }

    # new in 0.6-16
    # we assume the dummy lv's have already been removed
    if (transform) {
      VETA <- computeVETA(lavmodel = lavmodel, remove.dummy.lv = TRUE)
      EETA <- computeEETA(
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats, remove.dummy.lv = TRUE
      )
      out <- lapply(seq_len(lavdata@ngroups), function(g) {
        # determine block
        if (lavdata@nlevels == 1L) {
          bb <- g
        } else {
          bb <- (g - 1) * lavdata@nlevels + level
        }

        FS.centered <- scale(out[[g]],
          center = TRUE,
          scale = FALSE
        )
        FS.cov <- crossprod(FS.centered) / nrow(FS.centered)
        FS.cov.inv <- try(solve(FS.cov), silent = TRUE)
        if (inherits(FS.cov.inv, "try-error")) {
          lav_msg_warn(
            gettext("could not invert (co)variance matrix of factor scores;
                    returning original factor scores."))
          return(out[[g]])
        }
        fs.inv.sqrt <- lav_matrix_symmetric_sqrt(FS.cov.inv)
        veta.sqrt <- lav_matrix_symmetric_sqrt(VETA[[g]])
        if (fsm) {
          # change FSM
          FSM[[g]] <<- veta.sqrt %*% fs.inv.sqrt %*% FSM[[g]]
        }
        tmp <- FS.centered %*% fs.inv.sqrt %*% veta.sqrt
        ret <- t(t(tmp) + drop(EETA[[g]]))

        ret
      })
    }

    # new in 0.6-17
    if (mdist) {
      VETA <- computeVETA(lavmodel = lavmodel, remove.dummy.lv = TRUE)
      EETA <- computeEETA(
        lavmodel = lavmodel,
        lavsamplestats = lavsamplestats, remove.dummy.lv = TRUE
      )
      MDIST <- lapply(seq_len(lavdata@ngroups), function(g) {
        A <- FSM[[g]]
        Sigma <- lavimplied$cov[[g]]
        if (transform) {
          fs.cov <- VETA[[g]]
        } else {
          fs.cov <- A %*% Sigma %*% t(A)
        }
        fs.cov.inv <- solve(fs.cov)
        # Mahalobis distance
        fs.c <- t(t(out[[g]]) - EETA[[g]]) # center
        df.squared <- rowSums((fs.c %*% fs.cov.inv) * fs.c)
        ret <- df.squared # squared!
        ret
      })
    }

    # append original/new data? (also remove attr)
    if (append.data && level == 1L) {
      out <- lapply(seq_len(lavdata@ngroups), function(g) {
        ret <- cbind(out[[g]], data.obs[[g]])
        ret
      })
    }

    if (se != "none") {
      SE <- lapply(seq_len(lavdata@ngroups), function(g) {
        # determine block
        if (lavdata@nlevels == 1L) {
          bb <- g
        } else {
          bb <- (g - 1) * lavdata@nlevels + level
        }
        lv.idx <- c(
          lavmodel@ov.y.dummy.lv.idx[[bb]],
          lavmodel@ov.x.dummy.lv.idx[[bb]]
        )
        ret <- SE[[g]]
        if (length(lv.idx) > 0L) {
          ret <- SE[[g]][, -lv.idx, drop = FALSE]
        }
        ret
      })
      if (acov != "none") {
        ACOV <- lapply(seq_len(lavdata@ngroups), function(g) {
          # determine block
          if (lavdata@nlevels == 1L) {
            bb <- g
          } else {
            bb <- (g - 1) * lavdata@nlevels + level
          }
          lv.idx <- c(
            lavmodel@ov.y.dummy.lv.idx[[bb]],
            lavmodel@ov.x.dummy.lv.idx[[bb]]
          )
          ret <- ACOV[[g]]
          if (length(lv.idx) > 0L) {
            if (is.matrix(ACOV[[g]])) {
              ret <- ACOV[[g]][-lv.idx, -lv.idx, drop = FALSE]
            } else if (is.list(ACOV[[g]])) {
              ret <- lapply(ACOV[[g]], function(x) {
                ret <- x[-lv.idx, -lv.idx, drop = FALSE]
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

        if (append.data) {
          colnames(out[[g]]) <- c(
            lavpta$vnames$lv[[gg]],
            ov.names[[g]]
          ) # !not gg
        } else {
          colnames(out[[g]]) <- lavpta$vnames$lv[[gg]]
        }

        if (fsm) {
          if (is.null(FSM[[g]])) {
            # skip
          } else if (is.matrix(FSM[[g]])) {
            dimnames(FSM[[g]]) <- list(
              lavpta$vnames$lv[[gg]],
              # ov.names[[g]]) # !not gg
              lavpta$vnames$ov.ind[[gg]]
            )
          } else if (is.list(FSM[[g]])) {
            FSM[[g]] <- lapply(FSM[[g]], function(x) {
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
          if (!is.null(SE[[g]])) {
            colnames(SE[[g]]) <- lavpta$vnames$lv[[gg]]
          }
        }

        if (acov != "none") {
          if (is.null(ACOV[[g]])) {
            # skip
          } else if (is.matrix(ACOV[[g]])) {
            dimnames(ACOV[[g]]) <- list(
              lavpta$vnames$lv[[gg]],
              lavpta$vnames$lv[[gg]]
            )
          } else if (is.list(ACOV[[g]])) {
            ACOV[[g]] <- lapply(ACOV[[g]], function(x) {
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
          names(SE) <- lavdata@group.label
        }
        if (acov != "none") {
          names(ACOV) <- lavdata@group.label
        }
      }
    } # label

    # yhat: estimated value for the observed indicators, given (estimated)
    # factor scores
    # resid: y - yhat
  } else if (type %in% c("yhat", "resid")) {
    resid.flag <- type == "resid"
    out <- lav_predict_yhat(
      lavobject = NULL, lavmodel = lavmodel,
      lavdata = lavdata, lavsamplestats = lavsamplestats,
      lavimplied = lavimplied,
      data.obs = data.obs, eXo = eXo,
      ETA = ETA, method = method, optim.method = optim.method,
      fsm = fsm,
      resid.flag = resid.flag
    )
    if (fsm) {
      FSM <- attr(out, "fsm")
    }

    # label?
    if (label) {
      for (g in seq_len(lavdata@ngroups)) {
        colnames(out[[g]]) <- lavpta$vnames$ov[[g]]
      }
    }

    # mdist
    if (mdist) {
      LAMBDA <- computeLAMBDA(
        lavmodel = lavmodel,
        remove.dummy.lv = FALSE
      )
      MDIST <- lapply(seq_len(lavdata@ngroups), function(g) {
        Sigma <- lavimplied$cov[[g]]
        LA <- LAMBDA[[g]]
        if (type == "resid") {
          ILA <- diag(ncol(Sigma)) - LA %*% FSM[[g]]
          Omega.e <- ILA %*% Sigma %*% t(ILA)
          eig <- eigen(Omega.e, symmetric = TRUE)
          A <- eig$vectors[, seq_len(nrow(LA) - ncol(LA)),
            drop = FALSE
          ]
        } else if (type == "yhat") {
          LAA <- LA %*% FSM[[g]]
          Omega.e <- LAA %*% Sigma %*% t(LAA)
          eig <- eigen(Omega.e, symmetric = TRUE)
          A <- eig$vectors[, seq_len(ncol(LA)), drop = FALSE]
        }
        outA <- apply(out[[g]], 1L, function(x) {
          colSums(A * x, na.rm = TRUE)
        })
        if (is.matrix(outA)) {
          outA <- t(outA)
        } else {
          outA <- as.matrix(outA)
        }
        # if(lavmodel@meanstructure) {
        #     est.mean <- drop(t(lavimplied$mean[[g]]) %*% A)
        #     if(type == "resid") {
        #         obs.mean <- drop(lavh1$implied$mean[[g]] %*% A)
        #         est.mean <- drop(t(lavimplied$mean[[g]]) %*% A)
        #         outA.mean <- obs.mean - est.mean
        #     } else if(type == "yhat") {
        #         outA.mean <- est.mean
        #     }
        # } else {
        #     outA.mean <- colMeans(outA)
        # }
        outA.cov <- t(A) %*% Omega.e %*% A
        outA.cov.inv <- solve(outA.cov)
        # Mahalobis distance
        # outA.c <- t( t(outA) - outA.mean ) # center
        outA.c <- outA
        df.squared <- rowSums((outA.c %*% outA.cov.inv) * outA.c)
        ret <- df.squared # squared!
        ret
      })
    }


    # density for each observed item, given (estimated) factor scores
  } else if (type == "fy") {
    out <- lav_predict_fy(
      lavobject = NULL, lavmodel = lavmodel,
      lavdata = lavdata, lavsamplestats = lavsamplestats,
      lavimplied = lavimplied,
      data.obs = data.obs, eXo = eXo,
      ETA = ETA, method = method, optim.method = optim.method
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

  if (lavdata@ngroups == 1L && drop.list.single.group) {
    res <- out[[1L]]
  } else {
    res <- out
  }

  # assemble multiple groups into a single data.frame? (new in 0.6-4)
  if (lavdata@ngroups > 1L && assemble) {
    if (!is.null(newdata)) {
      lavdata <- newData
    }
    DATA <- matrix(as.numeric(NA),
      nrow = sum(unlist(lavdata@norig)),
      ncol = ncol(out[[1L]])
    ) # assume == per g
    colnames(DATA) <- colnames(out[[1L]])
    for (g in seq_len(lavdata@ngroups)) {
      DATA[lavdata@case.idx[[g]], ] <- out[[g]]
    }
    DATA <- as.data.frame(DATA, stringsAsFactors = FALSE)

    if (!is.null(newdata)) {
      DATA[, lavdata@group] <- newdata[, lavdata@group]
    } else {
      # add group
      DATA[, lavdata@group] <- rep(as.character(NA), nrow(DATA))
      if (lavdata@missing == "listwise") {
        # we will loose the group label of omitted variables!
        DATA[unlist(lavdata@case.idx), lavdata@group] <-
          rep(lavdata@group.label, unlist(lavdata@nobs))
      } else {
        DATA[unlist(lavdata@case.idx), lavdata@group] <-
          rep(lavdata@group.label, unlist(lavdata@norig))
      }
    }

    res <- DATA
  }

  if (fsm && type == "lv") {
    attr(res, "fsm") <- FSM
  }

  if (mdist) {
    attr(res, "mdist") <- MDIST
  }

  if (se != "none") {
    attr(res, "se") <- SE
    # return full sampling covariance matrix?
    if (acov == "standard") {
      attr(res, "acov") <- ACOV
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
                            data.obs = NULL, eXo = NULL,
                            # options
                            method = "EBM",
                            fsm = FALSE,
                            se = "none", acov = "none",
                            level = 1L,
                            optim.method = "bfgs") {
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
        data.obs = data.obs, eXo = eXo, fsm = fsm
      )
    } else if (method == "ml") {
      out <- lav_predict_eta_bartlett(
        lavobject = lavobject,
        lavmodel = lavmodel, lavdata = lavdata,
        lavimplied = lavimplied, se = se, acov = acov,
        level = level, lavsamplestats = lavsamplestats,
        data.obs = data.obs, eXo = eXo, fsm = fsm
      )
    } else {
      lav_msg_stop(gettextf("unkown method: %s.", method))
    }
  } else {
    if (method == "ebm") {
      out <- lav_predict_eta_ebm_ml(
        lavobject = lavobject,
        lavmodel = lavmodel, lavdata = lavdata,
        lavsamplestats = lavsamplestats, se = se, acov = acov,
        level = level, data.obs = data.obs, eXo = eXo,
        ML = FALSE, optim.method = optim.method
      )
    } else if (method == "ml") {
      out <- lav_predict_eta_ebm_ml(
        lavobject = lavobject,
        lavmodel = lavmodel, lavdata = lavdata,
        lavsamplestats = lavsamplestats, se = se, acov = acov,
        level = level, data.obs = data.obs, eXo = eXo,
        ML = TRUE, optim.method = optim.method
      )
    } else {
      lav_msg_stop(gettextf("unkown method: %s.", method))
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
                                   data.obs = NULL, eXo = NULL,
                                   se = "none", acov = "none", level = 1L,
                                   fsm = FALSE) {
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

  if (is.null(data.obs)) {
    data.obs <- lavdata@X
    newdata.flag <- FALSE
  } else {
    newdata.flag <- TRUE
  }
  # eXo not needed

  # missings? and missing = "ml"?
  if (lavdata@missing %in% c("ml", "ml.x")) {
    if (newdata.flag) {
      MP <- vector("list", lavdata@ngroups)
      for (g in seq_len(lavdata@ngroups)) {
        MP[[g]] <- lav_data_missing_patterns(data.obs[[g]])
      }
    } else {
      MP <- lavdata@Mp
    }
  }

  LAMBDA <- computeLAMBDA(lavmodel = lavmodel, remove.dummy.lv = FALSE)
  Sigma.hat <- lavimplied$cov
  Sigma.inv <- lapply(Sigma.hat, MASS::ginv)
  VETA <- computeVETA(lavmodel = lavmodel)
  EETA <- computeEETA(lavmodel = lavmodel, lavsamplestats = lavsamplestats)
  EY <- computeEY(lavmodel = lavmodel, lavsamplestats = lavsamplestats)

  FS <- vector("list", length = lavdata@ngroups)
  if (fsm) {
    FSM <- vector("list", length = lavdata@ngroups)
  }

  if (acov != "none") {
    se <- acov # ACOV implies SE
  }
  if (se != "none") {
    SE <- vector("list", length = lavdata@ngroups)
    # return full sampling covariance matrix?
    if (acov != "none") {
      ACOV <- vector("list", length = lavdata@ngroups)
    }
  }

  for (g in 1:lavdata@ngroups) {
    if (lavdata@nlevels > 1L) {
      Lp <- lavdata@Lp[[g]]
      YLp <- lavsamplestats@YLp[[g]]

      # implied for this group
      group.idx <- (g - 1) * lavdata@nlevels + seq_len(lavdata@nlevels)
      implied.group <- lapply(lavimplied, function(x) x[group.idx])

      # random effects (=random intercepts or cluster means)
      out <- lav_mvnorm_cluster_implied22l(
        Lp = Lp,
        implied = implied.group
      )
      MB.j <- lav_mvnorm_cluster_em_estep_ranef(
        YLp = YLp, Lp = Lp,
        sigma.w = out$sigma.w, sigma.b = out$sigma.b,
        sigma.zz = out$sigma.zz, sigma.yz = out$sigma.yz,
        mu.z = out$mu.z, mu.w = out$mu.w, mu.b = out$mu.b,
        se = FALSE
      )

      ov.idx <- Lp$ov.idx

      if (level == 1L) {
        data.W <- data.obs[[g]][, ov.idx[[1]]]
        data.B <- MB.j[Lp$cluster.idx[[2]], , drop = FALSE]

        # center
        data.obs.g <- data.W - data.B
      } else if (level == 2L) {
        Data.B <- matrix(0,
          nrow = nrow(MB.j),
          ncol = ncol(data.obs[[g]])
        )
        Data.B[, ov.idx[[1]]] <- MB.j
        between.idx <- Lp$between.idx[[2 * g]]
        if (length(between.idx) > 0L) {
          Data.B[, between.idx] <- data.obs[[g]][
            !duplicated(Lp$cluster.idx[[2]]),
            between.idx
          ]
        }
        data.obs.g <- Data.B[, ov.idx[[2]]]
      } else {
        lav_msg_stop(gettext("only 2 levels are supported"))
      }

      gg <- (g - 1) * lavdata@nlevels + level
      VETA.g <- VETA[[gg]]
      EETA.g <- EETA[[gg]]
      LAMBDA.g <- LAMBDA[[gg]]
      EY.g <- EY[[gg]]
      Sigma.inv.g <- Sigma.inv[[gg]]
    } else {
      data.obs.g <- data.obs[[g]]
      VETA.g <- VETA[[g]]
      EETA.g <- EETA[[g]]
      LAMBDA.g <- LAMBDA[[g]]
      EY.g <- EY[[g]]
      Sigma.inv.g <- Sigma.inv[[g]]
    }

    nfac <- ncol(VETA[[g]])
    if (nfac == 0L) {
      FS[[g]] <- matrix(0, lavdata@nobs[[g]], nfac)
      next
    }

    # center data
    Yc <- t(t(data.obs.g) - EY.g)

    # global factor score coefficient matrix 'C'
    FSC <- VETA.g %*% t(LAMBDA.g) %*% Sigma.inv.g

    # store fsm?
    if (fsm) {
      FSM.g <- FSC
    }

    # compute factor scores
    if (lavdata@missing %in% c("ml", "ml.x")) {
      # missing patterns for this group
      Mp <- MP[[g]]

      # factor scores container
      FS.g <- matrix(as.numeric(NA), nrow(Yc), ncol = length(EETA.g))

      # if(fsm) {
      #    FSM.g <- vector("list", length = Mp$npatterns)
      # }

      if (se == "standard") {
        SE.g <- matrix(as.numeric(NA), nrow(Yc), ncol = length(EETA.g))
      }

      if (acov == "standard") {
        ACOV.g <- vector("list", length = Mp$npatterns)
      }

      # compute FSC per pattern
      for (p in seq_len(Mp$npatterns)) {
        var.idx <- Mp$pat[p, ] # observed
        na.idx <- which(!var.idx) # missing

        # extract observed data for these (centered) cases
        Oc <- Yc[Mp$case.idx[[p]], Mp$pat[p, ], drop = FALSE]

        # invert Sigma (Sigma_22, observed part only) for this pattern
        Sigma_22.inv <- try(lav_matrix_symmetric_inverse_update(
          S.inv =
            Sigma.inv.g, rm.idx = na.idx,
          logdet = FALSE
        ), silent = TRUE)
        if (inherits(Sigma_22.inv, "try-error")) {
          lav_msg_stop(gettext("Sigma_22.inv cannot be inverted"))
        }

        lambda <- LAMBDA.g[var.idx, , drop = FALSE]
        FSC <- VETA.g %*% t(lambda) %*% Sigma_22.inv

        # FSM?
        # if(fsm) {
        #    tmp <- matrix(as.numeric(NA), nrow = ncol(lambda),
        #                  ncol = ncol(Yc))
        #    tmp[,var.idx] <- FSC
        #    FSM.g[[p]] <- tmp
        # }

        # factor score for this pattern
        FS.g[Mp$case.idx[[p]], ] <- t(FSC %*% t(Oc) + EETA.g)

        # SE?
        if (se == "standard") {
          tmp <- (VETA.g - VETA.g %*% t(lambda) %*%
            Sigma_22.inv %*% lambda %*% VETA.g)
          tmp.d <- diag(tmp)
          tmp.d[tmp.d < 1e-05] <- as.numeric(NA)

          # all cases in this pattern get the same SEs
          SE.g[Mp$case.idx[[p]], ] <- matrix(sqrt(tmp.d),
            nrow = length(Mp$case.idx[[p]]),
            ncol = ncol(SE.g), byrow = TRUE
          )
        }

        # ACOV?
        if (acov == "standard") {
          ACOV.g[[p]] <- tmp # for this pattern
        }
      } # p
    } else {
      # compute factor scores
      FS.g <- t(FSC %*% t(Yc) + EETA.g)
    }

    # replace values in dummy lv's by their observed counterpart
    if (length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L && level == 1L) {
      FS.g[, lavmodel@ov.y.dummy.lv.idx[[g]]] <-
        data.obs.g[, lavmodel@ov.y.dummy.ov.idx[[g]], drop = FALSE]
    }
    if (length(lavmodel@ov.x.dummy.lv.idx[[g]]) > 0L && level == 1L) {
      FS.g[, lavmodel@ov.x.dummy.lv.idx[[g]]] <-
        data.obs.g[, lavmodel@ov.x.dummy.ov.idx[[g]], drop = FALSE]
    }

    FS[[g]] <- FS.g

    # FSM
    if (fsm) {
      FSM[[g]] <- FSM.g
    }

    # standard error
    if (se == "standard") {
      if (lavdata@missing %in% c("ml", "ml.x")) {
        SE[[g]] <- SE.g
        if (acov == "standard") {
          ACOV[[g]] <- ACOV.g
        }
      } else { # complete data
        tmp <- (VETA.g - VETA.g %*% t(LAMBDA.g) %*%
          Sigma.inv.g %*% LAMBDA.g %*% VETA.g)
        tmp.d <- diag(tmp)
        tmp.d[tmp.d < 1e-05] <- as.numeric(NA)
        SE[[g]] <- matrix(sqrt(tmp.d), nrow = 1L)

        # return full sampling covariance matrix?
        if (acov == "standard") {
          ACOV[[g]] <- tmp
        }
      }
    } # se = "standard"
  } # g

  if (fsm) {
    attr(FS, "fsm") <- FSM
  }
  if (se != "none") {
    attr(FS, "se") <- SE
    # return full sampling covariance matrix?
    if (acov == "standard") {
      attr(FS, "acov") <- ACOV
    }
  }

  FS
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
                                     data.obs = NULL, eXo = NULL,
                                     se = "none", acov = "none", level = 1L,
                                     fsm = FALSE) {
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

  if (is.null(data.obs)) {
    data.obs <- lavdata@X
    newdata.flag <- FALSE
  } else {
    newdata.flag <- TRUE
  }
  # eXo not needed

  # missings? and missing = "ml"?
  if (lavdata@missing %in% c("ml", "ml.x")) {
    if (newdata.flag) {
      MP <- vector("list", lavdata@ngroups)
      for (g in seq_len(lavdata@ngroups)) {
        MP[[g]] <- lav_data_missing_patterns(data.obs[[g]])
      }
    } else {
      MP <- lavdata@Mp
    }
  }

  LAMBDA <- computeLAMBDA(lavmodel = lavmodel, remove.dummy.lv = FALSE)
  Sigma <- lavimplied$cov
  Sigma.inv <- lapply(lavimplied$cov, MASS::ginv)
  VETA <- computeVETA(lavmodel = lavmodel) # for se only
  EETA <- computeEETA(lavmodel = lavmodel, lavsamplestats = lavsamplestats)
  EY <- computeEY(lavmodel = lavmodel, lavsamplestats = lavsamplestats)

  FS <- vector("list", length = lavdata@ngroups)
  if (fsm) {
    FSM <- vector("list", length = lavdata@ngroups)
  }

  if (acov != "none") se <- acov # ACOV implies SE
  if (se != "none") {
    SE <- vector("list", length = lavdata@ngroups)
    # return full sampling covariance matrix?
    if (acov != "none") {
      ACOV <- vector("list", length = lavdata@ngroups)
    }
  }

  for (g in 1:lavdata@ngroups) {
    if (lavdata@nlevels > 1L) {
      Lp <- lavdata@Lp[[g]]
      YLp <- lavsamplestats@YLp[[g]]

      # implied for this group
      group.idx <- (g - 1) * lavdata@nlevels + seq_len(lavdata@nlevels)
      implied.group <- lapply(lavimplied, function(x) x[group.idx])

      # random effects (=random intercepts or cluster means)
      # NOTE: is the 'ML' way not simply using the observed cluster
      #       means?
      out <- lav_mvnorm_cluster_implied22l(
        Lp = Lp,
        implied = implied.group
      )
      MB.j <- lav_mvnorm_cluster_em_estep_ranef(
        YLp = YLp, Lp = Lp,
        sigma.w = out$sigma.w, sigma.b = out$sigma.b,
        sigma.zz = out$sigma.zz, sigma.yz = out$sigma.yz,
        mu.z = out$mu.z, mu.w = out$mu.w, mu.b = out$mu.b,
        se = FALSE
      )

      ov.idx <- Lp$ov.idx

      if (level == 1L) {
        data.W <- data.obs[[g]][, ov.idx[[1]]]
        data.B <- MB.j[Lp$cluster.idx[[2]], , drop = FALSE]

        # center
        data.obs.g <- data.W - data.B
      } else if (level == 2L) {
        Data.B <- matrix(0,
          nrow = nrow(MB.j),
          ncol = ncol(data.obs[[g]])
        )
        Data.B[, ov.idx[[1]]] <- MB.j
        between.idx <- Lp$between.idx[[2 * g]]
        if (length(between.idx) > 0L) {
          Data.B[, between.idx] <- data.obs[[g]][
            !duplicated(Lp$cluster.idx[[2]]),
            between.idx
          ]
        }
        data.obs.g <- Data.B[, ov.idx[[2]]]
      } else {
        lav_msg_stop(gettext("only 2 levels are supported"))
      }

      gg <- (g - 1) * lavdata@nlevels + level
      VETA.g <- VETA[[gg]]
      EETA.g <- EETA[[gg]]
      LAMBDA.g <- LAMBDA[[gg]]
      EY.g <- EY[[gg]]
      Sigma.inv.g <- Sigma.inv[[gg]]
    } else {
      data.obs.g <- data.obs[[g]]
      VETA.g <- VETA[[g]]
      EETA.g <- EETA[[g]]
      LAMBDA.g <- LAMBDA[[g]]
      EY.g <- EY[[g]]
      Sigma.inv.g <- Sigma.inv[[g]]
    }

    nfac <- length(EETA[[g]])
    if (nfac == 0L) {
      FS[[g]] <- matrix(0, lavdata@nobs[[g]], nfac)
      next
    }

    # center data
    Yc <- t(t(data.obs.g) - EY.g)

    # global factor score coefficient matrix 'C'
    FSC <- (MASS::ginv(t(LAMBDA.g) %*% Sigma.inv.g %*% LAMBDA.g)
    %*% t(LAMBDA.g) %*% Sigma.inv.g)

    # store fsm?
    if (fsm) {
      # store fsm?
      FSM.g <- FSC
    }

    # compute factor scores
    if (lavdata@missing %in% c("ml", "ml.x")) {
      # missing patterns for this group
      Mp <- MP[[g]]

      # factor scores container
      FS.g <- matrix(as.numeric(NA), nrow(Yc), ncol = length(EETA.g))

      # if(fsm) {
      #    FSM.g <- vector("list", length = Mp$npatterns)
      # }

      if (se == "standard") {
        SE.g <- matrix(as.numeric(NA), nrow(Yc), ncol = length(EETA.g))
      }

      if (acov == "standard") {
        ACOV.g <- vector("list", length = Mp$npatterns)
      }

      # compute FSC per pattern
      for (p in seq_len(Mp$npatterns)) {
        var.idx <- Mp$pat[p, ] # observed
        na.idx <- which(!var.idx) # missing

        # extract observed data for these (centered) cases
        Oc <- Yc[Mp$case.idx[[p]], Mp$pat[p, ], drop = FALSE]

        # invert Sigma (Sigma_22, observed part only) for this pattern
        Sigma_22.inv <- try(lav_matrix_symmetric_inverse_update(
          S.inv =
            Sigma.inv.g, rm.idx = na.idx,
          logdet = FALSE
        ), silent = TRUE)
        if (inherits(Sigma_22.inv, "try-error")) {
          lav_msg_stop(gettext("Sigma_22.inv cannot be inverted"))
        }

        lambda <- LAMBDA.g[var.idx, , drop = FALSE]
        FSC <- (MASS::ginv(t(lambda) %*% Sigma_22.inv %*% lambda)
        %*% t(lambda) %*% Sigma_22.inv)

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
        zero.idx <- which(apply(FSC, 1L, function(x) all(x == 0)))
        if (length(zero.idx) > 0L) {
          FSC[zero.idx, ] <- NA
        }

        # FSM?
        # if(fsm) {
        #    tmp <- matrix(as.numeric(NA), nrow = ncol(lambda),
        #                  ncol = ncol(Yc))
        #    tmp[,var.idx] <- FSC
        #    FSM.g[[p]] <- tmp
        # }

        # factor scores for this pattern
        FS.g[Mp$case.idx[[p]], ] <- t(FSC %*% t(Oc) + EETA.g)

        # SE?
        if (se == "standard") {
          tmp <- (MASS::ginv(t(lambda) %*% Sigma_22.inv %*% lambda)
          - VETA.g)
          tmp.d <- diag(tmp)
          tmp.d[tmp.d < 1e-05] <- as.numeric(NA)

          # all cases in this pattern get the same SEs
          SE.g[Mp$case.idx[[p]], ] <- matrix(sqrt(tmp.d),
            nrow = length(Mp$case.idx[[p]]),
            ncol = ncol(SE.g), byrow = TRUE
          )
        }

        # ACOV?
        if (acov == "standard") {
          ACOV.g[[p]] <- tmp # for this pattern
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
      FS.g <- t(FSC %*% t(Yc) + EETA.g)
    }

    # replace values in dummy lv's by their observed counterpart
    if (length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L && level == 1L) {
      FS.g[, lavmodel@ov.y.dummy.lv.idx[[g]]] <-
        data.obs[[g]][, lavmodel@ov.y.dummy.ov.idx[[g]], drop = FALSE]
    }
    if (length(lavmodel@ov.x.dummy.lv.idx[[g]]) > 0L && level == 1L) {
      FS.g[, lavmodel@ov.x.dummy.lv.idx[[g]]] <-
        data.obs[[g]][, lavmodel@ov.x.dummy.ov.idx[[g]], drop = FALSE]
    }

    FS[[g]] <- FS.g

    # FSM
    if (fsm) {
      FSM[[g]] <- FSM.g
    }

    # standard error
    if (se == "standard") {
      if (lavdata@missing %in% c("ml", "ml.x")) {
        SE[[g]] <- SE.g
        if (acov == "standard") {
          ACOV[[g]] <- ACOV.g
        }
      } else { # complete data

        # the traditional formula is:
        #     solve(t(lambda) %*% solve(theta) %*% lambda)
        # but we replace it by
        #     solve( t(lambda) %*% solve(sigma) %*% lambda ) - psi
        # to handle negative variances
        # in addition, we use ginv
        tmp <- (MASS::ginv(t(LAMBDA.g) %*% Sigma.inv.g %*% LAMBDA.g)
        - VETA.g)
        tmp.d <- diag(tmp)
        tmp.d[tmp.d < 1e-05] <- as.numeric(NA)
        SE[[g]] <- matrix(sqrt(tmp.d), nrow = 1L)

        # return full sampling covariance matrix?
        if (acov == "standard") {
          ACOV[[g]] <- tmp
        }
      }
    } # se
  } # g

  if (fsm) {
    attr(FS, "fsm") <- FSM
  }
  if (se != "none") {
    attr(FS, "se") <- SE
    # return full sampling covariance matrix?
    if (acov == "standard") {
      attr(FS, "acov") <- ACOV
    }
  }

  FS
}

# factor scores - EBM or ML
lav_predict_eta_ebm_ml <- function(lavobject = NULL, # for convenience
                                   # sub objects
                                   lavmodel = NULL, lavdata = NULL,
                                   lavsamplestats = NULL,
                                   # optional new data
                                   data.obs = NULL, eXo = NULL,
                                   se = "none", acov = "none", level = 1L,
                                   ML = FALSE,
                                   optim.method = "bfgs") {
  optim.method <- tolower(optim.method)

  stopifnot(optim.method %in% c("nlminb", "bfgs"))

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
  if (is.null(data.obs)) {
    data.obs <- lavdata@X
  }
  if (is.null(eXo)) {
    eXo <- lavdata@eXo
  }

  # se?
  if (acov != "none") {
    se <- acov # ACOV implies SE
  }
  # if(se != "none") {
  #    warning("lavaan WARNING: standard errors are not available (yet) for the non-normal case")
  # }

  VETAx <- computeVETAx(lavmodel = lavmodel)
  VETAx.inv <- VETAx
  for (g in seq_len(lavdata@ngroups)) {
    if (nrow(VETAx[[g]]) > 0L) {
      VETAx.inv[[g]] <- solve(VETAx[[g]])
    }
  }
  EETAx <- computeEETAx(
    lavmodel = lavmodel, lavsamplestats = lavsamplestats,
    eXo = eXo, nobs = lapply(data.obs, NROW),
    remove.dummy.lv = TRUE
  ) ## FIXME?
  TH <- computeTH(lavmodel = lavmodel, delta = FALSE)
  THETA <- computeTHETA(lavmodel = lavmodel)

  # check for zero entries in THETA (new in 0.6-4)
  for (g in seq_len(lavdata@ngroups)) {
    if (any(diag(THETA[[g]]) == 0)) {
      lav_msg_stop(gettext(
        "(residual) variance matrix THETA contains zero elements
        on the diagonal."))
    }
  }

  # local objective function: x = lv values
  f.eta.i <- function(x, y.i, x.i, mu.i) {
    # add 'dummy' values (if any) for ov.y
    if (length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L) {
      x2 <- c(x, data.obs[[g]][i,
        lavmodel@ov.y.dummy.ov.idx[[g]],
        drop = FALSE
      ])
    } else {
      x2 <- x
    }

    # conditional density of y, given eta.i(=x)
    log.fy <- lav_predict_fy_eta.i(
      lavmodel = lavmodel,
      lavdata = lavdata,
      lavsamplestats = lavsamplestats,
      y.i = y.i,
      x.i = x.i,
      eta.i = matrix(x2, nrow = 1L), # <---- eta!
      theta.sd = theta.sd,
      th = th,
      th.idx = th.idx,
      log = TRUE
    )

    if (ML) {
      # NOTE: 'true' ML is simply  -1*sum(log.fy)
      #  - but there is no upper/lower bound for the extrema:
      #    a pattern of all (in)correct drives the 'theta' parameter
      #    towards +/- Inf
      # - therefore, we add a vague prior, just to stabilize
      #
      diff <- t(x) - mu.i
      V <- diag(length(x)) * 1e-05
      tmp <- as.numeric(0.5 * diff %*% V %*% t(diff))
      out <- 1 + tmp - sum(log.fy, na.rm = TRUE)
    } else {
      diff <- t(x) - mu.i
      V <- VETAx.inv[[g]]
      tmp <- as.numeric(0.5 * diff %*% V %*% t(diff))
      out <- tmp - sum(log.fy, na.rm = TRUE)
    }
    out
  }

  FS <- vector("list", length = lavdata@ngroups)
  for (g in seq_len(lavdata@ngroups)) {
    nfac <- ncol(VETAx[[g]])
    nfac2 <- nfac
    if (length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L) {
      nfac2 <- nfac2 + length(lavmodel@ov.y.dummy.lv.idx[[g]])
    }
    FS[[g]] <- matrix(as.numeric(NA), nrow(data.obs[[g]]), nfac2)

    # special case: no regular lv's
    if (nfac == 0) {
      # impute dummy ov.y (if any)
      FS[[g]][, lavmodel@ov.y.dummy.ov.idx[[g]]] <-
        data.obs[[g]][, lavmodel@ov.y.dummy.ov.idx[[g]], drop = FALSE]
      next
    }

    ## FIXME: factor scores not identical (but close) to Mplus
    #         if delta elements not equal to 1??
    mm.in.group <- 1:lavmodel@nmat[g] + cumsum(c(0, lavmodel@nmat))[g]
    MLIST <- lavmodel@GLIST[mm.in.group]

    # check for negative values
    neg.var.idx <- which(diag(THETA[[g]]) < 0)
    if (length(neg.var.idx) > 0) {
      lav_msg_warn(
        gettext("factor scores could not be computed due to at least
                one negative (residual) variance"))
      next
    }

    # common values
    theta.sd <- sqrt(diag(THETA[[g]]))
    th <- TH[[g]]
    th.idx <- lavmodel@th.idx[[g]]

    # casewise for now
    N <- nrow(data.obs[[g]])
    for (i in 1:N) {
      # eXo?
      if (!is.null(eXo[[g]])) {
        x.i <- eXo[[g]][i, , drop = FALSE]
      } else {
        x.i <- NULL
      }
      mu.i <- EETAx[[g]][i, , drop = FALSE]
      y.i <- data.obs[[g]][i, , drop = FALSE]

      ### DEBUG ONLY:
      # cat("i = ", i, "mu.i = ", mu.i, "\n")

      START <- numeric(nfac) # initial values for eta

      if (!all(is.na(y.i))) {
        # find best values for eta.i
        if (optim.method == "nlminb") {
          out <- nlminb(
            start = START, objective = f.eta.i,
            gradient = NULL, # for now
            control = list(rel.tol = 1e-8),
            y.i = y.i, x.i = x.i, mu.i = mu.i
          )
        } else if (optim.method == "bfgs") {
          out <- optim(
            par = START, fn = f.eta.i,
            gr = NULL,
            control = list(reltol = 1e-8, fnscale = 1.1),
            method = "BFGS",
            y.i = y.i, x.i = x.i, mu.i = mu.i
          )
        }
        if (out$convergence == 0L) {
          eta.i <- out$par
        } else {
          eta.i <- rep(as.numeric(NA), nfac)
        }
      } else {
        eta.i <- rep(as.numeric(NA), nfac)
      }

      # add dummy ov.y lv values
      if (length(lavmodel@ov.y.dummy.lv.idx[[g]]) > 0L) {
        eta.i <- c(eta.i, data.obs[[g]][i,
          lavmodel@ov.y.dummy.ov.idx[[g]],
          drop = FALSE
        ])
      }

      FS[[g]][i, ] <- eta.i
    }
  }

  FS
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
                             data.obs = NULL, eXo = NULL,
                             # ETA values
                             ETA = NULL,
                             # options
                             method = "EBM",
                             duplicate = FALSE,
                             optim.method = "bfgs",
                             fsm = FALSE,
                             resid.flag = FALSE) {
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
  if (is.null(data.obs)) {
    data.obs <- lavdata@X
  }
  if (is.null(eXo)) {
    eXo <- lavdata@eXo
  }

  # do we get values for ETA? If not, use `predict' to get plausible values
  if (is.null(ETA)) {
    ETA <- lav_predict_eta(
      lavobject = NULL, lavmodel = lavmodel,
      lavdata = lavdata, lavsamplestats = lavsamplestats,
      lavimplied = lavimplied,
      data.obs = data.obs, eXo = eXo, method = method,
      optim.method = optim.method, fsm = fsm
    )
    FSM <- attr(ETA, "fsm")
  } else {
    # matrix
    if (is.matrix(ETA)) { # user-specified?
      if (nrow(ETA) == 1L) {
        tmp <- matrix(ETA, lavsamplestats@ntotal, length(ETA),
          byrow = TRUE
        )
      } else if (nrow(ETA) != lavsamplestats@ntotal) {
        lav_msg_stop(gettext("nrow(ETA) != lavsamplestats@ntotal"))
      } else {
        tmp <- ETA
      }
      ETA <- lapply(1:lavdata@ngroups, function(i) tmp[lavdata@case.idx[[i]], ])
      # vector: just 1 row of factor-scores
    } else if (is.numeric(ETA)) {
      # convert to matrix
      tmp <- matrix(ETA, lavsamplestats@ntotal, length(ETA), byrow = TRUE)
      ETA <- lapply(1:lavdata@ngroups, function(i) tmp[lavdata@case.idx[[i]], ])
    } else if (is.list(ETA)) {
      stopifnot(lavdata@ngroups == length(ETA))
    }
  }

  YHAT <- computeYHAT(
    lavmodel = lavmodel, GLIST = NULL,
    lavsamplestats = lavsamplestats, eXo = eXo,
    nobs = lapply(data.obs, NROW),
    ETA = ETA, duplicate = duplicate
  )

  # if conditional.x, paste eXo
  if (lavmodel@categorical && !is.null(eXo)) {
    YHAT <- lapply(seq_len(lavdata@ngroups), function(g) {
      ret <- cbind(YHAT[[g]], eXo[[g]])
      ret
    })
  }

  # residuals? compute y - yhat
  if (resid.flag) {
    RES <- lapply(seq_len(lavdata@ngroups), function(g) {
      ret <- data.obs[[g]] - YHAT[[g]]
      ret
    })
  } else {
    RES <- YHAT
  }

  # fsm?
  if (fsm) {
    attr(RES, "fsm") <- FSM
  }

  RES
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
                           data.obs = NULL, eXo = NULL,
                           # ETA values
                           ETA = NULL,
                           # options
                           method = "EBM",
                           log. = FALSE,
                           optim.method = "bfgs") {
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
  if (is.null(data.obs)) {
    data.obs <- lavdata@X
  }
  if (is.null(eXo)) {
    eXo <- lavdata@eXo
  }

  # we need the YHATs (per group)
  YHAT <- lav_predict_yhat(
    lavobject = NULL, lavmodel = lavmodel,
    lavdata = lavdata, lavsamplestats = lavsamplestats,
    lavimplied = lavimplied,
    data.obs = data.obs, eXo = eXo, ETA = ETA, method = method,
    duplicate = FALSE, optim.method = optim.method
  )

  THETA <- computeTHETA(lavmodel = lavmodel)
  TH <- computeTH(lavmodel = lavmodel, delta = FALSE)

  FY <- vector("list", length = lavdata@ngroups)
  for (g in seq_len(lavdata@ngroups)) {
    FY[[g]] <- lav_predict_fy_internal(
      X = data.obs[[g]], yhat = YHAT[[g]],
      TH = TH[[g]], THETA = THETA[[g]],
      num.idx = lavmodel@num.idx[[g]],
      th.idx = lavmodel@th.idx[[g]],
      link = lavmodel@link, log. = log.
    )
  }

  FY
}


# single group, internal function
lav_predict_fy_internal <- function(X = NULL, yhat = NULL,
                                    TH = NULL, THETA = NULL,
                                    num.idx = NULL, th.idx = NULL,
                                    link = NULL, log. = FALSE) {
  # shortcuts
  theta.var <- diag(THETA)

  # check size YHAT (either 1L or Nobs rows)
  if (!(nrow(yhat) == 1L || nrow(yhat) == nrow(X))) {
    lav_msg_stop(gettext("nrow(YHAT[[g]]) not 1L and not nrow(X))"))
  }

  FY.group <- matrix(0, nrow(X), ncol(X))
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

  ord.idx <- unique(th.idx[th.idx > 0L])

  # first, NUMERIC variables
  if (length(num.idx) > 0L) {
    for (v in num.idx) {
      FY.group[, v] <- dnorm(X[, v],
        # YHAT may change or not per case
        mean = yhat[, v],
        sd   = sqrt(theta.var[v]),
        log  = log.
      )
    }
  }

  # second, ORDERED variables
  for (v in ord.idx) {
    th.y <- TH[th.idx == v]
    TH.Y <- c(-Inf, th.y, Inf)
    ncat <- length(th.y) + 1L
    fy <- numeric(ncat)
    theta.v <- sqrt(theta.var[v])
    yhat.v <- yhat[, v]

    # two cases: yhat.v is a scalar, or has length = nobs
    fy <- matrix(0, nrow = length(yhat.v), ncol = ncat)

    # for each category
    for (k in seq_len(ncat)) {
      if (link == "probit") {
        fy[, k] <- pnorm((TH.Y[k + 1] - yhat.v) / theta.v) -
          pnorm((TH.Y[k] - yhat.v) / theta.v)
      } else if (link == "logit") {
        fy[, k] <- plogis((TH.Y[k + 1] - yhat.v) / theta.v) -
          plogis((TH.Y[k] - yhat.v) / theta.v)
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
    if (log.) {
      fy <- log(fy)
    }

    # case-wise expansion/selection
    if (length(yhat.v) == 1L) {
      # expand category probabilities for all observations
      FY.group[, v] <- fy[1L, X[, v]]
    } else {
      # select correct category probability per observation
      FY.group[, v] <- fy[cbind(seq_len(nrow(fy)), X[, v])]
    }
  } # ord

  FY.group
}



# conditional density y -- assuming independence!!
# f(y_i | eta_i, x_i)
#
# but for a SINGLE observation y_i (and x_i), for given values of eta_i
#
lav_predict_fy_eta.i <- function(lavmodel = NULL, lavdata = NULL,
                                 lavsamplestats = NULL,
                                 y.i = NULL, x.i = NULL,
                                 eta.i = NULL, theta.sd = NULL, g = 1L,
                                 th = NULL, th.idx = NULL, log = TRUE) {
  mm.in.group <- 1:lavmodel@nmat[g] + cumsum(c(0, lavmodel@nmat))[g]
  MLIST <- lavmodel@GLIST[mm.in.group]

  # linear predictor for all items
  YHAT <-
    computeEYetax.LISREL(
      MLIST = MLIST,
      eXo = x.i,
      ETA = eta.i,
      sample.mean = lavsamplestats@mean[[g]],
      ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
      ov.x.dummy.ov.idx = lavmodel@ov.x.dummy.ov.idx[[g]],
      ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]],
      ov.x.dummy.lv.idx = lavmodel@ov.x.dummy.lv.idx[[g]],
      delta = FALSE
    )

  # P(y_i | eta_i, x_i) for all items
  if (all(lavdata@ov$type == "numeric")) {
    # NORMAL case
    FY <- dnorm(y.i, mean = YHAT, sd = theta.sd, log = log)
  } else {
    FY <- numeric(lavmodel@nvar[g])
    for (v in seq_len(lavmodel@nvar[g])) {
      if (lavdata@ov$type[v] == "numeric") {
        ### FIXME!!! we can do all numeric vars at once!!
        FY[v] <- dnorm(y.i[v],
          mean = YHAT[v], sd = theta.sd[v],
          log = log
        )
      } else if (lavdata@ov$type[v] == "ordered") {
        # handle missing value
        if (is.na(y.i[v])) {
          FY[v] <- as.numeric(NA)
        } else {
          th.y <- th[th.idx == v]
          TH.Y <- c(-Inf, th.y, Inf)
          k <- y.i[v]
          p1 <- pnorm((TH.Y[k + 1] - YHAT[v]) / theta.sd[v])
          p2 <- pnorm((TH.Y[k] - YHAT[v]) / theta.sd[v])
          prob <- (p1 - p2)
          if (prob < .Machine$double.eps) {
            prob <- .Machine$double.eps
          }
          if (log) {
            FY[v] <- log(prob)
          } else {
            FY[v] <- prob
          }
        }
      } else {
        lav_msg_stop(gettextf("unknown type: `%1$s' for variable: %2$s",
          lavdata@ov$type[v], lavdata@ov$name[v])
        )
      }
    }
  }

  FY
}
