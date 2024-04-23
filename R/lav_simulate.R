# new version of lavSimulateData (replaced simulateData)
# from lavaan 0.6-1
# YR 23 March 2018
#
# - calls lavaan directly to get model-implied statistics
# - allows for groups with different sets of variables
# -


lavSimulateData <- function(model = NULL,
                            cmd.pop = "sem",
                            ...,
                            # data properties
                            sample.nobs = 1000L,
                            cluster.idx = NULL,
                            # control
                            empirical = FALSE,
                            # output
                            add.labels = TRUE,
                            return.fit = FALSE,
                            output = "data.frame") {
  # dotdotdot
  dotdotdot <- list(...)
  dotdotdot.orig <- dotdotdot

  # remove/override some options
  dotdotdot$verbose <- FALSE
  dotdotdot$debug <- FALSE
  dotdotdot$data <- NULL
  dotdotdot$sample.cov <- NULL


  # add sample.nobs/group.label to lavaan call
  dotdotdot$sample.nobs <- sample.nobs

  # always use meanstructure = TRUE
  dotdotdot$meanstructure <- TRUE


  # remove 'ordered' argument: we will first pretend we generate
  # continuous data only
  dotdotdot$ordered <- NULL

  # 'fit' population model
  fit.pop <- do.call(cmd.pop, args = c(list(model = model), dotdotdot))

  # categorical?
  if (fit.pop@Model@categorical) {
    # refit, as if continuous only
    dotdotdot$ordered <- NULL
    fit.con <- do.call(cmd.pop, args = c(list(model = model), dotdotdot))
    # restore
    dotdotdot$ordered <- dotdotdot.orig$ordered
  } else {
    fit.con <- fit.pop
  }

  # extract model implied statistics and data slot
  lavimplied <- fit.con@implied # take continuous mean/cov
  lavdata <- fit.pop@Data
  lavmodel <- fit.pop@Model
  lavpartable <- fit.pop@ParTable
  lavoptions <- fit.pop@Options

  # number of groups/levels
  ngroups <- lav_partable_ngroups(lavpartable)
  nblocks <- lav_partable_nblocks(lavpartable)

  # check sample.nobs argument
  if (lavdata@nlevels > 1L) {
    # multilevel
    if (is.null(cluster.idx)) {
      # default? -> 1000 per block
      if (is.null(sample.nobs)) {
        sample.nobs <- rep.int(
          c(
            1000L,
            rep.int(100L, lavdata@nlevels - 1L)
          ),
          times = ngroups
        )
      } else {
        # we assume sample.nobs only contains a single number
        sample.nobs <- rep.int(
          c(
            sample.nobs,
            rep.int(100L, lavdata@nlevels - 1L)
          ),
          times = ngroups
        )
      }
    } else {
      # we got a cluster.idx argument
      if (!is.list(cluster.idx)) {
        cluster.idx <- rep(list(cluster.idx), ngroups)
      }

      if (!is.null(sample.nobs) && (length(sample.nobs) > 1L ||
        sample.nobs != 1000L)) {
        lav_msg_warn(gettext(
          "sample.nobs will be ignored if cluster.idx is provided"))
      }
      sample.nobs <- numeric(nblocks)
      for (g in seq_len(ngroups)) {
        gg <- (g - 1) * lavdata@nlevels + 1L
        sample.nobs[gg] <- length(cluster.idx[[g]])
        sample.nobs[gg + 1] <- length(unique(cluster.idx[[g]]))
      }
    }
  } else {
    # single level
    if (length(sample.nobs) == ngroups) {
      # nothing to do
    } else if (ngroups > 1L && length(sample.nobs) == 1L) {
      sample.nobs <- rep.int(sample.nobs, ngroups)
    } else {
      lav_msg_stop(gettextf(
        "ngroups = %1$s but sample.nobs has length = %2$s",
        ngroups, length(sample.nobs)))
    }
  }

  # check if ov.names are the same for each group
  if (ngroups > 1L) {
    N1 <- lavdata@ov.names[[1]]
    if (!all(sapply(
      lavdata@ov.names,
      function(x) all(x %in% N1)
    ))) {
      if (output == "data.frame") {
        output <- "matrix"
        lav_msg_warn(gettext(
          "groups do not contain the same set of variables;
          changing output= argument to \"matrix\""))
      }
    }
  }

  # prepare data containers
  X <- vector("list", length = nblocks)

  # generate data per BLOCK
  for (b in seq_len(nblocks)) {
    if (lavoptions$conditional.x) {
      lav_msg_stop(gettext("conditional.x is not ready yet"))
    } else {
      COV <- lavimplied$cov[[b]]
      MU <- lavimplied$mean[[b]]
    }

    # if empirical = TRUE, rescale by N/(N-1), so that estimator=ML
    # returns exact results
    if (empirical) {
      # check if sample.nobs is large enough
      if (sample.nobs[b] < NCOL(COV)) {
        lav_msg_stop(gettextf(
          "empirical = TRUE requires sample.nobs = %1$s to be larger than the
          number of variables = %2$s in block = %3$s",
          sample.nobs[b], NCOL(COV), b
        ))
      }
      if (lavdata@nlevels > 1L && (b %% lavdata@nlevels == 1L)) {
        COV <- COV * sample.nobs[b] / (sample.nobs[b] - sample.nobs[b + 1])
      } else {
        COV <- COV * sample.nobs[b] / (sample.nobs[b] - 1)
      }
    }

    # generate normal data
    tmp <- try(
      MASS::mvrnorm(
        n = sample.nobs[b],
        mu = MU, Sigma = COV, empirical = empirical
      ),
      silent = TRUE
    )

    if (inherits(tmp, "try-error")) {
      # something went wrong; most likely: non-positive COV?
      ev <- eigen(COV, symmetric = TRUE, only.values = TRUE)$values
      if (any(ev < 0)) {
        lav_msg_stop(gettextf(
          "model-implied covariance matrix is not positive-definite in block
          = %1$s; smallest eigen value = %2$s; change the model parameters.",
          b, round(min(ev), 5)))
      } else {
        lav_msg_stop(gettextf("data generation failed for block = %s", b))
      }
    } else {
      X[[b]] <- unname(tmp)
    }
  } # block

  if (output == "block") {
    return(X)
  }

  # if multilevel, make a copy, and create X[[g]] per group
  if (lavdata@nlevels > 1L) {
    X.block <- X
    X <- vector("list", length = ngroups)
  }

  # assemble data per group
  group.values <- lav_partable_group_values(lavpartable)
  for (g in 1:ngroups) {
    # multilevel?
    if (lavdata@nlevels > 1L) {
      # which block?
      bb <- (g - 1) * lavdata@nlevels + 1L

      Lp <- lavdata@Lp[[g]]
      p.tilde <- length(lavdata@ov.names[[g]])
      tmp1 <- matrix(0, nrow(X.block[[bb]]), p.tilde + 1L) # one extra for
      tmp2 <- matrix(0, nrow(X.block[[bb]]), p.tilde + 1L) # the clus id

      # level 1
      # if(empirical) {
      if (FALSE) {
        # force the within-cluster means to be zero (for both.idx vars)
        Y2 <- unname(as.matrix(aggregate(X.block[[bb]],
          # NOTE: cluster.idx becomes a factor
          # should be 111122223333...
          by = list(cluster.idx[[g]]), FUN = mean, na.rm = TRUE
        )[, -1]))
        # don't touch within-only variables
        w.idx <- match(Lp$within.idx[[2]], Lp$ov.idx[[1]])
        Y2[, w.idx] <- 0

        # center cluster-wise
        Y1c <- X.block[[bb]] - Y2[cluster.idx[[g]], , drop = FALSE]

        # this destroys the within covariance matrix
        sigma.sqrt <- lav_matrix_symmetric_sqrt(lavimplied$cov[[bb]])
        NY <- NROW(Y1c)
        S <- cov(Y1c) * (NY - 1) / NY
        S.inv <- solve(S)
        S.inv.sqrt <- lav_matrix_symmetric_sqrt(S.inv)

        # transform
        X.block[[bb]] <- Y1c %*% S.inv.sqrt %*% sigma.sqrt
      }
      tmp1[, Lp$ov.idx[[1]]] <- X.block[[bb]]

      # level 2
      tmp2[, Lp$ov.idx[[2]]] <- X.block[[bb + 1L]][cluster.idx[[g]], ,
        drop = FALSE
      ]
      # final
      X[[g]] <- tmp1 + tmp2

      # cluster id
      X[[g]][, p.tilde + 1L] <- cluster.idx[[g]]
    }

    # add variable names?
    if (add.labels) {
      if (lavdata@nlevels > 1L) {
        colnames(X[[g]]) <- c(lavdata@ov.names[[g]], "cluster")
      } else {
        colnames(X[[g]]) <- lavdata@ov.names[[g]]
      }
    }

    # any categorical variables?
    ov.ord <- lavNames(fit.pop, "ov.ord", group = group.values[g])
    if (is.list(ov.ord)) {
      # multilvel -> use within level only
      ov.ord <- ov.ord[[1L]]
    }
    if (length(ov.ord) > 0L) {
      ov.names <- lavdata@ov.names[[g]]

      # which block?
      bb <- (g - 1) * lavdata@nlevels + 1L

      # th/names
      TH.VAL <- as.numeric(fit.pop@implied$th[[bb]])
      if (length(lavmodel@num.idx[[bb]]) > 0L) {
        NUM.idx <- which(lavmodel@th.idx[[bb]] == 0)
        TH.VAL <- TH.VAL[-NUM.idx]
      }
      th.names <- fit.pop@pta$vnames$th[[bb]]
      TH.NAMES <- sapply(strsplit(th.names,
        split = "|",
        fixed = TRUE
      ), "[[", 1L)

      # use thresholds to cut
      for (o in ov.ord) {
        o.idx <- which(o == ov.names)
        th.idx <- which(o == TH.NAMES)
        th.val <- c(-Inf, sort(TH.VAL[th.idx]), +Inf)
        # center (because model-implied 'mean' may be nonzero)
        tmp <- X[[g]][, o.idx]
        tmp <- tmp - mean(tmp, na.rm = TRUE)
        X[[g]][, o.idx] <- cut(tmp, th.val, labels = FALSE)
      }
    }
  }


  # output
  if (output == "matrix") {
    if (ngroups == 1L) {
      out <- X[[1L]]
    } else {
      out <- X
    }
  } else if (output == "data.frame") {
    if (ngroups == 1L) {
      # convert to data.frame
      out <- as.data.frame(X[[1L]], stringsAsFactors = FALSE)
    } else if (ngroups > 1L) {
      # rbind
      out <- do.call("rbind", X)

      # add group column
      group <- rep.int(1:ngroups, times = sapply(X, NROW))
      out <- cbind(out, group)

      # convert to data.frame
      out <- as.data.frame(out, stringsAsFactors = FALSE)
    }
  } else if (output == "cov") {
    if (ngroups == 1L) {
      out <- cov(X[[1L]])
    } else {
      out <- lapply(X, cov)
    }
  } else {
    lav_msg_stop(gettextf("unknown option for argument output: %s", output))
  }

  if (return.fit) {
    attr(out, "fit") <- fit.pop
  }

  out
}
