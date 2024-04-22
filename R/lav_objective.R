# fitting function for standard ML
estimator.ML <- function(Sigma.hat = NULL, Mu.hat = NULL,
                         data.cov = NULL, data.mean = NULL,
                         data.cov.log.det = NULL,
                         meanstructure = FALSE) {
  # FIXME: WHAT IS THE BEST THING TO DO HERE??
  # CURRENTLY: return Inf  (at least for nlminb, this works well)
  if (!attr(Sigma.hat, "po")) {
    return(Inf)
  }


  Sigma.hat.inv <- attr(Sigma.hat, "inv")
  Sigma.hat.log.det <- attr(Sigma.hat, "log.det")
  nvar <- ncol(Sigma.hat)

  if (!meanstructure) {
    fx <- (Sigma.hat.log.det + sum(data.cov * Sigma.hat.inv) -
      data.cov.log.det - nvar)
  } else {
    W.tilde <- data.cov + tcrossprod(data.mean - Mu.hat)
    fx <- (Sigma.hat.log.det + sum(W.tilde * Sigma.hat.inv) -
      data.cov.log.det - nvar)
  }

  # no negative values
  if (is.finite(fx) && fx < 0.0) fx <- 0.0

  fx
}

# fitting function for standard ML
estimator.ML_res <- function(Sigma.hat = NULL, Mu.hat = NULL, PI = NULL,
                             res.cov = NULL, res.int = NULL, res.slopes = NULL,
                             res.cov.log.det = NULL,
                             cov.x = NULL, mean.x = NULL) {
  if (!attr(Sigma.hat, "po")) {
    return(Inf)
  }

  # augmented mean.x + cov.x matrix
  C3 <- rbind(
    c(1, mean.x),
    cbind(mean.x, cov.x + tcrossprod(mean.x))
  )

  Sigma.hat.inv <- attr(Sigma.hat, "inv")
  Sigma.hat.log.det <- attr(Sigma.hat, "log.det")
  nvar <- ncol(Sigma.hat)

  # sigma
  objective.sigma <- (Sigma.hat.log.det + sum(res.cov * Sigma.hat.inv) -
    res.cov.log.det - nvar)
  # beta
  OBS <- t(cbind(res.int, res.slopes))
  EST <- t(cbind(Mu.hat, PI))
  Diff <- OBS - EST
  objective.beta <- sum(Sigma.hat.inv * crossprod(Diff, C3) %*% Diff)

  fx <- objective.sigma + objective.beta

  # no negative values
  if (is.finite(fx) && fx < 0.0) fx <- 0.0

  fx
}


# fitting function for restricted ML
estimator.REML <- function(Sigma.hat = NULL, Mu.hat = NULL,
                           data.cov = NULL, data.mean = NULL,
                           data.cov.log.det = NULL,
                           meanstructure = FALSE,
                           group = 1L, lavmodel = NULL,
                           lavsamplestats = NULL, lavdata = NULL) {
  if (!attr(Sigma.hat, "po")) {
    return(Inf)
  }

  Sigma.hat.inv <- attr(Sigma.hat, "inv")
  Sigma.hat.log.det <- attr(Sigma.hat, "log.det")
  nvar <- ncol(Sigma.hat)

  if (!meanstructure) {
    fx <- (Sigma.hat.log.det + sum(data.cov * Sigma.hat.inv) -
      data.cov.log.det - nvar)
  } else {
    W.tilde <- data.cov + tcrossprod(data.mean - Mu.hat)
    fx <- (Sigma.hat.log.det + sum(W.tilde * Sigma.hat.inv) -
      data.cov.log.det - nvar)
  }

  lambda.idx <- which(names(lavmodel@GLIST) == "lambda")
  LAMBDA <- lavmodel@GLIST[[lambda.idx[group]]]
  data.cov.inv <- lavsamplestats@icov[[group]]
  reml.h0 <- log(det(t(LAMBDA) %*% Sigma.hat.inv %*% LAMBDA))
  reml.h1 <- log(det(t(LAMBDA) %*% data.cov.inv %*% LAMBDA))
  nobs <- lavsamplestats@nobs[[group]]

  # fx <- (Sigma.hat.log.det + tmp - data.cov.log.det - nvar) + 1/Ng * (reml.h0  - reml.h1)
  fx <- fx + (1 / nobs * (reml.h0 - reml.h1))

  # no negative values
  if (is.finite(fx) && fx < 0.0) fx <- 0.0

  fx
}

# 'classic' fitting function for GLS
# used again since 0.6-10 (we used the much slower estimator.WLS before)
estimator.GLS <- function(Sigma.hat = NULL, Mu.hat = NULL,
                          data.cov = NULL, data.cov.inv = NULL, data.mean = NULL,
                          meanstructure = FALSE) {
  tmp <- data.cov.inv %*% (data.cov - Sigma.hat)
  # tmp is not perfectly symmetric, so we use t(tmp) on the next line
  # to obtain the same value as estimator.WLS
  fx <- 0.5 * sum(tmp * t(tmp))

  if (meanstructure) {
    tmp2 <- sum(data.cov.inv * tcrossprod(data.mean - Mu.hat))
    fx <- fx + tmp2
  }

  # no negative values
  if (is.finite(fx) && fx < 0.0) fx <- 0.0

  fx
}

# general WLS estimator (Muthen, Appendix 4, eq 99 single group)
# full weight (WLS.V) matrix
estimator.WLS <- function(WLS.est = NULL, WLS.obs = NULL, WLS.V = NULL) {
  # diff <- as.matrix(WLS.obs - WLS.est)
  # fx <- as.numeric( t(diff) %*% WLS.V %*% diff )

  # since 0.5-17, we use crossprod twice
  diff <- WLS.obs - WLS.est
  fx <- as.numeric(crossprod(crossprod(WLS.V, diff), diff))
  # todo alternative: using chol(WLS.V)

  # no negative values
  if (is.finite(fx) && fx < 0.0) fx <- 0.0

  fx
}

# diagonally weighted LS (DWLS)
estimator.DWLS <- function(WLS.est = NULL, WLS.obs = NULL, WLS.VD = NULL) {
  diff <- WLS.obs - WLS.est
  fx <- sum(diff * diff * WLS.VD)

  # no negative values
  if (is.finite(fx) && fx < 0.0) fx <- 0.0

  fx
}

# Full Information ML estimator (FIML) handling the missing values
estimator.FIML <- function(Sigma.hat = NULL, Mu.hat = NULL, Yp = NULL,
                           h1 = NULL, N = NULL) {
  if (is.null(N)) {
    N <- sum(sapply(Yp, "[[", "freq"))
  }

  # Note: we ignore x.idx (if any)
  fx <- lav_mvnorm_missing_loglik_samplestats(
    Yp = Yp,
    Mu = Mu.hat, Sigma = Sigma.hat,
    log2pi = FALSE,
    minus.two = TRUE
  ) / N

  # ajust for h1
  if (!is.null(h1)) {
    fx <- fx - h1

    # no negative values
    if (is.finite(fx) && fx < 0.0) fx <- 0.0
  }

  fx
}

# pairwise maximum likelihood
# this is adapted from code written by Myrsini Katsikatsou
#
# some changes:
# - no distinction between x/y (ksi/eta)
# - 29/03/2016: adapt for exogenous covariates
# - 21/09/2016: added code for missing = doubly.robust (contributed by
#   Myrsini Katsikatsou)
# - HJ 18/10/2023: For sampling weights the lavcache$bifreq are weighted
estimator.PML <- function(Sigma.hat = NULL, # model-based var/cov/cor
                          Mu.hat = NULL, # model-based means
                          TH = NULL, # model-based thresholds + means
                          PI = NULL, # slopes
                          th.idx = NULL, # threshold idx per variable
                          num.idx = NULL, # which variables are numeric
                          X = NULL, # raw data
                          eXo = NULL, # eXo data
                          wt = NULL, # case weights
                          lavcache = NULL, # housekeeping stuff
                          missing = NULL) { # how to deal with missings?

  # YR 3 okt 2012
  # - the idea is to compute for each pair of variables, the model-based
  #   probability (or likelihood in mixed case) (that we observe the data
  #   for this pair under the model)
  # - if we have exogenous variables + conditional.x, do this for each case
  # - after taking logs, the sum over the cases gives the
  #   log probablity/likelihood for this pair
  # - the sum over all pairs gives the final PL based logl

  # first of all: check if all correlations are within [-1,1]
  # if not, return Inf; (at least with nlminb, this works well)

  # diagonal of Sigma.hat is not necessarily 1, even for categorical vars
  Sigma.hat2 <- Sigma.hat
  if (length(num.idx) > 0L) {
    diag(Sigma.hat2)[-num.idx] <- 1
  } else {
    diag(Sigma.hat2) <- 1
  }
  # all positive variances? (for continuous variables)
  if (any(diag(Sigma.hat2) < 0)) {
    OUT <- +Inf
    attr(OUT, "logl") <- as.numeric(NA)
    return(OUT)
  }
  Cor.hat <- cov2cor(Sigma.hat2) # to get correlations (rho!)
  cors <- lav_matrix_vech(Cor.hat, diagonal = FALSE)

  if (length(cors) > 0L && (any(abs(cors) > 1) ||
    any(is.na(cors)))) {
    # question: what is the best approach here??
    OUT <- +Inf
    attr(OUT, "logl") <- as.numeric(NA)
    return(OUT)
  }

  nvar <- nrow(Sigma.hat)
  if (is.null(eXo)) {
    nexo <- 0L
  } else {
    nexo <- NCOL(eXo)
  }
  pstar <- nvar * (nvar - 1) / 2
  ov.types <- rep("ordered", nvar)
  if (length(num.idx) > 0L) {
    ov.types[num.idx] <- "numeric"
  }


  ##### Three cases:
  ##### 1) all ordered, no exogenous (fast!)
  ##### 2) mixed ordered + continuous, no exogenous
  ##### 3) mixed ordered + continuous, exogenous (conditional.x = TRUE)




  ##### Case 1:
  #####  all ordered
  #####  no exogenous covariates
  #####
  if (all(ov.types == "ordered") && nexo == 0L) {
    # prepare for Myrsini's vectorization scheme
    long2 <- LongVecTH.Rho(
      no.x = nvar,
      all.thres = TH,
      index.var.of.thres = th.idx,
      rho.xixj = cors
    )
    # get expected probability per table, per pair
    pairwisePI <- pairwiseExpProbVec(
      ind.vec = lavcache$long,
      th.rho.vec = long2
    )
    pairwisePI_orig <- pairwisePI # for doubly.robust

    # get frequency per table, per pair
    logl <- sum(lavcache$bifreq * log(pairwisePI))

    # >>>>>>>> HJ/MK PML CODE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # FYI the bifreq are already weighted so this will work. Alternatively:
    if (!is.null(wt)) {
      logl <- sum(lavcache$sum_obs_weights_xixj_ab_vec * log(pairwisePI))
    }

    # more convenient fit function
    prop <- lavcache$bifreq / lavcache$nobs
    freq <- lavcache$bifreq
    if (!is.null(wt)) {
      prop <- lavcache$sum_obs_weights_xixj_ab_vec / sum(wt)
      freq <- lavcache$sum_obs_weights_xixj_ab_vec
    }

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # remove zero props # FIXME!!! or add 0.5???
    # zero.idx <- which(prop == 0.0)
    zero.idx <- which((prop == 0.0) | !is.finite(prop))
    if (length(zero.idx) > 0L) {
      freq <- freq[-zero.idx]
      prop <- prop[-zero.idx]
      pairwisePI <- pairwisePI[-zero.idx]
    }
    ## Fmin <- sum( prop*log(prop/pairwisePI) )
    Fmin <- sum(freq * log(prop / pairwisePI)) # to avoid 'N'

    if (missing == "available.cases" || missing == "doubly.robust") {
      uniPI <- univariateExpProbVec(TH = TH, th.idx = th.idx)

      # shortcuts
      unifreq <- lavcache$unifreq
      uninobs <- lavcache$uninobs
      uniweights <- lavcache$uniweights

      logl <- logl + sum(uniweights * log(uniPI))

      uniprop <- unifreq / uninobs

      # remove zero props
      # uni.zero.idx <- which(uniprop == 0.0)
      uni.zero.idx <- which((uniprop == 0.0) | !is.finite(uniprop))
      if (length(uni.zero.idx) > 0L) {
        uniprop <- uniprop[-uni.zero.idx]
        uniPI <- uniPI[-uni.zero.idx]
        uniweights <- uniweights[-uni.zero.idx]
      }

      Fmin <- Fmin + sum(uniweights * log(uniprop / uniPI))
    }

    if (missing == "doubly.robust") {
      # COMPUTE THE SUM OF THE EXPECTED BIVARIATE CONDITIONAL LIKELIHOODS
      # SUM_{i,j} [ E_{Yi,Yj|y^o}(lnf(Yi,Yj))) ]

      # First compute the terms of the summand. Since the cells of
      # pairwiseProbGivObs are zero for the pairs of variables that at least
      # one of the variables is observed (hence not contributing to the summand)
      # there is no need to construct an index vector for summing appropriately
      # within each individual.
      log_pairwisePI_orig <- log(pairwisePI_orig)
      pairwiseProbGivObs <- lavcache$pairwiseProbGivObs
      tmp_prod <- t(t(pairwiseProbGivObs) * log_pairwisePI_orig)

      SumElnfijCasewise <- apply(tmp_prod, 1, sum)
      SumElnfij <- sum(SumElnfijCasewise)
      logl <- logl + SumElnfij
      Fmin <- Fmin - SumElnfij

      # COMPUTE THE THE SUM OF THE EXPECTED UNIVARIATE CONDITIONAL LIKELIHOODS
      # SUM_{i,j} [ E_{Yj|y^o}(lnf(Yj|yi))) ]

      # First compute the model-implied conditional univariate probabilities
      # p(y_i=a|y_j=b). Let ModProbY1Gy2 be the vector of these
      # probabilities. The order the probabilities
      # are listed in the vector ModProbY1Gy2 is as follows:
      # y1|y2, y1|y3, ..., y1|yp, y2|y1, y2|y3, ..., y2|yp,
      # ..., yp|y1, yp|y2, ..., yp|y(p-1). Within each pair of variables the
      # index "a" which represents the response category of variable yi runs faster than
      # "b" which represents the response category of the given variable yj.
      # The computation of these probabilities are based on the model-implied
      # bivariate probabilities p(y_i=a,y_j=b). To do the appropriate summations
      # and divisions we need some index vectors to keep track of the index i, j,
      # a, and b, as well as the pair index. These index vectors should be
      # computed once and stored in lavcache. About where in the lavaan code
      # we will add the computations and how they will be done please see the
      # file "new objects in lavcache for DR-PL.r"

      idx.pairs <- lavcache$idx.pairs
      idx.cat.y2.split <- lavcache$idx.cat.y2.split
      idx.cat.y1.split <- lavcache$idx.cat.y1.split
      idx.Y1 <- lavcache$idx.Y1
      idx.Gy2 <- lavcache$idx.Gy2
      idx.cat.Y1 <- lavcache$idx.cat.Y1
      idx.cat.Gy2 <- lavcache$idx.cat.Gy2
      id.uniPrGivObs <- lavcache$id.uniPrGivObs
      # the latter keeps track which variable each column of the matrix
      # univariateProbGivObs refers to

      # For the function compute_uniCondProb_based_on_bivProb see the .r file
      # with the same name.
      ModProbY1Gy2 <- compute_uniCondProb_based_on_bivProb(
        bivProb = pairwisePI_orig,
        nvar = nvar,
        idx.pairs = idx.pairs,
        idx.Y1 = idx.Y1,
        idx.Gy2 = idx.Gy2,
        idx.cat.y1.split = idx.cat.y1.split,
        idx.cat.y2.split = idx.cat.y2.split
      )

      log_ModProbY1Gy2 <- log(ModProbY1Gy2)

      # Let univariateProbGivObs be the matrix of the conditional univariate
      # probabilities Pr(y_i=a|y^o) that has been computed in advance and are
      # fed to the DR-PL function. The rows represent different individuals,
      # i.e. nrow=nobs, and the columns different probabilities. The columns
      # are listed as follows: a runs faster than i.

      # Note that the number of columns of univariateProbGivObs is not the
      # same with the length(log_ModProbY1Gy2), actually
      # ncol(univariateProbGivObs) < length(log_ModProbY1Gy2).
      # For this we use the following commands in order to multiply correctly.

      # Compute for each case the product  Pr(y_i=a|y^o) * log[ p(y_i=a|y_j=b) ]
      # i.e. univariateProbGivObs * log_ModProbY1Gy2
      univariateProbGivObs <- lavcache$univariateProbGivObs
      nobs <- nrow(X)
      uniweights.casewise <- lavcache$uniweights.casewise
      id.cases.with.missing <- which(uniweights.casewise > 0)
      no.cases.with.missing <- length(id.cases.with.missing)
      no.obs.casewise <- nvar - uniweights.casewise
      idx.missing.var <- apply(X, 1, function(x) {
        which(is.na(x))
      })
      idx.observed.var <- lapply(idx.missing.var, function(x) {
        c(1:nvar)[-x]
      })
      idx.cat.observed.var <- sapply(1:nobs, function(i) {
        X[i, idx.observed.var[[i]]]
      })
      ElnyiGivyjbCasewise <- sapply(1:no.cases.with.missing, function(i) {
        tmp.id.case <- id.cases.with.missing[i]
        tmp.no.mis <- uniweights.casewise[tmp.id.case]
        tmp.idx.mis <- idx.missing.var[[tmp.id.case]]
        tmp.idx.obs <- idx.observed.var[[tmp.id.case]]
        tmp.no.obs <- no.obs.casewise[tmp.id.case]
        tmp.idx.cat.obs <- idx.cat.observed.var[[tmp.id.case]]
        tmp.uniProbGivObs.i <- univariateProbGivObs[tmp.id.case, ]
        sapply(1:tmp.no.mis, function(k) {
          tmp.idx.mis.var <- tmp.idx.mis[k]
          tmp.uniProbGivObs.ik <-
            tmp.uniProbGivObs.i[id.uniPrGivObs == tmp.idx.mis.var]
          tmp.log_ModProbY1Gy2 <- sapply(1:tmp.no.obs, function(z) {
            log_ModProbY1Gy2[idx.Y1 == tmp.idx.mis.var &
              idx.Gy2 == tmp.idx.obs[z] &
              idx.cat.Gy2 == tmp.idx.cat.obs[z]]
          })
          sum(tmp.log_ModProbY1Gy2 * tmp.uniProbGivObs.ik)
        })
      })
      ElnyiGivyjb <- sum(unlist(ElnyiGivyjbCasewise))
      logl <- logl + ElnyiGivyjb
      # for the Fmin function
      Fmin <- Fmin - ElnyiGivyjb
    } # end of if (missing =="doubly.robust")


    ##### Case 2:
    #####  mixed ordered + numeric
    #####  no exogenous covariates
    #####
  } else if (nexo == 0L) {
    # mixed ordered/numeric variables, but no exogenous covariates
    # - no need to compute 'casewise' (log)likelihoods

    PSTAR <- matrix(0, nvar, nvar) # utility matrix, to get indices
    PSTAR[lav_matrix_vech_idx(nvar, diagonal = FALSE)] <- 1:pstar
    N <- NROW(X)

    logLikPair <- numeric(pstar) # logl per pair (summed over cases)
    for (j in seq_len(nvar - 1L)) {
      for (i in (j + 1L):nvar) {
        pstar.idx <- PSTAR[i, j]
        if (ov.types[i] == "numeric" &&
          ov.types[j] == "numeric") {
          logLIK <- lav_mvnorm_loglik_data(
            Y = X[, c(i, j)], wt = wt, Mu = Mu.hat[c(i, j)],
            Sigma = Sigma.hat[c(i, j), c(i, j)], casewise = TRUE
          )
          logLikPair[pstar.idx] <- sum(logLIK, na.rm = TRUE)
        } else if (ov.types[i] == "numeric" &&
          ov.types[j] == "ordered") {
          # polyserial correlation
          logLIK <- lav_bvmix_lik(
            Y1 = X[, i], Y2 = X[, j],
            wt = wt,
            evar.y1 = Sigma.hat[i, i],
            beta.y1 = Mu.hat[i],
            th.y2 = TH[th.idx == j],
            rho = Cor.hat[i, j], .log = TRUE
          )
          logLikPair[pstar.idx] <- sum(logLIK, na.rm = TRUE)
        } else if (ov.types[j] == "numeric" &&
          ov.types[i] == "ordered") {
          # polyserial correlation
          logLIK <- lav_bvmix_lik(
            Y1 = X[, j], Y2 = X[, i],
            wt = wt,
            evar.y1 = Sigma.hat[j, j],
            beta.y1 = Mu.hat[j],
            th.y2 = TH[th.idx == i],
            rho = Cor.hat[i, j], .log = TRUE
          )
          logLikPair[pstar.idx] <- sum(logLIK, na.rm = TRUE)
        } else if (ov.types[i] == "ordered" &&
          ov.types[j] == "ordered") {
          # polychoric correlation
          pairwisePI <- lav_bvord_noexo_pi(
            rho = Cor.hat[i, j],
            th.y1 = TH[th.idx == i],
            th.y2 = TH[th.idx == j]
          )
          # avoid zeroes
          pairwisePI[pairwisePI < .Machine$double.eps] <-
            .Machine$double.eps
          # note: missing values are just not counted
          FREQ <- lav_bvord_freq(X[, i], X[, j], wt = wt)
          logLikPair[pstar.idx] <- sum(FREQ * log(pairwisePI))
        }
      }
    } # all pairs

    na.idx <- which(is.na(logLikPair))
    if (length(na.idx) > 0L) {
      lav_msg_warn(gettext("some pairs produces NA values for logl:"),
        lav_msg_view(round(logLikPair, 3), "none")
      )
    }

    # sum over pairs
    logl <- sum(logLikPair)

    # Fmin
    Fmin <- (-1) * logl


    ##### Case 3:
    #####  mixed ordered + numeric
    #####  exogenous covariates
    #####  (conditional.x = TRUE)
  } else {
    LIK <- matrix(0, nrow(X), pstar) # likelihood per case, per pair
    PSTAR <- matrix(0, nvar, nvar) # utility matrix, to get indices
    PSTAR[lav_matrix_vech_idx(nvar, diagonal = FALSE)] <- 1:pstar
    N <- NROW(X)

    for (j in seq_len(nvar - 1L)) {
      for (i in (j + 1L):nvar) {
        pstar.idx <- PSTAR[i, j]
        # cat("pstar.idx =", pstar.idx, "i = ", i, " j = ", j, "\n")
        if (ov.types[i] == "numeric" &&
          ov.types[j] == "numeric") {
          # ordinary pearson correlation
          LIK[, pstar.idx] <-
            lav_bvreg_lik(
              Y1 = X[, i], Y2 = X[, j], eXo = eXo,
              wt = wt,
              evar.y1 = Sigma.hat[i, i],
              beta.y1 = c(Mu.hat[i], PI[i, ]),
              evar.y2 = Sigma.hat[j, j],
              beta.y2 = c(Mu.hat[j], PI[j, ]),
              rho = Cor.hat[i, j]
            )
        } else if (ov.types[i] == "numeric" &&
          ov.types[j] == "ordered") {
          # polyserial correlation
          ### FIXME: th.y2 should go into ps_lik!!!
          LIK[, pstar.idx] <-
            lav_bvmix_lik(
              Y1 = X[, i], Y2 = X[, j], eXo = eXo,
              wt = wt,
              evar.y1 = Sigma.hat[i, i],
              beta.y1 = c(Mu.hat[i], PI[i, ]),
              th.y2 = TH[th.idx == j],
              sl.y2 = PI[j, ],
              rho = Cor.hat[i, j]
            )
        } else if (ov.types[j] == "numeric" &&
          ov.types[i] == "ordered") {
          # polyserial correlation
          ### FIXME: th.y1 should go into ps_lik!!!
          LIK[, pstar.idx] <-
            lav_bvmix_lik(
              Y1 = X[, j], Y2 = X[, i], eXo = eXo,
              wt = wt,
              evar.y1 = Sigma.hat[j, j],
              beta.y1 = c(Mu.hat[j], PI[j, ]),
              th.y2 = TH[th.idx == i],
              sl.y2 = PI[i, ],
              rho = Cor.hat[i, j]
            )
        } else if (ov.types[i] == "ordered" &&
          ov.types[j] == "ordered") {
          LIK[, pstar.idx] <-
            pc_lik_PL_with_cov(
              Y1 = X[, i],
              Y2 = X[, j],
              Rho = Sigma.hat[i, j],
              th.y1 = TH[th.idx == i],
              th.y2 = TH[th.idx == j],
              eXo = eXo,
              PI.y1 = PI[i, ],
              PI.y2 = PI[j, ],
              missing.ind = missing
            )
        }
      }
    } # all pairs

    # check for zero likelihoods/probabilities
    # FIXME: or should we replace them with a tiny number?
    if (any(LIK == 0.0, na.rm = TRUE)) {
      OUT <- +Inf
      attr(OUT, "logl") <- as.numeric(NA)
      return(OUT)
    }

    # loglikelihood
    LogLIK.cases <- log(LIK)

    # sum over cases
    LogLIK.pairs <- colSums(LogLIK.cases, na.rm = TRUE)

    # sum over pairs
    logl <- logl_pairs <- sum(LogLIK.pairs)

    if (missing == "available.cases" && all(ov.types == "ordered") &&
      nexo != 0L) {
      uni_LIK <- matrix(0, nrow(X), ncol(X))
      for (i in seq_len(nvar)) {
        uni_LIK[, i] <- uni_lik(
          Y1 = X[, i],
          th.y1 = TH[th.idx == i],
          eXo = eXo,
          PI.y1 = PI[i, ]
        )
      }

      if (any(uni_LIK == 0.0, na.rm = TRUE)) {
        OUT <- +Inf
        attr(OUT, "logl") <- as.numeric(NA)
        return(OUT)
      }

      uni_logLIK_cases <- log(uni_LIK) * lavcache$uniweights.casewise

      # sum over cases
      uni_logLIK_varwise <- colSums(uni_logLIK_cases)

      # sum over variables
      uni_logLIK <- sum(uni_logLIK_varwise)

      # add with the pairwise part of LogLik
      logl <- logl_pairs + uni_logLIK
    }

    # we minimise
    Fmin <- (-1) * logl
  }


  # here, we should have two quantities: logl and Fmin

  # function value as returned to the minimizer
  fx <- Fmin

  # attach 'loglikelihood'
  attr(fx, "logl") <- logl

  fx
}

# full information maximum likelihood
# underlying multivariate normal approach (see Joreskog & Moustaki, 2001)
#
estimator.FML <- function(Sigma.hat = NULL, # model-based var/cov/cor
                          TH = NULL, # model-based thresholds + means
                          th.idx = NULL, # threshold idx per variable
                          num.idx = NULL, # which variables are numeric
                          X = NULL, # raw data
                          lavcache = NULL) { # patterns

  # YR 27 aug 2013
  # just for fun, and to compare with PML for small models

  # first of all: check if all correlations are within [-1,1]
  # if not, return Inf; (at least with nlminb, this works well)
  cors <- Sigma.hat[lower.tri(Sigma.hat)]

  if (any(abs(cors) > 1)) {
    return(+Inf)
  }

  nvar <- nrow(Sigma.hat)
  pstar <- nvar * (nvar - 1) / 2
  ov.types <- rep("ordered", nvar)
  if (length(num.idx) > 0L) ov.types[num.idx] <- "numeric"
  MEAN <- rep(0, nvar)

  # shortcut for all ordered - per pattern
  if (all(ov.types == "ordered")) {
    PAT <- lavcache$pat
    npatterns <- nrow(PAT)
    freq <- as.numeric(rownames(PAT))
    PI <- numeric(npatterns)
    TH.VAR <- lapply(1:nvar, function(x) c(-Inf, TH[th.idx == x], +Inf))
    # FIXME!!! ok to set diagonal to 1.0?
    diag(Sigma.hat) <- 1.0
    for (r in 1:npatterns) {
      # compute probability for each pattern
      lower <- sapply(1:nvar, function(x) TH.VAR[[x]][PAT[r, x]])
      upper <- sapply(1:nvar, function(x) TH.VAR[[x]][PAT[r, x] + 1L])


      # how accurate must we be here???
      PI[r] <- sadmvn(lower, upper,
        mean = MEAN, varcov = Sigma.hat,
        maxpts = 10000 * nvar, abseps = 1e-07
      )
    }
    # sum (log)likelihood over all patterns
    # LogLik <- sum(log(PI) * freq)

    # more convenient fit function
    prop <- freq / sum(freq)
    # remove zero props # FIXME!!! or add 0.5???
    zero.idx <- which(prop == 0.0)
    if (length(zero.idx) > 0L) {
      prop <- prop[-zero.idx]
      PI <- PI[-zero.idx]
    }
    Fmin <- sum(prop * log(prop / PI))
  } else { # case-wise
    PI <- numeric(nobs)
    for (i in 1:nobs) {
      # compute probability for each case
      PI[i] <- lav_msg_stop(gettext("not implemented"))
    }
    # sum (log)likelihood over all observations
    LogLik <- sum(log(PI))
    lav_msg_stop(gettext("not implemented"))
  }

  # function value as returned to the minimizer
  # fx <- -1 * LogLik
  fx <- Fmin

  fx
}

estimator.MML <- function(lavmodel = NULL,
                          THETA = NULL,
                          TH = NULL,
                          GLIST = NULL,
                          group = 1L,
                          lavdata = NULL,
                          sample.mean = NULL,
                          sample.mean.x = NULL,
                          lavcache = NULL) {
  # compute case-wise likelihoods
  lik <- lav_model_lik_mml(
    lavmodel = lavmodel, THETA = THETA, TH = TH,
    GLIST = GLIST, group = group, lavdata = lavdata,
    sample.mean = sample.mean, sample.mean.x = sample.mean.x,
    lavcache = lavcache
  )

  # log + sum over observations
  logl <- sum(log(lik))

  # function value as returned to the minimizer
  fx <- -logl

  fx
}

estimator.2L <- function(lavmodel = NULL,
                         GLIST = NULL,
                         Y1 = NULL, # only for missing
                         Lp = NULL,
                         Mp = NULL,
                         lavsamplestats = NULL,
                         group = 1L) {
  # compute model-implied statistics for all blocks
  implied <- lav_model_implied(lavmodel, GLIST = GLIST)

  # here, we assume only 2!!! levels, at [[1]] and [[2]]
  if (lavmodel@conditional.x) {
    Res.Sigma.W <- implied$res.cov[[(group - 1) * 2 + 1]]
    Res.Int.W <- implied$res.int[[(group - 1) * 2 + 1]]
    Res.Pi.W <- implied$res.slopes[[(group - 1) * 2 + 1]]

    Res.Sigma.B <- implied$res.cov[[(group - 1) * 2 + 2]]
    Res.Int.B <- implied$res.int[[(group - 1) * 2 + 2]]
    Res.Pi.B <- implied$res.slopes[[(group - 1) * 2 + 2]]
  } else {
    Sigma.W <- implied$cov[[(group - 1) * 2 + 1]]
    Mu.W <- implied$mean[[(group - 1) * 2 + 1]]
    Sigma.B <- implied$cov[[(group - 1) * 2 + 2]]
    Mu.B <- implied$mean[[(group - 1) * 2 + 2]]
  }

  if (lavsamplestats@missing.flag) {
    if (lavmodel@conditional.x) {
      lav_msg_stop(gettext("multilevel + conditional.x is not ready yet for
                           fiml; rerun with conditional.x = FALSE"))
    }
    # SIGMA.B <- Sigma.B[Lp$both.idx[[2]], Lp$both.idx[[2]], drop = FALSE]
    # if(any(diag(SIGMA.B) < 0)) {
    #    return(+Inf)
    # }
    # COR.B <- cov2cor(SIGMA.B)
    # if(any(abs(lav_matrix_vech(COR.B, diagonal = FALSE)) > 1)) {
    #   return(+Inf)
    # }

    Y2 <- lavsamplestats@YLp[[group]][[2]]$Y2
    Yp <- lavsamplestats@missing[[group]]
    loglik <- lav_mvnorm_cluster_missing_loglik_samplestats_2l(
      Y1 = Y1,
      Y2 = Y2, Lp = Lp, Mp = Mp,
      Mu.W = Mu.W, Sigma.W = Sigma.W,
      Mu.B = Mu.B, Sigma.B = Sigma.B,
      log2pi = FALSE, minus.two = TRUE
    )
  } else {
    YLp <- lavsamplestats@YLp[[group]]
    if (lavmodel@conditional.x) {
      loglik <- lav_mvreg_cluster_loglik_samplestats_2l(
        YLp = YLp, Lp = Lp,
        Res.Sigma.W = Res.Sigma.W,
        Res.Int.W = Res.Int.W, Res.Pi.W = Res.Pi.W,
        Res.Sigma.B = Res.Sigma.B,
        Res.Int.B = Res.Int.B, Res.Pi.B = Res.Pi.B,
        log2pi = FALSE, minus.two = TRUE
      )
    } else {
      loglik <- lav_mvnorm_cluster_loglik_samplestats_2l(
        YLp = YLp, Lp = Lp,
        Mu.W = Mu.W, Sigma.W = Sigma.W,
        Mu.B = Mu.B, Sigma.B = Sigma.B,
        log2pi = FALSE, minus.two = TRUE
      )
    }
  }

  # minimize
  objective <- 1 * loglik

  # divide by (N*2)
  objective <- objective / (lavsamplestats@ntotal * 2)

  # should be strictly positive
  # if(objective < 0) {
  #   objective <- +Inf
  # }

  objective
}
