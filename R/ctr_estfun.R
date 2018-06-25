# contributed by Ed Merkle (17 Jan 2013)


# YR 12 Feb 2013: small changes to match the results of lav_model_gradient
#                 in the multiple group case
# YR 30 May 2014: handle 1-variable case (fixing apply in lines 56, 62, 108)
# YR 05 Nov 2015: add remove.duplicated = TRUE, to cope with strucchange in
#                 case of simple equality constraints
# YR 19 Nov 2015: if constraints have been used, compute case-wise Lagrange
#                 multipliers, and define the scores as: SC + (t(R) lambda)
# YR 05 Feb 2016: catch conditional.x = TRUE: no support (for now), until
#                 we can use the generic 0.6 infrastructure for scores,
#                 including the missing-values case
# YR 16 Feb 2016: adapt to changed @Mp slot elements; add remove.empty.cases=
#                 argument

estfun.lavaan <- lavScores <- function(object, scaling = FALSE,
                                       ignore.constraints = FALSE,
                                       remove.duplicated = TRUE,
                                       remove.empty.cases = TRUE) {

  stopifnot(inherits(object, "lavaan"))

  # what if estimator != ML?
  # avoid hard error (using stop); throw a warning, and return an empty matrix
  if(object@Options$estimator != "ML") {
      warning("lavaan WARNING: scores only availalbe if estimator is ML")
      return(matrix(0,0,0))
  }

  # check if conditional.x = TRUE
  if(object@Model@conditional.x) {
      stop("lavaan ERROR: scores not available (yet) if conditional.x = TRUE")
  }

  # shortcuts
  lavdata        <- object@Data
  lavmodel       <- object@Model
  lavsamplestats <- object@SampleStats
  lavoptions     <- object@Options

  ## number variables/sample size
  #ntab <- unlist(lavsamplestats@nobs)
  ## change in 0.5-17: we keep the 'empty cases'
  ##                   and 'fill' in the scores at their 'case.idx'
  ##                   later, we remove the 'empty rows'
  #ntot <- max( object@Data@case.idx[[ object@Data@ngroups ]] )
  ntab <- unlist(lavdata@norig)
  ntot <- sum(ntab)

  npar <- lav_object_inspect_npar(object)
  #if(object@Model@eq.constraints) {
  #   npar <- NCOL(object@Model@eq.constraints.K)
  #}
  Score.mat <- matrix(NA, ntot, npar)

  for(g in 1:lavsamplestats@ngroups) {
    if (lavsamplestats@ngroups > 1){
      moments <- fitted(object)[[g]]
    } else {
      moments <- fitted(object)
    }
    Sigma.hat <- moments$cov

    if(lavoptions$likelihood == "wishart") {
        N1 <- lavsamplestats@nobs[[g]]/(lavsamplestats@nobs[[g]] - 1)
    } else {
        N1 <- 1
    }

    if(!lavsamplestats@missing.flag) { # complete data
      #if(lavmodel@meanstructure) { # mean structure
        nvar <- ncol(lavsamplestats@cov[[g]])
        Mu.hat <- moments$mean
        X <- lavdata@X[[g]]
        Sigma.inv <- inv.chol(Sigma.hat, logdet=FALSE)
        group.w <- (unlist(lavsamplestats@nobs)/lavsamplestats@ntotal)

        J <- matrix(1, 1L, ntab[g]) ## FIXME: needed? better maybe rowSums/colSums?
        J2 <- matrix(1, nvar, nvar)
        diag(J2) <- 0.5

        if(lavmodel@meanstructure) {
            ## scores.H1 (H1 = saturated model)
            mean.diff <- t(t(X) - Mu.hat %*% J)
            dx.Mu <- -1 * mean.diff %*% Sigma.inv
            dx.Sigma <- t(matrix(apply(mean.diff, 1L,
               function(x) lav_matrix_vech(- J2 * (Sigma.inv %*% (tcrossprod(x)*N1 - Sigma.hat) %*% Sigma.inv))), ncol=nrow(mean.diff)))

            scores.H1 <- cbind(dx.Mu, dx.Sigma)
        } else {
            mean.diff <- t(t(X) - lavsamplestats@mean[[g]] %*% J)
            dx.Sigma <- t(matrix(apply(mean.diff, 1L,
               function(x) lav_matrix_vech(- J2 * (Sigma.inv %*% (tcrossprod(x)*N1 - Sigma.hat) %*% Sigma.inv))), ncol=nrow(mean.diff)))
            scores.H1 <- dx.Sigma
        }
        ## FIXME? Seems like we would need group.w even in the
        ##        complete-data case:
        ##if(scaling){
        ##  scores.H1 <- group.w[g] * scores.H1
        ##}

      #} else {
      #  ## no mean structure
      #  stop("Score calculation with no mean structure is not implemented.")
      #}
    } else { # incomplete data
      nsub <- ntab[g]
      M <- lavsamplestats@missing[[g]]
      Mp <- lavdata@Mp[[g]]
      #pat.idx <- match(MP1$id, MP1$order)
      group.w <- (unlist(lavsamplestats@nobs)/lavsamplestats@ntotal)

      Mu.hat <- moments$mean
      nvar <- ncol(lavsamplestats@cov[[g]])
      score.sigma   <- matrix(0, nsub, nvar*(nvar+1)/2)
      score.mu <- matrix(0, nsub, nvar)

      for(p in 1:length(M)) {
        ## Data
        #X <- M[[p]][["X"]]
        case.idx <- Mp$case.idx[[p]]
        var.idx <- M[[p]][["var.idx"]]
        X <- lavdata@X[[g]][case.idx,var.idx,drop = FALSE]
        nobs <- M[[p]][["freq"]]
        ## Which unique entries of covariance matrix are estimated?
        ## (Used to keep track of scores in score.sigma)
        var.idx.mat <- tcrossprod(var.idx)
        Sigma.idx <- which(var.idx.mat[lower.tri(var.idx.mat, diag=T)]==1)

        J <- matrix(1, 1L, nobs) #[var.idx]
        J2 <- matrix(1, nvar, nvar)[var.idx, var.idx, drop = FALSE]
        diag(J2) <- 0.5
        Sigma.inv <- inv.chol(Sigma.hat[var.idx, var.idx, drop = FALSE],
                                       logdet=FALSE)
        Mu <- Mu.hat[var.idx]
        mean.diff <- t(t(X) - Mu %*% J)

        ## Scores for missing pattern p within group g
        score.mu[case.idx,var.idx] <- -1 * mean.diff %*% Sigma.inv
        score.sigma[case.idx,Sigma.idx] <- t(matrix(apply(mean.diff, 1L,
          function(x) lav_matrix_vech(- J2 * (Sigma.inv %*% (tcrossprod(x) - Sigma.hat[var.idx,var.idx,drop = FALSE]) %*% Sigma.inv)) ), ncol=nrow(mean.diff)) )

      }

      scores.H1 <- cbind(score.mu, score.sigma)
      if(scaling){
        scores.H1 <- group.w[g] * scores.H1
      }
    } # missing

    Delta <- computeDelta(lavmodel = lavmodel)[[g]]
    #if(lavmodel@eq.constraints) {
    #    Delta <- Delta %*% lavmodel@eq.constraints.K # + lavmodel@eq.constraints.k0
    #    #x <- as.numeric(lavmodel@eq.constraints.K %*% x) +
    #    #                lavmodel@eq.constraints.k0
    #}
    wi <- lavdata@case.idx[[g]]
    Score.mat[wi,] <- -scores.H1 %*% Delta
    if(scaling){
      Score.mat[wi,] <- (-1/ntot) * Score.mat[wi,]
    }

  } # g

  # handle empty rows
  if(remove.empty.cases) {
      #empty.idx <- which( apply(Score.mat, 1L,
      #                        function(x) sum(is.na(x))) == ncol(Score.mat) )
      empty.idx <- unlist(lapply(lavdata@Mp, "[[", "empty.idx"))
      if(length(empty.idx) > 0L) {
          Score.mat <- Score.mat[-empty.idx,,drop=FALSE]
      }
  }

  # provide column names
  colnames(Score.mat) <- names(lav_object_inspect_coef(object,
                               type = "free", add.labels = TRUE))

  # handle general constraints, so that the sum of the columns equals zero
  if(!ignore.constraints &&
     sum(lavmodel@ceq.linear.idx, lavmodel@ceq.nonlinear.idx,
         lavmodel@cin.linear.idx, lavmodel@cin.nonlinear.idx) > 0) {

      R <- object@Model@con.jac[,]
      PRE <- lav_constraints_lambda_pre(object)
      #LAMBDA <- -1 * t(PRE %*% t(Score.mat))
      #RLAMBDA <- t(t(R) %*% t(LAMBDA))
      Score.mat <- Score.mat - t( t(R) %*% PRE %*% t(Score.mat) )

  }

  # handle simple equality constraints
  if(remove.duplicated && lavmodel@eq.constraints) {
      simple.flag <- lav_constraints_check_simple(lavmodel)
      if(simple.flag) {
          K <- lav_constraints_R2K(lavmodel)
          Score.mat <- Score.mat %*% K
      } else {
          warning("lavaan WARNING: remove.duplicated is TRUE, but equality constraints do not appear to be simple; returning full scores")
      }
  }

  Score.mat
}

