# contributed by Ed Merkle (17 Jan 2013)

# small changes by YR (12 Feb 2013) to match the results
# of lav_model_gradient in the multiple group case

# YR 30 May 2014: handle 1-variable case (fixing apply in lines 56, 62, 108)

estfun.lavaan <- lavScores <- function(object, scaling=FALSE) {

  stopifnot(inherits(object, "lavaan"))

  # what if estimator != ML? 
  # avoid hard error (using stop); throw a warning, and return an empty matrix
  if(object@Options$estimator != "ML") {
      warning("lavaan WARNING: scores only availalbe if estimator is ML")
      return(matrix(0,0,0))
  }


  # shortcuts
  lavdata        <- object@Data
  lavmodel       <- object@Model
  lavsamplestats <- object@SampleStats
  lavoptions     <- object@Options

  ## number variables/sample size
  ntab <- unlist(lavsamplestats@nobs)
  ## change in 0.5-17: we keep the 'empty cases'
  ntot <- ( lavsamplestats@ntotal + 
           length(sapply(object@Data@Mp, "[[", "empty.idx")) )

  Score.mat <- matrix(NA, ntot, length(coef(object)))
  
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
               function(x) vech(- J2 * (Sigma.inv %*% (tcrossprod(x)*N1 - Sigma.hat) %*% Sigma.inv))), ncol=nrow(mean.diff)))

            scores.H1 <- cbind(dx.Mu, dx.Sigma)
        } else {
            mean.diff <- t(t(X) - lavsamplestats@mean[[g]] %*% J)
            dx.Sigma <- t(matrix(apply(mean.diff, 1L,
               function(x) vech(- J2 * (Sigma.inv %*% (tcrossprod(x)*N1 - Sigma.hat) %*% Sigma.inv))), ncol=nrow(mean.diff)))
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
      MP1 <- lavdata@Mp[[g]]
      pat.idx <- match(MP1$id, MP1$order)
      group.w <- (unlist(lavsamplestats@nobs)/lavsamplestats@ntotal)

      Mu.hat <- moments$mean
      nvar <- ncol(lavsamplestats@cov[[g]])
      score.sigma   <- matrix(0, nsub, nvar*(nvar+1)/2)
      score.mu <- matrix(0, nsub, nvar)
    
      for(p in 1:length(M)) {
        ## Data
        X <- M[[p]][["X"]]
        nobs <- M[[p]][["nobs"]]
        var.idx <- M[[p]][["var.idx"]]
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
        score.mu[pat.idx==p,var.idx] <- -1 * mean.diff %*% Sigma.inv
        score.sigma[pat.idx==p,Sigma.idx] <- t(matrix(apply(mean.diff, 1L,
          function(x) vech(- J2 * (Sigma.inv %*% (tcrossprod(x) - Sigma.hat[var.idx,var.idx,drop = FALSE]) %*% Sigma.inv)) ), ncol=nrow(mean.diff)) )

      }

      scores.H1 <- cbind(score.mu, score.sigma)
      if(scaling){
        scores.H1 <- group.w[g] * scores.H1
      }
    } # missing
    
    Delta <- computeDelta(lavmodel = lavmodel)[[g]]
    wi <- lavdata@case.idx[[g]]
    Score.mat[wi,] <- -scores.H1 %*% Delta
    if(scaling){
      Score.mat[wi,] <- (-1/ntot) * Score.mat[wi,]
    }
  } # g
  
  # provide column names
  colnames(Score.mat) <- names(coef(object))

  Score.mat
}
