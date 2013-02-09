# contributed by Ed Merkle (17 Jan 2013)
estfun <- lavScores <- function(object)
{
  ## number variables/sample size
  samplestats <- object@SampleStats
  ntab <- unlist(samplestats@nobs)
  nobs <- samplestats@ntotal

  Score.mat <- matrix(NA, nobs, length(coef(object)))
  
  for(g in 1:samplestats@ngroups) {
    if (samplestats@ngroups > 1){
      moments <- fitted(object)[[g]]
    } else {
      moments <- fitted(object)
    }
    Sigma.hat <- moments$cov
  
    if(!samplestats@missing.flag) { # complete data
      if(object@Model@meanstructure) { # mean structure
        nvar <- ncol(samplestats@cov[[g]])
        Mu.hat <- moments$mean
        X <- object@Data@X[[g]]
        Sigma.inv <- lavaan:::inv.chol(Sigma.hat, logdet=FALSE)

        J <- matrix(1, 1L, ntab[g]) ## FIXME: needed? better maybe rowSums/colSums?
        J2 <- matrix(1, nvar[g], nvar[g])
        diag(J2) <- 0.5

        ## scores.H1 (H1 = saturated model)
        mean.diff <- t(t(X) - Mu.hat %*% J)

        dx.Mu <- -1 * mean.diff %*% Sigma.inv

        dx.Sigma <- t(apply(mean.diff, 1L,
           function(x) lavaan:::vech(- J2 * (Sigma.inv %*% (tcrossprod(x) - Sigma.hat) %*% Sigma.inv))))

        scores.H1 <- cbind(dx.Mu, dx.Sigma)
    
      } else {
        ## no mean structure
        stop("Score calculation with no mean structure is not implemented.")
      }
    } else { # incomplete data
      nsub <- ntab[g]
      M <- samplestats@missing[[g]]
      MP1 <- object@Data@Mp[[g]]
      pat.idx <- match(MP1$id, MP1$order)

      Mu.hat <- moments$mean
      nvar <- ncol(samplestats@cov[[g]])
      score.sigma   <- matrix(0, nsub, nvar*(nvar+1)/2)
      score.mu <- matrix(0, nsub, nvar)
    
      for(p in 1:length(M)) {
        ## Data
        X <- M[[p]][["X"]]
        ## Cov matrix of data in pattern p (divide by n, not n-1)
        SX <- M[[p]][["SX"]]
        ## Mean vector of data in pattern p
        MX <- M[[p]][["MX"]]
        nobs <- M[[p]][["nobs"]]
        var.idx <- M[[p]][["var.idx"]]
        ## Which unique entries of covariance matrix are estimated?
        ## (Used to keep track of scores in score.sigma)
        var.idx.mat <- tcrossprod(var.idx)
        Sigma.idx <- which(var.idx.mat[lower.tri(var.idx.mat, diag=T)]==1)
      
        J <- matrix(1, 1L, nobs) #[var.idx]
        J2 <- matrix(1, nvar, nvar)[var.idx, var.idx]
        diag(J2) <- 0.5
        Sigma.inv <- lavaan:::inv.chol(Sigma.hat[var.idx, var.idx],
                                       logdet=FALSE)
        Mu <- Mu.hat[var.idx]
        mean.diff <- t(t(X) - Mu %*% J)

        ## Scores for group g
        score.mu[pat.idx==p,var.idx] <- -1 * mean.diff %*% Sigma.inv
        score.sigma[pat.idx==p,Sigma.idx] <- t(apply(mean.diff, 1L,
          function(x) lavaan:::vech(- J2 * (Sigma.inv %*% (tcrossprod(x) - Sigma.hat[var.idx,var.idx]) %*% Sigma.inv)) ) )

        scores.H1 <- cbind(score.mu, score.sigma)

        ## FIXME?: On the above line in computeOmega, Yves had W.tilde
        ##        (see below) instead of
        ##        tcrossprod(x).  But I believe this is only needed
        ##        when we take a gradient across a sample, vs for an
        ##        individual observation.

        ##W.tilde <- SX + tcrossprod(MX - Mu)
        ##OMEGA[var.idx, var.idx] <-
        ##  ( OMEGA[var.idx, var.idx] + nobs/samplestats@ntotal *
        ##   (Sigma.inv %*%
        ##    (W.tilde - Sigma.hat[[g]][var.idx,var.idx]) %*%
        ##    Sigma.inv ) )
      }
    } # missing
    
    Delta <- lavaan:::computeDelta(object@Model)[[g]]
    wi <- object@Data@case.idx[[g]]
    Score.mat[wi,] <- -scores.H1 %*% Delta
  } # g

  Score.mat
}
