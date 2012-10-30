# pairwise maximum likelihood
# this is adapted from code written by Myrsini Katsikatsou
#
# the first derivative of the pairwise logLik function with respect to the
# thresholds/slopes/var/correlations; together with DELTA, we can use the
# chain rule to get the gradient
#
# first attempt - YR 5 okt 2012

pml_deriv1 <- function(Sigma.hat = NULL,    # model-based var/cov/cor
                       TH        = NULL,    # model-based thresholds + means
                       th.idx    = NULL,    # threshold idx per variable
                       num.idx   = NULL,    # which variables are numeric
                       X         = NULL,    # data
                       eXo       = NULL,    # external covariates
                       scores    = FALSE,   # return case-wise scores
                       negative  = TRUE) {  # multiply by -1

    nvar <- nrow(Sigma.hat)
    pstar <- nvar*(nvar-1)/2
    ov.types <- rep("ordered", nvar)
    if(length(num.idx) > 0L) ov.types[num.idx] <- "numeric"
    if(!is.null(eXo)) {
        nexo <- ncol(eXo)
    } else {
        nexo <- 0
    }
    N.TH  <- length(th.idx)
    N.SL  <- nvar * nexo
    N.VAR <- length(num.idx)
    N.COR <- pstar

    #print(Sigma.hat); print(TH); print(th.idx); print(num.idx); print(str(X))

    # in this order: TH/MEANS + SLOPES + VAR + COR
    GRAD.size <- N.TH + N.SL + N.VAR + N.COR

    # scores or gradient?
    if(scores) {
        SCORES <- matrix(0, nrow(X), GRAD.size) # we will sum up over all pairs 
    } else {
        GRAD <- matrix(0, pstar, GRAD.size) # each pair is a row
    }
    PSTAR <- matrix(0, nvar, nvar)   # utility matrix, to get indices
    PSTAR[lavaan:::vech.idx(nvar, diag=FALSE)] <- 1:pstar

    for(j in seq_len(nvar-1L)) {
        for(i in (j+1L):nvar) {
            # cat(" i = ", i, " j = ", j, "\n") # debug only
            pstar.idx <- PSTAR[i,j]
            cor.idx <- N.TH + N.SL + N.VAR + PSTAR[i,j]
            th.idx_i <- which(th.idx == i)
            th.idx_j <- which(th.idx == j)
            if(nexo > 0L) {
                sl.idx_i <- N.TH + seq(i, by=nvar, length.out=nexo)
                sl.idx_j <- N.TH + seq(j, by=nvar, length.out=nexo)

                var.idx_i <- N.TH + N.SL + match(i, num.idx)
                var.idx_j <- N.TH + N.SL + match(j, num.idx)
            } else {
                var.idx_i <- N.TH + match(i, num.idx)
                var.idx_j <- N.TH + match(j, num.idx)
            }
            if(ov.types[i] == "numeric" && ov.types[j] == "numeric") {
                # ordinary pearson correlation
                stop("not done yet")
            } else if(ov.types[i] == "numeric" && ov.types[j] == "ordered") {
                # polyserial correlation
                stop("not done yet")
            } else if(ov.types[j] == "numeric" && ov.types[i] == "ordered") {
                # polyserial correlation
                stop("not done yet")
            } else if(ov.types[i] == "ordered" && ov.types[j] == "ordered") {
                # polychoric correlation
                SC.COR.UNI <- pc_cor_scores(Y1  = X[,i],
                                            Y2  = X[,j],
                                            eXo = NULL,
                                            rho = Sigma.hat[i,j],
                                            fit.y1 = NULL, # fixme
                                            fit.y2 = NULL, # fixme
                                            th.y1 = TH[ th.idx == i ],
                                            th.y2 = TH[ th.idx == j ],
                                            sl.y1 = NULL,
                                            sl.y2 = NULL)
                
                if(scores) {
                    # TH
                    SCORES[,th.idx_i] <- SCORES[,th.idx_i] + SC.COR.UNI$dx.th.y1
                    SCORES[,th.idx_j] <- SCORES[,th.idx_j] + SC.COR.UNI$dx.th.y2

                    # SL
                    if(nexo > 0L) {
                        SCORES[,sl.idx_i] <- SCORES[,sl.idx_i] + SC.COR.UNI$dx.sl.y1
                        SCORES[,sl.idx_j] <- SCORES[,sl.idx_j] + SC.COR.UNI$dx.sl.y2 
                    }
                    # NO VAR
                    # RHO
                    SCORES[,cor.idx] <- SCORES[,cor.idx] + SC.COR.UNI$dx.rho
                } else {
                    # TH
                    if(length(th.idx_i) > 1L) {
                        GRAD[pstar.idx, th.idx_i] <- colSums(SC.COR.UNI$dx.th.y1)
                    } else {
                        GRAD[pstar.idx, th.idx_i] <- sum(SC.COR.UNI$dx.th.y1)
                    }
                    if(length(th.idx_j) > 1L) {
                         GRAD[pstar.idx, th.idx_j] <- colSums(SC.COR.UNI$dx.th.y2)
                    } else {
                         GRAD[pstar.idx, th.idx_j] <- sum(SC.COR.UNI$dx.th.y2)
                    }
    
                    # SL
                    if(nexo > 0L) {
                        if(length(sl.idx_i) > 1L) {
                            GRAD[pstar.idx, sl.idx_i] <- colSums(SC.COR.UNI$dx.sl.y1)
                        } else {
                            GRAD[pstar.idx, sl.idx_i] <- sum(SC.COR.UNI$dx.sl.y1)
                        }
                        if(length(sl.idx_j) > 1L) {
                            GRAD[pstar.idx, sl.idx_j] <- colSums(SC.COR.UNI$dx.sl.y2)
                        } else {
                            GRAD[pstar.idx, sl.idx_j] <- sum(SC.COR.UNI$dx.sl.y2)
                        }
                    }
                    # NO VAR

                    # RHO
                    GRAD[pstar.idx,cor.idx] <- sum(SC.COR.UNI$dx.rho)
                }
            }
        }
    }

    # do we need scores?
    if(scores) return(SCORES)

    # gradient is sum over all pairs
    gradient <- colSums(GRAD)

    # we multiply by -1 because we minimize
    if(negative) {
        gradient <- -1 * gradient
    }

    gradient
}

