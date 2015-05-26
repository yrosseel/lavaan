# univariate modification indices
#

modindices <- function(object, 
                       standardized = TRUE, 

                       # power statistics?
                       power = FALSE, 
                       delta = 0.1, 
                       alpha = 0.05, 
                       high.power = 0.75,

                       # customize output
                       sort. = FALSE, 
                       minimum.value = 0.0, 
                       maximum.number = nrow(LIST),
                       na.remove = FALSE, 
                       op = NULL) {

    # check if model has converged
    if(object@Fit@npar > 0L && !object@Fit@converged) {
        warning("lavaan WARNING: model did not converge")
    }

    # not ready for estimator = "PML"
    if(object@Options$estimator == "PML") {
        stop("lavaan WARNING: modification indices for estimator PML are not implemented yet.")
    }
 
    # sanity check
    if(power) standardized <- TRUE

    # user-specified model parameters
    partable <- object@ParTable[c("lhs","op","rhs","group","free","exo")]
    # replace 'start' column, since lav_model will fill these in in GLIST
    partable$start <- parameterEstimates(object, 
                          remove.eq = FALSE, remove.ineq = FALSE)$est 

    # extended list (fixed-to-zero parameters)
    strict.exo <- FALSE
    if(object@Model@fixed.x && object@Model@categorical) {
        strict.exo <- TRUE ## truly conditional.x
    }
    FULL <- lav_partable_full(object@ParTable, free = TRUE, start = TRUE,
                              strict.exo = strict.exo)

    # merge
    LIST <- lav_partable_merge(partable, FULL, remove.duplicated = TRUE,
                               warn = FALSE)

    # remove  ==, <, :=, > rows from partable
    nonpar.idx <- which(LIST$op %in% c("==", ":=", "<", ">"))
    if(length(nonpar.idx) > 0L) {
        LIST <- LIST[-nonpar.idx,]   
    }

    # create lavmodel object for this 'full' LIST
    LIST2 <- LIST; LIST2$free <- 1:nrow(LIST)

    # reconstruct th.idx
    th.idx <- vector("list", length=object@Data@ngroups)
    for(g in 1:object@Data@ngroups) {
        th.idx[[g]] <- lav_partable_ov_idx(LIST2, type="th", group = g)
    }
    lavmodelFULL <- lav_model(lavpartable = LIST2,
                              representation = object@Model@representation,
                              th.idx = th.idx,
                              parameterization = object@Model@parameterization,
                              link = object@Model@link,
                              debug = FALSE)
    LIST$start <- NULL
    LIST$exo <- NULL
                              
    # compute information matrix 'full'
    # ALWAYS use *expected* information (for now)
    E <- 
       lav_model_information_expected(lavmodel       = lavmodelFULL,
                                      lavsamplestats = object@SampleStats,
                                      estimator      = object@Options$estimator)
    Q <- (1/object@SampleStats@ntotal) * E

    # compute gradient 'full'
    dx <- lav_model_gradient(lavmodel       = lavmodelFULL,
                             GLIST          = NULL, 
                             lavsamplestats = object@SampleStats,
                             lavdata        = object@Data,
                             lavcache       = object@Cache,
                             type           = "free",
                             estimator      = object@Options$estimator,
                             group.weight   = TRUE)

    # Saris, Satorra & Sorbom 1987
    # partition Q into Q_11, Q_22 and Q_12/Q_21
    # which elements of Q correspond with 'free' and 'nonfree' parameters?

    # NOTE: since 0.5-18, we do not make a distinction anymore between
    #       equality constrained parameters (using labels), and general
    #       linear constraints (using ==)
    #       As a result, the 'free' parameters are not necessarily a clean
    #       subset of the total set of parameters
    #
    #       But we can compute what would happen if we release a constraint,
    #       one at a time
    #       FIXME!!!

    model.idx <- which(LIST$free  > 0L)
    extra.idx <- which(LIST$free == 0L)

    # partition Q
    Q11 <- Q[extra.idx, extra.idx, drop = FALSE]
    Q12 <- Q[extra.idx, model.idx, drop = FALSE]
    Q21 <- Q[model.idx, extra.idx, drop = FALSE]
    Q22 <- Q[model.idx, model.idx, drop = FALSE]

    #Q22.inv <- vcov(object) * (nobs(object)) * (nobs(object))
    # ALWAYS use *expected* information (for now)
    Q22.inv <- 
        lavTech(object, "inverted.information.expected") * nobs(object)

    V <- Q11 - Q12 %*% Q22.inv %*% Q21
    V.diag <- diag(V)
    # dirty hack: catch very small or negative values in diag(V)
    # this is needed eg when parameters are not identified if freed-up;
    idx <- which(V.diag < 1.0e-15); V.diag[idx] <- as.numeric(NA)

    # create and fill in mi
    mi <- numeric( length(dx) )
    mi[extra.idx] <- dx[extra.idx]*dx[extra.idx] / V.diag
    if(length(model.idx) > 0L) {
        mi[model.idx]    <- dx[model.idx]*dx[model.idx] / diag(Q22)
    }


    # EPC
    d <- (-1 * object@SampleStats@ntotal) * dx
    # needed?
    d[which(abs(d) < 1e-15)] <- 1.0
    epc <- mi/d

    # FIXME: epc for equality constraints must be adapted!?
    # we need a reference for this!!!!
    #if(object@Model@eq.constraints) {
    #    for(i in 1:length(eq.idx)) {
    #        # this is just a temporary solution
    #        # until we get it right
    #        neq <- length(eq.idx[eq.id == eq.id[i]])
    #        epc[eq.idx[i]] <- epc[eq.idx[i]] / neq * (neq - 1)
    #    }
    #}

    LIST$mi <- mi
    if(length(object@Fit@test) > 1L) {
        LIST$mi.scaled <- mi / object@Fit@test[[2]]$scaling.factor
    }
    LIST$epc <- epc

    # remove some rows
    #idx <- which(LIST$free > 0L & !duplicated(LIST$free) & !LIST$eq.id > 0L)
    #LIST <- LIST[-idx,]
    #if(length(model.idx) > 0L) {
    #    LIST <- LIST[-model.idx,]
    #}

    # standardize?
    if(standardized) {
        # two problems: 
        #   - EPC of variances can be negative, and that is
        #     perfectly legal
        #   - EPC (of variances) can be tiny (near-zero), and we should 
        #     not divide by tiny variables
 
        EPC <- LIST$epc
        small.idx <- which(LIST$op == "~~" & 
                           LIST$lhs == LIST$rhs &
                           abs(EPC) < sqrt( .Machine$double.eps ) )
        if(length(small.idx) > 0L) {
            EPC[ small.idx ] <- as.numeric(NA)
        }

        # get the sign
        EPC.sign <- sign(LIST$epc)

        LIST$sepc.lv <- EPC.sign * standardize.est.lv(object, 
                                                      partable = LIST, 
                                                      est = abs(EPC),
                                                      cov.std = FALSE)
        if(length(small.idx) > 0L) {
            LIST$sepc.lv[small.idx] <- 0
        }
        LIST$sepc.all <- EPC.sign * standardize.est.all(object, 
                                                        partable = LIST, 
                                                        est = abs(EPC),
                                                        cov.std = FALSE)
        if(length(small.idx) > 0L) {
            LIST$sepc.all[small.idx] <- 0
        }
        LIST$sepc.nox <- EPC.sign * standardize.est.all.nox(object, 
                                                            partable = LIST,
                                                            est = abs(EPC),
                                                            cov.std = FALSE)
        if(length(small.idx) > 0L) {
            LIST$sepc.nox[small.idx] <- 0
        }
 
    }

    # power?
    if(power) {
        LIST$delta <- delta
        # FIXME: this is using epc in unstandardized metric
        #        this would be much more useful in standardized metric
        #        we need a standardize.est.all.reverse function...
        LIST$ncp <- (LIST$mi / LIST$epc*LIST$epc) * (delta*delta)
        LIST$power <- 1 - pchisq(qchisq((1.0 - alpha), df=1), 
                                 df=1, ncp=LIST$ncp)
        LIST$decision <- character( length(LIST$power) )

        # four possibilities (Table 6 in Saris, Satorra, van der Veld, 2009)
        mi.significant <- ifelse( 1 - pchisq(LIST$mi, df=1) < alpha,
                                  TRUE, FALSE )
        high.power <- LIST$power > high.power

        LIST$decision[ which(mi.significant &  high.power) ] <- "epc"
        LIST$decision[ which(mi.significant & !high.power) ] <- "***"
        LIST$decision[ which(!mi.significant & !high.power) ] <- "(i)"
    }


    # remove some columns
    LIST$free <- NULL
    LIST$mat <- LIST$row <- LIST$col <- LIST$id <- NULL
    #if(power) {
    #    LIST$epc <- NULL
    #    LIST$sepc.lv <- NULL
    #}

    class(LIST) <- c("lavaan.data.frame", "data.frame")

    # add eq constraints (if any)
    eq.idx <- which(partable$op == "==")
    if(length(eq.idx) > 0L) {
        warning("lavaan WARNING: modification indices do not reflect (yet) what happens if we release an equality constraint")
        # FIXME!!!!
        N <- length(eq.idx)
        TMP <- data.frame(lhs = partable$lhs[eq.idx],
                           op = partable$op[eq.idx],
                          rhs = partable$rhs[eq.idx],
                          group = partable$group[eq.idx],
                          mi = rep(as.numeric(NA), N),
                          epc = rep(as.numeric(NA), N),
                          sepc.lv = rep(as.numeric(NA), N),
                          sepc.all = rep(as.numeric(NA), N),
                          sepc.nox = rep(as.numeric(NA), N) )
        if(!standardized) {
            TMP$sepc.lv <- TMP$sepc.all <- TMP$sepc.nox <- NULL
        }
        if(!is.null(LIST$mi.scaled)) {
            TMP$mi.scaled <- rep(as.numeric(NA), N)
        }
        LIST <- rbind(LIST, TMP)
    }

    # remove 'existing' parameters
    #TMP1 <- data.frame(lhs = partable$lhs,
    #                   op = partable$op,
    #                   rhs = partable$rhs,
    #                   group = partable$group)
    #TMP2 <- data.frame(lhs = LIST$lhs, 
    #                   op = LIST$op,
    #                   rhs = LIST$rhs,
    #                   group = LIST$group)
    #idx <- which(duplicated(rbind(TMP2, TMP1), fromLast = TRUE))
    #if(length(idx) > 0L) {
    #    LIST <- LIST[-idx,]
    #}
    if(max(LIST$group) == 1) LIST$group <- NULL
    
    # sort?
    if(sort.) {
        LIST <- LIST[order(LIST$mi, decreasing = TRUE),]
    }
    if(minimum.value > 0.0) {
        LIST <- LIST[!is.na(LIST$mi) & LIST$mi > minimum.value,]  
    }
    if(maximum.number < nrow(LIST)) {
        LIST <- LIST[seq_len(maximum.number),]
    }
    if(na.remove) {
        idx <- which(is.na(LIST$mi))
        if(length(idx) > 0) {
            LIST <- LIST[-idx,]
        }
    }
    if(!is.null(op)) {
        idx <- LIST$op %in% op
        if(length(idx) > 0) {
            LIST <- LIST[idx,]
        }
    }

    LIST
}

# aliases
modificationIndices <- modificationindices <- modindices
