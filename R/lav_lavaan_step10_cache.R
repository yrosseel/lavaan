lav_lavaan_step10_cache <- function(slotCache, lavdata, lavmodel, lavpta, lavoptions, sampling.weights) {
  # # # # # # # # # # # 
  # #  10. lavcache # # 
  # # # # # # # # # # # 
  # if slotCache not NULL
  #   copy to lavcache
  # else
  #   lavcache = list of length lavdata@ngroups
  #   set tmp.ov.types = lavdata$ov$types
  #   if lavmodel@conditional.x and sum(lavmodel@nexo) > 0L remove elements lavpta$vidx$ov.x from tmp.ov.types
  #   if lavoptions$estimator == "PML" and all tmp.ov.types are "ordered"
  #     TH = computeTH(lavmodel)
  #     BI = lav_tables_pairwise_freq_cells(lavdata)
  #     if lavoptions$missing is "available.cases" or "doubly.robust" 
  #       UNI = lav_tables_univariate_freq_cell(lavdata)
  #       if lavoptions$missing is "doubly.robust"
  #         if lavoptions$control$pairwiseProbGivObs NULL: *** error ***
  #         if lavoptions$control$univariateProbGivObs NULL: *** error ***
  #     for all groups (1:lavdata@ngroups)
  #       set tmp.idx = 1:length(BI$ibs.freq)
  #       if BI$group not NULL and max(BI$group) > 1L set tmp.idx = indexes for this group in BI
  #       set bifreq = BI$obs.freq[tmp.idx]
  #       set binobs = BI$nobs[tmp.idx]
  #       set LONG = LongVecInd(no.x = ncol(lavdata@X[[g]]), all.thres = TH[[g]],
  #                             index.var.of.thres = lavmodel@th.idx[[g]])
  #       set lavcache[[g]] = list(bifreq = bifreq, nobs = binobs, LONG = LONG)
  #       if sampling.weights not NULL
  #         compute (for group g) lavcache[[g]]$sum_obs_weights_xixj_ab_vec (*)
  #       if lavoptions$missing is "available.cases" or "doubly.robust" 
  #         set tmp.idx = 1:length(BI$ibs.freq)
  #         if BI$group not NULL and max(BI$group) > 1L set tmp.idx = indexes for this group in BI
  #         set lavcache[[g]]$unifreq = unifreq = UNI$obs.freq[tmp.idx]
  #         set lavcache[[g]]$uninobs = uninobs = UNI$nobs[tmp.idx]
  #         set lavcache[[g]]$uniweights.casewise = uniweights.casewise = rowSums(is.na(lavdata@X[[g]]))
  #         compute lavcache[[g]]$uniweights (*)
  #       if lavoptions$missing is "doubly.robust" 
  #         lavcache[[g]]$pairwiseProbGivObs = lavoptions$control$pairwiseProbGivObs[[g]]
  #         lavcache[[g]]$univariateProbGivObs = lavoptions$control$univariateProbGivObs[[g]]
  #         compute members idx.Y1, idx.Gy2, idx.cat.Y1, idx.cat.Gy2 and id.uniPrGivObs from
  #           lavchache[[g]] (*)
  #   if lavdata$data.type is "full" and lavdata@Rp[[1L]] not NULL
  #     copy lavdata@Rp[[g]]$pat to lavcache[[g]]$pat for all groups g
  # if lavoptions$estimator is "MML"
  #   compute for all groups g lavcache[[g]]$GH via lav_integration_gauss_hermite
  #         
  # (*) !!! computations too complicated to summarize here !!!
  if (!is.null(slotCache)) {
    lavcache <- slotCache
  } else {
    # prepare cache -- stuff needed for estimation, but also post-estimation
    lavcache <- vector("list", length = lavdata@ngroups)
    
    # ov.types? (for PML check)
    tmp.ov.types <- lavdata@ov$type
    if (lavmodel@conditional.x && sum(lavmodel@nexo) > 0L) {
      # remove ov.x
      tmp.ov.x.idx <- unlist(lavpta$vidx$ov.x)
      tmp.ov.types <- tmp.ov.types[-tmp.ov.x.idx]
    }
    
    if (lavoptions$estimator == "PML" && all(tmp.ov.types == "ordered")) {
      TH <- computeTH(lavmodel)
      BI <- lav_tables_pairwise_freq_cell(lavdata)
      
      # handle option missing = "available.cases" or "doubly.robust"
      if (lavoptions$missing == "available.cases" ||
          lavoptions$missing == "doubly.robust") {
        UNI <- lav_tables_univariate_freq_cell(lavdata)
        # checks for missing = "double.robust"
        if (lavoptions$missing == "doubly.robust") {
          # check whether the probabilities pairwiseProbGivObs and
          # univariateProbGivObs are given by the user
          if (is.null(lavoptions$control$pairwiseProbGivObs)) {
            stop("lavaan ERROR: could not find `pairwiseProbGivObs' in control() list")
          }
          if (is.null(lavoptions$control$univariateProbGivObs)) {
            stop("lavaan ERROR: could not find `univariateProbGivObs' in control() list")
          }
        }
      }

      for (g in 1:lavdata@ngroups) {
        if (is.null(BI$group) || max(BI$group) == 1L) {
          bifreq <- BI$obs.freq
          binobs  <- BI$nobs
        } else {
          idx <- which(BI$group == g)
          bifreq <- BI$obs.freq[idx]
          binobs  <- BI$nobs[idx]
        }
        LONG  <- LongVecInd(no.x               = ncol(lavdata@X[[g]]),
                            all.thres          = TH[[g]],
                            index.var.of.thres = lavmodel@th.idx[[g]])
        lavcache[[g]] <- list(bifreq = bifreq,
                              nobs   = binobs,
                              LONG   = LONG)
        
        # >>>>>>>> HJ/MK PML CODE >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        # I need to add something that splits weights into g groups so
        # adjust what follows in the new code also compute the sum of
        # weights within a group, this will substitute n_g (group size)
        # of simple random sampling (SRS) and also compute the total the
        # total sum of weights over all observation over all groups,
        # this substitutes the total sample size of SRS.
        
        if (!is.null(sampling.weights)) {
          # Keep track of indices of the response categories (a,b) of a
          # pair of ordinal variables (xi,xj) appearing in the data as
          # well as the index of the pair.
          idx_ab_of_xixj_ab <- lapply(LONG[c(1:2, 5)], function(x) {
            x[(LONG$index.thres.var1.of.pair != 0) &
                (LONG$index.thres.var2.of.pair != 0)]
          })
          names(idx_ab_of_xixj_ab) <- c("idx_a", "idx_b", "idx_pairs")
          lavcache[[g]]$idx_ab_of_xixj_ab <- idx_ab_of_xixj_ab
          
          # Raw data for group g
          X_g <- lavdata@X[[g]]
          
          # I assume that X_g includes only the ordinal indicators nvar
          # gives the number of ordinal indicators
          nvar <- ncol(X_g)
          
          # pstar gives the number of pairs formed by the nvar ordinal
          # indicators
          pstar <- nvar * (nvar - 1) / 2
          
          # Keep track of the indices of variables forming each pair
          idx_vars_in_pair <- combn(nvar, 2)
          
          # The output of sapply below provides the sum of weights for
          # all bivariate response pattern for all pairs of indicators.
          
          # If all indicators have the same number of response
          # categories, the output of sapply function below is a matrix.
          # Each column refers to a different pair of indicators (i,j)
          # with j running faster than i, e.g. (1,2) (1,3) (2,3). Within
          # each column, each element (i.e. each row of the matrix)
          # refers to a different combination of response categories
          # (a,b) with a, the category index of indicator i, running
          # faster than b, the category index of indicator j, e.g.
          # (1,1), (2,1) (3,1) (1,2) (2,2) (3,2)
          
          # If the indicators have different number of response
          # categories, the output of sapply function below is a list.
          # Each element of the list refers to a different pair of
          # indicators (i,j) with j running faster than i and it is a
          # matrix with number of rows the number of response categories
          # of indicator i and ncol =  the number of response categories
          # of indicator j.
          
          sum_obs_weights_xixj_ab <- sapply(1:pstar, function(x) {
            tmp_idx_ab <- lapply(idx_ab_of_xixj_ab, function(y) {
              y[idx_ab_of_xixj_ab$idx_pairs == x]})
            tmp_idx_cols <- idx_vars_in_pair[, x]
            tmp_var1 <- factor(X_g[, tmp_idx_cols[1]], levels =
                                 as.character(unique(tmp_idx_ab$idx_a)))
            tmp_var2 <- factor(X_g[, tmp_idx_cols[2]], levels =
                                 as.character(unique(tmp_idx_ab$idx_b)))
            tapply(X = lavdata@weights[[g]],
                   INDEX = list(tmp_var1, tmp_var2),
                   FUN = sum)
          })
          
          # We need to transform the output of sapply into a vector
          # where the sum of weights (for all bivariate response
          # patterns for all pairs of indicators) are listed in the same
          # order as in pairwisePI vector, i.e. a runs the fastest,
          # followed by b, then by j and lastly by i.
          
          if (is.matrix(sum_obs_weights_xixj_ab)) {
            sum_obs_weights_xixj_ab_vec <- c(sum_obs_weights_xixj_ab)
          } else if (is.list(sum_obs_weights_xixj_ab)) {
            sum_obs_weights_xixj_ab_vec <-
              do.call(c, sum_obs_weights_xixj_ab)
          }
          
          # Note that sapply gives NA for these bivariate response
          # patterns which are not observed at all. Substitute NA with
          # 0.
          idx_na_sowxav <- is.na(sum_obs_weights_xixj_ab_vec)
          if (any(idx_na_sowxav)) {
            sum_obs_weights_xixj_ab_vec[idx_na_sowxav] <- 0
          }
          
          lavcache[[g]]$sum_obs_weights_xixj_ab_vec <-
            sum_obs_weights_xixj_ab_vec
        }
        
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        # available cases
        if (lavoptions$missing == "available.cases" ||
            lavoptions$missing == "doubly.robust") {
          if (is.null(UNI$group) || max(UNI$group) == 1L) {
            unifreq <- UNI$obs.freq
            uninobs <- UNI$nobs
          } else {
            idx <- which(UNI$group == g)
            unifreq <- UNI$obs.freq[idx]
            uninobs <- UNI$nobs[idx]
          }
          lavcache[[g]]$unifreq <- unifreq
          lavcache[[g]]$uninobs <- uninobs
          
          uniweights.casewise <- rowSums(is.na(lavdata@X[[g]]))
          lavcache[[g]]$uniweights.casewise <- uniweights.casewise
          
          #weights per response category per variable in the same
          # order as unifreq; i.e. w_ia, i = 1,...,p, (p variables),
          # a = 1,...,Ci, (Ci response categories for variable i),
          # a running faster than i
          tmp.uniweights <- apply(lavdata@X[[g]], 2,
                                  function(x) {
                                    tapply(uniweights.casewise, as.factor(x), sum,
                                           na.rm = TRUE) })
          if (is.matrix(tmp.uniweights)) {
            lavcache[[g]]$uniweights <- c(tmp.uniweights)
          }
          if (is.list(tmp.uniweights)) {
            lavcache[[g]]$uniweights <- unlist(tmp.uniweights)
          }
        } # "available.cases" or "double.robust"
        
        # doubly.robust only
        if (lavoptions$missing == "doubly.robust") {
          
          # add the provided by the user probabilities
          # pairwiseProbGivObs and univariateProbGivObs in Cache
          lavcache[[g]]$pairwiseProbGivObs <-
            lavoptions$control$pairwiseProbGivObs[[g]]
          lavcache[[g]]$univariateProbGivObs <-
            lavoptions$control$univariateProbGivObs[[g]]
          # compute different indices vectors that will help to do
          # calculations
          ind.vec <- as.data.frame(LONG[1:5])
          ind.vec <-
            ind.vec[((ind.vec$index.thres.var1.of.pair != 0) &
                       (ind.vec$index.thres.var2.of.pair != 0)) , ]
          idx.cat.y1 <- ind.vec$index.thres.var1.of.pair
          idx.cat.y2 <- ind.vec$index.thres.var2.of.pair
          idx.pairs  <- ind.vec$index.pairs.extended
          lavcache[[g]]$idx.pairs <- idx.pairs
          
          idx.cat.y1.split <- split(idx.cat.y1, idx.pairs)
          idx.cat.y2.split <- split(idx.cat.y2, idx.pairs)
          lavcache[[g]]$idx.cat.y1.split <- idx.cat.y1.split
          lavcache[[g]]$idx.cat.y2.split <- idx.cat.y2.split
          
          # generate the variables, categories indices vector which
          # keep track to which variables and categories the
          # elements of vector probY1Gy2 refer to
          nlev <- lavdata@ov$nlev
          nvar <- length(nlev)
          
          idx.var.matrix <- matrix(1:nvar, nrow = nvar, ncol = nvar)
          idx.diag <- diag(matrix(1:(nvar * nvar), nrow = nvar,
                                  ncol = nvar))
          idx.Y1Gy2.matrix <- rbind(t(idx.var.matrix)[-idx.diag],
                                    idx.var.matrix [-idx.diag])
          no.pairs.Y1Gy2 <- ncol(idx.Y1Gy2.matrix)
          idx.cat.Y1 <- unlist(lapply(1:no.pairs.Y1Gy2, function(x) {
            rep(1:nlev[idx.Y1Gy2.matrix[1, x]],
                times = nlev[idx.Y1Gy2.matrix[2, x]])}))
          idx.cat.Gy2 <- unlist(lapply(1:no.pairs.Y1Gy2, function(x) {
            rep(1:nlev[idx.Y1Gy2.matrix[2, x]],
                each = nlev[idx.Y1Gy2.matrix[1, x]])}))
          dim.pairs <- unlist(lapply(1:no.pairs.Y1Gy2, function(x) {
            nlev[idx.Y1Gy2.matrix[1, x]] *
              nlev[idx.Y1Gy2.matrix[2, x]] }))
          idx.Y1 <- unlist(mapply(rep, idx.Y1Gy2.matrix[1, ],
                                  each = dim.pairs))
          idx.Gy2 <- unlist(mapply(rep, idx.Y1Gy2.matrix[2, ],
                                   each = dim.pairs))
          
          lavcache[[g]]$idx.Y1      <- idx.Y1
          lavcache[[g]]$idx.Gy2     <- idx.Gy2
          lavcache[[g]]$idx.cat.Y1  <- idx.cat.Y1
          lavcache[[g]]$idx.cat.Gy2 <- idx.cat.Gy2
          
          # the vector below keeps track of the variable each column
          # of the matrix univariateProbGivObs refers to
          lavcache[[g]]$id.uniPrGivObs <-
            sort(c(unique(lavmodel@th.idx[[g]]),
                   lavmodel@th.idx[[g]]))
        } # doubly.robust

      } # g
    }
    # copy response patterns to cache -- FIXME!! (data not included
    # in Model only functions)
    if (lavdata@data.type == "full" && !is.null(lavdata@Rp[[1L]])) {
      for (g in 1:lavdata@ngroups) {
        lavcache[[g]]$pat <- lavdata@Rp[[g]]$pat
      }
    }
  }
  
  # If estimator = MML, store Gauss-Hermite nodes/weights
  if (lavoptions$estimator == "MML") {
    for (g in 1:lavdata@ngroups) {
      # count only the ones with non-normal indicators
      #nfac <- lavpta$nfac.nonnormal[[g]]
      nfac <- lavpta$nfac[[g]]
      lavcache[[g]]$GH <-
        lav_integration_gauss_hermite(n = lavoptions$integration.ngh,
                                      dnorm = TRUE,
                                      mean = 0, sd = 1,
                                      ndim = nfac)
      #lavcache[[g]]$DD <- lav_model_gradient_DD(lavmodel, group = g)
    }
  }
  return(lavcache)
}
