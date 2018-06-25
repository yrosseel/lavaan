lav_samplestats_wls_obs <- function(mean.g, cov.g, var.g,
                                    th.g, th.idx.g,
                                    res.int.g, res.cov.g, res.var.g, res.th.g,
                                    res.slopes.g,
                                    group.w.g,
                                    categorical    = FALSE,
                                    conditional.x  = FALSE,
                                    meanstructure  = FALSE,
                                    slopestructure = FALSE,
                                    group.w.free   = FALSE) {

    # WLS.obs
    if(categorical) {
        # order of elements is important here:
        # 1. thresholds + (negative) means (interleaved)
        # 2. slopes (if any)
        # 3. variances (if any)
        # 4. covariance matrix (no diagonal!)

        # NOTE: prior to 0.5-17, we had this:
        # TH[ov.types == "numeric"] <- -1*TH[ov.types == "numeric"]
        # which is WRONG if we have more than one threshold per variable
        # (thanks to Sacha Epskamp for spotting this!)
        if(conditional.x) {
            TH <- res.th.g
            TH[ th.idx.g == 0 ] <- -1*TH[ th.idx.g == 0 ]

            nvar <- length(res.var.g)
            num.idx <- which(!seq_len(nvar) %in% th.idx.g)

            WLS.obs <- c(TH,
                         lav_matrix_vec(res.slopes.g),
                         res.var.g[ num.idx ],
                         lav_matrix_vech(res.cov.g, diagonal = FALSE) )
        } else {
            TH <- th.g
            TH[ th.idx.g == 0 ] <- -1*TH[ th.idx.g == 0 ]

            nvar <- length(var.g)
            num.idx <- which(!seq_len(nvar) %in% th.idx.g)

            WLS.obs <- c(TH,
                         var.g[ num.idx ],
                         lav_matrix_vech(cov.g, diagonal = FALSE) )
        }
    } else {

        if(conditional.x) {
            if(meanstructure) {
                if(slopestructure) {
                    # order = vec(Beta), where first row are intercepts
                    # cbind(res.int, res.slopes) is t(Beta)
                    # so we need vecr
                    WLS.obs <- c( lav_matrix_vecr( cbind(res.int.g,
                                                         res.slopes.g)),
                                  lav_matrix_vech( res.cov.g )
                                )
                    #WLS.obs <- c(                 res.int.g,
                    #             lav_matrix_vec(  res.slopes.g  ),
                    #             lav_matrix_vech( res.cov.g     ))
                } else {
                    WLS.obs <- c(                 res.int.g,
                                 lav_matrix_vech( res.cov.g     ))
                }
            } else {
                if(slopestructure) {
                    WLS.obs <- c(lav_matrix_vecr( res.slopes.g ),
                                 lav_matrix_vech( res.cov.g    ))
                } else {
                    WLS.obs <-   lav_matrix_vech( res.cov.g    )
                }
            }

        } else {
            if(meanstructure) {
                 WLS.obs <- c(                 mean.g,
                              lav_matrix_vech(  cov.g  ))
            } else {
                 WLS.obs <-   lav_matrix_vech(  cov.g  )
            }
        }
    }

    # group.w.free?
    if(group.w.free) {
        #group.w.last <- nobs[[ngroups]] / sum(unlist(nobs))
        #WLS.obs[[g]] <- c(log(group.w[[g]]/group.w.last), WLS.obs[[g]])
        WLS.obs <- c(group.w.g, WLS.obs)
    }

    WLS.obs
}
