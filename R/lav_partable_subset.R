# YR 11 feb 2017: initial version

# given a parameter table (PT), extract a part of the model:
# eg.:
# - only the measurement model (with saturated latent variables
# - only the stuctural part
# - a single measurement model (1 factor only)
# ...

lav_partable_subset_measurement_model <- function(PT = NULL,
                                                  lavpta = NULL,
                                                  lv.names = NULL) {

    # PT
    PT <- as.data.frame(PT, stringsAsFactors = FALSE)

    # lavpta
    if(is.null(lavpta)) {
        lavpta <- lav_partable_attributes(PT)
    }

    # ngroups
    ngroups <- lavpta$ngroups

    # lv.names: list with element per group
    if(is.null(lv.names)) {
        lv.names <- lavpta$vnames$lv.regular
    } else if(!is.list(lv.names)) {
        lv.names <- list(lv.names)
    }
    
    # which latent variables should we remove?
    lv.names.rm <- lapply(1:ngroups, function(g) {
                       lavpta$vnames$lv.regular[[g]][ 
                        !lavpta$vnames$lv.regular[[g]] %in% lv.names[[g]] ] 
                   })

    # remove rows idx
    rm.idx <- integer(0L)

    # remove not-needed measurement models
    for(g in 1:ngroups) {
        # indicators for not-needed latent variables
        IND.idx <- which(  PT$op == "=~"              &
                          !PT$lhs %in% lv.names[[g]]  &
                           PT$group == g )
        IND <- PT$rhs[ IND.idx ]

        # remove =~
        rm.idx <- c(rm.idx, IND.idx)

        # remove ~~
        VAR.idx <- which( PT$op == "~~"                     &
                          ( PT$lhs %in% IND                 |
                            PT$rhs %in% IND                 |
                            PT$lhs %in% lv.names.rm[[g]]    |
                            PT$rhs %in% lv.names.rm[[g]] )  &
                          PT$group == g )
        rm.idx <- c(rm.idx, VAR.idx)

        # regressions, involving a latent variable
        LV.EQS.idx <- which( PT$op == "~"                               &
                             ( PT$lhs %in% lavpta$vnames$lv.regular[[g]]   |
                               PT$rhs %in% lavpta$vnames$lv.regular[[g]] ) &
                             PT$group == g )
        rm.idx <- c(rm.idx, LV.EQS.idx)

        # regressions, involving indicators
        OV.EQS.idx <- which( PT$op == "~"       &
                            ( PT$lhs %in% IND   |
                              PT$rhs %in% IND ) &
                            PT$group == g )
        rm.idx <- c(rm.idx, OV.EQS.idx)

        # intercepts indicators
        OV.INT.idx <- which( PT$op == "~1"    &
                             PT$lhs %in% IND  &
                             PT$group == g )
        rm.idx <- c(rm.idx, OV.INT.idx)

        # intercepts latent variables
        LV.INT.idx <- which( PT$op == "~1"                 &
                             PT$lhs %in% lv.names.rm[[g]]  &
                             PT$group == g )
        rm.idx <- c(rm.idx, LV.INT.idx)

        # thresholds
        TH.idx <- which( PT$op == "|"    &
                         PT$lhs %in% IND &
                         PT$group == g )
        rm.idx <- c(rm.idx, TH.idx)

        # scaling factors
        SC.idx <- which( PT$op == "~*~"  &
                         PT$lhs %in% IND &
                         PT$group == g )
        rm.idx <- c(rm.idx, SC.idx)

        # FIXME: ==, :=, <, >, == involving IND...
    }
    
    if(length(rm.idx) > 0L) {
        PT <- PT[-rm.idx,,drop = FALSE]
    }

    # clean up
    PT <- lav_partable_complete(PT)

    # check if we have enough indicators?
    # TODO
    
    PT
}

lav_partable_subset_structural_model <- function(PT = NULL,
                                                 lavpta = NULL) {

    # PT
    PT <- as.data.frame(PT, stringsAsFactors = FALSE)

    # lavpta
    if(is.null(lavpta)) {
        lavpta <- lav_partable_attributes(PT)
    }

    # ngroups
    ngroups <- lavpta$ngroups

    # eqs.names
    eqs.x.names <- lavpta$vnames$eqs.x
    eqs.y.names <- lavpta$vnames$eqs.y

    # keep rows idx
    keep.idx <- integer(0L)

    # remove not-needed measurement models
    for(g in 1:ngroups) {

        # eqs.names
        eqs.names <- unique( c(lavpta$vnames$eqs.x[[g]],
                               lavpta$vnames$eqs.y[[g]]) )

        # regressions
        reg.idx <- which(PT$op == "~" & PT$group == g &
                         PT$lhs %in% eqs.names &
                         PT$rhs %in% eqs.names)

        # the variances
        var.idx <- which(PT$op == "~~" & PT$group == g &
                         PT$lhs %in% eqs.names &
                         PT$rhs %in% eqs.names &
                         PT$lhs == PT$rhs)

        # optionally covariances (exo!)
        cov.idx <- which(PT$op == "~~" & PT$group == g &
                         PT$lhs %in% eqs.names &
                         PT$rhs %in% eqs.names &
                         PT$lhs != PT$rhs)

        # means/intercepts
        int.idx <- which(PT$op == "~1" & PT$group == g &
                         PT$lhs %in% eqs.names)

        keep.idx <- c(keep.idx, reg.idx, var.idx, cov.idx, int.idx)
    }

    PT <- PT[keep.idx, , drop = FALSE]

    # clean up
    PT <- lav_partable_complete(PT)

    PT
}
