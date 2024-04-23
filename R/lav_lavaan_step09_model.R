lav_lavaan_step09_model <- function(slotModel = NULL, # nolint
                                    lavoptions = NULL,
                                    lavpartable = NULL,
                                    lavsamplestats = NULL,
                                    lavdata = NULL) {
  # # # # # # # # # # #
  # #  9. lavmodel # #
  # # # # # # # # # # #

  # if slotModel not NULL
  #   copy to lavmodel
  # else
  #   compute lavmodel via lav_model
  #   if lavdata@data.type == "none" and categorical mode
  #     set parameters in lavmodel via lav_model_set_parameters and
  #       re-adjust start column in lavpartable
  #     if differences between start and ustart column (in lavpartable)
  #       if lavmodel$parameterization == "delta"
  #         if user specified delta values : ** warning **
  #       if lavmodel$parameterization == "theta"
  #         if user specified theta values : ** warning **

  if (!is.null(slotModel)) {
    lavmodel <- slotModel
  } else {
    if (lavoptions$verbose) {
      cat("lavmodel           ...")
    }
    lavmodel <- lav_model(
      lavpartable = lavpartable,
      lavoptions = lavoptions,
      th.idx = lavsamplestats@th.idx
    )
    # no longer needed: x values are in start
    # cov.x            = lavsamplestats@cov.x,
    # mean.x           = lavsamplestats@mean.x)

    # if no data, call lav_model_set_parameters once (for categorical case)
    if (lavdata@data.type == "none" && lavmodel@categorical) {
      lavmodel <- lav_model_set_parameters(
        lavmodel = lavmodel,
        x = lav_model_get_parameters(lavmodel)
      )
      # re-adjust parameter table
      lavpartable$start <- lav_model_get_parameters(lavmodel, type = "user")

      # check/warn if theta/delta values make sense
      if (!all(lavpartable$start == lavpartable$ustart)) {
        if (lavmodel@parameterization == "delta") {
          # did the user specify theta values?
          user.var.idx <- which(lavpartable$op == "~~" &
            lavpartable$lhs == lavpartable$rhs &
            lavpartable$lhs %in% unlist(attr(lavpartable, "vnames")$ov.ord) &
            lavpartable$user == 1L)
          if (length(user.var.idx)) {
            lav_msg_warn(
              gettextf("variance (theta) values for categorical variables
                       are ignored if parameterization = %s!",
                       "'delta'")
            )
          }
        } else if (lavmodel@parameterization == "theta") {
          # did the user specify theta values?
          user.delta.idx <- which(lavpartable$op == "~*~" &
            lavpartable$lhs == lavpartable$rhs &
            lavpartable$lhs %in% unlist(attr(lavpartable, "vnames")$ov.ord) &
            lavpartable$user == 1L)
          if (length(user.delta.idx)) {
            lav_msg_warn(
              gettextf("scaling (~*~) values for categorical variables
                       are ignored if parameterization = %s!",
                       "'theta'")
            )
          }
        }
      }
    }
    if (lavoptions$verbose) {
      cat(" done.\n")
    }
  }

  list(
    lavpartable = lavpartable,
    lavmodel    = lavmodel
  )
}

# 9b. bounds for EFA -- to force diag(LAMBDA) to be positive (new in 0.6-7)
# if((.hasSlot(lavmodel, "nefa")) && (lavmodel@nefa > 0L) &&
#    (lavoptions$rotation != "none")) {
#
#    # add lower column
#    if (is.null(lavpartable$lower)) {
#        lavpartable$lower <- rep(-Inf, length(lavpartable$lhs))
#    }
#    efa.values <- lav_partable_efa_values(lavpartable)
#    group.values <- lav_partable_group_values(lavpartable)
#    for (g in seq_len(lavdata@ngroups)) {
#        for (set in seq_len(lavmodel@nefa)) {
#            lv.efa <-
#                unique(lavpartable$lhs[lavpartable$op == "=~" &
#                                       lavpartable$block == g &
#                                       lavpartable$efa == efa.values[set] ])
#            for (f in seq_len(length(lv.efa))) {
#                lambda.idx <- which(lavpartable$lhs == lv.efa[f] &
#                                     lavpartable$op == "=~" &
#                                     lavpartable$group == group.values[g])
#                # get diagonal element of LAMBDA
#                midx <- lambda.idx[f] # diagonal element of LAMBDA
#                lavpartable$lower[midx] <- 0
#             } # factors
#        } # sets
#    } # groups
# }
