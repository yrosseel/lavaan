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
    if (lav_verbose()) {
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

    # same for composites: call lav_model_set_parameters once to set
    # total/residual variances of composites in PSI
    } else if (lavmodel@composites) {
      lavmodel <- lav_model_set_parameters(
        lavmodel = lavmodel,
        x = lav_model_get_parameters(lavmodel)
      )
      # re-adjust parameter table
      lavpartable$start <- lav_model_get_parameters(lavmodel, type = "user")
    }
    if (lav_verbose()) {
      cat(" done.\n")
    }
  }

  # if parameterization = "delta" and categorical/correlation: check if
  # we have an observed mediator (new in 0.6-19)
  # but fixed (for recursive models) in 0.6-20
#   check.flag <- FALSE
#   if (!is.null(lavoptions$check.delta.cat.mediator) && # lavaan >= 0.6.19
#       lavoptions$check.delta.cat.mediator) {
#     check.flag <- TRUE
#   }
#   if (check.flag &&
#       (lavmodel@categorical || lavmodel@correlation) &&
#       lavmodel@representation == "LISREL" &&
#       lavmodel@parameterization == "delta") {
#     # get idx BETA matrices
# 	beta.idx   <- which(names(lavmodel@GLIST) == "beta")
# 	# for every block
# 	for (i in seq_len(length(beta.idx))) {
# 	  this.beta <- abs(lavmodel@GLIST[[beta.idx[i]]])
#       this.beta[lavmodel@m.free.idx[[beta.idx[i]]]] <- 1.0
# 	  # exogenous variables:       have zero rows
# 	  # endogenous variables:      have non-zero rows
# 	  # endogenous-only variables: have non-zero rows and zero columns
# 	  # mediators:                 have non-zero rows and non-zero columns
# 	  m.idx <- which(apply(this.beta, 1, sum) != 0 &
#                      apply(this.beta, 2, sum) != 0)
# 	  if (length(m.idx) > 0L) {
# 	    # we have (at least) one mediator
# 	    # is one of them an observed variable? -> warn
#         m.names <- lavmodel@dimNames[[beta.idx[i]]][[1]][m.idx]
# 		ov.names   <- unlist(lavdata@ov.names)
#         ov.ordered <- unlist(lavdata@ordered)
#         # correlation model
# 		if (lavmodel@correlation && any(m.names %in% ov.names)) {
# 		  bad.names <- m.names[m.names %in% ov.names]
# 		  if (length(beta.idx) == 1L) {
#             lav_msg_warn(gettextf("model contains at least one observed mediator: [%s]; consider switching to parameterization = \"theta\"",
#             paste(bad.names, collapse = " ")))
# 		  } else {
#             lav_msg_warn(gettextf("model contains at least one observed mediator in block %i: [%s]; consider switching to parameterization = \"theta\"", i,
#             paste(bad.names, collapse = " ")))
# 		  }
#         # categorical mode
# 		} else if (lavmodel@categorical && any(m.names %in% ov.ordered)) {
#           bad.names <- m.names[m.names %in% ov.ordered]
#           if (length(beta.idx) == 1L) {
#             lav_msg_warn(gettextf("model contains at least one observed categorical mediator: [%s]; consider switching to parameterization = \"theta\"",
#             paste(bad.names, collapse = " ")))
#           } else {
#             lav_msg_warn(gettextf("model contains at least one observed categorical mediator in block %i: [%s]; consider switching to parameterization = \"theta\"", i,
#             paste(bad.names, collapse = " ")))
#           }
#         }
# 	  } # mediators
# 	} # block
#   } # delta parameterization

  list(
    lavpartable = lavpartable,
    lavmodel    = lavmodel
  )
}

