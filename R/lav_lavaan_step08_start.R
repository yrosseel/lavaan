lav_lavaan_step08_start <- function(slotModel = NULL, # nolint
                                    lavoptions = NULL,
                                    lavpartable = NULL,
                                    lavsamplestats = NULL,
                                    lavh1 = NULL) {
  # # # # # # # # # # #
  # #  8. lavstart # #
  # # # # # # # # # # #

  # if slotModel is NULL
  #   if lavpartable$est not NULL and lavoptions$start == "default"
  #     if there are free variances with est==0 or there are NA's in est
  #       compute start column in lavpartable via lav_start
  #     else
  #       set start column in lavpartable equal to est column
  #   else
  #     compute start column via lav_start and
  #       check via lav_start_check_cov if demanded (lavoptions$check.start)

  samplestats.flag <- TRUE
  if (!is.null(lavoptions$samplestats) &&
      !lavoptions$samplestats) {
    samplestats.flag <- FALSE
  }

  if (is.null(slotModel)) {
    # check if we have provided a full parameter table as model = input
    if (!is.null(lavpartable$est) && is.character(lavoptions$start) &&
      lavoptions$start == "default") {
      if (lav_verbose()) {
        cat("lavstart           ...")
      }
      # check if all 'est' values look ok
      # this is not the case, eg, if partables have been merged eg, as
      # in semTools' auxiliary() function

      # check for zero free variances and NA values
      if (is.null(lavpartable$lower)) {
        zero.idx <- which(lavpartable$free > 0L &
          lavpartable$op == "~~" &
          lavpartable$lhs == lavpartable$rhs &
          lavpartable$est == 0)
      } else {
        # ignore zero variances on the boundary; new in 0.6-21
        zero.idx <- which(lavpartable$free > 0L &
          lavpartable$op == "~~" &
          lavpartable$lhs == lavpartable$rhs &
          lavpartable$est == 0 &
          lavpartable$lower != 0)
      }

      if (length(zero.idx) > 0L || any(is.na(lavpartable$est))) {
        lavpartable$start <- lav_start(
          start.method = lavoptions$start,
          lavpartable = lavpartable,
          lavsamplestats = lavsamplestats,
          model.type = lavoptions$model.type,
          reflect = FALSE,
		      samplestats.flag = samplestats.flag,
          # order.lv.by  = lavoptions$rotation.args$order.lv.by,
          order.lv.by = "none"
          )
      } else {
        lavpartable$start <- lavpartable$est
      }

      # check for exogenous parameters: if the dataset changed, we must
      # update them! (new in 0.6-16)
      # ... or not? (not compatible with how we bootstrap under fixed.x = T)
      # we really need to think about this more carefully...
      #
      # if (any(lavpartable$exo == 1L)) {
      #     # FIXME: there should be an easier way just to
      #     # (re)initialize the the exogenous part of the model
      #     tmp <- lav_start(start.method   = "lavaan", # not "simple"
      #                                                 # if fixed.x = TRUE
      #                    lavpartable    = lavpartable,
      #                    lavsamplestats = lavsamplestats,
      #                    lavh1          = lavh1,
      #                    model.type     = lavoptions$model.type,
      #                    reflect      = FALSE,
      #                    #order.lv.by  = lavoptions$rotation.args$order.lv.by,
      #                    order.lv.by  = "none",
      #                    debug          = lav_debug())
      #     exo.idx <- which(lavpartable$exo == 1L)
      #     lavpartable$start[exo.idx] <- tmp[exo.idx]
      # }

      if (lav_verbose()) {
        cat(" done.\n")
      }
    } else {
      if (lav_verbose()) {
        cat("lavstart           ...")
      }
      start.values <- lav_start(
        start.method = lavoptions$start,
        lavpartable = lavpartable,
        lavsamplestats = lavsamplestats,
        lavh1 = lavh1,
        model.type = lavoptions$model.type,
        reflect = FALSE,
        samplestats.flag = samplestats.flag,
        # order.lv.by  = lavoptions$rotation.args$order.lv.by,
        order.lv.by = "none"
      )

      # sanity check
      if (!is.null(lavoptions$check.start) && lavoptions$check.start) {
        start.values <- lav_start_check_cov(
          lavpartable = lavpartable,
          start = start.values
        )
      }

      lavpartable$start <- start.values
      if (lav_verbose()) {
        cat(" done.\n")
      }
    }
  }

  lavpartable
}
