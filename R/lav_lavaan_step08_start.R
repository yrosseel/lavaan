lav_lavaan_step08_start <- function(slotModel, lavoptions, lavpartable, lavsamplestats, lavh1) {
  # # # # # # # # # # #
  # #  8. lavstart # # 
  # # # # # # # # # # #
  if (is.null(slotModel)) {
    # check if we have provided a full parameter table as model = input
    if (!is.null(lavpartable$est) && is.character(lavoptions$start) &&
        lavoptions$start == "default") {
      if (lavoptions$verbose) {
        cat("lavstart           ...")
      }
      # check if all 'est' values look ok
      # this is not the case, eg, if partables have been merged eg, as
      # in semTools' auxiliary() function
      
      # check for zero free variances and NA values
      zero.idx <- which(lavpartable$free > 0L &
                          lavpartable$op == "~~" &
                          lavpartable$lhs == lavpartable$rhs &
                          lavpartable$est == 0)
      
      if (length(zero.idx) > 0L || any(is.na(lavpartable$est))) {
        lavpartable$start <- lav_start(start.method = lavoptions$start,
                                       lavpartable     = lavpartable,
                                       lavsamplestats  = lavsamplestats,
                                       model.type   = lavoptions$model.type,
                                       reflect      = FALSE,
                                       #order.lv.by  = lavoptions$rotation.args$order.lv.by,
                                       order.lv.by  = "none",
                                       mimic        = lavoptions$mimic,
                                       debug        = lavoptions$debug)
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
      #                    mimic          = lavoptions$mimic,
      #                    debug          = lavoptions$debug)
      #     exo.idx <- which(lavpartable$exo == 1L)
      #     lavpartable$start[exo.idx] <- tmp[exo.idx]
      # }
      
      if (lavoptions$verbose) {
        cat(" done.\n")
      }
    } else {
      if (lavoptions$verbose) {
        cat("lavstart           ...")
      }
      start.values <- lav_start(start.method   = lavoptions$start,
                                lavpartable    = lavpartable,
                                lavsamplestats = lavsamplestats,
                                lavh1          = lavh1,
                                model.type     = lavoptions$model.type,
                                reflect      = FALSE,
                                #order.lv.by  = lavoptions$rotation.args$order.lv.by,
                                order.lv.by  = "none",
                                mimic          = lavoptions$mimic,
                                debug          = lavoptions$debug)
      
      # sanity check
      if (!is.null(lavoptions$check.start) && lavoptions$check.start) {
        start.values <- lav_start_check_cov(lavpartable = lavpartable,
                                            start       = start.values)
      }
      
      lavpartable$start <- start.values
      if (lavoptions$verbose) {
        cat(" done.\n")
      }
    }
  }
  return(lavpartable)
}

# 8b. EFA -- change user == 7 elements if columns
#     have been reordered
#        if (!is.null(lavpartable$efa) && any(lavpartable$user == 7L) &&
#            lavoptions$rotation != "none") {
#
#            # 7 to free
#            idx <- which(lavpartable$user == 7 &
#                         lavpartable$op == "=~" &
#                         abs(lavpartable$start) > sqrt(.Machine$double.eps))
#            if (length(idx) > 0L) {
#                lavpartable$user[idx] <- 1L
#                lavpartable$free[idx] <- 1L
#                lavpartable$ustart[idx] <- as.numeric(NA)
#            }
#
#            # free to 7
#            idx <- which(lavpartable$user != 7  &
#                         lavpartable$op == "=~" &
#                         lavpartable$free > 0L &
#                         nchar(lavpartable$efa) > 0L &
#                         abs(lavpartable$start) < sqrt(.Machine$double.eps))
#            if (length(idx) > 0L) {
#                lavpartable$user[idx] <- 7L
#                lavpartable$free[idx] <- 0L
#                lavpartable$ustart[idx] <- as.numeric(0)
#            }
#
#            # recount free parameters
#            idx <- which(lavpartable$free > 0L)
#            if (length(idx) > 0L) {
#                lavpartable$free[idx] <- seq_len(length(idx))
#            }
#        } # EFA
#
#