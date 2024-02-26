lav_lavaan_step04_partable <- function(slotParTable, model, flat.model, lavoptions, lavdata,
                                        constraints) {
  # # # # # # # # # # # # 
  # #  4. lavpartable # # 
  # # # # # # # # # # # # 
  if (!is.null(slotParTable)) {
    lavpartable <- slotParTable
  } else if (is.character(model) ||
             inherits(model, "formula")) {
    if (lavoptions$verbose) {
      cat("lavpartable        ...")
    }
    # check flat.model before we proceed
    if (lavoptions$debug) {
      print(as.data.frame(flat.model))
    }
    # catch ~~ of fixed.x covariates if fixed.x = TRUE
    # --> done inside lavaanify!
    
    #if(lavoptions$fixed.x) {
    #    tmp <- lav_partable_vnames(flat.model, type = "ov.x",
    #                               ov.x.fatal = FALSE, warn = TRUE)
    #tmp <- try(vnames(flat.model, type = "ov.x", ov.x.fatal = TRUE),
    #           silent = TRUE)
    #if(inherits(tmp, "try-error")) {
    #    warning("lavaan WARNING: syntax contains parameters involving exogenous covariates;",
    #    "switching to fixed.x = FALSE")
    #    lavoptions$fixed.x <- FALSE
    #}
    #}
    #if(lavoptions$conditional.x) {
    #    tmp <- vnames(flat.model, type = "ov.x", ov.x.fatal = TRUE)
    #}
    
    if (lavoptions$estimator == "catML") {
      lavoptions$meanstructure <- FALSE
      DataOV <- lavdata@ov
      DataOV$type <- rep("numeric", length(DataOV$type))
    } else {
      DataOV <- lavdata@ov
    }
    
    lavpartable <-
      lavaanify(model            = flat.model,
                constraints      = constraints,
                varTable         = DataOV,
                ngroups          = lavdata@ngroups,
                
                meanstructure    = lavoptions$meanstructure,
                int.ov.free      = lavoptions$int.ov.free,
                int.lv.free      = lavoptions$int.lv.free,
                marker.int.zero  = lavoptions$marker.int.zero,
                orthogonal       = lavoptions$orthogonal,
                orthogonal.x     = lavoptions$orthogonal.x,
                orthogonal.y     = lavoptions$orthogonal.y,
                orthogonal.efa   = lavoptions$rotation.args$orthogonal,
                conditional.x    = lavoptions$conditional.x,
                fixed.x          = lavoptions$fixed.x,
                std.lv           = lavoptions$std.lv,
                correlation      = lavoptions$correlation,
                effect.coding    = lavoptions$effect.coding,
                ceq.simple       = lavoptions$ceq.simple,
                parameterization = lavoptions$parameterization,
                auto.fix.first   = lavoptions$auto.fix.first,
                auto.fix.single  = lavoptions$auto.fix.single,
                auto.var         = lavoptions$auto.var,
                auto.cov.lv.x    = lavoptions$auto.cov.lv.x,
                auto.cov.y       = lavoptions$auto.cov.y,
                auto.th          = lavoptions$auto.th,
                auto.delta       = lavoptions$auto.delta,
                auto.efa         = lavoptions$auto.efa,
                group.equal      = lavoptions$group.equal,
                group.partial    = lavoptions$group.partial,
                group.w.free     = lavoptions$group.w.free,
                debug            = lavoptions$debug,
                warn             = lavoptions$warn,
                
                as.data.frame.   = FALSE)
    
    if (lavoptions$verbose) {
      cat(" done.\n")
    }
  } else if (inherits(model, "lavaan")) {
    lavpartable <- as.list(parTable(model))
  } else if (is.list(model)) {
    # we already checked this when creating flat.model
    # but we may need to complete it
    lavpartable <- as.list(flat.model) # in case model is a data.frame
    # complete table
    lavpartable <- as.list(lav_partable_complete(lavpartable))
  } else {
    stop("lavaan ERROR: model [type = ", class(model),
         "] is not of type character or list")
  }
  if (lavoptions$debug) {
    print(as.data.frame(lavpartable))
  }
  
  # at this point, we should check if the partable is complete
  # or not; this is especially relevant if the lavaan() function
  # was used, but the user has forgotten some variances/intercepts...
  if (is.null(slotParTable)) {
    junk <- lav_partable_check(lavpartable,
                               categorical = lavoptions$.categorical,
                               warn = TRUE)
  }
  
  # for EM only (for now), force fixed-to-zero (residual) variances
  # to be slightly larger than zero
  if (lavoptions$optim.method == "em") {
    zero.var.idx <- which(lavpartable$op == "~~" &
                            lavpartable$lhs == lavpartable$rhs &
                            lavpartable$free == 0L &
                            lavpartable$ustart == 0)
    if (length(zero.var.idx) > 0L) {
      lavpartable$ustart[zero.var.idx] <- lavoptions$em.zerovar.offset
    }
  }
  return (list(
    lavoptions = lavoptions,
    lavpartable = lavpartable
  ))
}

lav_lavaan_step04_pta <- function(lavpartable, lavoptions) {
  # # # # # # # # # # # # # # # # #
  # #  4b. parameter attributes # #
  # # # # # # # # # # # # # # # # #
  if (lavoptions$verbose) {
    cat("lavpta             ...")
  }
  lavpta <- lav_partable_attributes(lavpartable)
  if (lavoptions$verbose) {
    cat(" done.\n")
  }
  return (lavpta)
}

