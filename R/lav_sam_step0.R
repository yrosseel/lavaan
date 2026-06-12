# STEP 0: process full model, without fitting
lav_sam_step0 <- function(cmd = "sem", model = NULL, data = NULL,
                          se = "twostep", sam_method = "local",
                          dotdotdot = NULL) {

  # create dotdotdot0 for dummy fit
  dotdotdot0 <- dotdotdot

  # parse model, so we can inspect a few features
  flat_model <- lavParseModelString(model)

  # remove do.fit option if present
  dotdotdot0$do.fit <- NULL

  if (sam_method %in% c("local", "fsr", "cfsr")) {
    dotdotdot0$sample.icov <- FALSE # if N < nvar
  }
  dotdotdot0$se                   <- "none"
  dotdotdot0$test                 <- "none"
  dotdotdot0$verbose              <- FALSE # no output for this 'dummy' FIT
  dotdotdot0$ceq.simple           <- TRUE # if not the default yet
  dotdotdot0$check.lv.interaction <- FALSE # we allow for it
  # note: setting cat.wls.w = FALSE (no weight matrix if categorical)
  # would break the computation of twostep standard errors
  if (se %in% c("local", "ij", "twostep.robust")) {
    dotdotdot0$sample.icov <- TRUE
    dotdotdot0$NACOV <- TRUE
    dotdotdot0$fixed.x <- FALSE
    dotdotdot0$ov_order <- "force.model" # avoid data ordering...
  }

  # any lv interaction terms?
  if (length(lav_object_vnames(flat_model, "lv.interaction")) > 0L) {
    dotdotdot0$meanstructure <- TRUE
  }

  # initial processing of the model, no fitting
  fit <- do.call(cmd,
    args = c(list(
      model  = flat_model,
      data   = data,
      do.fit = FALSE
    ), dotdotdot0)
  )

  # restore options

  # do.fit
  fit@Options$do.fit <- TRUE

  # sample.icov
  if (sam_method %in% c("local", "fsr", "cfsr")) {
    fit@Options$sample.icov <- TRUE
  }

  # se
  if (fit@Model@categorical && se == "twostep") {
    # FIXME!
    # should do this for global too, but we need the 'P' matrix, which
    # we only have for local (for now)
    if (sam_method == "local") {
      se <- "twostep.robust"
    }
  }
  fit@Options$se <- se

  # test
  if (!is.null(dotdotdot$test)) {
    fit@Options$test <- dotdotdot$test
  } else {
    fit@Options$test <- "standard"
  }

  # adjust parameter table:
  pt_1 <- fit@ParTable

  # check parameter table
  pt_1$est <- pt_1$se <- NULL
  # est equals ustart by default (except exo values)
  pt_1$est <- pt_1$ustart
  if (any(pt_1$exo > 0L)) {
    pt_1$est[pt_1$exo > 0L] <- pt_1$start[pt_1$exo > 0L]
  }

  # clear se values (needed here?) only for global approach to compute SE
  pt_1$se <- rep(as.numeric(NA), length(pt_1$lhs))
  pt_1$se[pt_1$free == 0L & !is.na(pt_1$ustart)] <- 0.0

  fit@ParTable <- pt_1


  fit
}
