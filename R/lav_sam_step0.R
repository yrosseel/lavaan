# STEP 0: process full model, without fitting
lav_sam_step0 <- function(cmd = "sem", model = NULL, data = NULL,
                          se = "twostep", sam.method = "local",
                          dotdotdot = NULL) {

  # create dotdotdot0 for dummy fit
  dotdotdot0 <- dotdotdot

  # parse model, so we can inspect a few features
  flat.model <- lavParseModelString(model)

  # remove do.fit option if present
  dotdotdot0$do.fit <- NULL

  #
  if (sam.method %in% c("local", "fsr")) {
    dotdotdot0$sample.icov <- FALSE # if N < nvar
  }
  dotdotdot0$se                   <- "none"
  dotdotdot0$test                 <- "none"
  dotdotdot0$verbose              <- FALSE # no output for this 'dummy' FIT
  dotdotdot0$conditional.x        <- FALSE
  dotdotdot0$fixed.x              <- TRUE
  dotdotdot0$ceq.simple           <- TRUE # if not the default yet
  dotdotdot0$check.lv.interaction <- FALSE # we allow for it
  # dotdotdot0$cat.wls.w            <- FALSE # no weight matrix if categorical
  # note: this break the computation of twostep standard errors...

  # any lv interaction terms?
  if (length(lavNames(flat.model, "lv.interaction")) > 0L) {
    dotdotdot0$meanstructure   <- TRUE
    dotdotdot0$marker.int.zero <- TRUE
  }

  # initial processing of the model, no fitting
  FIT <- do.call(cmd,
    args = c(list(
      model  = flat.model,
      data   = data,
      do.fit = FALSE
    ), dotdotdot0)
  )

  # restore options

  # do.fit
  FIT@Options$do.fit <- TRUE
  # FIT@Options$cat.wls.w <- TRUE

  # sample.icov
  if (sam.method %in% c("local", "fsr")) {
    FIT@Options$sample.icov <- TRUE
  }

  # se
  FIT@Options$se <- se

  # test
  if (!is.null(dotdotdot$test)) {
    FIT@Options$test <- dotdotdot$test
  } else {
    FIT@Options$test <- "standard"
  }

  # verbose
  if (!is.null(dotdotdot$verbose)) {
    FIT@Options$verbose <- dotdotdot$verbose
  }

  # adjust parameter table:
  PT <- FIT@ParTable

  # check parameter table
  PT$est <- PT$se <- NULL
  # est equals ustart by default (except exo values)
  PT$est <- PT$ustart
  if (any(PT$exo > 0L)) {
    PT$est[PT$exo > 0L] <- PT$start[PT$exo > 0L]
  }

  # clear se values (needed here?) only for global approach to compute SE
  PT$se <- rep(as.numeric(NA), length(PT$lhs))
  PT$se[PT$free == 0L & !is.na(PT$ustart)] <- 0.0

  FIT@ParTable <- PT


  FIT
}
