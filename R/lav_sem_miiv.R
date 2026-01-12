# place-holder for MIIV estimation

# internal function to be used inside lav_optim_noniter
# return 'x', the estimated vector of free parameters
lav_sem_miiv_internal <- function(lavmodel = NULL, lavsamplestats = NULL,
                                  lavpartable = NULL,
                                  lavdata = NULL, lavoptions = NULL) {
  # this is the entry-point for MIIV estimation
  lav_msg_note(gettext("** Estimator MIIV is still under development! **\n"))

  lavpta <- lav_partable_attributes(lavpartable)
  nblocks <- lavmodel@nblocks

  # find ALL model-implied instrumental variables (miivs) per equation
  eqs <- lav_model_find_iv(lavmodel = lavmodel, lavpta = lavpta)

  # 2SLS per equation (for now)
  x <- lav_model_get_parameters(lavmodel) * as.numeric(NA)
  for(b in seq_len(nblocks)) {
    ov.names <- lavpta$vnames$ov[[b]]
    XY <- lavdata@X[[b]]
    for (j in seq_along(eqs[[b]])) {
      eq <- eqs[[b]][[j]]

      # variable indices
      y.idx <- match(eq$lhs_new, ov.names)
      x.idx <- match(eq$rhs_new, ov.names)
      i.idx <- match(eq$miiv, ov.names)

      # 0. check
      if (length(eq$miiv) < length(eq$rhs)) {
        next
      }

      # 1. regress x on instruments
      fit_x_on_z <- lm.fit(x = cbind(1, XY[, i.idx, drop = FALSE]),
                           y = XY[, x.idx, drop = FALSE])
      xhat <- as.matrix(fit_x_on_z$fitted.values)

      # 2. regress y on xhat
      fit_y_on_xres <- lm.fit(x = cbind(1, xhat),
                              y = XY[, y.idx, drop = FALSE])
      int <- fit_y_on_xres$coefficients[1]
      coefs <- fit_y_on_xres$coefficients[-1]

      # 3. fill estimates in x
      x[lavpartable$free[eq$pt]] <- coefs
    } # eqs
  } # nblocks

  x
}
