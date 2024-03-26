# place-holder for MIIV estimation

# internal function to be used inside lav_optim_noniter
# return 'x', the estimated vector of free parameters
lav_sem_miiv_internal <- function(lavmodel = NULL, lavsamplestats = NULL,
                                  lavpartable = NULL,
                                  lavdata = NULL, lavoptions = NULL) {
  # this is the entry-point for MIIV estimation
  lav_msg_note(gettext("** Estimator MIIV is still under development! **\n"))

  # return error (for now)
  x <- as.numeric(NA)
  class(x) <- "try-error"
  x
}
