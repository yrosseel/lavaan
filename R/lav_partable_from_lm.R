# build a bare-bones parameter table from a fitted lm object
#
# YR: this function was broken since Mar 3, 2017, but nobody noticed this!
#     fixed again Apr 29, 2025.

lav_partable_from_lm <- function(object, est = FALSE, label = FALSE,
                                 as.data.frame. = FALSE) {           # nolint
  # sanity check
  if (!inherits(object, "lm")) {
    lav_msg_stop(gettext("object must be of class lm"))
  }

  object_terms <- terms(object)

  response_index <- attr(object_terms, "response")
  var_names <- as.character(attr(object_terms, "variables"))[-1]
  response_name <- var_names[response_index]

  pred_coef <- coef(object)
  pred_names <- names(pred_coef)

  lhs <- rep(response_name, length(pred_names))
  op <- rep("~", length(pred_names))
  rhs <- pred_names

  # intercept?
  if (attr(object_terms, "intercept")) {
    int_idx <- which(rhs == "(Intercept)")
    op[int_idx] <- "~1"
    rhs[int_idx] <- ""
  }

  # always add residual variance?
  # lhs <- c(lhs, responseName)
  # op <- c(op, "~~")
  # rhs <- c(rhs, responseName)

  # construct minimal partable
  partable <- list(lhs = lhs, op = op, rhs = rhs)

  # include 'est' column?
  if (est) {
    # partable$est <- c(as.numeric(predCoef),
    #                  sum(resid(object)^2) / object$df.residual)
    partable$est <- as.numeric(pred_coef)
  }

  # include 'label' column?
  if (label) {
    # partable$label <- c(predNames, responseName)
    partable$label <- pred_names

    # convert all ':' to '.'
    partable$label <- gsub("[:()]", ".", partable$label)
  }

  # convert to data.frame?
  if (as.data.frame.) {
    partable <- as.data.frame(partable, stringsAsFactors = FALSE)
  }

  partable
}
