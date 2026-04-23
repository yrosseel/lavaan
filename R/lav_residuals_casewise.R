# casewise residuals

lav_residuals_casewise <- function(object, labels = labels) {

  # check object
  object <- lav_object_check_version(object)

  # check if we have full data
  if (object@Data@data.type != "full") {
    lav_msg_stop(gettext("casewise residuals not available if sample statistics
                         were used for fitting the model"))
  }
  # check if we have categorical data
  if (object@Model@categorical) {
    lav_msg_stop(gettext(
      "casewise residuals not available if data is categorical"))
  }

  g_1 <- object@Data@ngroups

  x_1 <- object@Data@X
  if (object@Model@categorical) {
    # add 'eXo' columns to X
    x_1 <- lapply(seq_len(object@Data@ngroups), function(g) {
      ret <- cbind(x_1[[g]], object@Data@eXo[[g]])
      ret
    })
  }
  m <- lav_predict_yhat(object)
  # Note: if M has already class lavaan.matrix, print goes crazy
  # with Error: C stack usage is too close to the limit
  out_1 <- lapply(seq_len(g_1), function(x) {
    out <- x_1[[x]] - m[[x]]
    class(out) <- c("lavaan.matrix", "matrix")
    out
  })

  if (labels) {
    for (g in 1:g_1) {
      colnames(out_1[[g]]) <- object@pta$vnames$ov[[g]]
    }
  }

  if (g_1 == 1) {
    out_1 <- out_1[[1]]
  } else {
    names(out_1) <- unlist(object@Data@group.label)
  }

  out_1
}
