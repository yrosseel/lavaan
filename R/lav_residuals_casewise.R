# casewise residuals

lav_residuals_casewise <- function(object, labels = labels) {
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

  G <- object@Data@ngroups
  ov.names <- object@Data@ov.names

  X <- object@Data@X
  if (object@Model@categorical) {
    # add 'eXo' columns to X
    X <- lapply(seq_len(object@Data@ngroups), function(g) {
      ret <- cbind(X[[g]], object@Data@eXo[[g]])
      ret
    })
  }
  M <- lav_predict_yhat(object)
  # Note: if M has already class lavaan.matrix, print goes crazy
  # with Error: C stack usage is too close to the limit
  OUT <- lapply(seq_len(G), function(x) {
    out <- X[[x]] - M[[x]]
    class(out) <- c("lavaan.matrix", "matrix")
    out
  })

  if (labels) {
    for (g in 1:G) {
      colnames(OUT[[g]]) <- object@pta$vnames$ov[[g]]
    }
  }

  if (G == 1) {
    OUT <- OUT[[1]]
  } else {
    names(OUT) <- unlist(object@Data@group.label)
  }

  OUT
}
