lav_samplestats_step2 <- function(UNI = NULL,
                                  wt = NULL,
                                  ov.names = NULL, # error message only
                                  # polychoric and empty cells
                                  zero.add = c(0.5, 0.0),
                                  zero.keep.margins = TRUE,
                                  zero.cell.warn = FALSE,
                                  # keep track of tables with zero cells?
                                  zero.cell.tables = TRUE) {
  nvar <- length(UNI)
  COR <- diag(nvar)

  if (zero.cell.tables) {
    zero.var1 <- character(0L)
    zero.var2 <- character(0L)
  }

  # one-by-one (for now)
  for (j in seq_len(nvar - 1L)) {
    for (i in (j + 1L):nvar) {
      if (is.null(UNI[[i]]$th.idx) &&
        is.null(UNI[[j]]$th.idx)) {
        rho <- lav_bvreg_cor_twostep_fit(
          fit.y1 = UNI[[i]], # linear
          fit.y2 = UNI[[j]], # linear
          wt = wt,
          Y1.name = ov.names[i],
          Y2.name = ov.names[j]
        )
        COR[i, j] <- COR[j, i] <- rho
      } else if (is.null(UNI[[i]]$th.idx) &&
        !is.null(UNI[[j]]$th.idx)) {
        # polyserial
        rho <- lav_bvmix_cor_twostep_fit(
          fit.y1 = UNI[[i]], # linear
          fit.y2 = UNI[[j]], # ordinal
          wt = wt,
          Y1.name = ov.names[i],
          Y2.name = ov.names[j]
        )
        COR[i, j] <- COR[j, i] <- rho
      } else if (is.null(UNI[[j]]$th.idx) &&
        !is.null(UNI[[i]]$th.idx)) {
        # polyserial
        rho <- lav_bvmix_cor_twostep_fit(
          fit.y1 = UNI[[j]], # linear
          fit.y2 = UNI[[i]], # ordinal
          wt = wt,
          Y1.name = ov.names[j],
          Y2.name = ov.names[i]
        )
        COR[i, j] <- COR[j, i] <- rho
      } else if (!is.null(UNI[[i]]$th.idx) &&
        !is.null(UNI[[j]]$th.idx)) {
        # polychoric correlation
        rho <- lav_bvord_cor_twostep_fit(
          fit.y1 = UNI[[j]], # ordinal
          fit.y2 = UNI[[i]], # ordinal
          wt = wt,
          zero.add = zero.add,
          zero.keep.margins = zero.keep.margins,
          zero.cell.warn = zero.cell.warn,
          zero.cell.flag = zero.cell.tables,
          Y1.name = ov.names[i],
          Y2.name = ov.names[j]
        )
        if (zero.cell.tables) {
          if (attr(rho, "zero.cell.flag")) {
            zero.var1 <- c(zero.var1, ov.names[j])
            zero.var2 <- c(zero.var2, ov.names[i])
          }
          attr(rho, "zero.cell.flag") <- NULL
        }
        COR[i, j] <- COR[j, i] <- rho
      }
      # check for near 1.0 correlations
      if (abs(COR[i, j]) > 0.99) {
        lav_msg_warn(gettextf(
          "correlation between variables %1$s and %2$s is (nearly) 1.0",
          ov.names[i], ov.names[j]))
      }
    }
  }

  # keep track of tables with zero cells
  if (zero.cell.tables) {
    zero.cell.tables <- cbind(zero.var1, zero.var2)
    attr(COR, "zero.cell.tables") <- zero.cell.tables
  }

  COR
}
