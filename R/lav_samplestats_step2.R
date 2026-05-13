lav_samplestats_step2 <- function(uni = NULL,
                                  wt = NULL,
                                  ov_names = NULL, # error message only
                                  # polychoric and empty cells
                                  zero_add = c(0.5, 0.0),
                                  zero_keep_margins = TRUE,
                                  zero_cell_warn = FALSE,
                                  # keep track of tables with zero cells?
                                  zero_cell_tables = TRUE) {
  nvar <- length(uni)
  cor_1 <- diag(nvar)

  if (zero_cell_tables) {
    zero_var1 <- character(0L)
    zero_var2 <- character(0L)
  }

  # one-by-one (for now)
  for (j in seq_len(nvar - 1L)) {
    for (i in (j + 1L):nvar) {
      if (is.null(uni[[i]]$th_idx) &&
        is.null(uni[[j]]$th_idx)) {
        rho <- lav_bvreg_cor_twostep_fit(
          fit_y1 = uni[[i]], # linear
          fit_y2 = uni[[j]], # linear
          wt = wt,
          y1_name = ov_names[i],
          y2_name = ov_names[j]
        )
        cor_1[i, j] <- cor_1[j, i] <- rho
      } else if (is.null(uni[[i]]$th_idx) &&
        !is.null(uni[[j]]$th_idx)) {
        # polyserial
        rho <- lav_bvmix_cor_twostep_fit(
          fit_y1 = uni[[i]], # linear
          fit_y2 = uni[[j]], # ordinal
          wt = wt,
          y1_name = ov_names[i],
          y2_name = ov_names[j]
        )
        cor_1[i, j] <- cor_1[j, i] <- rho
      } else if (is.null(uni[[j]]$th_idx) &&
        !is.null(uni[[i]]$th_idx)) {
        # polyserial
        rho <- lav_bvmix_cor_twostep_fit(
          fit_y1 = uni[[j]], # linear
          fit_y2 = uni[[i]], # ordinal
          wt = wt,
          y1_name = ov_names[j],
          y2_name = ov_names[i]
        )
        cor_1[i, j] <- cor_1[j, i] <- rho
      } else if (!is.null(uni[[i]]$th_idx) &&
        !is.null(uni[[j]]$th_idx)) {
        # polychoric correlation
        rho <- lav_bvord_cor_twostep_fit(
          fit_y1 = uni[[j]], # ordinal
          fit_y2 = uni[[i]], # ordinal
          wt = wt,
          zero_add = zero_add,
          zero_keep_margins = zero_keep_margins,
          zero_cell_warn = zero_cell_warn,
          zero_cell_flag = zero_cell_tables,
          y1_name = ov_names[i],
          y2_name = ov_names[j]
        )
        if (zero_cell_tables) {
          if (attr(rho, "zero.cell.flag")) {
            zero_var1 <- c(zero_var1, ov_names[j])
            zero_var2 <- c(zero_var2, ov_names[i])
          }
          attr(rho, "zero.cell.flag") <- NULL
        }
        cor_1[i, j] <- cor_1[j, i] <- rho
      }
      # check for near 1.0 correlations
      if (abs(cor_1[i, j]) > 0.99) {
        lav_msg_warn(gettextf(
          "correlation between variables %1$s and %2$s is (nearly) 1.0",
          ov_names[i], ov_names[j]))
      }
    }
  }

  # keep track of tables with zero cells
  if (zero_cell_tables) {
    zero_cell_tables <- cbind(zero_var1, zero_var2)
    attr(cor_1, "zero.cell.tables") <- zero_cell_tables
  }

  cor_1
}
