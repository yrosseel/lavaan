# SAM step 2: estimate structural part

lav_sam_step2 <- function(step1 = NULL, fit = NULL,
                          sam_method = "local", struc_args = list()) {
  lavoptions <- fit@Options
  lavpta <- fit@pta
  nlevels <- lavpta$nlevels
  pt_1 <- step1$PT

  # Gamma available?
  gamma_flag <- FALSE
  if (sam_method %in% c("local", "fsr", "cfsr") &&
      !is.null(step1$Gamma.eta[[1]])) {
    gamma_flag <- TRUE
  }

  lv_names <- unique(unlist(fit@pta$vnames$lv.regular))

  # adjust options
  lavoptions_pa <- lavoptions
  if (lavoptions_pa$se == "naive") {
    lavoptions_pa$se <- "standard"
  } else if (gamma_flag) {
    lavoptions_pa$se <- "robust.sem"
    lavoptions_pa$test <- "satorra.bentler"
  } else {
    # twostep or none -> none
    lavoptions_pa$se <- "none"
  }
  # lavoptions.PA$fixed.x <- TRUE # may be false if indicator is predictor
  if (!lavoptions_pa$conditional.x) {
    lavoptions_pa$fixed.x <- FALSE # until we fix this...
  }
  lavoptions_pa$categorical <- FALSE
  lavoptions_pa$.categorical <- FALSE
  lavoptions_pa$rotation <- "none"
  lavoptions_pa <- modifyList(lavoptions_pa, struc_args)

  if (gamma_flag) {
    lavoptions_pa$check.vcov <- FALSE # always non-pd
                                      # if interactions + fixed.x = FALSE
  }

  # override, no matter what
  lavoptions_pa$do.fit <- TRUE

  if (sam_method %in% c("local", "fsr", "cfsr")) {
    lavoptions_pa$missing <- "listwise"
    lavoptions_pa$sample.cov.rescale <- FALSE
    # lavoptions.PA$baseline <- FALSE
    # lavoptions.PA$h1 <- FALSE
    # lavoptions.PA$implied <- FALSE
    lavoptions_pa$loglik <- FALSE
  } else {
    lavoptions_pa$h1 <- FALSE
    # lavoptions.PA$implied <- FALSE
    lavoptions_pa$loglik <- FALSE
  }


  # construct PTS
  if (sam_method %in% c("local", "fsr", "cfsr")) {
    # extract structural part
    pts <- lav_pt_subset_structural_model(pt_1,
      add_idx = TRUE,
      add_exo_cov = TRUE,
      fixed_x = lavoptions_pa$fixed.x,
      conditional_x = lavoptions_pa$conditional.x,
      free_fixed_var = TRUE,
      meanstructure = lavoptions_pa$meanstructure
    )

    # any 'extra' parameters: not (free) in PT, but free in PTS (user == 3)
    #  - fixed.x in PT, but fixed.x = FALSE is PTS
    #  - fixed-to-zero intercepts in PT, but free in PTS
    #  - add.exo.cov: absent/fixed-to-zero in PT, but add/free in PTS
    extra_id <- which(pts$user == 3L)

    # remove est/se/start columns
    pts$est <- NULL
    pts$se <- NULL
    pts$start <- NULL

    if (nlevels > 1L) {
      pts$level <- NULL
      pts$group <- NULL
      pts$group <- pts$block
      nobs_1 <- fit@Data@Lp[[1]]$nclusters
    } else {
      nobs_1 <- fit@Data@nobs
    }

    # if meanstructure, 'free' user=0 intercepts?
    # if (lavoptions.PA$meanstructure) {
    #   extra.int.idx <- which(PTS$op == "~1" & PTS$user == 0L &
    #     PTS$free == 0L &
    #     PTS$exo == 0L) # needed?
    #   if (length(extra.int.idx) > 0L) {
    #     PTS$free[extra.int.idx] <- 1L
    #     PTS$ustart[extra.int.idx] <- as.numeric(NA)
    #     PTS$free[PTS$free > 0L] <-
    #       seq_len(length(PTS$free[PTS$free > 0L]))
    #     PTS$user[extra.int.idx] <- 3L
    #   }
    # } else {
    #   extra.int.idx <- integer(0L)
    # }
    # extra.id <- c(extra.id, extra.int.idx)

    reg_idx <- attr(pts, "idx")
    attr(pts, "idx") <- NULL
  } else {
    # global SAM

    # the measurement model parameters now become fixed ustart values
    pt_1$ustart[pt_1$free > 0] <- pt_1$est[pt_1$free > 0]

    reg_idx <- lav_pt_subset_structural_model(
      pt_1 = pt_1,
      idx_only = TRUE
    )

    # remove 'exogenous' factor variances (if any) from reg.idx
    lv_names_x <- lv_names[lv_names %in% unlist(lavpta$vnames$eqs.x) &
      !lv_names %in% unlist(lavpta$vnames$eqs.y)]
    if ((lavoptions_pa$fixed.x || lavoptions_pa$std.lv) &&
        length(lv_names_x) > 0L) {
      var_idx <- which(pt_1$lhs %in% lv_names_x &
        pt_1$op == "~~" &
        pt_1$lhs == pt_1$rhs)
      rm_idx <- which(reg_idx %in% var_idx)
      if (length(rm_idx) > 0L) {
        reg_idx <- reg_idx[-rm_idx]
      }
    }

    # adapt parameter table for structural part
    pts <- pt_1

    # remove constraints we don't need
    con_idx <- which(pts$op %in% c("==", "<", ">", ":="))
    if (length(con_idx) > 0L) {
      needed_idx <- which(con_idx %in% reg_idx)
      if (length(needed_idx) > 0L) {
        con_idx <- con_idx[-needed_idx]
      }
      if (length(con_idx) > 0L) {
        pts <- as.data.frame(pts, stringsAsFactors = FALSE)
        pts <- pts[-con_idx, ]
      }
    }
    pts$est <- NULL
    pts$se <- NULL

    # 'fix' step 1 parameters
    pts$free[!pts$id %in% reg_idx & pts$free > 0L] <- 0L

    # but free up residual variances if fixed (eg std.lv = TRUE) (new in 0.6-20)
    var_idx <- reg_idx[which(pt_1$free[reg_idx] == 0L &
                             pt_1$user[reg_idx] != 1L &
                             pt_1$op[reg_idx] == "~~")] # FIXME: more?
    pts$free[var_idx] <- max(pts$free) + seq_along(var_idx)

    # set 'ustart' values for free FIT.PA parameter to NA
    pts$ustart[pts$free > 0L] <- as.numeric(NA)

    pts <- lav_pt_complete(pts)

    extra_id <- integer(0L)
  } # global

  # fit structural model
  if (lav_verbose()) {
    cat("Fitting the structural part ... \n")
  }
  if (sam_method %in% c("local", "fsr", "cfsr")) {
    if (gamma_flag) {
      nacov <- step1$Gamma.eta
      # ov_order <- "data"
    } else {
      nacov <- NULL
      # ov_order <- "model"
    }
    fit_pa <- lavaan::lavaan(pts,
      sample.cov  = step1$VETA,
      sample.mean = step1$EETA,
      sample.nobs = nobs_1,
      NACOV       = nacov,
      slotOptions = lavoptions_pa,
      verbose     = FALSE
    )
  } else {
    fit_pa <- lavaan::lavaan(
      model = pts,
      slotData = fit@Data,
      slotSampleStats = fit@SampleStats,
      slotOptions = lavoptions_pa,
      verbose = FALSE
    )
  }
  if (lav_verbose()) {
    cat("Fitting the structural part ... done.\n")
  }

  # which parameters from PTS do we wish to fill in:
  # - all 'free' parameters
  # - :=, <, > (if any)
  # - and NOT element with user=3 (add.exo.cov = TRUE, extra.int.idx)
  pts_idx <- which((pts$free > 0L | (pts$op %in% c(":=", "<", ">"))) &
    !pts$user == 3L)

  # find corresponding rows in PT
  pts2 <- as.data.frame(pts, stringsAsFactors = FALSE)
  pt_idx <- lav_pt_map_id_p1_in_p2(pts2[pts_idx, ], pt_1,
    exclude_nonpar = FALSE
  )
  # fill in
  pt_1$est[pt_idx] <- fit_pa@ParTable$est[pts_idx]

  # create step2.free.idx
  p2_idx <- seq_along(pt_1$lhs) %in% pt_idx & pt_1$free > 0 # no def!
  step2_free_idx <- step1$PT.free[p2_idx]

  # add 'step' column in PT
  pt_1$step <- rep(1L, length(pt_1$lhs))
  pt_1$step[seq_along(pt_1$lhs) %in% reg_idx] <- 2L

  step2 <- list(
    FIT.PA = fit_pa, PT = pt_1, reg.idx = reg_idx,
    step2.free.idx = step2_free_idx, extra.id = extra_id,
    pt.idx = pt_idx, pts.idx = pts_idx
  )

  step2
}
