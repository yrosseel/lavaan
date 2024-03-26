# export go BUGS syntax

# we assume that N1, N2, ... are in data
lav2bugs <- function(partable, as.function. = FALSE) {
  # get parameter table attributes
  pta <- lav_partable_attributes(partable = partable)
  vnames <- pta$vnames
  nblocks <- pta$nblocks
  nvar <- pta$nvar
  nfac <- pta$nfac

  # sanity check
  partable <- lav2check(partable)

  # tabs
  t1 <- paste(rep(" ", 2L), collapse = "")
  t2 <- paste(rep(" ", 4L), collapse = "")
  t3 <- paste(rep(" ", 6L), collapse = "")
  t4 <- paste(rep(" ", 8L), collapse = "")

  # TXT header
  if (as.function.) {
    TXT <- paste("{\n", sep = "")
  } else {
    TXT <- paste("model {\n", sep = "")
  }

  # model for every i
  for (g in 1:nblocks) {
    ov.names <- vnames$ov[[g]]
    lv.names <- vnames$lv[[g]]
    yname <- paste("y", g, sep = "")
    if (nblocks > 1L) {
      TXT <- paste(TXT, t1,
        "# block ", g, "\n",
        sep = ""
      )
    } else {
      TXT <- paste(TXT, "\n")
    }
    TXT <- paste(TXT, t1,
      "for(i in 1:N", g, ") {\n",
      sep = ""
    )

    # ov.nox - all observed variables (except exogenous ones)
    ov.names.nox <- vnames$ov.nox[[g]]
    nov <- length(ov.names.nox)
    TXT <- paste(TXT, "\n", t2,
      "# ov.nox",
      sep = ""
    )
    for (i in 1:nov) {
      ov.idx <- match(ov.names.nox[i], ov.names)
      theta.free.idx <- which(partable$block == g &
        partable$op == "~~" &
        partable$lhs == partable$rhs &
        partable$lhs == ov.names.nox[i])
      if (length(theta.free.idx) != 1L) {
        lav_msg_stop(gettextf(
          "parameter for residual variance %s not found",
          ov.names.nox[i])
        )
      } else {
        theta.idx <- partable$free[theta.free.idx]
      }
      TXT <- paste(TXT, "\n", t2,
        yname, "[i,", ov.idx, "] ~ dnorm(mu", g, "[i,", ov.idx,
        "], itheta[", theta.idx, "])",
        sep = ""
      )
    }

    TXT <- paste(TXT, "\n", t2, sep = "")
    for (i in 1:nov) {
      ov.idx <- match(ov.names.nox[i], ov.names)
      TXT <- paste(TXT, "\n", t2,
        "mu", g, "[i,", ov.idx, "] <- ",
        sep = ""
      )

      # find rhs for this observed variable

      # 1. intercept?
      int.idx <- which(partable$block == g &
        partable$op == "~1" &
        partable$lhs == ov.names.nox[i])
      if (length(int.idx) == 1L) {
        # fixed or free?
        if (partable$free[int.idx] == 0L) {
          TXT <- paste(TXT,
            partable$ustart[int.idx],
            sep = ""
          )
        } else {
          TXT <- paste(TXT,
            "theta[", partable$free[int.idx], "]",
            sep = ""
          )
        }
      } else { # no intercept, say '0', so we always have rhs
        TXT <- paste(TXT, "0", sep = "")
      }

      # 2. factor loading?
      lam.idx <- which(partable$block == g &
        partable$op == "=~" &
        partable$rhs == ov.names.nox[i])
      for (j in lam.idx) {
        # fixed or free?
        if (partable$free[j] == 0L) {
          TXT <- paste(TXT, " + ",
            partable$ustart[j], "*eta", g, "[i,",
            match(partable$lhs[j], lv.names),
            "]",
            sep = ""
          )
        } else {
          TXT <- paste(TXT, " + ",
            "theta[", partable$free[j], "]*eta", g, "[i,",
            match(partable$lhs[j], lv.names),
            "]",
            sep = ""
          )
        }
      }

      # 3. regression?
      r.idx <- which(partable$block == g &
        partable$op == "~" &
        partable$lhs == ov.names.nox[i])
      for (j in r.idx) {
        # what is the rhs?
        rhs <- partable$rhs[j]
        if (rhs %in% lv.names) {
          RHS <- paste("eta", g, "[i,",
            match(rhs, lv.names), "]",
            sep = ""
          )
        } else if (rhs %in% vnames$ov[[g]]) {
          RHS <- paste("y", g, "[i,",
            match(rhs, ov.names), "]",
            sep = ""
          )
        }

        # fixed or free?
        if (partable$free[j] == 0L) {
          TXT <- paste(TXT, " + ",
            partable$ustart[j], "*", RHS,
            sep = ""
          )
        } else {
          TXT <- paste(TXT, " + ",
            "theta[", partable$free[j], "]*", RHS,
            sep = ""
          )
        }
      }
    }


    # lv.y
    # var(lv.y) = PSI (lisrel style)
    lv.y <- vnames$lv.y[[g]]
    if (length(lv.y) > 0L) {
      TXT <- paste(TXT, "\n\n", t2,
        "# lv.y",
        sep = ""
      )
      lv.y.idx <- match(lv.y, lv.names)
      ny <- length(lv.y.idx)
      for (j in 1:ny) {
        theta.free.idx <- which(partable$block == g &
          partable$op == "~~" &
          partable$lhs == partable$rhs &
          partable$lhs == lv.y[j])
        if (length(theta.free.idx) != 1L) {
          lav_msg_stop(gettextf(
            "parameter for residual variance %s not found",
            lv.y[j])
          )
        } else {
          theta.idx <- partable$free[theta.free.idx]
        }
        TXT <- paste(TXT, "\n", t2,
          # dnorm for now
          "eta", g, "[i,", lv.y.idx[j], "] ~ dnorm(mu.eta", g, "[i,",
          lv.y.idx[j], "], itheta[", theta.idx, "])",
          sep = ""
        )
      }
      for (j in 1:ny) {
        TXT <- paste(TXT, "\n", t2,
          # dnorm for now
          "mu.eta", g, "[i,", lv.y.idx[j], "] <- ",
          sep = ""
        )

        # lhs elements regression
        # 1. intercept?
        int.idx <- which(partable$block == g &
          partable$op == "~1" &
          partable$lhs == lv.y[j])
        if (length(int.idx) == 1L) {
          # fixed or free?
          if (partable$free[int.idx] == 0L) {
            TXT <- paste(TXT,
              partable$ustart[int.idx],
              sep = ""
            )
          } else {
            TXT <- paste(TXT,
              "theta[", partable$free[int.idx], "]",
              sep = ""
            )
          }
        } else { # no intercept, say '0', so we always have rhs
          TXT <- paste(TXT, "0", sep = "")
        }

        rhs.idx <- which(partable$block == g &
          partable$op == "~" &
          partable$lhs == lv.y[j])
        np <- length(rhs.idx)
        for (p in 1:np) {
          TXT <- paste(TXT, " + ",
            "theta[", partable$free[rhs.idx[p]],
            "]*eta", g, "[i,",
            match(partable$rhs[rhs.idx[p]], lv.names),
            "]",
            sep = ""
          )
        }
      }
    }

    # exogenous lv -- FIXME: we assume the lv.x array is continous
    # (eg 3,4,5, but NOT 3,5,6)
    # var(lv.x) = PHI (lisrel style)
    lv.x <- vnames$lv.x[[g]]
    if (length(lv.x) > 0L) {
      TXT <- paste(TXT, "\n\n", t2,
        "# lv.x",
        sep = ""
      )
      lv.x.idx <- match(lv.x, lv.names)
      nx <- length(lv.x.idx)
      TXT <- paste(TXT, "\n", t2,
        # dmnorm for now
        "eta", g, "[i,", min(lv.x.idx), ":", max(lv.x.idx),
        "] ~ dmnorm(mu.eta", g, "[i,", min(lv.x.idx), ":",
        max(lv.x.idx), "], iphi", g, "[1:", nx, ",1:", nx, "])",
        sep = ""
      )
      for (j in 1:nx) {
        TXT <- paste(TXT, "\n", t2,
          "mu.eta", g, "[i,", lv.x.idx[j], "] <- 0",
          sep = ""
        )
      }
    }


    # exogenous ov ??? (what to do here?)

    # end of this block
    TXT <- paste(TXT, "\n\n", t1,
      "} # end of block ", g, "\n",
      sep = ""
    )
  }

  # priors (both fixed and free)
  TXT <- paste(TXT, "\n", t1,
    "# Priors free parameters (univariate):",
    sep = ""
  )
  npt <- length(partable$lhs)
  for (i in seq_len(npt)) {
    if (partable$free[i] == 0L) next # skip non-free parameters
    lhs <- partable$lhs[i]
    op <- partable$op[i]
    rhs <- partable$rhs[i]
    free.idx <- partable$free[i]
    g <- partable$block[i]
    if (op == "=~") {
      # factor loading
      TXT <- paste(TXT, "\n", t1,
        "theta[", free.idx, "] ~ dnorm(0.8, 1)",
        sep = ""
      )
    } else if (op == "~") {
      # regression
      TXT <- paste(TXT, "\n", t1,
        "theta[", free.idx, "] ~ dnorm(0, 1)",
        sep = ""
      )
    } else if (op == "~~" && lhs == rhs) {
      # variance
      # 1. lv.x + lv.x (skip -> multivariate)
      # 2. lv.y + lv.y
      # 3. observed + observed
      # 4. else -> fix (upgrade to latent?)
      if (lhs %in% vnames$lv.x[[g]] && rhs %in% vnames$lv.x[[g]]) {
        # lv.x: move to multivariate... (dwish)
        next
      } else if (lhs %in% vnames$lv.y[[g]] && rhs %in% vnames$lv.y[[g]]) {
        # lv.y
        TXT <- paste(TXT, "\n", t1,
          "itheta[", free.idx, "] ~ dgamma(9, 4)",
          sep = ""
        )
        TXT <- paste(TXT, "\n", t1,
          "theta[", free.idx, "] <- 1/itheta[", free.idx, "]",
          sep = ""
        )
      } else if (lhs %in% vnames$ov[[g]] && rhs %in% vnames$ov[[g]]) {
        TXT <- paste(TXT, "\n", t1,
          "itheta[", free.idx, "] ~ dgamma(9, 4)",
          sep = ""
        )
        TXT <- paste(TXT, "\n", t1,
          "theta[", free.idx, "] <- 1/itheta[", free.idx, "]",
          sep = ""
        )
      } else {
        lav_msg_stop(gettextf("FIXME!! parameter %s", i))
      }
    } else if (op == "~~" && lhs != rhs) {
      # covariance
      # 1. lv.x + lv.x (skip -> multivariate)
      # 2. lv.y + lv.y
      # 3. observed + observed
      # 4. else -> fix (upgrade to latent?)
      if (lhs %in% vnames$lv.x[[g]] && rhs %in% vnames$lv.x[[g]]) {
        # exo lv covariance
        next
      } else if (lhs %in% vnames$lv.y[[g]] && rhs %in% vnames$lv.y[[g]]) {
        # lv.y
        lav_msg_stop(gettextf("FIXME!! parameter ", i))
      } else if (lhs %in% vnames$ov[[g]] && rhs %in% vnames$ov[[g]]) {
        TXT <- paste(TXT, "\n", t1,
          "itheta[", free.idx, "] ~ dgamma(9, 4)",
          sep = ""
        )
        TXT <- paste(TXT, "\n", t1,
          "theta[", free.idx, "] <- 1/itheta[", free.idx, "]",
          sep = ""
        )
      } else {
        lav_msg_stop(gettextf("FIXME!! parameter ", i))
      }
    } else if (op == "~1") {
      # intercept
      TXT <- paste(TXT, "\n", t1,
        "theta[", free.idx, "] ~ dnorm(0, 1)",
        sep = ""
      )
    } else {
      lav_msg_stop(gettextf("op not supported yet for parameter ", i))
    }
  }

  TXT <- paste(TXT, "\n\n", t1,
    "# Priors free parameters (multivariate):",
    sep = ""
  )
  for (g in 1:nblocks) {
    lv.phi.idx <- which(partable$block == g &
      partable$op == "~~" &
      partable$lhs %in% vnames$lv.x[[g]] &
      partable$rhs %in% vnames$lv.x[[g]])
    nx <- length(vnames$lv.x[[g]])
    if (length(nx) > 0L) {
      TXT <- paste(TXT, "\n", t1,
        "iphi", g, "[1:", nx, ",1:", nx, "] ~ dwish(R", g, "[1:",
        nx, ",1:", nx, "], 5)",
        sep = ""
      )
      TXT <- paste(TXT, "\n", t1,
        "phi", g, "[1:", nx, ",1:", nx, "] <- inverse(iphi", g, "[1:",
        nx, ",1:", nx, "])",
        sep = ""
      )
      for (idx in lv.phi.idx) {
        TXT <- paste(TXT, "\n", t1,
          "theta[", partable$free[idx], "] <- phi", g, "[",
          match(partable$lhs[idx], vnames$lv.x[[g]]), ",",
          match(partable$rhs[idx], vnames$lv.x[[g]]), "]",
          sep = ""
        )
      }
    }
  }

  # end of model
  TXT <- paste(TXT, "\n\n", "} # End of model\n", sep = "")

  # end of model
  if (as.function.) {
    out <- function() NULL
    formals(out) <- alist()
    body(out) <- parse(file = "", text = TXT)
  } else {
    out <- TXT
    class(out) <- c("lavaan.character", "character")
  }

  out
}
