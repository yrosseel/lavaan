# tools for the multivariate Bernoulli distribution
#
# see:
#
# Maydeu-Olivares & Joe (2005). Limited- and Full-Information Estimation and
# Goodness-of-Fit Testing in 2^n Contingency Tables: A Unified Framework.
# Journal of the American Statistical Association, 100, 1009--1020.

# YR. 15 April 2014 -- first version

# compute higher-order joint moments (Teugels 1991)
# prop must be an array, with dim = rep(2L, nitems)
lav_tables_mvb_prop_pidot <- function(prop, n_order = nitems) {
  # number of items/dimensions
  nitems <- length(dim(prop))

  # compute 'pi dot' up to order = n_order
  pidot <- unlist(
    lapply(1:n_order, function(order1) {
      idx <- utils::combn(1:nitems, order1)
      tmp <- apply(idx, 2L, function(idx1) {
        as.numeric(apply(prop, idx1, sum))[1L]
      })
      tmp
    })
  )

  pidot
}

# compute 'T' matrix, so that pidot = T %*% prop
lav_tables_mvb_gett <- function(nitems = 3L, n_order = nitems,
                                l_rbind = FALSE) {
  # index matrix
  index <- array(1:(2 ^ nitems), dim = rep(2L, nitems))

  t_r <- lapply(1:n_order, function(order1) {
    idx <- utils::combn(1:nitems, order1)
    m_tt <- matrix(0L, ncol(idx), 2 ^ nitems)
    m_tt <- do.call(
      "rbind",
      lapply(seq_len(ncol(idx)), function(i) {
        l_true <- as.list(rep(TRUE, nitems))
        l_true[idx[, i]] <- 1L
        v_args <- c(list(index), l_true)
        t1 <- integer(2 ^ nitems)
        t1[as.vector(do.call("[", v_args))] <- 1L
        t1
      })
    )
    m_tt
  })

  if (l_rbind) {
    t_r <- do.call("rbind", t_r)
  }

  t_r
}

# simple test function to check that  pidot = T %*% prop
lav_tables_mvb_test <- function(nitems = 3L) {
  freq <- sample(5:50, 2^nitems, replace = TRUE)
  prop <- freq / sum(freq)
  TABLE <- array(freq, dim = rep(2, nitems))
  PROP <- array(prop, dim = rep(2, nitems))
  # note: freq is always as.numeric(TABLE)
  #       prop is always as.numeric(PROP)

  pidot <- lav_tables_mvb_prop_pidot(PROP)
  T.r <- lav_tables_mvb_gett(nitems = nitems, n_order = nitems, l_rbind = TRUE)

  if (lav_verbose()) {
    out <- cbind(as.numeric(T.r %*% prop), pidot)
    colnames(out) <- c("T * prop", "pidot")
    print(out)
  }

  all.equal(pidot, as.numeric(T.r %*% prop))
}

# L_r test of Maydeu-Olivares & Joe (2005) eq (4)
lav_tables_mvb_Lr <- function(nitems = 0L,
                              obs.prop = NULL, est.prop = NULL, nobs = 0L,
                              n_order = 2L) {
  # recreate tables
  obs.PROP <- array(obs.prop, dim = rep(2L, nitems))
  est.PROP <- array(est.prop, dim = rep(2L, nitems))

  # compute {obs,est}.prop.dot
  obs.prop.dot <- lav_tables_mvb_prop_pidot(obs.PROP, n_order = n_order)
  est.prop.dot <- lav_tables_mvb_prop_pidot(est.PROP, n_order = n_order)

  # compute T.r
  T.r <- lav_tables_mvb_gett(nitems = nitems, n_order = n_order, l_rbind = TRUE)

  # compute GAMMA based on est.prop
  GAMMA <- diag(est.prop) - tcrossprod(est.prop)

  # compute XI
  XI <- T.r %*% GAMMA %*% t(T.r)

  # compute Lr
  diff.dot <- obs.prop.dot - est.prop.dot
  Lr <- as.numeric(nobs * t(diff.dot) %*% solve(XI) %*% diff.dot)
  df <- 2^nitems - 1L
  p.value <- 1 - pchisq(Lr, df = df)

  # return list
  list(Lr = Lr, df = df, p.value = p.value)
}
