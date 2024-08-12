lav_rls <- \(object) {
  s_inv <- chol2inv(chol(lavaan::lavInspect(object, "sigma.hat")))
  residuals <- lavaan::lavInspect(object, "residuals")
  mm <- residuals[[1]] %*% s_inv
  n <- lavaan::lavInspect(object, "nobs")
  .5 * n * sum(mm * t(mm))
}

lav_test_fmg <- \(object, input) {
  make_chisqs <- \(chisq) {
    if (chisq == "ml") lavaan::fitmeasures(object, "chisq") else rls(object)
  }

  name <- input$name
  param <- input$param
  #unbiased <- input$unbiased
  #chisq <- if (input$chisq == "ml") lavaan::fitmeasures(object, "chisq") else lav_rls(object)
  chisq <- fitmeasures(object, "chisq")
  df <- fitmeasures(object, "df")
  ug <- lav_object_gamma(object)[[1]] # Only for one group.
  #ug <- lav_ugamma_no_groups(object, unbiased)

  lambdas_list <- Re(eigen(ug)$values)[seq(df)]

  lambdas <- lambdas_list
  list(
    test = lav_fmg_reconstruct_label(input),
    pvalue = if (name == "fmg") {
    lav_fmg(chisq, lambdas, param)
  } else if (name == "EBA") {
    lav_fmg_eba(chisq, lambdas, param)
  } else if (name == "fmgols") {
    lav_fmgols(chisq, lambdas, param)
  },
  label = lav_fmg_reconstruct_label(input))
}

lav_fmg_parse_input <- \(string) {
  na_to_null <- \(x) if (is.na(x)) NULL else x
  default <- \(x) if (x != "") as.numeric(x) else 4

  tryCatch ({
    string <- tolower(string)
    splitted <- strsplit(string, "_")[[1]]
    name <- na_to_null(splitted[1])
    out <- NULL

    if (name != "fmg" && name != "fmgols") {
      return(NULL)
    }

    param <- if (length(splitted) == 2) {
      as.numeric(splitted[2])
    } else if (name == "fmg") {
      4
    } else {
      0.5
    }

    list(
      param = param,
      name = name
    )
  }, error = \(e){
    cat(e)
  })
}

lav_fmg_reconstruct_label <- \(input) paste0(input$name, "_", input$param)

lav_fmg_construct_label <- \(input) {
  if (input$name == "fmg") {
    paste0("FMG with ", input$param, " blocks.")
  } else {
    paste0("FMGOLS with weighting parameter ", input$param, ".")
  }
}

lav_fmg <- \(chisq, lambdas, j) {
  m <- length(lambdas)
  k <- ceiling(m / j)
  eig <- lambdas
  eig <- c(eig, rep(NA, k * j - length(eig)))
  dim(eig) <- c(k, j)
  eig_means <- colMeans(eig, na.rm = TRUE)
  eig_mean <- mean(lambdas)
  repeated <- rep(eig_means, each = k)[seq(m)]
  lav_fmg_imhof(chisq, (repeated + eig_mean) / 2)
}

lav_fmgols <- \(chisq, lambdas, gamma) {
  x <- seq_along(lambdas)
  beta1_hat <- 1 / gamma * stats::cov(x, lambdas) / stats::var(x)
  beta0_hat <- mean(lambdas) - beta1_hat * mean(x)
  lambda_hat <- pmax(beta0_hat + beta1_hat * x, 0)
  lav_fmg_imhof(chisq, lambda_hat)
}

lav_fmg_imhof <- \(x, lambda) {
  theta <- \(u, x, lambda) 0.5 * (sum(atan(lambda * u)) - x * u)
  rho <- \(u, lambda) exp(1 / 4 * sum(log(1 + lambda^2 * u^2)))
  integrand <- Vectorize(\(u) {
    sin(theta(u, x, lambda)) / (u * rho(u, lambda))
  })
  z <- integrate(integrand, lower = 0, upper = Inf)$value
  0.5 + z / pi
}


#' Calculate non-nested ugamma for multiple groups.
#' @param object A `lavaan` object.
#' @param unbiased If `TRUE`, uses the unbiased gamma estimate.
#' @keywords internal
#' @return Ugamma for non-nested object.
ugamma_non_nested <- \(object) {
  lavmodel <- object@Model

  ceq_idx <- attr(lavmodel@con.jac, "ceq.idx")
  if (length(ceq_idx) > 0L) {
    stop("Testing of models with groups and equality constraints not supported.")
  }

  test <- list()
  lavsamplestats <- object@SampleStats
  lavmodel <- object@Model
  lavoptions <- object@Options
  lavimplied <- object@implied
  lavdata <- object@Data
  test$standard <- object@test[[1]]

  if (test$standard$df == 0L || test$standard$df < 0) {
    stop("Df must be > 0.")
  }

  e <- lavaan:::lav_model_information(
    lavmodel = lavmodel,
    lavimplied = lavimplied,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    extra = TRUE
  )

  delta <- attr(e, "Delta")
  wls_v <- attr(e, "WLS.V")

  gamma <- lavaan:::lav_object_gamma(object)
  if (is.null(gamma[[1]])) {
    gamma <- lapply(lavaan::lavInspect(object, "gamma"), \(x) {
      class(x) <- "matrix"
      x
    })
  }

  gamma_global <- as.matrix(Matrix::bdiag(gamma))
  delta_global <- do.call(rbind, delta)
  v_global <- as.matrix(Matrix::bdiag(wls_v))
  x <- v_global %*% delta_global
  u_global <- v_global - crossprod(t(x), solve(t(delta_global) %*% x, t(x)))
  u_global %*% gamma_global
}

#' Calculate nested ugamma.
#'
#' This can also be used with restrictions.
#'
#' @param m0,m1 Two nested `lavaan` objects.
#' @param a The `A` matrix. If if `NULL`, gets calculated by
#'    `lavaan:::lav_test_diff_A` with `method = method`.
#' @param method Method passed to `lavaan:::lav_test_diff_A`.
#' @param unbiased If `TRUE`, uses the unbiased gamma estimate.
#' @keywords internal
#' @return Ugamma for nested object.
lav_ugamma_nested <- \(m0, m1, a = NULL, method = "delta", unbiased = FALSE) {
  m0@Options$gamma.unbiased <- unbiased
  m1@Options$gamma.unbiased <- unbiased

  # extract information from m1 and m2
  t1 <- m1@test[[1]]$stat
  r1 <- m1@test[[1]]$df

  t0 <- m0@test[[1]]$stat
  r0 <- m0@test[[1]]$df

  # m = difference between the df's
  m <- r0 - r1

  # check for identical df setting
  if (m == 0L) {
    return(list(
      T.delta = (t0 - t1), scaling.factor = as.numeric(NA),
      df.delta = m, a = as.numeric(NA), b = as.numeric(NA)
    ))
  }

  gamma <- lavaan:::lav_object_gamma(m0) # the same for m1 and m0

  if (is.null(gamma)) {
    stop("lavaan error: Can not compute gamma matrix; perhaps missing \"ml\"?")
  }

  wls_v <- lavaan::lavTech(m1, "WLS.V")
  pi <- lavaan::lavInspect(m1, "delta")

  p_inv <- lavaan::lavInspect(m1, what = "inverted.information")

  if (is.null(a)) {
    a <- lavaan:::lav_test_diff_A(m1, m0, method = method, reference = "H1")
    if (m1@Model@eq.constraints) {
      a <- a %*% t(m1@Model@eq.constraints.K)
    }
  }

  paapaap <- p_inv %*% t(a) %*% MASS::ginv(a %*% p_inv %*% t(a)) %*% a %*% p_inv

  # compute scaling factor
  fg <- unlist(m1@SampleStats@nobs) / m1@SampleStats@ntotal

  # We need the global gamma, cf. eq.~(10) in Satorra (2000).
  gamma_rescaled <- gamma
  for (i in (seq_along(gamma))) {
    gamma_rescaled[[i]] <- fg[i] * gamma_rescaled[[i]]
  }
  gamma_global <- as.matrix(Matrix::bdiag(gamma_rescaled))
  # Also the global V:
  v_global <- as.matrix(Matrix::bdiag(wls_v))
  pi_global <- do.call(rbind, pi)
  # U global version, eq.~(22) in Satorra (2000).
  u_global <- v_global %*% pi_global %*% paapaap %*% t(pi_global) %*% v_global
  return(u_global %*% gamma_global)
}
