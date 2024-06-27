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
  degree <- input$degree
  unbiased <- input$unbiased
  chisq <- if (input$chisq == "ml") lavaan::fitmeasures(object, "chisq") else lav_rls(object)
  df <- lavaan::fitmeasures(object, "df")
  ug <- lav_ugamma_no_groups(object, unbiased)
  lambdas_list <- Re(eigen(ug)$values)[seq(df)]

  lambdas <- lambdas_list
  list(
    test = lav_fmg_reconstruct_input(input),
    pvalue = if (name == "pEBA") {
    lav_fmg_peba(chisq, lambdas, degree)
  } else if (name == "EBA") {
    lav_fmg_eba(chisq, lambdas, degree)
  } else if (name == "pOLS") {
    lav_fmg_pols(chisq, lambdas, degree)
  },
  label = lav_fmg_reconstruct_label(input))
}

lav_fmg_parse_input <- \(string) {
  na_to_null <- \(x) if (is.na(x)) NULL else x
  default <- \(x) if (x != "") as.numeric(x) else 2

  tryCatch ({
    string <- tolower(string)
    splitted <- strsplit(string, "_")[[1]]
    unbiased <- FALSE
    chisq <- "ml"
    type <- na_to_null(splitted[1])

    if (length(splitted) == 3) {
      unbiased <- if (splitted[2] == "ug") TRUE else FALSE
      chisq <- splitted[3]
    } else if (length(splitted) == 2) {
      if (splitted[2] == "rls" || splitted[2] == "ml") {
        chisq <- splitted[2]
      } else if (splitted[2] == "ug") {
        unbiased <- TRUE
      }
    }

    out <- NULL
    name <- NULL
    if (startsWith(type, "peba")) {
      degree <- default(substring(type, 5))
      name <- "pEBA"
    } else if (startsWith(type, "eba")) {
      degree <- default(substring(type, 4))
      name <- "EBA"
    } else if (startsWith(type, "pols")) {
      degree <- default(substring(type, 5))
      name <- "pOLS"
    }

    if (is.null(degree) || (chisq != "ml" && chisq != "rls")) {
      return(NULL)
    }

    list(
      degree = degree,
      name = name,
      unbiased = unbiased,
      chisq = chisq
    )
  }, error = \(e){
    NULL
  })
}

lav_fmg_reconstruct_input <- \(input) {
  name <- input$name
  degree <- input$degree
  unbiased <- input$unbiased
  chisq <- input$chisq

  out <- paste0(name, degree)
  if (unbiased) out <- paste0(out, "_ug")
  if (chisq == "rls") out <- paste0(out, "_rls")
  out
}

lav_fmg_reconstruct_label <- \(input) {
  name <- tolower(input$name)
  degree <- input$degree
  unbiased <- input$unbiased
  chisq <- input$chisq

  out <- if (name == "peba") {
    paste0("Penalized EBA with ", degree, " blocks, ", if (unbiased) "unbiased" else "biased", " gamma, based on ",
    chisq, ".")
  } else if (name == "eba" && degree == 1) {
    paste0("Satorra-Bentler with ", if (unbiased) "unbiased" else "biased", " gamma, based on ",
           chisq, ".")
  } else if (name == "eba") {
    paste0("EBA with ", degree, " blocks, ", if (unbiased) "unbiased" else "biased", " gamma, based on ",
           chisq, ".")
  } else {
    paste0("Penalized OLS with weighting parameter ", degree, ", ", if (unbiased) "unbiased" else "biased", " gamma, based on ",
           chisq, ".")
  }
  out
}

lav_fmg_eba <- \(chisq, lambdas, j) {
  m <- length(lambdas)
  k <- ceiling(m / j)
  eig <- lambdas
  eig <- c(eig, rep(NA, k * j - length(eig)))
  dim(eig) <- c(k, j)
  eig_means <- colMeans(eig, na.rm = TRUE)
  repeated <- rep(eig_means, each = k)[seq(m)]
  lav_fmg_imhof(chisq, repeated)
}

lav_fmg_peba <- \(chisq, lambdas, j) {
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

lav_fmg_pols <- \(chisq, lambdas, gamma) {
  x <- seq_along(lambdas)
  beta1_hat <- 1 / gamma * stats::cov(x, lambdas) / stats::var(x)
  beta0_hat <- mean(lambdas) - beta1_hat * mean(x)
  lambda_hat <- pmax(beta0_hat + beta1_hat * x, 0)
  lav_fmg_imhof(chisq, lambda_hat)
}

lav_fmg_gamma_unbiased <- \(obj, gamma) {
  gamma_est_nt <- \(sigma) {
    lower_vec_indices <- \(n = 1L, diagonal = TRUE) {
      rows <- matrix(seq_len(n), n, n)
      cols <- matrix(seq_len(n), n, n, byrow = TRUE)
      if (diagonal) which(rows >= cols) else which(rows > cols)
    }

    upper_vec_indices <- \(n = 1L, diagonal = TRUE) {
      rows <- matrix(seq_len(n), n, n)
      cols <- matrix(seq_len(n), n, n, byrow = TRUE)
      tmp <- matrix(seq_len(n * n), n, n, byrow = TRUE)
      if (diagonal) tmp[rows >= cols] else tmp[rows > cols]
    }

    n <- ncol(sigma)
    lower <- lower_vec_indices(n)
    upper <- upper_vec_indices(n)
    y <- sigma %x% sigma
    out <- (y[lower, , drop = FALSE] + y[upper, , drop = FALSE]) / 2
    out[, lower, drop = FALSE] + out[, upper, drop = FALSE]
  }

  gamma_est_unbiased <- \(x, n = NULL, sigma = NULL, gamma_adf = NULL, gamma_nt = NULL) {
    vech <- \(x) x[row(x) >= col(x)]
    if (!missing(x)) n <- nrow(x)
    sigma <- if (is.null(sigma)) stats::cov(x) * (n - 1) / n else sigma
    gamma_adf <- if (is.null(gamma_adf)) gamma_est_adf(x) else gamma_adf
    gamma_nt <- if (is.null(gamma_nt)) gamma_est_nt(sigma) else gamma_nt
    gamma_rem <- tcrossprod(vech(sigma))
    mult <- n / ((n - 2) * (n - 3))
    mult * ((n - 1) * gamma_adf - (gamma_nt - 2 / (n - 1) * gamma_rem))
  }

  gamma_est_unbiased(
    n = lavaan::lavInspect(obj, "nobs"),
    sigma = obj@SampleStats@cov[[1]],
    gamma_adf = gamma,
    gamma_nt = NULL
  )
}

lav_ugamma_no_groups <- \(object, unbiased) {
  u <- lavaan::lavInspect(object, "U")
  bias <- object@Options$gamma.unbiased
  object@Options$gamma.unbiased <- FALSE
  gamma <- lavaan::lavInspect(object, "gamma")
  object@Options$gamma.unbiased <- bias
  if (unbiased) gamma <- lav_fmg_gamma_unbiased(object, gamma)
  u %*% gamma
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
