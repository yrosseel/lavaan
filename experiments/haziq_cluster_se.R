# pak::pak("haziqj/lavaan.bingof")
library(lavaan.bingof)
library(tidyverse)

# Function to obtain clustered SE ----------------------------------------------
get_clust_se <- function(object, dat) {
  
  # Obtain scores
  Hinv <- lavaan.bingof:::get_sensitivity_inv_mat(object)  # H^{-1}
  Delta <- lavaan:::computeDelta(lavmodel = object@Model)[[1]]
  SC <- lavaan:::pml_deriv1(
    Sigma.hat = object@implied$cov[[1]],
    Mu.hat = object@implied$mean[[1]],
    TH = object@implied$th[[1]],
    th.idx = object@Model@th.idx[[1]],
    num.idx = object@Model@num.idx[[1]],
    X = object@Data@X[[1]],
    eXo = NULL,
    wt = NULL,
    PI = NULL,
    lavcache = object@Cache[[1]],
    missing = object@Data@missing,
    scores = TRUE,
    negative = FALSE
  )
  
  # Combine with weights information from dat
  SCwt <-
    dat |>
    select(school, wt) |>
    bind_cols(as_tibble(SC))
  
  # Get the clustered J matrix and Godambe matrix
  clusters <- unique(SCwt$school)
  nclust <- length(clusters)
  zb <- list()
  for (b in seq_along(clusters)) {
    dat_b <- 
      SCwt |>
      filter(school == clusters[b])
    
    tmp <- select(dat_b, -school, -wt) * dat_b$wt
    tmp <- apply(tmp, 2, sum)
    zb[[b]] <- tmp
  }
  zbar <- apply(do.call(cbind, zb), 1, mean)
  
  B1c <-
    lapply(zb, \(z) tcrossprod(z - zbar)) |>
    Reduce(f = `+`)
  Jc <- nclust / (nclust - 1) * (t(Delta) %*% B1c %*% Delta) / nrow(dat)
  Gc <- Hinv %*% Jc %*% Hinv  # inverse Godambe matrix
  
  # Return standard errors
  out <- sqrt(diag(Gc) / nrow(dat)) 
  names(out) <- names(coef(object))
  out
}

# Test the function and compare with lavaan ------------------------------------
mod_no <- 5
dat <- gen_data_bin_clust(population = make_population(mod_no), n = 1000)
fit <- sem(
  model = txt_mod(mod_no), 
  data = dat, 
  estimator = "PML",
  std.lv = TRUE,
  sampling.weights = "wt"
  # cluster = "school"
)
se_manual <- get_clust_se(fit, dat)

fitc <- sem(
  model = txt_mod(mod_no), 
  data = dat, 
  estimator = "PML",
  std.lv = TRUE,
  sampling.weights = "wt",
  cluster = "school"
)
se_lav <- partable(fitc)$se[partable(fitc)$free != 0]

# Compare
sqrt(sum((se_manual - se_lav) ^ 2))  # should be zero

plot(se_lav, se_manual)
abline(a = 0, b = 1, col = "red")

