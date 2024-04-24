library(lavaan.bingof)
library(tidyverse)

dat <- gen_data_bin_clust(n = 1000)
fit <- sem(
  model = txt_mod(1), 
  data = dat, 
  estimator = "PML",
  std.lv = TRUE,
  sampling.weights = "wt"
  # cluster = "school"
)

# Normal non-clustered, but possibly weighted SE -------------------------------
Delta <- computeDelta(lavmodel = fit@Model)[[1]]
SC <- pml_deriv1(
  Sigma.hat = fit@implied$cov[[1]],
  Mu.hat = fit@implied$mean[[1]],
  TH = fit@implied$th[[1]],
  th.idx = fit@Model@th.idx[[1]],
  num.idx = fit@Model@num.idx[[1]],
  X = fit@Data@X[[1]],
  eXo = NULL,
  wt = NULL,
  PI = NULL,
  lavcache = fit@Cache[[1]],
  missing = fit@Data@missing,
  scores = TRUE,
  negative = FALSE
)
B1 <- crossprod(SC)
J <- t(Delta) %*% B1 %*% Delta / nrow(dat)
Hinv <- lavaan.bingof:::get_sensitivity_inv_mat(fit)
G <- Hinv %*% J %*% Hinv

sqrt(diag(G) / nrow(dat)) 
partable(fit)$se[partable(fit)$free != 0]

# Clustered SE -----------------------------------------------------------------
SCwt <-
  dat |>
  select(school, wt) |>
  bind_cols(as_tibble(SC))

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
Gc <- Hinv %*% Jc %*% Hinv

sqrt(diag(Gc) / nrow(dat)) 
partable(fit)$se[partable(fit)$free != 0]

# Export to ask Irini ----------------------------------------------------------
# SC_ <- SC %*% Delta
# colnames(SC_) <- names(get_true_values(1))
# bind_cols(
#   h = 1:nrow(dat),
#   select(dat, school, wt),
#   as_tibble(SC_)
# ) |>
#   write_csv("SC.csv")

# Function to obtain clustered SE ----------------------------------------------
get_clust_se <- function(object, dat) {
  
  # Obtain scores
  Hinv <- lavaan.bingof:::get_sensitivity_inv_mat(object)  # H^{-1}
  Delta <- computeDelta(lavmodel = object@Model)[[1]]
  SC <- pml_deriv1(
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

# test function
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
get_clust_se(fit, dat)
