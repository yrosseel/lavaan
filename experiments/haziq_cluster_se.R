library(lavaan.bingof)

dat <- gen_data_bin_clust(n = 1000)
fit <- sem(
  model = txt_mod(1), 
  data = dat, 
  estimator = "PML",
  std.lv = TRUE
  # sampling.weights = "wt"
  # cluster = "school"
)

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
for (b in 1) {
  dat_b <- 
    SCwt |>
    filter(school == clusters[b])
  
  zab <- select(dat_b, -school, -wt)
  zab <- zab * dat_b$wt
}

# stuck!!!

# Export to ask Irini ----------------------------------------------------------
SC_ <- SC %*% Delta
colnames(SC_) <- names(get_true_values(1))
bind_cols(
  h = 1:nrow(dat),
  select(dat, school, wt),
  as_tibble(SC_)
) |>
  write_csv("SC.csv")
