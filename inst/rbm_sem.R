# Reduced Bias M-estimation for SEM.

mod <- "
  # latent variables 
  ind60 =~ x1 + x2 + x3 
  dem60 =~ y1 + y2 + y3 + y4 
  dem65 =~ y5 + y6 + y7 + y8 
   
  # regressions
  dem60 ~ ind60 
  dem65 ~ ind60 + dem60 
  
  # residual covariances 
  y1 ~~ y5
  y2 ~~ y4 + y6 
  y3 ~~ y7 
  y4 ~~ y8
  y6 ~~ y8
"

# Maximum likelihood fit
fit_ML   <- sem(model = mod, data = PoliticalDemocracy)

# Implicit RBM fit (default)
fit_iRBM <- sem(model = mod, data = PoliticalDemocracy,
                estimator = list(estimator = "rbm", rbm.method = "implicit")) 

# Explicit RBM fit
fit_eRBM <- sem(model = mod, data = PoliticalDemocracy,
                estimator = list(estimator = "rbm", rbm.method = "explicit"))

# Compare
tab <- data.frame(
    ML = coef(fit_ML),
    iRBM = coef(fit_iRBM),
    eRBM = coef(fit_eRBM)
)
rownames(tab) <- names(coef(fit_ML))
print(round(tab, 3))
