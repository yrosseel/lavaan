n = 500
data = psych::bfi

model = "A =~ A1+A2+A3+A4+A5;
         C =~ C1+C2+C3+C4+C5
         "
object <- sem(model, data[1:n, ], test = "sb")
sum(diag(ugamma_non_nested(object))) / 34

dim(lavaan::lavInspect(object, "delta"))

model = "A =~ A1+A2+A3+A4+A5;
         C =~ C1+C2+C3+C4+C5
         "
object <- sem(model, data[1:n, ], test = "sb")
dim(lavaan::lavInspect(object, "delta"))

object_chisq <- sem(model, data[1:n, ], test = "sb")

#lavGamma(object_chisq, gamma.unbiased = TRUE)
gamma_biased <- lav_samplestats_Gamma(data[1:n, 1:10])
gamma_yves <- lav_samplestats_Gamma(data[1:n, 1:10], unbiased = TRUE)
gamma_fmg <- lav_fmg_gamma_unbiased(object_chisq, gamma_biased)
sum(abs(gamma_yves-gamma_fmg))

cov_yves = object_chisq@SampleStats@cov[[1]]
cov_fmg = cov(data[1:n, 1:10], use = "pairwise") * (n-1) / n
sum(abs(cov_yves-cov_fmg))


#object_rls <- lavaan::sem(model, psych::bfi[1:n, 1:10], test = "sb", scaled.test = "Browne.residual.nt")


## UNBIASED

GammaNT.cov <- 2 * lav_matrix_duplication_ginv_pre_post(COV %x% COV)

Gamma.u <- (N * (N - 1) / (N - 2) / (N - 3) * Gamma.cov -
              N / (N - 2) / (N - 3) * (GammaNT.cov -
                                         2 / (N - 1) * tcrossprod(cov.vech)))





















lavaan::lavTest(object_chisq, test = "sb", scaled.test = "Browne.residual.nt")


lavTest(object_chisq, scaled.test = "browne.residual.nt", test = "sb")

lav_model_test(lavobject = object_chisq)
lav_model_test(lavobject = object_rls)

test.orig <- object_chisq@Options$test
