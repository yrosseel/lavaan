# keep 'old' names for some function names that have been used
# (or are still being used) by external packages

lavJacobianD <- lav_func_jacobian_simple
lavJacobianC <- lav_func_jacobian_complex
lavGradientC <- lav_func_gradient_complex

# Myrsini
getHessian <- lav_object_inspect_hessian
getVariability <- lav_object_inspect_firstorder

# rsem
computeExpectedInformation <- lav_model_information_expected

# only for simsem ....
getParameterLabels <- lav_partable_labels

# standardize function names in lav_utils.R / 31 Oct 2025
getCov   <- lav_get_cov
char2num <- lav_char2num
cor2cov  <- lav_cor2cov

# standardize function names in lav_mplus_lavaan / 7 Nov 2025
mplus2lavaan.modelSyntax <- lav_mplus_syntax_model
mplus2lavaan             <- lav_mplus_lavaan

# standardize function names in ctr_estfun.R / 13 Nov 2025
estfun.lavaan <- lav_scores
lavScores <- lav_scores

# standardize function names in lav_export.R
lavExport <- lav_export

# standardize function names in lav_cor.R / 9 December 2025
lavCor <- lav_object_cor

# standardize function names in lav_partable_vnames.R / 9 December 2025
lavNames <- lav_object_vnames
lavaanNames <- lav_object_vnames

# standardize function names in lav_partable.R / 9 December 2025
lavParTable <- lav_model_partable
lavaanify <- lav_model_partable

# standardize function names in lav_simulate_old.R / 9 December 2025
simulateData <- lav_data_simulate_old

# standardize function names in lav_bootstrap.R / 9 December 2025
bootstrapLavaan <- lav_lavaan_bootstrap
