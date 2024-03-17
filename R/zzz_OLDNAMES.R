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
