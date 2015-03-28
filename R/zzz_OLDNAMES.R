# keep 'old' names for some function names that have been used
# (or are still being used) by external packages

lavJacobianD <- lav_func_jacobian_simple
lavJacobianC <- lav_func_jacobian_complex
lavGradientC <- lav_func_gradient_complex

getHessian     <- lav_object_inspect_hessian
getVariability <- lav_object_inspect_firstorder

vech   <- lav_matrix_vech
vechr  <- lav_matrix_vechr
vechu  <- lav_matrix_vechu
vechru <- lav_matrix_vechru
vech.reverse <- vechru.reverse <- upper2full <- lav_matrix_vech_reverse
vechr.reverse <- vechu.reverse <- lower2full <- lav_matrix_vechr_reverse
duplicationMatrix <- lav_matrix_duplication
commutationMatrix <- lav_matrix_commutation
sqrtSymmetricMatrix <- lav_matrix_symmetric_sqrt

# rsem
computeExpectedInformation <- lav_model_information_expected
