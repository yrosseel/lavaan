#' lav_export_estimation
#'
#' lavaan provides a range of optimization methods with the optim.method argument
#' (nlminb, BFGS, L-BFGS-B, GN, and nlminb.constr). `lav_export_estimation`
#' allows exporting objects and functions necessary to pass a lavaan model into
#' any optimizer that takes a combination of (1) starting values, (2) fit-function,
#' (3) gradient-function, and (4) upper and lower bounds. This allows testing new
#' optimization frameworks.
#'
#' @param lavaan_model a fitted lavaan model
#' @returns List with:
#' \itemize{
#'   \item get_coef - When working with equality constraints, lavaan internally
#'   uses some transformations. get_coef is a functions that recreates the coef
#'   function for the parameters.
#'   \item starting_values - starting_values to be used in the optimization
#'   \item objective_function - objective function, expecting the current parameter
#'   values and the lavaan model
#'   \item gradient_function - gradient function, expecting the current parameter
#'   values and the lavaan model
#'   \item lower - lower bounds for parameters
#'   \item upper - upper bound for parameters
#' }
#' @export
#' @examples
#' library(lavaan)
#' model <- "
#'   # latent variable definitions
#'      ind60 =~ x1 + x2 + x3
#'      dem60 =~ y1 + y2 + y3 + y4
#'      dem65 =~ y5 + a*y6 + y7 + y8
#'
#'   # regressions
#'     dem60 ~ ind60
#'     dem65 ~ ind60 + dem60
#' "
#'
#' fit <- sem(model,
#'   data = PoliticalDemocracy,
#'   do.fit = FALSE
#' )
#'
#' est <- lav_export_estimation(lavaan_model = fit)
#'
#' # The starting values are:
#' est$starting_values
#' # Note that these do not have labels (and may also differ from coef(fit)
#' # in case of equality constraints):
#' coef(fit)
#' # To get the same parameters, use:
#' est$get_coef(
#'   parameter_values = est$starting_values,
#'   lavaan_model = fit
#' )
#'
#' # The objective function can be used to compute the fit at the current estimates:
#' est$objective_function(
#'   parameter_values = est$starting_values,
#'   lavaan_model = fit
#' )
#'
#' # The gradient function can be used to compute the gradients at the current estimates:
#' est$gradient_function(
#'   parameter_values = est$starting_values,
#'   lavaan_model = fit
#' )
#'
#' # Together, these elements provide the means to estimate the parameters with a large
#' # range of optimizers. For simplicity, here is an example using optim:
#' est_fit <- optim(
#'   par = est$starting_values,
#'   fn = est$objective_function,
#'   gr = est$gradient_function,
#'   lavaan_model = fit,
#'   method = "BFGS"
#' )
#' est$get_coef(
#'   parameter_values = est_fit$par,
#'   lavaan_model = fit
#' )
#'
#' # This is identical to
#' coef(sem(model,
#'   data = PoliticalDemocracy
#' ))
#'
#' # Example using ridge regularization for parameter a
#' fn_ridge <- function(parameter_values, lavaan_model, est, lambda) {
#'   return(est$objective_function(
#'     parameter_values = parameter_values,
#'     lavaan_model = lavaan_model
#'   ) + lambda * parameter_values[6]^2)
#' }
#' ridge_fit <- optim(
#'   par = est$get_coef(est$starting_values,
#'     lavaan_model = fit
#'   ),
#'   fn = fn_ridge,
#'   lavaan_model = fit,
#'   est = est,
#'   lambda = 10
#' )
#' est$get_coef(
#'   parameter_values = ridge_fit$par,
#'   lavaan_model = fit
#' )
lav_export_estimation <- function(lavaan_model) {
  # define objective function
  objective_function <- function(parameter_values,
                                 lavaan_model) {
    if (lavaan_model@Model@eq.constraints) {
      parameter_values <- as.numeric(lavaan_model@Model@eq.constraints.K %*% parameter_values) +
        lavaan_model@Model@eq.constraints.k0
    }

    # create group list
    GLIST <- lav_model_x2GLIST(lavaan_model@Model, x = parameter_values)
    # get objective function **value**
    fx <- lav_model_objective(
      lavmodel = lavaan_model@Model,
      GLIST = GLIST,
      lavsamplestats = lavaan_model@SampleStats,
      lavdata = lavaan_model@Data,
      lavcache = list(),
      verbose = FALSE
    )
    if (lavaan_model@Options$estimator == "PML") {
      # rescale objective function value
      fx <- fx / lavaan_model@SampleStats@ntotal
    }

    if (!is.finite(fx)) {
      fx.group <- attr(fx, "fx.group")
      fx <- 1e+20
      attr(fx, "fx.group") <- fx.group
    }
    return(fx)
  }

  # define gradient function
  gradient_function <- function(parameter_values,
                                lavaan_model) {
    if (lavaan_model@Model@eq.constraints) {
      parameter_values <- as.numeric(lavaan_model@Model@eq.constraints.K %*% parameter_values) +
        lavaan_model@Model@eq.constraints.k0
    }

    GLIST <- lav_model_x2GLIST(lavaan_model@Model,
      x = parameter_values
    )
    dx <- lav_model_gradient(
      lavmodel = lavaan_model@Model,
      GLIST = GLIST,
      lavsamplestats = lavaan_model@SampleStats,
      lavdata = lavaan_model@Data,
      lavcache = list(),
      type = "free",
      group.weight = !(lavaan_model@SampleStats@missing.flag || lavaan_model@Options$estimator == "PML"),
      verbose = FALSE,
      ceq.simple = lavaan_model@Model@ceq.simple.only
    )

    if (lavaan_model@Model@eq.constraints) {
      dx <- as.numeric(dx %*% lavaan_model@Model@eq.constraints.K)
    }
    if (lavaan_model@Options$estimator == "PML") {
      dx <- dx / lavaan_model@SampleStats@ntotal
    }
    return(dx)
  }

  # extract bounds
  lower <- lavaan_model@ParTable$lower[lavaan_model@ParTable$free > 0L]
  upper <- lavaan_model@ParTable$upper[lavaan_model@ParTable$free > 0L]

  # get starting values
  starting_values <- lav_model_get_parameters(lavaan_model@Model)
  if (lavaan_model@Model@eq.constraints) {
    starting_values <- as.numeric((starting_values - lavaan_model@Model@eq.constraints.k0) %*%
      lavaan_model@Model@eq.constraints.K)
  }

  # lavaan internally uses transformations when there are equality constraints.
  # As a result, the parameters are not necessarily those one would expect when
  # fitting the model. The parameters can be translated with the following function:
  get_coef <- function(parameter_values,
                       lavaan_model) {
    if (lavaan_model@Model@eq.constraints) {
      parameter_values <- as.numeric(lavaan_model@Model@eq.constraints.K %*% parameter_values) +
        lavaan_model@Model@eq.constraints.k0
    }
    names(parameter_values) <- lav_partable_labels(lavaan_model@ParTable,
      type = "free"
    )
    return(parameter_values)
  }

  # Now we just return everything so that the user can use their own optimizer
  return(
    list(
      get_coef = get_coef,
      starting_values = starting_values,
      objective_function = objective_function,
      gradient_function = gradient_function,
      lower = lower,
      upper = upper
    )
  )
}
