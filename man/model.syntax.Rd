\name{model.syntax}
\alias{model.syntax}
\alias{lavaanify}
\alias{lavParTable}
\alias{parseModelString}
\alias{lavParseModelString}
\title{The Lavaan Model Syntax}
\description{
The lavaan model syntax describes a latent variable model. The
function \code{lavParTable} turns it into a table that represents the full
model as specified by the user. We refer to this table as the parameter table.}
\usage{
lavaanify(model = NULL, meanstructure = FALSE, int_ov_free = FALSE,
    int_lv_free = FALSE, marker_int_zero = FALSE,
    orthogonal = FALSE, orthogonal_y = FALSE,
    orthogonal_x = FALSE, orthogonal_efa = FALSE, std_lv = FALSE,
    correlation = FALSE, correlation_delta_free = FALSE,
    composites = TRUE, composites_cov_free = FALSE,
    effect_coding = "", conditional_x = FALSE,
    fixed_x = FALSE, parameterization = "delta", constraints = NULL,
    ceq_simple = FALSE, auto = FALSE, model_type = "sem",
    auto_fix_first = FALSE, marker = NULL,
    auto_fix_single = FALSE, auto_var = FALSE,
    auto_cov_lv_x = FALSE, auto_cov_x = FALSE, auto_cov_y = FALSE,
    auto_th = FALSE, auto_delta = FALSE, auto_efa = FALSE,
    var_table = NULL, ngroups = 1L, nthresholds = NULL,
    group_equal = NULL, group_partial = NULL, group_w_free = FALSE,
    debug = FALSE, warn = TRUE, as_data_frame = TRUE, ...)

lavParTable(model = NULL, meanstructure = FALSE, int_ov_free = FALSE,
    int_lv_free = FALSE, marker_int_zero = FALSE,
    orthogonal = FALSE, orthogonal_y = FALSE,
    orthogonal_x = FALSE, orthogonal_efa = FALSE, std_lv = FALSE,
    correlation = FALSE, correlation_delta_free = FALSE,
    composites = TRUE, composites_cov_free = FALSE,
    effect_coding = "", conditional_x = FALSE,
    fixed_x = FALSE, parameterization = "delta", constraints = NULL,
    ceq_simple = FALSE, auto = FALSE, model_type = "sem",
    auto_fix_first = FALSE, marker = NULL,
    auto_fix_single = FALSE, auto_var = FALSE,
    auto_cov_lv_x = FALSE, auto_cov_x = FALSE, auto_cov_y = FALSE,
    auto_th = FALSE, auto_delta = FALSE, auto_efa = FALSE,
    var_table = NULL, ngroups = 1L, nthresholds = NULL,
    group_equal = NULL, group_partial = NULL, group_w_free = FALSE,
    debug = FALSE, warn = TRUE, as_data_frame = TRUE, ...)

lavParseModelString(model_syntax = '', as_data_frame = FALSE,
                    parser = "open", warn = TRUE, debug = FALSE, ...)
}
\arguments{
\item{model}{A description of the user-specified model. Typically, the model
    is described using the lavaan model syntax; see details for more
    information. Alternatively, a parameter table (e.g., the output of
    \code{lavParseModelString}) is also accepted.}
\item{model_syntax}{The model syntax specifying the model. Must be a literal
string.}
\item{meanstructure}{If \code{TRUE}, intercepts/means will be added to
the model for both observed and latent variables.}
\item{int_ov_free}{If \code{FALSE}, the intercepts of the observed variables
    are fixed to zero.}
\item{int_lv_free}{If \code{FALSE}, the intercepts of the latent variables
    are fixed to zero.}
\item{marker_int_zero}{Logical. Only relevant if the metric of each latent
    variable is set by fixing the first factor loading to unity.
    If \code{TRUE}, it implies \code{meanstructure = TRUE} and
    \code{std_lv = FALSE}, and it fixes the intercepts of the marker
    indicators to zero, while freeing the means/intercepts of the latent
    variables. This only works correctly for single-group, single-level
    models.}
\item{orthogonal}{If \code{TRUE}, all covariances among
   latent variables are set to zero.}
\item{orthogonal_y}{If \code{TRUE}, all covariances among
   endogenous latent variables only are set to zero.}
\item{orthogonal_x}{If \code{TRUE}, all covariances among
   exogenous latent variables only are set to zero.}
\item{orthogonal_efa}{If \code{TRUE}, all covariances among
   latent variables involved in rotation only are set to zero.}
\item{std_lv}{If \code{TRUE}, the metric of each latent variable is
    determined by fixing their variances to 1.0. If \code{FALSE}, the metric
    of each latent variable is determined by fixing the factor loading of the
    first indicator to 1.0. If there are multiple
      groups, \code{std_lv = TRUE} and \code{"loadings"} is included in
      the \code{group.label} argument, then only the latent variances
      of the first group will be fixed to 1.0, while the latent
      variances of the other groups are freely estimated.}
\item{correlation}{If \code{TRUE}, a correlation structure is fitted. For
    continuous data, this implies that the (residual) variances are no
    longer parameters of the model.}
\item{correlation_delta_free}{Logical. Only used if
    \code{correlation = TRUE} (or a character vector). If \code{TRUE},
    the scaling (\code{~*~}) parameters are free scale (standard
    deviation) parameters, and the model is the D-augmented covariance
    structure \code{Sigma = Delta P(theta) Delta} with unit-diagonal
    \code{P} (used internally when an ML-family estimator is combined
    with a correlation structure).}
\item{composites}{Logical. If \code{TRUE}, use the new (0.6-20) approach
     to handle composites.}
\item{composites_cov_free}{Logical. Only relevant if the model contains
     composites. If \code{TRUE}, the composite-indicator (co)variances (the
     elements of the \code{T} matrix) are estimated as free parameters. If
     \code{FALSE} (the default), they are fixed to their sample values.}
\item{effect_coding}{Can be logical or character string. If
     logical and \code{TRUE}, this implies
     \code{effect_coding = c("loadings", "intercepts")}. If logical and
     \code{FALSE}, it is set equal to the empty string.
     If \code{"loadings"} is included, equality
     constraints are used so that the average of the factor loadings (per
     latent variable) equals 1. Note that this should not be used
     together with \code{std_lv = TRUE}. If \code{"intercepts"} is
     included, equality constraints are used so that the sum of the
     intercepts (belonging to the indicators of a single latent variable)
     equals zero.
     As a result, the latent mean will be freely estimated and will usually
     equal the average of the means of the indicators involved.}
\item{conditional_x}{If \code{TRUE}, we set up the model conditional on
    the exogenous `x' covariates; the model-implied sample statistics
    only include the non-x variables. If \code{FALSE}, the exogenous `x'
    variables are modeled jointly with the other variables, and the
    model-implied statistics reflect both sets of variables.}
\item{fixed_x}{If \code{TRUE}, the exogenous `x' covariates are considered
    fixed variables and the means, variances and covariances of these variables
    are fixed to their sample values. If \code{FALSE}, they are considered
    random, and the means, variances and covariances are free parameters.}
\item{parameterization}{Currently only used if data is categorical. If
    \code{"delta"}, the delta parameterization is used. If \code{"theta"},
    the theta parameterization is used.}
\item{constraints}{Additional (in)equality constraints. See details for
    more information.}
\item{ceq_simple}{If \code{TRUE}, and no other general constraints are
    used in the model, simple equality constraints are represented in the
    parameter table as duplicated free parameters (instead of extra rows
    with \code{op = "=="}).}
\item{auto}{If \code{TRUE}, the default values are used for the auto.*
    arguments, depending on the value of \code{model_type}.}
\item{model_type}{Either \code{"sem"} or \code{"growth"}; only used if
    \code{auto=TRUE}.}
\item{auto_fix_first}{If \code{TRUE}, the factor loading of the first indicator
    is set to 1.0 for every latent variable.}
\item{marker}{Optional named character vector mapping a latent variable
    (name) to the observed indicator (name) whose loading should be fixed to
    1.0 (instead of the first indicator) when \code{auto_fix_first = TRUE}.
    Latent variables that are not named in this vector keep the first
    indicator as their marker. This is used internally by the
    \code{bad.marker.crit} mechanism.}
\item{auto_fix_single}{If \code{TRUE}, the residual variance (if included)
    of an observed indicator is set to zero if it is the only indicator of a
    latent variable.}
\item{auto_var}{If \code{TRUE}, the (residual) variances of both observed
    and latent variables are set free.}
\item{auto_cov_lv_x}{If \code{TRUE}, the covariances of exogenous latent
    variables are included in the model and set free.}
\item{auto_cov_x}{If \code{TRUE}, the covariances between exogenous latent
    variables and observed exogenous covariates are also included in the
    model and set free (this implies \code{auto_cov_lv_x = TRUE}). Ignored
    (with a warning) if \code{conditional_x = TRUE}.}
\item{auto_cov_y}{If \code{TRUE}, the covariances of dependent variables
    (both observed and latent) are included in the model and set free.}
\item{auto_th}{If \code{TRUE}, thresholds for limited dependent variables
    are included in the model and set free.}
\item{auto_delta}{If \code{TRUE}, response scaling parameters for limited
    dependent variables are included in the model and set free.}
\item{auto_efa}{If \code{TRUE}, the necessary constraints are
    imposed to make the (unrotated) exploratory factor analysis blocks
    identifiable: for each block, factor variances are set to 1, factor
    covariances are constrained to be zero, and factor loadings are
    constrained to follow an echelon pattern.}
\item{var_table}{The variable table containing information about the
    observed variables in the model.}
\item{ngroups}{The number of (independent) groups.}
\item{nthresholds}{Either a single integer or a named vector of integers.
    If \code{nthresholds} is a single integer, all endogenous
    variables are assumed to be ordered with \code{nthresholds} indicating
    the number of thresholds needed in the model. If \code{nthresholds} is a
    named vector, it indicates the number of thresholds for these ordered
    variables only. This argument should not be used in combination with
    var_table.}
\item{group_equal}{A vector of character strings. Only used in
    a multiple group analysis. Can be one or more of the following:
    \code{"loadings"}, \code{"intercepts"},
    \code{"means"}, \code{"regressions"}, \code{"residuals"} or
    \code{"covariances"}, specifying the pattern of equality
    constraints across multiple groups. When (in the model syntax) a vector
    of labels is used as a modifier for a certain parameter, this will
    override the group_equal setting if it applies to this parameter.
    See also the Multiple groups section below for using modifiers in multiple
    groups.}
\item{group_partial}{A vector of character strings containing the labels
    of the parameters which should be free in all groups (thereby
    overriding the group_equal argument for some specific parameters).}
\item{group_w_free}{Logical. If \code{TRUE}, the group frequencies are
    considered to be free parameters in the model. In this case, a
    Poisson model is fitted to estimate the group frequencies. If
    \code{FALSE} (the default), the group frequencies are fixed to their
    observed values.}
\item{as_data_frame}{If \code{TRUE}, return the list of model parameters
    as a \code{data.frame}.}
\item{parser}{Character. If \code{"old"}, use the original/classic parser.
    If \code{"new"}, use the new/ldw parser.
    If \code{"open"}, use the newest, extensible parser.
    The default (as of version 0.7-1) is \code{"open"}.}
\item{warn}{If \code{TRUE}, some (possibly harmless) warnings are printed
    out.}
\item{debug}{If \code{TRUE}, debugging information is printed out.}
\item{...}{To accept old argument names with dots. No other arguments are accepted.}
}
\details{
The model syntax consists of one or more formula-like expressions, each one
describing a specific part of the model. The model syntax can be read from
a file (using \code{\link{readLines}}), or can be specified as a literal
string enclosed by single quotes as in the example below.
\preformatted{myModel <- '
  # 1. latent variable definitions
    f1 =~ y1 + y2 + y3
    f2 =~ y4 + y5 + y6
    f3 =~ y7 + y8 +
          y9 + y10
    f4 =~ y11 + y12 + y13

  ! this is also a comment

  # 2. regressions
    f1 ~ f3 + f4
    f2 ~ f4
    y1 + y2 ~ x1 + x2 + x3

  # 3. (co)variances
    y1 ~~ y1
    y2 ~~ y4 + y5
    f1 ~~ f2

  # 4. intercepts
    f1 ~ 1; y5 ~ 1

  # 5. thresholds
    y11 | t1 + t2 + t3
    y12 | t1
    y13 | t1 + t2

  # 6. scaling factors
    y11 ~*~ y11
    y12 ~*~ y12
    y13 ~*~ y13

  # 7. composites
    f5 <~ z1 + z2 + z3 + z4
'
}
Blank lines and comments can be used in between the formulas, and formulas can
be split over multiple lines. Both the sharp (#) and the exclamation (!)
characters can be used to start a comment. Multiple formulas can be placed
on a single line if they are separated by a semicolon (;).

There can be seven types of formula-like expressions in the model syntax:
\enumerate{

    \item Latent variable definitions: The \code{"=~"} operator can be
      used to define (continuous) latent variables. The name of the latent
      variable is on the left of the \code{"=~"} operator, while the terms
      on the right, separated by \code{"+"} operators, are the indicators
      of the latent variable.

      The operator \code{"=~"} can be read as ``is manifested by''.

    \item Regressions: The \code{"~"} operator specifies a regression.
      The dependent variable is on the left of a \code{"~"} operator and the
      independent variables, separated by \code{"+"} operators, are on the right.
      These regression formulas are similar to the way ordinary linear regression
      formulas are used in R, but they may include latent variables. Interaction
      terms are currently not supported.

    \item Variance-covariances: The \code{"~~"} (`double tilde') operator specifies
     (residual) variances of an observed or latent variable, or a set of
     covariances between one variable and several other variables (either
     observed or latent). Several variables, separated by \code{"+"}
     operators, can appear on the right. In this way, several pairwise
     (co)variances involving the same left-hand variable can be expressed in a
     single expression. The distinction between variances and residual variances
     is made automatically.

    \item Intercepts: A special case of a regression formula can be used to
      specify an intercept (or a mean) of either an observed or a latent variable.
     The variable name is on the left of a \code{"~"} operator, and the only
     term on the right is the number \code{"1"}, representing the intercept.
     Including an intercept
     formula in the model automatically implies \code{meanstructure = TRUE}. The
     distinction between intercepts and means is made automatically.

    \item Thresholds: The \code{"|"} operator can be used to define the
      thresholds of categorical endogenous variables (on the left-hand side
      of the operator). By convention, the
      thresholds (on the right-hand side, separated by the \code{"+"} operator)
      are named \code{"t1"}, \code{"t2"}, and so on.

    \item Scaling factors: The \code{"~*~"} operator defines a scale factor.
      The variable name on the left hand side must be the same as the variable
      name on the right hand side. Scale factors are used in the Delta
      parameterization, in a multiple group analysis when factor indicators
      are categorical.

    \item Composites: The \code{"<~"} operator can be used to define
      a composite (on the right hand side of the operator).
      A composite is a weighted linear combination of its composite indicators.
      The name of the composite variable is on the left of the \code{"<~"}
      operator, while the terms on the right, separated by \code{"+"}
      operators, are the indicators of the composite variable.
}
There are 4 additional operators, also with left- and right-hand sides, that can
be included in model syntax.  Three of them are used to specify (in)equality
constraints on estimated parameters (\code{==}, \code{>}, and \code{<}), and
those are demonstrated in a later section about
\bold{(In)equality constraints}.
The final additional operator (\code{:=}) can be used to define ``new'' parameters
that are functions of one or more other estimated parameters. The \code{:=}
operator is demonstrated in a section about \bold{User-defined parameters}.

Usually, only a single variable name appears on the left side of an
operator. However, if multiple variable names are specified,
separated by the \code{"+"} operator, the formula is repeated for each
element on the left side (as for example in the third regression
formula in the example above). The only exception is scaling factors, where
only a single element is allowed on the left-hand side.

In the right-hand side of these formula-like expressions, each element can be
modified (using the \code{"*"} operator) by either a numeric constant,
an expression resulting in a numeric constant, an expression resulting
in a character vector, or one
of three special functions: \code{start()}, \code{label()} and \code{equal()}.
This provides the user with a mechanism to fix parameters, to provide
alternative starting values, to label the parameters, and to define equality
constraints among model parameters. All \code{"*"} expressions are
referred to as \emph{modifiers}. They are explained in more detail in the
following sections.
}

\section{Fixing parameters}{
It is often desirable to fix a model parameter that is otherwise (by default)
free. Any parameter in a model can be fixed by using a modifier resulting
in a numerical constant. Here are some examples:
\itemize{

    \item Fixing the regression coefficient of the predictor
       \code{x2}: \preformatted{y ~ x1 + 2.4*x2 + x3}

    \item Specifying an orthogonal (zero) covariance between two latent
      variables: \preformatted{f1 ~~ 0*f2}

    \item Specifying an intercept and a linear slope in a growth
    model: \preformatted{i =~ 1*y11 + 1*y12 + 1*y13 + 1*y14
s =~ 0*y11 + 1*y12 + 2*y13 + 3*y14}
}
Instead of a numeric constant, one can use a mathematical function that returns
a numeric constant, for example \code{sqrt(10)}. Multiplying with \code{NA}
will force the corresponding parameter to be free.

Additionally, the \code{==} operator can be used to set a \emph{labeled} parameter
equal to a specific numeric value. This will be demonstrated in the section below
about \bold{(In)equality constraints}.
}

\section{Starting values}{
User-provided starting values can be given by using the special function
\code{start()}, containing a numeric constant. For example:
\preformatted{y ~ x1 + start(1.0)*x2 + x3}
Note that if a starting value is provided, the parameter is not
automatically considered to be free.
}

\section{Parameter labels and equality constraints}{
Each free parameter in a model is automatically given a name (or label).
The name given to a model
parameter consists of three parts, coerced to a single character vector.
The first part is the name of the variable in the left-hand side of the
formula where the parameter was
implied. The middle part is based on the special `operator' used in the
formula. This can be either one of \code{"=~"}, \code{"~"} or \code{"~~"}. The
third part is the name of the variable in the right-hand side of the formula
where the parameter was implied, or \code{"1"} if it is an intercept. The three
parts are pasted together in a single string. For example, the name of the
fixed regression coefficient in the regression formula
\code{y ~ x1 + 2.4*x2 + x3} is the string \code{"y~x2"}.
The name of the parameter
corresponding to the covariance between two latent variables in the
formula \code{f1 ~~ f2} is the string \code{"f1~~f2"}.

Although this automatic labeling of parameters is convenient, the user may
specify their own labels for specific parameters simply by pre-multiplying
the corresponding term (on the right hand side of the operator only) by
a character string (starting with a letter).
For example, in the formula \code{f1 =~ x1 + x2 + mylabel*x3}, the parameter
corresponding with the factor loading of
\code{x3} will be named \code{"mylabel"}.
An alternative way to specify the label is as follows:
\code{f1 =~ x1 + x2 + label("mylabel")*x3},
where the label is the argument of special function \code{label()};
this can be useful if the label contains a space, or an operator (like "~").

There are two ways to constrain a parameter
to be equal to another target parameter. If you
have specified your own labels, you can use the fact that
\emph{equal labels imply equal parameter values}.
If you rely on automatic parameter labels, you
can use the special function \code{equal()}. The argument of
\code{equal()} is the (automatic or user-specified) name of the target
parameter. For example, in the confirmatory factor analysis example below, the
intercepts of the three indicators of each latent variable are constrained to
be equal to each other. For the first three, we have used the default
names. For the last three, we have provided a custom label for the
\code{y2a} intercept.
\preformatted{model <- '
  # two latent variables with fixed loadings
    f1 =~ 1*y1a + 1*y1b + 1*y1c
    f2 =~ 1*y2a + 1*y2b + 1*y2c

  # intercepts constrained to be equal
  # using the default names
    y1a ~ 1
    y1b ~ equal("y1a~1") * 1
    y1c ~ equal("y1a~1") * 1

  # intercepts constrained to be equal
  # using a custom label
    y2a ~ int2*1
    y2b ~ int2*1
    y2c ~ int2*1
'
}
}

\section{Multiple groups}{
In a multiple group analysis, modifiers that contain a single element
should be replaced by a vector, having the same length as the number
of groups. If you provide a single element, it will be recycled
for all the groups. This may be dangerous, in particular when the modifier
is a label. In that case, the (same) label is copied across all groups,
and this would imply an equality constraint across groups.
Therefore, when using modifiers in a multiple group setting, it is always
safer (and cleaner) to specify the same number of
elements as the number of groups. Consider this example with two groups:
\preformatted{
HS.model <- ' visual  =~ x1 + 0.5*x2 + c(0.6, 0.8)*x3
              textual =~ x4 + start(c(1.2, 0.6))*x5 + x6
              speed   =~ x7 + x8 + c(x9.group1, x9.group2)*x9 '
}
In this example, the factor loading of the `x2' indicator is fixed to the
value 0.5 for both groups. However, the factor loadings of the `x3' indicator
are fixed to 0.6 and 0.8 for group 1 and group 2 respectively. The same
logic is used for all modifiers. Note that character vectors can contain
unquoted strings.
}

\section{Multiple modifiers}{
In the model syntax, you can specify a variable more than once on the right hand
side of an operator; therefore, several `modifiers' can be applied
simultaneously; for example, if you want to fix the value of a parameter and
also label that parameter, you can use something like:
\preformatted{ f1 =~ x1 + x2 + 4*x3 + x3.loading*x3}
}

\section{Random slopes (two-level models)}{
In a two-level model (with a \code{cluster=} argument), the \code{rv()}
modifier can be used in the level-1 (within) part of the model to turn a
regression coefficient into a \emph{random slope}: a latent variable at
the second (between) level. The (quoted) argument of \code{rv()} is the
name of this latent variable. At the between level, the random slope can
then be given a (residual) variance, covariances, regressions on
between-level covariates, and an intercept (which is the `fixed effect'
of the slope). For example (cfr. ex9.8 in the Mplus User's Guide):
\preformatted{model <- '
  level: 1
    fw =~ y1 + y2 + y3 + y4
    fw ~ rv("s1")*x1 + rv("s2")*x2
  level: 2
    fb =~ y1 + y2 + y3 + y4
    s1 + s2 + fb ~ w
'
}
Here, the within-level regression coefficients of \code{fw} on the
observed covariates \code{x1} and \code{x2} vary across clusters; the
random slopes \code{s1} and \code{s2} are regressed on the between-level
covariate \code{w}.

The covariate that carries a random slope may also be a \emph{latent}
(within-level) variable, as in
\preformatted{model <- '
  level: 1
    fxw =~ x1 + x2 + x3
    fyw =~ y1 + y2 + y3
    fyw ~ rv("s1")*fxw
  level: 2
    fxb =~ x1 + x2 + x3
    fyb =~ y1 + y2 + y3
    fyb ~ fxb
    fyb ~~ s1
'
}
Such a slope multiplies two random variables, and the marginal
loglikelihood is no longer available in closed form: it is computed by
(Gauss-Hermite) quadrature over the latent-covariate slopes only,
while all other random effects are still integrated out analytically
(Rockwood, 2020). The number of quadrature points (per dimension) can
be set via the \code{integration.ngh} option (default: 21).

The same route is taken when the covariate is a `split' both-level
observed variable (an observed variable that appears in both the
level-1 and the level-2 part of the model): the random slope then
multiplies its \emph{latent} within-cluster component (latent
centering), and the covariate is jointly modeled:
\preformatted{model <- '
  level: 1
    y ~ rv("s1")*x
  level: 2
    y ~ x
    y ~~ s1
'
}
This is the latent-covariate contextual model of Rockwood (2020,
section 4). Note the contrast with a within-only covariate (a
variable that does \emph{not} appear in the level-2 part): there, the
slope multiplies the observed values, which are conditioned upon.

The current implementation is limited to: two-level models, a single
group, continuous data, \code{estimator = "ML"} and \code{fixed.x =
TRUE}; the observed covariates that carry a random slope must be
within-only (level-1) variables. Missing data can be handled by
full-information maximum likelihood (\code{missing = "ml"}); the
exogenous covariates themselves must be complete (cases with missing
values on exogenous variables are removed). The (observed-data)
loglikelihood is \emph{conditional} on \emph{all} exogenous
covariates: the covariates carrying random slopes, any other
within-only exogenous covariates (as in Mplus
\code{TYPE = TWOLEVEL RANDOM}), and the between-level exogenous
covariates. Between-only \emph{endogenous} observed variables
(level-2 outcomes, e.g., indicators of a between-level factor) are
supported; they must be pure indicators (not involved in
regressions), with a free (or fixed to a nonzero value) residual
variance, and without missing values. No
chi-square test statistic (or fit indices) are available for models
with random slopes; \code{fitMeasures()} provides the
loglikelihood-based measures only (\code{npar}, \code{ntotal},
\code{logl}, \code{aic}, \code{bic}, \code{bic2}), plus (with
\code{estimator = "MLR"}) the H0 scaling correction factor
(\code{scaling.factor.h0}), which can be used for scaled
(loglikelihood) difference tests. Standard errors can be
\code{"standard"} (observed
information) or \code{"robust.huber.white"} (sandwich; this is also
what \code{estimator = "MLR"} implies). Both \code{optim.method =
"nlminb"} and \code{optim.method = "em"} are supported; with complete
data, \code{"nlminb"} is the default, while with \code{missing =
"ml"}, the EM algorithm is the default optimizer.
\code{lavPredict()} provides empirical Bayes (posterior mean)
predictions: with \code{level = 2}, the
between-level latent variables including the random slopes (and, with
\code{se = "standard"}, their posterior standard deviations); with
\code{level = 1} (the default), the within-level factor scores. For
latent-covariate models, these are EAP (expected a posteriori)
estimates: the posterior is a mixture over the quadrature nodes.

For random slopes with latent covariates, some additional restrictions
apply (for now): a random-slope label may not be attached to both an
observed and a latent covariate, and the optimizer is always
\code{"nlminb"}. The latent-covariate slopes may be regressed on
between-level covariates (e.g., \code{s1 ~ w} in the level-2 part),
as in the PISA example of Rockwood (2020).
}

\section{(In)equality constraints}{
The \code{==} operator can be used either to fix a parameter to a specific value,
or to set an estimated parameter equal to another parameter. Adapting the
example in the \bold{Parameter labels and equality constraints} section, we
could have used different labels for the second factor's intercepts:
\preformatted{
    y2a ~ int1*1
    y2b ~ int2*1
    y2c ~ int3*1
}
Then, we could fix the first intercept to zero by including in the syntax an
operation that indicates the parameter's label equals that value:
\preformatted{
    int1 == 0
}
Whereas we could still estimate the other two intercepts under an equality
constraint by setting their different labels equal to each other:
\preformatted{
    int2 == int3
}
Optimization can be less efficient when constraining parameters this way (see
the documentation linked under \bold{See also} for more information). But the
flexibility might be advantageous.  For example, the constraints could be
specified in a separate character-string object, which can be passed to the
\code{lavaan(..., constraints=)} argument, enabling users to compare results
with(out) the constraints.

Inequality constraints work in much the same way, using the \code{<} or \code{>}
operator to indicate which estimated parameter is hypothesized to be greater or
less than either a specific value or another estimated parameter. For example, a
variance can be constrained to be nonnegative:
\preformatted{
    y1a ~~ var1a*y1a
    ## hypothesized constraint:
    var1a > 0
}
Or the factor loading of a particular indicator might be expected to exceed
other indicators' loadings:
\preformatted{
    f1 =~ L1*y1a + L2*y1b + L3*y1c
    ## hypothesized constraints:
    L1 > L2
    L3 < L1
}
}

\section{User-defined parameters}{
Functions of parameters can be useful to test particular hypotheses. Following
from the \code{Multiple groups} example, we might be interested in which group's
factor loading is larger (i.e., an estimate of differential item functioning
(DIF) when the latent scales are linked by anchor items with equal loadings).
\preformatted{
  speed =~ c(L7, L7)*x7 + c(L8, L8)*x8 + c(L9.group1, L9.group2)*x9 '
  ## user-defined parameter:
  DIF_L9 := L9.group1 - L9.group2
}
Note that this hypothesis is easily tested without a user-defined parameter by
using the \code{lavTestWald()} function. However, a user-defined parameter
additionally provides an estimate of the parameter being tested.

User-defined parameters are particularly useful for specifying indirect effects
in models of mediation. For example:
\preformatted{
model <- ' # direct effect
             Y ~ c*X
           # mediator
             M ~ a*X
             Y ~ b*M

# user defined parameters:

    # indirect effect (a*b)
             ab := a*b
    # total effect (defined using another user-defined parameter)
             total := ab + c
'
}
}

\references{
Rockwood, N. J. (2020). Maximum likelihood estimation of multilevel
structural equation models with random slopes for latent covariates.
\emph{Psychometrika, 85}(2), 275--300.
\doi{10.1007/s11336-020-09702-9}

Rosseel, Y. (2012). \code{lavaan}: An R package for structural equation
modeling. \emph{Journal of Statistical Software, 48}(2), 1--36.
\doi{10.18637/jss.v048.i02}
}

