\name{model.syntax}
\alias{model.syntax}
\alias{lavaanify}
\alias{lavParTable}
\alias{lavpartable}
\alias{lavPartable}
\alias{parseModelString}
\alias{lavParseModelString}
\title{The Lavaan Model Syntax}
\description{
The lavaan model syntax describes a latent variable model. The 
function \code{lavaanify} turns it into a table that represents the full
model as specified by the user. We refer to this table as the parameter table.} 
\usage{
lavaanify(model = NULL, meanstructure = FALSE, int.ov.free = FALSE, 
    int.lv.free = FALSE, orthogonal = FALSE, std.lv = FALSE, 
    conditional.x = FALSE, fixed.x = TRUE, parameterization = "delta", 
    constraints = NULL, auto = FALSE, model.type = "sem", 
    auto.fix.first = FALSE, auto.fix.single = FALSE, auto.var = FALSE, 
    auto.cov.lv.x = FALSE, auto.cov.y = FALSE, auto.th = FALSE, 
    auto.delta = FALSE, varTable = NULL, ngroups = 1L, group.equal = NULL, 
    group.partial = NULL, group.w.free = FALSE,
    debug = FALSE, warn = TRUE, as.data.frame. = TRUE)

lavParTable(model = NULL, meanstructure = FALSE, int.ov.free = FALSE,
    int.lv.free = FALSE, orthogonal = FALSE, std.lv = FALSE,
    conditional.x = FALSE, fixed.x = TRUE, parameterization = "delta", 
    constraints = NULL, auto = FALSE, model.type = "sem",
    auto.fix.first = FALSE, auto.fix.single = FALSE, auto.var = FALSE,
    auto.cov.lv.x = FALSE, auto.cov.y = FALSE, auto.th = FALSE,
    auto.delta = FALSE, varTable = NULL, ngroups = 1L, group.equal = NULL,
    group.partial = NULL, group.w.free = FALSE,
    debug = FALSE, warn = TRUE, as.data.frame. = TRUE)

lavParseModelString(model.syntax = '', as.data.frame. = FALSE, 
                    warn = TRUE, debug = FALSE)

}
\arguments{
\item{model}{A description of the user-specified model. Typically, the model
    is described using the lavaan model syntax; see details for more 
    information. Alternatively, a parameter table (e.g., the output of 
    \code{lavParseModelString} is also accepted.}
\item{model.syntax}{The model syntax specifying the model. Must be a literal
string.}
\item{meanstructure}{If \code{TRUE}, intercepts/means will be added to
the model both for both observed and latent variables.}
\item{int.ov.free}{If \code{FALSE}, the intercepts of the observed variables
    are fixed to zero.}
\item{int.lv.free}{If \code{FALSE}, the intercepts of the latent variables
    are fixed to zero.}
\item{orthogonal}{If \code{TRUE}, the exogenous latent variables
    are assumed to be uncorrelated.}
\item{std.lv}{If \code{TRUE}, the metric of each latent variable is
    determined by fixing their variances to 1.0. If \code{FALSE}, the metric
    of each latent variable is determined by fixing the factor loading of the
    first indicator to 1.0.}
\item{conditional.x}{If \code{TRUE}, we set up the model conditional on
    the exogenous `x' covariates; the model-implied sample statistics
    only include the non-x variables. If \code{FALSE}, the exogenous `x'
    variables are modeled jointly with the other variables, and the
    model-implied statistics refect both sets of variables.}
\item{fixed.x}{If \code{TRUE}, the exogenous `x' covariates are considered
    fixed variables and the means, variances and covariances of these variables
    are fixed to their sample values. If \code{FALSE}, they are considered
    random, and the means, variances and covariances are free parameters.}
\item{parameterization}{Currently only used if data is categorical. If
    \code{"delta"}, the delta parameterization is used. If \code{"theta"},
    the theta parameterization is used.}
\item{constraints}{Additional (in)equality constraints. See details for
    more information.}
\item{auto}{If \code{TRUE}, the default values are used for the auto.*
    arguments, depending on the value of \code{model.type}.}
\item{model.type}{Either \code{"sem"} or \code{"growth"}; only used if
    \code{auto=TRUE}.}
\item{auto.fix.first}{If \code{TRUE}, the factor loading of the first indicator
    is set to 1.0 for every latent variable.}
\item{auto.fix.single}{If \code{TRUE}, the residual variance (if included)
    of an observed indicator is set to zero if it is the only indicator of a
    latent variable.}
\item{auto.var}{If \code{TRUE}, the residual variances and the variances
    of exogenous latent variables are included in the model and set free.}
\item{auto.cov.lv.x}{If \code{TRUE}, the covariances of exogenous latent
    variables are included in the model and set free.}
\item{auto.cov.y}{If \code{TRUE}, the covariances of dependent variables
    (both observed and latent) are included in the model and set free.}
\item{auto.th}{If \code{TRUE}, thresholds for limited dependent variables
    are included in the model and set free.}
\item{auto.delta}{If \code{TRUE}, response scaling parameters for limited
    dependent variables are included in the model and set free.}
\item{varTable}{The variable table containing information about the
    observed variables in the model.}
\item{ngroups}{The number of (independent) groups.}
\item{group.equal}{A vector of character strings. Only used in
    a multiple group analysis. Can be one or more of the following:
    \code{"loadings"}, \code{"intercepts"},
    \code{"means"}, \code{"regressions"}, \code{"residuals"} or
    \code{"covariances"}, specifying the pattern of equality
    constraints across multiple groups.}
\item{group.partial}{A vector of character strings containing the labels
    of the parameters which should be free in all groups (thereby
    overriding the group.equal argument for some specific parameters).}
\item{group.w.free}{Logical. If \code{TRUE}, the group frequencies are
    considered to be free parameters in the model. In this case, a
    Poisson model is fitted to estimate the group frequencies. If
    \code{FALSE} (the default), the group frequencies are fixed to their
    observed values.}
\item{warn}{If \code{TRUE}, some (possibly harmless) warnings are printed
    out.}
\item{as.data.frame.}{If \code{TRUE}, return the list of model parameters
as a \code{data.frame}.}
\item{debug}{If \code{TRUE}, debugging information is printed out.}
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

  # 7. formative factors
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
     covariances between one variable, and several other variables (either
     observed or latent). Several variables, separated by \code{"+"} 
     operators can appear on the right. This way, several pairwise 
     (co)variances involving the same left-hand variable can be expressed in a
     single expression. The distinction between variances and residual variances
     is made automatically.

    \item Intercepts: A special case of a regression formula can be used to
      specify an intercept (or a mean) of either an observed or a latent variable.
     The variable name is on the left of a \code{"~"} operator. On the right is
     only the number \code{"1"} representing the intercept. Including an intercept
     formula in the model automatically implies \code{meanstructure = TRUE}. The
     distinction between intercepts and means is made automatically.

    \item Thresholds: The \code{"|"} operator can be used to define the 
      thresholds of categorical endogenous variables (on the left hand side
      of the operator). By convention, the
      thresholds (on the right hand sided, separated by the \code{"+"} operator,
      are named \code{"t1"}, \code{"t2"}, etcetera.

    \item Scaling factors: The \code{"~*~"} operator defines a scale factor.
      The variable name on the left hand side must be the same as the variable
      name on the right hand side. Scale factors are used in the Delta
      parameterization, in a multiple group analysis when factor indicators
      are categorical.

    \item Formative factors: The \code{"<~"} operator can be used to define
      a formative factor (on the right hand side of the operator), in a
      similar why as a reflexive factor is defined (using the \code{"=~"}
      operator). This is just syntax sugar to define a phantom latent
      variable (equivalent to using \code{"f =~ 0"}). And in addition, the
      (residual) variance of the formative factor is fixed to zero.
}
Usually, only a single variable name appears on the left side of an
operator. However, if multiple variable names are specified, 
separated by the \code{"+"} operator, the formula is repeated for each
element on the left side (as for example in the third regression 
formula in the example above). The only exception are scaling factors, where
only a single element is allowed on the left hand side. 

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
in a numerical constaint. Here are some examples:
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
specify its own labels for specific parameters simply by pre-multiplying
the corresponding term (on the right hand side of the operator only) by
a character string (starting with a letter).
For example, in the formula \code{f1 =~ x1 + x2 + mylabel*x3}, the parameter 
corresponding with the factor loading of
\code{x3} will be named \code{"mylabel"}.
An alternative way to specify the label is as follows: 
\code{f1 =~ x1 + x2 + label("mylabel")*x3}, 
where the label is the argument of special function \code{label()};
this can be useful if the label contains a space, or an operator (like "~").

To constrain a parameter
to be equal to another target parameter, there are two ways. If you
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
In a multiple group analysis, modifiers that contain a single constant
must be replaced by a vector, having the same length as the number
of groups. The only exception are numerical constants (for fixing values): if
you provide only a single number, the same number will be used for all groups.
However, it is safer (and cleaner) to specify the same number of elements as
the number of groups. For example, if there are two groups:
\preformatted{
HS.model <- ' visual  =~ x1 + 0.5*x2 + c(0.6, 0.8)*x3
              textual =~ x4 + start(c(1.2, 0.6))*x5 + x6
              speed   =~ x7 + x8 + c(x9.group1, x9.group2)*x9 '
}
In this example, the factor loading of the `x2' indicator is fixed to the
value 0.5 for all groups. However, the factor loadings of the `x3' indicator 
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
\references{
Yves Rosseel (2012). lavaan: An R Package for Structural Equation
Modeling. Journal of Statistical Software, 48(2), 1-36. URL
http://www.jstatsoft.org/v48/i02/.}
