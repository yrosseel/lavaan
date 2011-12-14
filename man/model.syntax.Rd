\name{model.syntax}
\alias{model.syntax}
\alias{lavaanify}
\alias{lavaanNames}
\title{The Lavaan Model Syntax}
\description{
The lavaan model syntax describes a latent variable model. The 
function \code{lavaanify} turns it into a list that represents the full
model as specified by the user.} 
\usage{
lavaanify(model.syntax = NULL, meanstructure = FALSE, int.ov.free = FALSE, 
    int.lv.free = FALSE, orthogonal = FALSE, std.lv = FALSE, fixed.x = TRUE, 
    constraints = NULL, auto = FALSE, model.type = "sem",
    auto.fix.first = FALSE, auto.fix.single = FALSE,  auto.var = FALSE, 
    auto.cov.lv.x = FALSE, auto.cov.y = FALSE, ngroups = 1L, group.equal = NULL,
    group.partial = NULL, debug = FALSE, warn = TRUE, as.data.frame. = TRUE)

lavaanNames(object, type = "ov", group = NULL)
}
\arguments{
\item{model.syntax}{The model syntax specifying the model. See
    details for more information.}
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
\item{fixed.x}{If \code{TRUE}, the exogenous `x' covariates are considered
    fixed variables and the means, variances and covariances of these variables
    are fixed to their sample values. If \code{FALSE}, they are considered
    random, and the means, variances and covariances are free parameters.}
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
\item{warn}{If \code{TRUE}, some (possibly harmless) warnings are printed
    out during the iterations.}
\item{as.data.frame.}{If \code{TRUE}, return the list of model parameters
as a data.frame.}
\item{debug}{If \code{TRUE}, debugging information is printed out.}
\item{object}{Either a list containing the user-specificed model, as 
returned by \code{lavaanify}, or an object of class 
\code{\linkS4class{lavaan}}.}
\item{type}{Only used in the function \code{lavaanNames}. If type contains 
\code{"ov"}, only observed variable names are returned. If type contains
\code{"lv"}, only latent variable names are returned. The \code{"ov.x"}
and \code{"lv.x"} types return exogenous variables only. The \code{"ov.y"}
and \code{"lv.y"} types return dependent variables only (in the regression
sense, excluding the indicators of latent variables). The
\code{"ov.nox"} type returns all observed variables, except the exogenous
ones.}
\item{group}{Only used in the function \code{lavaanNames}. If \code{NULL},
all groups (if any) are used. If an integer (vector), only names from those
groups are extracted. The group numbers are found in the \code{group}
column of the parameter table.}
}
\details{
The model syntax consists of one or more formula-like expressions, each one
describing a specific part of the model. The model syntax can be read from
a file (using \code{\link{readLines}}), or can be specified as a literal
string enclosed by single quotes as in the example below. 
\preformatted{myModel <- '
  # latent variable definitions
    f1 =~ y1 + y2 + y3
    f2 =~ y4 + y5 + y6
    f3 =~ y7 + y8 + 
          y9 + y10 
    f4 =~ y11 + y12 + y13

  ! this is also a comment
 
  # regressions
    f1 ~ f3 + f4
    f2 ~ f4
    y1 + y2 ~ x1 + x2 + x3

  # (co)variances
    y1 ~~ y1
    y2 ~~ y4 + y5 
    f1 ~~ f2

  # intercepts
    f1 ~ 1; y5 ~ 1
'
}
Blank lines and comments can be used in between the formulas, and formulas can
be split over multiple lines. Both the sharp (#) and the exclamation (!) 
characters can be used to start a comment. Multiple formulas can be placed
on a single line if they are separated by a semicolon (;).

There can be four types of formula-like expressions in the model syntax:
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
}
Usually, only a single variable name appears on the left side of an
operator. However, if multiple variable names are specified, 
separated by the \code{"+"} operator, the formula is repeated for each
element on the left side (as for example in the third regression 
formula in the example above).

In the right-hand side of these formula-like expressions, each element can be
multiplied (using the \code{"*"} operator) by either a numeric constant,
an expression resulting in a numeric constant, or one
of three special functions: \code{start()}, \code{label()} and \code{equal()}.
This provides the user with a mechanism to fix parameters, to provide
alternative starting values, to label the parameters, and to define equality
constraints among model parameters. This is explained in more detail in the
following sections.
}

\section{Fixing parameters}{
It is often desirable to fix a model parameter that is otherwise (by default)
free. Any parameter in a model can be fixed by using
the pre-multiplication mechanism. Here are some examples:
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
User-provided starting values can be given by using the same pre-multiplication
mechanism. But the numeric constant is now the argument
of a special function \code{start()}. For example:
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
specify its own labels for specific parameters by using the \code{label()}
modifier in a formula. For example, in the formula \code{f1 =~ x1 + x2 +
label("mylabel")*x3}, the parameter corresponding with the factor loading of
\code{x3} will be named \code{"mylabel"} instead of the default name
\code{"f1=~x3"}. 

To constrain a parameter
to be equal to another target parameter, the pre-multiplication mechanism can
be used using the special function \code{equal()}. The argument of
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
    y2a ~ label("int2") * 1
    y2b ~ equal("int2") * 1
    y2c ~ equal("int2") * 1
'
}
}

\section{Multiple groups}{
In a multiple group analysis, you can still use the pre-multiplication 
techniques, but the single argument is now replaced by a vector of 
arguments, one for each group. If it is a single value, the same value
is used for all groups. For example, if there are two groups:
\preformatted{HS.model <- ' visual  =~             x1 +
                                 0.5*x2 +
                         c(0.6, 0.8)*x3

              textual =~                    x4 +
                         start(c(1.2, 0.6))*x5 +
                                            x6

              speed   =~                        x7 +
                                                x8 +
                         label(c("x9.group1",
                                 "x9.group2"))*x9 '
}
In this example, the factor loading of the `x2' indicator is fixed to the
value 0.5 for all groups. However, the factor loadings of the `x3' indicator 
are fixed to 0.6 and 0.8 for group 1 and group 2 respectively. The same
logic is used for all modifiers.
}

\section{Multiple modifiers}{
In the model syntax, you can specify a variable more than once on the right hand
side of an operator; therefore, several `modifiers' can be applied
simultaneously; for example, if you want to fix the value of a parameter and
also label that parameter, you can use something like:
\preformatted{ f1 =~ x1 + x2 + 4*x3 + label("x3.loading")*x3}
}
