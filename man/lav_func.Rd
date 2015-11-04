\name{lav_func}
\alias{lav_func_gradient_complex}
\alias{lav_func_gradient_simple}
\alias{lav_func_jacobian_complex}
\alias{lav_func_jacobian_simple}
\title{Utility Functions: Gradient and Jacobian}
\description{Utility functions for computing the gradient of a scalar-valued 
function or the Jacobian of a vector-valued function by numerical approximation.}
\usage{
lav_func_gradient_complex(func, x, h = .Machine$double.eps, ..., 
                          check.scalar = TRUE, fallback.simple = TRUE)
lav_func_jacobian_complex(func, x, h = .Machine$double.eps, ...,
                          fallback.simple = TRUE)

lav_func_gradient_simple(func, x, h = sqrt(.Machine$double.eps), ...,
                         check.scalar = TRUE)
lav_func_jacobian_simple(func, x, h = sqrt(.Machine$double.eps), ...)
}
\arguments{
\item{func}{A real-valued function returning a numeric scalar or a 
numeric vector.}
\item{x}{A numeric vector: the point(s) at which the gradient/Jacobian
of the function should be computed.}
\item{h}{Numeric value representing a small change in `x' when
computing the gradient/Jacobian.}
\item{...}{Additional arguments to be passed to the function `func'.}
\item{check.scalar}{Logical. If TRUE, check if the function is scalar-valued.}
\item{fallback.simple}{Logical. If TRUE, and the function evaluation fails,
we call the corresponding simple (non-complex) method instead.}
}
\details{
The complex versions use complex numbers to gain more precision, while
retaining the simplicity (and speed) of the simple forward method (see
references). These
functions were added to lavaan (around 2012) when the complex functionality
was not part of the numDeriv package. They were used internally, and made
public in 0.5-17 per request of other package developers.
}
\references{
Squire, W. and Trapp, G. (1998). Using Complex Variables to Estimate 
Derivatives of Real Functions. SIAM Review, 40(1), 110-112.
}
\examples{
# very accurate complex method
lav_func_gradient_complex(func = exp, x = 1) - exp(1)

# less accurate forward method
lav_func_gradient_simple(func = exp, x = 1) - exp(1)

# very accurate complex method
diag(lav_func_jacobian_complex(func = exp, x = c(1,2,3))) - exp(c(1,2,3))

# less accurate forward method
diag(lav_func_jacobian_simple(func = exp, x = c(1,2,3))) - exp(c(1,2,3))
}
