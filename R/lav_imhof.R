imhof <- \(x, lambda) {
  theta <- \(u, x, lambda) 0.5 * (sum(atan(lambda * u)) - x * u)
  rho <- \(u, lambda) exp(1 / 4 * sum(log(1 + lambda^2 * u^2)))
  integrand <- Vectorize(\(u) {
    sin(theta(u, x, lambda)) / (u * rho(u, lambda))
  })
  z <- integrate(integrand, lower = 0, upper = Inf)$value
  0.5 + z / pi
}

# k <- 50
# lambda <- runif(k)
# x <- k/2
# imhof(x, lambda)
# CompQuadForm::imhof(x, lambda)

# microbenchmark::microbenchmark(imhof(x, lambda), CompQuadForm::imhof(x, lambda))
