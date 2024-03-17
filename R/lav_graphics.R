# small functions to do something useful with the common
# plot commands

# suggested by JEB
pairs.lavaan <- function(x, group = 1L, ...) {
  X <- x@Data@X[[group]]
  colnames(X) <- x@Data@ov.names[[group]]
  pairs(X, ...)
}
