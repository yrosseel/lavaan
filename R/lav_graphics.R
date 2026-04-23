# small functions to do something useful with the common
# plot commands

# suggested by JEB
lav_lavaan_pairs <- function(x, group = 1L, ...) {
  m_x <- x@Data@X[[group]]
  colnames(m_x) <- x@Data@ov.names[[group]]
  pairs(m_x, ...)
}
