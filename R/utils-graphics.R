# small functions to do something useful with the common
# plot commands

# suggested by JEB
pairs.lavaan <- function(x, group=1L, ...) {
        pairs(x@Data[[group]])
}
