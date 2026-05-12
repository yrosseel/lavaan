# Test that lavInspect() with what="icc" works for clustered data
# This test verifies the fix for issue #502

library(lavaan)

# Demo from ?Demo.twolevel
model <- '
  level: 1
    fw =~ y1 + y2 + y3
    fw ~ x1 + x2 + x3
  level: 2
    fb =~ y1 + y2 + y3
    fb ~ w1 + w2
'
fit <- sem(model, data = Demo.twolevel, cluster = "cluster")

# This should work (and it does in current lavaan)
icc_values <- lavInspect(fit, "icc")
print("ICC values computed successfully:")
print(icc_values)

# Verify it's documented
help_text <- help("lavInspect", package = "lavaan")
if (length(grep("icc", capture.output(help_text), ignore.case = TRUE)) == 0) {
  stop("'icc' is not documented in lavInspect help page")
} else {
  cat("\n'icc' is now documented in lavInspect\n")
}
