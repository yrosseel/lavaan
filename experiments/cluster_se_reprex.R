library(lavaan)
set.seed(1234)
LETTERS702 <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))

# Generate a clustered population of 5 binary items
Nclust <- 500
clust_sizes <- sample(10:30, Nclust, replace = TRUE)
N <- sum(clust_sizes)
clust <- LETTERS702[rep(1:Nclust, clust_sizes)]

# Underlying variable approach
pop.model <- '
eta1 =~ 0.8 * y1 + 0.7 * y2 + 0.47 * y3 + 0.38 * y4 + 0.34 * y5
'
tau <- c(1.43, 0.55, 0.13, 0.72, 1.13)
ystar <- simulateData(pop.model, sample.nobs = N, standardized = TRUE)

# Generate binary data based on thresholds
Pop <- 1 * (ystar > matrix(tau, nrow = N, ncol = 5, byrow = TRUE))
colnames(Pop) <- paste0("y", 1:5)
Pop <- lapply(as.data.frame(Pop), \(x) ordered(x, levels = c(0, 1)))
Pop$clust <- clust
Pop <- as.data.frame(Pop)

# Perform a clustered sample (nc = 50)
samp_prob <- clust_sizes / sum(clust_sizes)  # PPS
names(samp_prob) <- LETTERS702[1:Nclust]
clust_samp <- sample(unique(clust), 50, prob = samp_prob)
Data <- Pop[Pop$clust %in% clust_samp, ]
Data$wt <- 1 / samp_prob[Data$clust]

# Fit model
mod <- '
eta1 =~ y1 + y2 + y3 + y4 + y5
'
fit <- sem(model = mod, data = Data, estimator = "PML", std.lv = TRUE,
           sampling.weights = "wt", cluster = "clust")
summary(fit)
