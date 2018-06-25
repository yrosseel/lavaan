# prelav
#
# a program not unlike prelis (unfinished)
#
# YR 19 Sept 2013

prelav <- function(object = NULL, ordered = NULL, ov.names.x = NULL,
                   group = NULL, missing = "pairwise",
                   output = list(MA ="none", # moment matrix
                                 SR = FALSE, # transformed raw data
                                 RA = FALSE, # transformed raw data
                                 SA = FALSE, # asymptotic covariance matrix
                                 AC = FALSE, # asymptotic covariance matrix
                                 SV = FALSE, # asymptotic variances
                                 TH = FALSE, # thresholds
                                 ME = FALSE, # means
                                 ND = 3L,    # number of decimals
                                 PK = FALSE, # Mardia's (1970) mult kurtosis
                                 WP = FALSE, # wide print
                                 XB = FALSE, # omit bivariate freq tables
                                 XT = FALSE, # omit omit test statistics
                                 XM = FALSE  # omit tests of mult normality
                                ),
                   mimic="LISREL") {

    # empty output?
    if(length(output) == 0L) {
        return( list() )
    } else {
        OU <- toupper(substr(names(output), 1, 2))
    }

    # parse data
    NAMES <- names(object)
    if(!is.null(group)) {
        NAMES <- NAMES[- match(group, NAMES)]
    }
    lav.data <- lavData(data = object, group = group,
                        ov.names = NAMES, ordered = ordered,
                        ov.names.x = ov.names.x,
                        lavoptions = list(missing = missing))

    lav.stats <- lav_samplestats_from_data(lavdata       = lav.data,
                                           missing       = missing,
                                           rescale       = FALSE,
                                           estimator     = "ML",
                                           mimic         = mimic,
                                           meanstructure = TRUE,
                                           conditional.x = FALSE,
                                           group.w.free  = FALSE,
                                           missing.h1    = FALSE,
                                           WLS.V         = NULL,
                                           NACOV         = NULL,
                                           ridge         = 1e-5,
                                           debug         = FALSE,
                                           verbose       = FALSE)

    out <- list(lav.data=lav.data, lav.stats=lav.stats, OU=OU)


    # output
    class(out) <- c("prelav", "list")
    out
}

# S3 method
print.prelav <- function(x, ..., nd=3) {

    # shorcuts
    lav.data <- x$lav.data
    lav.stats <- x$lav.stats
    ngroups <- lav.data@ngroups

    # header
    version <- read.dcf(file=system.file("DESCRIPTION", package="lavaan"),
                        fields="Version")
    cat("This is prelav ", version, ".\n", sep="")

    # data information + missingness
    cat("\n")
    print(lav.data)

    # univariate information
    cat("\n")
    cat("Univariate information:\n")
    # varTable!
    print(as.data.frame(lav.data@ov))

    invisible(x)
}
