# generate parameter table for an independence model
# YR - 12 Sep 2017: special case of lav_partable_unrestricted()
lav_partable_independence <- function(lavobject      = NULL,
                                      # if no object is available,
                                      lavdata        = NULL,
                                      lavpta         = NULL,
                                      lavoptions     = NULL,
                                      lavsamplestats = NULL,
                                      # optional user-provided sample stats
                                      sample.cov     = NULL,
                                      sample.mean    = NULL,
                                      sample.slopes  = NULL,
                                      sample.th      = NULL,
                                      sample.th.idx  = NULL,
                                      sample.cov.x   = NULL,
                                      sample.mean.x  = NULL) {

    lav_partable_indep_or_unrestricted(lavobject = lavobject,
        lavdata = lavdata, lavpta = lavpta, lavoptions = lavoptions,
        lavsamplestats = lavsamplestats, sample.cov = sample.cov,
        sample.mean = sample.mean , sample.slopes = sample.slopes,
        sample.th = sample.th, sample.th.idx = sample.th.idx,
        independent = TRUE)
}
