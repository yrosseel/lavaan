# factor score regression


# this function will call the sem/lavaan function several times
fsr <-  function(model = NULL, data = NULL, model.type = "sem", ...) {
   
    # we need full data
    if(is.null(data)) {
        stop("lavaan ERROR: full data is required for factor score regression")
    }

    # fit model without data, to get parameter table
    FIT <- do.call(model.type, data = NULL, do.fit = FALSE, ...)

    # check parameter table

    FIT
}

