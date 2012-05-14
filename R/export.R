# export `lavaan' lav model description to third-party software
# 

export <- function(object, target="mplus", file=NULL) {

    target <- tolower(target)

    if(class(object) == "lavaan") {
        lav <- object@ParTable
    } else if(is.list(object)) {
        lav <- object
    } else {
        stop("lavaan ERROR: object must be of class `lavaan' or a parTable")
    }

    if(target == "lavaan") {
        model <- export_lavaan(lav)
    } else if(target == "mplus") {
        model <- export_Mplus(lav)
    } else if(target == "lisrel") {
        model <- export_lisrel(lav)
    } else if(target == "eqs") {
        model <- export_eqs(lav)
    } else if(target == "sem") {
        model <- export_sem(lav)
    } else if(target == "openmx") {
        model <- export_openmx(lav)
    } else {
        stop("lavaan ERROR: target", target, "has not been implemented yet")
    }
}

## FIXME: this is completely UNFINISHED (just  used to quickly get something)
export_lavaan(lav) {
    N <- length(lav$lhs)
    out <- character(N)
    for(i in 1:length(lav$lhs)) {
        out[i]
    }
    
}

export_mplus <- function(lav) {
    stop("this function needs revision")
}

export_lisrel <- function(lav) {
    stop("this function needs revision")
}

export.eqs <- function(lav) {
    stop("this function needs revision")
}

export.sem <- function(lav) {
    stop("this function needs revision")
}

export.openmx <- function(lav) {
    stop("this function needs revision")
}

