# export `lavaan' user model description to third-party software
# 

export <- function(object, target="mplus", file=NULL) {

    target <- tolower(target)

    if(class(object) == "lavaan") {
        user <- object@ParTable
    } else if(is.list(object)) {
        user <- object
    } else {
        stop("lavaan ERROR: object must be a `lavaan' object or a lavaan user object")
    }

    if(target == "mplus") {
        export.mplus(user, file=file)
    } else if(target == "lisrel") {
        export.lisrel(user, file=file)
    } else if(target == "eqs") {
        export.eqs(user, file=file)
    } else if(target == "sem") {
        export.sem(user, file=file)
    } else if(target == "openmx") {
        export.openmx(user, file=file)
    } else {
        stop("lavaan ERROR: target", target, "has not been implemented yet")
    }
}

export.mplus <- function(user, file=NULL) {
    stop("this function needs revision")
}

export.lisrel <- function(user, file=NULL) {
    stop("this function needs revision")
}

export.eqs <- function(user, file=NULL) {
    stop("this function needs revision")
}

export.sem <- function(user, file=NULL) {
    stop("this function needs revision")
}

export.openmx <- function(object, file=NULL) {
    stop("this function needs revision")
}

