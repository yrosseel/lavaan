# Model-methods

#
# initial version: YR 25/03/2009: `methods' for the Model class

lav_model_get_parameters <- function(object, GLIST=NULL, type="free",
                               extra=TRUE) {

    # type == "free": only non-redundant free parameters (x)
    # type == "unco": all free parameters (including constrained ones)
    # type == "user": all parameters listed in User model

    # state or final?
    if(is.null(GLIST)) GLIST <- object@GLIST

    if(type == "free") {
        N <- object@nx.free
    } else if(type == "unco") {
        N <- object@nx.unco
    } else if(type == "user") {
        N <- object@nx.user
    }
    x <- numeric(N)

    for(mm in 1:length(object@GLIST)) {
        if(type == "free") {
            m.idx <- object@m.free.idx[[mm]]
            x.idx <- object@x.free.idx[[mm]]
        } else if(type == "unco") { 
            m.idx <- object@m.unco.idx[[mm]]
            x.idx <- object@x.unco.idx[[mm]]
        } else if(type == "user") {
            m.idx <- object@m.user.idx[[mm]]
            x.idx <- object@x.user.idx[[mm]]
        }
        x[x.idx] <- GLIST[[mm]][m.idx]
    }

    if(type == "user" && extra && sum(object@x.def.idx,
                                      object@x.ceq.idx, 
                                      object@x.cin.idx) > 0L) {
        # we need 'free' x
        x.free <- lav_model_get_parameters(object, GLIST=GLIST, type="free")
        if(length(object@x.def.idx) > 0L) {
            x[object@x.def.idx] <- object@def.function(x.free)
        }
        if(length(object@x.ceq.idx) > 0L) {
            x[object@x.ceq.idx] <- object@ceq.function(x.free)
        }
        if(length(object@x.cin.idx) > 0L) {
            x[object@x.cin.idx] <- object@cin.function(x.free)
        }
    }

    x
}

# warning: this will make a copy of object
lav_model_set_parameters <- function(object, x=NULL, estimator = "ML") {

    tmp <- object@GLIST
    for(mm in 1:length(object@GLIST)) {
        m.free.idx <- object@m.free.idx[[mm]]
        x.free.idx <- object@x.free.idx[[mm]]
        tmp[[mm]][m.free.idx] <- x[x.free.idx]
    }

    # categorical? set categorical theta elements (if any)
    if(object@categorical) {
        nmat <- object@nmat
        if(object@representation == "LISREL") {
            for(g in 1:object@ngroups) {
                # which mm belong to group g?
                mm.in.group <- 1:nmat[g] + cumsum(c(0L,nmat))[g]

                if(estimator %in% c("WLS","DWLS","ULS","PML")) {
                    if(object@parameterization == "delta") {
                        tmp[mm.in.group] <- 
                        setResidualElements.LISREL(MLIST = tmp[mm.in.group],
                            num.idx = object@num.idx[[g]],
                            ov.y.dummy.ov.idx = object@ov.y.dummy.ov.idx[[g]],
                            ov.y.dummy.lv.idx = object@ov.y.dummy.lv.idx[[g]])
                    } else if(object@parameterization == "theta") {
                        tmp[mm.in.group] <-
                        setDeltaElements.LISREL(MLIST = tmp[mm.in.group],
                            num.idx = object@num.idx[[g]])
                    }
                } else if(estimator %in% c("MML", "FML")) {
                    ttt <- diag(tmp[mm.in.group]$theta)
                    diag(tmp[mm.in.group]$theta) <- as.numeric(NA)
                    if(length(object@num.idx[[g]]) > 0L) {
                        diag(tmp[mm.in.group]$theta)[ object@num.idx[[g]] ] <-
                            ttt[ object@num.idx[[g]] ]
                    }
                }
            }
        } else {
            cat("FIXME: deal with theta elements in the categorical case")
        }
    }

    object@GLIST <- tmp

    object
}

# create a standalone GLIST, filled with (new) x values
# (avoiding a copy of object)
lav_model_x2GLIST <- function(object, x=NULL, type="free", setDelta = TRUE) {

    GLIST <- object@GLIST
    for(mm in 1:length(GLIST)) {
        if(type == "free") {
            m.el.idx <- object@m.free.idx[[mm]]
            x.el.idx <- object@x.free.idx[[mm]]
        } else if(type == "full") {
            if(object@isSymmetric[mm]) {
                N <- ncol(GLIST[[mm]])
                m.el.idx <- vech.idx(N)
            } else {
                m.el.idx <- 1:length(GLIST[[mm]])
            }
            x.el.idx <- 1:length(m.el.idx)
            if(mm > 1) x.el.idx <- x.el.idx + sum(object@mmSize[1:(mm-1)])
        }

        # assign
        GLIST[[mm]][m.el.idx] <- x[x.el.idx]

        # make symmetric (if full)
        if(type == "full" && object@isSymmetric[mm]) {
            T <- t(GLIST[[mm]])
            GLIST[[mm]][upper.tri(GLIST[[mm]])] <- T[upper.tri(T)]
        }
    }

    # theta parameterization: delta must be reset!
    if(setDelta && object@parameterization == "theta") {
        nmat <- object@nmat
        for(g in 1:object@ngroups) {
            # which mm belong to group g?
            mm.in.group <- 1:nmat[g] + cumsum(c(0L,nmat))[g]
            GLIST[mm.in.group] <-
                setDeltaElements.LISREL(MLIST = GLIST[mm.in.group],
                    num.idx = object@num.idx[[g]])
        }
    }

    GLIST
}

