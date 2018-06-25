# lav_model utility functions

# initial version: YR 25/03/2009: `methods' for the Model class
# - YR 14 Jan 2014: rename object -> lavmodel, all functions as lav_model_*

lav_model_get_parameters <- function(lavmodel = NULL, GLIST = NULL,
                                     type = "free", extra = TRUE) {

    # type == "free": only non-redundant free parameters (x)
    # type == "user": all parameters listed in User model

    # state or final?
    if(is.null(GLIST)) GLIST <- lavmodel@GLIST

    if(type == "free") {
        N <- lavmodel@nx.free
    #} else if(type == "unco") {
    #    N <- lavmodel@nx.unco
    } else if(type == "user") {
        N <- lavmodel@nx.user
    }
    x <- numeric(N)

    for(mm in 1:length(lavmodel@GLIST)) {
        if(type == "free") {
            m.idx <- lavmodel@m.free.idx[[mm]]
            x.idx <- lavmodel@x.free.idx[[mm]]
        #} else if(type == "unco") {
        #    m.idx <- lavmodel@m.unco.idx[[mm]]
        #    x.idx <- lavmodel@x.unco.idx[[mm]]
        } else if(type == "user") {
            m.idx <- lavmodel@m.user.idx[[mm]]
            x.idx <- lavmodel@x.user.idx[[mm]]
        }
        x[x.idx] <- GLIST[[mm]][m.idx]
    }

    if(type == "user" && extra && sum(lavmodel@x.def.idx,
                                      lavmodel@x.ceq.idx,
                                      lavmodel@x.cin.idx) > 0L) {
        # we need 'free' x
        x.free <- lav_model_get_parameters(lavmodel = lavmodel, GLIST = GLIST,
                                           type = "free")
        if(length(lavmodel@x.def.idx) > 0L) {
            x[lavmodel@x.def.idx] <- lavmodel@def.function(x.free)
        }
        if(length(lavmodel@x.ceq.idx) > 0L) {
            x[lavmodel@x.ceq.idx] <- lavmodel@ceq.function(x.free)
        }
        if(length(lavmodel@x.cin.idx) > 0L) {
            x[lavmodel@x.cin.idx] <- lavmodel@cin.function(x.free)
        }
    }

    x
}

# warning: this will make a copy of lavmodel
lav_model_set_parameters <- function(lavmodel = NULL, x = NULL) {

    tmp <- lavmodel@GLIST
    for(mm in 1:length(lavmodel@GLIST)) {
        m.free.idx <- lavmodel@m.free.idx[[mm]]
        x.free.idx <- lavmodel@x.free.idx[[mm]]
        tmp[[mm]][m.free.idx] <- x[x.free.idx]
    }

    # categorical? set categorical theta elements (if any)
    if(lavmodel@categorical) {
        nmat <- lavmodel@nmat
        if(lavmodel@representation == "LISREL") {
            for(g in 1:lavmodel@nblocks) {
                # which mm belong to group g?
                mm.in.group <- 1:nmat[g] + cumsum(c(0L,nmat))[g]

                if(lavmodel@estimator %in% c("WLS","DWLS","ULS","PML")) {
                    if(lavmodel@parameterization == "delta") {
                        tmp[mm.in.group] <-
                        setResidualElements.LISREL(MLIST = tmp[mm.in.group],
                            num.idx = lavmodel@num.idx[[g]],
                            ov.y.dummy.ov.idx = lavmodel@ov.y.dummy.ov.idx[[g]],
                            ov.y.dummy.lv.idx = lavmodel@ov.y.dummy.lv.idx[[g]])
                    } else if(lavmodel@parameterization == "theta") {
                        tmp[mm.in.group] <-
                        setDeltaElements.LISREL(MLIST = tmp[mm.in.group],
                            num.idx = lavmodel@num.idx[[g]])
                    }
                } else if(lavmodel@estimator %in% c("MML", "FML")) {
                  #  ttt <- diag(tmp[mm.in.group]$theta)
                  #  diag(tmp[mm.in.group]$theta) <- as.numeric(NA)
                  #  if(length(lavmodel@num.idx[[g]]) > 0L) {
                  #      diag(tmp[mm.in.group]$theta)[ lavmodel@num.idx[[g]] ] <-
                  #          ttt[ lavmodel@num.idx[[g]] ]
                  #  }
                }
            }
        } else {
            cat("FIXME: deal with theta elements in the categorical case")
        }
    }

    lavmodel@GLIST <- tmp

    lavmodel
}

# create a standalone GLIST, filled with (new) x values
# (avoiding a copy of lavmodel)
lav_model_x2GLIST <- function(lavmodel = NULL, x = NULL,
                              type = "free", setDelta = TRUE,
                              m.el.idx = NULL, x.el.idx = NULL) {

    GLIST <- lavmodel@GLIST
    for(mm in 1:length(GLIST)) {
        # skip empty matrix
        if(nrow(GLIST[[mm]]) == 0L)
            next
        if(type == "free") {
            M.EL.IDX <- lavmodel@m.free.idx[[mm]]
            X.EL.IDX <- lavmodel@x.free.idx[[mm]]
        } else if(type == "full") {
            if(lavmodel@isSymmetric[mm]) {
                N <- ncol(GLIST[[mm]])
                M.EL.IDX <- lav_matrix_vech_idx(N)
            } else {
                M.EL.IDX <- seq_len(length(GLIST[[mm]]))
            }
            X.EL.IDX <- seq_len(length(m.el.idx))
            if(mm > 1) X.EL.IDX <- X.EL.IDX + sum(lavmodel@mmSize[1:(mm-1)])
        } else if(type == "custom") {
            # nothing to do, m.el.idx and x.el.idx should be given
            M.EL.IDX <- m.el.idx[[mm]]
            X.EL.IDX <- x.el.idx[[mm]]
        }

        # assign
        GLIST[[mm]][M.EL.IDX] <- x[X.EL.IDX]

        # make symmetric (if full)
        if(type == "full" && lavmodel@isSymmetric[mm]) {
            T <- t(GLIST[[mm]])
            GLIST[[mm]][upper.tri(GLIST[[mm]])] <- T[upper.tri(T)]
        }
    }

    # theta parameterization: delta must be reset!
    if(lavmodel@categorical && setDelta &&
       lavmodel@parameterization == "theta") {
        nmat <- lavmodel@nmat
        for(g in 1:lavmodel@nblocks) {
            # which mm belong to group g?
            mm.in.group <- 1:nmat[g] + cumsum(c(0L,nmat))[g]
            GLIST[mm.in.group] <-
                setDeltaElements.LISREL(MLIST = GLIST[mm.in.group],
                    num.idx = lavmodel@num.idx[[g]])
        }
    }

    GLIST
}

# backwards compatibility
# getModelParameters <- lav_model_get_parameters
# setModelParameters <- lav_model_set_parameters
# x2GLIST            <- lav_model_x2GLIST
