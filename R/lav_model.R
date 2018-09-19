# constructor of the matrix lavoptions$representation
#
# initial version: YR 22/11/2010
# - YR 14 Jan 2014: moved to lav_model.R
# - YR 18 Nov 2014: more efficient handling of linear equality constraints
# - YR 02 Dec 2014: allow for bare-minimum parameter tables
# - YR 25 Jan 2017: collect options in lavoptions

# construct MATRIX lavoptions$representation of the model
lav_model <- function(lavpartable      = NULL,
                      lavoptions       = NULL,
                      th.idx           = list(),
                      cov.x            = list(),
                      mean.x           = list()) { # for conditional.x only
                                                   # (not really needed,
                                                   #  but as a failsafe)

    # handle bare-minimum partables
    lavpartable <- lav_partable_complete(lavpartable)

    # global info from user model
    nblocks <- lav_partable_nblocks(lavpartable)
    ngroups <- lav_partable_ngroups(lavpartable)
    meanstructure <- any(lavpartable$op == "~1")
    categorical <- any(lavpartable$op == "|")
    if(categorical) {
        meanstructure <- TRUE

        # handle th.idx if length(th.idx) != nblocks
        if(nblocks != length(th.idx)) {
            th.idx <- rep(th.idx, each = nblocks)
        }

    }
    group.w.free <- any(lavpartable$lhs == "group" & lavpartable$op == "%")
    multilevel <- FALSE
    if(!is.null(lavpartable$level)) {
        nlevels <- lav_partable_nlevels(lavpartable)
        if(nlevels > 1L) {
            multilevel <- TRUE
        }
    }

    # handle variable definitions and (in)equality constraints
    CON <- lav_constraints_parse(partable = lavpartable,
                                 constraints = NULL,
                                 debug = lavoptions$debug)

    # handle *linear* equality constraints special
    if(CON$ceq.linear.only.flag) {
        con.jac <- CON$ceq.JAC
        con.lambda <- numeric(nrow(CON$ceq.JAC))
        attr(con.jac, "inactive.idx") <- integer(0L)
        attr(con.jac, "ceq.idx") <- seq_len( nrow(CON$ceq.JAC) )
    } else {
        con.jac <- matrix(0,0,0)
        con.lambda <- numeric(0)
    }

    # select model matrices
    if(lavoptions$representation == "LISREL") {
        REP <- representation.LISREL(lavpartable, target = NULL, extra = TRUE)
    } else {
        stop("lavaan ERROR: only representation \"LISREL\" has been implemented.")
    }
    if(lavoptions$debug) print(REP)

    # FIXME: check for non-existing parameters
    bad.idx <- which(REP$mat == "" &
                     !lavpartable$op %in% c("==","<",">",":="))

    if(length(bad.idx) > 0L) {

        label <- paste(lavpartable$lhs[bad.idx[1]],
                       lavpartable$op[bad.idx[1]],
                       lavpartable$rhs[bad.idx[1]], sep = " ")
        stop("lavaan ERROR: parameter is not defined: ", label)
    }

    # prepare nG-sized slots
    nG <- sum(unlist(attr(REP, "mmNumber")))
    GLIST <- vector(mode="list", nG)
    names(GLIST) <- unlist(attr(REP, "mmNames"))
    dimNames    <- vector(mode="list", length=nG)
    isSymmetric <- logical(nG)
    mmSize      <- integer(nG)

    m.free.idx <- m.user.idx <- vector(mode="list", length=nG)
    x.free.idx <- x.user.idx <- vector(mode="list", length=nG)

    # prepare nblocks-sized slots
    nvar <- integer(nblocks)
    nmat <- unlist(attr(REP, "mmNumber"))
    num.idx <- vector("list", length=nblocks)
    nexo <- integer(nblocks)
    ov.x.dummy.ov.idx <- vector(mode="list", length=nblocks)
    ov.x.dummy.lv.idx <- vector(mode="list", length=nblocks)
    ov.y.dummy.ov.idx <- vector(mode="list", length=nblocks)
    ov.y.dummy.lv.idx <- vector(mode="list", length=nblocks)

    offset <- 0L
    for(g in 1:nblocks) {

        # observed and latent variables for this block
        ov.names <-     lav_partable_vnames(lavpartable, "ov",     block = g)
        ov.names.nox <- lav_partable_vnames(lavpartable, "ov.nox", block = g)
        ov.names.x <-   lav_partable_vnames(lavpartable, "ov.x",   block = g)
        nexo[g] <- length(ov.names.x)
        ov.num <-       lav_partable_vnames(lavpartable, "ov.num", block = g)
        if(lavoptions$conditional.x) {
            nvar[g] <- length(ov.names.nox)
            num.idx[[g]] <- which(ov.names.nox %in% ov.num)
        } else {
            nvar[g] <- length(ov.names)
            num.idx[[g]] <- which(ov.names %in% ov.num)
        }

        # model matrices for this block
        mmNumber    <- attr(REP, "mmNumber")[[g]]
        mmNames     <- attr(REP, "mmNames")[[g]]
        mmSymmetric <- attr(REP, "mmSymmetric")[[g]]
        mmDimNames  <- attr(REP, "mmDimNames")[[g]]
        mmRows      <- attr(REP, "mmRows")[[g]]
        mmCols      <- attr(REP, "mmCols")[[g]]

        for(mm in 1:mmNumber) {

            # offset in GLIST
            offset <- offset + 1L

            # matrix size, symmetric, dimNames
            if(mmSymmetric[mm]) {
                N <- mmRows[mm]
                mm.size <- as.integer(N*(N+1)/2)
            } else {
                mm.size <- as.integer(mmRows[mm] * mmCols[mm])
            }
            mmSize[offset] <- mm.size
            isSymmetric[offset] <- mmSymmetric[mm]
            dimNames[[offset]] <- mmDimNames[[mm]]

            # select elements for this matrix
            idx <- which(lavpartable$block == g & REP$mat == mmNames[mm])

            # create empty `pattern' matrix
            # FIXME: one day, we may want to use sparse matrices...
            #        but they should not slow things down!
            tmp <- matrix(0L, nrow=mmRows[mm],
                              ncol=mmCols[mm])

            # 1. first assign free values only, to get vector index
            #    -> to be used in lav_model_objective
            tmp[ cbind(REP$row[idx], REP$col[idx]) ] <- lavpartable$free[idx]
            if(mmSymmetric[mm]) {
                # NOTE: we assume everything is in the UPPER tri!
                T <- t(tmp); tmp[lower.tri(tmp)] <- T[lower.tri(T)]
            }
            m.free.idx[[offset]] <-     which(tmp > 0)
            x.free.idx[[offset]] <- tmp[which(tmp > 0)]

            # 2. if equality constraints, unconstrained free parameters
            #    -> to be used in lav_model_gradient
            #if(CON$ceq.linear.only.flag) {
            #    tmp[ cbind(REP$row[idx],
            #               REP$col[idx]) ] <- lavpartable$unco[idx]
            #    if(mmSymmetric[mm]) {
            #        # NOTE: we assume everything is in the UPPER tri!
            #        T <- t(tmp); tmp[lower.tri(tmp)] <- T[lower.tri(T)]
            #    }
            #    m.unco.idx[[offset]] <-     which(tmp > 0)
            #    x.unco.idx[[offset]] <- tmp[which(tmp > 0)]
            #} else {
            #    m.unco.idx[[offset]] <- m.free.idx[[offset]]
            #    x.unco.idx[[offset]] <- x.free.idx[[offset]]
            #}

            # 3. general mapping between user and GLIST
            tmp[ cbind(REP$row[idx], REP$col[idx]) ] <- lavpartable$id[idx]
            if(mmSymmetric[mm]) {
                T <- t(tmp); tmp[lower.tri(tmp)] <- T[lower.tri(T)]
            }
            m.user.idx[[offset]] <-     which(tmp > 0)
            x.user.idx[[offset]] <- tmp[which(tmp > 0)]

            # 4. now assign starting/fixed values
            # create empty matrix
            # FIXME: again, we may want to use sparse matrices here...
            tmp <- matrix(0.0, nrow=mmRows[mm],
                               ncol=mmCols[mm])
            tmp[ cbind(REP$row[idx], REP$col[idx]) ] <- lavpartable$start[idx]
            if(mmSymmetric[mm]) {
                T <- t(tmp); tmp[lower.tri(tmp)] <- T[lower.tri(T)]
            }

            # 4b. override with cov.x (if conditional.x = TRUE)
            # new in 0.6-1
            # shouldn't be needed, if lavpartable$start contains cov.x values
            if(mmNames[mm] == "cov.x") {
                tmp <- cov.x[[g]]
            }
            # 4c. override with mean.x (if conditional.x = TRUE)
            # new in 0.6-1
            # shouldn't be needed, if lavpartable$start contains mean.x values
            if(mmNames[mm] == "mean.x") {
                tmp <- as.matrix(mean.x[[g]])
            }

            # representation specific stuff
            if(lavoptions$representation == "LISREL" && mmNames[mm] == "lambda") {
                ov.dummy.names.nox <- attr(REP, "ov.dummy.names.nox")[[g]]
                ov.dummy.names.x   <- attr(REP, "ov.dummy.names.x")[[g]]
                ov.dummy.names <- c(ov.dummy.names.nox, ov.dummy.names.x)
                # define dummy latent variables
                if(length(ov.dummy.names)) {
                    # in this case, lv.names will be extended with the dummys
                    LV.names <- mmDimNames$psi[[1]]
                    row..idx <- match(ov.dummy.names, ov.names)
                    col..idx <- match(ov.dummy.names, LV.names)
                    # Fix lambda values to 1.0
                    tmp[ cbind(row..idx, col..idx)] <- 1.0

                    ov.x.dummy.ov.idx[[g]] <- match(ov.dummy.names.x,ov.names)
                    ov.x.dummy.lv.idx[[g]] <- match(ov.dummy.names.x,LV.names)
                    ov.y.dummy.ov.idx[[g]] <- match(ov.dummy.names.nox,ov.names)
                    ov.y.dummy.lv.idx[[g]] <- match(ov.dummy.names.nox,LV.names)
                }
            }

            # representation specific
            if(lavoptions$representation == "LISREL" && mmNames[mm] == "delta") {
                # only categorical values are listed in the lavpartable
                # but all remaining values should be 1.0
                idx <- which(tmp[,1L] == 0.0)
                tmp[idx,1L] <- 1.0
            }

            # assign matrix to GLIST
            GLIST[[offset]] <- tmp
        } # mm
    } # g

    # fixed.x parameters?
    #fixed.x <- any(lavpartable$exo > 0L & lavpartable$free == 0L)
    #if(categorical) {
    #    fixed.x <- TRUE
    #}

    # dirty hack to mimic MUML
    if(!is.null(lavoptions$tech.muml.scale)) {
        warning("lavaan WARNING: using muml scale in group 2")

        # find matrix
        lambda.idx <- which(names(GLIST) == "lambda")[2L]

        # find rows/cols
        B.names <- paste0("b", ov.names) ## ad-hoc assumption!!!
        COLS <-  match(B.names, LV.names)
        ROWS <- seq_len(nvar[2])
        stopifnot(length(COLS) == length(ROWS))
        GLIST[[ lambda.idx ]][ cbind(ROWS, COLS) ] <- lavoptions$tech.muml.scale
    }

    # which free parameters are observed variances?
    ov.names <- vnames(lavpartable, "ov")
    x.free.var.idx <- lavpartable$free[ lavpartable$free &
                                        #!duplicated(lavpartable$free) &
                                        lavpartable$lhs %in% ov.names &
                                        lavpartable$op == "~~" &
                                        lavpartable$lhs == lavpartable$rhs ]

    Model <- new("lavModel",
                 GLIST=GLIST,
                 dimNames=dimNames,
                 isSymmetric=isSymmetric,
                 mmSize=mmSize,
                 representation=lavoptions$representation,
                 meanstructure=meanstructure,
                 categorical=categorical,
                 multilevel=multilevel,
                 link=lavoptions$link,
                 nblocks=nblocks,
                 ngroups=ngroups, # breaks rsem????
                 group.w.free=group.w.free,
                 nmat=nmat,
                 nvar=nvar,
                 num.idx=num.idx,
                 th.idx=th.idx,
                 nx.free=max(lavpartable$free),
                 #nx.unco=max(lavpartable$unco),
                 nx.user=max(lavpartable$id),
                 m.free.idx=m.free.idx,
                 x.free.idx=x.free.idx,
                 x.free.var.idx=x.free.var.idx,
                 #m.unco.idx=m.unco.idx,
                 #x.unco.idx=x.unco.idx,
                 m.user.idx=m.user.idx,
                 x.user.idx=x.user.idx,
                 x.def.idx=which(lavpartable$op == ":="),
                 x.ceq.idx=which(lavpartable$op == "=="),
                 x.cin.idx=which(lavpartable$op == ">" | lavpartable$op == "<"),

                 eq.constraints      = CON$ceq.linear.only.flag,
                 eq.constraints.K    = CON$ceq.JAC.NULL,
                 eq.constraints.k0   = CON$ceq.rhs.NULL,

                 def.function        = CON$def.function,
                 ceq.function        = CON$ceq.function,
                 ceq.JAC             = CON$ceq.JAC,
                 ceq.rhs             = CON$ceq.rhs,
                 ceq.jacobian        = CON$ceq.jacobian,
                 ceq.linear.idx      = CON$ceq.linear.idx,
                 ceq.nonlinear.idx   = CON$ceq.nonlinear.idx,

                 cin.function        = CON$cin.function,
                 cin.JAC             = CON$cin.JAC,
                 cin.rhs             = CON$cin.rhs,
                 cin.jacobian        = CON$cin.jacobian,
                 cin.linear.idx      = CON$cin.linear.idx,
                 cin.nonlinear.idx   = CON$cin.nonlinear.idx,

                 con.jac             = con.jac,
                 con.lambda          = con.lambda,

                 nexo                = nexo,
                 fixed.x             = lavoptions$fixed.x,
                 conditional.x       = lavoptions$conditional.x,
                 parameterization    = lavoptions$parameterization,

                 ov.x.dummy.ov.idx   = ov.x.dummy.ov.idx,
                 ov.x.dummy.lv.idx   = ov.x.dummy.lv.idx,
                 ov.y.dummy.ov.idx   = ov.y.dummy.ov.idx,
                 ov.y.dummy.lv.idx   = ov.y.dummy.lv.idx,

                 estimator           = lavoptions$estimator)

    if(lavoptions$debug) {
         cat("lavaan lavoptions$debug: lavaanModel\n")
         print( str(Model) )
         print( Model@GLIST )
    }

    Model
}

# for backwards compatibility
# Model <- lav_model
