# various ways to compute a (scaled) difference chi-square test statistic

lav_test_diff_Satorra2000 <- function(m1, m0, H1 = TRUE, A.method = "delta",
                                      A = NULL,
                                      Satterthwaite = FALSE,
                                      scaled.shifted = FALSE,
                                      debug = FALSE) {

    if(scaled.shifted) {
        Satterthwaite <- TRUE
    }

    # extract information from m1 and m2
    T1 <- m1@test[[1]]$stat
    r1 <- m1@test[[1]]$df

    T0 <- m0@test[[1]]$stat
    r0 <- m0@test[[1]]$df

    # m = difference between the df's
    m <- r0 - r1

    Gamma <- lavTech(m1, "Gamma") # the same for m1 and m0
    # check for NULL
    if(is.null(Gamma)) {
        stop("lavaan ERROR: can not compute Gamma matrix; perhaps missing = \"ml\"?")
    }

    if(H1) {
        WLS.V <- lavTech(m1, "WLS.V")
        PI <- computeDelta(m1@Model)
        P <- lavTech(m1, "information")
        # needed? (yes, if H1 already has eq constraints)
        P.inv <- lav_model_information_augment_invert(m1@Model,
                                                      information = P,
                                                      inverted = TRUE)
        # compute 'A' matrix
        # NOTE: order of parameters may change between H1 and H0, so be
        # careful!
        if(is.null(A)) {
            A <- lav_test_diff_A(m1, m0, method = A.method, reference = "H1")
            # take into account equality constraints m1
            if(m1@Model@eq.constraints) {
                A <- A %*% t(m1@Model@eq.constraints.K)
            }
            if(debug) print(A)
        }
    } else {
        stop("not ready yet")

        WLS.V <- lavTech(m0, "WLS.V")
        PI <- computeDelta(m0@Model)
        P <- lavTech(m0, "information")
        # needed?
        P.inv <- lav_model_information_augment_invert(m0@Model,
                                                      information = P,
                                                      inverted = TRUE)

        # compute 'A' matrix
        # NOTE: order of parameters may change between H1 and H0, so be
        # careful!
        if(is.null(A)) {
            # m1, m0 OR m0, m1 (works for delta, but not for exact)
            A <- lav_test_diff_A(m1, m0, method = A.method, reference = "H0")
            # take into account equality constraints m1
            if(m0@Model@eq.constraints) {
                A <- A %*% t(m0@Model@eq.constraints.K)
            }
            if(debug) print(A)
        }
    }

    # compute tr UG per group
    ngroups <- m1@SampleStats@ngroups
    UG.group  <- vector("list", length=ngroups)

    # safety check: A %*% P.inv %*% t(A) should NOT contain all-zero
    # rows/columns
    # FIXME: is this really needed? As we use ginv later on
    APA <- A %*% P.inv %*% t(A)
    cSums <- colSums(APA)
    rSums <- rowSums(APA)
    empty.idx <- which( abs(cSums) < .Machine$double.eps^0.5 &
                        abs(rSums) < .Machine$double.eps^0.5 )
    if(length(empty.idx) > 0) {
        A <- A[-empty.idx,, drop = FALSE]
    }

    # PAAPAAP
    PAAPAAP <- P.inv %*% t(A) %*% MASS::ginv(A %*% P.inv %*% t(A)) %*% A %*% P.inv

    trace.UGamma  <- numeric(ngroups)
    trace.UGamma2 <- numeric(ngroups)
    for(g in 1:ngroups) {
        UG.group <- WLS.V[[g]] %*% Gamma[[g]] %*% WLS.V[[g]] %*%
                    PI[[g]] %*% PAAPAAP %*% t(PI[[g]])
        trace.UGamma[g]  <- sum(diag(UG.group))
        if(Satterthwaite) {
            trace.UGamma2[g] <- sum(diag(UG.group %*% UG.group))
        }
    }

    # compute scaling factor
    fg <- unlist(m1@SampleStats@nobs)/m1@SampleStats@ntotal

    trace.UGamma <- sum(fg * trace.UGamma)
    if(Satterthwaite) {
        trace.UGamma2 <- sum(fg * trace.UGamma2)
    }

    if(Satterthwaite && !scaled.shifted) {
        cd <- trace.UGamma2 / trace.UGamma
        df.delta <- trace.UGamma^2 / trace.UGamma2
        T.delta <- (T0 - T1)/cd
        a <- as.numeric(NA); b <- as.numeric(NA)
    } else if(Satterthwaite && scaled.shifted) {
        a <- sqrt(m/trace.UGamma2)
        #b <- m - sqrt(m * trace.UGamma^2 / trace.UGamma2)
        b <- m - a * trace.UGamma
        df.delta <- m
        T.delta <- (T0 - T1)*a + b
        cd <- as.numeric(NA)
    } else {
        cd <- 1/m * trace.UGamma
        df.delta <- m
        T.delta <- (T0 - T1)/cd
        a <- as.numeric(NA); b <- as.numeric(NA)
    }

    list(T.delta = T.delta, scaling.factor = cd, df.delta = df.delta,
         a = a, b = b)
}

lav_test_diff_SatorraBentler2001 <- function(m1, m0) {

    # extract information from m1 and m2
    T1 <- m1@test[[1]]$stat
    r1 <- m1@test[[1]]$df
    c1 <- m1@test[[2]]$scaling.factor
    if(r1 == 0) { # saturated model
        c1 <- 1
    }

    T0 <- m0@test[[1]]$stat
    r0 <- m0@test[[1]]$df
    c0 <- m0@test[[2]]$scaling.factor

    # m = difference between the df's
    m = r0 - r1

    # compute c_d
    cd <- (r0 * c0 - r1 * c1) / m

    # warn if cd is negative
    if(cd < 0) {
        warning("lavaan WARNING: scaling factor is negative")
        cd <- as.numeric(NA)
    }

    # compute scaled difference test
    T.delta <- (T0 - T1)/cd

    list(T.delta = T.delta, scaling.factor = cd, df.delta = m)
}

lav_test_diff_SatorraBentler2010 <- function(m1, m0, H1 = FALSE) {

    ### FIXME: check if models are nested at the parameter level!!!

    # extract information from m1 and m2
    T1 <- m1@test[[1]]$stat
    r1 <- m1@test[[1]]$df
    c1 <- m1@test[[2]]$scaling.factor
    if(r1 == 0) { # saturated model
        c1 <- 1
    }

    T0 <- m0@test[[1]]$stat
    r0 <- m0@test[[1]]$df
    c0 <- m0@test[[2]]$scaling.factor
    if(r0 == 0) { # should never happen
        c0 <- 1
    }

    # m = difference between the df's
    m = r0 - r1

    # generate `M10' model
    if(H1) {
        # M0 with M1 parameters
        M01 <- lav_test_diff_m10(m0, m1, test = TRUE)
        c01 <- M01@test[[2]]$scaling.factor

        # check if vcov is positive definite (new in 0.6)
        # if not, we may get negative values
        eigvals <- eigen(lavTech(M10, "information"),
                         symmetric=TRUE, only.values=TRUE)$values
        if(any(eigvals < -1 * .Machine$double.eps^(3/4))) {
            warning(
  "lavaan WARNING: information matrix of the M10 model is not positive definite.\n",
"                  As a result, the scale-factor can not be computed.")
            cd <- as.numeric(NA)
        } else {
            # compute c_d
            # cd.01 <- (r0 * c01 - r1 * c0) / m ???
            cd <- (r0 * c0 - r1 * c01) / m
        }

    } else {
        # M1 with M0 parameters (as in Satorra & Bentler 2010)
        M10 <- lav_test_diff_m10(m1, m0, test = TRUE)
        c10 <- M10@test[[2]]$scaling.factor

        # check if vcov is positive definite (new in 0.6)
        # if not, we may get negative values
        eigvals <- eigen(lavTech(M10, "information"),
                         symmetric=TRUE, only.values=TRUE)$values
        if(any(eigvals < -1 * .Machine$double.eps^(3/4))) {
            warning(
  "lavaan WARNING: information matrix of the M10 model is not positive definite.\n",
"                  As a result, the scale-factor can not be computed.")
            cd <- as.numeric(NA)
        } else {
            # compute c_d
            cd <- (r0 * c0 - r1 * c10) / m
        }
    }

    # compute scaled difference test
    T.delta <- (T0 - T1)/cd

    list(T.delta = T.delta, scaling.factor = cd, df.delta = m,
         T.delta.unscaled = (T0 - T1))
}

# create a new model 'm10', where we use model 'm1', but we
# inject it with the values of 'm0'
lav_test_diff_m10 <- function(m1, m0, test = FALSE) {

    # switch of verbose/se/test
    Options <- m1@Options
    Options$verbose <- FALSE

    # should we compute se/test statistics?
    if(!test) {
        Options$se <- "none"; Options$test <- "none"
    }

    PT.M0 <- m0@ParTable
    PT.M1 <- m1@ParTable

    # `extend' PT.M1 partable to include all `fixed-to-zero parameters'
    PT.M1.FULL <- lav_partable_full(partable = PT.M1, lavpta = m1@pta,
                                    free = TRUE, start = TRUE)
    PT.M1.extended <- lav_partable_merge(PT.M1, PT.M1.FULL,
                                         remove.duplicated = TRUE, warn = FALSE)

    # `extend' PT.M0 partable to include all `fixed-to-zero parameters'
    PT.M0.FULL <- lav_partable_full(partable = PT.M0, lavpta = m0@pta,
                                    free = TRUE, start = TRUE)
    PT.M0.extended <- lav_partable_merge(PT.M0, PT.M0.FULL,
                                         remove.duplicated = TRUE, warn = FALSE)

    # `extend' PE of M0 to include all `fixed-to-zero parameters'
    PE.M0 <- parameterEstimates(m0, remove.eq = FALSE, remove.ineq = FALSE,
                                remove.system.eq =  FALSE, remove.def = FALSE)
    PE.M0.FULL <- lav_partable_full(PE.M0)
    PE.M0.extended <- lav_partable_merge(PE.M0, PE.M0.FULL,
                                         remove.duplicated = TRUE, warn = FALSE)

    # FIXME:
    # - check if H0 does not contain additional parameters...

    Options$optim.method          = "none"
    Options$optim.force.converged = TRUE
    Options$start                 = PE.M0.extended # new in 0.6!
    m10 <- lavaan(model = PT.M1.extended,
                  #start = PE.M0.extended,
                  slotOptions     = Options,
                  slotSampleStats = m1@SampleStats,
                  slotData        = m1@Data,
                  slotCache       = m1@Cache)

    m10
}

# compute the `A' matrix: the jacobian of the constraint function a(\delta)
# (see Satorra 2000)
#
#
#
lav_test_diff_A <- function(m1, m0, method = "delta", reference = "H1") {

    # FIXME!!!!

    if(method == "exact") {
        if(reference == "H1") {
            af <- lav_test_diff_af_h1(m1 = m1, m0 = m0)
            xx <- m1@optim$x
        } else { # evaluate under H0
            stop("not ready yet")
            #af <- .test_compute_partable_A_diff_h0(m1 = m1, m0 = m0)
            xx <- m0@optim$x
        }
        A <- try(lav_func_jacobian_complex(func = af, x = xx), silent = TRUE)
        if(inherits(A, "try-error")) {
            A <- lav_func_jacobian_simple(func = af, x = xx)
        }
    } else if(method == "delta") {
        # use a numeric approximation of `A'
        Delta1.list <- computeDelta(m1@Model)
        Delta0.list <- computeDelta(m0@Model)
        Delta1 <- do.call(rbind, Delta1.list)
        Delta0 <- do.call(rbind, Delta0.list)

        # take into account equality constraints m0
        if(m0@Model@eq.constraints) {
            Delta0 <- Delta0 %*% m0@Model@eq.constraints.K
        }

        # take into account equality constraints m1
        if(m1@Model@eq.constraints) {
            Delta1 <- Delta1 %*% m1@Model@eq.constraints.K
        }

        #H <- solve(t(Delta1) %*% Delta1) %*% t(Delta1) %*% Delta0
        H <- MASS::ginv(Delta1) %*% Delta0
        A <- t(lav_matrix_orthogonal_complement(H))
    }

    A
}


# for each parameter in H1 (m1), see if we have somehow constrained
# this parameter under H0 (m0)
#
# since we work 'under H0', we need to use the labels/constraints/def
# as they appear in H0. Unfortunately, the order of the parameters, and
# even the (p)labels may be different in the two models...
#
# Therefore, we will attempt to:
#   - change the 'order' of the 'free' column in m0, so that they map to
#     to the 'x' that we will provide from H1
#   - the plabels used in "==" constraints must be renamed, if necessary
#
lav_test_diff_af_h1 <- function(m1, m0) {

    PT.M0 <- parTable(m0)
    PT.M1 <- parTable(m1)

    # select .p*. parameters only
    M0.p.idx <- which(grepl("\\.p", PT.M0$plabel)); np0 <- length(M0.p.idx)
    M1.p.idx <- which(grepl("\\.p", PT.M1$plabel)); np1 <- length(M1.p.idx)

    # check if parameter space is the same
    if(np0 != np1) {
        stop("lavaan ERROR: unconstrained parameter set is not the same in m0 and m1")
    }

    # split partable in 'parameter' and 'constraints' section
    PT.M0.part1 <- PT.M0[ M0.p.idx,]
    PT.M0.part2 <- PT.M0[-M0.p.idx,]

    PT.M1.part1 <- PT.M1[ M1.p.idx,]
    PT.M1.part2 <- PT.M1[-M1.p.idx,]

    #figure out relationship between m0 and m1
    p1.id <- lav_partable_map_id_p1_in_p2(PT.M0.part1, PT.M1.part1)
    p0.free.idx <- which(PT.M0.part1$free > 0)

    # change 'free' order in m0
    # NOTE: this only works all the free parameters in h0 are also free
    # in h1 (and if not, they will become fixed in h0)
    PT.M0.part1$free[p0.free.idx] <-
    PT.M1.part1$free[ PT.M0.part1$id[p1.id][p0.free.idx] ]

    # paste back
    PT.M0 <- rbind(PT.M0.part1, PT.M0.part2)
    PT.M1 <- rbind(PT.M1.part1, PT.M1.part2)

    # `extend' PT.M1 partable to include all `fixed-to-zero parameters'
    PT.M1.FULL <- lav_partable_full(partable = PT.M1, lavpta = m1@pta,
                                    free = TRUE, start = TRUE)
    PT.M1.extended <- lav_partable_merge(PT.M1, PT.M1.FULL,
                                         remove.duplicated = TRUE, warn = FALSE)

    # `extend' PT.M0 partable to include all `fixed-to-zero parameters'
    PT.M0.FULL <- lav_partable_full(partable = PT.M0, lavpta = m0@pta,
                                    free = TRUE, start = TRUE)
    PT.M0.extended <- lav_partable_merge(PT.M0, PT.M0.FULL,
                                         remove.duplicated = TRUE, warn = FALSE)

    p1 <- PT.M1.extended; np1 <- length(p1$lhs)
    p0 <- PT.M0.extended; np0 <- length(p0$lhs)

    con.function <- function() NULL
    formals(con.function) <- alist(.x.=, ...=)
    BODY.txt <- paste("{\nout <- numeric(0L)\n", sep = "")


    # first handle def + == constraints
    # but FIRST, remove == constraints that also appear in H1!!!

    # remove equivalent eq constraints from p0
    P0 <- p0

    p0.eq.idx <- which(p0$op == "==")
    p1.eq.idx <- which(p1$op == "==")
    p0.remove.idx <- integer(0L)
    if(length(p0.eq.idx) > 0L) {
        for(i in seq_along(p0.eq.idx)) {
            # e0 in p0
            e0 <- p0.eq.idx[i]
            lhs <- p0$lhs[e0]; rhs <- p0$rhs[e0]

            # do we have an equivalent constraint in H1?
            # NOTE!! the (p)labels may differ

            # SO, we will use an 'empirical' approach: if we fill in (random)
            # values, and work out the constraint, do we get identical values?
            # if yes, constraint is equivalent, and we should NOT add it here

            if(length(p1.eq.idx) > 0) {
                # generate random parameter values
                xx1 <- rnorm( length(M1.p.idx) )
                xx0 <- xx1[ p1.id ]

                con.h0.value <- m0@Model@ceq.function(xx0)[i]
                con.h1.values <- m1@Model@ceq.function(xx1)

                if(con.h0.value %in% con.h1.values) {
                    p0.remove.idx <- c(p0.remove.idx, e0)
                }
            }
        }
    }
    if(length(p0.remove.idx) > 0L) {
        P0 <- P0[-p0.remove.idx,]
    }

    # only for the UNIQUE equality constraints in H0, generate syntax
    DEFCON.txt <- lav_partable_constraints_ceq(P0, txtOnly=TRUE)
    BODY.txt <- paste(BODY.txt, DEFCON.txt, "\n", sep="")


    # for each parameter in p1, we 'check' is it is fixed to a constant in p0
    ncon <- length( which(P0$op == "==") )
    for(i in seq_len(np1)) {
        # p in p1
        lhs <- p1$lhs[i]; op <- p1$op[i]; rhs <- p1$rhs[i]; group <- p1$group[i]

        # ignore '==', '<', '>' and ':=' for now
        if(op == "==" || op == ">" || op == "<" || op == ":=") next

        # search for corresponding parameter in p0
        p0.idx <- which(p0$lhs == lhs & p0$op == op & p0$rhs == rhs &
                        p0$group == group)
        if(length(p0.idx) == 0L) {
            stop("lavaan ERROR: parameter in H1 not found in H0: ",
                 paste(lhs, op, rhs, "(group = ", group, ")", sep=" "))
        }

        # 4 possibilities: p is free/fixed in p1, p is free/fixed in p0
        if(p1$free[i] == 0L) {
            if(p0$free[p0.idx] == 0L) {
                # match, nothing to do
            } else {
                warning("lavaan WARNING: fixed parameter in H1 is free in H0: ",
                     paste("\"", lhs, " ", op, " ", rhs,
                           "\" (group = ", group, ")", sep=""))
            }
        } else {
            if(p0$free[p0.idx] == 0L) {
                # match, this is a contrained parameter in H0
                ncon <- ncon + 1L
                BODY.txt <- paste(BODY.txt,
                    "out[", ncon, "] = .x.[", p1$free[i], "] - ",
                        p0$ustart[p0.idx], "\n", sep="")
                next
            } else {
                # match, nothing to do
            }
        }
    }


    # wrap function
    BODY.txt <- paste(BODY.txt, "return(out)\n}\n", sep="")
    body(con.function) <- parse(file="", text=BODY.txt)

    con.function
}

