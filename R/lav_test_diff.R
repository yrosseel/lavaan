# various ways to compute a (scaled) difference chi-square test statistic

lav_test_diff_Satorra2000 <- function(m1, m0, H1 = TRUE, A.method = "delta",
                                      A = NULL,
                                      Satterthwaite = FALSE,
                                      debug = FALSE) {

    # extract information from m1 and m2
    T1 <- m1@Fit@test[[1]]$stat
    r1 <- m1@Fit@test[[1]]$df

    T0 <- m0@Fit@test[[1]]$stat
    r0 <- m0@Fit@test[[1]]$df

    # m = difference between the df's
    m <- r0 - r1
     
    Gamma <- lavTech(m1, "Gamma") # the same for m1 and m0
    if(H1) {
        WLS.V <- lavTech(m1, "WLS.V")
        PI <- computeDelta(m1@Model)
        P <- lavTech(m1, "information")
        # needed? (yes, if H1 already has eq constraints)
        P.inv <- lav_model_information_augment_invert(m1@Model,
                                                      information = P,
                                                      inverted = TRUE)
        #P.inv <- solve(P)
    
        # compute 'A' matrix 
        # NOTE: order of parameters may change between H1 and H0, so be
        # careful!
        if(is.null(A)) {
            A <- lav_test_diff_A(m1, m0, method = A.method)
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
        #P.inv <- solve(P)

        # compute 'A' matrix 
        # NOTE: order of parameters may change between H1 and H0, so be
        # careful!
        if(is.null(A)) {
            # m1, m0 OR m0, m1 (works for delta, but not for exact)
            A <- lav_test_diff_A(m1, m0, method = A.method)
            if(debug) print(A)
        }
    }

    # compute tr UG per group
    ngroups <- m1@SampleStats@ngroups
    UG.group  <- vector("list", length=ngroups)

    # PAAPAAP
    PAAPAAP <- P.inv %*% t(A) %*% solve(A %*% P.inv %*% t(A)) %*% A %*% P.inv

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
    fg <- unlist(m1@SampleStats@nobs[[g]])/m1@SampleStats@ntotal

    if(Satterthwaite) {
        cd <- sum(fg * trace.UGamma2) / sum(fg * trace.UGamma)
        df.delta <- (sum(fg * trace.UGamma))^2 / sum(fg * trace.UGamma2)
    } else {
        cd <- 1/m * sum(fg * trace.UGamma)
        df.delta <- m
    }

    # compute scaled difference test      
    T.delta <- (T0 - T1)/cd

    list(T.delta = T.delta, scaling.factor = cd, df.delta = df.delta)
}

lav_test_diff_SatorraBentler2001 <- function(m1, m0) {
    
    # extract information from m1 and m2
    T1 <- m1@Fit@test[[1]]$stat
    r1 <- m1@Fit@test[[1]]$df
    c1 <- m1@Fit@test[[2]]$scaling.factor
    if(r1 == 0) { # saturated model
        c1 <- 1
    }

    T0 <- m0@Fit@test[[1]]$stat
    r0 <- m0@Fit@test[[1]]$df
    c0 <- m0@Fit@test[[2]]$scaling.factor

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

    # extract information from m1 and m2
    T1 <- m1@Fit@test[[1]]$stat
    r1 <- m1@Fit@test[[1]]$df
    c1 <- m1@Fit@test[[2]]$scaling.factor
    if(r1 == 0) { # saturated model
        c1 <- 1
    }

    T0 <- m0@Fit@test[[1]]$stat
    r0 <- m0@Fit@test[[1]]$df
    c0 <- m0@Fit@test[[2]]$scaling.factor
    if(r0 == 0) { # should never happen
        c0 <- 1
    }

    # m = difference between the df's
    m = r0 - r1

    # generate `M10' model
    if(H1) {
        # M0 with M1 parameters
        M01 <- lav_test_diff_m10(m0, m1, test = TRUE)
        c01 <- M01@Fit@test[[2]]$scaling.factor

        # compute c_d
        # cd.01 <- (r0 * c01 - r1 * c0) / m ???
        cd <- (r0 * c0 - r1 * c01) / m
    } else {
        # M1 with M0 parameters (as in Satorra & Bentler 2010)
        M10 <- lav_test_diff_m10(m1, m0, test = TRUE)
        c10 <- M10@Fit@test[[2]]$scaling.factor

        # compute c_d
        cd <- (r0 * c0 - r1 * c10) / m
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
    PT.M1.FULL <- lav_partable_full(PT.M1, free = TRUE, start = TRUE)
    PT.M1.extended <- lav_partable_merge(PT.M1, PT.M1.FULL,
                                         remove.duplicated = TRUE, warn = FALSE)

    # `extend' PT.M0 partable to include all `fixed-to-zero parameters'
    PT.M0.FULL <- lav_partable_full(PT.M0, free = TRUE, start = TRUE)
    PT.M0.extended <- lav_partable_merge(PT.M0, PT.M0.FULL,
                                         remove.duplicated = TRUE, warn = FALSE)

    # `extend' PE of M0 to include all `fixed-to-zero parameters'
    PE.M0 <- parameterEstimates(m0, remove.eq = FALSE, remove.ineq = FALSE)
    PE.M0.FULL <- lav_partable_full(PE.M0)
    PE.M0.extended <- lav_partable_merge(PE.M0, PE.M0.FULL,
                                         remove.duplicated = TRUE, warn = FALSE)

    # FIXME:
    # - check if H0 does not contain additional parameters...

    m10 <- lavaan(model = PT.M1.extended,
                  start = PE.M0.extended,
                  control=list(optim.method          = "none",
                               optim.force.converged = TRUE) ,
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
lav_test_diff_A <- function(m1, m0, method = "exact") {

    # FIXME!!!!

    if(method == "exact") {
        # not ready yet, just for testing
        af <- .test_compute_partable_A_diff(m1 = m1, m0 = m0)
        A <- lav_func_jacobian_complex(func = af, x = m1@Fit@x)
    } else if(method == "delta") {
        # use a numeric approximation of `A'
        Delta1 <- do.call(rbind, computeDelta(m1@Model))
        Delta0 <- do.call(rbind, computeDelta(m0@Model))

        # take into account equality constraints m0
        if(m0@Model@eq.constraints) {
            # the normalization creates a lot of distortion...
            Delta0 <- Delta0 %*% m0@Model@eq.constraints.K
        }

        # take into account equality constraints m1
        if(m1@Model@eq.constraints) {
            # we need a better solution here...
            warning("lavaan WARNING: H1 contains equality constraints; this routine can not handle this (yet)")
        }

        # take into account equality constraints m1
        #tDelta1Delta1 <- crossprod(Delta1)
        #tDelta1Delta1.inv <- 
        #    lav_model_information_augment_invert(m1@Model,
        #                                         information = tDelta1Delta1,
        #                                         inverted = TRUE)
        #H <- solve(t(Delta1) %*% Delta1) %*% t(Delta1) %*% Delta0
        #H <- tDelta1Delta1.inv %*% t(Delta1) %*% Delta0 ## still wrong?
        #                                                ## Delta1 not corrected

        H <- solve(t(Delta1) %*% Delta1) %*% t(Delta1) %*% Delta0
        A <- t(lav_matrix_orthogonal_complement(H))
    }

    A
}

### 
### af()
###
### - problem: we need the plabels, but plabels in h1 need not be
###             the same plabels in h0
###             so we need to somehow 'match' the plabels as they appear
###             in == constraints, etc...
###
.test_compute_partable_A_diff  <- function(m1, m0) {

    PT.M0 <- m0@ParTable
    PT.M1 <- m1@ParTable

    # `extend' PT.M1 partable to include all `fixed-to-zero parameters'
    PT.M1.FULL <- lav_partable_full(PT.M1, free = TRUE, start = TRUE)
    PT.M1.extended <- lav_partable_merge(PT.M1, PT.M1.FULL,
                                         remove.duplicated = TRUE, warn = FALSE)

    # `extend' PT.M0 partable to include all `fixed-to-zero parameters'
    PT.M0.FULL <- lav_partable_full(PT.M0, free = TRUE, start = TRUE)
    PT.M0.extended <- lav_partable_merge(PT.M0, PT.M0.FULL,
                                         remove.duplicated = TRUE, warn = FALSE)

    # for each parameter in P0, see if we have constrained this parameter
    # somehow
    p1 <- PT.M1.extended; np1 <- length(p1$lhs)
    p0 <- PT.M0.extended; np0 <- length(p0$lhs)

    con.function <- function() NULL
    formals(con.function) <- alist(x=, ...=)
    BODY.txt <- paste("{\nout <- numeric(0L)\n", sep="")

    # first come the variable definitions, from H0 only
    DEF.txt <- lav_partable_constraints_def(p0, defTxtOnly=TRUE)
    def.idx <- which(p0$op == ":=")
    BODY.txt <- paste(BODY.txt, DEF.txt, "\n", sep="")

    # then we need all the labels, as they correspond to the parameters
    # FREE only! in H0
    p0.free.idx <- which(p0$free > 0L)
    BODY.txt <- paste(BODY.txt, "# parameter labels\n",
                      paste(  p0$plabel[p0.free.idx],
                              " <- ",
                              paste("x[", p0$free[p0.free.idx], "]",sep=""), 
                              collapse = "\n"),
                      "\n", sep="")

    # we also need any user-specified labels, in case we have
    # explicit == constraints
    user.label.idx <- which(nchar(p0$label) > 0L & p0$free > 0L)
    if(length(user.label.idx) > 0L) {
        BODY.txt <- paste(BODY.txt, "# parameter labels\n",
                      paste(  p0$label[user.label.idx],
                              " <- ",
                              paste("x[", p0$free[user.label.idx], "]",sep=""),
                              collapse = "\n"),
                      "\n", sep="")
    }
   
    # for each parameter in p1, we 'check' is it is fixed to a constant in p0
    ncon <- 0L
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
                    "out[", ncon, "] = x[", p1$free[i], "] - ",
                        p0$ustart[p0.idx], "\n", sep="")
                next
            } else {
                # match, nothing to do
            }
        }
    }

    # add 'new' equality constraints only
    eq.idx <- which(p0$op == "==")
    if(length(eq.idx) > 0L) {
        for(e in eq.idx) {
            # e in p0
            lhs <- p0$lhs[e]; rhs <- p0$rhs[e]

            # do we have the same == in h1? if so, ignore
            p1.idx <- which(p1$lhs == lhs & p1$op == "==" & p1$rhs == rhs)

            if(length(p1.idx) > 0) {
                # ignore
            } else {
                # new == constraint!
                ncon <- ncon + 1L
                eq.string <- paste(p0$lhs[e], " - (", p0$rhs[e], ")", sep="")
                BODY.txt <- paste(BODY.txt, 
                             "out[", ncon, "] <- ", eq.string, "\n", sep="")
            }
        }
    }    

    # wrap function
    BODY.txt <- paste(BODY.txt, "return(out)\n}\n", sep="")
    body(con.function) <- parse(file="", text=BODY.txt)

    con.function
}

