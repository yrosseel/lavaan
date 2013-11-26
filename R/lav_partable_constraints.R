
# given two models, M1 and M0, where M0 is nested in M1,
# create a function 'af(x)' where 'x' is the full parameter vector of M1
# and af(x) returns the evaluated restrictions under M0).
# The 'jacobian' of this function 'A' will be used in the anova
# anova() function, and elsewhere
lav_partable_constraints_function <- function(p1, p0) {

    # check for inequality constraints
    if(any(c(p0$op,p1$op) %in% c(">","<"))) 
        stop("lavaan ERROR: anova() can not handle inequality constraints; use InformativeTesting() instead")

    npar.p1 <- max(p1$free)
    npar.p0 <- max(p0$free)
    npar.diff <- npar.p1 - npar.p0

    con.function <- function() NULL
    formals(con.function) <- alist(x=, ...=)
    BODY.txt <- paste("{\nout <- rep(NA, ", npar.diff, ")\n", sep="")

    # for each free parameter in p1, we 'check' is it is somehow 
    # restricted in p0
    ncon <- 0L; EQ.ID <- integer(); EQ.P1 <- integer()
    for(i in seq_len(npar.p1)) {
        idx <- which(p1$free == i)[1L]
        lhs <- p1$lhs[idx]; op <- p1$op[idx]; rhs <- p1$rhs[idx]
        group <- p1$group[idx]
        p0.idx <- which(p0$lhs == lhs & p0$op == op & p0$rhs == rhs &
                        p0$group == group)
        if(length(p0.idx) == 0L) {
            # this parameter is not to be found in M0, we will
            # assume that is fixed to zero
            ncon <- ncon + 1L 
            BODY.txt <- paste(BODY.txt,
                "out[", ncon, "] = x[", i, "] - 0\n", sep="")
            next
        }
        if(p0$free[p0.idx] == 0L) {
            ncon <- ncon + 1L
            # simple fixed value
            BODY.txt <- paste(BODY.txt,
                "out[", ncon, "] = x[", i, "] - ",
                                   p0$ustart[p0.idx], "\n", sep="")
        } 
        # how to deal with *new* equality constraints?
        # check if p0 has an eq.id not yet in EQ.ID (while p1 eq.id is empty)
        if(p0$eq.id[p0.idx] != 0 && p1$eq.id[idx] == 0) {
            # if not in EQ.ID, put it there and continue
            EQ <- p0$eq.id[p0.idx]
            if(EQ %in% EQ.ID) {
                # add constraint
                ncon <- ncon + 1L
                BODY.txt <- paste(BODY.txt,
                                  "out[", ncon, "] = x[", EQ.P1[EQ.ID == EQ], 
                                                "] - x[", i, "]\n", sep="")
            } else {
                EQ.ID <- c(EQ.ID, EQ)
                EQ.P1 <- c(EQ.P1,  i)
            }
        }
    }

    # add all NEW constraints using "==" (not already in p1)
    p0.con.idx <- which(p0$op == "==")
    if(length(p0.con.idx) > 0L) {
        # first, remove those also present in p1
        del.idx <- integer(0L)
        for(con in 1:length(p0.con.idx)) {
            p0.idx <- p0.con.idx[con]
            p1.idx <- which(p1$op == "==" &
                            p1$lhs == p0$lhs[p0.idx] &
                            p1$rhs == p0$rhs[p0.idx])
            if(length(p1.idx) > 0L)
                del.idx <- c(del.idx, con)
        }
        if(length(del.idx) > 0L)
            p0.con.idx <- p0.con.idx[-del.idx]   
    }
    eq.idx <- p0.con.idx
    if(length(eq.idx) > 0L) {
        def.idx <- which(p0$op == ":=")
        # first come the variable definitions
        if(length(def.idx) > 0L) {
            for(i in 1:length(def.idx)) {
                lhs <- p0$lhs[ def.idx[i] ]
                rhs <- p0$rhs[ def.idx[i] ]
                def.string <- rhs
                # coerce to expression to extract variable names
                def.labels <- all.vars( parse(file="", text=def.string) )
                # get corresponding 'x' indices
                def.x.idx  <- p0$free[match(def.labels, p0$label)]
                def.x.lab  <- paste("x[", def.x.idx, "]",sep="")
                # put both the labels and the expression in the function BODY
                BODY.txt <- paste(BODY.txt,
                    paste(def.labels, "=",def.x.lab, collapse=";"),"\n",
                    lhs, " = ", def.string, "\n", sep="")
            }
        }

        for(i in 1:length(eq.idx)) {
            lhs <- p0$lhs[ eq.idx[i] ]
            rhs <- p0$rhs[ eq.idx[i] ]
            if(rhs == "0") {
                eq.string <- lhs
            } else {
                eq.string <- paste(lhs, "- (", rhs, ")", sep="")
            }
            # coerce to expression to extract variable names
            eq.labels <- all.vars( parse(file="", text=eq.string) )
            # get corresponding 'x' indices
            if(length(def.idx) > 0L) {
                # remove def.names from ineq.labels
                def.names <- as.character(p0$lhs[def.idx])
                d.idx <- which(eq.labels %in% def.names)
                if(length(d.idx) > 0) eq.labels <- eq.labels[-d.idx]
            }
            if(length(eq.labels) > 0L) {
                eq.x.idx  <- p0$free[match(eq.labels, p0$label)]
                if(any(is.na(eq.x.idx))) {
                    stop("lavaan ERROR: unknown label(s) in equality constraint: ",
                         paste(eq.labels[which(is.na(eq.x.idx))], collapse=" "))
                }
                if(any(eq.x.idx == 0)) {
                    stop("lavaan ERROR: non-free parameter(s) in inequality constraint: ",
                        paste(eq.labels[which(eq.x.idx == 0)], collapse=" "))
                }
                eq.x.lab  <- paste("x[", eq.x.idx, "]",sep="")
                # put both the labels and the expression in the function BODY
                BODY.txt <- paste(BODY.txt,
                    paste(eq.labels, "=", eq.x.lab, collapse=";"),"\n",
                    "out[", i, "] = ", eq.string, "\n", sep="")
            } else {
                BODY.txt <- paste(BODY.txt,
                    "out[", i, "] = ", eq.string, "\n", sep="")
            }
        }
               
    }

    # wrap function
    BODY.txt <- paste(BODY.txt, "return(out)\n}\n", sep="")
    body(con.function) <- parse(file="", text=BODY.txt)

    con.function
}

