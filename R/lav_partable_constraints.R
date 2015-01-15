# build def function from partable
lav_partable_constraints_def <- function(partable, con = NULL, debug = FALSE,
                                         defTxtOnly = FALSE) {

    # empty function
    def.function <- function() NULL

    # if 'con', merge partable + con
    if(!is.null(con)) {
        partable$lhs <- c(partable$lhs, con$lhs)
        partable$op  <- c(partable$op,  con$op )
        partable$rhs <- c(partable$rhs, con$rhs)
    }

    # get := definitions
    def.idx <- which(partable$op == ":=")
    
    # catch empty def
    if(length(def.idx) == 0L) {
        if(defTxtOnly) {
            return(character(0L))
        } else {
            return(def.function)
        }
    }

    # create function
    formals(def.function) <- alist(x=, ...=)
    if(defTxtOnly) {
        BODY.txt <- ""
    } else {
        BODY.txt <- paste("{\n# parameter definitions\n\n")
    }

    lhs.names <- partable$lhs[def.idx]
    def.labels <- all.vars( parse(file="", text=partable$rhs[def.idx]) )
    # remove the ones in lhs.names
    idx <- which(def.labels %in% lhs.names)
    if(length(idx) > 0L) def.labels <- def.labels[-idx]

    # get corresponding 'x' indices
    def.x.idx  <- partable$free[match(def.labels, partable$label)]
    if(any(is.na(def.x.idx))) {
        stop("lavaan ERROR: unknown label(s) in variable definition(s): ",
         paste(def.labels[which(is.na(def.x.idx))], collapse=" "))
    }
    if(any(def.x.idx == 0)) {
        stop("lavaan ERROR: non-free parameter(s) in variable definition(s): ",
            paste(def.labels[which(def.x.idx == 0)], collapse=" "))
    }
    def.x.lab  <- paste("x[", def.x.idx, "]",sep="")
    # put both the labels the function BODY
    if(length(def.x.idx) > 0L) {
        BODY.txt <- paste(BODY.txt, "# parameter labels\n",
            paste(def.labels, " <- ",def.x.lab, collapse="\n"),
            "\n", sep="")
    }

    # write the definitions literally
    BODY.txt <- paste(BODY.txt, "\n# parameter definitions\n", sep="")
    for(i in 1:length(def.idx)) {
        BODY.txt <- paste(BODY.txt,
            lhs.names[i], " <- ", partable$rhs[def.idx[i]], "\n", sep="")
    }

    if(defTxtOnly) return(BODY.txt)

    # put the results in 'out'
    BODY.txt <- paste(BODY.txt, "\nout <- ",
        paste("c(", paste(lhs.names, collapse=","),")\n", sep=""), sep="")
    # what to do with NA values? -> return +Inf???
    BODY.txt <- paste(BODY.txt, "out[is.na(out)] <- Inf\n", sep="")
    BODY.txt <- paste(BODY.txt, "names(out) <- ",
        paste("c(\"", paste(lhs.names, collapse="\",\""), "\")\n", sep=""),
        sep="")
    BODY.txt <- paste(BODY.txt, "return(out)\n}\n", sep="")

    body(def.function) <- parse(file="", text=BODY.txt)
    if(debug) { cat("def.function = \n"); print(def.function); cat("\n") }

    def.function
}

# build ceq function from partable
#     non-trivial equality constraints (linear or nonlinear)
#     convert to 'ceq(x)' function where 'x' is the (free) parameter vector
#     and ceq(x) returns the evaluated equality constraints
#
#     eg. if b1 + b2 == 2 (and b1 correspond to, say,  x[10] and x[17])
#         ceq <- function(x) {
#             out <- rep(NA, 1)
#             b1 = x[10]; b2 = x[17] 
#             out[1] <- b1 + b2 - 2
#         }
lav_partable_constraints_ceq <- function(partable, con = NULL, debug = FALSE) {

    # empty function
    ceq.function <- function() NULL

    # if 'con', merge partable + con
    if(!is.null(con)) {
        partable$lhs <- c(partable$lhs, con$lhs)
        partable$op  <- c(partable$op,  con$op )
        partable$rhs <- c(partable$rhs, con$rhs)
    }
    
    # get equality constraints
    eq.idx <- which(partable$op == "==")

    # catch empty ceq
    if(length(eq.idx) == 0L) {
        return(ceq.function)
    }

    # create function
    formals(ceq.function) <- alist(x=, ...=)
    BODY.txt <- paste("{\nout <- rep(NA, ", length(eq.idx), ")\n", sep="")

    # first come the variable definitions
    DEF.txt <- lav_partable_constraints_def(partable, defTxtOnly=TRUE)
    def.idx <- which(partable$op == ":=")
    BODY.txt <- paste(BODY.txt, DEF.txt, "\n", sep="")


    # extract labels
    lhs.labels <- all.vars( parse(file="", text=partable$lhs[eq.idx]) )
    rhs.labels <- all.vars( parse(file="", text=partable$rhs[eq.idx]) )
    eq.labels <- unique(c(lhs.labels, rhs.labels))
    # remove def.names from eq.labels
    if(length(def.idx) > 0L) {
        def.names <- as.character(partable$lhs[def.idx])
        d.idx <- which(eq.labels %in% def.names)
        if(length(d.idx) > 0) eq.labels <- eq.labels[-d.idx]
    }
    eq.x.idx <- rep(as.integer(NA), length(eq.labels))
    # get user-labels ids
    ulab.idx <- which(eq.labels %in% partable$label)
    if(length(ulab.idx) > 0L) {
        eq.x.idx[ ulab.idx] <- partable$free[match(eq.labels[ulab.idx], 
                                                   partable$label)]
    }
    # get plabels ids
    plab.idx <- which(eq.labels %in% partable$plabel)
    if(length(plab.idx) > 0L) {
        eq.x.idx[ plab.idx] <- partable$free[match(eq.labels[plab.idx],  
                                                   partable$plabel)]
    }

    # check if we have found the label
    if(any(is.na(eq.x.idx))) {
        stop("lavaan ERROR: unknown label(s) in equality constraint(s): ",
         paste(eq.labels[which(is.na(eq.x.idx))], collapse=" "))
    }
    # check if they are all 'free'
    if(any(eq.x.idx == 0)) {
        fixed.eq.idx <- which(eq.x.idx == 0)
        # FIXME: what should we do here? we used to stop with an error
        # from 0.5.18, we give a warning, and replace the non-free label
        # with its fixed value in ustart
        warning("lavaan WARNING: non-free parameter(s) in equality constraint(s): ",
            paste(eq.labels[fixed.eq.idx], collapse=" "))

        fixed.lab.lhs <- eq.labels[fixed.eq.idx]
        fixed.lab.rhs <- partable$ustart[match(fixed.lab.lhs, partable$label)]
        BODY.txt <- paste(BODY.txt, "# non-free parameter labels\n",
            paste(fixed.lab.lhs, "<-", fixed.lab.rhs, collapse="\n"),
            "\n", sep="")

        eq.x.idx <- eq.x.idx[-fixed.eq.idx]
        eq.labels <- eq.labels[-fixed.eq.idx]
    }

    # put the labels the function BODY
    eq.x.lab  <- paste("x[", eq.x.idx, "]",sep="")
    if(length(eq.x.idx) > 0L) {
        BODY.txt <- paste(BODY.txt, "# parameter labels\n",
            paste(eq.labels, "<-", eq.x.lab, collapse="\n"),
            "\n", sep="")
    }

    # write the equality constraints literally
    BODY.txt <- paste(BODY.txt, "\n# equality constraints\n", sep="")
    for(i in 1:length(eq.idx)) {
        lhs <- partable$lhs[ eq.idx[i] ]
        rhs <- partable$rhs[ eq.idx[i] ]
        if(rhs == "0") { 
            eq.string <- lhs
        } else {
            eq.string <- paste(lhs, " - (", rhs, ")", sep="")
        }
        BODY.txt <- paste(BODY.txt, "out[", i, "] <- ", eq.string, "\n", sep="")
    }
    # put the results in 'out'
    #BODY.txt <- paste(BODY.txt, "\nout <- ",
    #    paste("c(", paste(lhs.names, collapse=","),")\n", sep=""), sep="")

    # what to do with NA values? -> return +Inf???
    BODY.txt <- paste(BODY.txt, "\n", "out[is.na(out)] <- Inf\n", sep="")
    BODY.txt <- paste(BODY.txt, "return(out)\n}\n", sep="")
    body(ceq.function) <- parse(file="", text=BODY.txt)
    if(debug) { cat("ceq.function = \n"); print(ceq.function); cat("\n") }

    ceq.function
}


# build ciq function from partable
#     non-trivial inequality constraints (linear or nonlinear)
#     convert to 'cin(x)' function where 'x' is the (free) parameter vector
#     and cin(x) returns the evaluated inequality constraints
#
#     eg. if b1 + b2 > 2 (and b1 correspond to, say,  x[10] and x[17])
#         cin <- function(x) {
#             out <- rep(NA, 1)
#             b1 = x[10]; b2 = x[17] 
#             out[1] <- b1 + b2 - 2
#         }
#
# NOTE: very similar, but not identitical to ceq, because we need to take
#       care of the difference between '<' and '>'
lav_partable_constraints_ciq <- function(partable, con = NULL, debug = FALSE) {


    # empty function
    cin.function <- function() NULL

    # if 'con', merge partable + con
    if(!is.null(con)) {
        partable$lhs <- c(partable$lhs, con$lhs)
        partable$op  <- c(partable$op,  con$op )
        partable$rhs <- c(partable$rhs, con$rhs)
    }
    
    # get inequality constraints
    ineq.idx <- which(partable$op == ">" | partable$op == "<")

    # catch empty ceq
    if(length(ineq.idx) == 0L) {
        return(cin.function)
    }

    # create function
    formals(cin.function) <- alist(x=, ...=)
    BODY.txt <- paste("{\nout <- rep(NA, ", length(ineq.idx), ")\n", sep="")

    # first come the variable definitions
    DEF.txt <- lav_partable_constraints_def(partable, defTxtOnly=TRUE)
    def.idx <- which(partable$op == ":=")
    BODY.txt <- paste(BODY.txt, DEF.txt, "\n", sep="")

    # extract labels
    lhs.labels <- all.vars( parse(file="", text=partable$lhs[ineq.idx]) )
    rhs.labels <- all.vars( parse(file="", text=partable$rhs[ineq.idx]) )
    ineq.labels <- unique(c(lhs.labels, rhs.labels))
    # remove def.names from ineq.labels
    if(length(def.idx) > 0L) {
        def.names <- as.character(partable$lhs[def.idx])
        d.idx <- which(ineq.labels %in% def.names)
        if(length(d.idx) > 0) ineq.labels <- ineq.labels[-d.idx]
    }
    ineq.x.idx <- rep(as.integer(NA), length(ineq.labels))
    # get user-labels ids
    ulab.idx <- which(ineq.labels %in% partable$label)
    if(length(ulab.idx) > 0L) {
        ineq.x.idx[ ulab.idx] <- partable$free[match(ineq.labels[ulab.idx], 
                                                   partable$label)]
    }
    # get plabels ids
    plab.idx <- which(ineq.labels %in% partable$plabel)
    if(length(plab.idx) > 0L) {
        ineq.x.idx[ plab.idx] <- partable$free[match(ineq.labels[plab.idx],  
                                                   partable$plabel)]
    }

    # check if we have found the label
    if(any(is.na(ineq.x.idx))) {
        stop("lavaan ERROR: unknown label(s) in inequality constraint(s): ",
         paste(ineq.labels[which(is.na(ineq.x.idx))], collapse=" "))
    }
    # check if they are all 'free'
    if(any(ineq.x.idx == 0)) {
        fixed.ineq.idx <- which(ineq.x.idx == 0)
        # FIXME: what should we do here? we used to stop with an error
        # from 0.5.18, we give a warning, and replace the non-free label
        # with its fixed value in ustart
        warning("lavaan WARNING: non-free parameter(s) in inequality constraint(s): ",
            paste(ineq.labels[fixed.ineq.idx], collapse=" "))

        fixed.lab.lhs <- ineq.labels[fixed.ineq.idx]
        fixed.lab.rhs <- partable$ustart[match(fixed.lab.lhs, partable$label)]
        BODY.txt <- paste(BODY.txt, "# non-free parameter labels\n",
            paste(fixed.lab.lhs, "<-", fixed.lab.rhs, collapse="\n"),
            "\n", sep="")

        ineq.x.idx <- ineq.x.idx[-fixed.ineq.idx]
        ineq.labels <- ineq.labels[-fixed.ineq.idx]
    }

    # put the labels the function BODY
    ineq.x.lab  <- paste("x[", ineq.x.idx, "]",sep="")
    if(length(ineq.x.idx) > 0L) {
        BODY.txt <- paste(BODY.txt, "# parameter labels\n",
            paste(ineq.labels, "<-", ineq.x.lab, collapse="\n"),
            "\n", sep="")
    }

    # write the constraints literally
    BODY.txt <- paste(BODY.txt, "\n# inequality constraints\n", sep="")
    for(i in 1:length(ineq.idx)) {
        lhs <- partable$lhs[ ineq.idx[i] ]
         op <- partable$op[  ineq.idx[i] ]
        rhs <- partable$rhs[ ineq.idx[i] ]

        # note,this is different from ==, because we have < AND >
        if(rhs == "0" && op == ">") {
            ineq.string <- lhs
        } else if(rhs == "0" && op == "<") {
            ineq.string <- paste(rhs, " - (", lhs, ")", sep="")
        } else if(rhs != "0" && op == ">") {
            ineq.string <- paste(lhs, " - (", rhs, ")", sep="")
        } else if(rhs != "0" && op == "<") {
            ineq.string <- paste(rhs, " - (", lhs, ")", sep="")
        }

        BODY.txt <- paste(BODY.txt, "out[", i, "] <- ", ineq.string, "\n", sep="")
    }
    # put the results in 'out'
    #BODY.txt <- paste(BODY.txt, "\nout <- ",
    #    paste("c(", paste(lhs.names, collapse=","),")\n", sep=""), sep="")

    # what to do with NA values? -> return +Inf???
    BODY.txt <- paste(BODY.txt, "\n", "out[is.na(out)] <- Inf\n", sep="")
    BODY.txt <- paste(BODY.txt, "return(out)\n}\n", sep="")
    body(cin.function) <- parse(file="", text=BODY.txt)
    if(debug) { cat("cin.function = \n"); print(cin.function); cat("\n") }

    cin.function
}


# given two models, M1 and M0, where M0 is nested in M1,
# create a function 'af(x)' where 'x' is the full parameter vector of M1
# and af(x) returns the evaluated restrictions under M0).
# The 'jacobian' of this function 'A' will be used in the anova
# anova() function, and elsewhere
lav_partable_constraints_function <- function(p1, p0) {

    stop("lavaan DEBUG: this function needs revision; use SB.classic = TRUE for now")

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
                    stop("lavaan ERROR: non-free parameter(s) in inequality constraint(s): ",
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

