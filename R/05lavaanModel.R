# constructor of the matrix representation
#
# initial version: YR 22/11/2010

# construct MATRIX representation of the model
Model <- function(partable       = NULL,
                  start          = NULL, 
                  representation = "LISREL",
                  th.idx         = list(),
                  debug          = FALSE) {

    # global info from user model
    ngroups <- max(partable$group)
    meanstructure <- any(partable$op == "~1")
    categorical <- any(partable$op == "|")
    if(categorical) meanstructure <- TRUE


    # what if no starting values are provided? 
    if(is.null(start)) 
        startValues <- StartingValues(start.method="simple", partable=partable)
    else
        startValues <- start
 
    # check start length
    stopifnot(length(startValues) == nrow(partable))

    # only representation = "LISREL" for now
    stopifnot(representation == "LISREL")


    # capture equality constraints in 'K' matrix
    # rows are the unconstrained parameters, cols are the unique parameters
    n.unco <- max(partable$unco)
    n.free       <- max(partable$free)
    if(n.free == n.unco) {
        eq.constraints <- FALSE
        K <- matrix(0, 0, 0)
    } else {
        K <- matrix(0, nrow=n.unco, ncol=n.free)
        #####
        #####     FIXME !
        #####
        idx.free <- partable$free[ partable$free > 0 ]
        for(k in 1:n.unco) {
            c <- idx.free[k]
            K[k, c] <- 1
        }
        eq.constraints <- TRUE
    }

    # Ku matrix (relation th and ov.ord)
    # FIXME (not for mixed!)
    #if(categorical) {
    #    th <- vnames(partable, "th")
    #    ov <- vnames(partable, "ov")
    #    Ku <- t(sapply(ov, grepl, th) + 0L)
    #} else {
    #    Ku <- matrix(0,0,0)
    #}



    # select model matrices
    if(representation == "LISREL") {
        REP <- representation.LISREL(partable, target=NULL, extra=TRUE)
    } else {
        stop("lavaan ERROR: only representation \"LISREL\" has been implemented.")
    }
    if(debug) print(REP)

    # prepare nG-sized slots
    nG <- sum(unlist(attr(REP, "mmNumber")))
    GLIST <- vector(mode="list", nG)
    names(GLIST) <- unlist(attr(REP, "mmNames"))
    dimNames    <- vector(mode="list", length=nG)
    isSymmetric <- logical(nG)
    mmSize      <- integer(nG)

    m.free.idx <- m.unco.idx <- m.user.idx <- 
        vector(mode="list", length=nG)
    x.free.idx <- x.unco.idx <- x.user.idx <-
        vector(mode="list", length=nG)

    # prepare ngroups-sized slots
    nvar <- integer(ngroups)
    nmat <- unlist(attr(REP, "mmNumber"))
    num.idx <- vector("list", length=ngroups)
    nexo <- integer(ngroups)

    offset <- 0L
    for(g in 1:ngroups) {

        # observed and latent variables for this group
        ov.names <- vnames(partable, "ov", group=g)
        ov.names.nox <- vnames(partable, "ov.nox", group=g)
        ov.names.x <- vnames(partable, "ov.x", group=g)
        nexo[g] <- length(ov.names.x)
        ov.num <- vnames(partable, "ov.num", group=g)
        if(categorical) {
            nvar[g] <- length(ov.names.nox)
        } else {
            nvar[g] <- length(ov.names)
        }
        num.idx[[g]] <- match(ov.num, ov.names.nox)

        # model matrices for this group
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
            idx <- which(partable$group == g & REP$mat == mmNames[mm]) 

            # create empty `pattern' matrix
            # FIXME: one day, we may want to use sparse matrices...
            #        but they should not slow things down!
            tmp <- matrix(0L, nrow=mmRows[mm],
                              ncol=mmCols[mm])

            # 1. first assign free values only, to get vector index
            #    -> to be used in computeObjective
            tmp[ cbind(REP$row[idx], REP$col[idx]) ] <- partable$free[idx]
            if(mmSymmetric[mm]) {
                # NOTE: we assume everything is in the UPPER tri!
                T <- t(tmp); tmp[lower.tri(tmp)] <- T[lower.tri(T)]
            }
            m.free.idx[[offset]] <-     which(tmp > 0)
            x.free.idx[[offset]] <- tmp[which(tmp > 0)]

            # 2. if equality constraints, unconstrained free parameters
            #    -> to be used in computeGradient
            if(eq.constraints) {
                tmp[ cbind(REP$row[idx], 
                           REP$col[idx]) ] <- partable$unco[idx]
                if(mmSymmetric[mm]) {
                    # NOTE: we assume everything is in the UPPER tri!
                    T <- t(tmp); tmp[lower.tri(tmp)] <- T[lower.tri(T)]
                }
                m.unco.idx[[offset]] <-     which(tmp > 0)
                x.unco.idx[[offset]] <- tmp[which(tmp > 0)]
            } else {
                m.unco.idx[[offset]] <- m.free.idx[[offset]]
                x.unco.idx[[offset]] <- x.free.idx[[offset]]
            }

            # 3. general mapping between user and GLIST
            tmp[ cbind(REP$row[idx], REP$col[idx]) ] <- partable$id[idx]
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
            tmp[ cbind(REP$row[idx], REP$col[idx]) ] <- startValues[idx]
            if(mmSymmetric[mm]) {
                T <- t(tmp); tmp[lower.tri(tmp)] <- T[lower.tri(T)]
            }

            # representation specific stuff
            if(representation == "LISREL" && mmNames[mm] == "lambda") { 
                ov.dummy.names <- attr(REP, "ov.dummy.names")[[g]]
                # define dummy latent variables
                if(length(ov.dummy.names)) {
                    # in this case, lv.names will be extended with the dummys
                    LV.names <- mmDimNames$psi[[1]]
                    row..idx <- match(ov.dummy.names, ov.names)
                    col..idx <- match(ov.dummy.names, LV.names)
                    tmp[ cbind(row..idx, col..idx)] <- 1.0
                }
            }

            # representation specific stuff (part 2)
            #if(representation == "LISREL" && mmNames[mm] == "psi") {
            #    # this only seems to happen in the categorical case
            #    idx <- which(diag(tmp) == 0.0)
            #    diag(tmp)[idx] <- 1.0
            #}
            
            # assign matrix to GLIST
            GLIST[[offset]] <- tmp
        } # mm
    } # g

    # fixed.x parameters?
    fixed.x <- any(partable$exo > 0L & partable$free == 0L) 
    
    # second check (categorical)
    if(categorical)
        fixed.x <- TRUE



    # constraints

    # 1. simple equality constraints (eg b1 == b2)
    #    capture equality constraints in 'K' matrix
    #    rows are the unconstrained parameters, cols are the unique parameters
    n.unco <- max(partable$unco)
    n.free       <- max(partable$free)
    if(n.free == n.unco) {
        eq.constraints <- FALSE
        K <- matrix(0, 0, 0)
    } else {
        K <- matrix(0, nrow=n.unco, ncol=n.free)
        idx.free <- partable$free[ partable$free > 0 ]
        for(k in 1:n.unco) {
            c <- idx.free[k]
            K[k, c] <- 1
        }
        eq.constraints <- TRUE
    }

    # 2. variable definitions
    #    define a new variable as a arbitrary expression of free parameters
    #
    def.function <- function() NULL
    def.idx <- which(partable$op == ":=")
    if(length(def.idx) > 0L) {
        formals(def.function) <- alist(x=, ...=)
        BODY.txt <- paste("{\n# parameter definitions\n\n")

        lhs.names <- partable$lhs[def.idx]
        def.labels <- all.vars( parse(file="", text=partable$rhs[def.idx]) )
        # remove the ones in lhs.names
        idx <- which(def.labels %in% lhs.names)
        if(length(idx) > 0L) def.labels <- def.labels[-idx]

        # get corresponding 'x' indices
        def.x.idx  <- partable$free[match(def.labels, partable$label)]
        if(any(is.na(def.x.idx))) {
            stop("lavaan ERROR: unknown label(s) in variable definitions: ",
                 paste(def.labels[which(is.na(def.x.idx))], collapse=" "))
        }
        if(any(def.x.idx == 0)) {
            stop("lavaan ERROR: non-free parameter(s) in variable definitions: ",
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
    }

    # 3a. non-trivial equality constraints (linear or nonlinear)
    #     convert to 'ceq(x)' function where 'x' is the (free) parameter vector
    #     and ceq(x) returns the evaluated equality constraints
    #
    #     eg. if b1 + b2 == 2 (and b1 correspond to, say,  x[10] and x[17])
    #         ceq <- function(x) {
    #             out <- rep(NA, 1)
    #             b1 = x[10]; b2 = x[17] 
    #             out[1] <- b1 + b2 - 2
    #         }
    
    ceq.function <- function() NULL
    eq.idx <- which(partable$op == "==")
    if(length(eq.idx) > 0L) {
        formals(ceq.function) <- alist(x=, ...=)
        BODY.txt <- paste("{\nout <- rep(NA, ", length(eq.idx), ")\n", sep="")

        # first come the variable definitions
        if(length(def.idx) > 0L) {
            for(i in 1:length(def.idx)) {
                lhs <- partable$lhs[ def.idx[i] ]
                rhs <- partable$rhs[ def.idx[i] ]
                def.string <- rhs
                # coerce to expression to extract variable names
                def.labels <- all.vars( parse(file="", text=def.string) )
                # get corresponding 'x' indices
                def.x.idx  <- partable$free[match(def.labels, partable$label)]
                def.x.lab  <- paste("x[", def.x.idx, "]",sep="")
                # put both the labels and the expression in the function BODY
                BODY.txt <- paste(BODY.txt,
                    paste(def.labels, "=",def.x.lab, collapse=";"),"\n",
                    lhs, " = ", def.string, "\n", sep="")
            }
        }

        for(i in 1:length(eq.idx)) {
            lhs <- partable$lhs[ eq.idx[i] ]
            rhs <- partable$rhs[ eq.idx[i] ]
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
                def.names <- as.character(partable$lhs[def.idx])
                d.idx <- which(eq.labels %in% def.names)
                if(length(d.idx) > 0) eq.labels <- eq.labels[-d.idx]
            }
            if(length(eq.labels) > 0L) {
                eq.x.idx  <- partable$free[match(eq.labels, partable$label)]
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
        # what to do with NA values? -> return +Inf???
        BODY.txt <- paste(BODY.txt, "out[is.na(out)] <- Inf\n", sep="")
        BODY.txt <- paste(BODY.txt, "return(out)\n}\n", sep="")
        body(ceq.function) <- parse(file="", text=BODY.txt)
        if(debug) { cat("ceq.function = \n"); print(ceq.function); cat("\n") }
    }
    
    # 3b. construct jacobian function 
    #     corresponding with the ceq.function constraints
    ceq.jacobian <- function() NULL
    # TODO!!!
   
    # 4a. non-trivial inequality constraints (linear or nonlinear)
    #     convert to 'cin(x)' function where 'x' is the (free) parameter vector
    #     and cin(x) returns the evaluated inequality constraints
    #
    #     eg. if b1 + b2 > 2 (and b1 correspond to, say,  x[10] and x[17])
    #         cin <- function(x) {
    #             out <- rep(NA, 1)
    #             b1 = x[10]; b2 = x[17] 
    #             out[1] <- b1 + b2 - 2
    #         }
    cin.function <- function() NULL
    ineq.idx <- which(partable$op == ">" | partable$op == "<")
    if(length(ineq.idx) > 0L) {
        formals(cin.function) <- alist(x=, ...=)
        BODY.txt <- paste("{\nout <- rep(NA, ", length(ineq.idx), ")\n", sep="")

        # first come the variable definitions
        if(length(def.idx) > 0L) {
            for(i in 1:length(def.idx)) {
                lhs <- partable$lhs[ def.idx[i] ]
                rhs <- partable$rhs[ def.idx[i] ]
                def.string <- rhs 
                # coerce to expression to extract variable names
                def.labels <- all.vars( parse(file="", text=def.string) )
                # get corresponding 'x' indices
                def.x.idx  <- partable$free[match(def.labels, partable$label)]
                def.x.lab  <- paste("x[", def.x.idx, "]",sep="")
                # put both the labels and the expression in the function BODY
                BODY.txt <- paste(BODY.txt,
                    paste(def.labels, "=",def.x.lab, collapse=";"),"\n",
                    lhs, " = ", def.string, "\n", sep="")
            }
        }

        for(i in 1:length(ineq.idx)) {
            lhs <- partable$lhs[ ineq.idx[i] ]
             op <- partable$op[  ineq.idx[i] ]
            rhs <- partable$rhs[ ineq.idx[i] ]
            if(rhs == "0" && op == ">") {
                ineq.string <- lhs
            } else if(rhs == "0" && op == "<") {
                ineq.string <- paste(rhs, " - (", lhs, ")", sep="")   
            } else if(rhs != "0" && op == ">") {
                ineq.string <- paste(lhs, " - (", rhs, ")", sep="")
            } else if(rhs != "0" && op == "<") {
                ineq.string <- paste(rhs, " - (", lhs, ")", sep="")
            }
            # coerce to expression to extract variable names
            ineq.labels <- all.vars( parse(file="", text=ineq.string) )
            # get corresponding 'x' indices
            if(length(def.idx) > 0L) {
                # remove def.names from ineq.labels
                def.names <- as.character(partable$lhs[def.idx])
                d.idx <- which(ineq.labels %in% def.names)   
                if(length(d.idx) > 0) ineq.labels <- ineq.labels[-d.idx]
            } 
            if(length(ineq.labels) > 0L) {
                ineq.x.idx  <- partable$free[match(ineq.labels, partable$label)]
                if(any(is.na(ineq.x.idx))) {
                   stop("lavaan ERROR: unknown label(s) in inequality constraint: ",
                        paste(ineq.labels[which(is.na(ineq.x.idx))], collapse=" "))
                }
                if(any(ineq.x.idx == 0)) {
                    stop("lavaan ERROR: non-free parameter(s) in inequality constraint: ",
                        paste(ineq.labels[which(ineq.x.idx == 0)], collapse=" "))
                }
                ineq.x.lab  <- paste("x[", ineq.x.idx, "]",sep="")
                # put both the labels and the expression in the function BODY
                BODY.txt <- paste(BODY.txt,
                    paste(ineq.labels, "=", ineq.x.lab, collapse=";"),"\n",
                    "out[", i, "] = ", ineq.string, "\n", sep="")
            } else {
                BODY.txt <- paste(BODY.txt, 
                    "out[", i, "] = ", ineq.string, "\n", sep="")
            }
        }
        # what to do with NA values? -> return +Inf???
        BODY.txt <- paste(BODY.txt, "out[is.na(out)] <- Inf\n", sep="")   
        BODY.txt <- paste(BODY.txt, "return(out)\n}\n", sep="")
        body(cin.function) <- parse(file="", text=BODY.txt)
        if(debug) { cat("cin.function = \n"); print(cin.function); cat("\n") }
    }

    # 4b. construct jacobian function 
    #     corresponding with the cin.function constraints
    cin.jacobian <- function() NULL
    # TODO!!!


    # which free parameters are observed variances?
    #ov.names <- vnames(partable, "ov")
    #x.free.var.idx <-  partable$free[ partable$free & !duplicated(partable$free) &
    #                              partable$lhs %in% ov.names &
    #                              partable$op == "~~" & partable$lhs == partable$rhs ]

    Model <- new("Model",
                 GLIST=GLIST,
                 dimNames=dimNames,
                 isSymmetric=isSymmetric,
                 mmSize=mmSize,
                 representation=representation,
                 meanstructure=meanstructure,
                 categorical=categorical,
                 ngroups=ngroups,
                 nmat=nmat,
                 nvar=nvar,
                 num.idx=num.idx,
                 th.idx=th.idx,
                 nx.free=max(partable$free),
                 nx.unco=max(partable$unco),
                 nx.user=max(partable$id),
                 m.free.idx=m.free.idx,
                 x.free.idx=x.free.idx,
                 m.unco.idx=m.unco.idx,
                 x.unco.idx=x.unco.idx,
                 m.user.idx=m.user.idx,
                 x.user.idx=x.user.idx,
                 x.def.idx=def.idx,
                 x.ceq.idx=eq.idx,
                 x.cin.idx=ineq.idx,
                 #x.free.var.idx=x.free.var.idx,
                 eq.constraints=eq.constraints,
                 eq.constraints.K=K,
                 def.function=def.function,
                 ceq.function=ceq.function,
                 ceq.jacobian=ceq.jacobian,
                 cin.function=cin.function,
                 cin.jacobian=cin.jacobian, 
                 nexo=nexo,
                 fixed.x=fixed.x)

    if(debug) {
         cat("lavaan DEBUG: lavaanModel\n")
         print( str(Model) )
         print( Model@GLIST )
    }

    Model
}
