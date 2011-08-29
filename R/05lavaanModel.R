# constructor of the matrix representation
#
# initial version: YR 22/11/2010

# construct MATRIX representation of the model
Model <- function(user           = NULL,
                  start          = NULL, 
                  representation = "LISREL",
                  debug          = FALSE) {

    # info from user model
    ov.names    <- vnames(user, "ov"); nvar <- length(ov.names)
    lv.names    <- vnames(user, "lv")
    ngroups <- max(user$group)
    meanstructure <- any(user$op == "~1")

    # what if no starting values are provided? 
    if(is.null(start)) {
        start <- user$ustart

        # 0. set everything to 0.0
        start[ is.na(user$ustart) ] <- 0.0

        # 1. factor loadings: set to 1.0
        start[ is.na(user$ustart) & user$op == "=~" ] <- 1.0

        # 2. all lv variances: set to 0.05
        start[ is.na(user$ustart) & user$op == "~~"
                                  & user$lhs %in% lv.names 
                                  & user$lhs == user$rhs ] <- 0.05
 
        # 3. residual ov variances
        start[ is.na(user$ustart) & user$op == "~~"
                                  & user$lhs %in% ov.names
                                  & user$lhs == user$rhs ] <- 1.0
    }
 
    stopifnot(length(start) == nrow(user))

    # only representation = "LISREL" for now
    stopifnot(representation == "LISREL")


    # capture equality constraints in 'K' matrix
    # rows are the unconstrained parameters, cols are the unique parameters
    n.free.uncon <- max(user$free.uncon)
    n.free       <- max(user$free)
    if(n.free == n.free.uncon) {
        eq.constraints <- FALSE
        K <- matrix(0, 0, 0)
    } else {
        K <- matrix(0, nrow=n.free.uncon, ncol=n.free)
        #####
        #####     FIXME !
        #####
        idx.free <- user$free[ user$free > 0 ]
        for(k in 1:n.free.uncon) {
            c <- idx.free[k]
            K[k, c] <- 1
        }
        eq.constraints <- TRUE
    }

    # select model matrices
    if(representation == "LISREL") {
        REP <- representation.LISREL(user, target=NULL, extra=TRUE)
    } else {
        stop("lavaan ERROR: only representation \"LISREL\" has been implemented.")
    }
    REP$mat <- REP$mat
    REP$row <- REP$row
    REP$col <- REP$col

    # summarize representation
    mmNumber    <- attr(REP, "mmNumber") 
    mmNames     <- attr(REP, "mmNames")
    mmSymmetric <- attr(REP, "mmSymmetric")
    mmDimNames  <- attr(REP, "mmDimNames")
    mmRows      <- attr(REP, "mmRows")
    mmCols      <- attr(REP, "mmCols")

    # matrix size (unique elements only)
    mm.size <- integer(mmNumber)
    for(mm in 1:mmNumber) {
        if(mmSymmetric[mm]) {
            N <- mmRows[mm]
            mm.size[mm] <- as.integer(N*(N+1)/2)
        } else {
            mm.size[mm] <- as.integer(mmRows[mm] * mmCols[mm])
        }
    }


    # fill in GLIST
    GLIST <- vector(mode="list", length=mmNumber*ngroups)
    names(GLIST) <- rep(mmNames, times=ngroups)

    m.free.idx <- m.unco.idx <- m.user.idx <- 
        vector(mode="list", length=mmNumber*ngroups)
    x.free.idx <- x.unco.idx <- x.user.idx <-
        vector(mode="list", length=mmNumber*ngroups)
    

    for(g in 1:ngroups) {
        for(mm in 1:mmNumber) {

            # offset in GLIST
            offset <- mmNumber*(g - 1L) + mm

            # select elements for this matrix
            idx <- which(user$group == g & REP$mat == mmNames[mm]) 

            # create empty `pattern' matrix
            # FIXME: one day, we may want to use sparse matrices...
            #        but they may not slow things down!
            tmp <- matrix(0L, nrow=mmRows[mm],
                              ncol=mmCols[mm])

            # 1. first assign free values only, to get vector index
            #    -> to be used in computeObjective
            tmp[ cbind(REP$row[idx], REP$col[idx]) ] <- user$free[idx]
            if(mmSymmetric[mm]) {
                # NOTE: we assume everything is in the UPPER tri!
                T <- t(tmp); tmp[lower.tri(tmp)] <- T[lower.tri(T)]
            }
            m.free.idx[[offset]] <-     which(tmp > 0)
            x.free.idx[[offset]] <- tmp[which(tmp > 0)]

            # 2. if equality constraints, unconstrained free parameters
            #    -> to be used in computeGradient
            if(eq.constraints) {
                tmp[ cbind(REP$row[idx], REP$col[idx]) ] <- user$free.uncon[idx]
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
            tmp[ cbind(REP$row[idx], REP$col[idx]) ] <- user$id[idx]
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
            tmp[ cbind(REP$row[idx], REP$col[idx]) ] <- start[idx]
            if(mmSymmetric[mm]) {
                T <- t(tmp); tmp[lower.tri(tmp)] <- T[lower.tri(T)]
            }

            # representation specific stuff
            if(representation == "LISREL" && mmNames[mm] == "lambda") { 
                ov.dummy.names <- attr(REP, "ov.dummy.names")
                # define dummy latent variables
                if(length(ov.dummy.names)) {
                    # in this case, lv.names will be extended with the dummys
                    LV.names <- mmDimNames$psi[[1]]
                    row..idx <- match(ov.dummy.names, ov.names)
                    col..idx <- match(ov.dummy.names, LV.names)
                    tmp[ cbind(row..idx, col..idx)] <- 1.0
                }
            }
            
            # assign matrix to GLIST
            GLIST[[offset]] <- tmp
        } # mm
    } # g

    # fixed.x parameters?
    fixed.x <- any(user$fixed.x > 1 & user$free == 0)



    # constraints

    # 1. simple equality constraints (eg b1 == b2)
    #    capture equality constraints in 'K' matrix
    #    rows are the unconstrained parameters, cols are the unique parameters
    n.free.uncon <- max(user$free.uncon)
    n.free       <- max(user$free)
    if(n.free == n.free.uncon) {
        eq.constraints <- FALSE
        K <- matrix(0, 0, 0)
    } else {
        K <- matrix(0, nrow=n.free.uncon, ncol=n.free)
        idx.free <- user$free[ user$free > 0 ]
        for(k in 1:n.free.uncon) {
            c <- idx.free[k]
            K[k, c] <- 1
        }
        eq.constraints <- TRUE
    }

    # 2. variable definitions
    #    define a new variable as a arbitrary expression of free parameters
    #
    def.function <- function(x, ...) NULL
    def.idx <- which(user$op == ":=")
    if(length(def.idx) > 0L) {
        lhs.names <- character(length(def.idx))
        BODY.txt <- paste("{\nout <- rep(NA, ", length(def.idx), ")\n", sep="")
        for(i in 1:length(def.idx)) {
            lhs <- user$lhs[ def.idx[i] ]
            rhs <- user$rhs[ def.idx[i] ]
            lhs.names[i] <- as.character(lhs)
            def.string <- rhs # expression must be on the right-hand side
            # coerce to expression to extract variable names
            def.labels <- all.vars( parse(file="", text=def.string) )
            # get corresponding 'x' indices
            def.x.idx  <- user$free[match(def.labels, user$label)]
            if(any(is.na(def.x.idx))) {
                stop("lavaan ERROR: unknown label(s) in variable definition: ",
                     paste(def.labels[which(is.na(def.x.idx))], collapse=" "))
            }
            if(any(def.x.idx == 0)) {
                stop("lavaan ERROR: non-free parameter(s) in variable definition: ",
                    paste(def.labels[which(def.x.idx == 0)], collapse=" "))
            }
            def.x.lab  <- paste("x[", def.x.idx, "]",sep="")
            # put both the labels and the expression in the function BODY
            BODY.txt <- paste(BODY.txt,
                paste(def.labels, "=",def.x.lab, collapse=";"),"\n",
                "out[", i, "] = ", def.string, "\n", sep="")
        }
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
    
    ceq.function <- function(x, ...) NULL
    eq.idx <- which(user$op == "==")
    if(length(eq.idx) > 0L) {
        BODY.txt <- paste("{\nout <- rep(NA, ", length(eq.idx), ")\n", sep="")

        # first come the variable definitions
        if(length(def.idx) > 0L) {
            for(i in 1:length(def.idx)) {
                lhs <- user$lhs[ def.idx[i] ]
                rhs <- user$rhs[ def.idx[i] ]
                def.string <- rhs
                # coerce to expression to extract variable names
                def.labels <- all.vars( parse(file="", text=def.string) )
                # get corresponding 'x' indices
                def.x.idx  <- user$free[match(def.labels, user$label)]
                def.x.lab  <- paste("x[", def.x.idx, "]",sep="")
                # put both the labels and the expression in the function BODY
                BODY.txt <- paste(BODY.txt,
                    paste(def.labels, "=",def.x.lab, collapse=";"),"\n",
                    lhs, " = ", def.string, "\n", sep="")
            }
        }

        for(i in 1:length(eq.idx)) {
            lhs <- user$lhs[ eq.idx[i] ]
            rhs <- user$rhs[ eq.idx[i] ]
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
                def.names <- as.character(user$lhs[def.idx])
                d.idx <- which(eq.labels %in% def.names)
                if(length(d.idx) > 0) eq.labels <- eq.labels[-d.idx]
            }
            if(length(eq.labels) > 0L) {
                eq.x.idx  <- user$free[match(eq.labels, user$label)]
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
    ceq.jacobian <- function(x, ...) NULL
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
    cin.function <- function(x, ...) NULL
    ineq.idx <- which(user$op == ">" | user$op == "<")
    if(length(ineq.idx) > 0L) {
        BODY.txt <- paste("{\nout <- rep(NA, ", length(ineq.idx), ")\n", sep="")

        # first come the variable definitions
        if(length(def.idx) > 0L) {
            for(i in 1:length(def.idx)) {
                lhs <- user$lhs[ def.idx[i] ]
                rhs <- user$rhs[ def.idx[i] ]
                def.string <- rhs 
                # coerce to expression to extract variable names
                def.labels <- all.vars( parse(file="", text=def.string) )
                # get corresponding 'x' indices
                def.x.idx  <- user$free[match(def.labels, user$label)]
                def.x.lab  <- paste("x[", def.x.idx, "]",sep="")
                # put both the labels and the expression in the function BODY
                BODY.txt <- paste(BODY.txt,
                    paste(def.labels, "=",def.x.lab, collapse=";"),"\n",
                    lhs, " = ", def.string, "\n", sep="")
            }
        }

        for(i in 1:length(ineq.idx)) {
            lhs <- user$lhs[ ineq.idx[i] ]
             op <- user$op[  ineq.idx[i] ]
            rhs <- user$rhs[ ineq.idx[i] ]
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
                def.names <- as.character(user$lhs[def.idx])
                d.idx <- which(ineq.labels %in% def.names)   
                if(length(d.idx) > 0) ineq.labels <- ineq.labels[-d.idx]
            } 
            if(length(ineq.labels) > 0L) {
                ineq.x.idx  <- user$free[match(ineq.labels, user$label)]
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
    cin.jacobian <- function(x, ...) NULL
    # TODO!!!


    # which free parameters are observed variances?
    x.free.var.idx <-  user$free[ user$free & !duplicated(user$free) &
                                  user$lhs %in% ov.names &
                                  user$op == "~~" & user$lhs == user$rhs ]

    Model <- new("Model",
                 GLIST=GLIST,
                 dimNames=rep(mmDimNames, ngroups),
                 isSymmetric=rep(mmSymmetric, ngroups),
                 mmSize=rep(mm.size, ngroups),
                 representation=representation,
                 meanstructure=meanstructure,
                 ngroups=ngroups,
                 nmat=mmNumber,
                 nvar=nvar,
                 nx.free=max(user$free),
                 nx.unco=max(user$free.uncon),
                 nx.user=max(user$id),
                 m.free.idx=m.free.idx,
                 x.free.idx=x.free.idx,
                 m.unco.idx=m.unco.idx,
                 x.unco.idx=x.unco.idx,
                 m.user.idx=m.user.idx,
                 x.user.idx=x.user.idx,
                 x.def.idx=def.idx,
                 x.ceq.idx=eq.idx,
                 x.cin.idx=ineq.idx,
                 x.free.var.idx=x.free.var.idx,
                 eq.constraints=eq.constraints,
                 eq.constraints.K=K,
                 def.function=def.function,
                 ceq.function=ceq.function,
                 ceq.jacobian=ceq.jacobian,
                 cin.function=cin.function,
                 cin.jacobian=cin.jacobian, 
                 fixed.x=fixed.x)

    if(debug) {
         cat("lavaan DEBUG: lavaanModel\n")
         print( str(Model) )
    }

    Model
}
