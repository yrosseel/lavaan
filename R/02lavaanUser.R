# constructor for the lavParTable model description
#
# initial version: YR 22/05/2009
#  major revision: YR 02/11/2010: - FLATTEN the model syntax and turn it into a
#                                   data.frame, with a "modifiers" attribute
#                                 - add default elements here
#                                 - check for duplicate elements
#                                 - allow for every possible model...
#                                 - since 0.4-5
#                                 - the end result is a full description of
#                                   a model (but no matrix representation)


lavaanify <- lavParTable <- function(

                      model           = NULL, 
                      meanstructure   = FALSE,
                      int.ov.free     = FALSE,
                      int.lv.free     = FALSE,
                      orthogonal      = FALSE, 
                      std.lv          = FALSE,
                      fixed.x         = TRUE,
                      constraints     = NULL,

                      auto            = FALSE,
                      model.type      = "sem",
                      auto.fix.first  = FALSE,
                      auto.fix.single = FALSE,
                      auto.var        = FALSE,
                      auto.cov.lv.x   = FALSE,
                      auto.cov.y      = FALSE,
                      auto.th         = FALSE,
                      auto.delta      = FALSE,

                      varTable        = NULL,
                      ngroups         = 1L,
                      group.equal     = NULL,
                      group.partial   = NULL,
                      group.w.free    = FALSE,
                      debug           = FALSE,
                      warn            = TRUE,
                      
                      as.data.frame.   = TRUE) {


    # check if model is already FLAT or a full parameter table
    if(is.list(model) && !is.null(model$lhs)) {
        if(is.null(model$mod.idx)) {
            warning("lavaan WARNING: input already looks like a parameter table; returning as is")
            return(model)
        } else {
            FLAT <- model
        }
    } else {
        # parse the model syntax and flatten the user-specified model
        # return a data.frame, where each line is a model element (rhs, op, lhs)
        FLAT <- lavParseModelString(model.syntax=model, warn=warn, debug=FALSE)
    }
    # user-specified *modifiers* are returned as an attribute
    MOD  <- attr(FLAT, "modifiers"); attr(FLAT, "modifiers") <- NULL
    # user-specified *constraints* are returned as an attribute
    CON  <- attr(FLAT, "constraints"); attr(FLAT, "constraints") <- NULL

    # extra constraints?
    if(!is.null(constraints) && nchar(constraints) > 0L) {
        FLAT2 <- lavParseModelString(model.syntax=constraints, warn=warn)
        CON2 <- attr(FLAT2, "constraints"); rm(FLAT2)
        CON <- c(CON, CON2)
    }
   
    if(debug) {
        cat("[lavaan DEBUG]: FLAT (flattened user model):\n")
        print(FLAT)
        cat("[lavaan DEBUG]: MOD (modifiers):\n")
        print( str(MOD) )
        cat("[lavaan DEBUG]: CON (constraints):\n")
        print( str(CON) )
    }

    # check for wrongly specified variances/covariances/intercepts
    # of exogenous variables in model syntax (if fixed.x=TRUE)
    if(fixed.x) { # we ignore the groups here!
        # we only call this function for the warning message
        tmp <- vnames(FLAT, "ov.x", warn=TRUE); rm(tmp)
    }

    # auto=TRUE?
    if(auto && model.type == "sem") { # mimic sem/cfa auto behavior
        if(model.type == "sem") {
            int.ov.free     = TRUE
            int.lv.free     = FALSE
            auto.fix.first  = !std.lv
            auto.fix.single = TRUE
            auto.var        = TRUE
            auto.cov.lv.x   = TRUE
            auto.cov.y      = TRUE
            auto.th         = TRUE
            auto.delta      = TRUE
        } else 
        
        if(model.type == "growth") {
            model.type      = "growth"
            int.ov.free     = FALSE
            int.lv.free     = TRUE
            auto.fix.first  = !std.lv
            auto.fix.single = TRUE
            auto.var        = TRUE
            auto.cov.lv.x   = TRUE
            auto.cov.y      = TRUE
            auto.th         = TRUE    
            auto.delta      = TRUE
        }
    }

    # check for meanstructure
    if(any(FLAT$op == "~1")) meanstructure <- TRUE

    # check for group identifiers in the syntax (op = ":")
    n.group.flat <- length(which(FLAT$op == ":"))
    if(n.group.flat > 0L && n.group.flat != ngroups) {
        stop("lavaan ERROR: syntax defines ", n.group.flat, " groups; ",
             "data suggests ", ngroups, " groups")
    }

    # for each `group' in n.group.flat, produce a USER
    if(n.group.flat > 0L) {
        # split the FLAT data.frame per `group', create LIST
        # for each `group', and bind them together
        FLAT <- as.data.frame(FLAT, stringsAsFactors=FALSE)
        group.op.idx <- c(which(FLAT$op == ":"), nrow(FLAT)+1L)
        for(g in 1:n.group.flat) {
            FLAT.group <- FLAT[(group.op.idx[g]+1L):(group.op.idx[g+1]-1L),]    
            LIST.group <- getLIST(FLAT.group, meanstructure = meanstructure, 
                int.ov.free = int.ov.free, int.lv.free = int.lv.free,
                orthogonal = orthogonal, std.lv = std.lv, fixed.x = fixed.x,
                auto.fix.first = auto.fix.first, 
                auto.fix.single = auto.fix.single,
                auto.var = auto.var, auto.cov.lv.x = auto.cov.lv.x,
                auto.cov.y = auto.cov.y, auto.th = auto.th, 
                auto.delta = auto.delta, 
                varTable = varTable, group.equal = NULL, 
                group.w.free = group.w.free, ngroups = 1L)
            LIST.group <- as.data.frame(LIST.group, stringsAsFactors=FALSE)
            if(g == 1L) {
                LIST <- LIST.group
            } else {
                LIST.group$group <- rep(g, length(LIST.group$lhs))
                LIST.group$id <- LIST.group$id + max(LIST$id)
                LIST <- rbind(LIST, LIST.group)
            }
        }
        LIST <- as.list(LIST)
    } else {
        LIST <- getLIST(FLAT, meanstructure = meanstructure, 
            int.ov.free = int.ov.free, int.lv.free = int.lv.free,
            orthogonal = orthogonal, std.lv = std.lv, fixed.x = fixed.x,
            auto.fix.first = auto.fix.first, auto.fix.single = auto.fix.single,
            auto.var = auto.var, auto.cov.lv.x = auto.cov.lv.x,
            auto.cov.y = auto.cov.y, auto.th = auto.th, 
            auto.delta = auto.delta,
            varTable = varTable, group.equal = group.equal, 
            group.w.free = group.w.free,
            ngroups = ngroups)
    }        
    if(debug) {
        cat("[lavaan DEBUG]: parameter LIST:\n")
        print( as.data.frame(LIST, stringsAsFactors=FALSE) )
    }

    # apply user-specified modifiers
    if(length(MOD)) {
        for(el in 1:length(MOD)) {
            idx <- which(LIST$mod.idx == el) # for each group
            MOD.fixed <- MOD[[el]]$fixed  
            MOD.start <- MOD[[el]]$start
            MOD.label <- MOD[[el]]$label 

            # check for single argument if multiple groups
            if(ngroups > 1L && length(idx) > 1L) {
                # Ok, this is not very consistent:
                # A) here we force same behavior across groups
                if(length(MOD.fixed) == 1L) MOD.fixed <- rep(MOD.fixed, ngroups)
                if(length(MOD.start) == 1L) MOD.start <- rep(MOD.start, ngroups)
                # B) here we do NOT! otherwise, it would imply an equality 
                #                    constraint...
                #    except if group.equal="loadings"!
                if(length(MOD.label) == 1L) {
                    if("loadings" %in% group.equal) {
                        MOD.label <- rep(MOD.label, ngroups)
                    } else {
                        MOD.label <- c(MOD.label, rep("", (ngroups-1L)) )
                    }
                }
            }

            # check for wrong number of arguments if multiple groups
            nidx <- length(idx)
            if( (!is.null(MOD.fixed) && nidx != length(MOD.fixed)) ||
                (!is.null(MOD.start) && nidx != length(MOD.start)) ||
                (!is.null(MOD.label) && nidx != length(MOD.label)) ) {
                el.idx <- which(LIST$mod.idx == el)[1L]
                stop("lavaan ERROR: wrong number of arguments in modifier (",
                    paste(MOD.label, collapse=","), ") of element ", 
                    LIST$lhs[el.idx], LIST$op[el.idx], LIST$rhs[el.idx])
            }

            # apply modifiers
            if(!is.null(MOD.fixed)) {
                na.idx <- which(is.na(MOD.fixed))
                not.na.idx <- which(!is.na(MOD.fixed))
                LIST$ustart[idx][not.na.idx] <- MOD.fixed[not.na.idx]
                LIST$free[  idx][not.na.idx] <- 0L
                LIST$free[  idx][na.idx]     <- 1L # eg factor loading
            }
            if(!is.null(MOD.start)) {
                LIST$ustart[idx] <- MOD.start
            }
            if(!is.null(MOD.label)) {
                LIST$label[idx] <- MOD.label
            }
        }
    }
    # remove mod.idx column
    LIST$mod.idx <- NULL

    # check if CON contains *simple* equality constraints (eg b1 == b2)
    # FIXME!!!
    # b1 == b3
    # b2 == b3 does not work
    # better a general approach for linear constraints!
    #if(length(CON) > 0L) {
    #    el.idx <- integer(0L)
    #    for(el in 1:length(CON)) {
    #        if(CON[[el]]$op == "==") {
    #            lhs.idx <- which(LIST$label == CON[[el]]$lhs)
    #            rhs.idx <- which(LIST$label == CON[[el]]$rhs)
    #            if(length(lhs.idx) && length(rhs.idx)) {
    #                # flag this constraint (to be removed)
    #                el.idx <- c(el.idx, el)
    #                # fill in equal and fix
    #                if(LIST$label[lhs.idx] == "") {
    #                    LIST$label[rhs.idx] <- LIST$label[lhs.idx]
    #                } else {
    #                    LIST$label[rhs.idx] <- LIST$label[lhs.idx]
    #                }
    #                LIST$free[ rhs.idx] <- 0L
                     # FIXME: what is needed is to replace all occurences
                     #        of rhs.idx by lhs.idx in CON!!!
    #            }
    #        }
    #    }
    #    if(length(el.idx) > 0L) CON <- CON[-el.idx]
    #}

    # get 'virtual' parameter labels
    LABEL <- getParameterLabels(partable=LIST, group.equal=group.equal,
                                group.partial=group.partial)
    #cat("DEBUG: label after getParameterLabels:\n"); print(LABEL); cat("\n")
    #cat("DEBUG: eq.id after group.equal:\n"); print(LIST$eq.id); cat("\n")
    #cat("DEBUG: LIST$label:\n"); print(LIST$label); cat("\n")

    # handle user-specified equality constraints
    # insert 'target/reference' id's in eq.id/label columns
    # so that the eq.id column reflects sets of equal parameters
    # the 'reference' parameter has free > 0; the 'equal to' pars free == 0
    idx.eq.label <- which(duplicated(LABEL))
    if(length(idx.eq.label) > 0L) {
        for(idx in idx.eq.label) {
            eq.label <- LABEL[idx]
            ref.idx <- which(LABEL == eq.label)[1L] # the first one only
            # set eq.id equal
            LIST$eq.id[ref.idx] <- LIST$eq.id[idx] <- ref.idx
            # fix target
            LIST$free[idx] <- 0L

            # special case: check if there are any more instances 
            # of idx in  LIST$eq.id (perhaps due to group.equal)
            idx.all <- which(LIST$eq.id == idx)
            if(length(idx.all) > 0L) {
                ref.idx.all <- which(LIST$eq.id == ref.idx)
                LIST$label[ref.idx.all] <- eq.label
                LIST$eq.id[idx.all] <- ref.idx
                LIST$free[idx.all] <- 0L  # fix!
                LIST$label[idx.all] <- eq.label
            }
            # special case: ref.idx is a fixed parameter
            if(LIST$free[ref.idx] == 0L) {
                ref.idx.all <- which(LIST$eq.id == ref.idx)
                LIST$ustart[ref.idx.all] <- LIST$ustart[ref.idx]
                LIST$eq.id[ref.idx.all] <- 0L
            }
        }
    }
 
    #cat("DEBUG: eq.id after eq.id:\n"); print(LIST$eq.id); cat("\n")    
    #cat("DEBUG: LIST$label:\n"); print(LIST$label); cat("\n")

    # count free parameters
    idx.free <- which(LIST$free > 0)
    LIST$free[idx.free] <- 1:length(idx.free)

    # 2. add free counter to this element
    idx.equal <- which(LIST$eq.id > 0)
    LIST$free[idx.equal] <- LIST$free[ LIST$eq.id[idx.equal] ]

    # 3. which parameters would be free without equality constraints?
    idx.unco <- which(LIST$free > 0)
    LIST$unco[idx.unco] <- 1:length(idx.unco)


    # handle constraints (if any) (NOT per group, but overall - 0.4-11)
    if(length(CON) > 0L) {
        #cat("DEBUG:\n"); print(CON)
        lhs = unlist(lapply(CON, "[[", "lhs"))
         op = unlist(lapply(CON, "[[",  "op"))
        rhs = unlist(lapply(CON, "[[", "rhs"))
        LIST$id         <- c(LIST$id,         max(LIST$id) + 1:length(lhs) )
        LIST$lhs        <- c(LIST$lhs,        lhs)
        LIST$op         <- c(LIST$op,         op)
        LIST$rhs        <- c(LIST$rhs,        rhs)
        LIST$user       <- c(LIST$user,       rep(1L, length(lhs)) )
        LIST$group      <- c(LIST$group,      rep(0L, each=length(CON)))
        LIST$free       <- c(LIST$free,       rep(0L, length(lhs)) )
        LIST$ustart     <- c(LIST$ustart,     rep(as.numeric(NA), length(lhs)))
        LIST$exo        <- c(LIST$exo,        rep(0L, length(lhs)) )
        LIST$label      <- c(LIST$label,      rep("",  length(lhs)) )
        LIST$eq.id      <- c(LIST$eq.id,      rep(0L,  length(lhs)) )
        LIST$unco       <- c(LIST$unco,       rep(0L,  length(lhs)) )
    }

    # put lhs of := elements in label column
    def.idx <- which(LIST$op == ":=")
    LIST$label[def.idx] <- LIST$lhs[def.idx]


    if(debug) { 
        cat("[lavaan DEBUG] lavParTable\n")
        print( as.data.frame(LIST) )
    }

    # data.frame?
    if(as.data.frame.) LIST <- as.data.frame(LIST, stringsAsFactors = FALSE)

    LIST
}


lavParseModelString <- function(model.syntax = '', as.data.frame. = FALSE,
                             warn = TRUE, debug = FALSE) {
  
    # check for empty syntax
    if(length(model.syntax) == 0) {
        stop("lavaan ERROR: empty model syntax")
    }
  
    # replace semicolons with newlines prior to split
    model.syntax <- gsub(";", "\n", model.syntax, fixed=TRUE)
  
    # remove comments prior to split. 
    # Match from comment character to newline, but don't eliminate newline
    model.syntax <- gsub("[#!].*(?=\n)","", model.syntax, perl=TRUE)
  
    #remove whitespace prior to split
    model.syntax <- gsub("[ \t]+", "", model.syntax, perl=TRUE)
    # remove any occurrence of >= 2 consecutive newlines to eliminate \
    # blank statements; this retains a blank newline at the beginning, 
    # if such exists, but parser will not choke because of start.idx
    model.syntax <- gsub("\n{2,}", "\n", model.syntax, perl=TRUE)
  
    # break up in lines 
    model <- unlist( strsplit(model.syntax, "\n") )
    
    # check for multi-line formulas: they contain no "~" or "=" character
    # but before we do that, we remove all modifiers
    # to avoid confusion with for example equal("f1=~x1") statements
    model.simple <- gsub("\\(.*\\)\\*", "MODIFIER*", model)
  
    start.idx <- grep("[~=<>:|]", model.simple)
    end.idx <- c( start.idx[-1]-1, length(model) )
    model.orig    <- model
    model <- character( length(start.idx) )
    for(i in 1:length(start.idx)) {
        model[i] <- paste(model.orig[start.idx[i]:end.idx[i]], collapse="")
    }
  
    # ok, in all remaining lines, we should have a '~' operator
    # OR one of '=', '<', '>', '|' outside the ""
    model.simple <- gsub("\\\".[^\\\"]*\\\"", "LABEL", model)
    idx.wrong <- which(!grepl("[~=<>:|]", model.simple))
    if(length(idx.wrong) > 0) {
        cat("lavaan: missing operator in formula(s):\n")
        print(model[idx.wrong])
        stop("lavaan ERROR: syntax error in lavaan model syntax")
    }

    # but perhaps we have a '+' as the first character?
    idx.wrong <- which(grepl("^\\+", model))
    if(length(idx.wrong) > 0) {
        cat("lavaan: some formula(s) start with a plus (+) sign:\n")
        print(model[idx.wrong])
        stop("lavaan ERROR: syntax error in lavaan model syntax")
    }
  
  
    # main operation: flatten formulas into single bivariate pieces
    # with a left-hand-side (lhs), an operator (eg "=~"), and a 
    # right-hand-side (rhs)
    # both lhs and rhs can have a modifier 
    # (but we ignore the lhs modifier for now)
    FLAT.lhs         <- character(0)
    #FLAT.lhs.mod    <- character(0)
    FLAT.op          <- character(0)
    FLAT.rhs         <- character(0)
    FLAT.rhs.mod.idx <- integer(0)
    FLAT.group       <- integer(0)    # keep track of groups using ":" operator
  
    FLAT.fixed       <- character(0)  # only for display purposes! 
    FLAT.start       <- character(0)  # only for display purposes!
    FLAT.label       <- character(0)  # only for display purposes!
    FLAT.idx <- 0L
    MOD.idx  <- 0L
    CON.idx  <- 0L
    MOD <- vector("list", length=0L)
    CON <- vector("list", length=0L)
    GRP <- 1L
    GRP_OP <- FALSE
    for(i in 1:length(model)) {
        x <- model[i]
        if(debug) {
           cat("formula to parse:\n"); print(x); cat("\n")
        }
    
        # 1. which operator is used?
        line.simple <- gsub("\\\".[^\\\"]*\\\"", "LABEL", x)
        # "=~" operator?
        if(grepl("=~", line.simple, fixed=TRUE)) {
            op <- "=~"
        # "<~" operator?
        } else if(grepl("<~", line.simple, fixed=TRUE)) {
            op <- "<~"
        } else if(grepl("~*~", line.simple, fixed=TRUE)) {
            op <- "~*~"
        # "~~" operator?
        } else if(grepl("~~", line.simple, fixed=TRUE)) {
            op <- "~~"
        # "~" operator?
        } else if(grepl("~", line.simple, fixed=TRUE)) {
            op <- "~"           
        # "==" operator?
        } else if(grepl("==", line.simple, fixed=TRUE)) {
            op <- "=="  
        # "<" operator?
        } else if(grepl("<", line.simple, fixed=TRUE)) {
            op <- "<"
        # ">" operator?
        } else if(grepl(">", line.simple, fixed=TRUE)) {
            op <- ">"
        # ":=" operator?
        } else if(grepl(":=", line.simple, fixed=TRUE)) {
            op <- ":="
        # ":" operator?
        } else if(grepl(":", line.simple, fixed=TRUE)) {
            op <- ":"
        } else if(grepl("|", line.simple, fixed=TRUE)) {
            op <- "|"
        } else {
            stop("unknown operator in ", model[i])
        }
    
        # 2. split by operator (only the *first* occurence!)
        # check first if equal/label modifier has been used on the LEFT!
        if(substr(x,1,5) == "label") 
            stop("label modifier can not be used on the left-hand side of the operator")
        if(op == "|") {
            op.idx <- regexpr("\\|", x)
        } else if(op == "~*~") {
            op.idx <- regexpr("~\\*~", x)    
        } else {
            op.idx <- regexpr(op, x)
        }
        lhs <- substr(x, 1L, op.idx-1L)
        # fix for 'NA' names in lhs; not likely to happen to ov.names
        # since 'NA' is not a valid name for list elements/data.frame columns
        if(lhs == "NA") lhs <- "NA."
        rhs <- substr(x, op.idx+attr(op.idx, "match.length"), nchar(x))
    
        # 2b. if operator is "==" or "<" or ">" or ":=", put it in CON
        if(op == "==" || op == "<" || op == ">" || op == ":=") {
            # remove quotes, if any
            lhs <- gsub("\\\"", "", lhs)
            rhs <- gsub("\\\"", "", rhs)
            CON.idx <- CON.idx + 1L
            CON[[CON.idx]] <- list(op=op, lhs=lhs, rhs=rhs)
            next
        }
    
        # 2c if operator is ":", put it in GRP
        if(op == ":") {
            FLAT.idx <- FLAT.idx + 1L
            FLAT.lhs[FLAT.idx] <- lhs
            FLAT.op[ FLAT.idx] <- op
            FLAT.rhs[FLAT.idx] <- ""
            FLAT.fixed[FLAT.idx] <- ""
            FLAT.start[FLAT.idx] <- ""
            FLAT.label[FLAT.idx] <- ""
            FLAT.rhs.mod.idx[FLAT.idx] <- 0L
            if(GRP_OP) {
                GRP <- GRP + 1L
            }
            FLAT.group[FLAT.idx] <- GRP
            GRP_OP <- TRUE
            next
        }
    
        # 3. parse left hand
        #    lhs modifiers will be ignored for now
        lhs.formula <- as.formula(paste("~",lhs))
        out <- parse.rhs(rhs=lhs.formula[[2L]])
        lhs.names <- names(out)
        # check if we have modifiers
        if(sum(sapply(out, length)) > 0L) {
            warning("lavaan WARNING: left-hand side of formula below contains modifier:\n", x,"\n")
        }
    
        # 4. parse rhs (as rhs of a single-sided formula)

        # new 0.5-12: before we do this, replace '0.2?' by 'start(0.2)*'
        # requested by the simsem folks
        rhs <- gsub('([0-9]*\\.*[0-9]*)\\?',"start(\\1)\\*",rhs)
        rhs.formula <- as.formula(paste("~",rhs))
        out <- parse.rhs(rhs=rhs.formula[[2L]],op=op)

        if(debug) print(out)
    
        # for each lhs element
        for(l in 1:length(lhs.names)) {
      
            # for each rhs element
            for(j in 1:length(out)) {
        
                # catch intercepts
                if(names(out)[j] == "intercept") {
                    if(op == "~") {
                        rhs.name <- ""
                    } else {
                        stop("lavaan ERROR: right-hand side of formula contains an intercept, but operator is \"", op, "\" in: ", x)
                    }
                } else {
                    rhs.name <- names(out)[j]
                }
        
                # check if we not already have this combination (in this group)
                # 1. asymmetric (=~, ~, ~1)
                if(op != "~~") {
                    idx <- which(FLAT.lhs == lhs.names[l] &
                                 FLAT.op  == op &
                                 FLAT.group == GRP &
                                 FLAT.rhs == rhs.name)
                    if(length(idx) > 0L) {
                        stop("lavaan ERROR: duplicate model element in: ", model[i])
                    }
                } else {
                    # 2. symmetric (~~)
                    idx <- which(FLAT.lhs == rhs.name &
                                 FLAT.op  == "~~" &
                                 FLAT.group == GRP &
                                 FLAT.rhs == lhs.names[l])
                    if(length(idx) > 0L) {
                        stop("lavaan ERROR: duplicate model element in: ", model[i])
                    }
                }
                FLAT.idx <- FLAT.idx + 1L
                FLAT.lhs[FLAT.idx] <- lhs.names[l]
                FLAT.op[ FLAT.idx] <- op
                FLAT.rhs[FLAT.idx] <- rhs.name
                FLAT.group[FLAT.idx] <- GRP
                FLAT.fixed[FLAT.idx] <- ""
                FLAT.start[FLAT.idx] <- ""
                FLAT.label[FLAT.idx] <- ""
        
                mod <- list()
                rhs.mod <- 0L
                if(length(out[[j]]$fixed) > 0L) {
                    mod$fixed <- out[[j]]$fixed
                    FLAT.fixed[FLAT.idx] <- paste(mod$fixed, collapse=";")
                    rhs.mod <- 1L
                }
                if(length(out[[j]]$start) > 0L) {
                    mod$start <- out[[j]]$start
                    FLAT.start[FLAT.idx] <- paste(mod$start, collapse=";")
                    rhs.mod <- 1L
                }
                if(length(out[[j]]$label) > 0L) {
                    mod$label <- out[[j]]$label
                    FLAT.label[FLAT.idx] <- paste(mod$label, collapse=";")
                    rhs.mod <- 1L
                }
                if(op == "~1" && rhs == "0") {
                    mod$fixed <- 0
                    FLAT.fixed[FLAT.idx] <- paste(mod$fixed, collapse=";")
                    rhs.mod <- 1L
                }
                if(op == "=~" && rhs == "0") {
                    mod$fixed <- 0
                    FLAT.rhs[FLAT.idx] <- FLAT.lhs[FLAT.idx]
                    FLAT.fixed[FLAT.idx] <- paste(mod$fixed, collapse=";")
                    rhs.mod <- 1L
                }
        
                FLAT.rhs.mod.idx[FLAT.idx] <- rhs.mod
        
                if(rhs.mod > 0L) {
                    MOD.idx <- MOD.idx + 1L
                    MOD[[MOD.idx]] <- mod
                }
           } # rhs elements
        } # lhs elements
    } # model elements
  
    # enumerate modifier indices
    mod.idx <- which(FLAT.rhs.mod.idx > 0L)
    FLAT.rhs.mod.idx[ mod.idx ] <- 1:length(mod.idx)
  
    FLAT <- list(lhs=FLAT.lhs, op=FLAT.op, rhs=FLAT.rhs,
                 mod.idx=FLAT.rhs.mod.idx, group=FLAT.group,
                 fixed=FLAT.fixed, start=FLAT.start,
                 label=FLAT.label)

    # change op for intercepts (for convenience only)
    int.idx <- which(FLAT$op == "~" & FLAT$rhs == "")
    if(length(int.idx) > 0L)
        FLAT$op[int.idx] <- "~1"
  
    if(as.data.frame.) FLAT <- as.data.frame(FLAT, stringsAsFactors=FALSE)
  
    attr(FLAT, "modifiers") <- MOD
    attr(FLAT, "constraints") <- CON
  
    FLAT
}

parse.rhs <- function(rhs, op="") {

    # new version YR 15 dec 2011!
    # - no 'equal' field anymore (only labels!)
    # - every modifier is evaluated
    # - unquoted labels are allowed (eg. x1 + x2 + c(v1,v2,v3)*x3)

    getModifier <- function(mod) {
        if(length(mod) == 1L) {
            # three possibilites: 1) numeric, 2) NA, or 3) quoted character
            if( is.numeric(mod) ) 
                return( list(fixed=mod) )
            if( is.na(mod) ) 
                return( list(fixed=as.numeric(NA)) )
            if( is.character(mod) )
                return( list(label=mod) )
        } else if(mod[[1L]] == "start") {
            cof <- unlist(lapply(as.list(mod)[-1], 
                                 eval, envir=NULL, enclos=NULL))
            return( list(start=cof) )
        } else if(mod[[1L]] == "equal") {
            label <- unlist(lapply(as.list(mod)[-1],    
                            eval, envir=NULL, enclos=NULL))
            return( list(label=label) )
        } else if(mod[[1L]] == "label") {
            label <- unlist(lapply(as.list(mod)[-1],
                            eval, envir=NULL, enclos=NULL))
            label[is.na(label)] <- "" # catch 'NA' elements in a label
            return( list(label=label) )
        } else if(mod[[1L]] == "c") {
            # vector: we allow numeric and character only!
            cof <- unlist(lapply(as.list(mod)[-1],    
                                 eval, envir=NULL, enclos=NULL))
            if(is.numeric(cof)) 
                 return( list(fixed=cof) )
            else if(is.character(cof)) {
                 cof[is.na(cof)] <- "" # catch 'NA' elements in a label
                 return( list(label=cof) )
            } else {
                stop("lavaan ERROR: can not parse modifier:", mod, "\n")
            }
        } else {
            # unknown expression
            # as a final attempt, we will evaluate it and coerce it
            # to either a numeric or character (vector)
            cof <- try( eval(mod, envir=NULL, enclos=NULL), silent=TRUE)
            if(is.numeric(cof))
                 return( list(fixed=cof) )
            else if(is.character(cof))
                 return( list(label=cof) )
            else {
                stop("lavaan ERROR: can not parse modifier:", mod, "\n")
            }
        }
    }

    # fill in rhs list
    out <- list()
    repeat {
        if(length(rhs) == 1L) { # last one and only a single element
            out <- c(vector("list", 1L), out)
            NAME <- all.vars(rhs)
            if(length(NAME) > 0L) {
                names(out)[1L] <- NAME
            } else { # intercept or zero?
                if(as.character(rhs) == "1") {
                    names(out)[1L] <- "intercept"
                } else if(as.character(rhs) == "0") {
                    names(out)[1L] <- "zero"
                    out[[1L]]$fixed <- 0
                } else {
                    names(out)[1L] <- "constant"
                    out[[1L]]$fixed <- 0 
                }
            }
            break
        } else if(rhs[[1L]] == "*") { # last one, but with modifier
            out <- c(vector("list", 1L), out)
            NAME <- all.vars(rhs[[3L]])
            if(length(NAME) > 0L) {
                names(out)[1L] <- NAME
            } else { # intercept
                names(out)[1L] <- "intercept"
            }
            i.var <- all.vars(rhs[[2L]], unique=FALSE)
            if(length(i.var) > 0L) {
                # modifier are unquoted labels
                out[[1L]]$label <- i.var
            } else {
                # modifer is something else
                out[[1L]] <- getModifier(rhs[[2L]])
            }
            break
        } else if(rhs[[1L]] == "+") { # not last one!
            i.var <- all.vars(rhs[[3L]], unique=FALSE)
            n.var <- length(i.var)
            out <- c(vector("list", 1L), out)
            if(length(i.var) > 0L) {
                names(out)[1L] <- i.var[n.var]
            } else {
                names(out)[1L] <- "intercept"
            }
            if(n.var > 1L) { 
                # modifier are unquoted labels
                out[[1L]]$label <- i.var[-n.var]
            } else if(length(rhs[[3]]) == 3L) {
                # modifiers!!
                out[[1L]] <- getModifier(rhs[[3L]][[2L]])
            }

            # next element
            rhs <- rhs[[2L]]
        } else {
            stop("lavaan ERROR: I'm confused parsing this line: ", rhs, "\n")
        }
    }

    # if multiple elements, check for duplicated elements and merge if found
    if(length(out) > 1L) {
        rhs.names <- names(out)
        while( !is.na(idx <- which(duplicated(rhs.names))[1L]) ) {
            dup.name <- rhs.names[ idx ]
            orig.idx <- match(dup.name, rhs.names)
            merged <- c( out[[orig.idx]], out[[idx]] )
            if(!is.null(merged)) # be careful, NULL will delete element
                out[[orig.idx]] <- merged
            out <- out[-idx]
            rhs.names <- names(out)
        }
    }

    # if thresholds, check order and reorder if necessary
    #if(op == "|") {
    #    t.names <- names(out)
    #    idx <- match(sort(t.names), t.names)
    #    out <- out[idx]
    #}

    out
}

