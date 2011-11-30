# class definitions
# 
# initial version: YR 25/03/2009
# added ModelSyntax: YR 02/08/2010
# deleted ModelSyntax: YR 01/11/2010 (using flattened model syntax now)

setClass("SampleStats",            # sample moments
    representation(
        cov="list",                # observed var/cov matrix (per group)
        mean="list",               # observed mean vector (per group)

        nobs="list",               # effective number of obs (per group)
        norig="list",              # original number of obs (per group)
        ntotal="integer",          # total number of obs (all groups)
        ngroups="integer",         # number of groups
        ov.names="list",           # variable names (per group)

        missing="list"             # missingness information
    )
)

setClass("SampleStatsExtra",       # additional sample statistics 
    representation(
        icov="list",
        cov.log.det="list",
        cov.vecs="list",
        WLS.V="list"
    )
)

setClass("Model",          # MATRIX representation of the sem model
    representation(
        GLIST="list",              # list of all model matrices (for all groups)
        dimNames="list",           # dim names for the model matrices
        isSymmetric="logical",     # model matrix symmetric?
        mmSize="integer",          # model matrix size (unique only)

        representation="character",  # stub, until we define more classes
        meanstructure="logical",

        ngroups="integer",
        nmat="integer",
        nvar="integer",

        nx.free="integer",
        nx.unco="integer",
        nx.user="integer",

        m.free.idx="list",
        x.free.idx="list",
        m.unco.idx="list",
        x.unco.idx="list",
        m.user.idx="list",
        x.user.idx="list",
        x.def.idx="integer",
        x.ceq.idx="integer",
        x.cin.idx="integer",
        x.free.var.idx="integer",

        eq.constraints="logical",
        eq.constraints.K="matrix",

        def.function="function",
        ceq.function="function",
        ceq.jacobian="function",
        cin.function="function",
        cin.jacobian="function",
        con.jac="matrix",

        fixed.x="logical"

    )
)

setClass("Fit",
    representation(
        npar="integer",            # number of free parameters
        #ndat="integer",
        #df="integer",
        x="numeric",               # x
        start="numeric",           # starting values
        est="numeric",             # estimated values
        se="numeric",              # standard errors
        fx="numeric",
        fx.group="numeric",
        iterations="integer",      # number of iterations
        converged="logical",
        Sigma.hat="list",
        Mu.hat="list",
        test="list"
    )
)

setClass("lavaan",
    representation(
        call    = "call",            # matched call
        timing  = "list",            # timing information
        Options = "list",            # lavaanOptions
        User    = "list",            # user specified model
        Data    = "list",            # full data 
        Sample  = "SampleStats",     # sample statistics
        Extra   = "SampleStatsExtra",# extra sample statistics
        Model   = "Model",           # internal matrix representation
        Fit     = "Fit"              # optimization info
    ) 
)






                      
