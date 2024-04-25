# class definitions
#
# initial version: YR 25/03/2009
# added ModelSyntax: YR 02/08/2010
# deleted ModelSyntax: YR 01/11/2010 (using flattened model syntax now)
# ldw 20/11/2023: replace 'representation()' by 'slots='

setClass("lavData",
  slots = c(
    data.type = "character", # "full", "moment" or "none"
    group = "character", # group variable
    ngroups = "integer", # number of groups
    group.label = "character", # group labels
    block.label = "character", # block labels
    cluster = "character", # cluster variable(s)
    nlevels = "integer", # number of levels
    level.label = "character", # level labels
    std.ov = "logical", # standardize observed variables?
    nobs = "list", # effective number of observations
    norig = "list", # original number of observations
    ov.names = "list", # variable names (per group)
    ov.names.x = "list", # exo variable names (per group)
    ov.names.l = "list", # names per level
    # ov.types          = "list",           # variable types (per group)
    # ov.idx            = "list",           # column indices (all observed variables)
    ordered = "character", # ordered variables
    weights = "list", # sampling weights (per group)
    sampling.weights = "character", # sampling weights variable
    ov = "list", # variable table
    case.idx = "list", # case indices per group
    missing = "character", # "listwise" or not?
    Mp = "list", # if not complete, missing patterns
    # we need this here, to get nobs right!
    Rp = "list", # response patterns (categorical only)
    Lp = "list", # level patterns
    eXo = "list", # local copy exo only
    X = "list" # local copy
  )
)


setClass("lavSampleStats", # sample moments
  slots = c(
    var                = "list", # observed variances (per group)
    cov                = "list", # observed var/cov matrix (per group)
    mean               = "list", # observed mean vector (per group)
    th                 = "list", # thresholds for non-numeric var (per group)
    th.idx             = "list", # th index (0 for numeric)
    th.names           = "list", # threshold names

    res.cov            = "list", # residual var/cov matrix (if conditional.x)
    res.var            = "list", # residual variances
    res.th             = "list", # residual thresholds
    res.th.nox         = "list", # residual thresholds ignoring x
    res.slopes         = "list", # slopes exo (if conditional.x)
    res.int            = "list", # intercepts (if conditional.x)

    mean.x             = "list", # mean exo
    cov.x              = "list", # variance/covariance exo
    bifreq             = "list", # bivariate frequency tables
    group.w            = "list", # group weight

    nobs               = "list", # effective number of obs (per group)
    ntotal             = "numeric", # total number of obs (all groups)
    ngroups            = "integer", # number of groups
    x.idx              = "list", # x.idx if fixed.x = TRUE

    icov               = "list", # inverse of observed cov (per group)
    cov.log.det        = "list", # log det of observed cov (per group)
    res.icov           = "list",
    res.cov.log.det    = "list",
    # ridge.constant   = "numeric",        # ridge constant (per group)
    # ridge.constant.x = "numeric",        # ridge constant (per group) for eXo
    ridge              = "numeric",
    WLS.obs            = "list", # all relevant observed stats in a vector
    WLS.V              = "list", # weight matrix for GLS/WLS
    WLS.VD             = "list", # diagonal of weight matrix only
    NACOV              = "list", # N times the asymptotic covariance matrix
    NACOV.user         = "logical", # user-specified NACOV?

    missing.flag       = "logical", # missing patterns?
    missing            = "list", # missingness information
    missing.h1         = "list", # h1 model

    YLp                = "list", # cluster/level information

    zero.cell.tables   = "list" # bivariate tables with empty cells
  )
)

setClass("lavModel", # MATRIX representation of the sem model
  slots = c(
    GLIST              = "list", # list of all model matrices (for all groups)
    dimNames           = "list", # dim names for the model matrices
    isSymmetric        = "logical", # model matrix symmetric?
    mmSize             = "integer", # model matrix size (unique only)

    representation     = "character", # stub, until we define more classes
    modprop            = "list", # model properties
    meanstructure      = "logical",
    correlation        = "logical",
    categorical        = "logical",
    multilevel         = "logical",
    group.w.free       = "logical",
    link               = "character",
    nblocks            = "integer",
    ngroups            = "integer", # only for rsem!! (which uses rsem:::computeDelta)
    nefa               = "integer",
    nmat               = "integer",
    nvar               = "integer",
    num.idx            = "list",
    th.idx             = "list",
    nx.free            = "integer",
    nx.unco            = "integer",
    nx.user            = "integer",
    m.free.idx         = "list",
    x.free.idx         = "list",
    # m.unco.idx        = "list",           # always the same as m.free.idx
    x.unco.idx         = "list",
    m.user.idx         = "list",
    x.user.idx         = "list",
    x.def.idx          = "integer",
    x.ceq.idx          = "integer",
    x.cin.idx          = "integer",
    x.free.var.idx     = "integer",
    ceq.simple.only    = "logical",
    ceq.simple.K       = "matrix",
    eq.constraints     = "logical",
    eq.constraints.K   = "matrix",
    eq.constraints.k0  = "numeric",
    def.function       = "function",
    ceq.function       = "function",
    ceq.jacobian       = "function",
    ceq.JAC            = "matrix",
    ceq.rhs            = "numeric",
    ceq.linear.idx     = "integer",
    ceq.nonlinear.idx  = "integer",
    cin.function       = "function",
    cin.jacobian       = "function",
    cin.JAC            = "matrix",
    cin.rhs            = "numeric",
    cin.linear.idx     = "integer",
    cin.nonlinear.idx  = "integer",
    ceq.efa.JAC        = "matrix",
    con.jac            = "matrix",
    con.lambda         = "numeric",
    nexo               = "integer",
    conditional.x      = "logical",
    fixed.x            = "logical",
    parameterization   = "character",
    ov.x.dummy.ov.idx  = "list",
    ov.x.dummy.lv.idx  = "list",
    ov.y.dummy.ov.idx  = "list",
    ov.y.dummy.lv.idx  = "list",
    ov.efa.idx         = "list",
    lv.efa.idx         = "list",
    rv.ov              = "list",
    rv.lv              = "list",
    H                  = "list",
    lv.order           = "list",
    estimator          = "character",
    estimator.args     = "list"
  )
)

setClass("Fit",
  slots = c(
    npar               = "integer", # number of free parameters
    # ndat              = "integer",
    # df                = "integer",
    x                  = "numeric", # x
    partrace           = "matrix", # parameter trace
    start              = "numeric", # starting values (only for other packages)
    est                = "numeric", # estimated values (only for other packages)
    se                 = "numeric", # standard errors
    fx                 = "numeric",
    fx.group           = "numeric",
    logl               = "numeric",
    logl.group         = "numeric",
    iterations         = "integer", # number of iterations
    converged          = "logical",
    control            = "list",
    Sigma.hat          = "list",
    Mu.hat             = "list",
    TH                 = "list",
    test               = "list"
  )
)

setClass("lavaan",
  slots = c(
    version            = "character", # lavaan version
    call               = "call", # matched call
    timing             = "list", # timing information
    Options            = "list", # lavOptions
    ParTable           = "list", # parameter table user-specified model
    pta                = "list", # parameter table attributes
    Data               = "lavData", # full data
    SampleStats        = "lavSampleStats", # sample statistics
    Model              = "lavModel", # internal matrix representation
    Cache              = "list", # housekeeping stuff
    Fit                = "Fit", # fitted results
    boot               = "list", # bootstrap results
    optim              = "list", # optimizer results
    loglik             = "list", # loglik values and info
    implied            = "list", # model implied moments
    vcov               = "list", # vcov
    test               = "list", # test
    h1                 = "list", # unrestricted model results
    baseline           = "list", # baseline model results
    internal           = "list", # optional slot, for internal use
    external           = "list" # optional slot, for add-on packages
  )
)

setClass("lavaanList",
  slots = c(
    version            = "character", # lavaan version
    call               = "call", # matched call
    Options            = "list", # lavOptions
    ParTable           = "list",
    pta                = "list",
    Data               = "lavData", # from first dataset (ngroups!)
    Model              = "lavModel", # based on first dataset
    meta               = "list",
    timingList         = "list",
    ParTableList       = "list",
    DataList           = "list",
    SampleStatsList    = "list",
    CacheList          = "list",
    vcovList           = "list",
    testList           = "list",
    optimList          = "list",
    impliedList        = "list",
    h1List             = "list",
    loglikList         = "list",
    baselineList       = "list",
    internalList       = "list",
    funList            = "list",
    external           = "list" # optional slot, for add-on packages
  )
)
