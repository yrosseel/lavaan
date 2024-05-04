# LDW 26 Mar 2024: use option settings and store in cache environment

lavaan_cache_env <- new.env(parent = emptyenv())

# set the default options (including unspecified values "default")
lav_options_default <- function() {
  if (exists("opt.default", lavaan_cache_env)) {
    opt <- get("opt.default", lavaan_cache_env)
    return(opt)
  }
  # ---------------- preparation -----------------
  opt.default <- list()
  opt.check <- list()
  elm <- function(
    name = NULL,          # name of option, if length 2 first is sublist name
    dflt = NULL,          # default value
    chr = NULL,           # valid strings (names) and replacement values
    nm = NULL,            # valid numeric interval
    bl = FALSE,           # logical OK?
    oklen = c(1L, 1L),    # lengte > 1 OK, second negative to have a warning
    #                and not an error when length greater then abs(oklen[2])
    num2int = FALSE       # should numerical values be converted to int?
  ) {
    stopifnot(any(length(name) == 1:2))
    stopifnot(is.null(nm) || (length(nm) == 1 && grepl("^[][].*,.*[][]$", nm)))
    stopifnot(length(bl) == 1, is.logical(bl))
    stopifnot(is.null(chr) || is.character(chr))
    stopifnot(length(oklen) == 2, oklen[1] <= abs(oklen[2]))
    stopifnot(length(num2int) == 1, is.logical(num2int))

    # prepare list to store for checking option
    list2store <- list(oklen = oklen)
    if (!is.null(chr)) {
      if (is.null(names(chr))) names(chr) <- chr
      list2store$chr <- chr
    }
    if (!is.null(nm)) {
      first.in <- grepl("^\\[", nm)
      last.in <- grepl("\\]$", nm)
      elems <- as.numeric(strsplit(gsub("[][ ]", "", nm), ",")[[1]])
      if (num2int) {
        elems[elems == -Inf] <- -2e9
        elems[elems == Inf] <- 2e9
        elems <- as.integer(elems)
      }
      list2store$nm <- list(bounds = elems,
                            first.in = first.in, last.in = last.in)
    }
    if (bl) list2store$bl <- TRUE
    if (num2int) list2store$num2int <- TRUE

    # store default and list for checking
    if (length(name) == 1) name <- c("", name)
    if (name[1] != "") {
      if (is.null(opt.default[[name[1]]])) { # make sure sublists exist
        opt.default[[name[1]]] <<- list()
        sublist <- list()
        attr(sublist, "SUB") <- TRUE # indicate as sublist
        opt.check[[name[1]]] <<- sublist
      }
      opt.default[[name[1]]][[name[2]]] <<- dflt
      if (is.null(dflt)) opt.default[[name[1]]][name[2]] <<- list(NULL)
      opt.check[[name[1]]][[name[2]]] <<- list2store
    } else {
      opt.default[[name[2]]] <<- dflt
      if (is.null(dflt)) opt.default[name[2]] <<- list(NULL)
      opt.check[[name[2]]] <<- list2store
    }
    NULL
  }
  elmdup <- function(
    name = NULL, # name(s) of option
    from = NULL  # name(s) of option to duplicatie
  ) {
    if (length(name) == 1) name <- c("", name)
    if (length(from) == 1) from <- c("", from)
    if (from[1] != "") {
      from.default <- opt.default[[from[1]]][[from[2]]]
      from.check <- opt.check[[from[1]]][[from[2]]]
    } else {
      from.default <- opt.default[[from[2]]]
      from.check <- opt.check[[from[2]]]
    }
    if (name[1] != "") {
      if (is.null(opt.default[[name[1]]])) { # make sure sublists exist
        opt.default[[name[1]]] <<- list()
        sublist <- list()
        attr(sublist, "SUB") <- TRUE # indicate as sublist
        opt.check[[name[1]]] <<- sublist
      }
      opt.default[[name[1]]][[name[2]]] <<- from.default
      opt.check[[name[1]]][[name[2]]] <<- from.check
    } else {
      opt.default[[name[2]]] <<- from.default
      opt.check[[name[2]]] <<- from.check
    }
  }
  # ------------------------- store options --------------------------
  elm("model.type", "sem", chr = c(lavaan = "lavaan", cfa = "cfa",
            growth = "growth", sem = "sem", efa = "efa", path = "path",
            unrestricted = "unrestricted"))

  # global
  elm("mimic", "lavaan", chr = c(default = "lavaan", lavaan = "lavaan",
                                 regression = "lm", lisrel = "EQS",
                                 eqs = "EQS", lm = "lm", mplus = "Mplus"
  ))

  # model modifiers
  elm("meanstructure", "default", chr = "default", bl = TRUE)
  elm("int.ov.free", FALSE, bl = TRUE)
  elm("int.lv.free", FALSE, bl = TRUE)
  elm("marker.int.zero", FALSE, bl = TRUE) # fix maker intercepts
  # free lv means
  elm("conditional.x", "default", chr = "default", bl = TRUE)
  elm("fixed.x", "default", chr = "default", bl = TRUE)
  elm("orthogonal", FALSE, bl = TRUE)
  elm("orthogonal.x", FALSE, bl = TRUE)
  elm("orthogonal.y", FALSE, bl = TRUE)
  elm("std.lv", FALSE, bl = TRUE)
  elm("correlation", FALSE, bl = TRUE) # correlation structure
  elm("effect.coding", FALSE, chr = c("",
                                      "loadings", "intercepts",
                                      "mg.lv.efa.variances", "mg.lv.variances",
                                      "mg.lv.means", "mg.lv.intercepts"),
      bl = TRUE, oklen = c(0L, 6L))
  elm("ceq.simple", FALSE, bl = TRUE) # treat simple eq cons special?
  elm("parameterization", "default", c(
    "default", "mml", "delta", "theta"))
  elm("auto.fix.first", FALSE, bl = TRUE)
  elm("auto.fix.single", FALSE, bl = TRUE)
  elm("auto.var", FALSE, bl = TRUE)
  elm("auto.cov.lv.x", FALSE, bl = TRUE)
  elm("auto.cov.y", FALSE, bl = TRUE)
  elm("auto.th", FALSE, bl = TRUE)
  elm("auto.delta", FALSE, bl = TRUE)
  elm("auto.efa", FALSE, bl = TRUE)

  # rotation
  elm("rotation", "geomin", chr = c(crawfer = "cf", crawford.ferguson = "cf",
    crawfordferguson = "cf", cf = "cf",
    varimax = "varimax", quartimax = "quartimax", orthomax = "orthomax",
    oblimin = "oblimin", quartimin = "quartimin", geomin = "geomin",
    entropy = "entropy", mccammon = "mccammon", infomax = "infomax",
    tandem1 = "tandem1", tandem2 = "tandem2", none = "none", promax = "promax",
    oblimax = "oblimax", bentler = "bentler", simplimax = "simplimax",
    target = "target", pst = "pst", cf.quartimax = "cf-quartimax",
    cf.varimax = "cf-varimax", cf.equamax = "cf-equamax",
    cf.parsimax = "cf-parsimax", cf.facparsim = "cf-facparsim",
    bi.quartimin = "biquartimin",
    biquartimin = "biquartimin", bi.geomin = "bigeomin", bigeomin = "bigeomin"
  ))
  elm("rotation.se", "bordered", chr = c("delta", "bordered"))

  # rotation-args sublist
  elm(c("rotation.args", "orthogonal"), FALSE, bl = TRUE)
  elm(c("rotation.args", "row.weights"), "default", chr = c(
    default = "default", kaiser = "kaiser", none = "none",
    cureton.mulaik = "cm", cm = "cm"))
  elm(c("rotation.args", "std.ov"), TRUE, bl = TRUE)
  elm(c("rotation.args", "geomin.epsilon"), 0.001, nm = "]0, 0.01]")
  # was 0.01 < 0.6-10
  elm(c("rotation.args", "orthomax.gamma"), 1, nm = "[0, 1]")
  elm(c("rotation.args", "cf.gamma"), 0, nm = "[0, 1]")
  elm(c("rotation.args", "oblimin.gamma"), 0, nm = "[0, 1000]")
  elm(c("rotation.args", "promax.kappa"), 4, nm = "[0, 1000]")
  elm(c("rotation.args", "target"), matrix(0, 0, 0), oklen = c(0L, 1000L))
  elm(c("rotation.args", "target.mask"), matrix(0, 0, 0), oklen = c(0L, 1000L))
  elm(c("rotation.args", "rstarts"), 30L, nm = "[0, 1000000]")
  elm(c("rotation.args", "algorithm"), "gpa", chr = c("gpa", "pairwise"))
  elm(c("rotation.args", "reflect"), TRUE, bl = TRUE)
  elm(c("rotation.args", "order.lv.by"), "index",
      chr = c("sumofsquares", "index", "none"))
  elm(c("rotation.args", "gpa.tol"), 1e-05, nm = "]0, 0.01]")
  elm(c("rotation.args", "tol"), 1e-08, nm = "]0, 0.01]")
  elm(c("rotation.args", "warn"), FALSE, bl = TRUE)
  elm(c("rotation.args", "verbose"), FALSE, bl = TRUE)
  elm(c("rotation.args", "jac.init.rot"), TRUE, bl = TRUE)
  elm(c("rotation.args", "max.iter"), 10000L, nm = "[0, 1000000]")

  # full data
  elm("std.ov", FALSE, bl = TRUE)
  elm("missing", "default", chr = c(
    default = "default", ml = "ml", direct = "ml",
    ml.x = "ml.x", direct.x = "ml.x", fiml.x = "ml.x", fiml = "ml",
    two.stage = "two.stage", twostage = "two.stage", two.step = "two.stage",
    twostep = "two.stage", robust.two.stage = "robust.two.stage",
    robust.twostage = "robust.two.stage",
    robust.two.step = "robust.two.stage",
    robust.twostep = "robust.two.stage",
    two.stage.robust = "robust.two.stage",
    twostage.robust = "robust.two.stage",
    two.step.robust = "robust.two.stage",
    twostep.robust = "robust.two.stage",
    listwise = "listwise", pairwise = "pairwise",
    available.cases = "available.cases", doubly.robust = "doubly.robust"))
  elm("sampling.weights.normalization", "total", chr = c(
    "total", "group", "none"))
  elm("samplestats", TRUE, bl = TRUE)

  # summary data
  elm("sample.cov.rescale", "default", bl = TRUE)
  elm("sample.cov.robust", FALSE, bl = TRUE)
  elm("sample.icov", TRUE, bl = TRUE)
  elm("ridge", FALSE, bl = TRUE)
  elm("ridge.constant", "default", chr = "default", nm = "[0, Inf[")

  # multiple groups !!! group.label and group.partial capitals OK !!!
  elm("group.label", NULL, oklen = c(0L, 100L)) # no checks
  elm("group.equal", "", chr =
        c("", "none", "loadings", "intercepts", "means", "composite.loadings",
          "regressions", "residuals", "residual.covariances", "thresholds",
          "lv.variances", "lv.covariances"), oklen = c(0L, 100L))
  elm("group.partial", "", oklen = c(0L, 100L)) # no checks
  elm("group.w.free", FALSE, bl = TRUE)

  # clusters
  elm("level.label", NULL, oklen = c(0L, 100L)) # no checks

  # estimation
  elm("estimator", "default", chr = c(
    none = "none", default = "default", wlsmv = "wlsmv", ml = "ml", mlr = "mlr",
    mlf = "mlf", mlm = "mlm", mlmv = "mlmv", mlmvs = "mlmvs", gls = "gls",
    wls = "wls", wlsm = "wlsm", uls = "uls", ulsm = "ulsm", ulsmv = "ulsmv",
    pml = "pml", dls = "dls", ntrls = "ntrls", catml = "catml",
    dwls = "dwls", wlsmvs = "wlsmvs", ulsmvs = "ulsmvs", fml = "fml",
    umn = "fml", reml = "reml", mml = "mml", fabin = "fabin2",
    fabin2 = "fabin2", fabin3 = "fabin3", mgm = "mgm", guttman = "mgm",
    gutman = "mgm", gutmann = "mgm", guttman1952 = "mgm",
    js = "js", jsa = "jsa", james.stein = "js",
    james.stein.aggregated = "jsa", bentler = "bentler1982",
    bentler1982 = "bentler1982", miiv = "miiv", iv = "miiv",
    miiv.2sls = "miiv"
  ))
  elmdup("estimator.orig", "estimator")

  elm("estimator.args", list(), oklen = c(0L, 100L))
  elm("likelihood", "default", chr = c("default", "normal", "wishart"))
  elm("link", "default", chr = c("default", "logit", "probit"))
  elm("representation", "default", chr = c(
    default = "LISREL", lisrel = "LISREL", ram = "RAM"))
  elm("do.fit", TRUE, bl = TRUE)
  elm("bounds", "none", chr = c(
    "none", "default", "standard", "user", "wide", "wide.zerovar", "pos.var",
    "pos.ov.var", "pos.lv.var"))      # new in 0.6-6
  elm("rstarts", 0L, nm = "[0, 1000]", num2int = TRUE) # new in 0.6-18

  # inference
  elm("se", "default", chr = c(
    default = "default", none = "none", standard = "standard",
    robust.huber.white = "robust.huber.white", robust = "robust",
    robust.cluster = "robust.cluster",
    robust.cluster.sem = "robust.cluster.sem",
    sandwich = "robust.huber.white", robust.sem = "robust.sem",
    two.stage = "two.stage", robust.two.stage = "robust.two.stage",
    bootstrap = "bootstrap", boot = "bootstrap", first.order = "first.order",
    robust.mlm = "robust.sem", robust.mlr = "robust.huber.white",
    observed = "observed", expected = "expected"),
    oklen = c(1L, -1L)
  )
  elm("test", "default", oklen = c(1L, 100L))
                  # checks for 'test' are in lav_test_rename !!!

  # information (se + test)
  elm("information", c("default", "default"), chr = c(
    "default", "expected", "observed", "first.order"), oklen = c(1L, 2L))
  elm("h1.information", c("structured", "structured"), chr = c(
    "structured", "unstructured"), oklen =  c(1L, 2L))
  elm("observed.information", c("hessian", "default"), chr = c(
    "default", "hessian", "h1"), oklen =  c(1L, 2L))

  # information se only
  elm("information.meat", "default",
      chr = c(default = "first.order", first.order = "first.order"))
  elm("h1.information.meat", "default", chr = c(
    "default", "structured", "unstructured"))

  # information for 'Omega' (yuan-benter test only)
  elm("omega.information", "default", chr = c(
    "default", "expected", "observed"
  ))
  elm("omega.h1.information", "default", chr = c(
    "default", "structured", "unstructured"
  ))
  elm("omega.information.meat", "default", chr = c(
    default = "first.order", first.order = "first.order"
  ))
  elm("omega.h1.information.meat", "default", chr = c(
    "default", "structured", "unstructured"
  ))

  # test statistic for scaling
  elm("scaled.test", "standard", oklen = c(1L, 100L))

  # old approach trace.UGamma2
  elm("ug2.old.approach", FALSE, bl = TRUE)

  # bootstrap
  elm("bootstrap", 1000L, nm = "[1, Inf[", num2int = TRUE)

  # gamma
  elm("gamma.n.minus.one", FALSE, bl = TRUE)
  elm("gamma.unbiased", FALSE, bl = TRUE)

  # optimization
  elm("control", list(), oklen = c(0L, 100L))
  elm("optim.method", "default", chr = c(
    "nlminb", "gn", "default", "noniter", "none", "em"
  )) # gn for DLS, nlminb rest
  elm("optim.attempts", 4L, nm = "[1, 4]")
  elm("optim.force.converged", FALSE, bl = TRUE)
  elm("optim.gradient", "analytic", chr = c(
    analytic = "analytic", analytical = "analytic",
    numeric = "numerical", numerical = "numerical"
  ))
  elm("optim.init_nelder_mead", FALSE, bl = TRUE)
  elm("optim.var.transform", "none", chr = c(
    "none", "sqrt"
  ))
  elm("optim.parscale", "none", chr = c(
    none = "none", st = "standardized", stand = "standardized",
    standardize = "standardized", standardized = "standardized"
  ))
  elm("optim.partrace", FALSE, bl = TRUE)
  elm("optim.dx.tol", 1e-03, nm = "]0, 0.01]") # not too strict
  elm("optim.bounds", list(), oklen = c(0L, 100L))
  elm("em.iter.max", 10000L, nm = "[100, 1e8]", num2int = TRUE)
  elm("em.fx.tol", 1e-08, nm = "]0, 0.01]")
  elm("em.dx.tol", 1e-04, nm = "]0, 0.01]")
  elm("em.zerovar.offset", 0.0001, nm = "]0, 0.01]")
  elm("em.h1.iter.max", 500L, nm = "[10, 1e7]", num2int = TRUE)
  elm("em.h1.tol", 1e-05, nm = "]0, 0.01]") # was 1e-06 < 0.6-9
  elm("em.h1.warn", TRUE, bl = TRUE)
  elm("optim.gn.iter.max", 200L, nm = "[100, 1e8]", num2int = TRUE)
  elm("optim.gn.stephalf.max", 10L, nm = "[1, 1e8]", num2int = TRUE)
  elm("optim.gn.tol.x", 1e-05, nm = "]0, 0.01]")

  # numerical integration
  elm("integration.ngh", 21L, nm = "[1, 1000]", num2int = TRUE)

  # parallel
  elm("parallel", "no", chr = c(
    "no", "multicore", "snow"
  ))
  maxcpu <- parallel::detectCores() - 1L
  elm("ncpus", maxcpu, nm = paste0("[1,", maxcpu, "]"))
  elm("cl", NULL, oklen = c(0L, 1L))
  elm("iseed", NULL, oklen = c(0L, 1L))

  # categorical
  elm("zero.add", c(0.5, 0.0), chr = "default",
      nm = "[0, 1]", oklen = c(1L, -2L))
  elm("zero.keep.margins", "default", chr = "default", bl = TRUE)
  elm("zero.cell.warn", FALSE, bl = TRUE) # since 0.6-1
  elm("cat.wls.w", TRUE, bl = TRUE) # since 0.6-18

  # starting values (char values checked in lav_options_set())
  elm("start", "default", oklen = c(1L, 1000L))

  # sanity checks
  elm("check.start", TRUE, bl = TRUE)
  elm("check.post", TRUE, bl = TRUE)
  elm("check.gradient", TRUE, bl = TRUE)
  elm("check.vcov", TRUE, bl = TRUE)
  elm("check.lv.names", TRUE, bl = TRUE)
  elm("check.lv.interaction", TRUE, bl = TRUE)

  # more models/info
  elm("h1", TRUE, bl = TRUE)
  elm("baseline", TRUE, bl = TRUE)
  elm("baseline.conditional.x.free.slopes", TRUE, bl = TRUE)
  elm("implied", TRUE, bl = TRUE)
  elm("loglik", TRUE, bl = TRUE)

  # storage of information
  elm("store.vcov", "default", chr = "default", bl = TRUE)

  # internal
  elm("parser", "new", chr = c(old = "old", orig = "old", new = "new",
                               classic = "old"))

  # verbosity
  elm("verbose", FALSE, bl = TRUE)
  elm("warn", TRUE, bl = TRUE)
  elm("debug", FALSE, bl = TRUE)

  # categorical
  elm("categorical", "default", chr = "default", bl = TRUE)

  # ------------- store info in lavaan environment ---------------
  assign("opt.default", opt.default, lavaan_cache_env)
  assign("opt.check", opt.check, lavaan_cache_env)

  # return defaults
  return(opt.default)
}

# public function
lavOptions <- function(x = NULL, default = NULL, mimic = "lavaan") { # nolint
  lavoptions <- lav_options_default()

  # selection only
  if (!is.null(x)) {
    if (is.character(x)) {
      # lower case only
      x <- tolower(x)

      # check if x is in names(lavoptions)
      not.ok <- which(!x %in% names(lavoptions))
      if (length(not.ok) > 0L) {
        lav_msg_warn(gettextf(
          "option(s) %s not available", lav_msg_view(x[not.ok]))
        )
        x <- x[-not.ok]
      }

      # return requested option(s)
      if (length(x) == 0L) {
        return(default)
      } else {
        lavoptions[x]
      }
    } else {
      lav_msg_stop(gettext("`x' must be a character string"))
    }
  } else {
    lavoptions
  }
}
