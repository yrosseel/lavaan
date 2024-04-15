# initial version YR 02/08/2010

# YR 28 Jan 2017: add lavOptions(), lav_options_default()
# LDW 26 Mar 2024: use option settings and store in cache environment

lavaan_cache_env <- new.env(parent = emptyenv()) # perhaps better to put
# this in lav_utils.R ?

lav_options_settings <- function() {
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
      opt.check[[name[1]]][[name[2]]] <<- list2store
    } else {
      opt.default[[name[2]]] <<- dflt
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
  elm("model.type", "sem", chr = c("sem", "cfa", "growth", "unrestricted"))
  
  # global
  elm("mimic", "lavaan", chr = c(default = "lavaan",
                                 lavaan = "lavaan", regression = "lm", lisrel = "EQS",
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
                                      "loadings", "intercepts", "mg.lv.efa.variances", "mg.lv.variances",
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
  elm("rotation", "geomin", chr = c(
    crawfer = "cf", crawford.ferguson = "cf", crawfordferguson = "cf", cf = "cf",
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
  
  # summary data
  elm("sample.cov.rescale", "default", bl = TRUE)
  elm("sample.cov.robust", FALSE, bl = TRUE)
  elm("sample.icov", TRUE, bl = TRUE)
  elm("ridge", FALSE, bl = TRUE)
  elm("ridge.constant", "default", chr = "default", nm = "[0, Inf[")
  
  # multiple groups !!! group.label and group.partial capitals OK !!!
  elm("group.label", NULL, oklen = c(0L, 100L)) # no checks
  elm("group.equal", "", chr = c("",
                                 "none", "loadings", "intercepts", "means", "composite.loadings",
                                 "regressions", "residuals", "residual.covariances", "thresholds",
                                 "lv.variances", "lv.covariances"), oklen = c(0L, 100L))
  elm("group.partial", "", oklen = c(0L, 100L)) # no checks
  elm("group.w.free", FALSE, bl = TRUE)
  
  # clusters
  elm("level.label", NULL, oklen = c(0L, 100L)) # no checks
  
  # estimation
  elm("estimator", "default", chr = c(none = "none",
                                      default = "default", wlsmv = "wlsmv", ml = "ml", mlr = "mlr",
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
  elm("test", "default", oklen = c(1L, 100L)) # checks for 'test' are in lav_test_rename !!!
  
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
  elm("cl", NULL, oklen = c(0L, 1L))     # !!! don't know what this is used for
  elm("iseed", NULL, oklen = c(0L, 1L))  # !!! don't know what this is used for
  
  # zero values
  elm("zero.add", c(0.5, 0.0), chr = "default",
      nm = "[0, 1]", oklen = c(1L, -2L))
  elm("zero.keep.margins", "default", chr = "default", bl = TRUE)
  elm("zero.cell.warn", FALSE, bl = TRUE) # since 0.6-1
  
  # starting values
  elm("start", "default", chr = c(
    "default", "simple", "est"
  ))
  
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
  elm("parser", "old", chr = c(old = "old", orig = "old", new = "new",
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
}


# public function
lavOptions <- function(x = NULL, default = NULL, mimic = "lavaan") {
  lavoptions <- lav_options_default(mimic = mimic)
  
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

# set the default options (including unspecified values "default")
lav_options_default <- function(mimic = "lavaan") {
  if (!exists("opt.default", lavaan_cache_env)) lav_options_settings()
  opt <- get("opt.default", lavaan_cache_env)
  opt
}

# help functions ####
# help function to determine estimator 'group'
lav_options_estimatorgroup <- function(estimator) {
  goal <- switch(estimator,
                 ml = ,
                 mlf = ,
                 mlm = ,
                 mlmv = ,
                 mlmvs = ,
                 mlr = "ML",
                 catml = "catML",
                 dwls = ,
                 wlsm = ,
                 wlsmv = ,
                 wlsmvs = "DWLS",
                 uls = ,
                 ulsm = ,
                 ulsmv = ,
                 ulsmvs = "ULS",
                 none = "none",
                 toupper(estimator)
  )
  goal
}
lav_options_checkinterval <- function(x, nm, num2int) {
  if (num2int) x <- as.integer(x)
  oks <- vapply(x, function(x1) {
    (x1 > nm$bounds[1] || (x1 == nm$bounds[1] && nm$first.in)) &&
      (x1 < nm$bounds[2] || (x1 == nm$bounds[2] && nm$last.in))
  }, TRUE)
  all(oks)
}
lav_options_checkvalues <- function(optname, optvalue, chr) {
  optvalid <- names(chr)
  if (is.null(optvalid)) optvalid <- chr
  if (any(optvalid == "empty.string")) {
    optvalid[optvalid == "empty.string"] <- ""
  }
  optvals <- gsub("[_-]", ".", tolower(optvalue))
  optvalsok <- match(optvals, optvalid)
  if (any(is.na(optvalsok))) {
    lav_msg_notallowed_multi(optname, optvalue[is.na(optvalsok)])
  }
  as.vector(chr[optvalsok])
}
lav_options_check <- function(opts, opt.check, subname) {
  opt.names <- names(opts)
  hiddens <- startsWith(opt.names, ".")
  if (any(hiddens)) { # remove hidden options temporarily
    opts.hidden <- opts[hiddens]
    opts <- opts[!hiddens]
    opt.names <- opt.names[!hiddens]
  }
  check.names <- names(opt.check)
  match.opt <- match(opt.names, check.names)
  if (any(is.na(match.opt))) {
    lav_msg_stop(gettextf(
      "Some option(s) unknown: %s !",
      lav_msg_view(opt.names[is.na(match.opt)], log.sep = "none")))
  }
  for (j in seq_along(opts)) {
    opt.name <- opt.names[j]
    opt.value <- opts[[j]]
    opt.check1 <- opt.check[[match.opt[j]]]
    if (!is.null(attr(opt.check1, "SUB"))) {
      opts[[j]] <- lav_options_check(opt.value, opt.check1,
                                     paste0(opt.name, "$"))
      next
    }
    # check length of option value
    if (length(opt.value) < opt.check1$oklen[1]) {
      lav_msg_stop(gettextf(
        "Length of option '%1$s' value must be at least %2$s.",
        paste0(subname, opt.name), opt.check1$oklen[1]))
    }
    if (length(opt.value) > abs(opt.check1$oklen[2])) {
      if (opt.check1$oklen[2] > 0L) {
        lav_msg_stop(gettextf(
          "Length of option '%1$s' value must be maximum %2$s.",
          paste0(subname, opt.name), opt.check1$oklen[2]))
      } else {
        lav_msg_warn(gettextf(
          "Length of option '%1$s' value should be maximum %2$s;",
          paste0(subname, opt.name), -opt.check1$oklen[2]),
          gettextf("Only first %s elements used.", -opt.check1$oklen[2]))
      }
    }
    if (is.null(opt.check1$bl)) opt.check1$bl <- FALSE
    if (!is.null(opt.check1$chr) || !is.null(opt.check1$nm)          ||
        opt.check1$bl) {
      if (!opt.check1$bl || !is.logical(opt.value)) {
        if (!is.null(opt.check1$nm) && is.numeric(opt.value)) {
          num2int <- FALSE
          if (!is.null(opt.check1$num2int)) num2int <- opt.check1$num2int
          if (!lav_options_checkinterval(opt.value, opt.check1$nm, num2int)) {
            lav_msg_stop(gettextf(
              "Value(s) of option %1$s out of range (%2$s)!",
              paste0(subname, opt.name),
              paste0(opt.check1$nm$bounds[1],
                     if (opt.check1$nm$first.in) " <= " else " < ",
                     "x",
                     if (opt.check1$nm$last.in) " <= " else " < ",
                     opt.check1$nm$bounds[2])))
          }
        }
        if (!is.null(opt.check1$chr) && is.character(opt.value)) {
          opt.value <- lav_options_checkvalues(opt.name, opt.value, opt.check1$chr)
          opts[[j]] <- opt.value
        }
      }
    }
  }
  if (any(hiddens)) { # add hidden options
    opts <- modifyList(opts, opts.hidden)
  }
  opts
}
# this function collects and checks the user-provided options/arguments,
# and fills in the "default" values, or changes them in an attempt to
# produce a consistent set of values...
#
# returns a list with the named options
lav_options_set <- function(opt = NULL) {
  # check the presence of necessary hidden options ####
  if (is.null(opt$.categorical) || is.null(opt$.multilevel) || 
      is.null(opt$.clustered)) lav_msg_fixme(
        ".categorical, .multilevel and .clustered must be present")
  
  # get opt.default and opt.check ####
  if (!exists("opt.check", lavaan_cache_env)) lav_options_settings()
  opt.check <- get("opt.check", lavaan_cache_env)
  
  if (opt$debug) {
    cat("lavaan DEBUG: lavaanOptions IN\n")
    str(opt)
    opt$optim.partrace <- TRUE
  }
  
  # check options with definitions ####
  opt <- lav_options_check(opt, opt.check, "")
  
  # first of all: set estimator ####
  if (opt$estimator == "default") {
    if (opt$.categorical) {
      opt$estimator <- "wlsmv"
    } else {
      opt$estimator <- "ml"
    }
  }
  
  # defaults for opt$sample.cov.rescale
  if (opt$sample.cov.rescale == "default") {
    opt$sample.cov.rescale <- switch(
      opt$estimator,
      dls = TRUE,
      fabin2 = ,
      fabin3 = ,
      mgm = ,
      js = ,
      jsa = ,
      bentler1982 = TRUE,
      miiv = TRUE,
      "default"
    )
  }
  
  # option defaults specific for mimic=...
  opt <- lav_options_mimic(opt)
  
  # store opt$estimator as estimator.orig in upper case
  opt$estimator.orig <- toupper(opt$estimator)
  
  # rename names of test statistics if needed, check for invalid values ####
  opt$test <- lav_test_rename(opt$test, check = TRUE)
  
  # same for scaled.test
  opt$scaled.test <- lav_test_rename(opt$scaled.test, check = TRUE)
  
  # rename names of se values, check illegal combinations se/estimator ####
  # pass-through function: may change value of information
  # for backwards compatibility (eg if se = "expected")
  opt <- lav_options_check_se(opt)
  
  # do.fit FALSE implies se="none" and test="none" (unless not default) ####
  if (!opt$do.fit) {
    if (opt$se == "default") opt$se <- "none"
    if (opt$test[1] == "default") opt$test <- "none"
  }
  
  # marker.int.fixed ####
  if (opt$marker.int.zero) {
    opt$meanstructure <- TRUE
    opt$int.ov.free <- TRUE
    if ((is.logical(opt$effect.coding) && opt$effect.coding) ||
        (is.character(opt$effect.coding) && nchar(opt$effect.coding) > 0L)) {
      lav_msg_stop(gettext(
        "effect coding cannot be combined with marker.int.zero = TRUE option"))
    }
    if (opt$std.lv) {
      lav_msg_stop(gettext(
        "std.lv = TRUE cannot be combined with marker.int.zero = TRUE"))
    }
  }
  
  # group.equal and group.partial ####
  if (length(opt$group.equal) > 0L && opt$group.equal[1] == "none") {
    opt$group.equal <- character(0)
  } else if (is.null(opt$group.equal) || all(nchar(opt$group.equal) == 0L)) {
    opt$group.equal <- character(0)
  }
  
  if (is.null(opt$group.partial) || all(nchar(opt$group.partial) == 0L)) {
    opt$group.partial <- character(0)
  } else if (length(opt$group.partial) == 0) {
    # nothing to do
  } else {
    # strip white space
    opt$group.partial <- gsub("[[:space:]]+", "", opt$group.partial)
  }
  
  # if categorical, and group.equal contains "intercepts", also add
  # thresholds (and vice versa)
  if (opt$.categorical && any("intercepts" == opt$group.equal)) {
    opt$group.equal <- unique(c(opt$group.equal, "thresholds"))
  }
  if (opt$.categorical && any("thresholds" == opt$group.equal)) {
    opt$group.equal <- unique(c(opt$group.equal, "intercepts"))
  }
  
  # clustered ####
  # brute-force override (for now)
  if (opt$.clustered && !opt$.multilevel) {
    opt$meanstructure <- TRUE
    
    if (opt$estimator == "mlr") {
      opt$estimator <- "ml"
      opt$test <- "yuan.bentler.mplus"
      opt$se <- "robust.cluster"
    } else if (opt$estimator == "mlm") {
      opt$estimator <- "ml"
      opt$test <- "satorra.bentler"
      opt$se <- "robust.cluster.sem"
    } else if (opt$.categorical) {
      opt$test <- "satorra.bentler"
      opt$se <- "robust.cluster.sem"
    }
    
    # test ####
    if (length(opt$test) == 1L && opt$test == "default") {
      opt$test <- "yuan.bentler.mplus"
    } else if (all(opt$test %in% c(
      "none", "standard",
      "satorra.bentler",
      "yuan.bentler", "yuan.bentler.mplus"
    ))) {
      # nothing to do
    } else if (opt$se == "robust") {
      opt$test <- "yuan.bentler.mplus"
    } else {
      lav_msg_stop(
        gettextf("`test' argument must one of %s in the clustered case",
                 lav_msg_view(c("none", "yuan.bentler", "yuan.bentler.mplus",
                                "satorra.bentler"), log.sep = "or")))
    }
    
    # se ####
    if (opt$se == "default") {
      opt$se <- "robust.cluster"
    } else if (any(opt$se == c("none", "robust.cluster",
                               "robust.cluster.sem"))) {
      # nothing to do
    } else if (opt$se == "robust") {
      opt$se <- "robust.cluster"
    }
    
    # information ####
    if (opt$information[1] == "default") {
      if (opt$se == "robust.cluster" && opt$estimator == "ml") {
        opt$information[1] <- "observed"
      } else {
        opt$information[1] <- "expected"
      }
    }
    if (length(opt$information) > 1L && opt$information[2] == "default") {
      if (opt$se == "robust.cluster") {
        opt$information[2] <- "observed"
      } else {
        opt$information[2] <- "expected"
      }
    }
  }
  
  # multilevel ####
  # brute-force override (for now)
  if (opt$.multilevel) {
    opt$meanstructure <- TRUE
    
    # test
    if (length(opt$test) == 1L && opt$test == "default") {
      # ok, will be set later
    } else if (all(opt$test %in% c("none", "standard", "yuan.bentler"))) {
      # nothing to do
    } else {
      lav_msg_stop(gettextf(
        "`test' argument must one of %s in the multilevel case",
        lav_msg_view(c("none", "standard", "yuan.bentler"), log.sep = "or")))
    }
    
    # se
    if (opt$se == "default") {
      # ok, will be set later
    } else if (any(opt$se == c(
      "none", "standard", "robust.huber.white", "sandwich"))) {
      # nothing to do
    } else if (opt$se == "robust") {
      opt$se <- "robust.huber.white"
    } else {
      lav_msg_stop(gettextf(
        "`se' argument must one of %s  in the multilevel case",
        lav_msg_view(c("none", "standard", "robust.huber.white"),
                     log.sep = "or")))
    }
    
    # information
    if (opt$information[1] == "default") {
      opt$information[1] <- "observed"
    }
    if (length(opt$information) > 1L && opt$information[2] == "default") {
      opt$information[2] <- "observed"
    }
  }
  
  # missing ####
  if (opt$missing == "default") {
    opt$missing <- "listwise"
  } else if (opt$missing == "ml") {
    if (opt$.categorical) {
      lav_msg_stop(gettextf(
        "missing = %s not available in the categorical setting",
        dQuote(opt$missing)))
    }
    if (any(opt$estimator == c(
      "mlm", "mlmv", "gls", "wls", "wlsm", "wlsmv",
      "uls", "ulsm", "ulsmv", "pml", "dls"
    ))) {
      lav_msg_stop(gettextf(
        "missing=%1$s is not allowed for estimator %2$s",
        dQuote(opt$missing), dQuote(lav_options_estimatorgroup(opt$estimator))))
    }
  } else if (opt$missing == "ml.x") {
    if (opt$.categorical) {
      lav_msg_stop(gettextf(
        "missing = %s not available in the categorical setting",
        dQuote(opt$missing)))
    }
    if (any(opt$estimator == c(
      "mlm", "mlmv", "gls", "wls", "wlsm", "wlsmv",
      "uls", "ulsm", "ulsmv", "pml", "dls"
    ))) {
      lav_msg_stop(gettextf(
        "missing=%1$s is not allowed for estimator %2$s",
        dQuote(opt$missing), dQuote(lav_options_estimatorgroup(opt$estimator))))
    }
  } else if (opt$missing == "two.stage") {
    if (opt$.categorical) {
      lav_msg_stop(gettextf(
        "missing = %s not available in the categorical setting",
        dQuote(opt$missing)))
    }
    if (any(opt$estimator == c(
      "mlm", "mlmv", "gls", "wls", "wlsm", "wlsmv",
      "uls", "ulsm", "ulsmv", "pml", "mml", "dls"
    ))) {
      lav_msg_stop(gettextf(
        "missing=%1$s is not allowed for estimator %2$s",
        dQuote(opt$missing), dQuote(lav_options_estimatorgroup(opt$estimator))))
    }
  } else if (opt$missing == "robust.two.stage") {
    if (opt$.categorical) {
      lav_msg_stop(gettextf(
        "missing = %s not available in the categorical setting",
        dQuote(opt$missing)))
    }
    if (any(opt$estimator == c(
      "mlm", "mlmv", "gls", "wls", "wlsm", "wlsmv",
      "uls", "ulsm", "ulsmv", "pml", "mml", "dls"
    ))) {
      lav_msg_stop(gettextf(
        "missing=%1$s is not allowed for estimator %2$s",
        dQuote(opt$missing), dQuote(lav_options_estimatorgroup(opt$estimator))))
    }
  } else if (opt$missing == "doubly.robust") {
    if (opt$estimator != "pml") {
      lav_msg_stop(gettextf(
        "missing=%s option only available for estimator PML",
        dQuote(opt$missing)))
    }
  }
  
  # check missing ####
  if (any(opt$missing == c("ml", "ml.x")) && opt$se == "robust.sem") {
    lav_msg_warn(gettextf(
      "missing will be set to %1$s for se = %2$s.",
      dQuote("listwise"), dQuote(opt$se)))
    opt$missing <- "listwise"
  }
  if (any(opt$missing == c("ml", "ml.x")) &&
      any(opt$test %in% c(
        "satorra.bentler",
        "mean.var.adjusted", "scaled.shifted"
      ))) {
    lav_msg_warn(gettextf(
      "missing will be set to %s for satorra.bentler style test",
      dQuote("listwise")))
    opt$missing <- "listwise"
  }
  
  # checks if missing = "two.stage" or "robust.two.stage" #### 
  if (any(opt$missing == c("two.stage", "robust.two.stage"))) {
    opt$meanstructure <- TRUE
    # se
    if (opt$se == "default") {
      if (opt$missing == "two.stage") {
        opt$se <- "two.stage"
      } else {
        opt$se <- "robust.two.stage"
      }
    } else if (opt$missing == "two.stage" &&
               opt$se == "two.stage") {
      # nothing to do
    } else if (opt$missing == "robust.two.stage" &&
               opt$se == "robust.two.stage") {
      # nothing to do
    } else {
      lav_msg_warn(gettextf(
        "se will be set to %1$s if missing = %2$s",
        dQuote(opt$missing), dQuote(opt$missing)))
      opt$se <- opt$missing
    }
    # information 
    if (opt$information[1] == "default") {
      # for both two.stage and robust.two.stage
      opt$information[1] <- "observed"
    } else if (opt$information[1] == "first.order") {
      lav_msg_warn(gettextf(
        "information will be set to %1$s if missing = %2$s",
        dQuote("observed"), dQuote(opt$missing)
      ))
      opt$information[1] <- "observed"
    }
    
    # observed.information (ALWAYS "h1" for now)
    opt$observed.information[1] <- "h1"
    opt$observed.information[2] <- "h1"
    
    # new in 0.6-9: ALWAYS h1.information = "unstructured"
    opt$h1.information <- c("unstructured", "unstructured")
    
    if (length(opt$information) > 1L && opt$information[2] == "default") {
      # for both two.stage and robust.two.stage
      opt$information[2] <- "observed"
    }
    
    # test
    if (length(opt$test) > 1L) {
      lav_msg_warn(gettextf(
        "test= argument can only contain a single element if missing = %s ",
        dQuote(opt$missing)), gettext("(taking the first)"))
      opt$test <- opt$test[1]
    }
    
    if (length(opt$test) == 1L && opt$test == "default") {
      opt$test <- "satorra.bentler"
    } else if (length(opt$test) == 1L && any(opt$test ==
                                             c("satorra", "sb", "satorra.bentler", "satorra-bentler"))) {
      opt$test <- "satorra.bentler"
    } else {
      lav_msg_warn(gettextf(
        "test will be set to %1$s if missing = %2$s",
        dQuote("satorra.bentler"), dQuote(opt$missing)
      ))
      opt$test <- "satorra.bentler"
    }
  }
  
  # meanstructure ####
  if (is.logical(opt$meanstructure)) {
    if (opt$meanstructure == FALSE) {
      if (any(opt$missing == c("ml", "ml.x", "two.stage"))) {
        lav_msg_warn(gettextf(
          "missing argument %s forces meanstructure = TRUE"),
          opt$missing)
      }
    }
  } else if (opt$meanstructure == "default") {
    # by default: no meanstructure!
    if (opt$estimator == "pml") {
      opt$meanstructure <- TRUE
    } else {
      opt$meanstructure <- FALSE
    }

  }
  
  # bootstrap ####
  if (opt$se == "bootstrap") {
    opt$information[1] <- "observed"
    if (length(opt$information) > 1L && opt$information[2] == "default") {
      opt$information[2] <- "observed"
    }
    opt$bootstrap <- as.integer(opt$bootstrap)
  }
  
  # specific per estimator (group) ####
  opt <- switch(opt$estimator, 
                ml = ,
                mlf = ,
                mlm = ,
                mlmv = ,
                mlmvs = ,
                mlr = lav_options_est_ml(opt),
                gls = lav_options_est_gls(opt),
                ntrls = lav_options_est_ntrls(opt),
                catml = lav_options_est_catml(opt),
                wls = lav_options_est_wls(opt),
                dls = lav_options_est_dls(opt),
                dwls = ,
                wlsm = ,
                wlsmv = ,
                wlsmvs = lav_options_est_dwls(opt),
                uls = ,
                ulsm = ,
                ulsmv = ,
                ulsmvs = lav_options_est_uls(opt),
                pml = lav_options_est_pml(opt),
                fml = lav_options_est_fml(opt),
                reml = lav_options_est_reml(opt),
                mml = lav_options_est_mml(opt),
                fabin2 = ,
                fabin3 = ,
                mgm = ,
                js = ,
                jsa = ,
                bentler1982 = lav_options_est_fabin(opt),
                miiv = lav_options_est_miiv(opt),
                lav_options_est_none(opt)  # estimator = none
  )
  
  # after code specific to estimator types                          ####
  # optim.method - if still "default" at this point -> set to "nlminb"
  if (opt$optim.method == "default") {
    opt$optim.method <- "nlminb"
  }
  
  # special stuff for categorical
  if (opt$.categorical) {
    opt$meanstructure <- TRUE # Mplus style
    if (lav_options_estimatorgroup(opt$estimator) == "ML") {
      lav_msg_stop(gettext(
        "estimator ML for ordered data is not supported yet. Use WLSMV instead."
      ))
    }
  }
  
  # link
  if (opt$link == "logit") {
    if (opt$estimator != "mml") {
      lav_msg_warn(gettextf(
        "link will be set to %1$s for estimator = %2$s",
        dQuote("probit"), dQuote(opt$estimator)
      ))
    }
  }
  
  # likelihood approach (wishart or normal) + sample.cov.rescale
  if (!any(lav_options_estimatorgroup(opt$estimator) == 
           c("ML", "REML", "PML", "FML", "NTRLS", "catML"))) {
    # if(opt$likelihood != "default") {
    #    lav_msg_stop(gettext(
    #    "likelihood argument is only relevant if estimator = ML"))
    # }
    if (opt$sample.cov.rescale == "default") {
      opt$sample.cov.rescale <- FALSE
    } # else {
    #    lav_msg_warn(gettext(
    #    "sample.cov.rescale argument is only relevant if estimator = ML"))
    # }
  } else { # ml and friends
    if (any(lav_options_estimatorgroup(opt$estimator) == c("PML", "FML"))) {
      opt$likelihood <- "normal"
    } else if (opt$likelihood == "default") {
      opt$likelihood <- "normal"
    }
    
    if (opt$sample.cov.rescale == "default") {
      opt$sample.cov.rescale <- FALSE
      if (opt$likelihood == "normal") {
        opt$sample.cov.rescale <- TRUE
      }
    }
  }
  
  # se information
  if (opt$information[1] == "default") {
    if (any(opt$missing == c("ml", "ml.x")) ||
        any(opt$se == c("robust.huber.white", "first.order"))) {
      # nchar(opt$constraints) > 0L) {
      opt$information[1] <- "observed"
    } else {
      opt$information[1] <- "expected"
    }
  }
  
  # first.order information can not be used with robust
  if (opt$information[1] == "first.order" &&
      any(opt$se == c("robust.huber.white", "robust.sem"))) {
    lav_msg_stop(gettextf(
      "information must be either %s if robust standard errors are requested.",
      lav_msg_view(c("expected", "observed"), log.sep = "or")))
  }
  
  # test information
  if (length(opt$information) == 1L) {
    opt$information <- rep(opt$information, 2L)
  }
  if (opt$information[2] == "default") {
    if (any(opt$missing == c("ml", "ml.x")) ||
        any(opt$se == c("robust.huber.white", "first.order"))) {
      # nchar(opt$constraints) > 0L) {
      opt$information[2] <- "observed"
    } else {
      opt$information[2] <- "expected"
    }
  }
  
  # first.order information cannot be used with robust
  if (opt$information[2] == "first.order" &&
      any(opt$test %in% c("satorra.bentler", "yuan.bentler", "yuan.bentler.mplus",
                          "mean.var.adjusted", "scaled.shifted"))) {
    lav_msg_stop(gettextf(
      "information must be either %s if robust test statistics are requested.",
      lav_msg_view(c("expected", "observed"), log.sep = "or")))
  }
  
  
  if (length(opt$observed.information) == 1L) {
    opt$observed.information <- rep(opt$observed.information, 2L)
  }
  
  if (all(opt$observed.information[2] != c("hessian", "h1"))) {
    if (opt$observed.information[2] == "default") {
      if (any(opt$test %in% c(
        "satorra.bentler",
        "yuan.bentler",
        "yuan.bentler.mplus",
        "mean.var.adjusted",
        "scaled.shifted"
      ))) {
        if (length(opt$test) > 1L) {
          opt$observed.information[2] <- "h1" # CHANGED in 0.6-6!
          if (any(opt$test == "yuan.bentler.mplus")) {
            lav_msg_warn(gettext(
              "observed.information for ALL test statistics is set to h1."))
          }
        } else {
          if (opt$estimator == "PML" ||
              opt$test[1] == "yuan.bentler.mplus") {
            opt$observed.information[2] <- "hessian"
          } else {
            opt$observed.information[2] <- "h1" # CHANGED in 0.6-6!
          }
        }
      } else {
        # default is "hessian"
        opt$observed.information[2] <- "hessian"
      }
    }
  }
  
  if (length(opt$h1.information) == 1L) {
    opt$h1.information <- rep(opt$h1.information, 2L)
  }
  
  if (opt$h1.information.meat == "default") {
    opt$h1.information.meat <- opt$h1.information[1]
  } 
  
  # check information if estimator is uls/wls and friends
  if (any(lav_options_estimatorgroup(opt$estimator) == c("ULS", "WLS", "DWLS"))) {
    if (opt$information[1] != "expected") {
      lav_msg_warn(gettextf(
        "information will be set to %1$s for estimator = %2$s",
        dQuote("expected"), dQuote(opt$estimator))
      )
      opt$information <- rep.int("expected", 2L)
    }
    opt$h1.information <- rep.int("unstructured", 2L) #FIXME: allow option?
  }
  
  
  # omega information 
  if (opt$omega.information == "default") {
    opt$omega.information <- opt$information[2] # test version!
  }
  
  if (opt$omega.h1.information == "default") {
    # opt$omega.h1.information <- opt$h1.information[2] # test version!
    opt$omega.h1.information <- "unstructured"
  }
  
  if (opt$omega.h1.information.meat == "default") {
    opt$omega.h1.information.meat <- opt$omega.h1.information
  }
  
  # conditional.x
  if (is.character(opt$conditional.x)) { # = "default"
    if (opt$.categorical) {
      opt$conditional.x <- TRUE
    } else {
      opt$conditional.x <- FALSE
    }
  }
  
  # if conditional.x, always use a meanstructure
  if (opt$conditional.x) {
    opt$meanstructure <- TRUE
  }
  
  # fixed.x
  if (is.logical(opt$fixed.x)) {
    # if(opt$conditional.x && opt$fixed.x == FALSE && !opt$.multilevel) {
    if (opt$conditional.x && opt$fixed.x == FALSE) {
      lav_msg_stop(gettext(
        "fixed.x = FALSE is not supported when conditional.x = TRUE."))
    }
    if (opt$fixed.x && is.character(opt$start) && opt$start == "simple") {
      lav_msg_warn(gettextf(
        "start = %s implies fixed.x = FALSE", dQuote(opt$start)))
      opt$fixed.x <- FALSE
    }
  } else if (opt$fixed.x == "default") {
    if (opt$conditional.x) {
      opt$fixed.x <- TRUE
    } else {
      opt$fixed.x <- FALSE
    }
  }
  
  # meanstructure again
  if (any(opt$missing == c("ml", "ml.x")) || opt$model.type == "growth") {
    opt$meanstructure <- TRUE
  }
  if (any(c("intercepts", "means") %in% opt$group.equal)) {
    opt$meanstructure <- TRUE
  }
  # if(opt$se == "robust.huber.white" ||
  #   opt$se == "robust.sem" ||
  #   opt$test == "satorra.bentler" ||
  #   opt$test == "mean.var.adjusted" ||
  #   opt$test == "scaled.shifted" ||
  #   opt$test == "yuan.bentler") {
  #    opt$meanstructure <- TRUE
  # }
  if (!is.logical(opt$meanstructure)) {
    lav_msg_fixme("meanstructure must be logical at this point!")
  }
  
  if (opt$debug) {
    opt$verbose <- opt$warn <- TRUE
  }
  
  # zero cell frequencies
  if (is.character(opt$zero.add)) { # = "default"
    opt$zero.add <- c(0.5, 0.0)
    # FIXME: TODO: mimic EQS , LISREL (0.0, 0.0)
  } else {
    if (length(opt$zero.add) == 1L) {
      opt$zero.add <- c(opt$zero.add, opt$zero.add)
    }
  }
  
  if (is.character(opt$zero.keep.margins)) { # = "default"
    opt$zero.keep.margins <- FALSE
  }


  # parameterization
  if (opt$parameterization == "default") {
    # for now, default is always delta
    opt$parameterization <- "delta"
  }
  
  # std.lv vs auto.fix.first # new in 0.6-5 (used to be in sem/cfa/growth)
  if (opt$std.lv) {
    opt$auto.fix.first <- FALSE
  }
  
  # std.lv vs effect.coding # new in 0.6-4
  if (is.logical(opt$effect.coding)) {
    if (opt$effect.coding) {
      opt$effect.coding <- c("loadings", "intercepts")
    } else {
      opt$effect.coding <- ""
    }
  }
  
  # if we use effect coding for the factor loadings, we don't need/want
  # std.lv = TRUE
  if (any("loadings" == opt$effect.coding)) {
    if (opt$std.lv) {
      lav_msg_stop(gettextf(
        "std.lv is set to FALSE but effect.coding contains %s", 
        dQuote("loadings")))
    }
    # shut off auto.fix.first
    opt$auto.fix.first <- FALSE
  }
  
  # test again
  
  # unless test = "none", always add test = "standard" as the
  # first entry
  # NO: this breaks lavaan.survey pval.pFsum, which has the following check:
  # if (!lavInspect(lavaan.fit, "options")$test %in% c("satorra.bentler",
  #    "mean.var.adjusted", "Satterthwaite")) {
  #    lav_msg_stop(
  #    gettext("Please refit the model with Satorra-Bentler (MLM)"), 
  #    gettext(" or Satterthwaite (MLMVS) adjustment."))
  # }
  # if(! (length(opt$test) == 1L && opt$test == "none") ) {
  #    opt$test <- c("standard", opt$test)
  #    opt$test <- unique(opt$test)
  # }
  
  # add scaled.test to test (if not already there)
  if (opt$scaled.test != "standard") {
    if (length(opt$test) == 1L && opt$test[1] == "standard") {
      opt$test <- unique(c(opt$test, opt$scaled.test))
    } else {
      opt$test <- unique(c(opt$scaled.test, opt$test))
    }
    
    # make sure "standard" comes first
    standard.idx <- which(opt$test == "standard")[1]
    if (length(standard.idx) > 0L && standard.idx != 1L) {
      opt$test <- c("standard", opt$test[-standard.idx])
    }
  }
  
  
  # final check
  wrong.idx <- which(!opt$test %in% c(
    "none", "standard", "satorra.bentler",
    "yuan.bentler", "yuan.bentler.mplus",
    "mean.var.adjusted", "scaled.shifted",
    "browne.residual.adf", "browne.residual.nt",
    "browne.residual.nt.model",
    "browne.residual.adf.model",
    "bollen.stine"
  ))
  if (length(wrong.idx) > 0L) {
    lav_msg_stop(
      gettextf("invalid option(s) for test argument: %s.",
               paste(dQuote(opt$test[wrong.idx]), collapse = " ")),
      gettextf(" Possible options are: %s.",
               lav_msg_view(c("none", "standard", "browne.residual.adf", 
                              "browne.residual.nt", "browne.residual.adf.model",
                              "browne.residual.nt.model", "satorra.bentler", "yuan.bentler",
                              "yuan.bentler.mplus", "mean.var.adjusted", "scaled.shifted",
                              "bollen.stine"), log.sep = "or"))
    )
  }
  
  # bounds
  if (is.null(opt$bounds)) {
    if (length(opt$optim.bounds) > 0L) {
      opt$bounds <- "user"
    } else {
      opt$bounds <- "none" # for now
    }
  } else if (is.logical(opt$bounds)) {
    if (opt$bounds) {
      opt$bounds <- "wide" # default for most estimators
    } else {
      opt$bounds <- "none"
    }
  }
  
  # optim.bounds
  if (length(opt$optim.bounds) > 0L) {
    # opt$bounds should be "default", or "user" (or "none")
    if (any(opt$bounds == c("default", "none", "user"))) {
      opt$bounds <- "user"
    } else {
      lav_msg_stop(
        gettext("bounds and optim.bounds arguments can not be used together;"),
        gettext("remove the bounds= argument or set it to \"user\".")
      )
    }
  }
  
  # handle different 'profiles'
  if (opt$bounds == "none") {
    opt$optim.bounds <- list(
      lower = character(0L),
      upper = character(0L)
    )
  } else if (opt$bounds == "user") {
    if (length(opt$optim.bounds) == 0L) {
      lav_msg_stop(gettextf(
        "bounds= is %s but optim.bounds= argument is empty",  dQuote("user")
      ))
    }
  } else if (opt$bounds == "default" || opt$bounds == "wide") {
    opt$optim.bounds <- list(
      lower = c("ov.var", "lv.var", "loadings", "covariances"),
      upper = c("ov.var", "lv.var", "loadings", "covariances"),
      lower.factor = c(1.05, 1.0, 1.1, 1.0),
      upper.factor = c(1.20, 1.3, 1.1, 1.0),
      min.reliability.marker = 0.1,
      min.var.lv.endo = 0.005
    )
  } else if (opt$bounds == "wide.zerovar") {
    opt$optim.bounds <- list(
      lower = c("ov.var", "lv.var", "loadings", "covariances"),
      upper = c("ov.var", "lv.var", "loadings", "covariances"),
      lower.factor = c(1.00, 1.0, 1.1, 1.0),
      upper.factor = c(1.20, 1.3, 1.1, 1.0),
      min.reliability.marker = 0.1,
      min.var.lv.endo = 0.005
    )
  } else if (opt$bounds == "standard") {
    opt$optim.bounds <- list(
      lower = c("ov.var", "lv.var", "loadings", "covariances"),
      upper = c("ov.var", "lv.var", "loadings", "covariances"),
      lower.factor = c(1.0, 1.0, 1.0, 0.999),
      upper.factor = c(1.0, 1.0, 1.0, 0.999),
      min.reliability.marker = 0.1,
      min.var.lv.endo = 0.005
    )
  } else if (opt$bounds == "pos.var") {
    opt$optim.bounds <- list(
      lower = c("ov.var", "lv.var"),
      lower.factor = c(1, 1),
      min.reliability.marker = 0.0,
      min.var.lv.exo = 0.0,
      min.var.lv.endo = 0.0
    )
  } else if (opt$bounds == "pos.ov.var") {
    opt$optim.bounds <- list(
      lower = c("ov.var"),
      lower.factor = 1
    )
  } else if (opt$bounds == "pos.lv.var") {
    opt$optim.bounds <- list(
      lower = c("lv.var"),
      lower.factor = 1,
      min.reliability.marker = 0.0,
      min.var.lv.exo = 0.0,
      min.var.lv.endo = 0.0
    )
  } 
  
  # rotations.args
  if (!is.list(opt$rotation.args)) {
    lav_msg_stop(gettext("rotation.args should be be list."))
  }
  
  # force orthogonal for some rotation algorithms
  if (any(opt$rotation == c("varimax", "entropy", "mccammon",
                            "tandem1", "tandem2"))) {
    opt$rotation.args$orthogonal <- TRUE
  }
  
  # if target, check target matrix
  if (opt$rotation == "target" || opt$rotation == "pst") {
    target <- opt$rotation.args$target
    if (is.null(target) || !is.matrix(target)) {
      lav_msg_stop(gettext("rotation target matrix is NULL, or not a matrix"))
    }
  }
  if (opt$rotation == "pst") {
    target.mask <- opt$rotation.args$target.mask
    if (is.null(target.mask) || !is.matrix(target.mask)) {
      lav_msg_stop(gettext(
        "rotation target.mask matrix is NULL, or not a matrix"
      ))
    }
  }
  # if NAs, force opt$rotation to be 'pst' and create target.mask
  if (opt$rotation == "target" && anyNA(target)) {
    opt$rotation <- "pst"
    target.mask <- matrix(1, nrow = nrow(target), ncol = ncol(target))
    target.mask[is.na(target)] <- 0
    opt$rotation.args$target.mask <- target.mask
  }
  
  # set row.weights
  opt$rotation.args$row.weights <- tolower(opt$rotation.args$row.weights)
  if (opt$rotation.args$row.weights == "default") {
    # the default is "none", except for varimax and promax
    if (any(opt$rotation == c("varimax", "promax"))) {
      opt$rotation.args$row.weights <- "kaiser"
    } else {
      opt$rotation.args$row.weights <- "none"
    }
  }
  
  # override if bifactor
  if (any(opt$rotation == c("bi-geomin", "bigeomin", "bi-quartimin", 
                            "biquartimin"))) {
    opt$rotation.args$order.lv.by <- "none"
  }
  
  # no standard errors for promax (for now)...
  if (opt$rotation == "promax") {
    opt$se <- "none"
    opt$rotation.args$algorithm <- "promax"
    opt$rotation.args$rstarts <- 0L
  }
  
  # correlation
  if (opt$correlation) {
    if (opt$missing == "ml") {
      lav_msg_stop(gettext(
        "correlation structures only work for complete data (for now)."))
    }
    if (opt$.multilevel) {
      lav_msg_stop(gettext(
        "correlation structures only work for single-level data."))
    }
    if (opt$conditional.x) {
      lav_msg_stop(gettext(
        "correlation structures only work for conditional.x = FALSE (for now)."
      ))
    }
    if (opt$representation == "RAM") {
      lav_msg_stop(gettext(
        "correlation structures only work for representation = \"LISREL\"."))
    }
    if (opt$fixed.x) {
      # first fix eliminate.pstar.idx in lav_mvnorm_information_expected()
      lav_msg_stop(gettext(
        "correlation structures only work for fixed.x = FALSE (for now)."))
    }
  }
  
  # sample.cov.robust
  # new in 0.6-17
  # sample.cov.robust cannot be used if:
  # - data is missing (for now),
  # - sampling weights are used
  # - estimator is (D)WLS
  # - multilevel
  # - conditional.x
  if (opt$sample.cov.robust) {
    if (opt$missing != "listwise") {
      lav_msg_stop(gettext(
        "sample.cov.robust = TRUE does not work (yet) if data is missing."))
    }
    if (opt$.categorical) {
      lav_msg_stop(gettext(
        "sample.cov.robust = TRUE does not work (yet) if data is categorical"))
    }
    if (opt$.clustered || opt$.multilevel) {
      lav_msg_stop(gettext(
        "sample.cov.robust = TRUE does not work (yet) if data is clustered"))
    }
    if (opt$conditional.x) {
      lav_msg_stop(gettext(
        "sample.cov.robust = TRUE does not work (yet) if conditional.x = TRUE"))
    }
    if (all(lav_options_estimatorgroup(opt$estimator) != c("ML", "GLS"))) {
      lav_msg_stop(gettext(
        "sample.cov.robust = TRUE does not work (yet) if estimator is not GLS or ML"))
    }
  }
  
  opt$estimator <- lav_options_estimatorgroup(opt$estimator)
  
  # group.w.free
  # if(opt$group.w.free && opt$.categorical) {
  #    lav_msg_stop(gettext(
  #    "group.w.free = TRUE is not supported (yet) in the categorical setting."
  #    ))
  # }
  
  # in order not to break semTools and blavaan, we restore categorical:
  opt$categorical <- opt$.categorical
  
  if (opt$debug) {
    cat("lavaan DEBUG: lavaanOptions OUT\n")
    str(opt)
  }
  
  opt
}
