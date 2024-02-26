lav_lavaan_step00_checkdata <- function(data, dotdotdot, sample.cov, sample.nobs,
                                      sample.mean, sample.th, NACOV, WLS.V, ov.order) {
  # als data niet NULL is :
  #   als het een 'uitgebreide' data.frame is (bvb tibble), vereenvoudig naar gewone data.frame
  #   als het klasse 'lavMoments' heeft : 
  #             check of het sample.cov en sample.nobs bevat (***error*** indien niet)
  #             indien sample.mean, sample.th, NACOV en/of WLS.V aanwezig,
  #                       kopieren naar corr. argumenten vd functie
  #             indien lavOptions aanwezig, deze die niet meegegeven zijn in de functie
  #                   oproep (dotdotdot) er naartoe kopieren
  #             zet data = NULL
  #   als het een functie is --> ***error***
  #   TODO: andere tests gebeuren in lavData, misschien beter alles op één plaats ???
  # als NACOV of WLS.V niet NULL, zet order = "data" 
  # als dotdotdot$control aanwezig overhevelen naar dotdotdot$... van
  #   optim.method, 
  #   optim.force.converged, 
  #   gradient -> optim.gradient ! overschreven door dotdotdot$gradient indien aanwezig
  #   init_nelder_mead -> optim.init_nelder_mead
  #   TODO: laatste stap gebruikt enkel parameters en heeft niets te maken met data parameter, kan dus vroeger
  
  # check data
  if (!is.null(data)) {
    if (inherits(data, "data.frame")) {
      # just in case it is not a traditional data.frame
      data <- as.data.frame(data)
      
    } else if (inherits(data, "lavMoments")) {
      # This object must contain summary statistics
      # e.g., created by lavaan.mi::poolSat
      
      # set required-data arguments
      if ("sample.cov" %in% names(data)) {
        sample.cov  <- data$sample.cov
      } else {
        stop("When data= is of class lavMoments, it must contain sample.cov")
      }
      
      if ("sample.nobs" %in% names(data)) {
        sample.nobs <- data$sample.nobs
      } else {
        stop("When data= is of class lavMoments, it must contain sample.nobs")
      }
      
      # check for optional-data arguments
      if ("sample.mean" %in% names(data)) sample.mean <- data$sample.mean
      if ("sample.th"   %in% names(data)) sample.th   <- data$sample.th
      if ("NACOV"       %in% names(data)) NACOV       <- data$NACOV
      if ("WLS.V"       %in% names(data)) WLS.V       <- data$WLS.V
      
      # set other args not included in dotdotdot
      if (length(data$lavOptions)) {
        newdots <- setdiff(names(data$lavOptions), names(dotdotdot))
        if (length(newdots)) {
          for (dd in newdots) dotdotdot[[dd]] <- data$lavOptions[[dd]]
        }
      }
      
      #FIXME: Should WLS.V be an I(dentity) matrix when ULS is requested?
      #       Unused for point estimates, but still used to scale/shift test
      # if (!is.null(dotdotdot$estimator)) {
      #   if (grepl(pattern = "ULS", x = toupper(dotdotdot$estimator[1L])) &&
      #       !is.null(WLS.V)) {
      #     #  set to diagonal
      #     if (is.list(WLS.V)) {
      #       WLS.V <- lapply(WLS.V, function(w) {diag(w) <- 1 ; return(w) })
      #     } else diag(WLS.V) <- 1
      #   }
      # }
      
      # get rid of data= argument
      data <- NULL
    }
    
    if (is.function(data)) {
      stop("lavaan ERROR: data is a function; it should be a data.frame")
    }
  }
  # new in 0.6-14: if NACOV and/or WLS.V are provided, we force
  # ov.order="data" for now
  # until we have reliable code to re-arrange/select col/rows for
  # of NACOV/WLS.V based on the model-based ov.names
  if (!is.null(NACOV) || !is.null(WLS.V)) {
    ov.order <- "data"
  }
  # backwards compatibility, control= argument (<0.5-23)
  if (!is.null(dotdotdot$control)) {
    # optim.method
    if (!is.null(dotdotdot$control$optim.method)) {
      dotdotdot$optim.method <- dotdotdot$control$optim.method
    }
    # cor.optim.method
    if (!is.null(dotdotdot$control$cor.optim.method)) {
      # ignore it silently
    }
    # control$optim.force.converged
    if (!is.null(dotdotdot$control$optim.force.converged)) {
      dotdotdot$optim.force.converged <-
        dotdotdot$control$optim.force.converged
    }
    # gradient
    if (!is.null(dotdotdot$control$gradient)) {
      dotdotdot$optim.gradient <- dotdotdot$control$gradient
    }
    if (!is.null(dotdotdot$gradient)) {
      dotdotdot$optim.gradient <- dotdotdot$gradient
    }
    # init_nelder_mead
    if (!is.null(dotdotdot$control$init_nelder_mead)) {
      dotdotdot$optim.init_nelder_mead <-
        dotdotdot$control$init_nelder_mead
    }
  }
  return(list(data = data, dotdotdot = dotdotdot, sample.cov = sample.cov, sample.nobs = sample.nobs,
              sample.mean = sample.mean, sample.th = sample.th, NACOV = NACOV, WLS.V = WLS.V,
              ov.order = ov.order))
}
attr(lav_lavaan_step00_checkdata, "output") <- c("data", "dotdotdot", "sample.cov",
                                               "sample.nobs", "sample.mean", "sample.th",
                                               "NACOV", "WLS.V", "ov.order")
