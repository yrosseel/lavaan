lav_label_code <- function(label = "", value = "", show = FALSE,
                           idx.font.size = 20L, dy = 7L,
                           italic = TRUE, auto.subscript = TRUE) {
  latexsymbols <- c("varepsilon",
    "Alpha", "Beta", "Gamma", "Delta", "Epsilon", "Zeta", "Eta", "Theta",
    "Iota", "Kappa", "Lambda", "Mu", "Nu", "Xi", "Omicron", "Pi",
    "Rho", "Sigma", "Tau", "Upsilon", "Phi", "Chi", "Psi", "Omega",
    "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta", "theta",
    "iota", "kappa", "lambda", "mu", "nu", "xi", "omicron", "pi",
    "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega"
  )
  unicodes <- c("1013",
    "913", "914", "915", "916", "917", "918", "919", "920",
    "921", "922", "923", "924", "925", "926", "927", "928",
    "929", "931", "932", "933", "934", "935", "936", "937",
    "945", "946", "947", "948", "949", "950", "951", "952",
    "953", "954", "955", "956", "957", "958", "959", "960",
    "961", "963", "964", "965", "966", "967", "968", "969"
  )
  if (label == "" && value == "") return(list(svg = "", tikz = "", r = ""))
  if (auto.subscript) {
    label <- gsub("^([a-zA-Z]+)([0-9]+)", "\\1_\\2", label)
  }
  rexpression <- FALSE
  if (grepl("1van", label)) label <- "1"
  if (value == "") {
    splitted <- strsplit(label, "=", fixed =  TRUE)[[1L]]
    if (length(splitted) > 1L) {
      label <- splitted[1L]
      value <- splitted[2L]
    }
  }
  if (label == "") {
    label <- value
    value <- ""
  }
  svgpart <- tikzpart <- rpart <- character(3)
  rpart[3L] <- tikzpart[3L] <- svgpart[3L] <- value
  if (label != "") {
    # separate label and value by equal sign
    if (svgpart[3L] != "") svgpart[3L] <- paste0("=", svgpart[3L])
    if (tikzpart[3L] != "") tikzpart[3L] <- paste0("=", tikzpart[3L])
    if (rpart[3L] != "") {
      rpart[3L] <- paste0(" == ", rpart[3L])
      rexpression <- TRUE
    }
    # subscript and Greek character set handling
    splitunderscore <- strsplit(label, "_", TRUE)[[1L]]
    svgpart[1] <- tikzpart[1] <- rpart[1] <- splitunderscore[1L]
    if (rpart[1] %in% latexsymbols) {
      rexpression <- TRUE
      svgpart[1] <- paste0("&#", unicodes[latexsymbols == svgpart[1]], ";")
      tikzpart[1] <- paste0("\\", tikzpart[1])
      if (rpart[1] == "varepsilon") rpart[1] <- "epsilon"
      rpart[1] <- str2expression(rpart[1])
    }
    if (length(splitunderscore) > 1L) {
      rexpression <- TRUE
      svgpart[2L] <- paste0("<tspan dy=\"", dy  ,"\" font-size=\"",
                     idx.font.size, "\">", splitunderscore[2L], "</tspan>")
      if (svgpart[3L] != "") {
        svgpart[3L] <- paste0("<tspan dy=\"-", dy, "\">",
                               svgpart[3L], "</tspan>")
      }
      tikzpart[2L] <- paste0("_{", splitunderscore[2L], "}")
      rpart[2L] <- paste0("[\"", splitunderscore[2L], "\"]")
    }
    if (!italic) tikzpart <- c("\\mathrm{", tikzpart, "}")
  }
  if (show) {
    plot(c(0,3), c(0,2), type = "n")
    if (rexpression) {
      text(1, 1, str2expression(paste(rpart, collapse = "")))
    } else {
      text(1, 1, paste(rpart, collapse = ""), font = ifelse(italic, 3L, 1L))
    }
  }
  list(svg = paste(svgpart, collapse = ""),
       tikz = paste(c("$", tikzpart, "$"), collapse = ""),
       r = ifelse(rexpression,
         str2expression(paste(rpart, collapse = "")),
         paste(rpart, collapse = ""))
  )
}
