# print 'blocks' of test statistics
# - blocks with 'scaling.factors' come first (in 'two columns')
# - then come the 'single-column' test statistics (eg browne.residual.adf)
# - print additional information (eg information matrix, h1.information, ...)
#   if they deviate from what is used for the standard errors

# this is used by the summary() function and lavTest(, output = "text")

# YR 13 Feb 2026: if test[[1]] has stat = NA, skip if there are other tests

lav_test_print <- function(object, nd = 3L) {
  # object is list of tests
  test <- object

  # empty list?
  if (is.null(test) || length(test) == 0L || !is.list(test)) {
    return(character(0L))
  }

  # test = "none"?
  if (test[[1]]$test == "none") {
    return(character(0L))
  }

  # remove empty first test (stat = NA) if multiple tests are available)
  if (length(test) > 1L && is.na(test[[1]]$stat)) {
    test <- test[-1]
  }

  # meta data
  info <- attr(object, "info")
  ngroups <- info$ngroups
  group_label <- info$group.label
  information <- info$information
  h1_information <- info$h1.information
  observed_information <- info$observed.information

  # num format
  num_format <- paste("%", max(8L, nd + 5L), ".", nd, "f", sep = "")

  # header
  cat("Model Test User Model:\n")

  # locate 'robust' tests (here: having a scaling factor)
  has_no_scaling <- unname(sapply(
    (lapply(test, "[[", "scaling.factor")),
    is.null
  ))
  robust_idx <- which(!has_no_scaling)
  non_robust_idx <- which(has_no_scaling)
  scaled_idx <- 1L
  if (length(robust_idx) > 0L) {
    scaled_idx <- which(names(test) == test[[robust_idx[1]]]$scaled.test)
    if (length(scaled_idx) == 0L) {
      scaled_idx <- 1L
    }
    # remove 'scaled.test', because it is shown together with robust
    non_robust_idx <- non_robust_idx[-scaled_idx]
  }
  blocks <- c(robust_idx, non_robust_idx)
  # n_blocks <- length(blocks)

  # print out blocks
  for (block in blocks) {
    # one or two-columns for this block?
    if (length(robust_idx) > 0L && block %in% robust_idx) {
      twocolumn <- TRUE
    } else {
      twocolumn <- FALSE
    }

    if (!twocolumn) {
      # print label
      c1 <- c2 <- c3 <- character(0L)
      if (!is.null(test[[block]]$label)) {
        c1 <- c(c1, test[[block]]$label)
        c2 <- c(c2, "")
        c3 <- c(c3, "")
      }
      if (is.na(test[[block]]$df) || test[[block]]$df == 0L) {
        c1 <- c(c1, c("Test statistic", "Degrees of freedom"))
        c2 <- c(c2, c(
          sprintf(num_format, test[[block]]$stat),
          ifelse(test[[block]]$df %% 1 == 0, # integer
            test[[block]]$df,
            sprintf(num_format, test[[block]]$df)
          )
        ))
        c3 <- c(c3, c("", ""))
      } else {
        plabel <- "P-value"
        if (!is.null(test[[block]]$refdistr)) {
          if (test[[block]]$refdistr == "chisq") {
            plabel <- "P-value (Chi-square)"
          } else if (test[[block]]$refdistr == "unknown") {
            plabel <- "P-value (Unknown)"
          } else if (test[[block]]$refdistr == "bootstrap") {
            plabel <- "P-value (Bollen-Stine bootstrap)"
          }
        }
        c1 <- c(c1, c("Test statistic", "Degrees of freedom", plabel))
        c2 <- c(c2, c(
          sprintf(num_format, test[[block]]$stat),
          ifelse(test[[block]]$df %% 1 == 0, # integer
            test[[block]]$df,
            sprintf(num_format, test[[block]]$df)
          ),
          sprintf(num_format, test[[block]]$pvalue)
        ))
        c3 <- c(c3, c("", "", ""))
      }

      # two-column
    } else {
      # print label
      c1 <- c2 <- c3 <- character(0L)
      if (!is.null(test[[scaled_idx]]$label)) {
        c1 <- c(c1, test[[scaled_idx]]$label)
        c2 <- c(c2, "")
        c3 <- c(c3, "")
      }
      if (is.na(test[[block]]$df) || test[[block]]$df == 0L) {
        c1 <- c(c1, c("Test Statistic", "Degrees of freedom"))
        c2 <- c(
          c2,
          c(
            sprintf(num_format, test[[scaled_idx]]$stat),
            ifelse(test[[scaled_idx]]$df %% 1 == 0, # integer
              test[[scaled_idx]]$df,
              sprintf(num_format, test[[scaled_idx]]$df)
            )
          )
        )
        c3 <- c(
          c3,
          c(
            sprintf(num_format, test[[block]]$stat),
            ifelse(test[[block]]$df %% 1 == 0, # integer
              test[[block]]$df,
              sprintf(num_format, test[[block]]$df)
            )
          )
        )
      } else {
        if (!is.null(test[[scaled_idx]]$refdistr)) {
          if (test[[scaled_idx]]$refdistr == "chisq") {
            plabel <- "P-value (Chi-square)"
          } else if (test[[scaled_idx]]$refdistr == "unknown") {
            plabel <- "P-value (Unknown)"
          } else {
            plabel <- "P-value"
          }
        }
        c1 <- c(c1, c(
          "Test Statistic", "Degrees of freedom", plabel,
          "Scaling correction factor"
        ))
        c2 <- c(
          c2,
          c(
            sprintf(num_format, test[[scaled_idx]]$stat),
            ifelse(test[[scaled_idx]]$df %% 1 == 0, # integer
              test[[scaled_idx]]$df,
              sprintf(num_format, test[[scaled_idx]]$df)
            ),
            sprintf(num_format, test[[scaled_idx]]$pvalue), ""
          )
        )
        c3 <- c(
          c3,
          c(
            sprintf(num_format, test[[block]]$stat),
            ifelse(test[[block]]$df %% 1 == 0, # integer
              test[[block]]$df,
              sprintf(num_format, test[[block]]$df)
            ),
            sprintf(num_format, test[[block]]$pvalue),
            sprintf(num_format, test[[block]]$scaling.factor)
          )
        )

        if (test[[block]]$test == "scaled.shifted") {
          if (ngroups == 1L ||
            length(test[[block]]$shift.parameter) == 1L) {
            c1 <- c(c1, "Shift parameter")
            c2 <- c(c2, "")
            c3 <- c(
              c3,
              sprintf(num_format, test[[block]]$shift.parameter)
            )
          } else {
            c1 <- c(c1, "Shift parameter for each group:")
            c2 <- c(c2, "")
            c3 <- c(c3, "")
            for (g in 1:ngroups) {
              c1 <- c(c1, sprintf("    %-38s", group_label[[g]]))
              c2 <- c(c2, "")
              c3 <- c(c3, sprintf(
                num_format,
                test[[block]]$shift.parameter[g]
              ))
            }
          }
        } # shift

        # which correction factor?
        c1 <- c(c1, paste("  ", test[[block]]$label, sep = ""))
        c2 <- c(c2, "")
        c3 <- c(c3, "")
      }
    }

    # if twocolumn, label first row
    if (twocolumn && block == blocks[1]) {
      c1 <- c("", c1)
      c2 <- c("Standard", c2)
      c3 <- c("Scaled", c3)
    } else {
      # empty row
      c1 <- c("", c1)
      c2 <- c("", c2)
      c3 <- c("", c3)
    }

    # if information type is different from 'se', print it
    if (length(information) > 1L &&
      information[1] != information[2]) {
      c1 <- c(c1, "Information")
      tmp_txt <- information[2]
      c2 <- c(c2, paste(toupper(substring(tmp_txt, 1, 1)),
        substring(tmp_txt, 2),
        sep = ""
      ))
      c3 <- c(c3, "")
    }
    # if h1.information type is different from 'se', print it
    if (length(h1_information) > 1L &&
      h1_information[1] != h1_information[2]) {
      c1 <- c(c1, "Information saturated (h1) model")
      tmp_txt <- h1_information[2]
      c2 <- c(c2, paste(toupper(substring(tmp_txt, 1, 1)),
        substring(tmp_txt, 2),
        sep = ""
      ))
      c3 <- c(c3, "")
    }
    # if observed.information type is different from 'se', print it
    if (length(observed_information) > 1L &&
      information[2] == "observed" &&
      (observed_information[1] !=
        observed_information[2])) {
      c1 <- c(c1, "Observed information based on")
      tmp_txt <- observed_information[2]
      c2 <- c(c2, paste(toupper(substring(tmp_txt, 1, 1)),
        substring(tmp_txt, 2),
        sep = ""
      ))
      c3 <- c(c3, "")
    }


    # format c1/c2/c3 (note: fitMeasures uses 35/16/8)
    c1 <- format(c1, width = 43L)
    c2 <- format(c2, width = 8L + max(0, (nd - 3L)) * 4L, justify = "right")
    c3 <- format(c3, width = 8L + nd, justify = "right")

    # create character matrix
    if (twocolumn) {
      m <- cbind(c1, c2, c3, deparse.level = 0)
    } else {
      m <- cbind(c1, c2, deparse.level = 0)
    }
    colnames(m) <- rep("", ncol(m))
    rownames(m) <- rep(" ", nrow(m))

    # print
    write.table(m, row.names = TRUE, col.names = FALSE, quote = FALSE)

    # multiple groups?
    ngroups <- ngroups
    if (ngroups > 1L && !is.null(test[[block]]$stat.group)) {
      c1 <- c2 <- c3 <- character(ngroups)
      for (g in 1:ngroups) {
        tmp <- sprintf("  %-40s", group_label[[g]])
        c1[g] <- format(tmp, width = 43L)
        if (!twocolumn) {
          tmp <- sprintf(num_format, test[[block]]$stat.group[g])
          c2[g] <- format(tmp,
            width = 8L + max(0, (nd - 3L)) * 4L,
            justify = "right"
          )
        } else {
          tmp <- sprintf(num_format, test[[block]]$stat.group[g])
          c2[g] <- format(tmp,
            width = 8L + max(0, (nd - 3L)) * 4L,
            justify = "right"
          )
          tmp <- sprintf(num_format, test[[block]]$stat.group[g])
          c3[g] <- format(tmp, width = 8L + nd, justify = "right")
        }
      }
      if (twocolumn) {
        m <- cbind(c1, c2, c3, deparse.level = 0)
      } else {
        m <- cbind(c1, c2, deparse.level = 0)
      }
      colnames(m) <- rep("", ncol(m))
      rownames(m) <- rep(" ", nrow(m))
      cat("  Test statistic for each group:\n")
      write.table(m, row.names = TRUE, col.names = FALSE, quote = FALSE)
    }
  } # blocks

  # invisible(M)
}
