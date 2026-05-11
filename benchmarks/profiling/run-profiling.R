#!/usr/bin/env Rscript

required_packages <- c(
  "bench", "profmem", "pkgload", "withr", "openxlsx"
)

check_required_packages <- function(packages = required_packages) {
  missing_packages <- packages[
    !vapply(packages, requireNamespace, logical(1L), quietly = TRUE)
  ]
  if (length(missing_packages) > 0L) {
    stop(
      "Missing required profiling packages: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    description <- file.path(current, "DESCRIPTION")
    if (file.exists(description)) {
      desc <- read.dcf(description)
      if (identical(unname(desc[1L, "Package"]), "lavaan")) {
        return(current)
      }
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not find lavaan package root from ", start, call. = FALSE)
    }
    current <- parent
  }
}

normalize_candidate_path <- function(path, base) {
  expanded <- path.expand(path)
  is_windows_abs <- grepl("^[A-Za-z]:[\\/]", expanded)
  is_unix_abs <- startsWith(expanded, "/")
  candidate <- if (is_windows_abs || is_unix_abs) {
    expanded
  } else {
    file.path(base, expanded)
  }
  collapse_path(candidate)
}

collapse_path <- function(path) {
  path <- chartr("\\", "/", path)
  parts <- strsplit(path, "/", fixed = TRUE)[[1L]]
  prefix <- character(0L)
  if (length(parts) > 0L && grepl("^[A-Za-z]:$", parts[[1L]])) {
    prefix <- parts[[1L]]
    parts <- parts[-1L]
  } else if (startsWith(path, "/")) {
    prefix <- ""
  }

  stack <- character(0L)
  for (part in parts) {
    if (!nzchar(part) || identical(part, ".")) {
      next
    }
    if (identical(part, "..")) {
      if (length(stack) > 0L && !identical(stack[[length(stack)]], "..")) {
        stack <- stack[-length(stack)]
      } else {
        stack <- c(stack, part)
      }
    } else {
      stack <- c(stack, part)
    }
  }

  if (length(prefix) > 0L && nzchar(prefix)) {
    paste0(prefix, "/", paste(stack, collapse = "/"))
  } else if (length(prefix) > 0L) {
    paste0("/", paste(stack, collapse = "/"))
  } else {
    paste(stack, collapse = "/")
  }
}

path_is_within <- function(path, parent) {
  path_norm <- tolower(chartr("\\", "/", path))
  parent_norm <- tolower(chartr("\\", "/", parent))
  identical(path_norm, parent_norm) ||
    startsWith(paste0(path_norm, "/"), paste0(parent_norm, "/"))
}

windows_write_path <- function(path) {
  path
}

windows_shell_path <- function(path) {
  path <- sub("^//\\?/", "", path)
  chartr("/", "\\", path)
}

windows_mkdir <- function(path) {
  if (.Platform$OS.type != "windows") {
    return(FALSE)
  }
  cmd_path <- windows_shell_path(path)
  debug <- identical(Sys.getenv("LAVAAN_PROFILE_DEBUG"), "1")
  status <- suppressWarnings(system2(
    "cmd.exe",
    c("/c", "mkdir", shQuote(cmd_path)),
    stdout = if (debug) TRUE else FALSE,
    stderr = if (debug) TRUE else FALSE
  ))
  if (debug) {
    cat(
      "windows_mkdir status=", paste(status, collapse = " "),
      " cmd_path=", cmd_path,
      "\n",
      sep = ""
    )
  }
  identical(status, 0L) || dir.exists(path)
}

copy_file <- function(from, to) {
  if (file.copy(from, to, overwrite = TRUE)) {
    return(TRUE)
  }
  if (.Platform$OS.type != "windows") {
    return(FALSE)
  }
  status <- suppressWarnings(system2(
    "cmd.exe",
    c(
      "/c", "copy", "/Y",
      shQuote(windows_shell_path(from)),
      shQuote(windows_shell_path(to))
    ),
    stdout = FALSE,
    stderr = FALSE
  ))
  identical(status, 0L) && file.exists(to)
}

parse_args <- function(args, project_root) {
  run_label <- if (length(args) >= 1L && nzchar(args[[1L]])) {
    args[[1L]]
  } else {
    "current-branch"
  }

  mode <- if (length(args) >= 2L && nzchar(args[[2L]])) {
    args[[2L]]
  } else {
    "quick"
  }
  if (!mode %in% c("quick", "full")) {
    stop("mode must be 'quick' or 'full'; got '", mode, "'", call. = FALSE)
  }

  default_iterations <- if (identical(mode, "full")) 900L else 1L
  iterations <- if (length(args) >= 3L && nzchar(args[[3L]])) {
    as.integer(args[[3L]])
  } else {
    default_iterations
  }
  if (is.na(iterations) || iterations < 1L) {
    stop("iterations must be a positive integer", call. = FALSE)
  }

  output_root_arg <- if (length(args) >= 4L && nzchar(args[[4L]])) {
    args[[4L]]
  } else {
    "../profiling-output"
  }
  output_root_abs <- normalize_candidate_path(output_root_arg, project_root)
  output_dir_abs <- file.path(output_root_abs, run_label)
  if (path_is_within(output_dir_abs, project_root)) {
    stop(
      "Refusing to write profiling artifacts inside the lavaan directory: ",
      output_dir_abs,
      call. = FALSE
    )
  }

  scenario_filter <- if (length(args) >= 5L && nzchar(args[[5L]])) {
    strsplit(args[[5L]], ",", fixed = TRUE)[[1L]]
  } else {
    "all"
  }
  scenario_filter <- trimws(scenario_filter)

  log_root_arg <- if (length(args) >= 6L && nzchar(args[[6L]])) {
    args[[6L]]
  } else {
    "../logs"
  }
  log_root_abs <- normalize_candidate_path(log_root_arg, project_root)
  log_file <- file.path(log_root_abs, paste0(run_label, ".log"))
  if (path_is_within(log_file, project_root)) {
    stop(
      "Refusing to write profiling logs inside the lavaan directory: ",
      log_file,
      call. = FALSE
    )
  }

  list(
    run_label = run_label,
    mode = mode,
    iterations = iterations,
    output_root = output_root_arg,
    output_root_abs = output_root_abs,
    output_dir = windows_write_path(output_dir_abs),
    output_dir_abs = output_dir_abs,
    workbook_path = windows_write_path(
      file.path(output_dir_abs, "profiling-summary.xlsx")
    ),
    workbook_path_abs = file.path(output_dir_abs, "profiling-summary.xlsx"),
    log_root = log_root_arg,
    log_root_abs = log_root_abs,
    log_file = windows_write_path(log_file),
    log_file_abs = log_file,
    scenario_filter = scenario_filter
  )
}

ensure_dir <- function(path, attempts = 10L) {
  for (attempt in seq_len(attempts)) {
    created <- dir.create(path, recursive = TRUE, showWarnings = FALSE)
    exists <- dir.exists(path)
    if (!exists) {
      created <- windows_mkdir(path)
      exists <- dir.exists(path)
    }
    if (identical(Sys.getenv("LAVAAN_PROFILE_DEBUG"), "1")) {
      cat(
        "ensure_dir attempt=", attempt,
        " created=", created,
        " exists=", exists,
        " cwd=", getwd(),
        " path=", path,
        "\n",
        sep = ""
      )
    }
    if (exists) {
      return(invisible(path))
    }
    Sys.sleep(0.1 * attempt)
  }
  stop("Could not create directory: ", path, call. = FALSE)
}

write_with_retries <- function(path, writer, attempts = 6L) {
  ensure_dir(dirname(path))
  last_error <- NULL
  for (attempt in seq_len(attempts)) {
    ok <- tryCatch(
      {
        direct_ok <- tryCatch(
          {
            writer(path)
            TRUE
          },
          error = function(error) {
            last_error <<- error
            FALSE
          }
        )
        if (!direct_ok) {
          extension <- tools::file_ext(path)
          temp_path <- tempfile(
            pattern = paste0(basename(path), "-"),
            fileext = if (nzchar(extension)) paste0(".", extension) else ""
          )
          on.exit(unlink(temp_path), add = TRUE)
          writer(temp_path)
          if (!copy_file(temp_path, path)) {
            stop("file.copy returned FALSE")
          }
        }
        TRUE
      },
      error = function(error) {
        last_error <<- error
        FALSE
      }
    )
    if (ok && file.exists(path)) {
      return(invisible(path))
    }
    Sys.sleep(0.2 * attempt)
  }
  stop(
    "Could not write ", path, ": ", conditionMessage(last_error),
    call. = FALSE
  )
}

as_number <- function(x) {
  as.numeric(x)
}

sanitize_sheet_name <- function(name) {
  name <- gsub("[][*/?:\\\\]", "-", name)
  substr(name, 1L, 31L)
}

write_xlsx_workbook <- function(sheets, path) {
  write_with_retries(path, function(target) {
    workbook <- openxlsx::createWorkbook()
    used_names <- character(0L)
    for (sheet_name in names(sheets)) {
      safe_name <- sanitize_sheet_name(sheet_name)
      if (safe_name %in% used_names) {
        suffix <- length(used_names) + 1L
        safe_name <- substr(
          paste0(substr(safe_name, 1L, 26L), "-", suffix),
          1L,
          31L
        )
      }
      used_names <- c(used_names, safe_name)
      openxlsx::addWorksheet(workbook, safe_name)
      openxlsx::writeData(workbook, safe_name, sheets[[sheet_name]])
      openxlsx::freezePane(workbook, safe_name, firstRow = TRUE)
      openxlsx::setColWidths(workbook, safe_name, cols = 1:50, widths = "auto")
    }
    openxlsx::saveWorkbook(workbook, target, overwrite = TRUE)
  })
}

key_value_frame <- function(values) {
  data.frame(
    metric = names(values),
    value = as.character(unlist(values, use.names = FALSE)),
    stringsAsFactors = FALSE
  )
}

dependency_versions_frame <- function(packages) {
  data.frame(
    package = names(packages),
    version = unlist(packages, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

ansi_color <- function(text, code) {
  if (nzchar(Sys.getenv("NO_COLOR"))) {
    return(text)
  }
  paste0("\033[", code, "m", text, "\033[0m")
}

cat_section <- function(text) {
  cat("\n", ansi_color(paste0("== ", text, " =="), "36;1"), "\n", sep = "")
}

cat_key_value <- function(key, value) {
  cat(ansi_color(sprintf("%-24s", paste0(key, ":")), "33;1"), value, "\n",
      sep = "")
}

format_seconds <- function(seconds) {
  if (is.na(seconds)) {
    return("NA")
  }
  if (seconds < 60) {
    return(sprintf("%.2f sec", seconds))
  }
  sprintf("%.2f min", seconds / 60)
}

format_bytes <- function(bytes) {
  if (is.na(bytes)) {
    return("NA")
  }
  units <- c("B", "KB", "MB", "GB", "TB")
  value <- as.numeric(bytes)
  unit_id <- 1L
  while (abs(value) >= 1024 && unit_id < length(units)) {
    value <- value / 1024
    unit_id <- unit_id + 1L
  }
  sprintf("%.2f %s", value, units[[unit_id]])
}

format_time <- function(time) {
  format(time, "%Y-%m-%d %H:%M:%S %Z")
}

current_process_peak_memory_bytes <- function() {
  if (.Platform$OS.type == "windows") {
    value <- suppressWarnings(system2(
      "powershell.exe",
      c(
        "-NoProfile", "-Command",
        paste0("(Get-Process -Id ", Sys.getpid(), ").PeakWorkingSet64")
      ),
      stdout = TRUE,
      stderr = FALSE
    ))
    value <- suppressWarnings(as.numeric(value[[1L]]))
    if (!is.na(value)) {
      return(value)
    }
  }

  status_file <- "/proc/self/status"
  if (file.exists(status_file)) {
    status <- readLines(status_file, warn = FALSE)
    vm_hwm <- grep("^VmHWM:", status, value = TRUE)
    if (length(vm_hwm) > 0L) {
      kb <- suppressWarnings(as.numeric(gsub("[^0-9]", "", vm_hwm[[1L]])))
      if (!is.na(kb)) {
        return(kb * 1024)
      }
    }
  }

  NA_real_
}

start_log_capture <- function(log_file) {
  ensure_dir(dirname(log_file))
  connection <- file(log_file, open = "wt", encoding = "UTF-8")
  sink(connection, split = TRUE)
  connection
}

stop_log_capture <- function(connection) {
  if (sink.number() > 0L) {
    sink()
  }
  close(connection)
}

extract_gc_count <- function(bench_row) {
  if ("n_gc" %in% names(bench_row)) {
    return(as.integer(bench_row[["n_gc"]][[1L]]))
  }
  if (!"gc" %in% names(bench_row)) {
    return(NA_integer_)
  }
  gc_value <- bench_row[["gc"]][[1L]]
  gc_df <- tryCatch(as.data.frame(gc_value), error = function(e) NULL)
  if (is.null(gc_df)) {
    return(NA_integer_)
  }
  numeric_cols <- vapply(gc_df, is.numeric, logical(1L))
  if (!any(numeric_cols)) {
    return(NA_integer_)
  }
  as.integer(sum(as.matrix(gc_df[numeric_cols]), na.rm = TRUE))
}

fit_summary_row <- function(scenario, fit) {
  optim <- fit@optim
  data.frame(
    scenario = scenario$label,
    converged = isTRUE(optim$converged),
    iterations = if (is.null(optim$iterations)) NA_integer_ else
      as.integer(optim$iterations),
    npar = if (is.null(optim$npar)) NA_integer_ else as.integer(optim$npar),
    fx = if (is.null(optim$fx)) NA_real_ else as.numeric(optim$fx),
    estimator = fit@Options$estimator,
    se = fit@Options$se,
    test = paste(fit@Options$test, collapse = ";"),
    ngroups = fit@Data@ngroups,
    nlevels = fit@Data@nlevels,
    stringsAsFactors = FALSE
  )
}

stage_timing_rows <- function(scenario, fit) {
  timing <- fit@timing
  timing_names <- setdiff(names(timing), "start_time")
  if (length(timing_names) == 0L) {
    return(data.frame(
      scenario = character(0L),
      stage = character(0L),
      elapsed_sec = numeric(0L),
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    scenario = scenario$label,
    stage = timing_names,
    elapsed_sec = vapply(timing[timing_names], as_number, numeric(1L)),
    stringsAsFactors = FALSE
  )
}

benchmark_scenario <- function(scenario, iterations) {
  result <- bench::mark(
    {
      fit <- scenario$fit()
      invisible(isTRUE(fit@optim$converged))
    },
    iterations = iterations,
    check = FALSE,
    memory = TRUE,
    filter_gc = FALSE
  )
  timings <- as_number(result$time[[1L]])

  data.frame(
    scenario = scenario$label,
    iterations = iterations,
    min_sec = as_number(result$min[[1L]]),
    mean_sec = mean(timings),
    median_sec = as_number(result$median[[1L]]),
    itr_per_sec = as_number(result[["itr/sec"]][[1L]]),
    mem_alloc_bytes = as_number(result$mem_alloc[[1L]]),
    gc_count = extract_gc_count(result[1L, , drop = FALSE]),
    stringsAsFactors = FALSE
  )
}

hotspot_rows_from_summary <- function(scenario, summary, source_file,
                                      kind = c("by.self", "by.total"),
                                      top_n = 25L) {
  kind <- match.arg(kind)
  table <- summary[[kind]]
  if (is.null(table) || nrow(table) == 0L) {
    return(data.frame(
      scenario = character(0L),
      rank = integer(0L),
      kind = character(0L),
      function_name = character(0L),
      self_time_sec = numeric(0L),
      self_pct = numeric(0L),
      total_time_sec = numeric(0L),
      total_pct = numeric(0L),
      source_file = character(0L),
      stringsAsFactors = FALSE
    ))
  }
  n <- min(top_n, nrow(table))
  table <- table[seq_len(n), , drop = FALSE]
  data.frame(
    scenario = scenario$label,
    rank = seq_len(n),
    kind = kind,
    function_name = rownames(table),
    self_time_sec = table[, "self.time"],
    self_pct = table[, "self.pct"],
    total_time_sec = table[, "total.time"],
    total_pct = table[, "total.pct"],
    source_file = source_file,
    stringsAsFactors = FALSE
  )
}

profile_scenario_rprof <- function(scenario, top_n = 25L) {
  rprof_file <- tempfile(paste0(scenario$label, "-"), fileext = ".rprof")
  on.exit(unlink(rprof_file, force = TRUE), add = TRUE)
  utils::Rprof(rprof_file, memory.profiling = TRUE)
  on.exit(utils::Rprof(NULL), add = TRUE)
  fit <- scenario$fit()
  invisible(isTRUE(fit@optim$converged))
  utils::Rprof(NULL)

  summary <- utils::summaryRprof(rprof_file, memory = "both")
  rbind(
    hotspot_rows_from_summary(scenario, summary, NA_character_, "by.self",
                              top_n),
    hotspot_rows_from_summary(scenario, summary, NA_character_, "by.total",
                              top_n)
  )
}

profile_scenario_memory <- function(scenario) {
  result <- profmem::profmem({
    fit <- scenario$fit()
    invisible(isTRUE(fit@optim$converged))
  })
  bytes <- result$bytes[!is.na(result$bytes)]
  data.frame(
    scenario = scenario$label,
    allocations = length(bytes),
    total_alloc_bytes = sum(bytes),
    largest_alloc_bytes = if (length(bytes) == 0L) NA_real_ else max(bytes),
    stringsAsFactors = FALSE
  )
}

select_scenarios <- function(scenarios, filter) {
  if (identical(filter, "all")) {
    return(scenarios)
  }
  labels <- vapply(scenarios, `[[`, character(1L), "label")
  missing <- setdiff(filter, labels)
  if (length(missing) > 0L) {
    stop("Unknown scenario label(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }
  scenarios[match(filter, labels)]
}

package_version_list <- function(packages) {
  versions <- vapply(packages, function(package) {
    as.character(utils::packageVersion(package))
  }, character(1L))
  as.list(versions)
}

remove_stale_directory <- function(output_dir, name) {
  target <- file.path(output_dir, name)
  if (!dir.exists(target)) {
    return(invisible(FALSE))
  }
  if (!path_is_within(target, output_dir)) {
    stop("Refusing to remove directory outside run output: ", target,
         call. = FALSE)
  }
  unlink(target, recursive = TRUE, force = TRUE)
  invisible(TRUE)
}

remove_stale_files <- function(output_dir, patterns) {
  files <- unlist(lapply(patterns, function(pattern) {
    list.files(output_dir, pattern = pattern, full.names = TRUE,
               recursive = FALSE)
  }), use.names = FALSE)
  if (length(files) == 0L) {
    return(invisible(FALSE))
  }
  for (file in files) {
    if (!path_is_within(file, output_dir)) {
      stop("Refusing to remove file outside run output: ", file,
           call. = FALSE)
    }
  }
  unlink(files, force = TRUE)
  invisible(TRUE)
}

scenario_latency_summary_frame <- function(benchmark_summary) {
  data.frame(
    scenario = benchmark_summary$scenario,
    iterations = benchmark_summary$iterations,
    mean_sec = benchmark_summary$mean_sec,
    median_sec = benchmark_summary$median_sec,
    min_sec = benchmark_summary$min_sec,
    itr_per_sec = benchmark_summary$itr_per_sec,
    mem_alloc_bytes = benchmark_summary$mem_alloc_bytes,
    gc_count = benchmark_summary$gc_count,
    stringsAsFactors = FALSE
  )
}

print_latency_table <- function(latency_summary) {
  display <- data.frame(
    scenario = latency_summary$scenario,
    iterations = latency_summary$iterations,
    mean_ms = sprintf("%.2f", latency_summary$mean_sec * 1000),
    median_ms = sprintf("%.2f", latency_summary$median_sec * 1000),
    mem_alloc = vapply(
      latency_summary$mem_alloc_bytes,
      format_bytes,
      character(1L)
    ),
    gc_count = latency_summary$gc_count,
    stringsAsFactors = FALSE
  )
  print(display, row.names = FALSE, right = FALSE)
}

run_lavaan_profiling <- function(cli_args = commandArgs(trailingOnly = TRUE)) {
  run_start_time <- Sys.time()
  run_start_elapsed <- proc.time()[["elapsed"]]
  project_root <- find_project_root()
  args <- parse_args(cli_args, project_root)
  setwd(project_root)
  source(file.path(project_root, "benchmarks", "profiling", "scenarios.R"))

  all_scenarios <- lavaan_profiling_scenarios()
  scenarios <- select_scenarios(all_scenarios, args$scenario_filter)
  ensure_dir(args$output_dir)
  remove_stale_directory(args$output_dir, "profvis")
  remove_stale_directory(args$output_dir, "rprof")
  remove_stale_files(args$output_dir, c("\\.csv$", "\\.json$"))
  log_connection <- start_log_capture(args$log_file)
  on.exit(stop_log_capture(log_connection), add = TRUE)
  check_required_packages()
  pkgload::load_all(project_root, quiet = TRUE)
  setwd(project_root)

  cat_section("lavaan profiling run")
  cat_key_value("run label", args$run_label)
  cat_key_value("mode", args$mode)
  cat_key_value("iterations", args$iterations)
  cat_key_value("start time", format_time(run_start_time))
  cat_key_value("workbook", args$workbook_path_abs)
  cat_key_value("log file", args$log_file_abs)
  cat_key_value(
    "scenarios",
    paste(vapply(scenarios, `[[`, character(1L), "label"), collapse = ", ")
  )

  metadata <- list(
    run_label = args$run_label,
    mode = args$mode,
    iterations = args$iterations,
    output_dir = args$output_dir_abs,
    workbook_path = args$workbook_path_abs,
    log_file = args$log_file_abs,
    project_root = project_root,
    start_time = format(run_start_time, "%Y-%m-%dT%H:%M:%S%z"),
    r_version = as.character(getRversion()),
    platform = R.version$platform,
    lavaan_version = as.character(utils::packageVersion("lavaan")),
    git_head = tryCatch(
      system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
      error = function(e) NA_character_
    ),
    git_branch = tryCatch(
      system2("git", c("branch", "--show-current"), stdout = TRUE,
              stderr = FALSE),
      error = function(e) NA_character_
    ),
    dependency_versions = package_version_list(required_packages)
  )

  fit_rows <- list()
  timing_rows <- list()
  bench_rows <- list()
  hotspot_rows <- list()
  profmem_rows <- list()

  for (i in seq_along(scenarios)) {
    scenario <- scenarios[[i]]
    cat(
      "\n",
      ansi_color(
        paste0("[", i, "/", length(scenarios), "] ", scenario$label),
        "35;1"
      ),
      "\n",
      sep = ""
    )

    fit <- scenario$fit()
    fit_rows[[i]] <- fit_summary_row(scenario, fit)
    timing_rows[[i]] <- stage_timing_rows(scenario, fit)

    bench_rows[[i]] <- benchmark_scenario(scenario, args$iterations)
    profmem_rows[[i]] <- profile_scenario_memory(scenario)
    hotspot_rows[[i]] <- profile_scenario_rprof(scenario)
  }

  fit_summary <- do.call(rbind, fit_rows)
  stage_timing_summary <- do.call(rbind, timing_rows)
  benchmark_summary <- do.call(rbind, bench_rows)
  hotspot_summary <- do.call(rbind, hotspot_rows)
  profmem_summary <- do.call(rbind, profmem_rows)
  latency_summary <- scenario_latency_summary_frame(benchmark_summary)
  run_end_time <- Sys.time()
  total_run_time_sec <- proc.time()[["elapsed"]] - run_start_elapsed
  peak_process_memory_bytes <- current_process_peak_memory_bytes()
  max_scenario_alloc_bytes <- max(benchmark_summary$mem_alloc_bytes,
                                  na.rm = TRUE)
  max_scenario_alloc <- benchmark_summary[
    which.max(benchmark_summary$mem_alloc_bytes),
    c("scenario", "mem_alloc_bytes"),
    drop = FALSE
  ]

  metadata$end_time <- format(run_end_time, "%Y-%m-%dT%H:%M:%S%z")
  metadata$total_run_time_sec <- total_run_time_sec
  metadata$peak_process_memory_bytes <- peak_process_memory_bytes
  metadata$max_scenario_alloc_bytes <- max_scenario_alloc_bytes

  run_summary <- key_value_frame(c(
    run_label = args$run_label,
    mode = args$mode,
    iterations = args$iterations,
    scenario_count = length(scenarios),
    start_time = format_time(run_start_time),
    end_time = format_time(run_end_time),
    total_run_time_sec = sprintf("%.3f", total_run_time_sec),
    total_run_time = format_seconds(total_run_time_sec),
    peak_process_memory_bytes = peak_process_memory_bytes,
    peak_process_memory = format_bytes(peak_process_memory_bytes),
    max_scenario_alloc_bytes = max_scenario_alloc_bytes,
    max_scenario_alloc = format_bytes(max_scenario_alloc_bytes),
    max_scenario_alloc_scenario = max_scenario_alloc$scenario,
    workbook_path = args$workbook_path_abs,
    log_file = args$log_file_abs,
    project_root = project_root,
    git_head = metadata$git_head,
    git_branch = metadata$git_branch
  ))

  write_xlsx_workbook(
    list(
      "run-summary" = run_summary,
      "dependency-versions" =
        dependency_versions_frame(metadata$dependency_versions),
      "scenario-metadata" = profiling_scenario_metadata(scenarios),
      "scenario-latency" = latency_summary,
      "benchmark-summary" = benchmark_summary,
      "fit-summary" = fit_summary,
      "profmem-summary" = profmem_summary,
      "stage-timing-summary" = stage_timing_summary,
      "hotspot-summary" = hotspot_summary
    ),
    args$workbook_path
  )

  cat_section("profiling complete")
  cat_key_value("workbook", args$workbook_path_abs)
  cat_key_value("log file", args$log_file_abs)
  cat_key_value("start time", format_time(run_start_time))
  cat_key_value("end time", format_time(run_end_time))
  cat_key_value("total run time", format_seconds(total_run_time_sec))
  cat_key_value("peak process memory", format_bytes(peak_process_memory_bytes))
  cat_key_value(
    "max scenario alloc",
    paste0(
      format_bytes(max_scenario_alloc_bytes),
      " (",
      max_scenario_alloc$scenario,
      ")"
    )
  )

  cat_section("scenario latency summary")
  print_latency_table(latency_summary)

  invisible(list(
    output_dir = args$output_dir_abs,
    workbook_path = args$workbook_path_abs,
    log_file = args$log_file_abs,
    scenarios = vapply(scenarios, `[[`, character(1L), "label")
  ))
}

run_lavaan_profiling_cli <- function() {
  run_lavaan_profiling(commandArgs(trailingOnly = TRUE))
}

if (sys.nframe() == 0L) {
  run_lavaan_profiling_cli()
}
