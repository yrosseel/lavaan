# lavaan profiling suite

This directory contains the tracked profiling harness for evidence-first
lavaan performance work. It deliberately stores generated profiling artifacts
outside the package directory so that benchmark snapshots do not become part of
the package source tree.

## What it measures

The suite starts from public, user-facing workflows and then records enough
detail to justify later hotspot work:

- repeated runtime and allocation measurements from `bench::mark`
- per-stage lavaan timing from the fitted object's `@timing` slot
- convergence, optimizer iteration, estimator, group, and level metadata
- call-stack hotspots from `Rprof`/`summaryRprof`
- a single workbook with all run, scenario, benchmark, memory, timing, and
  hotspot summaries

The current scenario matrix covers CFA, SEM, growth, two-level SEM, ordinal
WLSMV, FIML missing data, multigroup equality constraints, sample-covariance
input, and a deterministic larger synthetic CFA.

## Run it

Run commands from the `lavaan/` package root.

Quick smoke run:

```sh
Rscript benchmarks/profiling/run-profiling.R smoke quick 1
```

Full current-branch baseline:

```sh
Rscript benchmarks/profiling/run-profiling.R current-branch full
```

By default, outputs are written to a sibling directory:

```text
../profiling-output/<run-label>/profiling-summary.xlsx
../logs/<run-label>.log
```

The runner refuses to write artifacts inside the `lavaan/` directory. To use a
different output root, pass it as the fourth argument:

```sh
Rscript benchmarks/profiling/run-profiling.R current-branch full 900 ../profiling-output
```

To run only selected scenarios, pass a comma-separated label list as the fifth
argument:

```sh
Rscript benchmarks/profiling/run-profiling.R smoke quick 1 ../profiling-output cfa_holzinger_ml,sem_political_democracy_ml
```

To write logs somewhere other than `../logs`, pass the log root as the sixth
argument:

```sh
Rscript benchmarks/profiling/run-profiling.R smoke quick 1 ../profiling-output all ../logs
```

Arguments are:

```text
1. run_label      default: current-branch
2. mode           quick or full; default: quick
3. iterations     bench::mark iterations; default: 1 for quick, 900 for full
4. output_root    default: ../profiling-output
5. scenarios      comma-separated labels or all; default: all
6. log_root       default: ../logs
```

Use `full` mode for branch baselines. It is intentionally sized for a
15-30 minute run on a typical development machine; pass a smaller explicit
iteration count only when you are debugging the harness itself.

## Outputs

Each run directory contains one workbook:

- `profiling-summary.xlsx`

The workbook includes these sheets:

- `run-summary`: start/end time, total run time, peak memory, paths, Git metadata
- `dependency-versions`: package versions used by the harness
- `scenario-metadata`: scenario labels and descriptions
- `scenario-latency`: mean and median latency per scenario
- `benchmark-summary`: repeated timing, allocation, and GC counts
- `fit-summary`: convergence, iterations, estimator, groups, and levels
- `profmem-summary`: allocation counts and total allocated bytes
- `stage-timing-summary`: lavaan `@timing` stage durations
- `hotspot-summary`: top `Rprof` functions by self and total time

The raw `Rprof` trace is temporary and deleted after summary extraction. No raw
`.out` files, `rprof/` directory, or HTML profiler directory are kept.

Each run also writes a console transcript to `../logs/<run-label>.log`.


## Evidence standard

Do not call a function a bottleneck from one number alone. Treat something as a
candidate optimization target only when repeated benchmark timing/allocation
evidence and profiler call-stack evidence point to the same workflow or
function. The workbook is the durable baseline; generated evidence stays
outside `lavaan/`.
