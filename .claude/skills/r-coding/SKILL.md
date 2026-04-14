---
name: r-coding
description: Senior R developer conventions, style guide, and best practices for scientific computing. Use when writing, reviewing, or refactoring any R code.
user-invocable: false
---

# R Coding Standards — Senior Developer Level

## Style & Naming

- Follow the tidyverse style guide: snake_case for variables and functions, PascalCase only for S4/R6 classes
- Maximum line length: 80 characters. Break long pipelines with one verb per line
- Use `<-` for assignment, never `=` at the top level
- Prefix Boolean variables with `is_`, `has_`, `should_`
- Name functions as verbs: `read_peaks()`, `annotate_regions()`, `compute_signal()`
- Use meaningful, domain-specific names: `peak_granges` not `x`, `signal_matrix` not `mat`

## Function Design

- Every function should do one thing. If a function exceeds 30 lines, consider splitting
- Always validate inputs at function boundaries with informative error messages:
  ```r
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  if (!is(granges, "GRanges")) stop("Expected GRanges, got ", class(granges)[1], call. = FALSE)
  ```
- Use explicit `return()` only for early returns; rely on implicit return for the final value
- Document non-obvious parameters and return types with roxygen2 `@param` and `@return`
- Avoid side effects: functions should not modify global state, call `setwd()`, or install packages
- Use `...` (dots) sparingly and only when genuinely forwarding arguments

## Vectorization & Performance

- Prefer vectorized operations over loops. Use `vapply()` over `sapply()` for type safety
- For grouped operations, use `dplyr::summarise()` or `data.table` — never row-wise loops
- Pre-allocate vectors when loops are unavoidable: `result <- vector("list", n)`
- Profile before optimizing: use `bench::mark()` or `system.time()`
- For large genomic data, prefer `data.table` or Bioconductor structures over base data.frames

## Error Handling

- Use `tryCatch()` / `withCallingHandlers()` for operations that may fail (file I/O, network)
- Use `cli::cli_abort()` or `rlang::abort()` over base `stop()` in package code
- Never use `suppressWarnings()` / `suppressMessages()` without documenting why
- Wrap external tool calls (system commands, API requests) in error handlers

## Dependencies & Namespace

- Always use explicit namespacing for non-base functions: `dplyr::filter()`, `purrr::map()`
- Avoid `library()` calls inside functions; use `requireNamespace()` for soft dependencies
- Group `library()` calls at the top of scripts
- Never use `require()` — use `library()` (fails loudly) or `requireNamespace()` (returns logical)

## Code Organization

- Structure scripts in sections: setup, data loading, processing, output
- Use roxygen2 section headers `# Section Name ----` for RStudio navigation
- Keep configuration (file paths, parameters) at the top as named constants
- Separate reusable logic into sourced utility files or packages

## Memory Management

- Remove large intermediate objects with `rm()` and call `gc()` after processing large datasets
- Use streaming/chunked reads for files > 1 GB
- Prefer in-place operations where possible (data.table `:=`)

## Output & Logging

- Use `message()` for progress info, `warning()` for recoverable issues, `stop()` for fatal errors
- Never use `print()` or `cat()` for logging — these go to stdout and break pipelines
- Include the object/file/step name in all messages for traceability
