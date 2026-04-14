---
name: r-reproducibility
description: Best practices for reproducible scientific computing in R — environment management, project structure, and provenance tracking. Use when setting up projects, managing dependencies, or ensuring reproducibility.
user-invocable: false
---

# Reproducible Scientific Computing in R

## Environment Management

- Use `renv` for dependency management — never install packages inside analysis scripts:
  ```r
  renv::init()       # initialize project
  renv::snapshot()   # lock current package versions
  renv::restore()    # recreate environment from lockfile
  ```
- Pin Bioconductor version alongside R version: they must be compatible
- Record `sessionInfo()` or `renv::diagnostics()` at the end of every analysis
- Never use `install.packages()` or `BiocManager::install()` inside analysis scripts — separate setup from analysis

## Project Structure

```
project/
  R/                  # reusable functions
  scripts/            # analysis scripts (numbered: 01_prep.R, 02_annotate.R)
  data/               # raw input data (read-only, never modified)
  output/             # generated outputs (figures, tables, models)
  renv/               # dependency lockfile and library
  renv.lock           # version-locked dependencies
  .Rprofile           # renv activation
  DESCRIPTION         # project metadata (optional but recommended)
```

- Keep raw data read-only — never overwrite input files
- Number scripts in execution order: `01_`, `02_`, etc.
- Separate reusable functions (in `R/`) from workflow scripts (in `scripts/`)

## File Paths & Portability

- Never use `setwd()` — it breaks reproducibility across machines
- Use `here::here()` for project-relative paths:
  ```r
  peak_file <- here::here("data", "tissue37_H3K27ac.bed.gz")
  ```
- Never use absolute paths like `"D:/chromatin_data"` or `"/home/user/project"`
- Use `fs::path()` for safe path construction across OS

## Parameterization

- Define all analysis parameters as named constants at the top of the script:
  ```r
  TSS_UPSTREAM   <- 1000
  TSS_DOWNSTREAM <- 100
  QVALUE_CUTOFF  <- 0.05
  GENOME         <- "hg38"
  ```
- Never hardcode thresholds inside functions — pass them as arguments
- Document the biological rationale for parameter choices in comments

## Provenance Tracking

- Log which input files were used and their checksums:
  ```r
  file_info <- data.frame(
    file = input_files,
    md5 = tools::md5sum(input_files),
    modified = file.mtime(input_files)
  )
  ```
- Save analysis metadata alongside results:
  ```r
  attr(result, "analysis_date") <- Sys.time()
  attr(result, "r_version") <- R.version.string
  attr(result, "params") <- list(tss_region = c(-1000, 100), genome = "hg38")
  ```
- Use git for version control of scripts — never version control large data files

## Workflow Orchestration

- For multi-step pipelines, use `targets` package (successor to `drake`):
  - Automatic dependency tracking between pipeline steps
  - Only re-runs steps whose inputs changed
  - Built-in caching and provenance
- For simpler workflows, use a Makefile or a numbered script runner
- Never rely on manual execution order — document and automate the full pipeline

## Seed Management

- Set `set.seed()` before every stochastic operation, not just once at the top
- Use different but documented seeds for different random operations
- Record the RNG kind: `RNGkind()` — defaults vary across R versions
