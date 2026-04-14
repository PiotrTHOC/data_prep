---
name: r-code-review
description: Review R scripts for correctness, performance, and scientific rigor. Use when asked to review, audit, or improve R code for data analysis.
---

# R Code Review — Scientific Data Analysis

When reviewing R code, check these areas systematically:

## 1. Correctness

- Verify join keys match across datasets (ENSEMBL vs Entrez vs Symbol)
- Check that NA handling is biologically appropriate (NA -> 0 vs. filtering)
- Verify aggregation logic: `sum()` vs `mean()` chosen for the right reasons
- Check for off-by-one errors in genomic coordinates (0-based BED vs 1-based GRanges)
- Ensure `group_by()` variables are correct and `.groups = "drop"` is set
- Watch for accidental column name collisions after joins

## 2. Data Integrity

- Are duplicated rows/genes detected and handled?
- Is gene ID versioning consistent (ENSG00000141510 vs ENSG00000141510.18)?
- Are chromosome naming styles consistent (chr1 vs 1)?
- Are genome assemblies consistent across all inputs (hg38 vs hg19)?
- Is the join type appropriate? `inner_join` drops unmatched rows silently

## 3. Performance

- Flag unnecessary repeated operations (rebuilding TxDb in loops)
- Identify row-wise operations that should be vectorized
- Check for `sapply()` calls that should be `vapply()` for type safety
- Flag `rbind()` in loops — use `dplyr::bind_rows()` or pre-allocate
- Look for unnecessary copies of large objects

## 4. Reproducibility

- Are random seeds set before stochastic operations?
- Are file paths relative and configurable, not hardcoded?
- Is `setwd()` avoided? Use project-relative paths instead
- Are package versions recorded?
- Are parameters (thresholds, window sizes) defined as named constants?

## 5. Scientific Rigor

- Is the normalization method appropriate for the data type?
- Are batch effects considered when merging samples?
- Is the experimental design reflected in the analysis (paired samples, replicates)?
- Are QC steps included (peak count checks, signal distributions)?
- Are biological controls used appropriately?

## Output Format

Structure your review as:
1. **Critical issues** — bugs or errors that produce wrong results
2. **Improvements** — changes that improve correctness or performance
3. **Suggestions** — optional enhancements for maintainability

For each finding, include the line number, the issue, and a concrete fix.
