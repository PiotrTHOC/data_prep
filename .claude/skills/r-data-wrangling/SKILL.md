---
name: r-data-wrangling
description: Best practices for data manipulation in R using tidyverse and data.table — joins, reshaping, aggregation, and feature engineering for scientific datasets. Use when transforming, merging, or reshaping data.
user-invocable: false
---

# R Data Wrangling — Scientific Computing

## Tidyverse Pipeline Design

- One verb per line in pipelines for readability:
  ```r
  result <- df %>%
    dplyr::filter(qvalue < 0.05) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(
      signal_sum = sum(signal_value, na.rm = TRUE),
      peak_count = dplyr::n(),
      .groups = "drop"
    )
  ```
- Always pass `.groups = "drop"` to `summarise()` to avoid unexpected grouping
- Use explicit `dplyr::` namespace to avoid conflicts with `stats::filter()`, `stats::lag()`
- Prefer `dplyr::across()` for applying the same operation to multiple columns

## Joins

- Choose the correct join for the analytical goal:
  - `inner_join()`: only genes present in ALL datasets — use for ML feature matrices
  - `full_join()`: retain all genes — use when NA = "no signal" is informative
  - `left_join()`: keep all from primary dataset — use when one dataset is the reference
- Always specify the `by` argument explicitly — never rely on auto-detection
- After joining, check for unexpected row count changes (many-to-many joins)
- When merging iteratively, use `purrr::reduce(list_of_dfs, dplyr::full_join, by = "key")`

## Missing Data

- Distinguish between structural zeros (no peak = no signal) and missing data (experiment not performed)
- Replace NA with 0 only when biologically justified (absent ChIP peaks = zero signal)
- Use `tidyr::replace_na()` or direct assignment for explicit NA handling
- Document the NA replacement strategy in comments
- Check for all-zero rows after NA replacement — these genes carry no information

## Feature Engineering

- Log-transform right-skewed distributions (expression, signal values): `log2(x + 1)`
  - Add pseudocount of 1 to handle zeros
  - Use log2 (not log10 or ln) — standard in genomics, interpretable as fold-changes
- Scale/center features for ML with `scale()` — apply to numeric columns only:
  ```r
  numeric_cols <- sapply(df, is.numeric)
  df[, numeric_cols] <- scale(df[, numeric_cols])
  ```
- Create ratio features when biologically meaningful: `H3K27ac_promoter / (H3K4me3_promoter + 1)`
- Binarize labels explicitly: use factors with defined levels for classification

## Data Validation

- After every major transformation, verify dimensions: `dim()`, `str()`, `summary()`
- Check for duplicated keys before and after joins: `sum(duplicated(df$gene))`
- Verify no unexpected NAs were introduced: `colSums(is.na(df))`
- Check value ranges make biological sense (TPM >= 0, signal >= 0)
- Use `stopifnot()` or `assertthat::assert_that()` for critical invariants

## I/O Best Practices

- Read TSV/CSV with explicit types: `read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)`
- For large files, use `data.table::fread()` or `readr::read_tsv()` — 10-100x faster
- Write outputs with `write.csv(row.names = FALSE)` or `readr::write_csv()`
- Save intermediate R objects as `.rds` (single object) not `.RData` (entire environment)
- Include metadata (date, parameters, source files) in output filenames or as attributes

## Reshaping

- Use `tidyr::pivot_longer()` and `tidyr::pivot_wider()` — never `reshape()` or `reshape2`
- When converting wide feature matrices to long format for plotting:
  ```r
  df_long <- tidyr::pivot_longer(df, cols = -gene, names_to = "feature", values_to = "value")
  ```
- Keep data in wide format for ML model input; long format for ggplot2 visualization
