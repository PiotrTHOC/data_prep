---
name: r-statistical-ml
description: Best practices for statistical analysis and machine learning in R — feature preparation, model building, evaluation, and biological interpretation. Use when preparing data for ML, building models, or performing statistical tests on genomic data.
user-invocable: false
---

# Statistical Analysis & Machine Learning in R

## Feature Matrix Preparation

- Final ML matrix structure: rows = observations (gene x sample), columns = features + label
- Required preprocessing pipeline:
  1. Filter zero-variance features: `nearZeroVar()` from caret
  2. Log-transform skewed features: `log2(x + 1)`
  3. Scale and center: `scale()` — compute on training set, apply to test set
  4. Handle multicollinearity: check `cor()` matrix, remove features with r > 0.95
- Label encoding:
  ```r
  df$label <- factor(df$label, levels = c("normal", "cancer"))
  ```
  - Set reference level explicitly — the first level is the negative class in most R packages

## Train/Test Split

- Use `caret::createDataPartition()` for stratified splits preserving class balance
- For genomic data, be aware of data leakage:
  - Same gene in different samples is NOT independent — consider gene-level splitting
  - Samples from the same tissue/cell line are correlated — consider sample-level CV
- Typical split: 70/30 or 80/20 for train/test
- For small datasets (< 500 samples), use repeated k-fold CV instead of a hold-out set

## Model Selection for Genomic Data

- Start simple: logistic regression or elastic net (`glmnet`) — interpretable, handles high-dimensional data
- Random forest (`ranger` or `randomForest`): good baseline, handles non-linear relationships, gives feature importance
- XGBoost (`xgboost`): best performance for tabular data, requires careful tuning
- SVM (`e1071::svm`): effective for high-dimensional, small-sample genomic data
- Avoid deep learning unless n > 10,000 and you have a clear architectural hypothesis

## Cross-Validation

- Use `caret::trainControl()` with `method = "repeatedcv"`, `number = 10`, `repeats = 3`
- For imbalanced classes: use `classProbs = TRUE`, `summaryFunction = twoClassSummary`
- Report AUC-ROC as primary metric for binary classification — accuracy is misleading with class imbalance
- Include confidence intervals: bootstrap or repeated CV standard deviation

## Feature Importance & Interpretation

- Extract variable importance from models: `caret::varImp()` or `randomForest::importance()`
- For genomic ML, always map feature importance back to biology:
  - Which chromatin marks are most predictive?
  - Are promoter or distal features more informative?
  - Does accessibility (ATAC) add value over histone marks alone?
- Use SHAP values (`shapviz` package) for model-agnostic interpretation
- Visualize top features with `ggplot2` — barplots for importance, boxplots for distributions

## Statistical Testing

- For differential expression: use DESeq2 or edgeR, never t-tests on raw counts
- For enrichment analysis: `clusterProfiler::enrichGO()`, `enrichKEGG()`
- Adjust for multiple testing: `p.adjust(method = "BH")` (Benjamini-Hochberg)
- Report effect sizes alongside p-values — statistical significance != biological significance
- For comparing distributions: Wilcoxon rank-sum (non-parametric) over t-test when normality is questionable

## Reproducibility in ML

- Set seeds before any random operation: `set.seed(42)`
- Record the full preprocessing pipeline — transformations must be reproducible on new data
- Save trained models: `saveRDS(model, "model_v1.rds")`
- Log hyperparameters, performance metrics, and data versions
- Use `sessionInfo()` to record package versions
