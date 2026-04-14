# load libraries
library(here)
library(ChIPseeker)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(txdbmaker)
library(GenomeInfoDb)
library(dplyr)
library(purrr)
library(org.Hs.eg.db)

# TxDb
txdb <- txdbmaker::makeTxDbFromGFF(here::here("input", "GENCODE_v32.gtf"))
seqlevelsStyle(txdb) <- "UCSC"
genome(txdb) <- "hg38"

# Build transcript (ENST) → gene (ENSG) mapping via org.Hs.eg.db
# The GTF uses ENST IDs for both gene_id and transcript_id (UCSC knownGene format)
all_enst <- keys(org.Hs.eg.db, keytype = "ENSEMBLTRANS")
tx2gene_clean <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = all_enst,
  keytype = "ENSEMBLTRANS",
  column = "ENSEMBL",
  multiVals = "first"
)

# Core function: annotate peaks and summarize per gene
annotate_and_summarize <- function(file, mark_name, txdb, tx2gene_clean) {

  if (!file.exists(file)) {
    stop(paste("Missing file:", file))
  }

  peak <- readPeakFile(file)

  # normalize signal column
  peak$signalValue <- peak$V7
  peak$pvalue <- peak$V8
  peak$qvalue <- peak$V9
  genome(peak) <- "hg38"
  seqlevelsStyle(peak) <- "UCSC"

  peakAnno <- annotatePeak(
    peak,
    TxDb = txdb,
    tssRegion = c(-1000, 100)
  )

  df <- as.data.frame(peakAnno)

  features <- df %>%
    dplyr::group_by(geneId) %>%
    dplyr::summarise(
      promoter = sum(signalValue[grepl("Promoter", annotation)], na.rm = TRUE),
      distal = sum(signalValue[annotation == "Distal Intergenic"], na.rm = TRUE),
      genebody = sum(signalValue[grepl("Exon|Intron", annotation)], na.rm = TRUE),
      peak_count = dplyr::n(),
      .groups = "drop"
    )

  # Map transcript IDs (ENST) to gene IDs (ENSG) using GTF-derived mapping
  colnames(features)[1] <- "tx_id"
  features$tx_id <- sub("\\..*", "", features$tx_id)
  features$gene <- tx2gene_clean[features$tx_id]
  features <- features[!is.na(features$gene) & grepl("^ENSG", features$gene), ]
  features$tx_id <- NULL

  # Aggregate to gene level (multiple transcripts may map to same gene)
  features <- features %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), sum), .groups = "drop")

  # prefix columns with mark name
  colnames(features)[colnames(features) != "gene"] <-
    paste0(mark_name, "_", colnames(features)[colnames(features) != "gene"])

  return(features)
}

# Sample definitions
samples <- list(
  tissue37 = list(
    peaks = c(
      "tissue37_H3K27ac" = here::here("input", "tissue37_H3K27ac.bed.gz"),
      "tissue37_H3K4me3" = here::here("input", "tissue37_H3K4me3.bed.gz"),
      "tissue37_CTCF"    = here::here("input", "tissue37_CTCF.bed.gz")
    ),
    expr = here::here("input", "tissue37_RNAseq.tsv")
  ),
  tissue54 = list(
    peaks = c(
      "tissue54_H3K27ac" = here::here("input", "tissue54_H3K27ac.bed.gz"),
      "tissue54_H3K4me3" = here::here("input", "tissue54_H3K4me3.bed.gz"),
      "tissue54_CTCF"    = here::here("input", "tissue54_CTCF.bed.gz"),
      "tissue54_ATAC"    = here::here("input", "tissue54_ATAC.bed.gz")
    ),
    expr = here::here("input", "tissue54_RNAseq.tsv")
  ),
  PC3_1 = list(
    peaks = c(
      "PC3_H3K27ac" = here::here("input", "PC3_H3K27ac.bed.gz"),
      "PC3_H3K4me3" = here::here("input", "PC3_H3K4me3.bed.gz"),
      "PC3_CTCF"    = here::here("input", "PC3_CTCF.bed.gz"),
      "PC3_ATAC"    = here::here("input", "PC3_ATAC.bed.gz")
    ),
    expr = here::here("input", "PC3_RNAseq1.tsv")
  ),
  PC3_2 = list(
    peaks = c(
      "PC3_H3K27ac" = here::here("input", "PC3_H3K27ac.bed.gz"),
      "PC3_H3K4me3" = here::here("input", "PC3_H3K4me3.bed.gz"),
      "PC3_CTCF"    = here::here("input", "PC3_CTCF.bed.gz"),
      "PC3_ATAC"    = here::here("input", "PC3_ATAC.bed.gz")
    ),
    expr = here::here("input", "PC3_RNAseq2.tsv")
  )
)

# Chromatin marks per sample
process_sample_chromatin <- function(peak_files, txdb, tx2gene_clean) {

  feature_list <- lapply(names(peak_files), function(mark) {
    annotate_and_summarize(peak_files[[mark]], mark, txdb, tx2gene_clean)
  })

  chromatin <- purrr::reduce(feature_list, dplyr::full_join, by = "gene")
  chromatin[is.na(chromatin)] <- 0

  return(chromatin)
}

# Read expression data
read_expr <- function(f) {
  message("Reading: ", f)

  expr_df <- read.table(f, header = TRUE, sep = "\t")

  # keep only ENSEMBL genes
  expr_df <- expr_df[grepl("^ENSG", expr_df$gene_id), ]

  # remove version numbers
  expr_df$gene <- sub("\\..*", "", expr_df$gene_id)

  # keep relevant columns
  expr_df <- expr_df[, c("gene", "TPM")]

  return(expr_df)
}

# Merge per sample: chromatin + expression
process_sample <- function(sample_name, sample_info, txdb, tx2gene_clean) {

  message("Processing: ", sample_name)

  chromatin <- process_sample_chromatin(sample_info$peaks, txdb, tx2gene_clean)
  expr <- read_expr(sample_info$expr)

  merged <- dplyr::inner_join(chromatin, expr, by = "gene")

  merged$sample <- sample_name

  return(merged)
}

# Build full dataset
ml_list <- lapply(names(samples), function(s) {
  process_sample(s, samples[[s]], txdb, tx2gene_clean)
})

ml_matrix <- dplyr::bind_rows(ml_list)

# Add labels
ml_matrix$label <- ifelse(
  grepl("PC3", ml_matrix$sample),
  "cancer",
  "normal"
)

# Log-transform expression to stabilize variance
ml_matrix$TPM <- log2(ml_matrix$TPM + 1)

# Scale/center all numeric features
numeric_cols <- sapply(ml_matrix, is.numeric)
ml_matrix[, numeric_cols] <- scale(ml_matrix[, numeric_cols])

cat("Dimensions:", nrow(ml_matrix), "rows x", ncol(ml_matrix), "columns\n")
cat("Samples:\n")
print(table(ml_matrix$sample))
cat("\nLabels:\n")
print(table(ml_matrix$label))
cat("\nFeature summary:\n")
summary(ml_matrix)

write.csv(ml_matrix, here::here("ml_matrix.csv"), row.names = FALSE)
