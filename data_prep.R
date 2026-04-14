# load libraries
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install("ChIPseeker")
BiocManager::install(c("rlang", "cli", "digest", "htmltools", "xfun", "yaml"), 
                     force = TRUE, ask = FALSE)
BiocManager::install("ChIPseeker", force = TRUE, ask = FALSE)

options(repos = BiocManager::repositories())

BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("txdbmaker")

library(ChIPseeker)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(txdbmaker)
library(GenomeInfoDb)
library(dplyr)
library(purrr)
library(org.Hs.eg.db)

setwd("D:/chromatin_data")

# Files
peak_files <- list(
  tissue37_H3K27ac = "tissue37_H3K27ac.bed.gz",
  tissue37_H3K4me3 = "tissue37_H3K4me3.bed.gz",
  tissue37_CTCF = "tissue37_CTCF.bed.gz",
  tissue54_H3K27ac = "tissue54_H3K27ac.bed.gz",
  tissue54_H3K4me3 = "tissue54_H3K4me3.bed.gz",
  tissue54_CTCF = "tissue54_CTCF.bed.gz",
  tissue54_ATAC = "tissue54_ATAC.bed.gz",
  PC3_H3K27ac = "PC3_H3K27ac.bed.gz",
  PC3_H3K4me3 = "PC3_H3K4me3.bed.gz",
  PC3_CTCF = "PC3_CTCF.bed.gz",
  PC3_ATAC = "PC3_ATAC.bed.gz"
)


# TxDb
txdb <- txdbmaker::makeTxDbFromGFF("GENCODE_v32.gtf")
seqlevelsStyle(txdb) <- "UCSC"
genome(txdb) <- "hg38"

# Function
annotate_and_summarize <- function(file, mark_name, txdb) {
  
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
    group_by(geneId) %>%
    summarise(
      promoter = sum(signalValue[grepl("Promoter", annotation)], na.rm = TRUE),
      distal = sum(signalValue[annotation == "Distal Intergenic"], na.rm = TRUE),
      genebody = sum(signalValue[grepl("Exon|Intron", annotation)], na.rm = TRUE),
      peak_count = n(),
      .groups = "drop"
    )
  
  # rename gene column
  colnames(features)[1] <- "gene"
  
  # refix columns
  colnames(features)[-1] <- paste0(mark_name, "_", colnames(features)[-1])
  
  return(features)
}

# Run loop
feature_list <- lapply(names(peak_files), function(mark) {
  annotate_and_summarize(peak_files[[mark]], mark, txdb)
})

# Merge
histone_features <- purrr::reduce(feature_list, full_join, by = "gene")

# Replace NA
histone_features[is.na(histone_features)] <- 0
summary(histone_features)
rowSums(histone_features[,-1]) == 0

# create peak-derived gene signals

# promoter signals

df <- df %>%
  mutate(
    signalValue = ifelse(
      "signalValue" %in% colnames(df),
      signalValue,
      1
    )
  )

promoter_features <- df %>%
  dplyr::filter(grepl("Promoter", annotation)) %>%
  dplyr::group_by(geneId) %>%
  dplyr::summarise(
    promoter_peak_count = dplyr::n(),
    promoter_signal_sum = sum(signalValue, na.rm = TRUE),
    promoter_signal_mean = mean(signalValue, na.rm = TRUE),
    .groups = "drop"
  )

# distal regulatory signal
df <- df %>%
  mutate(
    signalValue = ifelse(
      "signalValue" %in% colnames(df),
      signalValue,
      1
    )
  )

distal_features <- df %>%
  dplyr::filter(grepl("Distal Intergenic", annotation)) %>%
  dplyr::group_by(geneId) %>%
  dplyr::summarise(
    distal_peak_count = dplyr::n(),
    distal_signal_sum = sum(signalValue, na.rm = TRUE),
    distal_signal_mean = mean(signalValue, na.rm = TRUE),
    .groups = "drop"
  )

# gene body signal

df <- df %>%
  mutate(
    signalValue = ifelse(
      "signalValue" %in% colnames(df),
      signalValue,
      1
    )
  )

genebody_features <- df %>%
  dplyr::filter(grepl("Exon|Intron", annotation)) %>%
  dplyr::group_by(geneId) %>%
  dplyr::summarise(
    genebody_peak_count = dplyr::n(),
    genebody_signal_sum = sum(signalValue, na.rm = TRUE),
    genebody_signal_mean = mean(signalValue, na.rm = TRUE),
    .groups = "drop"
  )

# integrating ATACseq

df <- df %>%
  mutate(
    signalValue = ifelse(
      "signalValue" %in% colnames(df),
      signalValue,
      1
    )
  )

atac_distal_features <- df %>%
  dplyr::filter(grepl("Distal Intergenic", annotation)) %>%
  dplyr::group_by(geneId) %>%
  dplyr::summarise(
    atac_distal_peak_count = dplyr::n(),
    atac_distal_signal_sum = sum(signalValue, na.rm = TRUE),
    atac_distal_signal_mean = mean(signalValue, na.rm = TRUE),
    .groups = "drop"
  )

atac_promoter_features <- df %>%
  dplyr::filter(grepl("Promoter", annotation)) %>%
  dplyr::group_by(geneId) %>%
  dplyr::summarise(
    atac_promoter_peak_count = dplyr::n(),
    atac_promoter_signal_sum = sum(signalValue, na.rm = TRUE),
    atac_promoter_signal_mean = mean(signalValue, na.rm = TRUE),
    .groups = "drop"
  )

# integrating HiC
# hic_features <- data.frame(
#  SYMBOL = c("TP53", "MYC"),
#  hic_local_100kb = c(120, 300),
#  hic_local_1Mb = c(500, 1200),
#  compartment = c(0.8, -0.3),
#  tad_boundary_dist = c(20000, 5000)
# )

# merge all feature tables
features <- promoter_features %>%
  full_join(distal_features, by = "geneId") %>%
  full_join(genebody_features, by = "geneId") %>%
  full_join(atac_promoter_features, by = "geneId") %>%
  full_join(atac_distal_features, by = "geneId")

# replace NA
histone_features[is.na(histone_features)] <- 0

# remove duplicates if any
histone_features <- histone_features[!duplicated(histone_features$gene), ]

#expression data

files <- c("tissue37_RNAseq.tsv", "tissue54_RNAseq.tsv", "PC3_RNAseq1.tsv", "PC3_RNAseq2.tsv")

read_expr <- function(f) {
  message("Reading: ", f)
  
  expr_df <- read.table(f, header = TRUE, sep = "\t")
  
  # keep only ENSEMBL genes
  expr_df <- expr_df[grepl("^ENSG", expr_df$gene_id), ]
  
  # remove version numbers
  expr_df$gene <- sub("\\..*", "", expr_df$gene_id)
  
  # keep relevant columns
  expr_df <- expr_df[, c("gene", "TPM")]
  
  # rename TPM column to sample name
  colnames(expr_df)[2] <- tools::file_path_sans_ext(basename(f))
  
  return(expr_df)
}

expr_list <- lapply(files, read_expr)

str(expr_list[[1]])
head(expr_list[[1]])

# mapping
samples <- list(
  tissue37 = list(
    peaks = c("tissue37_H3K27ac.bed.gz", "tissue37_H3K4me3.bed.gz", "tissue37_CTCF.bed.gz", "tissue37_ATAC.bed.gz"),
    expr = "tissue37_RNAseq.tsv"
  ),
  tissue54 = list(
    peaks = c("tissue54_H3K27ac.bed.gz", "tissue54_H3K4me3.bed.gz", "tissue54_CTCF.bed.gz", "tissue54_ATAC.bed.gz"),
    expr = "tissue54_RNAseq.tsv"
  ),
  PC3_1 = list(
    peaks = c("PC3_H3K27ac.bed.gz", "PC3_H3K4me3.bed.gz", "PC3_CTCF.bed.gz", "PC3_ATAC.bed.gz"),
    expr = "PC3_RNAseq1.tsv"
  ),
  PC3_2 = list(
    peaks = c("PC3_H3K27ac.bed.gz", "PC3_H3K4me3.bed.gz", "PC3_CTCF.bed.gz", "PC3_ATAC.bed.gz"),
    expr = "PC3_RNAseq2.tsv"
  )
)

# chromatin marks per sample
process_sample_chromatin <- function(peak_files, txdb) {
  
  feature_list <- lapply(names(peak_files), function(mark) {
    annotate_and_summarize(peak_files[[mark]], mark, txdb)
  })
  
  chromatin <- purrr::reduce(feature_list, dplyr::full_join, by = "gene")
  chromatin[is.na(chromatin)] <- 0
  
  return(chromatin)
}

#expression per sample
expr_df <- read_expr(sample_info$expr)

#merge per sample
process_sample <- function(sample_name, sample_info, txdb) {
  
  message("Processing: ", sample_name)
  
  # build named peak list
  peak_files <- setNames(
    sample_info$peaks,
    tools::file_path_sans_ext(basename(sample_info$peaks))
  )
  
  chromatin <- process_sample_chromatin(peak_files, txdb)
  expr <- read_expr(sample_info$expr)
  
  merged <- dplyr::inner_join(chromatin, expr, by = "gene")
  
  merged$sample <- sample_name
  
  return(merged)
}

#full dataset
ml_list <- lapply(names(samples), function(s) {
  process_sample(s, samples[[s]], txdb)
})

ml_matrix <- dplyr::bind_rows(ml_list)

#adding labels
ml_matrix$label <- ifelse(
  grepl("PC3", ml_matrix$sample),
  "cancer",
  "normal"
)

#log transform to stabilize variance
expr_matrix[,-1] <- log2(expr_matrix[,-1] + 1)
numeric_cols <- sapply(ml_matrix, is.numeric)

ml_matrix[, numeric_cols] <- scale(ml_matrix[, numeric_cols])

#merge with chromatin data
ml_matrix <- dplyr::inner_join(
  histone_features,
  expr_matrix,
  by = "gene"
)