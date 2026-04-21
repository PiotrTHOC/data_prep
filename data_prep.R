if (!requireNamespace("TCGAbiolinks", quietly = TRUE))
  BiocManager::install("TCGAbiolinks")


## load libraries
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
library(TCGAbiolinks)
library(SummarizedExperiment)
library(tidyr)


## TxDb and gene mapping
txdb <- txdbmaker::makeTxDbFromGFF(here::here("GENCODE_v32.gtf"))
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


## Core functions
# 1. annotate peaks and summarize per gene
annotate_and_summarize <- function(file, mark_name, txdb, tx2gene_clean) {

  if (!file.exists(file)) {
    stop(paste("Missing file:", file))
  }

  peak <- readPeakFile(file)

  # normalize signal column
  peak$signalValue <- as.numeric(peak$V7)
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

# 2.Extract HiC distal peaks
extract_distal_peaks <- function(file) {
  
  peak <- readPeakFile(file)
  
  peak$signalValue <- as.numeric(peak$V7)
  genome(peak) <- "hg38"
  seqlevelsStyle(peak) <- "UCSC"
  
  return(peak)
}

# 3.Build HiC interactions
load_hic_interactions <- function(hic_file) {
  
  message("Processing Hi-C: ", hic_file)
  
  hic <- read.table(hic_file, header = FALSE, stringsAsFactors = FALSE)
  
  # ensure at least 6 columns
  if (ncol(hic) < 6) {
    stop("Hi-C file has fewer than 6 columns (invalid BEDPE)")
  }
  
  # assign basic columns
  colnames(hic)[1:6] <- c(
    "chr1","start1","end1",
    "chr2","start2","end2"
  )
  
  # build GRanges
  gr1 <- GRanges(hic$chr1, IRanges(hic$start1, hic$end1))
  gr2 <- GRanges(hic$chr2, IRanges(hic$start2, hic$end2))
  
  genome(gr1) <- "hg38"
  genome(gr2) <- "hg38"
  
  seqlevelsStyle(gr1) <- "UCSC"
  seqlevelsStyle(gr2) <- "UCSC"
  
  return(list(gr1 = gr1, gr2 = gr2))
}


# 4.Compute HiC-connected enhancer signals
compute_hic_connected_distal <- function(distal_peaks, hic_obj) {
  
  gr1 <- hic_obj$gr1
  gr2 <- hic_obj$gr2
  prom_links <- hic_obj$links
  
  # safety check
  if (nrow(prom_links) == 0 || length(distal_peaks) == 0) {
    return(data.frame(gene=character(), distal_connected=numeric()))
  }
  
  distal_hits1 <- findOverlaps(distal_peaks, gr1)
  distal_hits2 <- findOverlaps(distal_peaks, gr2)
  
  distal_map <- rbind(
    data.frame(
      peak = queryHits(distal_hits1),
      region = subjectHits(distal_hits1),
      side = "gr1"
    ),
    data.frame(
      peak = queryHits(distal_hits2),
      region = subjectHits(distal_hits2),
      side = "gr2"
    )
  )
  
  # SAFE EXIT
  if (!is.data.frame(distal_map) || nrow(distal_map) == 0) {
    return(data.frame(gene=character(), distal_connected=numeric()))
  }
  
  merged <- merge(prom_links, distal_map, by=c("region","side"))
  
  if (nrow(merged) == 0) {
    return(data.frame(gene=character(), distal_connected=numeric()))
  }
  
  merged$signal <- distal_peaks$signalValue[merged$peak]
  merged$signal[is.na(merged$signal)] <- 0
  
  merged <- merged %>%
    dplyr::distinct(gene, peak, .keep_all = TRUE)
  
  result <- merged %>%
    group_by(gene) %>%
    summarise(
      distal_connected = sum(signal),
      .groups = "drop"
    )
  
  return(result)
}


# 5.final merge
safe_full_join <- function(list_of_tables) {
  Reduce(function(x, y) {
    dplyr::full_join(x, y, by="gene")
  }, list_of_tables)
}


## Sample definitions
samples <- list(
  tissue37 = list(
    peaks = c(
      "tissue37_H3K27ac" = here::here("tissue37_H3K27ac.bed.gz"),
      "tissue37_H3K4me3" = here::here("tissue37_H3K4me3.bed.gz"),
      "tissue37_CTCF"    = here::here("tissue37_CTCF.bed.gz")
    ),
    hic = here::here("tissue54_HiC.bedpe.gz"),
    expr = here::here("tissue37_RNAseq.tsv")
  ),
  tissue54 = list(
    peaks = c(
      "tissue54_H3K27ac" = here::here("tissue54_H3K27ac.bed.gz"),
      "tissue54_H3K4me3" = here::here("tissue54_H3K4me3.bed.gz"),
      "tissue54_CTCF"    = here::here("tissue54_CTCF.bed.gz"),
      "tissue54_ATAC"    = here::here("tissue54_ATAC.bed.gz")
    ),
    hic = here::here("tissue54_HiC.bedpe.gz"),
    expr = here::here("tissue54_RNAseq.tsv")
  ),
  PC3_1 = list(
    peaks = c(
      "PC3_H3K27ac" = here::here("PC3_H3K27ac.bed.gz"),
      "PC3_H3K4me3" = here::here("PC3_H3K4me3.bed.gz"),
      "PC3_CTCF"    = here::here("PC3_CTCF.bed.gz"),
      "PC3_ATAC"    = here::here("PC3_ATAC.bed.gz")
    ),
    hic = here::here("PC3_HiC.bedpe.gz"),
    expr = here::here("PC3_RNAseq1.tsv")
  ),
  PC3_2 = list(
    peaks = c(
      "PC3_H3K27ac" = here::here("PC3_H3K27ac.bed.gz"),
      "PC3_H3K4me3" = here::here("PC3_H3K4me3.bed.gz"),
      "PC3_CTCF"    = here::here("PC3_CTCF.bed.gz"),
      "PC3_ATAC"    = here::here("PC3_ATAC.bed.gz")
    ),
    hic = here::here("PC3_HiC.bedpe.gz"),
    expr = here::here("PC3_RNAseq2.tsv")
  )
)

## Processing functions  
# Chromatin marks per sample
process_sample_chromatin <- function(peak_files, hic_file, txdb, tx2gene_clean) {
  
  hic_obj <- process_hic_sample(hic_file, txdb, tx2gene_clean)
  
  feature_list <- list()
  
  for (mark in names(peak_files)) {
    
    file <- peak_files[[mark]]
    
    base <- annotate_and_summarize(file, mark, txdb, tx2gene_clean)
    
    connected <- data.frame(
      gene = character(),
      distal_connected = numeric()
    )
    
    if (grepl("H3K27ac|ATAC", mark)) {

      distal_peaks <- readPeakFile(file)
      distal_peaks$signalValue <- as.numeric(distal_peaks$V7)

      connected <- compute_hic_connected_distal(distal_peaks, hic_obj)

      colnames(connected)[colnames(connected) != "gene"] <-
        paste0(mark, "_", colnames(connected)[colnames(connected) != "gene"])

      base <- dplyr::left_join(base, connected, by="gene")
    }
    
    feature_list[[mark]] <- base
  }
  
  chromatin <- safe_full_join(feature_list)
  chromatin[is.na(chromatin)] <- 0
  
  chromatin
}

# HiC marks
process_hic_sample <- function(hic_file, txdb, tx2gene_clean) {
  
  message("Processing Hi-C: ", hic_file)
  
  hic <- read.table(hic_file, header = FALSE)
  
  colnames(hic)[1:6] <- c("chr1","start1","end1","chr2","start2","end2")
  
  gr1 <- GRanges(hic$chr1, IRanges(hic$start1, hic$end1))
  gr2 <- GRanges(hic$chr2, IRanges(hic$start2, hic$end2))
  
  genome(gr1) <- "hg38"
  genome(gr2) <- "hg38"
  
  seqlevelsStyle(gr1) <- "UCSC"
  seqlevelsStyle(gr2) <- "UCSC"
  
  # promoters
  prom <- promoters(txdb, upstream = 1000, downstream = 100)
  
  prom_df <- as.data.frame(prom)
  prom_df$tx_id <- sub("\\..*", "", prom_df$tx_name)
  prom_df$gene <- tx2gene_clean[prom_df$tx_id]
  
  keep <- !is.na(prom_df$gene)
  prom <- prom[keep]
  prom$gene <- prom_df$gene[keep]
  
  hits1 <- findOverlaps(prom, gr1)
  hits2 <- findOverlaps(prom, gr2)
  
  # SAFE early exit
  if (length(hits1) == 0 && length(hits2) == 0) {
    return(list(
      links = data.frame(gene=character(), region=integer(), side=character()),
      gr1 = gr1,
      gr2 = gr2
    ))
  }
  
  prom_links <- rbind(
    data.frame(
      gene = prom$gene[queryHits(hits1)],
      region = subjectHits(hits1),
      side = "gr1"
    ),
    data.frame(
      gene = prom$gene[queryHits(hits2)],
      region = subjectHits(hits2),
      side = "gr2"
    )
  )
  
  prom_links <- prom_links[complete.cases(prom_links), ]
  
  list(
    links = prom_links,
    gr1 = gr1,
    gr2 = gr2
  )
}


## Execution
# Build per-sample chromatin features (no local RNA-seq)
process_sample <- function(sample_name, sample_info, txdb, tx2gene_clean) {

  message("Processing: ", sample_name)

  chromatin <- process_sample_chromatin(sample_info$peaks, sample_info$hic, txdb, tx2gene_clean)

  chromatin$sample <- sample_name

  return(chromatin)
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


# Scale/center all numeric features (applied after TCGA join below)
scale_per_feature <- function(df) {
  num_cols <- sapply(df, is.numeric)
  df[, num_cols] <- lapply(df[, num_cols], function(x) {
    if (all(is.na(x))) return(x)
    scale(x)
  })
  return(df)
}


# download TCGA prostate RNA-seq
query <- GDCquery(
  project = "TCGA-PRAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)
tcga_data <- GDCprepare(query)

# expression matrix (rows = genes, cols = samples), CPM-normalized
tcga_expr <- assay(tcga_data)
tcga_cpm <- apply(tcga_expr, 2, function(x) (x / sum(x)) * 1e6)

# strip ENSG version suffixes so gene keys match the chromatin matrix
rownames(tcga_cpm) <- sub("\\..*", "", rownames(tcga_cpm))


# label samples: TP = tumor (cancer), NT = solid tissue normal (adjacent tissue)
sample_types <- as.character(colData(tcga_data)$shortLetterCode)
tcga_labels <- ifelse(sample_types == "TP", "cancer",
                      ifelse(sample_types == "NT", "normal", NA))

valid <- !is.na(tcga_labels)
tcga_cpm    <- tcga_cpm[, valid]
tcga_labels <- tcga_labels[valid]

message("TCGA samples per label: ",
        paste(names(table(tcga_labels)), table(tcga_labels),
              sep = "=", collapse = ", "))


# log2(TPM+1) so variance reflects fold-change-scale dispersion
tcga_log <- log2(tcga_cpm + 1)

cancer_cols <- which(tcga_labels == "cancer")
normal_cols <- which(tcga_labels == "normal")

row_var <- function(mat, cols) {
  apply(mat[, cols, drop = FALSE], 1, stats::var, na.rm = TRUE)
}

tcga_agg <- data.frame(
  gene             = rownames(tcga_log),
  tcga_mean_cancer = rowMeans(tcga_log[, cancer_cols, drop = FALSE], na.rm = TRUE),
  tcga_var_cancer  = row_var (tcga_log, cancer_cols),
  tcga_mean_normal = rowMeans(tcga_log[, normal_cols, drop = FALSE], na.rm = TRUE),
  tcga_var_normal  = row_var (tcga_log, normal_cols),
  row.names = NULL,
  stringsAsFactors = FALSE
)

# stripping ENSG version suffixes can produce duplicate gene keys; collapse them
tcga_agg <- tcga_agg %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(dplyr::across(dplyr::everything(), mean), .groups = "drop")


# join TCGA aggregates onto chromatin matrix, then scale all numeric features
ml_matrix <- dplyr::inner_join(ml_matrix, tcga_agg, by = "gene")

ml_matrix <- scale_per_feature(ml_matrix)

summary(ml_matrix)

write.csv(ml_matrix, here::here("ml_matrix_tcga.csv"), row.names = FALSE)
