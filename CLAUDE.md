# Chromatin Data — ML Feature Engineering Pipeline

## Project Overview

This project processes multi-omics chromatin data (ChIP-seq, ATAC-seq, RNA-seq) from prostate tissue and cell line samples, annotates peaks to genes, and constructs a gene-level feature matrix for machine learning classification (cancer vs. normal).

## Data Sources

- **Genome**: Human hg38 (GRCh38), GENCODE v32 gene annotation (`GENCODE_v32.gtf`)
- **ChIP-seq marks**: H3K27ac (active enhancers/promoters), H3K4me3 (active promoters), CTCF (insulators)
- **ATAC-seq**: Chromatin accessibility
- **RNA-seq**: Gene expression (TPM values)
- **Samples**:
  - `tissue37` — normal prostate tissue
  - `tissue54` — normal prostate tissue
  - `PC3` — prostate cancer cell line (2 RNA-seq replicates)
  - Additional samples: LNCaP, VCaP, 22Rv1, RWPE1, RWPE2, epithelial (available but not all used in current pipeline)

## Pipeline Architecture

```
BED peak files ──> ChIPseeker annotation ──> per-gene signal summary ──> merged feature matrix
                        (via TxDb)               (promoter/distal/       (gene x features)
                                                   genebody signals)
                                                                              │
RNA-seq TSV ──────────────────────────────────────────────────────────────> inner_join
                                                                              │
                                                                         ML matrix + labels
```

### Processing Steps

1. **TxDb construction** — Build transcript database from GENCODE v32 GTF, set UCSC chromosome style
2. **Peak annotation** — For each ChIP-seq/ATAC-seq sample, annotate peaks to nearest gene using ChIPseeker with TSS region `c(-1000, 100)`
3. **Signal summarization** — Per gene, compute signal sums for promoter, distal intergenic, and gene body regions; count peaks
4. **Cross-mark merging** — `full_join` all marks for a sample into one wide table (NA → 0)
5. **Expression integration** — Read RNA-seq TPM values, strip ENSEMBL version suffixes, `inner_join` with chromatin features
6. **Sample stacking** — `bind_rows` across samples, add sample and label columns
7. **Normalization** — Log2 transform expression, scale/center all numeric features

## Key Functions

- `annotate_and_summarize(file, mark_name, txdb)` — Core function: reads a peak file, annotates with ChIPseeker, summarizes signal by gene and annotation category
- `read_expr(f)` — Reads RNA-seq TSV, filters to ENSEMBL genes, strips version suffixes
- `process_sample_chromatin(peak_files, txdb)` — Processes all chromatin marks for one sample
- `process_sample(sample_name, sample_info, txdb)` — Full per-sample pipeline: chromatin + expression + merge

## Feature Columns

Each mark produces 4 features per gene: `{mark}_promoter`, `{mark}_distal`, `{mark}_genebody`, `{mark}_peak_count`. Additional granular features: `promoter_signal_sum/mean`, `distal_signal_sum/mean`, `genebody_signal_sum/mean`, `atac_promoter_*`, `atac_distal_*`.

## Gene ID Convention

- Chromatin features use Entrez gene IDs (from ChIPseeker annotation via `geneId`)
- RNA-seq uses ENSEMBL IDs (stripped of version suffix)
- The join key between chromatin and expression data must be reconciled — currently uses `gene` column

## Classification Target

- `label`: binary — `"cancer"` (PC3 samples) vs `"normal"` (tissue samples)
- Determined by sample name: `grepl("PC3", sample)` → cancer

## R Dependencies

- **Bioconductor**: ChIPseeker, GenomicRanges, GenomicFeatures, rtracklayer, txdbmaker, GenomeInfoDb, TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db
- **Tidyverse**: dplyr, purrr, tidyr
- **Base R**: tools

## Conventions

- Use `dplyr::` explicit namespacing for all dplyr verbs
- Use `here::here()` for file paths — avoid `setwd()` and absolute paths
- All signal aggregations use `na.rm = TRUE`
- NA from `full_join` means "no peak detected" — replace with 0
- Column naming: `{sample}_{mark}_{feature}` pattern
- TSS region: `c(-1000, 100)` upstream/downstream
- Chromosome style: UCSC (`chr1`, `chr2`, ...)
- Genome assembly: hg38
