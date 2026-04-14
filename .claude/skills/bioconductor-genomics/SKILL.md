---
name: bioconductor-genomics
description: Best practices for Bioconductor genomic analysis — GRanges, ChIPseeker, GenomicFeatures, peak annotation, and multi-omics integration. Use when working with genomic intervals, peak files, GTF/GFF annotations, or Bioconductor objects.
user-invocable: false
paths: "*.R,*.Rmd,*.bed,*.bed.gz,*.gtf,*.gff,*.bigWig,*.bedpe.gz"
---

# Bioconductor & Genomic Data Analysis

## Core Data Structures

- Use `GRanges` as the canonical representation for genomic intervals — never plain data.frames with chr/start/end columns
- Always set genome and seqlevels style consistently:
  ```r
  seqlevelsStyle(gr) <- "UCSC"   # chr1, chr2, ...
  genome(gr) <- "hg38"
  ```
- Use `GRangesList` for grouped intervals (e.g., peaks per sample)
- Convert to data.frame only at the final output stage, not for intermediate processing

## Genome Annotations

- Build TxDb from GTF once and cache it; never rebuild inside loops:
  ```r
  txdb <- txdbmaker::makeTxDbFromGFF("GENCODE_v32.gtf")
  seqlevelsStyle(txdb) <- "UCSC"
  genome(txdb) <- "hg38"
  ```
- Use `org.Hs.eg.db` for gene ID mapping (Entrez <-> Symbol <-> ENSEMBL)
- Always specify `columns` in `AnnotationDbi::select()` to avoid pulling unnecessary data
- When working with GENCODE, strip version suffixes from ENSEMBL IDs: `sub("\\..*", "", gene_id)`

## Peak File Handling

- Read peak files with `rtracklayer::import()` (generic) or `ChIPseeker::readPeakFile()` (narrowPeak/broadPeak)
- Validate peak files on load: check column count, chromosome naming, coordinate ranges
- Always filter to standard chromosomes:
  ```r
  peak <- keepStandardChromosomes(peak, pruning.mode = "coarse")
  ```
- Document the peak caller and parameters used (MACS2 q-value, IDR threshold)

## Peak Annotation with ChIPseeker

- Set biologically meaningful TSS regions — default `c(-3000, 3000)` is wide; for promoter-focused analysis use `c(-1000, 100)` or `c(-2000, 200)`
- Always specify `TxDb` and optionally `annoDb` for gene symbol mapping:
  ```r
  peakAnno <- annotatePeak(peak, TxDb = txdb, annoDb = "org.Hs.eg.db",
                           tssRegion = c(-1000, 100))
  ```
- Understand annotation priority: ChIPseeker assigns each peak to ONE annotation category. Overlapping categories follow a priority order (Promoter > 5'UTR > 3'UTR > Exon > Intron > Downstream > Intergenic)
- Use `as.data.frame(peakAnno)` to extract results for downstream processing

## Signal Quantification

- When summarizing signal per gene, aggregate by annotation category:
  - **Promoter signal**: peaks overlapping TSS region
  - **Gene body signal**: peaks in exons + introns
  - **Distal signal**: intergenic peaks (enhancers, insulators)
- Use `signalValue` (column 7 in narrowPeak) as the primary signal metric; `pValue` and `qValue` for filtering
- Always use `na.rm = TRUE` in aggregation functions
- Count peaks per gene as a complementary feature to signal intensity

## Multi-Omics Integration

- Use gene-level identifiers (ENSEMBL or Entrez) as the join key across data types
- When merging ChIP-seq, ATAC-seq, and RNA-seq:
  1. Annotate peaks to genes separately per mark/assay
  2. Summarize per gene (signal sum, mean, peak count)
  3. Merge with `full_join()` or `inner_join()` depending on analysis goal
  4. Replace NA with 0 for marks with no peaks at a gene
- Prefix all feature columns with the assay/mark name: `H3K27ac_promoter`, `ATAC_distal`
- Normalize across samples before merging (quantile normalization, TPM for RNA-seq)

## ChIP-seq Specific

- H3K27ac: active enhancers and promoters — expect broad/sharp peaks
- H3K4me3: active promoters — expect sharp peaks near TSS
- CTCF: insulator elements, chromatin boundaries — expect sharp peaks
- Distinguish activating marks (H3K27ac, H3K4me3) from repressive marks (H3K27me3, H3K9me3)

## ATAC-seq Specific

- ATAC-seq measures chromatin accessibility, not a specific modification
- Peaks at promoters indicate open chromatin at active/poised genes
- Distal ATAC peaks often overlap enhancers — cross-reference with H3K27ac
- Consider nucleosome-free regions (NFR) vs. mono-nucleosome fragments

## Quality Control

- Check peak count per sample — flag samples with < 5,000 or > 200,000 peaks
- Verify FRiP (Fraction of Reads in Peaks) > 1% for ChIP, > 15% for ATAC
- Compare peak width distributions across samples — outliers indicate technical issues
- Use `plotAnnoBar()` and `plotDistToTSS()` from ChIPseeker for annotation QC
