# =============================================================================
# RE-MEND PROJECT - Task 2b: Paired DEA - Healthy Women
# Description: DEA between pre-partum (v38) and post-partum (pp) timepoints
#              using paired control samples only
# =============================================================================


# --- LIBRARIES ----------------------------------------------------------------
library(DESeq2)
library(dplyr)
library(ggplot2)
library(openxlsx)

set.seed(1234)


# --- CONFIGURATION ------------------------------------------------------------
paths <- list(
  counts      = "path/to/raw_count_data.txt",
  metadata    = "path/to/info.xlsx",
  paired      = "path/to/BASIC_Paired.xlsx",
  annotation  = "path/to/mart_Human.txt",
  results_dir = "results/task2b_paired_DEA/"
)

# DESeq2 filtering thresholds
filter <- list(
  min_counts  = 25,
  min_samples = 4    # lower threshold due to smaller paired sample size
)

# DEG thresholds
deg_thresholds <- list(
  padj   = 0.05,
  log2fc = 1
)


# --- LOAD DATA ----------------------------------------------------------------
load_data <- function(paths) {
  # Count matrix
  counts <- read.delim2(paths$counts)
  rownames(counts) <- counts$Ensembl
  counts <- counts[, -1]

  # Metadata
  metadata <- read.xlsx(paths$metadata)
  metadata <- metadata[order(metadata$LABCode), ]
  metadata <- metadata[!is.na(metadata$LABCode), ]

  # Paired sample list
  paired <- read.xlsx(paths$paired)

  # Keep only common samples between counts and metadata
  common_samples <- intersect(colnames(counts), metadata$LABCode)
  counts   <- counts[, common_samples]
  metadata <- metadata[metadata$LABCode %in% common_samples, ]
  metadata <- metadata[match(colnames(counts), metadata$LABCode), ]

  stopifnot("Sample order mismatch between counts and metadata" =
              all(metadata$LABCode == colnames(counts)))

  # Annotation
  annotation <- read.table(paths$annotation, header = TRUE, sep = "\t")
  annotation <- annotation[!duplicated(annotation$Gene.stable.ID), ]

  message(sprintf("Loaded %d genes x %d samples", nrow(counts), ncol(counts)))

  return(list(counts = counts, metadata = metadata,
              paired = paired, annotation = annotation))
}


# --- PREPROCESS METADATA ------------------------------------------------------
preprocess_metadata <- function(metadata, paired) {
  # Keep only paired control samples
  paired_controls <- paired[paired$Trajectory_EPDS == "Control", "Sample_Name"]
  metadata <- metadata[metadata$LABCode %in% paired_controls, ]

  # Timepoint as factor (reference: v38 pre-partum)
  metadata$Timepoint <- factor(metadata$Timepoint, levels = c("v38", "pp"))

  message(sprintf("Paired control samples: %d", nrow(metadata)))

  return(metadata)
}


# --- RUN DESEQ2 ---------------------------------------------------------------
run_deseq2 <- function(counts, metadata, filter) {
  counts_sub <- counts[, metadata$LABCode]

  stopifnot("Sample order mismatch before DESeq2" =
              all(metadata$LABCode == colnames(counts_sub)))

  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData   = metadata,
    design    = ~ Timepoint
  )

  # Low-count filtering
  keep <- rowSums(counts(dds) >= filter$min_counts) >= filter$min_samples
  dds  <- dds[keep, ]

  message(sprintf("Genes retained after filtering: %d", sum(keep)))

  dds <- DESeq(dds)
  return(dds)
}


# --- EXTRACT AND ANNOTATE RESULTS ---------------------------------------------
extract_results <- function(dds, annotation, thresholds) {
  res         <- results(dds, contrast = c("Timepoint", "pp", "v38"))
  res_ordered <- res[order(res$padj), ]

  # Annotate
  res_df <- as.data.frame(res_ordered)
  res_df$Gene.stable.ID <- sub("\\..*", "", rownames(res_df))
  res_df <- res_df %>% left_join(annotation, by = "Gene.stable.ID")
  res_df <- res_df %>%
    select(Gene.name, Gene.stable.ID, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

  # Filter DEGs
  res_complete <- res_df[complete.cases(res_df$padj), ]
  res_sig      <- res_complete[res_complete$padj < thresholds$padj, ]
  res_deg      <- res_sig[abs(res_sig$log2FoldChange) > thresholds$log2fc, ]

  up_genes   <- res_deg[res_deg$log2FoldChange > 0, "Gene.name"]
  down_genes <- res_deg[res_deg$log2FoldChange < 0, "Gene.name"]

  message(sprintf("pp vs v38 (paired): %d DEGs (%d up, %d down)",
                  nrow(res_deg), length(up_genes), length(down_genes)))

  return(list(
    full    = res_df,
    sig     = res_deg,
    raw_res = res
  ))
}


# --- VOLCANO PLOT -------------------------------------------------------------
plot_volcano <- function(raw_res, thresholds) {
  df <- as.data.frame(raw_res)
  df$log10padj   <- -log10(df$padj)
  df$significant <- abs(df$log2FoldChange) > thresholds$log2fc &
                    !is.na(df$padj) & df$padj < thresholds$padj

  ggplot(df, aes(x = log2FoldChange, y = log10padj, color = significant)) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_color_manual(
      values = c("FALSE" = "#181C14", "TRUE" = "#f8766d"),
      labels = c("Not significant", "Significant"),
      name   = NULL
    ) +
    geom_vline(xintercept = c(-thresholds$log2fc, thresholds$log2fc),
               linetype = "dashed", color = "blue") +
    geom_hline(yintercept = -log10(thresholds$padj),
               linetype = "dashed", color = "blue") +
    xlab("log2(Fold Change)") +
    ylab("-log10(adjusted p-value)") +
    ggtitle("Volcano Plot - Post-partum vs Pre-partum (Paired Controls)") +
    theme_bw()
}


# --- SAVE RESULTS -------------------------------------------------------------
save_results <- function(results_list, results_dir) {
  if (is.null(results_dir) || results_dir == "") {
    stop("results_dir is not set. Please update paths$results_dir in the configuration block.")
  }
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  # Full DEA table
  write.table(results_list$full,
              file.path(results_dir, "pp_vs_v38_paired_full.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  # Filtered DEGs
  write.table(results_list$sig,
              file.path(results_dir, "pp_vs_v38_paired_DEGs.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  message(sprintf("Results saved in: %s", results_dir))
}


# =============================================================================
# MAIN
# =============================================================================

# Load data
data_list <- load_data(paths)

# Preprocess metadata
meta <- preprocess_metadata(data_list$metadata, data_list$paired)

# Run DESeq2
dds <- run_deseq2(data_list$counts, meta, filter)

# Extract and annotate results
res <- extract_results(dds, data_list$annotation, deg_thresholds)

# Save results
save_results(res, paths$results_dir)

# Volcano plot
print(plot_volcano(res$raw_res, deg_thresholds))
