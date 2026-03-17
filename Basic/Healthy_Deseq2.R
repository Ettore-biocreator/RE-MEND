# =============================================================================
# RE-MEND PROJECT - Task 2: Differential Expression Analysis - Healthy Women
# Description: DEA between pre-partum (v38) and post-partum (pp) timepoints
#              in healthy (control) women, excluding SSRI users
# =============================================================================


# --- LIBRARIES ----------------------------------------------------------------
library(DESeq2)
library(dplyr)
library(ggplot2)
library(openxlsx)

set.seed(1234)


# --- CONFIGURATION ------------------------------------------------------------
paths <- list(
  counts      = "/home/tigem/e.aiello/Progetti/RE-MEND/RE_MEND_healthy/raw_count_data.txt",
  metadata    = "/home/tigem/e.aiello/Progetti/RE-MEND/RE_MEND_healthy/info.xlsx",
  annotation  = "/home/tigem/e.aiello/Progetti/RE-MEND/RE_MEND_healthy/mart_Human.txt",
  results_dir = "/home/tigem/e.aiello/Progetti/RE-MEND/RE_MEND_healthy/task1_DEA/"
)

# Outlier samples to exclude (identified during QC)
# Note: only one outlier here as DE25NGSUKBR129218 is a depression case,
#       already excluded by the Controls filter
outlier_samples <- c("DE79NGSUKBR129216")

# DESeq2 filtering thresholds
filter <- list(
  min_counts  = 25,
  min_samples = 10
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
  
  # Keep only common samples
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
  
  return(list(counts = counts, metadata = metadata, annotation = annotation))
}


# --- PREPROCESS METADATA ------------------------------------------------------
preprocess_metadata <- function(metadata, outliers) {
  # Remove rows with NA LABCode (can arise from Excel header misreading)
  metadata <- metadata[!is.na(metadata$LABCode), ]
  
  # Keep only healthy controls
  metadata <- metadata[metadata$RB_PPD_trajectory_EPDS == "Controls", ]
  
  # Remove outlier samples
  metadata <- metadata %>% filter(!LABCode %in% outliers)
  
  # Exclude SSRI users (pregnancy and post-partum)
  metadata <- metadata[complete.cases(metadata$NK_SSRI_pp),]
  metadata <- metadata[complete.cases(metadata$NK_SSRI_pregnancy_loose),]
  metadata <- metadata[metadata$NK_SSRI_pregnancy_loose == "No", ]
  metadata <- metadata[metadata$NK_SSRI_pp == "No", ]
  
  # BMI: convert to numeric, round, impute NA with mean, center and scale
  metadata$v17_BMI_innan_R <- as.numeric(metadata$v17_BMI_innan_R)
  metadata$v17_BMI_innan_R <- round(metadata$v17_BMI_innan_R)
  bmi_mean <- round(mean(metadata$v17_BMI_innan_R, na.rm = TRUE))
  metadata$v17_BMI_innan_R[is.na(metadata$v17_BMI_innan_R)] <- bmi_mean
  metadata$v17_BMI_innan_R <- scale(metadata$v17_BMI_innan_R)[, 1]
  
  # Age: convert to numeric, center and scale
  metadata$NK_Age_at_partus <- as.numeric(metadata$NK_Age_at_partus)
  metadata$NK_Age_at_partus <- scale(metadata$NK_Age_at_partus)[, 1]
  
  # Timepoint as factor (reference: v38 pre-partum)
  metadata$Timepoint <- factor(metadata$Timepoint, levels = c("v38", "pp"))
  
  message(sprintf("Samples after filtering: %d", nrow(metadata)))
  
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
    design    = ~ Timepoint + NK_Age_at_partus + v17_BMI_innan_R
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
  res        <- results(dds, contrast = c("Timepoint", "pp", "v38"))
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
  
  message(sprintf("pp vs v38: %d DEGs (%d up, %d down)",
                  nrow(res_deg), length(up_genes), length(down_genes)))
  
  return(list(
    full    = res_df,
    sig     = res_deg,
    up      = up_genes,
    down    = down_genes,
    raw_res = res
  ))
}


# --- PLOTS --------------------------------------------------------------------
plot_pca <- function(dds) {
  vsd      <- vst(dds, blind = FALSE)
  pca_data <- plotPCA(vsd, intgroup = "Timepoint", returnData = TRUE)
  pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)
  
  ggplot(pca_data, aes(x = PC1, y = PC2, color = Timepoint)) +
    geom_point(size = 3) +
    scale_color_manual(
      values = c("v38" = "#f8766d", "pp" = "#00bfc4"),
      labels = c("v38" = "Pre-partum", "pp" = "Post-partum"),
      name   = "Timepoint"
    ) +
    xlab(paste0("PC1: ", pct_var[1], "% variance")) +
    ylab(paste0("PC2: ", pct_var[2], "% variance")) +
    ggtitle("PCA - Healthy women: Pre vs Post-partum") +
    coord_fixed() +
    theme_bw()
}


plot_volcano <- function(raw_res, thresholds) {
  df <- as.data.frame(raw_res)
  df$log10padj  <- -log10(df$padj)
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
    ggtitle("Volcano Plot - Post-partum vs Pre-partum (Healthy)") +
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
              file.path(results_dir, "pp_vs_v38_healthy_full.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Filtered DEGs only
  write.table(results_list$sig,
              file.path(results_dir, "pp_vs_v38_healthy_DEGs.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  message(sprintf("Results saved in: %s", results_dir))
}


# =============================================================================
# MAIN
# =============================================================================

# Load data
data_list <- load_data(paths)

# Preprocess metadata
meta <- preprocess_metadata(data_list$metadata, outlier_samples)

# Run DESeq2
dds <- run_deseq2(data_list$counts, meta, filter)

# Extract and annotate results
res <- extract_results(dds, data_list$annotation, deg_thresholds)

# Save results
save_results(res, paths$results_dir)

# Plots
print(plot_pca(dds))
print(plot_volcano(res$raw_res, deg_thresholds))
