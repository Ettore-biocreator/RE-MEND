# =============================================================================
# RE-MEND PROJECT - Task 1: Differential Expression Analysis
# Description: DEA on depressed vs non-depressed women at two timepoints
#              (pre-partum: v38, post-partum: pp) with drug repurposing output
# =============================================================================


# --- LIBRARIES ----------------------------------------------------------------
library(openxlsx)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)


# --- CONFIGURATION ------------------------------------------------------------
# Modify these paths before running

paths <- list(
  counts      = "path/to/raw_count_data.txt",
  metadata    = "path/to/info.xlsx",
  annotation  = "path/to/annotation.txt",
  results_dir = "results/task1_DEA/"
)

# Outlier samples to exclude (identified during QC)
outlier_samples <- c("DE79NGSUKBR129216", "DE25NGSUKBR129218")

# DESeq2 filtering thresholds
filter <- list(
  min_counts      = 25,   # minimum counts per gene
  min_samples     = 10    # minimum number of samples with min_counts
)

# DEG thresholds
deg_thresholds <- list(
  padj     = 0.1,
  log2fc   = 1
)


# --- LOAD DATA ----------------------------------------------------------------
load_data <- function(paths) {
  # Count matrix
  counts <- read.delim2(paths$counts)
  rownames(counts) <- counts$Ensembl
  counts <- counts[, -1]  # remove Ensembl column

  # Metadata
  metadata <- read.xlsx(paths$metadata)
  metadata <- metadata[order(metadata$LABCode), ]

  # Keep only common samples between count matrix and metadata
  common_samples <- intersect(colnames(counts), metadata$LABCode)
  counts <- counts[, common_samples]
  metadata <- metadata[metadata$LABCode %in% common_samples, ]
  metadata <- metadata[match(colnames(counts), metadata$LABCode), ]

  stopifnot("Sample order mismatch between counts and metadata" =
              all(metadata$LABCode == colnames(counts)))

  # Annotation: load and remove duplicated Ensembl IDs
  annotation <- read.table(paths$annotation, header = TRUE, sep = "\t")
  annotation <- annotation[!duplicated(annotation$Gene.stable.ID), ]

  message(sprintf("Loaded %d genes x %d samples", nrow(counts), ncol(counts)))

  return(list(counts = counts, metadata = metadata, annotation = annotation))
}


# --- PREPROCESS METADATA ------------------------------------------------------
preprocess_metadata <- function(metadata, timepoint, outliers) {
  # Filter by timepoint: "pp" = post-partum, "v38" = pre-partum
  metadata <- metadata[metadata$Timepoint == timepoint, ]

  # Remove outlier samples
  metadata <- metadata %>% filter(!LABCode %in% outliers)

  # BMI: convert to numeric, round, impute NA with mean
  metadata$v17_BMI_innan_R <- as.numeric(metadata$v17_BMI_innan_R)
  metadata$v17_BMI_innan_R <- round(metadata$v17_BMI_innan_R)
  bmi_mean <- round(mean(metadata$v17_BMI_innan_R, na.rm = TRUE))
  metadata$v17_BMI_innan_R[is.na(metadata$v17_BMI_innan_R)] <- bmi_mean

  # Standardize trajectory labels (remove spaces)
  metadata$RB_PPD_trajectory <- gsub(" ", "_", metadata$RB_PPD_trajectory)
  metadata$RB_PPD_trajectory_EPDS <- gsub(" ", "_", metadata$RB_PPD_trajectory_EPDS)

  # Convert design variables to factors (required by DESeq2)
  metadata$RB_PPD_trajectory_EPDS <- factor(metadata$RB_PPD_trajectory_EPDS)
  metadata$NK_Age_at_partus        <- factor(metadata$NK_Age_at_partus)

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
    design    = ~ RB_PPD_trajectory_EPDS + NK_Age_at_partus + v17_BMI_innan_R
  )

  # Low-count filtering
  keep <- rowSums(counts(dds) >= filter$min_counts) >= filter$min_samples
  dds  <- dds[keep, ]

  message(sprintf("Genes retained after filtering: %d", sum(keep)))

  dds <- DESeq(dds)
  return(dds)
}


# --- EXTRACT AND ANNOTATE RESULTS ---------------------------------------------
extract_results <- function(dds, contrast_group, reference_group, annotation, thresholds) {
  res <- results(dds, contrast = c("RB_PPD_trajectory_EPDS", contrast_group, reference_group))
  res_ordered <- res[order(res$padj), ]

  # Annotate: strip Ensembl version suffix (e.g. ENSG00000001.12 -> ENSG00000001)
  res_df <- as.data.frame(res_ordered)
  res_df$Gene.stable.ID <- sub("\\..*", "", rownames(res_df))
  res_df <- res_df %>% left_join(annotation, by = "Gene.stable.ID")
  res_df <- res_df[, c("Gene.stable.ID", "Gene.name", "log2FoldChange", "padj", "pvalue")]

  # Filter DEGs
  res_complete  <- res_df[complete.cases(res_df$padj), ]
  res_sig       <- res_complete[res_complete$padj < thresholds$padj, ]
  res_deg       <- res_sig[abs(res_sig$log2FoldChange) > thresholds$log2fc, ]

  up_genes   <- res_deg[res_deg$log2FoldChange > 0, "Gene.name"]
  down_genes <- res_deg[res_deg$log2FoldChange < 0, "Gene.name"]

  message(sprintf("Contrast %s vs %s: %d DEGs (%d up, %d down)",
                  contrast_group, reference_group,
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
plot_pca <- function(dds, timepoint_label) {
  vsd      <- vst(dds, blind = FALSE)
  pca_data <- plotPCA(vsd, intgroup = "RB_PPD_trajectory", returnData = TRUE)

  pct_var <- round(100 * attr(pca_data, "percentVar"), 1)

  color_map <- c(
    "Both"             = "#f8766d",
    "Controls"         = "#7cae00",
    "Postpartum_only"  = "#00bfc4",
    "Pregnancy_only"   = "#c77cff"
  )

  ggplot(pca_data, aes(x = PC1, y = PC2, color = RB_PPD_trajectory)) +
    geom_point(size = 3) +
    scale_color_manual(values = color_map, name = "Depression period") +
    xlab(paste0("PC1: ", pct_var[1], "% variance")) +
    ylab(paste0("PC2: ", pct_var[2], "% variance")) +
    ggtitle(paste("PCA -", timepoint_label)) +
    coord_fixed() +
    theme_bw()
}


plot_volcano <- function(raw_res, contrast_label, thresholds) {
  df <- as.data.frame(raw_res)
  df$log10padj <- -log10(df$padj)
  df$significant <- abs(df$log2FoldChange) > thresholds$log2fc &
                    df$padj < thresholds$padj & !is.na(df$padj)

  ggplot(df, aes(x = log2FoldChange, y = log10padj, color = significant)) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_color_manual(values = c("FALSE" = "#181C14", "TRUE" = "#f8766d"),
                       labels = c("Not significant", "Significant"),
                       name   = NULL) +
    geom_vline(xintercept = c(-thresholds$log2fc, thresholds$log2fc),
               linetype = "dashed", color = "blue") +
    geom_hline(yintercept = -log10(thresholds$padj),
               linetype = "dashed", color = "blue") +
    xlab("log2(Fold Change)") +
    ylab("-log10(adjusted p-value)") +
    ggtitle(paste("Volcano Plot -", contrast_label)) +
    theme_bw()
}


# --- GO ENRICHMENT ------------------------------------------------------------
run_go_enrichment <- function(gene_list, direction_label, ont = "MF") {
  if (length(gene_list) == 0) {
    message(sprintf("No genes for GO analysis (%s)", direction_label))
    return(NULL)
  }

  enrich_res <- enrichGO(
    gene          = gene_list,
    OrgDb         = "org.Hs.eg.db",
    keyType       = "SYMBOL",
    ont           = ont,
    pvalueCutoff  = 0.25,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.5,
    minGSSize     = 10,
    maxGSSize     = 500
  )

  if (!is.null(enrich_res) && nrow(enrich_res@result) > 0) {
    print(dotplot(enrich_res) + ggtitle(paste("GO", ont, "-", direction_label)))
  } else {
    message(sprintf("No significant GO terms for %s", direction_label))
  }

  return(enrich_res)
}


# --- SAVE RESULTS -------------------------------------------------------------
save_results <- function(results_list, contrast_group, reference_group,
                         timepoint, results_dir) {
  if (is.null(results_dir) || results_dir == "") {
    stop("results_dir is not set. Please update PATHS$results_dir in the configuration block.")
  }
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  label <- paste0(contrast_group, "_vs_", reference_group, "_", timepoint)

  # Full DEA table MANTRA (drug repurposing)
  write.table(
    results_list$full,
    file      = file.path(results_dir, paste0(label, "_full.txt")),
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE
  )

  # Filtered DEGs only
  write.table(
    results_list$sig,
    file      = file.path(results_dir, paste0(label, "_DEGs.txt")),
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE
  )

  # Gene lists 
  write.table(
    data.frame(gene = results_list$up),
    file      = file.path(results_dir, paste0(label, "_upregulated.txt")),
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE
  )

  write.table(
    data.frame(gene = results_list$down),
    file      = file.path(results_dir, paste0(label, "_downregulated.txt")),
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE
  )

  message(sprintf("Results saved: %s", label))
}


# --- MAIN WRAPPER FUNCTION ----------------------------------------------------
#' Run full DEA pipeline for one contrast at one timepoint
#'
#' @param timepoint    Character: "pp" (post-partum) or "v38" (pre-partum)
#' @param contrast     Character: depression group ("Both", "Pregnancy_only", "Postpartum_only")
#' @param reference    Character: reference group (default "Controls")
#' @param data_list    List returned by load_data()
#' @param paths        PATHS config list
#' @param filter       FILTER config list
#' @param thresholds   DEG_THRESHOLDS config list
#' @param outliers     OUTLIER_SAMPLES vector

run_dea_pipeline <- function(timepoint,
                             contrast,
                             reference  = "Controls",
                             data_list,
                             cfg_paths      = paths,
                             cfg_filter     = filter,
                             cfg_thresholds = deg_thresholds,
                             cfg_outliers   = outlier_samples) {

  timepoint_label <- ifelse(timepoint == "pp", "Post-partum", "Pre-partum")
  contrast_label  <- paste(contrast, "vs", reference, "-", timepoint_label)
  message("\n", paste(rep("=", 60), collapse = ""))
  message("Running: ", contrast_label)
  message(paste(rep("=", 60), collapse = ""))

  # 1. Preprocess metadata for this timepoint
  meta_sub <- preprocess_metadata(data_list$metadata, timepoint, cfg_outliers)

  # 2. Run DESeq2
  dds <- run_deseq2(data_list$counts, meta_sub, cfg_filter)

  # 3. Extract results for this contrast
  res <- extract_results(dds, contrast, reference, data_list$annotation, cfg_thresholds)

  # 4. Save results
  save_results(res, contrast, reference, timepoint, cfg_paths$results_dir)

  # 5. Plots
  print(plot_pca(dds, timepoint_label))
  print(plot_volcano(res$raw_res, contrast_label, cfg_thresholds))

  # 6. GO enrichment
  go_up   <- run_go_enrichment(res$up,   paste("Upregulated -",   contrast_label))
  go_down <- run_go_enrichment(res$down, paste("Downregulated -", contrast_label))

  return(list(dds = dds, results = res, go_up = go_up, go_down = go_down))
}


# =============================================================================
# RUN ALL 6 CONTRASTS
# =============================================================================

# Load data once
data_list <- load_data(paths)

# Define all contrasts: list(timepoint, contrast_group)
contrasts_to_run <- list(
  list(timepoint = "v38", contrast = "Both"),
  list(timepoint = "v38", contrast = "Pregnancy_only"),
  list(timepoint = "v38", contrast = "Postpartum_only"),
  list(timepoint = "pp",  contrast = "Both"),
  list(timepoint = "pp",  contrast = "Pregnancy_only"),
  list(timepoint = "pp",  contrast = "Postpartum_only")
)

# Run all contrasts and store results
all_results <- lapply(contrasts_to_run, function(x) {
  run_dea_pipeline(
    timepoint = x$timepoint,
    contrast  = x$contrast,
    data_list = data_list
  )
})

# Name the results list for easy access
names(all_results) <- sapply(contrasts_to_run, function(x) {
  paste0(x$contrast, "_", x$timepoint)
})

