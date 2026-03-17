# =============================================================================
# RE-MEND PROJECT - Task 3: Organoids - Hormones Experiment
# Description: DEA between hormone treatments (agonist/inhibitor) and DMSO
#              control in female (CTL04) and male (CTL08) brain organoids.
#              Full ranked gene lists are produced as input for drug repurposing.
# =============================================================================


# --- LIBRARIES ----------------------------------------------------------------
library(DESeq2)
library(dplyr)
library(ggplot2)
library(openxlsx)


# --- CONFIGURATION ------------------------------------------------------------
# Each model has its own count matrix, sample sheet, and filtering threshold

models <- list(
  CTL04 = list(
    counts          = "path/to/Counts_ref_CTL04.txt",
    metadata        = "path/to/SampleSheet_ref_CTL04.xlsx",
    min_samples     = 4,
    results_dir     = "results/task3_hormones/CTL04/"
  ),
  CTL08 = list(
    counts          = "path/to/Counts_ref_CTL08.txt",
    metadata        = "path/to/SampleSheet_ref_CTL08.xlsx",
    min_samples     = 4,
    results_dir     = "results/task3_hormones/CTL08/"
  )
)

annotation_path <- NULL  # not needed: count matrix already contains gene names

# DESeq2 filtering threshold (shared across models)
min_counts <- 10

# Reference condition
reference <- "DMSO"

# Conditions to exclude from contrasts
exclude_conditions <- c("DMSO", "CTL")


# --- LOAD DATA ----------------------------------------------------------------
load_data <- function(model_cfg) {
  # Count matrix (first column = gene names, already annotated)
  counts <- read.delim2(model_cfg$counts)
  rownames(counts) <- counts[, 1]
  counts <- counts[, -1]

  # Metadata
  metadata <- read.xlsx(model_cfg$metadata)
  colnames(metadata)[1] <- "Samples"

  message(sprintf("Loaded %d genes x %d samples", nrow(counts), ncol(counts)))

  return(list(counts = counts, metadata = metadata))
}


# --- RUN DESEQ2 ---------------------------------------------------------------
run_deseq2 <- function(counts, metadata, min_counts, min_samples) {
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData   = metadata,
    design    = ~ Condition
  )

  # Low-count filtering
  keep <- rowSums(counts(dds) >= min_counts) >= min_samples
  dds  <- dds[keep, ]

  message(sprintf("Genes retained after filtering: %d", sum(keep)))

  dds <- DESeq(dds)
  return(dds)
}


# --- EXTRACT AND ANNOTATE RESULTS ---------------------------------------------
extract_results <- function(dds, contrast_cond, reference) {
  res         <- results(dds, contrast = c("Condition", contrast_cond, reference))
  res_ordered <- res[order(res$stat, decreasing = TRUE), ]  # rank by stat for MANTRA

  res_df <- as.data.frame(res_ordered)
  res_df$Gene.name <- rownames(res_df)
  res_df <- res_df %>%
    select(Gene.name, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

  return(res_df)
}


# --- VOLCANO PLOT -------------------------------------------------------------
plot_volcano <- function(res_df, contrast_label) {
  df <- res_df
  df$log10padj   <- -log10(df$padj)
  df$significant <- abs(df$log2FoldChange) > 1 &
                    !is.na(df$padj) & df$padj < 0.05

  ggplot(df, aes(x = log2FoldChange, y = log10padj, color = significant)) +
    geom_point(size = 0.8, alpha = 0.7) +
    scale_color_manual(
      values = c("FALSE" = "#181C14", "TRUE" = "#f8766d"),
      labels = c("Not significant", "Significant"),
      name   = NULL
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
    xlab("log2(Fold Change)") +
    ylab("-log10(adjusted p-value)") +
    ggtitle(paste("Volcano Plot -", contrast_label)) +
    theme_bw()
}


# --- PCA PLOT -----------------------------------------------------------------
plot_pca <- function(dds, model_name) {
  vsd      <- vst(dds, blind = FALSE)
  pca_data <- plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
  pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

  ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", pct_var[1], "% variance")) +
    ylab(paste0("PC2: ", pct_var[2], "% variance")) +
    ggtitle(paste("PCA -", model_name)) +
    theme_bw()
}


# --- SAVE RESULTS -------------------------------------------------------------
save_results <- function(res_df, contrast_cond, results_dir) {
  if (is.null(results_dir) || results_dir == "") {
    stop("results_dir is not set. Please update the model configuration.")
  }
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  write.table(
    res_df,
    file      = file.path(results_dir, paste0(contrast_cond, "_vs_DMSO.tsv")),
    sep       = "\t",
    row.names = FALSE,
    quote     = FALSE
  )
}


# --- MAIN WRAPPER FUNCTION ----------------------------------------------------
#' Run full DEA pipeline for one organoid model
#'
#' @param model_name  Character: "CTL04" or "CTL08"
#' @param model_cfg   List: model-specific configuration from models list

run_organoid_pipeline <- function(model_name, model_cfg) {
  message("\n", paste(rep("=", 60), collapse = ""))
  message("Model: ", model_name)
  message(paste(rep("=", 60), collapse = ""))

  # 1. Load data
  data_list <- load_data(model_cfg)

  # 2. Run DESeq2
  dds <- run_deseq2(data_list$counts, data_list$metadata,
                    min_counts, model_cfg$min_samples)

  # 3. PCA
  print(plot_pca(dds, model_name))

  # 4. Get conditions to test (exclude DMSO and CTL)
  conds <- unique(data_list$metadata$Condition)
  conds <- conds[!(conds %in% exclude_conditions)]
  message(sprintf("Conditions to test: %s", paste(conds, collapse = ", ")))

  # 5. Loop over conditions
  results_list <- lapply(conds, function(cond) {
    message(sprintf("  Contrast: %s vs %s", cond, reference))

    res_df <- extract_results(dds, cond, reference)
    save_results(res_df, cond, model_cfg$results_dir)
    print(plot_volcano(res_df, paste(cond, "vs", reference, "-", model_name)))

    return(res_df)
  })

  names(results_list) <- conds
  return(results_list)
}


# =============================================================================
# RUN ALL MODELS
# =============================================================================

all_results <- lapply(names(models), function(model_name) {
  run_organoid_pipeline(model_name, models[[model_name]])
})

names(all_results) <- names(models)
