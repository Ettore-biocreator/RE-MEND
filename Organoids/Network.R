# =============================================================================
# RE-MEND PROJECT - Task 3: MANTRA Results Network Visualization
# Description: Network plot summarizing MANTRA drug repurposing results for
#              Hormones experiment (CTL04 + CTL08) and EDC experiment (CTL04).
#              Edge weights are based on MANTRA similarity scores.
# =============================================================================


# --- LIBRARIES ----------------------------------------------------------------
library(igraph)
library(RColorBrewer)


# --- CONFIGURATION ------------------------------------------------------------
paths <- list(
  hormones    = "path/to/Combined_dataset.txt",
  edc_ctl04   = "path/to/Dataset_EDC_CTL04.txt",
  output      = "results/task3_network/network_plot.png"
)

# Edge weight transformation parameters
power_exp     <- 14       # exponent for power transformation of edge weights
layout_niter  <- 500      # number of iterations for Fruchterman-Reingold layout

# Visual parameters
node_size     <- 10
label_cex     <- 0.6
label_color   <- "darkblue"

# Color palette for Sex (CTL04 = female, CTL08 = male)
sex_palette <- c("CTL04" = "#FD8B51", "CTL08" = "#73EC8B")

# Shape palette for Experiment
shape_palette <- c("Hormones" = "circle", "EDCs" = "square")


# --- LOAD AND PREPARE DATA ----------------------------------------------------
prepare_data <- function(paths) {
  # --- Hormones dataset ---
  hormones <- read.delim2(paths$hormones)

  # Assign Agonist/Inhibitor class based on Model column
  hormones$Class <- "AG"
  hormones[hormones$Model == "B", ]$Class <- "INH"

  # Assign direction type based on Model column
  hormones$Type <- "Inv"
  hormones[hormones$Model %in% c("A", "C"), ]$Type <- "Dir"

  # Build Column label: Type_Column_Class
  hormones$Column <- paste(hormones$Type, hormones$Column, sep = "_")
  hormones$Column <- paste(hormones$Column, hormones$Class, sep = "_")

  # Keep only relevant columns and add experiment metadata
  hormones <- hormones[, c(1:3)]
  hormones$Experiment <- "Hormones"
  hormones$Sex        <- c(rep("CTL04", 11), rep("CTL08", 11))
  hormones$Column     <- paste(hormones$Column, hormones$Sex, sep = "_")

  # --- EDC dataset (CTL04 only) ---
  edc_ctl04 <- read.delim2(paths$edc_ctl04)
  edc_ctl04$Column     <- sub(".*_", "", edc_ctl04$Column)
  edc_ctl04$Experiment <- "EDCs"
  edc_ctl04$Sex        <- "CTL04"

  # --- Combine ---
  combined <- rbind(edc_ctl04, hormones)
  combined$Value <- as.numeric(combined$Value)

  return(combined)
}


# --- BUILD NETWORK ------------------------------------------------------------
build_network <- function(data, power_exp, layout_niter) {
  # Create graph from edge list
  network <- graph_from_data_frame(d = data, directed = FALSE)

  # Power-transform and normalize edge weights
  values        <- data$Value
  power_w       <- values ^ power_exp
  normalized    <- (power_w - min(power_w)) / (max(power_w) - min(power_w))
  normalized    <- normalized + 1e-6  # avoid exact zeros

  E(network)$weight <- normalized

  # Compute inverted normalized weights for edge width (high similarity = thinner edge)
  inv_w      <- max(normalized) - normalized + min(normalized)
  edge_widths <- (inv_w - min(inv_w)) / (max(inv_w) - min(inv_w)) * 2

  E(network)$edge_width <- edge_widths

  # Node metadata
  node_meta <- unique(data.frame(
    name       = data$Column,
    Experiment = data$Experiment,
    Sex        = data$Sex,
    stringsAsFactors = FALSE
  ))

  # Handle any nodes present in graph but missing from metadata
  all_nodes     <- V(network)$name
  missing_nodes <- setdiff(all_nodes, node_meta$name)
  if (length(missing_nodes) > 0) {
    node_meta <- rbind(node_meta, data.frame(
      name       = missing_nodes,
      Experiment = NA,
      Sex        = NA
    ))
  }

  node_meta <- node_meta[match(all_nodes, node_meta$name), ]

  # Compute layout
  layout_pos <- layout_with_fr(network, niter = layout_niter) * 100

  return(list(network = network, node_meta = node_meta, layout = layout_pos))
}


# --- PLOT NETWORK -------------------------------------------------------------
plot_network <- function(network_obj, sex_palette, shape_palette,
                         node_size, label_cex, label_color, output_path) {
  network   <- network_obj$network
  node_meta <- network_obj$node_meta
  layout    <- network_obj$layout

  # Node colors by Sex
  node_colors <- sex_palette[as.character(node_meta$Sex)]
  node_colors[is.na(node_colors)] <- "white"

  # Node shapes by Experiment
  node_shapes <- shape_palette[as.character(node_meta$Experiment)]
  node_shapes[is.na(node_shapes)] <- "circle"

  sex_levels <- names(sex_palette)
  exp_levels <- names(shape_palette)

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

  png(output_path, width = 2500, height = 2500, res = 300)

  plot(
    network,
    layout             = layout,
    vertex.color       = node_colors,
    vertex.shape       = node_shapes,
    vertex.size        = node_size,
    edge.width         = E(network)$edge_width,
    edge.color         = "black",
    vertex.label       = V(network)$name,
    vertex.label.cex   = label_cex,
    vertex.label.color = label_color,
    vertex.label.family = "sans"
  )

  legend("bottomleft",  legend = sex_levels, col = sex_palette,
         pch = 16, pt.cex = 2, title = "Model", bty = "n")
  legend("bottomright", legend = exp_levels, pch = c(21, 22),
         pt.bg = "grey", pt.cex = 2, title = "Experiment", bty = "n")

  dev.off()
  message(sprintf("Network plot saved: %s", output_path))
}


# =============================================================================
# MAIN
# =============================================================================

data    <- prepare_data(paths)
net_obj <- build_network(data, power_exp, layout_niter)
plot_network(net_obj, sex_palette, shape_palette,
             node_size, label_cex, label_color, paths$output)
