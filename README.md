# RE-MEND: Transcriptomic Analysis of Postpartum Depression

This repository contains the bioinformatic analysis scripts developed as part of the **RE-MEND European Project**, which investigates the molecular basis of postpartum depression (PPD) through transcriptomic profiling and drug repurposing.

> ⚠️ **Note:** These scripts represent a general framework of the analyses performed. File paths, sample names, and parameters may need to be adapted to match your specific dataset and directory structure before running.

---

## Repository Structure

```
RE-MEND/
│
├── Basic/
│   ├── Task1_DEA_Depression.R
│   ├── Task2_DEA_Healthy.R
│   └── Task2b_DEA_Healthy_Paired.R
│
├── Organoids/
│   ├── Task3_Hormones_DEA.R
│   └── Task3_Network_Visualization.R
│
└── README.md
```

---

## Analysis Overview

### Basic/

#### Task 1 — Differential Expression Analysis: Depressed vs Healthy Women
RNA-seq differential expression analysis (DESeq2) on ~340 blood samples collected at two timepoints (pre-partum: `v38`, post-partum: `pp`). Depression status is defined using EPDS scores and classified into four trajectories: **Both**, **Pregnancy_only**, **Postpartum_only**, and **Controls**. Six contrasts are tested (3 per timepoint), correcting for age and BMI. DEG lists are used as input for drug repurposing via MANTRA.

#### Task 2 — Differential Expression Analysis: Pre vs Post-partum in Healthy Women
DESeq2 analysis comparing pre-partum and post-partum timepoints exclusively in healthy control women, excluding SSRI users. Corrects for age and BMI.

#### Task 2b — Paired Analysis: Pre vs Post-partum in Healthy Women
DESeq2 analysis restricted to paired samples from healthy control women (both timepoints available for the same subject). Simplified design without covariate correction due to the paired structure.

---

### Organoids/

#### Task 3 — Hormones Experiment: DEA in Brain Organoids
DESeq2 analysis on cortical brain organoids (female: **CTL04**, male: **CTL08**) treated with agonists and inhibitors of different hormone classes (e.g. Androgens, Glucocorticoids). Each treatment is compared against DMSO control. Full ranked gene lists are produced as input for drug repurposing via MANTRA.

#### Task 3 — Network Visualization of MANTRA Results
Network visualization (igraph) summarizing drug repurposing results from MANTRA across the Hormones experiment (CTL04 + CTL08) and EDC experiment (CTL04 only). Nodes represent experimental conditions, edge weights reflect MANTRA similarity scores.

---

## Drug Repurposing: MANTRA

Drug repurposing analyses were performed using **MANTRA** (Mode of Action by NeTwoRk Analysis).  
The MANTRA pipeline is **not included** in this repository as it is a proprietary tool.  
MANTRA can be accessed and executed at: [https://mantra.tigem.it/](https://mantra.tigem.it/)
---

## Dependencies

All scripts are written in **R**. The following packages are required:

| Package | Version | Use |
|---|---|---|
| DESeq2 | ≥ 1.38 | Differential expression analysis |
| dplyr | ≥ 1.1 | Data manipulation |
| ggplot2 | ≥ 3.4 | Visualization |
| openxlsx | ≥ 4.2 | Reading Excel files |
| clusterProfiler | ≥ 4.6 | GO enrichment (Task 1) |
| org.Hs.eg.db | ≥ 3.16 | Human gene annotation (Task 1) |
| igraph | ≥ 1.4 | Network analysis and visualization |
| RColorBrewer | ≥ 1.1 | Color palettes |

Install all dependencies with:

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db"))
install.packages(c("dplyr", "ggplot2", "openxlsx", "igraph", "RColorBrewer"))
```

---

## Citation

If you use these scripts, please cite the RE-MEND project and the relevant tools:

- **DESeq2:** Love MI, Huber W, Anders S. *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biology, 2014.
- **MANTRA:** Carrella D, et al. *Mantra 2.0: an online collaborative resource for drug mode of action and repurposing by network analysis.* Bioinformatics, 2014.
- **clusterProfiler:** Wu T, et al. *clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.* Innovation, 2021.

---

## Contact

For questions regarding the analysis pipeline, please contact the RE-MEND bioinformatics team.
