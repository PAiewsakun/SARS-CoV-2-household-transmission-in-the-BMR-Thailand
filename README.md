# Integrating clinical, demographic, epidemiological, and viral genome data to investigate SARS-CoV-2 household transmission in the Bangkok metropolitan region, Thailand, 2022-2024

**Date:** 26/08/2026

## Description
This repository contains the R scripts, input data, and generated outputs used to produce the figures, tables, and statistical analyses presented in the manuscript: *"Integrating clinical, demographic, epidemiological, and viral genome data to investigate SARS-CoV-2 household transmission in the Bangkok metropolitan region, Thailand, 2022-2024"* (Aiewsakun et al., In Preparation).

## Repository contents
```
│   README.md
│   LICENSE
│
└───code
    │   1 Table S4 Figure S1 Symptom level clustering analysis.R
    │   2 Figure 1 Overview of the study cohort.R
    │   3 Figure 2 Major lineage temporal dynamics.R
    │   4 Table S5 and Table 1 Overview of the study cohort at the household lv.R
    │   5 Table 2 Iden factors assoc with endpoint COVID-19 status.R
    │   6 Table 3 Figure S2 Iden factors assoc with COVID-19 symp lv.R
    │   7 Figure 3 S3 S4 Table S9 Transmission clusters analysis.R
    │   8 Hidden comm acqu rates and attack rates.R
    │   miscellaneous.R
    │
    ├───data
    │       Th_map
    │       Data S1.ML phylogeny.treefile
    │       Data S2.Time-calibrated phylogeny.treefile
    │       Table S1.Participant metadata.txt
    │       Table S2.Household metadata.txt
    │       Table S3.Weekly counts of enrolled households COVID-19 cases and sequences.txt
    │       Table S6.Weekly sequence counts by major lineage.txt
    │       Table S7.Sequence metadata.txt
    │       Table S8.Root-to-tip regression analysis to detect outliers.xlsx
    │
    └───out
            clustering_robustness.txt
            clustering_score.txt
            Figure 1.Overview of the study cohort.raw.png
            Figure 1.Overview of the study cohort.raw.svg
            Figure 2.Genome lineage diversity and temporal dynamics of SARS-CoV-2 sequences.raw.png
            Figure 2.Genome lineage diversity and temporal dynamics of SARS-CoV-2 sequences.raw.svg
            Figure 3.Phylogenetic analysis and sequence clustering.raw.png
            Figure 3.Phylogenetic analysis and sequence clustering.raw.svg
            Figure S1.Symptom profile clustering analysis.raw.png
            Figure S1.Symptom profile clustering analysis.raw.svg
            Figure S2.Symptom level vs ORF1ab Ct value.raw.png
            Figure S2.Symptom level vs ORF1ab Ct value.raw.svg
            Figure S3.Soft-pruned time-calibrated phylogeny and sequence clustering.raw.png
            Figure S3.Soft-pruned time-calibrated phylogeny and sequence clustering.raw.svg
            Figure S4.Distribution of SARS-CoV-2 sequences across families and sequence clusters.raw.png
            Figure S4.Distribution of SARS-CoV-2 sequences across families and sequence clusters.raw.svg
            Table 1.Overview of the study cohort at HH level by HH type.2 lvs.raw.docx
            Table 2.Descriptive statistics and analyses of factors associated with endpoint positivity.raw.docx
            Table 3.Descriptive statistics and analyses of factors associated with symptom level.raw.docx
            Table S4.Severity clustering analysis results.raw.txt
            Table S5.Overview of the study cohort at HH level by HH type.3 lvs.raw.docx
            Table S9.Sequence cluster summary statistics.raw.txt
```

## Instructions for Use
1. Download or clone this repository.
2. Open R and set the working directory to `codes/`.
3. Run the R scripts individually to generate the results.
   * All required input data is located in the `codes/data/` directory.
   * All generated outputs (figures, tables, etc.) will be saved to the `codes/out/` directory.
   * `miscellaneous.R` contains global variables. It is sourced at the beginning of each of the main analysis scripts, and does not need to be executed independently.
   * Because some analyses involve stochastic procedures (e.g., simulation-based assignment of primary cases and clustering analyses), random seeds were set within the corresponding scripts to ensure reproducibility.


## System requirements & dependencies
The analyses were performed using **R version 4.4.2**. The following R packages are required to run the scripts. Package versions correspond to the environment used to generate the published results. We thus recommend installing these specific versions to ensure full reproducibility:

### Data manipulation & utilities

* `readxl` (v1.4.5)
* `tidyverse` (v2.0.0)
* `dplyr` (v1.2.0)
* `purrr` (v1.2.1)
* `lubridate` (v1.9.5)
* `fuzzyjoin` (v0.1.8)

### Dimension reduction & clustering analysis

* `umap` (v0.2.10.0)
* `cluster` (v2.1.6)
* `mclust` (v6.1.2)
* `aricode` (v1.0.3)

### Phylogenetics & network analysis

* `ape` (v5.8.1)
* `igraph` (v2.2.2)

### Statistics

* `lme4` (v1.1.38)
* `glmmTMB` (v1.1.14)
* `ordinal` (v2025.12.29)
* `gamlss` (v5.5.0)
* `buildmer` (v2.12)
* `emmeans` (v2.0.1)
* `performance` (v0.16.0)
* `gofcat` (v0.1.2)
* `gtsummary` (v2.5.0)
* `rstatix` (v0.7.3)

### Visualisation

* `ggplot2` (v4.0.2)
* `GGally` (v2.4.0)
* `ggpmisc` (v0.6.3)
* `ggforce` (v0.5.0)
* `ggtext` (v0.1.2)
* `ggh4x` (v0.3.1)
* `ggtree` (v3.14.0)
* `tidytree` (v0.4.7)
* `ggnewscale` (v0.5.2)
* `sf` (v1.1.0)
* `ComplexHeatmap` (v2.22.0)
* `cowplot` (v1.2.0)
* `patchwork` (v1.3.2)

## Data availability
All raw sequencing data are available from the National Center for Biotechnology Information (NCBI) under accession numbers PV890700-PV890817 and PX491811-PX491820. All other data required to reproduce the analyses are provided in the `codes/data/` directory.

## Citation
If you use these scripts or find them helpful for your research, please cite:
> Aiewsakun, P. et al. (In Preparation). Integrating clinical, demographic, epidemiological, and viral genome data to investigate SARS-CoV-2 household transmission in the Bangkok metropolitan region, Thailand, 2022-2024.

## License
This project is licensed under the GNU General Public License v3.0 - see the `LICENSE` file for the full legal text.