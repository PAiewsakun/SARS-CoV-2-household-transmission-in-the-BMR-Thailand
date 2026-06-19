# ###############################
# Load required libraries
# ###############################
library(dplyr)
library(readxl)
library(ggplot2) #for defining global plotting theme my_theme

# ###############################
# Paths to data
# ###############################
participant_metadata_filepath <- file.path("data/Table S1.Participant metadata.txt")
HH_metadata_filepath <- file.path("data/Table S2.Household metadata.txt")
weekly_HH_case_seq_count_dat_filepath <- file.path("data/Table S3.Weekly counts of enrolled households COVID-19 cases and sequences.txt")

weekly_seq_count_by_lineage_dat_filepath <- file.path("data/Table S6.Weekly sequence counts by major lineage.txt")

tree_filepath <- file.path("data/Data S2.Time-calibrated phylogeny.treefile")
seq_metadata_filepath <- file.path("data/Table S7.Sequence metadata.txt")

Th_prov_map_filepath <- file.path("data/Th_map/tha_admbnda_adm1_rtsd_20220121.shp")
Th_dist_map_filepath <- file.path("data/Th_map/tha_admbnda_adm2_rtsd_20220121.shp")

# ###############################
# Paths to output files
# ###############################
Table_1_file <- file.path("out/Table 1.Overview of the study cohort at HH level by HH type.2 lvs.raw.docx")
Table_2_file <- file.path("out/Table 2.Descriptive statistics and analyses of factors associated with endpoint positivity.raw.docx")
Table_3_file <- file.path("out/Table 3.Descriptive statistics and analyses of factors associated with symptom level.raw.docx")

Table_S4_file <- file.path("out/Table S4.Severity clustering analysis results.raw.txt")
Table_S5_file <- file.path("out/Table S5.Overview of the study cohort at HH level by HH type.3 lvs.raw.docx")
Table_S9_file <- file.path("out/Table S9.Sequence cluster summary statistics.raw.txt")

Figure_1_file_png <- file.path("out/Figure 1.Overview of the study cohort.raw.png")
Figure_1_file_svg <- file.path("out/Figure 1.Overview of the study cohort.raw.svg")

Figure_2_file_png <- file.path("out/Figure 2.Genome lineage diversity and temporal dynamics of SARS-CoV-2 sequences.raw.png")
Figure_2_file_svg <- file.path("out/Figure 2.Genome lineage diversity and temporal dynamics of SARS-CoV-2 sequences.raw.svg")

Figure_3_file_png <- file.path("out/Figure 3.Phylogenetic analysis and sequence clustering.raw.png")
Figure_3_file_svg <- file.path("out/Figure 3.Phylogenetic analysis and sequence clustering.raw.svg")

Figure_S1_file_png <- file.path("out/Figure S1.Symptom profile clustering analysis.raw.png")
Figure_S1_file_svg <- file.path("out/Figure S1.Symptom profile clustering analysis.raw.svg")

Figure_S2_file_png <- file.path("out/Figure S2.Symptom level vs ORF1ab Ct value.raw.png")
Figure_S2_file_svg <- file.path("out/Figure S2.Symptom level vs ORF1ab Ct value.raw.svg")

Figure_S3_file_png <- file.path("out/Figure S3.Soft-pruned time-calibrated phylogeny and sequence clustering.raw.png")
Figure_S3_file_svg <- file.path("out/Figure S3.Soft-pruned time-calibrated phylogeny and sequence clustering.raw.svg")

Figure_S4_file_png <- file.path("out/Figure S4.Distribution of SARS-CoV-2 sequences across families and sequence clusters.raw.png")
Figure_S4_file_svg <- file.path("out/Figure S4.Distribution of SARS-CoV-2 sequences across families and sequence clusters.raw.svg")

# ###############################
# Global vars
# ###############################
BMR_prov_list <- c("Nakhon Pathom", "Nonthaburi", "Pathum Thani", "Samut Sakhon", "Bangkok", "Samut Prakan")
BMR_prov_and_beyond_list <- c(BMR_prov_list, c("Samut Songkhram", "Ratchaburi", "Kanchanaburi", "Suphan Buri", "Phra Nakhon Si Ayutthaya", "Saraburi", "Nakhon Nayok", "Chachoengsao", "Chon Buri"))

# ###############################
# Global plotting theme
# ###############################
my_theme <- theme(
    plot.title = element_text(size = 8, face = "bold"),
    plot.margin = unit(c(0, 0.05, 0, 0.05), "cm"),
    
    axis.title = element_text(size = 6),
    axis.text = element_text(size = 6),
    
    legend.position = "inside", legend.justification = c(1, 1), legend.position.inside = c(1, 1),
    legend.background = element_blank(),
    
    legend.key = element_rect(fill = NA, colour = NA),
    legend.key.size = unit(0.25, 'cm'),
    legend.key.spacing.x = unit(0, 'cm'),
    
    legend.title = element_text(size = 6, face = "bold", margin = margin(t = 0, b = 1)),
    legend.text = element_text(size = 6, margin = margin(r = 0, l = 0)),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0), # Remove all margins
    
    strip.placement.x = 'outside',
    strip.text = element_text(size = 6),
    strip.background.x = element_blank(),
    
  )

# ###############################
# Col vars
# ###############################
BMR_cols <- c(
  "BMR" = "#78c5d6", 
  "non-BMR" = "#459ba8"
)

prov_cols <- c(
  "Nakhon Pathom" = "#1b9e77",
  "Nonthaburi" = "#d95f02",
  "Pathum Thani" = "#7570b3",
  "Samut Sakhon" = "#e7298a",
  "Bangkok" = "#66a61e",
  "Samut Prakan" = "#e6ab02"
)

overall_res_enrol_cols <- c(
  "Pos" = "#ff3b30", 
  "Neg" = "#43BF71", 
  "Unconfirmed" = "grey85"
)

overall_res_end_cols <- c(
  "Pos" = "#ff3b30", 
  "Neg" = "#43BF71", 
  "Unconfirmed" = "grey85"
)

severity_cols <- c(
  "Minimally/(a)symptomatic" = "#dffe00",
  "Moderately symptomatic" = "#ff9500",
  "Highly symptomatic" = "#ff3b30"
)

overall_res_end_with_symp_severity_cols <- c(
  "Pos: no symptom data" = "#FFC0CB",
  "Pos: highly symptomatic" = "#ff3b30",
  "Pos: moderately symptomatic" = "#ff9500",
  "Pos: minimally/(a)symptomatic" = "#dffe00", #"#fce205",
  "Neg" = "#43BF71",
  "Unconfirmed" = "grey85"
)

major_lineage_cols <- c(
  "BA.5" = "#BBDF27", #viridis(1, begin = 0.9, end = 1),
  "BA.2.75" = "#43BF71", #viridis(1, begin = 0.7, end = 0.8),
  "XBB" = "#21908C", #viridis(1, begin = 0.5, end = 0.6),
  "XBL" = "#414487", #viridis(1, begin = 0.2, end = 0.5),
  "XBB or XBL" = "#35608D", #viridis(1, begin = 0.3, end = 0.5),
  "BA.2.86" = "#440154", #viridis(1, begin = 0, end = 0.2),
  "Others" = "gray85"
)

age_cat_cols <- c(
  "Adult" = "#F39C6E",
  "Child" = "#F2DB92"
)

sex_cols <- c(
  "Male" = "#8CD3FF",
  "Female" = "#FFC8DC"
)

HH_type_cols <- c(
  "Adult-first household" = "#92809A",
  "Child-first household" = "#BC8BA1" # "#7A97B7"
)

