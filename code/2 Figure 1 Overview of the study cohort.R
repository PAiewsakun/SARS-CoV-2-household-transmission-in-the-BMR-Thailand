# ###############################
# Load global variables
# ###############################
source(file.path("miscellaneous.R"))

# ###############################
# Load required libraries
# ###############################
library(tidyverse)
library(lubridate)
library(fuzzyjoin)  # for non-equi join

library(sf) #st_read

library(ggplot2)
library(GGally)
library(ggpmisc)
library(cowplot)
library(patchwork)

# ###############################
# Load data
# ###############################
# Participant metadata
# -------------------------------
participant_metadata <- read.delim(participant_metadata_filepath, header = TRUE, sep = "\t") %>%
  mutate(
    overall_res_enrol = factor(overall_res_enrol, levels = names(overall_res_enrol_cols)),
    overall_res_enrol_strict = factor(overall_res_enrol_strict, levels = names(overall_res_enrol_cols)),
    overall_res_end = factor(overall_res_end, levels = names(overall_res_end_cols)),
    overall_res_end_with_symp_severity = factor(overall_res_end_with_symp_severity, levels = names(overall_res_end_with_symp_severity_cols)),
    
    group_wrt_index_case = factor(group_wrt_index_case, levels = c("Index case", "Co-index case", "Pre-index case / additional case identified at enrolment", "Post-index case / incident case", "Uninfected contact", "Unconfirmed status")),
    symp_cluster_res = factor(symp_cluster_res, levels = c("Minimally/(a)symptomatic", "Moderately symptomatic", "Highly symptomatic")), 
    time_to_last_vac_cat = factor(time_to_last_vac_cat, levels = c( "Never", ">1 yr", "0.5-1 yr", "<0.5 yr")),
    time_to_last_prev_COVID_cat = factor(time_to_last_prev_COVID_cat, levels = c( "Never", ">1 yr", "0.5-1 yr", "<0.5 yr")),
    
    major_lineage = factor(major_lineage, levels = names(major_lineage_cols))
  )

# Household metadata
# -------------------------------
HH_metadata <- read.delim(HH_metadata_filepath, header = TRUE, sep = "\t")

# Weekly counts of enrolled households, COVID-19 cases, GISAID sequences
# -------------------------------
weekly_HH_case_seq_count_dat <- read.delim(weekly_HH_case_seq_count_dat_filepath, header = TRUE, sep = "\t") 

# ###############################
# Plot Figure 1A: Geographic distribution of the recruited households
# ###############################
# Load map data
# -------------------------------
Th_prov_map <- st_read(dsn = Th_prov_map_filepath) %>% as.data.frame %>% 
  select(geometry, ADM1_EN) %>% st_as_sf(na.fail = FALSE)

BMR_prov_map <- Th_prov_map %>% filter(ADM1_EN %in% BMR_prov_list) 

BMR_and_surrounding_prov_map <- Th_prov_map %>%
  filter(ADM1_EN %in% BMR_prov_and_beyond_list)

BMR_dist_map <- st_read(dsn = Th_dist_map_filepath) %>% as.data.frame %>% 
  select(geometry, ADM2_EN, ADM1_EN) %>% filter(ADM1_EN %in% BMR_prov_list) %>%
  st_as_sf(na.fail = FALSE)

# Plot Thailand map (province lv)
# -------------------------------
# Main plot
Th_prov_map_plot <- ggplot() +
  geom_sf(data = Th_prov_map, lwd = 0.1, fill = NA, color = "black") +
  geom_sf(data = BMR_prov_map, aes(fill = ADM1_EN), linewidth = 0.1, color = "black", alpha = 0.7) +
  scale_fill_manual(values = prov_cols) +
  annotate(
    geom = "rect", 
    xmin = 99.7, xmax = 101, ymin = 13.4, ymax = 14.3,
    fill = NA, colour = "black", linewidth = 0.5
  ) + 
  theme_void() + theme(legend.position = "none")

# Zoom-in plot of the BMR region
Th_prov_map_plot_zoomed_in <- Th_prov_map_plot + 
  coord_sf(xlim = c(99.7, 101), ylim = c(13.4, 14.3), expand = FALSE) +
  guides(colour = "none", x = "none", y = "none")

# Combine the main Thailand map with the BMR zoom-in plot
Th_prov_map_plot <- Th_prov_map_plot + 
  geom_plot_npc(
    data = tibble(x = 100, y = 0, plot = list(Th_prov_map_plot_zoomed_in)),
    aes(npcx = x, npcy = y, label = plot, vp.width = 0.55, vp.height = 0.55)
  )

# Plot BMR map (dist lv)
# -------------------------------
# Calculate the number of recruited households per dist using HH_metadata
HH_counts_by_dist <- HH_metadata %>%
  group_by(district, province) %>%
  summarise(num_HHs = n(), .groups = "drop")

# Join dist counts to the BMR dist map
HH_counts_by_dist <- left_join(st_centroid(BMR_dist_map), HH_counts_by_dist, by = c("ADM2_EN" = "district", "ADM1_EN" = "province"))

# Plot BMR map (dist level) with circles representing household counts
BMR_and_surrounding_prov_map_plot <- ggplot() +
  geom_sf(data = BMR_and_surrounding_prov_map, fill = NA, lwd = 0.1, color = "black") +
  geom_sf(data = BMR_prov_map, aes(fill = ADM1_EN), lwd = 1, color = "black", alpha = 0.7) + 
    scale_fill_manual(values = prov_cols, guide = "none") +
  geom_sf(data = BMR_dist_map, fill = NA, lwd = 0.1, color = "black") +
  
  # Add circles for dist counts
  geom_sf(data = HH_counts_by_dist, 
          aes(size = num_HHs, color = num_HHs), 
          fill = NA, alpha = 0.7) +
  scale_size_area(max_size = 7, guide = "none") +
  scale_color_continuous(low = "#0000FF", high = "#FF0000") + 
  guides(
    colour = guide_colourbar(
      direction = "horizontal", 
      title.position = "left",
      title = "Number of\nenrolled households",
      theme = theme(
        legend.key.width  = unit(2, "cm"),
        legend.key.height = unit(0.25, "cm"),
        legend.text.position = "bottom"
      )
    )
  ) +
  
  # Add text labels for HH counts
  geom_sf_text(data = HH_counts_by_dist, 
               aes(label = num_HHs), 
               size = 6/(72.27/25.4), color = "black") +
  coord_sf(xlim = c(99.7, 101), ylim = c(13.4, 14.3), expand = FALSE) +
  theme_void() +  
  theme(
    legend.position = "inside", legend.justification = c(1, 0), legend.position.inside = c(1, 0),
    legend.background = element_blank(),
    legend.title = element_text(size = 6, face = "bold", margin = margin(r = 1), hjust = 1),
    legend.text = element_text(size = 6, margin = margin(t = 1)),
  )

# Combine Thailand map (province lv) and BMR map (dist level)
# -------------------------------
Figure_1A <- Th_prov_map_plot + BMR_and_surrounding_prov_map_plot + plot_layout(nrow = 1) +
  plot_annotation(
    title = "Geographic distribution of enrolled households",
    theme = theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5))
  )

# ###############################
# Plot Figure 1B: Distribution of household sizes
# ###############################
Figure_1B <- HH_metadata %>%
  ggplot(aes(x = HH_size)) +
  geom_histogram(binwidth = 1, color = NA, alpha = 0.7) +
  scale_x_continuous(breaks = seq(0, 10, 1)) +
  labs(
    x = "Household size",
    y = "# of households",
    title = "Household size"
  ) +
  theme_bw() + my_theme

# ###############################
# Plot Figure 1C: Distribution of child-to-adult ratio in the households
# ###############################
Figure_1C <- HH_metadata %>%
  ggplot(aes(x = child_to_adult_ratio)) +
  geom_histogram(binwidth = 0.1, color = NA, alpha = 0.7) +
  scale_x_continuous(breaks = seq(0.1, 2, 0.2)) +
  labs(
    x = "Child-to-adult ratio",
    y = "# of households",
    title = "Child-to-adult ratio"
  ) +
  theme_bw() + my_theme

# ###############################
# Plot Figure 1D: Distribution of male-to-female ratio in the households
# ###############################
Figure_1D <- HH_metadata %>%
  ggplot(aes(x = male_to_female_ratio)) +
  geom_histogram(binwidth = 0.1, color = NA, alpha = 0.7) +
  scale_x_continuous(breaks = seq(0, 3, 0.5)) +
  labs(
    x = "Male-to-female ratio",
    y = "# of households",
    title = "Male-to-female ratio"
  ) +
  theme_bw() + my_theme

# ###############################
# Plot Figure 1E1: Numbers of:
#   weekly recruited households, 
#   outpatient COVID-19 cases at the recruitment site, 
#   BMR hospitalised COVID-19 cases, and
#   SARS-CoV-2 genome sequences reported from Thailand on the GISAID database 
# ###############################
timeline_lines_df <- data.frame(
  epiyear           = c(  2022,   2023,   2023,   2023,   2023,   2024,   2024,  2024),
  epiweek_start     = c(    45,      0,     14,     20,     45,      0,     14,    20),
  epiweek_end       = c(  52.5,     13,     17,     40,   52.5,     13,     17,    32),
  y                 = c(     2,      2,    2.5,      2,      2,      2,    2.5,     2),
  col               = c( "red",  "red", "blue",  "red",  "red",  "red", "blue", "red"),
  dataset           = rep("weekly_HH_enrolled", 8)
) %>% 
  mutate(dataset = factor(dataset, c("weekly_HH_enrolled", "weekly_opd_cases", "weekly_BMR_hosp_cases", "weekly_GISAID_seqs")))

timeline_lines_df
#   epiyear epiweek_start epiweek_end   y  col            dataset
# 1    2022            45        52.5 2.0  red weekly_HH_enrolled # Term 2 (End of year)
# 2    2023             0        13.0 2.0  red weekly_HH_enrolled # Term 2 (Start of year)
# 3    2023            14        17.0 2.5 blue weekly_HH_enrolled # April: Major summer break / Songkran
# 4    2023            20        40.0 2.0  red weekly_HH_enrolled # Term 1
# 5    2023            45        52.5 2.0  red weekly_HH_enrolled # Term 2 (End of year)
# 6    2024             0        13.0 2.0  red weekly_HH_enrolled # Term 2 (Start of year)
# 7    2024            14        17.0 2.5 blue weekly_HH_enrolled # April: Major summer break / Songkran
# 8    2024            20        32.0 2.0  red weekly_HH_enrolled # Term 1

timeline_lab_df <- data.frame(
  epiyear           = c(            2023,               2023,              2023),
  epiweek           = c(               0,                 15,                30),
  y                 = c(               2,                2.5,                 2),
  lab               = c( "School term-2", "Songkarn festival",  "School term-1"),
  dataset           = rep("weekly_HH_enrolled", 3)
) %>% 
  mutate(dataset = factor(dataset, c("weekly_HH_enrolled", "weekly_opd_cases", "weekly_BMR_hosp_cases", "weekly_GISAID_seqs")))

Figure_1E1 <- weekly_HH_case_seq_count_dat %>%
  pivot_longer(
    cols = weekly_HH_enrolled:weekly_GISAID_seqs,
    names_to = "dataset",
    values_to = "count"
  ) %>% 
  mutate(dataset = factor(dataset, c("weekly_HH_enrolled", "weekly_opd_cases", "weekly_BMR_hosp_cases", "weekly_GISAID_seqs"))) %>%
  ggplot(aes(x = epiweek, y = count)) +
  geom_col(position = "stack") + 
  
  # Add horizontal lines to indicate school terms and the national major holiday: Songkran
  geom_segment(
    data = timeline_lines_df,
    aes(x = epiweek_start, xend = epiweek_end, y = y, yend = y, color = col),
    linewidth = 0.5,
    inherit.aes = FALSE 
  ) + scale_color_identity() + # Manual colors for the lines

  # Add timeline label
  geom_text(
    data = timeline_lab_df,
    aes(x = epiweek, y = y, label = lab),
    size = 5/.pt, hjust = 0, vjust = -0.5,
    inherit.aes = FALSE 
  ) + 

  facet_grid(dataset ~ epiyear, 
             space = "free_x", 
             scales = "free", 
             labeller = labeller(
               epiyear = 
                 c("2022" = "Year\n2022", 
                   "2023" = "Year\n2023", 
                   "2024" = "Year\n2024"),
               dataset = 
                 c("weekly_HH_enrolled" = "Households\nenrolled into\nthis study", 
                   "weekly_opd_cases" = "OPD\ncases", 
                   "weekly_BMR_hosp_cases" = "BMR\nhospitalised\ncases", 
                   "weekly_GISAID_seqs" = "GISAID\nsequence\ncounts")
             ),
  ) +
  scale_x_continuous(
    name = "Epidemiological week",
    breaks = seq(5, 52, by = 10),
    expand = c(0, 0), 
    position = "bottom",
    labels = function(x) sprintf("w.%s", x),
  ) +
  scale_y_continuous(
    name = "Weekly count",
    expand = c(0, 0)
  ) + 
  theme_bw() + my_theme + 
  theme(
    axis.text.y = element_text(size = 6, vjust = 0),
    
    panel.spacing.x = unit(0,"line"),
    panel.background = element_blank(),
  ) + 
  coord_cartesian(clip = "off")

# ###############################
# Plot Figure 1E2: Complete pairwise lm fitting between
#   weekly recruited households, 
#   outpatient COVID-19 cases at the recruitment site, 
#   BMR hospitalised COVID-19 cases, and
#   SARS-CoV-2 genome sequences reported from Thailand on the GISAID database 
# ###############################
# Upper triangle: display equation, R2, correlation coefficient, and p-value
# -------------------------------
upper_tri_plot_fun <- function(data, mapping) {
  # Extract x and y values from the data
  x <- eval_data_col(data, mapping$x)
  y <- eval_data_col(data, mapping$y)
  
  # Fit linear model forced through (0,0)
  mod <- lm(y ~ x + 0)
  
  # Compute Pearson correlation test
  ct <- cor.test(x, y)
  
  # Create the label string with equation, R2, correlation and p-value
  label <- sprintf("Y = %.3f * X\nadj R² = %.3f\nCorr = %.3f\np = %.3g", coef(mod)[1], summary(mod)$adj.r.squared, ct$estimate, ct$p.value)
  
  # Create the plot
  ggplot(data = data, mapping = mapping) +
    annotate("text", x = 0.5, y = 0.5, label = label, size = 6/(72.27/25.4), hjust = 0.5, vjust = 0.5) +
    theme_void() + 
    theme(
      panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5),
      panel.background = element_blank(),
    )
}

# Diagonal: display text labels
# -------------------------------
diag_plot_fun <- function(data, mapping) {
  var <- as_label(mapping$x)
  labels <- c(
    weekly_HH_enrolled = "Households\nenrolled into\nthis study",
    weekly_opd_cases = "OPD\ncases",
    weekly_BMR_hosp_cases = "BMR\nhospitalised\ncases",
    weekly_GISAID_seqs = "GISAID\nsequence\ncounts"
  )
  label <- labels[[var]]
  
  ggplot(data = data, mapping = mapping) +
    annotate("text", x = 0.5, y = 0.5, label = label, size = 6/(72.27/25.4), hjust = 0.5, vjust = 0.5) +
    theme_void() + 
    theme(
      panel.background = element_rect(fill = "grey85", colour = "black", linewidth = 0.5)
    )
}

# Lower triangle: scatter plot with trend line forced through (0,0)
# -------------------------------
lower_tri_plot_fun <- function(data, mapping) {
  ggplot(data = data, mapping = mapping) + 
    geom_point(size = 0.5) +
    geom_smooth(method = "lm", formula = y ~ 0 + x, se = F) +
    theme_bw() +
    theme(
      panel.background = element_blank(),
      axis.text = element_text(size = 6)
    )
}

# Plot everything together
# -------------------------------
Figure_1E2 <- weekly_HH_case_seq_count_dat %>%
  select(weekly_HH_enrolled, weekly_opd_cases, weekly_BMR_hosp_cases, weekly_GISAID_seqs) %>%
  ggpairs(
    upper = list(continuous = wrap(upper_tri_plot_fun)),
    diag = list(continuous = wrap(diag_plot_fun)),
    lower =list(continuous = wrap(lower_tri_plot_fun)), 
    title = "Correlations between weekly counts",
    xlab = "Count",
    columnLabels = c(
      "Households enrolled\ninto this study", 
      "OPD\ncases", 
      "BMR\nhospitalised cases", 
      "GISAID\nsequence count"
    )
  ) +
  theme(
    plot.title = element_text(size = 8, face = "bold"),
    plot.margin = unit(c(0, 0.05, 0, 0.05), "cm"),
    
    axis.title = element_text(size = 6),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    
    strip.text = element_blank(),
    
  ) + coord_cartesian(clip = "off")

# ###############################
# Plot Figure 1F: Participant age distribution with COVID-19 status at enrolment
# ###############################
Figure_1F <- participant_metadata %>%
  ggplot(aes(x = age, fill = overall_res_enrol)) +
  geom_histogram(breaks = seq(0,80,3), color = NA) +
  geom_vline(xintercept = 18, linetype = "solid") +
  geom_text(x = 18, y = 0, label="Children", size = 6/.pt, hjust = 1, vjust = -0.5) +
  geom_text(x = 80, y = 0, label="Adults", size = 6/.pt, hjust = 1, vjust = -0.5) +
  labs(
    x = "Age (years)",
    y = "# of participants",
    title = "Age\ndistribution"
  ) +
  scale_fill_manual(values = overall_res_enrol_cols,
                    breaks = names(overall_res_enrol_cols),
                    labels = c("Case", "Assumed at risk", "NA"),
                    name = "Case status at enrolment"
  ) +
  theme_bw() + my_theme + guides(fill = guide_legend(nrow = 2, byrow = FALSE))

# ###############################
# Plot Figure 1G: Distribution of Vaccine doses received with COVID-19 status at enrolment
# ###############################
Figure_1G <- participant_metadata %>%
  ggplot() +
  geom_bar(aes(x = round(num_vac), fill = overall_res_enrol), position = "stack", color = NA) +
  scale_x_continuous(breaks = seq(0, 6, 1)) +
  labs(
    x = "Vaccine doses\nreceived",
    y = "# of participants\nat enrolment",
    title = "Vaccine doses\nreceived"
  ) +
  scale_fill_manual(values = overall_res_enrol_cols) +
  theme_bw() + my_theme + theme(legend.position = "none")

# ###############################
# Plot Figure 1H: Distribution of time to last vaccination with COVID-19 status at enrolment
# ###############################
Figure_1H <- participant_metadata %>%
  ggplot() +
  geom_bar(aes(x = time_to_last_vac_cat, fill = overall_res_enrol), position = "stack", color = NA) +
  scale_y_continuous(breaks = seq(0, 140, 20)) +
  labs(
    x = "Time to last vaccination",
    y = "# of participants\nat enrolment",
    title = "Time to last\nvaccination"
  ) +
  scale_fill_manual(values = overall_res_enrol_cols) +
  theme_bw() + my_theme + theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
    legend.position = "none"
    )

# ###############################
# Plot Figure 1I: Distribution of previous infections with COVID-19 status at enrolment
# ###############################
Figure_1I <- participant_metadata %>%
  ggplot(aes(x = num_prev_COVID, fill = overall_res_enrol)) +
  geom_bar(position = "stack", color = NA) +
  scale_y_continuous(breaks = seq(0, 140, 20)) +
  labs(
    x = "Times",
    y = "# of participants\nat enrolment",
    title = "Previous\nCOVID-19"
  ) +
  scale_fill_manual(values = overall_res_enrol_cols) +
  theme_bw() + my_theme + theme(legend.position = "none")

# ###############################
# Plot Figure 1J: Distribution of time to last prev infection with COVID-19 status at enrolment
# ###############################
Figure_1J <- participant_metadata %>%
  ggplot() +
  geom_bar(aes(x = time_to_last_prev_COVID_cat , fill = overall_res_enrol), position = "stack", color = NA) +
  scale_y_continuous(breaks = seq(0, 140, 20)) +
  labs(
    x = "Time to last COVID-19",
    y = "# of participants\nat enrolment",
    title = "Time to last\nCOVID-19"
  ) +
  scale_fill_manual(values = overall_res_enrol_cols) +
  theme_bw() + my_theme + theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
    legend.position = "none"
  )

# ###############################
# Plot Figure 1K: Distribution of SARS-CoV-2 positive HH members at the time of enrolment
# ###############################
Figure_1K <- HH_metadata %>%
  ggplot(aes(x = prop_pos_enrol)) +
  geom_histogram(binwidth = 0.1, color = NA, alpha = 0.7) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(breaks = seq(0, 15, 5)) +
  labs(
    x = "Case proportion\nat enrolment",
    y = "# of households",
    title = "Case\nproprotion"
  ) +
  theme_bw() + my_theme

# ###############################
# Plot Figure 1L: Participant age distribution with COVID-19 status at endpoint
# ###############################
Figure_1L <- participant_metadata %>%
  ggplot(aes(x = age, fill = overall_res_end_with_symp_severity)) +
  geom_histogram(breaks = seq(0,80,3), color = NA) +
  geom_vline(xintercept = 18, linetype = "solid") +
  geom_text(x = 18, y = 0, label="Children", size = 6/.pt, hjust = 1, vjust = -0.5) +
  geom_text(x = 80, y = 0, label="Adults", size = 6/.pt, hjust = 1, vjust = -0.5) +
  labs(
    x = "Age (years)",
    y = "# of participants",
    title = "Age\ndistribution"
  ) +
  scale_fill_manual(values = overall_res_end_with_symp_severity_cols,
                    breaks = names(overall_res_end_with_symp_severity_cols),
                    labels = c(
                      "Case: no symp data",
                      "Case: highly symp",
                      "Case: mod symp",
                      "Case: min symp",
                      "Uninf. contact",
                      "No info."
                    ),
                    name = "Case status by the end of follow-up"
  ) +
  theme_bw() + my_theme + guides(fill = guide_legend(nrow = 4, byrow = FALSE))

# ###############################
# Plot Figure 1M: Distribution of vaccine doses received with COVID-19 status at endpoint
# ###############################
Figure_1M <- participant_metadata %>%
  ggplot() +
  geom_bar(aes(x = round(num_vac), fill = overall_res_end_with_symp_severity), position = "stack", color = NA) +
  scale_x_continuous(breaks = seq(0, 6, 1)) +
  labs(
    x = "Vaccine doses\nreceived",
    y = "# of participants\nat endpoint",
    title = "Vaccine doses\nreceived"
  ) +
  scale_fill_manual(values = overall_res_end_with_symp_severity_cols) +
  theme_bw() + my_theme + theme(legend.position = "none")

# ###############################
# Plot Figure 1N: Distribution of time to last vaccination with COVID-19 status at endpoint
# ###############################
Figure_1N <- participant_metadata %>%
  ggplot() +
  geom_bar(aes(x = time_to_last_vac_cat, fill = overall_res_end_with_symp_severity), position = "stack", color = NA) +
  scale_y_continuous(breaks = seq(0, 140, 20)) +
  labs(
    x = "Time to last vaccination",
    y = "# of participants\nat enrolment",
    title = "Time to last\nvaccination"
  ) +
  scale_fill_manual(values = overall_res_end_with_symp_severity_cols) +
  theme_bw() + my_theme + theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
    legend.position = "none"
  )

# ###############################
# Plot Figure 1O: Distribution of previous infections with COVID-19 status at endpoint
# ###############################
Figure_1O <- participant_metadata %>%
  ggplot(aes(x = num_prev_COVID, fill = overall_res_end_with_symp_severity)) +
  geom_bar(position = "stack", color = NA) +
  scale_y_continuous(breaks = seq(0, 140, 20)) +
  labs(
    x = "Times",
    y = "# of participants\nat endpoint",
    title = "Previous\nCOVID-19"
  ) +
  scale_fill_manual(values = overall_res_end_with_symp_severity_cols) +
  theme_bw() + my_theme + theme(legend.position = "none")

# ###############################
# Plot Figure 1P: Distribution of time to last prev infection with COVID-19 status at endpoint
# ###############################
Figure_1P <- participant_metadata %>%
  ggplot() +
  geom_bar(aes(x = time_to_last_prev_COVID_cat , fill = overall_res_end_with_symp_severity), position = "stack", color = NA) +
  scale_y_continuous(breaks = seq(0, 140, 20)) +
  labs(
    x = "Time to last COVID-19",
    y = "# of participants\nat enrolment",
    title = "Time to last\nCOVID-19"
  ) +
  scale_fill_manual(values = overall_res_end_with_symp_severity_cols) +
  theme_bw() + my_theme + theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1),
    legend.position = "none"
  )

# ###############################
# Plot Figure 1Q: Distribution of SARS-CoV-2 positive HH members at endpoint
# ###############################
Figure_1Q <- HH_metadata %>%
  ggplot(aes(x = prop_pos_end)) +
  geom_histogram(binwidth = 0.1, color = NA, alpha = 0.7) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(breaks = seq(0, 15, 5)) +
  labs(
    x = "Case proportion\nby endpoint",
    y = "# of households",
    title = "Case\nproprotion"
  ) +
  theme_bw() + my_theme

# ###############################
# Make Figure 1 and save it to file
# ###############################
left_aligned_plots <- cowplot::align_plots(
  Figure_1E1,
  Figure_1F,
  Figure_1L,
  align = "v", axis = "l")

Figure_1E1_aln <- left_aligned_plots[[1]]
Figure_1F_aln <- left_aligned_plots[[2]]
Figure_1L_aln <- left_aligned_plots[[3]]

right_aligned_plots <- cowplot::align_plots(
  Figure_1B, 
  Figure_1C,
  Figure_1D,
  ggmatrix_gtable(Figure_1E2),
  Figure_1K,
  Figure_1Q,
  align = "v", axis = "r")

Figure_1B_aln <- right_aligned_plots[[1]]
Figure_1C_aln <- right_aligned_plots[[2]]
Figure_1D_aln <- right_aligned_plots[[3]]
Figure_1E2_aln <- right_aligned_plots[[4]]
Figure_1K_aln <- right_aligned_plots[[5]]
Figure_1Q_aln <- right_aligned_plots[[6]]

top_aligned_plots <- cowplot::align_plots(
  Figure_1A, 
  Figure_1B_aln,
  align = "h", axis = "t")

Figure_1A_aln <- top_aligned_plots[[1]]
Figure_1B_aln <- top_aligned_plots[[2]]

btm_aligned_plots <- cowplot::align_plots(
  Figure_1A_aln, 
  Figure_1D_aln,
  align = "h", axis = "b")

Figure_1A_aln <- btm_aligned_plots[[1]]
Figure_1D_aln <- btm_aligned_plots[[2]]

Figure_1B_1C_1D <- plot_grid(
  Figure_1B_aln, 
  Figure_1C_aln, 
  Figure_1D_aln,
  labels = c("B)", "C)", "D)"), label_size = 8,
  label_x = 0, label_y = 1, hjust = 0, vjust = 1,
  nrow = 3, align = "v", axis = "lr", rel_heights = c(1, 1, 1)
)

Figure_1A_1B_1C_1D <- plot_grid(
  Figure_1A_aln, Figure_1B_1C_1D, 
  labels = c("A)", ""), label_size = 8,
  label_x = 0, label_y = 1, hjust = 0, vjust = 1,
  nrow = 1, rel_widths = c(3, 1)
) 

Figure_1E1_1E2 <- plot_grid(
  Figure_1E1_aln, Figure_1E2_aln,
  labels = c("E)", ""), label_size = 8,
  label_x = 0, label_y = 1, hjust = 0, vjust = 1,
  nrow = 1, align = "h", axis = "bt", rel_widths = c(1.4, 1)
)

Figure_1F_1G_1H_1I_1J_1K <- plot_grid(
  Figure_1F_aln,
  (Figure_1G + theme(axis.title.y = element_blank())),
  (Figure_1H + theme(axis.title.y = element_blank())),
  (Figure_1I + theme(axis.title.y = element_blank())),
  (Figure_1J + theme(axis.title.y = element_blank())),
  Figure_1K_aln, 
  labels = c("F)", "G)", "H)", "I)", "J)", "K)"), label_size = 8,
  label_x = 0, label_y = 1, hjust = 0, vjust = 1,
  nrow = 1, align = "h", axis = "bt", rel_widths = c(3.8, 2, 1.7, 1.5, 1.7, 1.8)
) 

Figure_1L_1M_1N_1O_1P_1Q <- plot_grid(
  Figure_1L_aln,
  (Figure_1M + theme(axis.title.y = element_blank())),
  (Figure_1N + theme(axis.title.y = element_blank())),
  (Figure_1O + theme(axis.title.y = element_blank())),
  (Figure_1P + theme(axis.title.y = element_blank())),
  Figure_1Q_aln, 
  labels = c("L)", "M)", "N)", "O)", "P)", "Q)"), label_size = 8,
  label_x = 0, label_y = 1, hjust = 0, vjust = 1,
  nrow = 1, align = "h", axis = "bt", rel_widths = c(3.8, 2, 1.7, 1.5, 1.7, 1.8)
) 

Figure_1 <- plot_grid(
  Figure_1A_1B_1C_1D,
  Figure_1E1_1E2,
  Figure_1F_1G_1H_1I_1J_1K,
  Figure_1L_1M_1N_1O_1P_1Q,
  labels = c("", "", "", ""),
  label_x = 0, label_y = 1, hjust = 0, vjust = 1,
  ncol = 1, rel_heights = c(1, 1, 0.7, 0.7)
) 

ggsave(filename = Figure_1_file_png, device = "png",
       plot = Figure_1, 
       width = 16, height = 23, units = "cm", dpi = 300, bg = "white")

ggsave(filename = Figure_1_file_svg, device = "svg",
       plot = Figure_1, 
       width = 16, height = 23, units = "cm", dpi = 300, bg = "white")

# ###############################
# Compute miscellaneous stats
# ###############################
# Main stats
# -------------------------------
# Number of child and adult participants by sex
participant_metadata %>% select(age_cat, sex) %>% 
  table(useNA = "ifany") %>% addmargins(FUN = sum) 
# Margins computed over dimensions
# in the following order:
# 1: age_cat
# 2: sex
#        sex
# age_cat Female Male sum
#   Adult    102   72 174
#   Child     47   45  92
#   sum      149  117 266

# Stats related to HH sizes
HH_metadata %>% 
  summarise(
    avg_HH_size = mean(HH_size, na.rm = TRUE), 
    min_HH_size = min(HH_size, na.rm = TRUE),
    max_HH_size = max(HH_size, na.rm = TRUE),
    sd_HH_size = sd(HH_size, na.rm = TRUE),
  )
#   avg_HH_size min_HH_size max_HH_size sd_HH_size
# 1    4.433333           2           9   1.454003

# Stats related to HH sizes: children
HH_metadata %>% 
  summarise(
    avg_num_children = mean(num_children, na.rm = TRUE), 
    min_num_children = min(num_children, na.rm = TRUE),
    max_num_children = max(num_children, na.rm = TRUE),
    sd_num_children = sd(num_children, na.rm = TRUE),
  )
#   avg_num_children min_num_children max_num_children sd_num_children
# 1         1.533333                1                4       0.7002824

# Stats related to HH sizes: adults
HH_metadata %>% 
  summarise(
    avg_num_adults = mean(num_adults, na.rm = TRUE), 
    min_num_adults = min(num_adults, na.rm = TRUE),
    max_num_adults = max(num_adults, na.rm = TRUE),
    sd_num_adults = sd(num_adults, na.rm = TRUE),
  )
#   avg_num_adults min_num_adults max_num_adults sd_num_adults
# 1            2.9              1              6      1.115378

# Stats related to HH sizes: children:adult ratio
HH_metadata %>%
  summarise(
    avg_child_to_adult_ratio = mean(child_to_adult_ratio, na.rm = TRUE), 
    min_child_to_adult_ratio = min(child_to_adult_ratio, na.rm = TRUE),
    max_child_to_adult_ratio = max(child_to_adult_ratio, na.rm = TRUE),
    sd_child_to_adult_ratio = sd(child_to_adult_ratio, na.rm = TRUE),
  )
#   avg_child_to_adult_ratio min_child_to_adult_ratio max_child_to_adult_ratio sd_child_to_adult_ratio
# 1                0.5844444                      0.2                      1.5               0.2946776

# Stats related to HH sizes: male:female ratio
HH_metadata %>%
  summarise(
    avg_male_to_female_ratio = mean(male_to_female_ratio, na.rm = TRUE), 
    min_male_to_female_ratio = min(male_to_female_ratio, na.rm = TRUE),
    max_male_to_female_ratio = max(male_to_female_ratio, na.rm = TRUE),
    sd_male_to_female_ratio = sd(male_to_female_ratio, na.rm = TRUE),
  )
#   avg_male_to_female_ratio min_male_to_female_ratio max_male_to_female_ratio sd_male_to_female_ratio
# 1                 1.003333                        0                        3               0.7434637

# Number of child and adult participants by baseline COVID-19 status at enrolment
# Pos == "Participants who tested positive for SARS-CoV-2 within 14 days before enrolment and/or at enrolment--i.e., had a positive test result on or before the enrolment date"
# Neg == "Otherwise--i.e., all remaining household members were considered household members at risk at enrolment"
participant_metadata %>% 
  select(age_cat, overall_res_enrol) %>% 
  table(useNA = "ifany") %>% addmargins(FUN = sum) 
# Margins computed over dimensions
# in the following order:
# 1: age_cat
# 2: overall_res_enrol
#        overall_res_enrol
# age_cat Pos Neg Unconfirmed sum
#   Adult  50 124           0 174
#   Child  68  24           0  92
#   sum   118 148           0 266

# Number of child and adult participants by baseline COVID-19 status at enrolment -- stricter classification
# Pos == "Participants who tested positive for SARS-CoV-2 within 14 days before enrolment and/or at enrolment--i.e., had a positive test result on or before the enrolment date"
# Neg == "Participants who had a negative test result within two days after enrolment and no prior positive result"
# Unconfirmed == "Otherwise"
participant_metadata %>% 
  select(age_cat, overall_res_enrol_strict) %>% 
  table(useNA = "ifany") %>% addmargins(FUN = sum) 
# Margins computed over dimensions
# in the following order:
# 1: age_cat
# 2: overall_res_enrol_strict
#        overall_res_enrol_strict
# age_cat Pos Neg Unconfirmed sum
#   Adult  50 110          14 174
#   Child  68  24           0  92
#   sum   118 134          14 266

# Stats related to baseline HH case rate
HH_metadata %>%
  summarise(
    avg_prop_pos_enrol = mean(prop_pos_enrol, na.rm = TRUE), 
    min_prop_pos_enrol = min(prop_pos_enrol, na.rm = TRUE),
    max_prop_pos_enrol = max(prop_pos_enrol, na.rm = TRUE),
    sd_prop_pos_enrol = sd(prop_pos_enrol, na.rm = TRUE),
  )
#   avg_prop_pos_enrol min_prop_pos_enrol max_prop_pos_enrol sd_prop_pos_enrol
# 1          0.4443651          0.1666667          0.8333333         0.1920013

# Stats related to baseline HH case number
HH_metadata %>%
  summarise(
    avg_num_pos_enrol = mean(num_pos_enrol, na.rm = TRUE), 
    min_num_pos_enrol = min(num_pos_enrol, na.rm = TRUE),
    max_num_pos_enrol = max(num_pos_enrol, na.rm = TRUE),
    sd_num_pos_enrol = sd(num_pos_enrol, na.rm = TRUE),
  )
#   avg_num_pos_enrol min_num_pos_enrol max_num_pos_enrol sd_num_pos_enrol
# 1          1.966667                 1                 5         1.104178

# Number of child and adult participants by participant group defined wrt index case
participant_metadata %>% 
  select(age_cat, group_wrt_index_case) %>% 
  table(useNA = "ifany") %>% addmargins(FUN = sum) 
# Margins computed over dimensions
# in the following order:
# 1: age_cat
# 2: group_wrt_index_case
#        group_wrt_index_case
# age_cat Index case Co-index case Pre-index case / additional case identified at enrolment Post-index case / incident case Uninfected contact Unconfirmed status sum
#   Adult          0            10                                                       40                              40                 83                  1 174
#   Child         64             1                                                        3                              12                 12                  0  92
#   sum           64            11                                                       43                              52                 95                  1 266

# Number of child and adult cases by symp lv
participant_metadata %>% 
  filter(overall_res_end == "Pos") %>%
  select(age_cat, symp_cluster_res) %>% 
  table(useNA = "ifany") %>% addmargins(FUN = sum) 
# Margins computed over dimensions
# in the following order:
# 1: age_cat
# 2: symp_cluster_res
#        symp_cluster_res
# age_cat Minimally/(a)symptomatic Moderately symptomatic Highly symptomatic <NA> sum
#   Adult                        5                     13                 25   47  90
#   Child                       11                     29                 40    0  80
#   sum                         16                     42                 65   47 170

# Stats related to endpoint HH case rate
HH_metadata %>%
  summarise(
    avg_prop_pos_end = mean(prop_pos_end, na.rm = TRUE), 
    min_prop_pos_end = min(prop_pos_end, na.rm = TRUE),
    max_prop_pos_end = max(prop_pos_end, na.rm = TRUE),
    sd_prop_pos_end = sd(prop_pos_end, na.rm = TRUE),
  )
#   avg_prop_pos_end min_prop_pos_end max_prop_pos_end sd_prop_pos_end
# 1        0.6286839        0.1666667                1        0.281717

# Stats related to endpoint HH case number
HH_metadata %>%
  summarise(
    avg_num_pos_end = mean(num_pos_end, na.rm = TRUE), 
    min_num_pos_end = min(num_pos_end, na.rm = TRUE),
    max_num_pos_end = max(num_pos_end, na.rm = TRUE),
    sd_num_pos_end = sd(num_pos_end, na.rm = TRUE),
  )
#   avg_num_pos_end min_num_pos_end max_num_pos_end sd_num_pos_end
# 1        2.833333               1               8       1.728786

# Number of child and adult cases by sequence availability
participant_metadata %>% 
  filter(overall_res_end == "Pos") %>%
  mutate(seq_aval = !is.na(seq_id)) %>%
  select(age_cat, seq_aval) %>% 
  table(useNA = "ifany") %>% addmargins(FUN = sum) 
# Margins computed over dimensions
# in the following order:
# 1: age_cat
# 2: seq_aval
#        seq_aval
# age_cat FALSE TRUE sum
#   Adult    25   65  90
#   Child    17   63  80
#   sum      42  128 170

# Number of HHs by sequence availability
HH_metadata %>% 
  mutate(seq_aval = prop_seq_avail!=0) %>% 
  select(seq_aval) %>% 
  table(useNA = "ifany") %>% addmargins()
# seq_aval
# FALSE  TRUE   Sum 
#     6    54    60 

# Stats related to HH sequence number
HH_metadata %>%
  filter(num_seq_avail != 0) %>%
  summarise(
    avg_num_seq_avail = mean(num_seq_avail, na.rm = TRUE), 
    min_num_seq_avail = min(num_seq_avail, na.rm = TRUE),
    max_num_seq_avail = max(num_seq_avail, na.rm = TRUE),
    sd_num_seq_avail = sd(num_seq_avail, na.rm = TRUE),
  )
#   avg_num_seq_avail min_num_seq_avail max_num_seq_avail sd_num_seq_avail
# 1           2.37037                 1                 7         1.483193

# Number of HHs with some sequences by complete sequence status
HH_metadata %>% 
  filter(prop_seq_avail!=0) %>%
  mutate(complete_seq_aval = prop_seq_avail==1) %>% 
  select(complete_seq_aval) %>% 
  table(useNA = "ifany") %>% addmargins()
# complete_seq_aval
# FALSE  TRUE   Sum 
#    25    29    54 

# Virus lineage dist
participant_metadata %>% 
  filter(!is.na(seq_id)) %>%
  select(major_lineage) %>% 
  table(useNA = "ifany") %>% addmargins() 
# major_lineage
#       BA.5    BA.2.75        XBB        XBL XBB or XBL    BA.2.86     Others        Sum 
#          5         30         57          2          0         34          0        128 

# Additional stats
# -------------------------------
# Stats related to vac doses
participant_metadata %>%
  summarise(
    avg_num_vac = mean(num_vac, na.rm = TRUE), 
    min_num_vac = min(num_vac, na.rm = TRUE),
    max_num_vac = max(num_vac, na.rm = TRUE),
    sd_num_vac = sd(num_vac, na.rm = TRUE),
  )
#   avg_num_vac min_num_vac max_num_vac sd_num_vac
# 1    2.853383           0           6   1.325059

# Stats related to the number of prev inf
participant_metadata %>%
  summarise(
    avg_num_prev_COVID = mean(num_prev_COVID, na.rm = TRUE), 
    min_num_prev_COVID = min(num_prev_COVID, na.rm = TRUE),
    max_num_prev_COVID = max(num_prev_COVID, na.rm = TRUE),
    sd_num_prev_COVID = sd(num_prev_COVID, na.rm = TRUE),
  )
#   avg_num_prev_COVID min_num_prev_COVID max_num_prev_COVID sd_num_prev_COVID
# 1          0.6165414                  0                  2         0.6108677

# Household type dist
HH_metadata %>% 
  select(first_HH_case_age_cat) %>% 
  table(useNA = "ifany") %>% addmargins(FUN = sum)
# first_HH_case_age_cat
#         Adult Adult & child         Child           sum 
#            22             6            32            60 

