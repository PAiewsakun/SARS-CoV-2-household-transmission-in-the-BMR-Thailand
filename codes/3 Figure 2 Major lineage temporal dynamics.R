# ###############################
# Load global variables
# ###############################
source(file.path("miscellaneous.R"))

# ###############################
# Load required libraries
# ###############################
library(tidyverse)

library(ape) # read.nexus
library(ggtree)
library(ggnewscale) # new_scale_fill

library(ggplot2)
library(ggforce) # facet_wrap_paginate
library(ggpmisc) # stat_poly_eq
library(cowplot) # plot_grid

# ###############################
# Load data
# ###############################
# Sequence metadata
# -------------------------------
seq_metadata <- read.delim(seq_metadata_filepath, header = TRUE, sep = "\t") %>%
  mutate(
    collection_date = as.Date(collection_date, format = "%d/%b/%Y"),
    major_lineage = factor(major_lineage, names(major_lineage_cols)),
    region = factor(region, levels = names(BMR_cols))
  ) %>%
  filter(source == "This study")

# Household metadata
# -------------------------------
HH_metadata <- read.delim(HH_metadata_filepath, header = TRUE, sep = "\t") %>%
  mutate(major_lineage = factor(major_lineage, names(major_lineage_cols)))

# Sequence metadata
# -------------------------------
weekly_seq_count_by_lineage <- read.delim(weekly_seq_count_by_lineage_dat_filepath, header = TRUE, sep = "\t") %>%
  mutate(
    major_lineage = case_when(
      major_lineage %in% c("XBB", "XBL") ~ "XBB or XBL", # Combine XBB and XBL together
      TRUE ~ major_lineage
    ),
    major_lineage = factor(major_lineage, names(major_lineage_cols))
  ) %>%
  group_by(epiyear, epiweek, major_lineage) %>%
  mutate(
    weekly_seq_count_this_study = sum(weekly_seq_count_this_study, na.rm = F),
    weekly_seq_count_GISAID = sum(weekly_seq_count_GISAID, na.rm = F),
  )

# ###############################
# Plot Figure 2A: HH sequence count
# ###############################
seq_count_by_HH <- data.frame(seq_count_by_HH = seq_metadata %>% pull(family_id) %>% table() %>% as.vector)

Figure_2A <- seq_count_by_HH %>% 
  ggplot() +
  geom_histogram(aes(x = seq_count_by_HH), binwidth = 1, color = NA, alpha = 0.7) +
  scale_x_continuous(breaks = seq(0, 8, 1)) +
  labs(
    x = "Sequence count",
    y = "# of households",
  ) +
  theme_bw() + my_theme

# ###############################
# Plot Figure 2B: HH sequence availability
# ###############################
# Main plot
# -------------------------------
Figure_2B <- HH_metadata %>% 
  select(prop_seq_avail, major_lineage) %>%
  filter(!is.na(major_lineage)) %>%
  ggplot() +
  geom_histogram(aes(x = prop_seq_avail), binwidth = 0.1, color = NA, alpha = 0.7) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(breaks = seq(0, 13, 2)) +
  facet_grid(~ major_lineage) +
  labs(
    x = "Sequence availability (%)",
    y = "# of households",
  ) + theme_bw() + my_theme

# Pie charts to show households with 100% sequence availability
# -------------------------------
Figure_2B_inset_HH_with_full_seq_avail <- lapply(
  seq(5), 
  function(i) {
    HH_metadata %>% 
      select(prop_seq_avail, major_lineage) %>%
      filter(!is.na(major_lineage)) %>%
      group_by(major_lineage) %>%
      summarise(
        total = n(),
        num_full = sum(prop_seq_avail == 1.0, na.rm = T),
        prop_full = num_full/total,
        
        num_partial = sum(prop_seq_avail != 1.0, na.rm = T),
        prop_partial = num_partial/total,
        .groups = "drop"
      ) %>% 
      pivot_longer(
        cols = c(num_full, prop_full, num_partial, prop_partial),
        names_to = c(".value", "data_avail"),
        names_pattern = "^(num|prop)_(.*)$",
        values_drop_na = TRUE
      ) %>%
      mutate(data_avail = factor(data_avail, c("partial", "full"))) %>%
      ggplot(aes(x = "", y = prop, fill = data_avail)) +
      # Pie charts
      geom_bar(stat = "identity", width = 1, color = NA) +
      coord_polar(theta = "y") +
      # count labeling
      geom_text(aes(label = num), position = position_stack(vjust = 0.5), size = 2) +
      
      scale_fill_manual(values = c("full" = "gray", "partial" = "white")) +
      facet_wrap_paginate( ~ major_lineage, nrow = 1, ncol = 1, page = i) +
      theme_void() +
      theme(
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        
        axis.title = element_blank(),
        
        plot.background = element_blank()
      )
  } # end function
) # end lapply

Figure_2B_inset_HH_with_full_seq_avail <- tibble(
  x = 0,
  y = 1,
  plot = Figure_2B_inset_HH_with_full_seq_avail,
  major_lineage = factor(
    c("BA.5", "BA.2.75", "XBB", "XBL", "BA.2.86"), 
    levels = c("BA.5", "BA.2.75", "XBB", "XBL", "BA.2.86")
    ),
) 

# Combine them together
# -------------------------------
Figure_2B <- Figure_2B +
  geom_plot_npc(
    data = Figure_2B_inset_HH_with_full_seq_avail, 
    aes(npcx = x, npcy = y, label = plot, vp.width = 0.3, vp.height = 0.3)
  ) 

# ###############################
# Plot Figure 2C: GISAID SARS-CoV-2 sequences by major variants
# ###############################
Figure_2C <- weekly_seq_count_by_lineage %>%
  pivot_longer(
    cols = c(weekly_seq_count_this_study, weekly_seq_count_GISAID),
    names_to = "dataset",
    values_to = "count"
  ) %>% 
  mutate(
    dataset = factor(dataset, c("weekly_seq_count_this_study", "weekly_seq_count_GISAID"))
  ) %>%
  ggplot(aes(x = epiweek, y = count, fill = major_lineage)) +
  geom_col(position = "stack") +
  facet_grid(dataset ~ epiyear, 
             space = "free_x", 
             scales = "free", 
             labeller = labeller(
               dataset = 
                 c("weekly_seq_count_this_study" = "This study", 
                   "weekly_seq_count_GISAID" = "GISAID"),
               epiyear = 
                 c("2022" = "Year\n2022", 
                   "2023" = "Year\n2023", 
                   "2024" = "Year\n2024")
             ),
  ) +
  scale_x_continuous(
    name = "Epidemiological week",
    breaks = seq(0, 52, by = 10),
    expand = c(0, 0), 
    position = "bottom",
    labels = function(x) sprintf("w.%s", x),
  ) +
  scale_y_continuous(
    name = "Weekly count",
    expand = c(0, 0)
  ) +
  scale_fill_manual(
    values = major_lineage_cols,
    name = "Lineage"
  ) + 
  theme_bw() + my_theme +
  theme(
    legend.title = element_blank(), legend.position.inside = c(1, 0.45),
    axis.text.y = element_text(size = 6, vjust = 0),
    panel.spacing.x = unit(0,"line"),
  ) + coord_cartesian(clip = "off")

# ###############################
# Plot Figure 2D: Correlation of weekly sequence counts by major lineage
# ###############################
Figure_2D <- weekly_seq_count_by_lineage %>%
  mutate(
    weekly_seq_count_this_study = replace_na(weekly_seq_count_this_study, 0),
    weekly_seq_count_GISAID = replace_na(weekly_seq_count_GISAID, 0)
  ) %>%
  filter(major_lineage %in% head(names(major_lineage_cols), -1)) %>%
  ggplot(aes(x = weekly_seq_count_GISAID, y = weekly_seq_count_this_study, col = major_lineage)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + 0, se = F) +
  facet_wrap(~ major_lineage, nrow = 2, scales = "free") +
  stat_poly_eq(
    use_label("eq"),
    formula = y ~ x + 0,
    parse = TRUE,
    label.x = 0.05, label.y = 0.9,
    size = 6/(72.27/25.4), col = "black"
  ) +
  stat_poly_eq(
    use_label("R2"),
    formula = y ~ x + 0,
    parse = TRUE,
    label.x = 0.05, label.y = 0.8,
    size = 6/(72.27/25.4), col = "black"
  ) +
  stat_poly_eq(
    use_label("P"),
    formula = y ~ x + 0,
    parse = TRUE,
    label.x = 0.05, label.y = 0.6,
    size = 6/(72.27/25.4), col = "black"
  ) +
  scale_x_continuous(
    name = "Weekly sequence count: GISAID",
    position = "bottom"
  ) +
  scale_y_continuous(
    name = "Weekly sequence count: Our cohort",
  ) +
  scale_color_manual(
    values = major_lineage_cols,
    name = "Lineage"
  ) + 
  theme_bw() + my_theme + 
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 6, vjust = 0),
    panel.spacing.x = unit(0,"line"),
  ) + coord_cartesian(clip = "off")

# ###############################
# Make Figure 2 and save it to file
# ###############################
left_aligned_plots <- cowplot::align_plots(
  Figure_2A, 
  Figure_2C,
  align = "v", axis = "l")

Figure_2A_aln <- left_aligned_plots[[1]]
Figure_2C_aln <- left_aligned_plots[[2]]

right_aligned_plots <- cowplot::align_plots(
  Figure_2B,
  Figure_2D,
  align = "v", axis = "r")

Figure_2B_aln <- right_aligned_plots[[1]]
Figure_2D_aln <- right_aligned_plots[[2]]

Figure_2A_2B <- plot_grid(
  Figure_2A_aln,
  Figure_2B_aln, 
  labels = c("A)", "B)"), label_size = 8, 
  label_x = 0, label_y = 1, hjust = 0, vjust = 1,
  nrow = 1, align = "h", axis = "bt", rel_widths = c(1, 5)
) 

Figure_2C_2D <- plot_grid(
  Figure_2C_aln,
  Figure_2D_aln, 
  labels = c("C)", "D)"), label_size = 8, 
  label_x = 0, label_y = 1, hjust = 0, vjust = 1,
  nrow = 1, align = "h", axis = "bt", rel_widths = c(2, 1)
) 

Figure_2 <- plot_grid(
  Figure_2A_2B,
  Figure_2C_2D, 
  labels = c("", ""),
  nrow = 2, 
  rel_heights = c(0.5, 1)
)

# Save it to file
ggsave(filename = Figure_2_file_png, device = "png",
       plot = Figure_2,
       width = 16, height = 8, units = "cm", dpi = 300, bg = "white")

ggsave(filename = Figure_2_file_svg, device = "svg",
       plot = Figure_2,
       width = 16, height = 8, units = "cm", dpi = 300, bg = "white")

# ###############################
# Compute some miscellaneous stats
# ###############################
# sequence count by HH
# -------------------------------
seq_count_by_HH %>% 
  table(useNA = "ifany") %>% addmargins()
# seq_count_by_HH
#   1   2   3   4   5   6   7 Sum 
#  20  12  13   4   2   2   1  54 

seq_count_by_HH %>%
  summarise(
    avg_seq_count_by_HH = mean(seq_count_by_HH, na.rm = TRUE), 
    min_seq_count_by_HH = min(seq_count_by_HH, na.rm = TRUE),
    max_seq_count_by_HH = max(seq_count_by_HH, na.rm = TRUE),
    sd_seq_count_by_HH = sd(seq_count_by_HH, na.rm = TRUE),
  )
#   avg_seq_count_by_HH min_seq_count_by_HH max_seq_count_by_HH sd_seq_count_by_HH
# 1             2.37037                   1                   7           1.483193
